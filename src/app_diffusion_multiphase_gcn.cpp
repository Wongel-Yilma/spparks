/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "solve.h"
#include "domain.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include "comm_lattice.h"
#include "app_diffusion_multiphase_gcn.h"
#include <vector>
#include "json.hpp"
#include <fstream>
#include <cmath>
#include <torch/torch.h>
#include <torch/script.h>
#include "matplotlibcpp.h"
using namespace SPPARKS_NS;
using json= nlohmann::json;
using std::map;
using std::set;

namespace plt = matplotlibcpp;

enum{LINEAR};
enum {ZERO, NI, H, V, T};   // $ types of sites Ni, H, V, Tetrahedral(Pinned vacancy)
enum {ZERO_TYPE, FCC, OCTA, TETRA};

#define DELTAEVENT 100000
#define NUHOP 1e13
#define TOLERANCE 1.0e-4

// This app is based on the diffusion app (for Kawasaki dynamics)
// These are the significant changes and simplifications
// 1. No Schowebel hops, only exchanges to site neighbors
// 2. Allow for an arbitrary number of phases (species for an atomic model)
// 3. Only the linear energy style is supported
// 4. One or more phases can be pinned, meaning they do not diffuse

/* ---------------------------------------------------------------------- */

// AppDiffusionMultiphaseGCN::AppDiffusionMultiphaseGCN(SPPARKS *spk, int narg, char **arg) : 
//   AppLattice(spk,narg,arg), phase_labels(), is_pinned(), weights(), device(torch::cuda::is_available() ? torch::kCUDA : torch::kCPU)
// {
AppDiffusionMultiphaseGCN::AppDiffusionMultiphaseGCN(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg), phase_labels(), is_pinned(), weights()
{
  // need to double check these values
  torch::NoGradGuard no_grad;
  ninteger = 2;
  ndouble = 3;
  delpropensity = 2;
  delevent = 1;
  allow_kmc = 1;
  allow_rejection = 1;
  allow_masking = 0;
  numrandom = 1;
  neigh_check = NULL;
  hopsite = NULL;
  // Check if CUDA is available
  // if (device == torch::kCPU) {
  //   std::cout << "CUDA is not available. Using CPU." << std::endl;
  // } else {
  //   std::cout << "CUDA is available! Running Eb calculation on GPU." << std::endl;
  // }

  torch::set_num_threads(1);
  torch::set_num_interop_threads(1);

  // no args for this app

  if (narg < 6) error->all(FLERR,"Illegal app_style command");
  std::ifstream in(arg[1]);
  json j;
  in >> j;
  in.close();
  edge_index_vec = j["edge_index"];
  first_shell_coords = j["coords"];
  destinations = j["octa_destinations"];
  gcn = torch::jit::load(arg[2]);
  // gcn.to(device);
  batch =  torch::zeros({42}, torch::dtype(torch::kLong));
  edge_index_linearized = linearize_int(edge_index_vec);
  edge_index = torch::from_blob(edge_index_linearized.data(), {2,478}, torch::dtype(torch::kLong)).clone();
  // edge_index = edge_index.to(device);
  // batch = batch.to(device);
  box_dims_md  = new double[3];
  box_dims_md[0] = std::stof(arg[3]);
  box_dims_md[1] = std::stof(arg[4]);
  box_dims_md[2] = std::stof(arg[5]);

  engstyle = LINEAR;
  num_evaluations = 0;
  num_events_removed= 0;
  create_arrays();
  esites = NULL;
  echeck = NULL;
  maxevent = 0;
  events = NULL;
  firstevent = NULL;

  allocated = 0;
}

/* ---------------------------------------------------------------------- */

AppDiffusionMultiphaseGCN::~AppDiffusionMultiphaseGCN()
{
  std::string file_name = "event_stats"+std::to_string(me)+".png";
  // std::string file_name_eb = "eb"+std::to_string(me)+".png";

  plt::plot(active_events, "b-");
  plt::plot(removed_events, "r-");
  plt::plot(all_events_generated, "g-");
  // plt::show();
  plt::save(file_name); 
  // plt::plot(eb_values, "b-");
  // plt::save(file_name_eb);
  delete [] esites;
  delete [] echeck;
  delete [] neigh_check;
  memory->sfree(events);
  memory->destroy(firstevent);
}

/* ----------------------------------------------------------------------
   input script commands unique to this app
------------------------------------------------------------------------- */

void AppDiffusionMultiphaseGCN::input_app(char *command, int narg, char **arg)
{
  if (sites_exist == 0) {
    char str[128];
    sprintf(str,"Cannot use %s command until sites exist",command);
    error->all(FLERR,str);
  }

  if (!allocated) allocate_data();
  allocated = 1;

  if (strcmp(command,"diffusion/multiphase") == 0) 
    parse_diffmultiphase(narg,arg);
  else error->all(FLERR,"Unrecognized command");

}

/* ---------------------------------------------------------------------- */

void AppDiffusionMultiphaseGCN::parse_diffmultiphase(int narg, char **arg)
{
   // 2 args: diffusion/multiphase phase <int value>
   // 2 args: diffusion/multiphase pin <int value>
   // 5 args: diffusion/multiphase weight <double> pair <int,int>

   if (narg < 2) 
     error->all(FLERR,"Illegal diffusion/multiphase command");

   if (strcmp(arg[0],"phase")==0){
      if (narg != 2) 
        error->all(FLERR,"Illegal diffusion/multiphase phase command: "
                   "num args != 2");
      int phase = std::atoi(arg[1]);
      phase_labels.insert(phase);
      // phases are not pinned by default
      is_pinned[phase] = false;
      if (phase < 1) 
        error->all(FLERR,"Illegal diffusion/multiphase phase value: "
                   "must be >= 1");
   } else if (strcmp(arg[0],"pin") == 0) {
      if (narg != 2) 
        error->all(FLERR,"Illegal diffusion/multiphase pin command: "
                   "num args != 2");
      int pin_phase = std::atoi(arg[1]);
      phase_labels.insert(pin_phase);
      is_pinned[pin_phase] = true;
      if (pin_phase < 1)
        error->all(FLERR,"Illegal diffusion/multiphase pin value: must be >= 1");
   } else if (strcmp(arg[0],"weight")==0){
      if (narg != 5) 
        error->all(FLERR,"Illegal diffusion/multiphase weight command: "
                   "num args != 5");
      double w = std::atof(arg[1]);
      if (strcmp(arg[2],"pair") == 0) {
         int p1 = std::atoi(arg[3]);
         int p2 = std::atoi(arg[4]);
         if (p1 == p2) 
           error->all(FLERR,"Cannot set diffusion/multiphase weight for "
                      "pair of identical phases");
         if (p1 < 1 || p2 < 1) 
           error->all(FLERR,"Illegal diffusion/multiphase weight command: "
                      "phases must be >= 1");
         if (w < 0.0) 
           error->all(FLERR,"Illegal diffusion/multiphase weight command: "
                      "weight must be >= 0");
         weights[{p1,p2}] = w;
         weights[{p2,p1}] = w;
      } else 
        error->all(FLERR,"Illegal diffusion/multiphase weight command: "
                   "expected keyword pair");
   } else error->all(FLERR,"Illegal diffusion/multiphase command: "
                     "expected phase, pin, or weight");
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppDiffusionMultiphaseGCN::grow_app()
{
  lattice = iarray[0];
  type = iarray[1];
  x_md = darray[0];
  y_md = darray[1];
  z_md = darray[2];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppDiffusionMultiphaseGCN::init_app()
{
   if (!allocated) allocate_data();
   allocated = 1;

   dimension = domain->dimension;
   dt_sweep = 1.0/maxneigh;

   {
     // insure all site values are in set of phase labels, otherwise error

     std::set<int>::iterator not_found=phase_labels.end();
     int flag = 0;
     for (int i = 0; i < nlocal; i++)
       if (not_found == phase_labels.find(lattice[i])) flag=1;
     int flagall;
     MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
     if (flagall) error->all(FLERR,"One or more sites have invalid values");
   }

   {
     // create default weights if they do not already exist
     // only unlike phases contribute to energy

     std::map<std::pair<int,int>,double>::iterator not_found=weights.end();
     for (auto p : phase_labels) {
       for (auto q : phase_labels) {
         if (p==q) continue;
         std::map<std::pair<int,int>,double>::iterator i = weights.find({p,q});
         // if pair not found, set default weight = 1.0
         if (not_found==i) {
           weights[{p,q}] = 1.0;
           weights[{q,p}] = 1.0;
         }
       }
     }
   }
   
}

/* ----------------------------------------------------------------------
   setup before each run
------------------------------------------------------------------------- */

void AppDiffusionMultiphaseGCN::setup_app()
{
  for (int i = 0; i < nlocal+nghost; i++) echeck[i] = 0;

  // clear event list
  nevents = 0;
  for (int i = 0; i < nlocal; i++) firstevent[i] = -1;
  for (int i = 0; i < maxevent; i++) events[i].next = i+1;
  freeevent = 0;


  comm->all();

  std::vector<double>V_coords (3,0.0);
  int m;
  int num_sites = nlocal+nghost;
  int zero_type_count =0;
  // Loop over vacancy sites and set their darray values to their neighboring Ni sites
  for (int k = 0; k<num_sites; k++){
  
    if (iarray[0][k]==V){
      int num_Ni_neighbors =0;
    //   // std::cout<<"Rank "<< me << " processing vacancy site " << k << std::endl;
      // for (int kk=0; kk<numneigh[k];kk++){
        for (int kk=0; kk< 6; kk++){
          m = neighbor[k][kk];
          if (iarray[0][m]==ZERO){
            zero_type_count++;
            continue;
          }
          V_coords[0] += darray[0][m];
          if (xyz[m][0]-xyz[k][0]>half_box_dims[0]){
            V_coords[0] -= box_dims_md[0];
          }
          else if (xyz[m][0]-xyz[k][0]<-half_box_dims[0]){
            V_coords[0] += box_dims_md[0];
          }
          V_coords[1] += darray[1][m];
          if (xyz[m][1]-xyz[k][1]>half_box_dims[1]){
            V_coords[1] -= box_dims_md[1];
          }
          else if (xyz[m][1]-xyz[k][1]<-half_box_dims[1]){
            V_coords[1] += box_dims_md[1];
          }
          V_coords[2] += darray[2][m];
          if (xyz[m][2]-xyz[k][2]>half_box_dims[2]){
            V_coords[2] -= box_dims_md[2];
          }
          else if (xyz[m][2]-xyz[k][2]<-half_box_dims[2]){
            V_coords[2] += box_dims_md[2];
          }
        // }
      }
      darray[0][k] = V_coords[0]/6.0;
      darray[1][k] = V_coords[1]/6.0;
      darray[2][k] = V_coords[2]/6.0;
      for (int o=0; o<3; o++) V_coords[o]=0.0; // Setting back to zero for next vacancy site
    }
  }
  std::cout<< "Rank " << me << " has " << zero_type_count << " ZERO type sites "<< std::endl;
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppDiffusionMultiphaseGCN::site_energy(int i)
{
  // energy of site = linear sum of bond weights
  // only unlike phase pairs contribute

   double energy(0.0);
   int ip = lattice[i];
   for (int j = 0; j < numneigh[i]; j++){
      int nj = neighbor[i][j];
      int jp = lattice[nj];
      if (ip == jp) continue;
      energy += weights[{ip,jp}];
   }

   // each site carries half the interaction energy
   // neighbor sites carry the other half 
   return 0.5*energy;
}

/* ----------------------------------------------------------------------
   rKMC method
   perform a site event with null bin rejection
   null bin extends to size maxneigh
------------------------------------------------------------------------- */

void AppDiffusionMultiphaseGCN::site_event_rejection(int i, RandomPark *random)
{
  double einitial,edelta;
  int i_old, j_old;

  // pinned sites can't exchange

  if (is_pinned[lattice[i]]) return;

  // need to double check neighborhood

  int iran = (int) (maxneigh*random->uniform());
  if (iran > maxneigh) iran = maxneigh-1;
  int j = neighbor[i][iran];

  // if site j pinned or if site i and site j have same spin, just return

  if (is_pinned[lattice[j]] || lattice[j] == lattice[i]) return;

  i_old = lattice[i];
  j_old = lattice[j];
  
  // accept or reject via energy model

  int hop = 0;
  einitial = site_energy(i)+site_energy(j);

  lattice[i] = j_old;
  lattice[j] = i_old;

  // compute energy difference from exchange

  edelta = site_energy(i)+site_energy(j) - einitial;

  // if edelta is negative, accept
  // otherwise if temperature is non-zero, can still accept

  if (edelta <= 0.0) hop = 1;
  else if (temperature > 0.0) {
    if (random->uniform() < exp(-1.0*edelta*t_inverse)) hop = 1;
  }
    
  if (hop) {
    naccept++;
  } else {
    lattice[i] = i_old;
    lattice[j] = j_old;
  }
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppDiffusionMultiphaseGCN::site_propensity(int i)
{
  return site_propensity_linear(i);
}

/* ---------------------------------------------------------------------- */

double AppDiffusionMultiphaseGCN::site_propensity_linear(int i)
{
  int j,k;
  double probone,proball;
  
  // add event if neighbor sites are Vacancy sites
  
  clear_events(i);
  // if (lattice[i]==V){
  //   int num_Ni_neighbors =0;
  //   for (k = 0; k < numneigh[i]; k++) {
  //     j = neighbor[i][k];
  //     if (lattice[j]==NI) num_Ni_neighbors++;
  //   }
  //   std::cout<<"Rank "<< me << " processing H site " << i << " num sites " << nlocal +nghost<< std::endl;
  //   assert (num_Ni_neighbors==6), "H site does not have 6 Ni neighbors";
  // }
  
  // if (is_pinned[lattice[i]]) return 0.0;
  if (lattice[i]!=H) return 0.0;    // only allow phase '2' or 'H' to undergo diffusion jumps

  int nhop1 {0};
  int num_occupied_neighbors {0};
  int l, ll;
  int m, mm,o, ihop;
  std::vector<int> local_neigh_check(nlocal + nghost, 0);
  std::vector<std::vector<double>> neighbor_coords(200,std::vector<double>(3,0.0));
  std::vector<std::vector<double>> md_coords(200,std::vector<double>(3,0.0));
  std::vector<int> neighbor_types(200);
  std::vector<double> xi = {xyz[i][0], xyz[i][1], xyz[i][2]};   // Coordinate of the H atom
  std::vector<double> xmdi = {x_md[i], y_md[i], z_md[i]};
  double dist {0.0};
  double dist_md {0.0};



  // loop over all possible hops, go through neighbor shell

  // einitial = site_energy(i);
  proball = 0.0;
  probone = 0.0;
  // The central atom itself is included in the neighbor list (To formulate complete local graph)
  // Both its local and MD coordinates are zero.
  neighbor_types[num_occupied_neighbors++] = lattice[i];
  // this is similar to the Potts approach
  neigh_check[i]=1;
  local_neigh_check[i]=1;
  for (k = 0; k < numneigh[i]; k++) {
    j = neighbor[i][k];
  
    if (lattice[j]==V){// only allow exchanges with Vacancy
      hopsite [nhop1++] = j;
    } 
    if (lattice[j]==T) continue; // Tetrahedral sites are ignored
    else{  // Loop through the second neighbor and 
      for (o=0; o<3; o++){
        neighbor_coords[num_occupied_neighbors][o] = xyz[j][o] - xi[o];
      }
      md_coords[num_occupied_neighbors][0] = x_md[j]-xmdi[0];
      md_coords[num_occupied_neighbors][1] = y_md[j]-xmdi[1];
      md_coords[num_occupied_neighbors][2] = z_md[j]-xmdi[2];

      neighbor_types[num_occupied_neighbors++] = lattice[j];
      local_neigh_check[j]=1;

      for (l=0; l<numneigh[j]; l++){
        m = neighbor[j][l];
        if (lattice[m]!=T && local_neigh_check[m]==0){
          for (o=0; o<3; o++){
            neighbor_coords[num_occupied_neighbors][o] = xyz[m][o] - xi[o];
          }
          md_coords[num_occupied_neighbors][0] = x_md[m]-xmdi[0];
          md_coords[num_occupied_neighbors][1] = y_md[m]-xmdi[1];
          md_coords[num_occupied_neighbors][2] = z_md[m]-xmdi[2];
          neighbor_types[num_occupied_neighbors++] = lattice[m];
          local_neigh_check[m]=1;
          for (ll=0; ll<numneigh[m]; ll++){
            mm = neighbor[m][ll];
            if (lattice[mm]==NI && local_neigh_check[mm]==0){
              for (o=0; o<3; o++){
                neighbor_coords[num_occupied_neighbors][o] = xyz[mm][o] - xi[o];
              }
              md_coords[num_occupied_neighbors][0] = x_md[mm]-xmdi[0];
              md_coords[num_occupied_neighbors][1] = y_md[mm]-xmdi[1];
              md_coords[num_occupied_neighbors][2] = z_md[mm]-xmdi[2];
              neighbor_types[num_occupied_neighbors++] = lattice[mm];
              local_neigh_check[mm]=1;
            }
          }
        }
      }
    }
  }

  for(l=0; l<num_occupied_neighbors; l++){
    for(o = 0; o<3; o++){
      dist = neighbor_coords[l][o];
      dist_md = md_coords[l][o];
      if (dist>half_box_dims[o]) {
        dist-=box_dims[o];
        dist_md -=box_dims_md[o];
      }
      else if (dist<-half_box_dims[o]) {
      dist+=box_dims[o];
      dist_md +=box_dims_md[o];
      }
      neighbor_coords[l][o] = dist;
      md_coords[l][o] = dist_md;
    }
  }

  for(ihop= 0; ihop<nhop1; ihop++){
    j =  hopsite[ihop];
    double eb = calculate_barrier_energy(i,j,neighbor_coords,md_coords,neighbor_types,num_occupied_neighbors);
    probone = (NUHOP)*exp(-eb*t_inverse);
    add_event(i,j,probone);
    proball += probone;
    num_evaluations++;
    eb_values.push_back(eb);
  }
 for ( l = 0; l < nlocal + nghost; l++) local_neigh_check[l] = 0;
  return proball;
}

/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppDiffusionMultiphaseGCN::site_event(int i, class RandomPark *random)
{
  return site_event_linear(i,random);
}

/* ---------------------------------------------------------------------- */

void AppDiffusionMultiphaseGCN::site_event_linear(int i, class RandomPark *random)
{
  int j,k,m,isite,i_old,j_old;

  // pick one event from total propensity by accumulating its probability
  // compare prob to threshhold, break when reach it to select event
  // perform event
  // std::cout<< "Number of GCN evaluations: " << num_evaluations << std::endl;
  // std::cout<< "Number of Events: " << nevents << std::endl;
  // std::cout<< "Number of Events decreased: " << num_events_removed << std::endl;

  removed_events.push_back(num_events_removed);
  active_events.push_back(nevents);
  all_events_generated.push_back(num_evaluations);


  double threshhold = random->uniform() * propensity[i2site[i]];
  double proball = 0.0;

  int ievent = firstevent[i];
  while (1) {
    proball += events[ievent].propensity;
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
    if (ievent < 0) error->one(FLERR,"Did not reach event propensity threshhold");
  }

  // exchange values

  j = events[ievent].destination;
  i_old = lattice[i];
  j_old = lattice[j];
  lattice[i] = j_old;
  lattice[j] = i_old;

  // compute propensity changes for self and swap site and their neighs
  // 1,2 neighs for NNHOP and 1,2,3 neighs for SCHWOEBEL
  // ignore update of sites with isite < 0
  // use echeck[] to avoid resetting propensity of same site

  int nsites = 0;

  isite = i2site[i];
  propensity[isite] = site_propensity(i);
  esites[nsites++] = isite;
  echeck[isite] = 1;
  
  // update site i's neighbors, this will include the exchanged site

  for (k = 0; k < numneigh[i]; k++) {
    m = neighbor[i][k];
    isite = i2site[m];
    // not quite sure what this does
    if (isite < 0) continue;
    // add to update list
    esites[nsites++] = isite;
    propensity[isite] = site_propensity(m);
    echeck[isite] = 1;
  }
  
  // update exchanged site's neighbors
  // avoid any that have already been found

  for (k = 0; k < numneigh[j]; k++) {
    m = neighbor[j][k];
    isite = i2site[m];
    // not quite sure what this does
    if (isite < 0) continue;
    // make sure site is not already updated
    if (echeck[isite] == 1) continue;
    // add to update list
    esites[nsites++] = isite;
    propensity[isite] = site_propensity(m);
    echeck[isite] = 1;
  }
  
  solve->update(nsites,esites,propensity);

  // clear echeck array

  for (k = 0; k < nsites; k++) echeck[esites[k]] = 0;
}

double AppDiffusionMultiphaseGCN::calculate_distance(const std::vector<double> &a,const std::vector<double> &b){
  return std::sqrt(std::pow((a[0]-b[0]),2)+std::pow((a[1]-b[1]),2)+std::pow((a[2]-b[2]),2));
}

double AppDiffusionMultiphaseGCN::vec_dot(const std::vector<double> &a,const std::vector<double> &b){
  double c = a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
  return c;
}

std::vector<double> AppDiffusionMultiphaseGCN::vec_cross(const std::vector<double> &a,const std::vector<double> &b){
  std::vector<double> c(3,0.0);
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
  return c;
}

double AppDiffusionMultiphaseGCN::vec_norm(const std::vector<double> &a){
  double c = std::sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
  return c;
}
std::vector<double> AppDiffusionMultiphaseGCN::mat_vecmul(const std::vector<std::vector<double>> &B,const std::vector<double> &a){
  
  int m = a.size();
  int p = B.size();
  int q = B[0].size();
  std::vector<double> c(q,0.0);
  for (size_t i =0; i<p;i++){
      for(size_t k=0; k<m; k++){
          c[i] +=  B[k][i]*a[k];
      }
  }
  return c;
}
std::vector<int64_t> AppDiffusionMultiphaseGCN::linearize_int (const std::vector<std::vector<int64_t>> &A){
  int64_t num_rows = A.size();
  int64_t num_cols = A[0].size();
  std::vector<int64_t> B;
  B.reserve(num_rows * num_cols);
  for (size_t col_idx = 0; col_idx < num_cols; col_idx++){
    for (size_t row_idx = 0; row_idx < num_rows; row_idx++){
        B.push_back(A[row_idx][col_idx]);
    }
  }
  return B;
}

double AppDiffusionMultiphaseGCN::calculate_barrier_energy(int i, int j, std::vector<std::vector<double>>& local_coords,std::vector<std::vector<double>>& md_coords, std::vector<int>& neighbor_types,int num_occupied_sites){
  int o,l, k, kk,flag;
  int src, dst;
  std::vector<int64_t> first_shell(42,0);
  std::vector<double> new_x(3);
  std::vector<double> new_y(3);
  std::vector<double> new_z(3);
  std::vector<std::vector<double>> original_cs ={
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
      };

  std::vector<std::vector<double>> rotated_cs(3, std::vector<double>(3,0.0));
  std::vector<std::vector<double>> R(3, std::vector<double>(3,0.0));
  std::vector<std::vector<double>> rotated_coords(num_occupied_sites,std::vector<double>(3,0.0));
  std::vector<std::vector<double>> rotated_coords_md(num_occupied_sites,std::vector<double>(3,0.0));

  std::vector<std::vector<double>> selected_coords(42, std::vector<double>(3,0.0));
  std::vector<float> edge_attr_vec(edge_index_vec.size(),0.0);
  for ( o=0; o<3; o++){
    new_x[o] = xyz[j][o] - xyz[i][o];
    if (new_x[o]>half_box_dims[o]){
      new_x[o] = new_x[o]- box_dims[o] ;
    }
    else if (new_x[o]<-half_box_dims[o]){
      new_x[o] = new_x[o] + box_dims[o];
    }
  }

  double dot_product = 1.0;
  for ( l =0;l<destinations.size();l++){
    dot_product = vec_dot(new_x, destinations[l]);
    if (dot_product==0){
      new_y = destinations[l];
      break;
    }
  }
  rotated_cs[0]= new_x;
  rotated_cs[1]= new_y;
  new_z = vec_cross(new_x, new_y);
  rotated_cs[2]= new_z;

  for ( k=0; k<3; k++){
    for ( kk=0; kk<3; kk++){
        R[k][kk] = vec_dot(original_cs[k],rotated_cs[kk])/(vec_norm(original_cs[k])*vec_norm(rotated_cs[kk]));
    }
  }
  // Rotating the local coordinates with R
  for ( k = 0; k<num_occupied_sites;k++){
    rotated_coords[k] = mat_vecmul(R,local_coords[k]);
    rotated_coords_md[k] = mat_vecmul(R,md_coords[k]);
  }
 //Checking the occupancy of the local coordinates
  for (k =0; k<num_occupied_sites;k++){
    for (kk=0;kk<first_shell.size();kk++){
      flag=1;
      for(o=0; o<3;o++){
        if (fabs(first_shell_coords[kk][o]-rotated_coords[k][o])>TOLERANCE){
          flag=0;
          break;
        }
      }
      if(flag==1){
        first_shell[kk] = neighbor_types[k]-1; // Mapping Ni=0, H=1, V=2 -- as defined in the GCN 
        for (o=0; o<3;o++){
          selected_coords[kk][o] = rotated_coords_md[k][o];
        }
        break;
      }                                                                                      
    }
  }
  for (l=0;l<edge_index_vec.size(); l++){
    src = edge_index_vec[l][0];
    dst = edge_index_vec[l][1];
    edge_attr_vec[l] = calculate_distance(selected_coords[src], selected_coords[dst] ); // 3 is the number of site types
  }
  // Call the torch GCN model here
  torch::Tensor atom_type = torch::from_blob(first_shell.data(),{42}, torch::dtype(torch::kLong)).clone();
  torch::Tensor edge_attr = torch::from_blob(edge_attr_vec.data(),{478}, torch::dtype(torch::kFloat)).clone();
  // atom_type= atom_type.to(device);
  // edge_attr = edge_attr.to(device);

  std::vector<torch::jit::IValue> inputs = { atom_type ,edge_attr, batch,edge_index};
  // auto eb = gcn.forward(inputs).toTensor().to(torch::kCPU).item<double>();
  // Mean, std: 0.39268717 0.03553854
  auto eb = gcn.forward(inputs).toTensor().item<double>()*0.03553854+0.39268717;
  // double eb = 0.40; // Dummy value for debugging
  return eb;
}
/* ----------------------------------------------------------------------
   clear all events out of list for site I
   add cleared events to free list
------------------------------------------------------------------------- */

void AppDiffusionMultiphaseGCN::clear_events(int i)
{
  int next;
  int index = firstevent[i];
  while (index >= 0) {
    next = events[index].next;
    events[index].next = freeevent;
    freeevent = index;
    nevents--;
    index = next;
    num_events_removed++;
  }
  firstevent[i] = -1;
}

/* ----------------------------------------------------------------------
   add an event to list for site I
   event = exchange with site J with probability = propensity
------------------------------------------------------------------------- */

void AppDiffusionMultiphaseGCN::add_event(int i, int destination, 
			      double propensity)
{
  // grow event list and setup free list

  if (nevents == maxevent) {
    maxevent += DELTAEVENT;
    events = 
      (Event *) memory->srealloc(events,maxevent*sizeof(Event),"app:events");
    for (int m = nevents; m < maxevent; m++) events[m].next = m+1;
    freeevent = nevents;
  }

  int next = events[freeevent].next;
  events[freeevent].propensity = propensity;
  events[freeevent].destination = destination;
  events[freeevent].next = firstevent[i];
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}

/* ----------------------------------------------------------------------
   allocate data structs that have to wait until sites exist
   so that nlocal,nghost,maxneigh are set
------------------------------------------------------------------------- */

void AppDiffusionMultiphaseGCN::allocate_data()
{
  // for linear:
  //   make esites large enough for 1 sites and their 1,2 neighbors

  if (engstyle == LINEAR) {
    int emax = 1 + maxneigh*2;
    esites = new int[2*emax];
  }

  echeck = new int[nlocal+nghost];
  neigh_check = new int[nlocal +nghost];
  hopsite = new int[maxneigh*maxneigh + maxneigh];

  box_dims = new double[3];
  box_dims[0] =  domain->boxxhi-domain->boxxlo;
  box_dims[1] = domain->boxyhi-domain->boxylo;
  box_dims[2] = domain->boxzhi-domain->boxzlo;
  half_box_dims = new double[3];
  half_box_dims[0] = box_dims[0] * 0.5;
  half_box_dims[1] = box_dims[1] * 0.5;
  half_box_dims[2] = box_dims[2] * 0.5;

  memory->create(firstevent,nlocal,"app:firstevent");
}
