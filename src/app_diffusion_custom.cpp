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
#include "app_diffusion_custom.h"
#include "comm_lattice.h"
#include "solve.h"
#include "domain.h"
#include "random_mars.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <bits/stdc++.h>
#include "json.hpp"
#include <fstream>
using json = nlohmann::json;


using namespace SPPARKS_NS;

enum{ZERO,VACANT,OCCUPIED,TOP};
enum{NO_ENERGY,LINEAR,NONLINEAR};
enum{DEPOSITION,NNHOP,SCHWOEBEL};
enum{NOSWEEP,RANDOM,RASTER,COLOR,COLOR_STRICT};  // from app_lattice.cpp
enum{DEP_NONE,DEP_EVENT,DEP_BATCH};

#define DELTAEVENT 100000
#define DELTABATCH 1024
#define NUHOP 1e13
#define TOLERANCE 1.0e-5

/* ---------------------------------------------------------------------- */

AppDiffusionCustom::AppDiffusionCustom(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  // these can be changed by model choice, see below

  ninteger = 1;
  ndouble = 0;
  delpropensity = 2;
  delevent = 1;
  allow_kmc = 1;
  allow_rejection = 0;
  allow_masking = 0;
  numrandom = 1;

  // parse arguments

  if (narg < 3) error->all(FLERR,"Illegal app_style command");
  // if (strcmp(arg[1],"off") == 0) engstyle = NO_ENERGY;
  // else if (strcmp(arg[1],"linear") == 0) engstyle = LINEAR;
  // else if (strcmp(arg[1],"nonlinear") == 0) engstyle = NONLINEAR;
  // else error->all(FLERR,"Illegal app_style command");
  nn = std::stoi(arg[1]);
  engstyle = LINEAR; // Only linear energy is supported for now
  if (strcmp(arg[2],"hop") == 0) {
    // if (narg != 3) error->all(FLERR,"Illegal app_style command");
    hopstyle = NNHOP;
  } else error->all(FLERR,"Illegal app_style command: only nearest-neighbor hop is supported");

  // increment delpropensity by 1 for nonlinear energy
  // increment delpropensity and delevent by 1 for Schwoebel hops
  // change allow_rejection to 1 for linear energy and non-Schwoebel hops

  // if (engstyle == NONLINEAR) delpropensity++;
  // if (hopstyle == SCHWOEBEL) delpropensity++;
  // if (hopstyle == SCHWOEBEL) delevent++;
  // if (engstyle == LINEAR && hopstyle == NNHOP) allow_rejection = 1;
  allow_rejection = 1;
  create_arrays();

  std::ifstream in2(arg[3]);
  json json_file_1;
  in2 >> json_file_1;
  in2.close();
  
  w1_1 = json_file_1["w1"][0];
  w1_2 = json_file_1["w1"][1];
  w2_1 = json_file_1["w2"][0];
  w2_2 = json_file_1["w2"][1];
  w3_1 = json_file_1["w3"][0];
  w3_2 = json_file_1["w3"][1];
  w4_1 = json_file_1["w4"][0];
  w4_2 = json_file_1["w4"][1];
  w5_1 = json_file_1["w5"][0];
  w5_2 = json_file_1["w5"][1];
  
  
  w6 = json_file_1["w6"][0];
  
  
  b1 = json_file_1["b1"];
  b2 = json_file_1["b2"];
  b3 = json_file_1["b3"];
  b4 = json_file_1["b4"];
  b5 = json_file_1["b5"];
  b6 = json_file_1["b6"];
  
  std::ifstream in(arg[4]);
  json j;
  in >> j;

  first = j["1nn"];
  second = j["2nn"];
  third = j["3nn"];
  fourth = j["4nn"];
  fifth = j["5nn"];
  destinations = j["fcc_destinations"];
  in.close();

  // csv_file.open(arg[5], std::ios::app); // Open file in append mode

  esites = psites = NULL;
  echeck = pcheck  =  NULL;
  maxevent = 0;
  events = NULL;
  firstevent = NULL;
  box_dims = NULL;
  half_box_dims = NULL;

  hbarrier = sbarrier = NULL;
  ecoord = NULL;
  hopsite = NULL;
  neigh_check = NULL;

  marklist = NULL;
  mark = NULL;

  allocated = 0;

  // default settings for app-specific commands

  depmode = DEP_NONE;
  barrierflag = 0;

  // batch deposition data

  ranbatch = NULL;

  maxbatch = 0;
  startpos = NULL;
  depinfo = NULL;
  depinfo_copy = NULL;
  elist = NULL;

  // statistics

  ndeposit = ndeposit_failed = 0;
  nfirst = nsecond = 0;
}

/* ---------------------------------------------------------------------- */

AppDiffusionCustom::~AppDiffusionCustom()
{
  delete [] esites;
  delete [] psites;
  delete [] echeck;
  delete [] pcheck;
  delete [] box_dims;
  delete [] half_box_dims;
  // Destroy vectors
  first.clear();
  second.clear();
  first.shrink_to_fit();
  second.shrink_to_fit();
  
  // if (csv_file.is_open()) {
  //   csv_file.close(); // Close the file
  // }

  memory->sfree(events);
  memory->destroy(firstevent);

  memory->destroy(hbarrier);
  memory->destroy(sbarrier);
  delete [] ecoord;

  delete [] hopsite;
  delete [] neigh_check;

  delete [] marklist;
  memory->destroy(mark);

  delete ranbatch;
  memory->destroy(startpos);
  memory->sfree(depinfo);
  memory->sfree(depinfo_copy);
  memory->destroy(elist);
}

/* ----------------------------------------------------------------------
   input script commands unique to this app
------------------------------------------------------------------------- */

void AppDiffusionCustom::input_app(char *command, int narg, char **arg)
{
  if (sites_exist == 0) {
    char str[128];
    sprintf(str,"Cannot use %s command until sites exist",command);
    error->all(FLERR,str);
  }

  if (!allocated) allocate_data();
  allocated = 1;

  if (strcmp(command,"ecoord") == 0) {
    if (engstyle != NONLINEAR)
      error->all(FLERR,"Can only use ecoord command with "
		 "app_style diffusion nonlinear");
    if (narg != 2) error->all(FLERR,"Illegal ecoord command");

    int lo,hi;
    bounds(arg[0],0,maxneigh,lo,hi);
    double value = atof(arg[1]);

    for (int i = lo; i <= hi; i++) ecoord[i] = value;

  } else if (strcmp(command,"deposition") == 0) {
    if (narg < 1) error->all(FLERR,"Illegal deposition command");

    if (strcmp(arg[0],"off") == 0) {
      if (narg != 1) error->all(FLERR,"Illegal deposition command");
      depmode = DEP_NONE;
      return;
    }

    if (narg != 8) error->all(FLERR,"Illegal deposition command");

    if (strcmp(arg[0],"event") == 0) depmode = DEP_EVENT;
    else if (strcmp(arg[0],"batch") == 0) depmode = DEP_BATCH;

    deprate_total = atof(arg[1]);
    dir[0] = atof(arg[2]);
    dir[1] = atof(arg[3]);
    dir[2] = atof(arg[4]);
    d0 = atof(arg[5]);
    coordlo = atoi(arg[6]);
    coordhi = atoi(arg[7]);

    if (deprate_total < 0.0) error->all(FLERR,"Illegal deposition command");
    if (domain->dimension == 2 && (dir[1] >= 0.0 || dir[2] != 0.0))
      error->all(FLERR,"Illegal deposition command");
    if (domain->dimension == 3 && dir[2] >= 0.0)
      error->all(FLERR,"Illegal deposition command");
    if (d0 < 0.0) error->all(FLERR,"Illegal deposition command");
    if (coordlo < 0 || coordhi > maxneigh || coordlo > coordhi)
      error->all(FLERR,"Illegal deposition command");

    double len = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
    dir[0] /= len;
    dir[1] /= len;
    dir[2] /= len;

    // ranbatch = RNG for batch deposition (same on all procs)

    if (depmode == DEP_BATCH) {
      delete ranbatch;
      ranbatch = new RandomPark(ranmaster->uniform());
    }

  } else if (strcmp(command,"barrier") == 0) {
    if (narg < 1) error->all(FLERR,"Illegal barrier command");
    barrierflag = 1;

    double **barrier;
    if (strcmp(arg[0],"none") == 0) {
      if (narg != 1) error->all(FLERR,"Illegal barrier command");
      barrierflag = 0;
      return;
    } else if (strcmp(arg[0],"hop") == 0) {
      barrier = hbarrier;
    } else if (strcmp(arg[0],"schwoebel") == 0) {
      barrier = sbarrier;
    } else error->all(FLERR,"Illegal barrier command");
    if (barrier == sbarrier && hopstyle != SCHWOEBEL)
      error->all(FLERR,
                 "Cannot define Schwoebel barrier without Schwoebel model");

    if (narg < 2 || narg > 4) error->all(FLERR,"Illegal barrier command");
    if (narg == 2) {
      double q = atof(arg[1]);
      int i,j;
      for (i = 0; i <= maxneigh; i++)
	for (j = 0; j <= maxneigh; j++)
	  barrier[i][j] = q;

    } else if (narg == 3) {
      int delta = atoi(arg[1]);
      double q = atof(arg[2]);
      int i,j;
      for (i = 0; i <= maxneigh; i++)
	for (j = 0; j <= maxneigh; j++)
	  if (j-i == delta) barrier[i][j] = q;

    } else {
      int ilo,ihi,jlo,jhi;
      bounds(arg[1],0,maxneigh,ilo,ihi);
      bounds(arg[2],0,maxneigh,jlo,jhi);
      double q = atof(arg[3]);

      for (int i = ilo; i <= ihi; i++)
	for (int j = jlo; j <= jhi; j++)
	  barrier[i][j] = q;
    }

  } else error->all(FLERR,"Unrecognized command");
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppDiffusionCustom::grow_app()
{
  lattice = iarray[0];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppDiffusionCustom::init_app()
{
  if (engstyle != LINEAR && !solve)
    error->all(FLERR,"App diffusion off or nonlinear requires KMC solver");
  if (depmode != DEP_NONE && !solve)
    error->all(FLERR,"App diffusion with deposition requires KMC solver");

  if (depmode == DEP_EVENT && nprocs > 1)
    error->all(FLERR,"Cannot use deposition event in parallel - use batch");
  if (depmode == DEP_EVENT && sectorflag)
    error->all(FLERR,"Cannot use deposition event with multiple sectors - "
               "use batch");
  if (depmode == DEP_BATCH && !sectorflag)
    error->all(FLERR,"Cannot use deposition batch without sectors - use event");

  // set per-processor deposition rate
  // deposition in batch mode uses callback at end of sweep

  if (depmode != DEP_NONE) {
    if (depmode == DEP_EVENT) deprate = deprate_total;
    else if (depmode == DEP_BATCH) {
      int flag = 1;
      if (nlocal == 0) flag = 0;
      int active;
      MPI_Allreduce(&flag,&active,1,MPI_INT,MPI_SUM,world);
      deprate = deprate_total / active;

      if (!elist) memory->create(elist,nlocal,"app:elist");
    }
  }

  nbatch = 0;
  if (depmode == DEP_BATCH) allow_app_update = 1;
  else allow_app_update = 0;

  // allocate data

  if (!allocated) allocate_data();
  allocated = 1;

  dimension = domain->dimension;
  dt_sweep = 1.0/maxneigh;

  // site validity

  int flag = 0;
  for (int i = 0; i < nlocal; i++)
    if (lattice[i] < VACANT || lattice[i] > TOP) flag = 1;
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
}

/* ----------------------------------------------------------------------
   setup before each run
------------------------------------------------------------------------- */

void AppDiffusionCustom::setup_app()
{
  for (int i = 0; i < nlocal+nghost; i++) echeck[i] = pcheck[i] = 0;

  // clear event list

  nevents = 0;
  for (int i = 0; i < nlocal; i++) firstevent[i] = -1;
  for (int i = 0; i < maxevent; i++) events[i].next = i+1;
  freeevent = 0;
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */


double AppDiffusionCustom::site_energy(int i)
{
  // energy only non-zero for OCCUPIED sites when energy included in model

  if (lattice[i] != OCCUPIED || engstyle == NO_ENERGY) return 0.0;

  // energy is a non-linear function of coordination number
  // calculate from user-specified tabulated values

  if (engstyle == NONLINEAR) {
    int n = 0;
    for (int j = 0; j < numneigh[i]; j++)
      if (lattice[neighbor[i][j]] == OCCUPIED) n++;
    return ecoord[n];
  }

  // energy is a linear function of coordination number, just count bonds

  int eng = 0;
  for (int j = 0; j < numneigh[i]; j++)
    if (lattice[neighbor[i][j]] == VACANT) eng++;
  return (double) eng;
}

/* ----------------------------------------------------------------------
   rKMC method
   perform a site event with null bin rejection
   null bin extends to size maxneigh
------------------------------------------------------------------------- */
void AppDiffusionCustom::site_event_rejection(int i, RandomPark *random)
{
  double einitial,edelta;

  // OCCUPIED site exchanges with random neighbor if VACANT

  if (lattice[i] != OCCUPIED) return;
  int iran = (int) (maxneigh*random->uniform());
  if (iran > maxneigh) iran = maxneigh-1;
  int j = neighbor[i][iran];
  if (lattice[j] != VACANT) return;

  // accept or reject via energy and barrier model
  // factor of 2 in edelta accounts for energy change of neighbors of I,J

  int hop = 0;
  if (engstyle != NO_ENERGY) einitial = site_energy(i);

  lattice[i] = VACANT;
  lattice[j] = OCCUPIED;

  if (engstyle == NO_ENERGY) {
    if (!barrierflag) hop = 1;
    else if (temperature > 0.0) {
      if (random->uniform() < exp(-hbarrier[ncoord(i)-1][ncoord(j)]*t_inverse))
	  hop = 1;
    }

  } else {
    edelta = site_energy(j) - einitial;

    if (!barrierflag) {
      if (edelta <= 0.0) hop = 1;
      else if (temperature > 0.0) {
	if (random->uniform() < exp(-2.0*edelta*t_inverse)) hop = 1;
      }
    } else if (temperature > 0.0) {
      if (edelta <= 0.0) {
	if (random->uniform() < 
	    exp(-hbarrier[ncoord(i)-1][ncoord(j)]*t_inverse)) hop = 1;
      } else {
	if (random->uniform() < 
	    exp((-2.0*edelta-hbarrier[ncoord(i)-1][ncoord(j)]) * t_inverse))
	  hop = 1;
      }
    }
  }

  if (hop) {
    naccept++;
    nfirst++;
  } else {
    lattice[i] = OCCUPIED;
    lattice[j] = VACANT;
  }
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppDiffusionCustom::site_propensity(int i)
{
  return site_propensity_linear(i);
}

/* ---------------------------------------------------------------------- */

double AppDiffusionCustom::site_propensity_linear(int i)
{
  int j,ihop,nhop1,nhop2,eflag;
  double probone,proball, eb;
  int l,ll,lll,llll,o;
  int m;
  double dist;
  int  mm,mmm,mmmm;
  int num_occupied_sites = 0;
  std::vector<std::vector<double>> occupied_coords(128,std::vector<double>(3,0.0));
  std::vector<double> xi = {xyz[i][0], xyz[i][1], xyz[i][2]};

  clear_events(i);

  if (lattice[i] != OCCUPIED) {
    if (i == 0 && depmode != DEP_NONE) {
      add_event(i,-1,deprate,DEPOSITION);
      return deprate;
    } else return 0.0;
  }
 
  nhop1 = 0;
  neigh_check[i] = 1; // Marking the site i as checked
  // Check for neighbors of site i that are occupied
  // and store their coordinates in occupied_coords
  // And identify possible vacancy sites and store them in hopsite
  for (l=0; l<numneigh[i];l++){
    m = neighbor[i][l];
    if (lattice[m]==OCCUPIED && neigh_check[m] == 0){
      for(o = 0; o<3; o++){
        occupied_coords[num_occupied_sites][o] = xyz[m][o] - xi[0];
      }
      
      num_occupied_sites++;
    }
    else if (lattice[m]==VACANT){
      hopsite[nhop1++] = m;
    }
    neigh_check[m] = 1; // Marking the neighbor as checked
    for (ll=0; ll<maxneigh; ll++){
        mm = neighbor[m][ll];

        if (lattice[mm] == OCCUPIED && neigh_check[mm] == 0){
          for(o = 0; o<3; o++){
            occupied_coords[num_occupied_sites][o] = xyz[mm][o] - xi[0];
          }
          num_occupied_sites++;
        }
        neigh_check[mm] = 1;  // Setting the neighbor check for the second layer
        for (lll=0; lll<maxneigh; lll++){
          mmm = neighbor[mm][lll];
          if (lattice[mmm]==OCCUPIED && neigh_check[mmm]==0){
            for(o = 0; o<3; o++){
              occupied_coords[num_occupied_sites][o] = xyz[mmm][o] - xi[0];
            }
            num_occupied_sites++;
          }
          neigh_check[mmm] = 1;
          for (llll=0; llll<maxneigh; llll++){
            mmmm = neighbor[mmm][llll];
            if (lattice[mmmm]==OCCUPIED && neigh_check[mmmm]==0){
              for(o = 0; o<3; o++){
                occupied_coords[num_occupied_sites][o] = xyz[mmmm][o] - xi[0];
              }
              num_occupied_sites++;
            }
            neigh_check[mmmm] = 1;
          }
        }
    }
  }
  // consider the periodicity, if the distance is greater than half the box dimension,
  // adjust the distance to account for periodic boundary conditions
  for(l=0; l<num_occupied_sites; l++){
    for(o = 0; o<3; o++){
      dist = occupied_coords[l][o];
      if (dist>half_box_dims[o]) dist-=box_dims[o];
      else if (dist<-half_box_dims[o]) dist+=box_dims[o];
      occupied_coords[l][o] = dist;
    }
  }
  // Calculate the probabilities of hopping from site i to j
  // and add them to the event list
  proball = 0.0;
  for(ihop= 0; ihop<nhop1; ihop++){
    j =  hopsite[ihop];
    eflag = NNHOP;
    eb = calculate_barrier_energy(i,j,occupied_coords,num_occupied_sites);
    probone = (NUHOP)*exp(-eb*t_inverse);
    add_event(i,j,probone,eflag);
    proball += probone;
 }
  if (i == 0 && depmode != DEP_NONE) {
    add_event(i,-1,deprate,DEPOSITION);
    proball += deprate;
  }
  // Reset the neigh_check array for the next site
  for ( l = 0; l < nlocal + nghost; l++) neigh_check[l] = 0;
  return proball;
}

/* ---------------------------------------------------------------------- */
//    Helper functions to perform vector and matrix operations
/* ---------------------------------------------------------------------- */


double AppDiffusionCustom::vec_dot(const std::vector<double> &a,const std::vector<double> &b){
  double c = a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
  return c;
}

std::vector<double> AppDiffusionCustom::vec_cross(const std::vector<double> &a,const std::vector<double> &b){
  std::vector<double> c(3,0.0);
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
  return c;
}

double AppDiffusionCustom::vec_norm(const std::vector<double> &a){
  double c = std::sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
  return c;
}

std::vector<double> AppDiffusionCustom::mat_vecmul(const std::vector<std::vector<double>> &B,const std::vector<double> &a){
  
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

std::vector<double> AppDiffusionCustom::relu(std::vector<double> &z){
  std::vector <double> a(z.size());
  int l;
  for (l = 0; l<z.size();l++){
    if (z[l] < 0){
      a[l] = 0;
    }
    else{
      a[l] = z[l];
    }
  }
  return a;
}

std::vector<double> AppDiffusionCustom::add_vec(const std::vector<double> &a, const std::vector<double> &b){
  if (a.size() != b.size()){
      throw std::invalid_argument("Vectors must be of the same size");
  }
  std::vector<double> c(a.size());
  for (int k = 0; k<a.size();k++){
      c[k] = a[k] + b[k];
  }
  return c;
}

std::vector<double> AppDiffusionCustom::vec_matmul(const std::vector<double> &a,const std::vector<std::vector<double>> &B){
  
  int m = a.size();
  int p = B.size();
  int q = B[0].size();
  int k, l;
  
  if (m !=p){
      throw std::invalid_argument("Matrix dimensions do not match for multiplication");
  }
  std::vector<double> c(q,0.0);
  for (int l =0; l<q;l++){
      for(int k=0; k<m; k++){
          c[l] += a[k] * B[k][l];
      }
  }
  return c;
}

/* ---------------------------------------------------------------------- */
//    Main function to calculate the barrier energy using the NN-model provided by the user
/* ---------------------------------------------------------------------- */
double AppDiffusionCustom::calculate_barrier_energy(int i, int j, std::vector<std::vector <double>> &local_coords, int num_occupied_sites){
  // Calculate the barrier energy for the hop from site i to site j
  // using the local coordinates of the occupied sites
  int k, kk, o, flag,l;
  double x , y,z, dist;
  double n_occupied=0.0;
  
  std::vector <double> first_shell(18,0.0);
  std::vector <double> second_shell(8,0.0);
  std::vector <double> third_shell(32,0.0);
  std::vector <double> fourth_shell(14,0.0);
  std::vector <double> fifth_shell(28,0.0);

  std::vector <double> new_x(3);
  std::vector <double> new_y(3);
  std::vector <double> new_z(3);

  std::vector<std::vector<double>> R(3,std::vector<double>(3,0));
  std::vector<std::vector<double>> original_cs = {{1,0,0}, {0,1,0}, {0,0,1}};
  std::vector<std::vector<double>> rotated_cs(3,std::vector<double>(3,0));

  std::vector<std::vector<double>> rotated_coords(num_occupied_sites, std::vector<double>(3,0));

  for ( o=0; o<3; o++){
    new_x[o] = xyz[j][o] - xyz[i][o];
    if (new_x[o]>4){
      new_x[o] = new_x[o]- box_dims[o] ;
    }
    else if (new_x[o]<-4){
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
  }

  //Checking the occupancy of the local coordinates

  for (k =0; k<num_occupied_sites;k++){
    for (kk=0;kk<first_shell.size();kk++){
      flag=1;
      for(o=0; o<3;o++){
        if (fabs(first[kk][o]-rotated_coords[k][o])>TOLERANCE){
          flag=0;
          break;
        }
      }
      if(flag==1){
        first_shell[kk] = 1.0/2;
        n_occupied++;
        break;
      }                                                                                      
    }
    // Moving to the second shell
    for (kk=0;kk<second_shell.size();kk++){
      flag=1;
      for(o=0; o<3;o++){
        if (fabs(second[kk][o]-rotated_coords[k][o])>TOLERANCE){
          flag=0;
          break;
        }
      }
      if(flag==1){
        second_shell[kk] = 1.0/4;
        n_occupied++;
        break;
      }
    }
    // Moving to the third shell
    for (kk=0;kk<third_shell.size();kk++){
      flag=1;
      for(o=0; o<3;o++){
        if (fabs(third[kk][o]-rotated_coords[k][o])>TOLERANCE){
          flag=0;
          break;
        }
      }
      if(flag==1){
        third_shell[kk] = 1.0/6;
        n_occupied++;
        break;
      }
    }
    // Moving to the fourth shell 
    for (kk=0;kk<fourth_shell.size();kk++){
      flag=1;
      for(o=0; o<3;o++){
        if (fabs(fourth[kk][o]-rotated_coords[k][o])>TOLERANCE){
          flag=0;
          break;
        }
      }
      if(flag==1){
        fourth_shell[kk] = 1.0/8;
        n_occupied++;
        break;
      }
    }
    // Moving to the fifth shell
    for (kk=0;kk<fifth_shell.size();kk++){
      flag=1;
      for(o=0; o<3;o++){
        if (fabs(fifth[kk][o]-rotated_coords[k][o])>TOLERANCE){
          flag=0;
          break;
        }
      }
      if(flag==1){
        fifth_shell[kk] = 1.0/10;
        n_occupied++;
        break;
      }
    }
  }

  // First layer matrix operations
    std::vector<double> z1_1 = add_vec( vec_matmul(first_shell, w1_1), b1[0]);
    std::vector<double> z2_1 = add_vec(vec_matmul(second_shell, w2_1), b2[0]);
    std::vector<double> z3_1 = add_vec( vec_matmul(third_shell, w3_1), b3[0]);
    std::vector<double> z4_1 = add_vec(vec_matmul(fourth_shell, w4_1), b4[0]);
    std::vector<double> z5_1 = add_vec(vec_matmul(fifth_shell, w5_1), b5[0]);


    std::vector<double> a1_1 = relu(z1_1);
    std::vector<double> a2_1 = relu(z2_1);
    std::vector<double> a3_1 = relu(z3_1);
    std::vector<double> a4_1 = relu(z4_1);
    std::vector<double> a5_1 = relu(z5_1);



  // Second layer matrix operations

    std::vector<double> z1_2 = add_vec(vec_matmul(a1_1, w1_2), b1[1]);
    std::vector<double> z2_2 = add_vec(vec_matmul(a2_1, w2_2), b2[1]);
    std::vector<double> z3_2 = add_vec(vec_matmul(a3_1, w3_2), b3[1]);
    std::vector<double> z4_2 = add_vec(vec_matmul(a4_1, w4_2), b4[1]);
    std::vector<double> z5_2 = add_vec(vec_matmul(a5_1, w5_2), b5[1]);

    std::vector<double> a1_2 = relu(z1_2);
    std::vector<double> a2_2 = relu(z2_2);
    std::vector<double> a3_2 = relu(z3_2);
    std::vector<double> a4_2 = relu(z4_2);
    std::vector<double> a5_2 = relu(z5_2);

  // Concatenating and performing final layer matrix operations
    a1_2.insert(a1_2.end(), a2_2.begin(), a2_2.end());
    a1_2.insert(a1_2.end(), a3_2.begin(), a3_2.end());
    a1_2.insert(a1_2.end(), a4_2.begin(), a4_2.end());
    a1_2.insert(a1_2.end(), a5_2.begin(), a5_2.end());
    std::vector<double> z3 = add_vec(vec_matmul(a1_2, w6), b6[0]);


  // Code for diagnosis purpose
  // double sum= 0.0;
  // for(l=0; l<11; l++){
  //   sum += first_shell[l]*2;
  // }
  
  // Write data to the CSV file for diagnosis
  // if (csv_file.is_open()) {
  //   for (l = 0; l < first_shell.size(); l++) {
  //     csv_file << first_shell[l] << ",";
  //   }
  //   for (l = 0; l < second_shell.size(); l++) {
  //     csv_file << second_shell[l] << ",";
  //   }
  //   for (l = 0; l < third_shell.size(); l++) {
  //     csv_file << third_shell[l] << ",";
  //   }
  //   for (l = 0; l < fourth_shell.size(); l++) {
  //     csv_file << fourth_shell[l] << ",";
  //   }
  //   for (l = 0; l < fifth_shell.size(); l++) {
  //     csv_file << fifth_shell[l] << ",";
  //   }

  //   csv_file << z3[0] << "\n"; 
  // }

  return z3[0];
}



/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppDiffusionCustom::site_event(int i, class RandomPark *random)
{
  return site_event_linear(i,random);
}

/* ---------------------------------------------------------------------- */

void AppDiffusionCustom::site_event_linear(int i, class RandomPark *random)
{
  int j,m;

  // pick one event from total propensity by accumulating its probability
  // compare prob to threshhold, break when reach it to select event
  // perform event

  double threshhold = random->uniform() * propensity[i2site[i]];
  double proball = 0.0;

  int ievent = firstevent[i];
  while (1) {
    proball += events[ievent].propensity;
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
    if (ievent < 0) error->one(FLERR,"Did not reach event propensity threshhold");
  }

  // deposition or hop event
  // for deposition event:
  //   if batch mode, tally for end-of-sweep deposition
  //   if event mode, find site to deposit on
  //   after deposition, reset i and j to that site,
  //   so propensity around it is updated correctly
  // for hop event: perform hop

  if (events[ievent].style == DEPOSITION) {
    if (depmode == DEP_BATCH) {
      nbatch++;
      return;
    }
    m = find_deposition_site(random);
    if (m < 0) return;
    lattice[m] = OCCUPIED;
    i = j = m;
  } else {
    j = events[ievent].destination;
    if (events[ievent].style == NNHOP) nfirst++;
    else nsecond++;
    lattice[i] = VACANT;
    lattice[j] = OCCUPIED;
  }

  // update propensities of all affected sites
  
  update_propensities(i,j);
}

/* ---------------------------------------------------------------------- */


/* ----------------------------------------------------------------------
   compute propensity changes for self and swap site and their neighs
   if i = j, then is a deposition event
   engstyle = NO_ENERGY or LINEAR:
     1,2 neighs for NNHOP and 1,2,3 neighs for SCHWOEBEL
   engstyle = NONLINEAR:
     1,2,3 neighs for NNHOP and 1,2,3,4 neighs for SCHWOEBEL
   ignore update of sites with isite < 0
   use echeck[] to avoid resetting propensity of same site
------------------------------------------------------------------------- */

void AppDiffusionCustom::update_propensities(int i, int j)
{
  int nsites = 0;

  int isite = i2site[i];
  propensity[isite] = site_propensity(i);
  esites[nsites++] = isite;
  echeck[isite] = 1;

  if (j != i) {
    isite = i2site[j];
    if (isite >= 0) {
      propensity[isite] = site_propensity(j);
      esites[nsites++] = isite;
      echeck[isite] = 1;
    }
  }

  if (engstyle != NONLINEAR) {
    if (hopstyle == NNHOP) {
      nsites += neighbor2(i,&esites[nsites]);
      if (j != i) nsites += neighbor2(j,&esites[nsites]);
    } else {
      nsites += neighbor3(i,&esites[nsites]);
      if (j != i) nsites += neighbor3(j,&esites[nsites]);
    }
  } else {
    if (hopstyle == NNHOP) {
      nsites += neighbor3(i,&esites[nsites]);
      if (j != i) nsites += neighbor3(j,&esites[nsites]);
    } else {
      nsites += neighbor4(i,&esites[nsites]);
      if (j != i) nsites += neighbor4(j,&esites[nsites]);
    }
  }

  // reset propensities of all affected sites within solver

  solve->update(nsites,esites,propensity);

  // clear echeck array

  for (int m = 0; m < nsites; m++) echeck[esites[m]] = 0;

  // sanity check on all propensity values

  /*
  printf("EVENT %d %d\n",i,j);
  for (m = 0; m < nlocal; m++) {
    if (fabs(propensity[m]-site_propensity(m)) > 1.0e-6) {
      printf("BAD PROP = %d %d %d %g %g\n",
	     id[i],id[j],id[m],propensity[m],site_propensity(m));
      error->one(FLERR,"BAD DONE");
    }
  }
  */
}

/* ----------------------------------------------------------------------
   perform deposition events as one batch per sweep
   only done when KMC and sectors are being used
   not done when rKMC b/c cannot do deposition
   not done when no sectors, b/c then app_update() is not called
------------------------------------------------------------------------- */

void AppDiffusionCustom::app_update(double stoptime)
{
  int i,j,m;

  // ntotal = total # of atoms to deposit

  int nbatch_total;
  MPI_Allreduce(&nbatch,&nbatch_total,1,MPI_INT,MPI_SUM,world);

  // data structs for all depositions, stored by all procs
  // initially has local info, at end will have global info

  if (nbatch_total > maxbatch) {
    memory->destroy(startpos);
    memory->sfree(depinfo);
    memory->sfree(depinfo_copy);
    while (nbatch_total > maxbatch) maxbatch += DELTABATCH;
    memory->create(startpos,maxbatch,3,"diffusion:startpos");
    depinfo = (DepInfo *) 
      memory->smalloc(maxbatch*sizeof(DepInfo),"diffusion:depinfo");
    depinfo_copy = (DepInfo *) 
      memory->smalloc(maxbatch*sizeof(DepInfo),"diffusion:depinfo_copy");
  }

  // pick random positions at top of box for all depositions

  for (i = 0; i < nbatch_total; i++) {
    startpos[i][0] = domain->boxxlo + domain->xprd*ranbatch->uniform();
    if (dimension == 2) {
      startpos[i][1] = domain->boxyhi;
      startpos[i][2] = 0.0;
    } else {
      startpos[i][1] = domain->boxylo + domain->yprd*ranbatch->uniform();
      startpos[i][2] = domain->boxzhi;
    }
  }

  // neligible = # of my sites eligible for a deposition atom
  // must be vacant site
  // must have neighbor count between coordlo and coordhi
  // elist = list of eligible site indices
  
  int ncount;
  int neligible = 0;
 
  for (i = 0; i < nlocal; i++) {
    if (lattice[i] != VACANT) continue;
    ncount = 0;
    for (j = 0; j < numneigh[i]; j++)
      if (lattice[neighbor[i][j]] == OCCUPIED) ncount++;
    if (ncount < coordlo || ncount > coordhi) continue;

    elist[neligible] = i;
    neligible++;
  }

  // double loop over depositions and my eligible sites
  // check 3rd condition: 
  //   must be close enough to deposition path, tested by exceeds_limit()
  // depinfo[I] = my site with projected distance closest to start point
  //   for Ith deposition in batch

  int closesite;
  double closedist,dist2start;

  for (i = 0; i < nbatch_total; i++) {
    closesite = -1;
    closedist = 1.0e20;
    for (j = 0; j < neligible; j++) {
      m = elist[j];
      if (exceeds_limit(m,startpos[i],dist2start)) continue;
      if (dist2start < closedist) {
        closedist = dist2start;
        closesite = m;
      }
    }
    depinfo[i].proc = me;
    depinfo[i].site = closesite;
    depinfo[i].distance = closedist;
  }

  // merge DepInfo data struct across all procs
  // logarithmic pairwise reduction
  // final depinfo = list of which procs own which sites to deposit onto
  // nsplit = proc ID of first proc in upper half at each iteration
  // procs in upper half send to lower half, allowing for non-power-of-two procs

  MPI_Request request;
  MPI_Status status;

  int nsplit = 1;
  while (nsplit < nprocs) nsplit *= 2;
  nsplit /= 2;

  int partner;

  while (nsplit >= 1) {
    if (me >= nsplit && me < 2*nsplit) partner = me - nsplit;
    else if (me < nsplit && me + nsplit < nprocs) partner = me + nsplit;
    else partner = -1;

    if (partner >= 0 && me < partner) {
      MPI_Irecv(depinfo_copy,nbatch_total*sizeof(DepInfo),MPI_CHAR,
                partner,0,world,&request);
      MPI_Wait(&request,&status);

      for (i = 0; i < nbatch_total; i++) {
        if (depinfo_copy[i].site < 0) continue;
        if (depinfo_copy[i].distance < depinfo[i].distance) {
          depinfo[i].proc = depinfo_copy[i].proc;
          depinfo[i].site = depinfo_copy[i].site;
          depinfo[i].distance = depinfo_copy[i].distance;
        }
      }
    }

    if (partner >= 0 && me > partner)
      MPI_Send(depinfo,nbatch_total*sizeof(DepInfo),MPI_CHAR,partner,0,world);

    nsplit /= 2;
  }

  // broadcast final merged DepInfo data struct to all procs

  MPI_Bcast(depinfo,nbatch_total*sizeof(DepInfo),MPI_CHAR,0,world);

  // loop over list of all batched depositions
  // tally ndeposit and ndeposit_failed
  // if I own the site, deposit the atom on my site
  // elist = list of my sites with new atoms

  int n = 0;
  for (i = 0; i < nbatch_total; i++) {
    if (depinfo[i].site < 0) {
      ndeposit_failed++;
      continue;
    }
    ndeposit++;

    if (depinfo[i].proc != me) continue;
    m = depinfo[i].site;
    lattice[m] = OCCUPIED;
    elist[n] = m;
    n++;
  }

  // comm to acquire all my ghost site values
  // since deposition may have updated them

  comm->all();

  // reset propensities for all sites affected by deposition events
  // this must be done within sectors
  //   b/c KMC solver for that sector must be updated
  // 2 kinds of affected sites:
  //   (a) all sites affected by a new atom in that sector
  //   (b) sites in owned border region of each sector,
  //       may be affected by a new atom in another sector (or other proc)

  int iset;

  for (i = 0; i < nbatch_total; i++) {
    if (depinfo[i].site < 0) continue;
    if (depinfo[i].proc != me) continue;

    // iset = the set which site M is in
    // set solve,propensity,i2site specific to that set
    // used by update_propensities()

    m = depinfo[i].site;
    iset = whichset(m);

    solve = set[iset].solve;
    propensity = set[iset].propensity;
    i2site = set[iset].i2site;

    update_propensities(m,m);
  }

  update_kmc_sector_border_propensities();

  // reset nbatch for next AppLattice loop

  nbatch = 0;
}

// ----------------------------------------------------------------------
// private functions
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   re-compute propensities out to 2nd neighbors of site I
------------------------------------------------------------------------- */

int AppDiffusionCustom::neighbor2(int i, int *sites)
{

  int k,kk,m,mm,isite;

  int nsites = 0;

  for (k = 0; k < numneigh[i]; k++) {
    m = neighbor[i][k];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
      propensity[isite] = site_propensity(m);
      sites[nsites++] = isite;
      echeck[isite] = 1;
    }
    for (kk = 0; kk < numneigh[m]; kk++) {
      mm = neighbor[m][kk];
      isite = i2site[mm];
      if (isite >= 0 && echeck[isite] == 0) {
	propensity[isite] = site_propensity(mm);
	sites[nsites++] = isite;
	echeck[isite] = 1;
      }
    }
  }

  return nsites;
}

/* ----------------------------------------------------------------------
   re-compute propensities out to 3rd neighbors of site I
------------------------------------------------------------------------- */

int AppDiffusionCustom::neighbor3(int i, int *sites)
{
  int k,kk,kkk,m,mm,mmm,isite;

  int nsites = 0;

  for (k = 0; k < numneigh[i]; k++) {
    m = neighbor[i][k];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
      propensity[isite] = site_propensity(m);
      sites[nsites++] = isite;
      echeck[isite] = 1;
    }
    for (kk = 0; kk < numneigh[m]; kk++) {
      mm = neighbor[m][kk];
      isite = i2site[mm];
      if (isite >= 0 && echeck[isite] == 0) {
	propensity[isite] = site_propensity(mm);
	sites[nsites++] = isite;
	echeck[isite] = 1;
      }
      for (kkk = 0; kkk < numneigh[mm]; kkk++) {
	mmm = neighbor[mm][kkk];
	isite = i2site[mmm];
	if (isite >= 0 && echeck[isite] == 0) {
	  propensity[isite] = site_propensity(mmm);
	  sites[nsites++] = isite;
	  echeck[isite] = 1;
	}
      }
    }
  }

  return nsites;
}

/* ----------------------------------------------------------------------
   re-compute propensities out to 4th neighbors of site I
------------------------------------------------------------------------- */

int AppDiffusionCustom::neighbor4(int i, int *sites)
{
  int k,kk,kkk,kkkk,m,mm,mmm,mmmm,isite;

  int nsites = 0;

  for (k = 0; k < numneigh[i]; k++) {
    m = neighbor[i][k];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
      propensity[isite] = site_propensity(m);
      sites[nsites++] = isite;
      echeck[isite] = 1;
    }
    for (kk = 0; kk < numneigh[m]; kk++) {
      mm = neighbor[m][kk];
      isite = i2site[mm];
      if (isite >= 0 && echeck[isite] == 0) {
	propensity[isite] = site_propensity(mm);
	sites[nsites++] = isite;
	echeck[isite] = 1;
      }
      for (kkk = 0; kkk < numneigh[mm]; kkk++) {
	mmm = neighbor[mm][kkk];
	isite = i2site[mmm];
	if (isite >= 0 && echeck[isite] == 0) {
	  propensity[isite] = site_propensity(mmm);
	  sites[nsites++] = isite;
	  echeck[isite] = 1;
	}
	for (kkkk = 0; kkkk < numneigh[mmm]; kkkk++) {
	  mmmm = neighbor[mmm][kkkk];
	  isite = i2site[mmmm];
	  if (isite >= 0 && echeck[isite] == 0) {
	    propensity[isite] = site_propensity(mmmm);
	    sites[nsites++] = isite;
	    echeck[isite] = 1;
	  }
	}
      }
    }
  }

  return nsites;
}

/* ---------------------------------------------------------------------- */

int AppDiffusionCustom::ncoord(int i)
{
  int count = 0;
  for (int j = 0; j < numneigh[i]; j++)
    if (lattice[neighbor[i][j]] == OCCUPIED) count++;
  return count;
}

/* ----------------------------------------------------------------------
   clear all events out of list for site I
   add cleared events to free list
------------------------------------------------------------------------- */

void AppDiffusionCustom::clear_events(int i)
{
  int next;
  int index = firstevent[i];
  while (index >= 0) {
    next = events[index].next;
    events[index].next = freeevent;
    freeevent = index;
    nevents--;
    index = next;
  }
  firstevent[i] = -1;
}

/* ----------------------------------------------------------------------
   add an event to list for site I
   hop event = exchange with site J with probability = propensity
   deposition event = add an atom
------------------------------------------------------------------------- */

void AppDiffusionCustom::add_event(int i, int destination, 
                             double propensity, int eventflag)
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
  events[freeevent].style = eventflag;
  events[freeevent].next = firstevent[i];
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}

/* ----------------------------------------------------------------------
   enumerate Schwoebel hop events centered around OCCUPIED site I
   assume mark array is currently cleared, use it, clear it when done
------------------------------------------------------------------------- */

int AppDiffusionCustom::schwoebel_enumerate(int i, int *site)
{
  int j,k,m,jneigh,kneigh,count;

  int nhop = 0;

  // if coord(I) > Nmax, no hops possible

  count = 0;
  for (j = 0; j < numneigh[i]; j++)
    if (lattice[neighbor[i][j]] == OCCUPIED) count++;
  if (count > nsmax) return nhop;

  // mark first neighbors of site I as vacant = 1, occupied = 2

  for (jneigh = 0; jneigh < numneigh[i]; jneigh++) {
    j = neighbor[i][jneigh];
    if (lattice[j] == VACANT) mark[j] = 1;
    else if (lattice[j] == OCCUPIED) mark[j] = 2;
  }

  // loop over 1st and 2nd neighbors of site I to find possible hops to K
  // if K not vacant, no hop
  // if mark(K) = 1 or 2, K is a 1st neigh, not a 2nd neigh
  // if mark(K) = 30, already seen as a possible hop
  // mark(K) = 10 or 20, it has a 1st neigh that is vacant or occupied
  // mark(K) = 30, has both a vacant/occupied 1st neigh, so consider it
  // if coord(K) < Nmin, no hop possible
  // if all criteria met, then it is a candidate hop, add to site[]
  
  int nlist = 0;
  for (jneigh = 0; jneigh < numneigh[i]; jneigh++) {
    j = neighbor[i][jneigh];
    for (kneigh = 0; kneigh < numneigh[j]; kneigh++) {
      k = neighbor[j][kneigh];
      if (lattice[k] != VACANT) continue;
      if (mark[k] == 1 || lattice[k] == 2) continue;
      if (mark[k] == 30) continue;
      if (mark[k] == 10*mark[j]) continue;
      if (mark[k] == 0) marklist[nlist++] = k;
      mark[k] += 10*mark[j];
      if (mark[k] != 30) continue;

      count = 0;
      for (m = 0; m < numneigh[k]; m++)
	if (lattice[neighbor[k][m]] == OCCUPIED) count++;
      if (count < nsmin) continue;

      site[nhop++] = k;
    }
  }

  // clear marked sites, 1st and 2nd neighbors

  for (j = 0; j < numneigh[i]; j++) mark[neighbor[i][j]] = 0;
  for (k = 0; k < nlist; k++) mark[marklist[k]] = 0;

  return nhop;
}

/* ----------------------------------------------------------------------
   identify a VACANT site to deposit an atom
   return -1 if could not find a suitable site
------------------------------------------------------------------------- */

int AppDiffusionCustom::find_deposition_site(RandomPark *random)
{
  // pick a random position at top of box

  double start[3];
  start[0] = domain->boxxlo + domain->xprd*random->uniform();
  if (dimension == 2) {
    start[1] = domain->boxyhi;
    start[2] = 0.0;
  } else {
    start[1] = domain->boxylo + domain->yprd*random->uniform();
    start[2] = domain->boxzhi;
  }

  // for each vacant site:
  // discard site if neighbor count not between coordlo and coordhi
  // find site whose projected distance is closest to start point

  int i,ncount;
  double dist2start;

  int closesite = -1;
  double closedist = 1.0e20;

  for (i = 0; i < nlocal; i++) {
    if (lattice[i] != VACANT) continue;
    ncount = 0;
    for (int j = 0; j < numneigh[i]; j++)
      if (lattice[neighbor[i][j]] == OCCUPIED) ncount++;
    if (ncount < coordlo || ncount > coordhi) continue;

    if (exceeds_limit(i,start,dist2start)) continue;
    if (dist2start < closedist) {
      closedist = dist2start;
      closesite = i;
    }
  }

  if (closesite < 0) ndeposit_failed++;
  else ndeposit++;

  return closesite;
}

/* ----------------------------------------------------------------------
   test if site M is within normal distance d0 from incident line
   if so, return 0 and dist2start, else return 1
   site M really becomes a periodic image in XY of M, adjusted via iprd/jprd
   dist2start = dist from site M to starting point of incident line
   dist2start is dist along incident line from start point to
     normal projection point of M
------------------------------------------------------------------------- */

int AppDiffusionCustom::exceeds_limit(int m, double *start, double &dist2start)
{
  int increment,iprd,jprd;

  iprd = jprd = 0;
  double d0sq = d0*d0;

  double distsq = distsq_to_line(m,start,iprd,jprd,dist2start);
  double newdistsq = distsq_to_line(m,start,iprd-1,jprd,dist2start);
  if (newdistsq < distsq) increment = -1;
  else increment = 1;

  iprd += increment;
  newdistsq = distsq_to_line(m,start,iprd,jprd,dist2start);
  while (newdistsq < distsq) {
    distsq = newdistsq;
    iprd += increment;
    newdistsq = distsq_to_line(m,start,iprd,jprd,dist2start);
  }
  iprd -= increment;

  if (dimension == 3) {
    newdistsq = distsq_to_line(m,start,iprd,jprd-1,dist2start);
    if (newdistsq < distsq) increment = -1;
    else increment = 1;

    jprd += increment;
    newdistsq = distsq_to_line(m,start,iprd,jprd,dist2start);
    while (newdistsq < distsq) {
      distsq = newdistsq;
      jprd += increment;
      newdistsq = distsq_to_line(m,start,iprd,jprd,dist2start);
    }
  }
  jprd -= increment;

  if (distsq > d0sq) return 1;
  distsq = distsq_to_line(m,start,iprd,jprd,dist2start);
  return 0;
}

/* ----------------------------------------------------------------------
   compute normal distsq from site M to incident line of deposition
   site M really becomes a periodic image in XY of M, adjusted via iprd/jprd
   also compute and return dist2start
   dist2start = dist from site M to starting point of incident line
   dist2start is dist along incident line from start point to
     normal projection point of M
------------------------------------------------------------------------- */

double AppDiffusionCustom::distsq_to_line(int m, double *start,
				    int iprd, int jprd, double &dist2start)
{
  double delta[3],projection[3],offset[3];

  delta[0] = xyz[m][0] + iprd*domain->xprd - start[0];
  delta[1] = xyz[m][1] + jprd*domain->yprd - start[1];
  delta[2] = xyz[m][2] - start[2];
    
  dist2start = dir[0]*delta[0] + dir[1]*delta[1] + dir[2]*delta[2];
  projection[0] = dist2start*dir[0];
  projection[1] = dist2start*dir[1];
  projection[2] = dist2start*dir[2];
  
  offset[0] = delta[0] - projection[0];
  offset[1] = delta[1] - projection[1];
  offset[2] = delta[2] - projection[2];
  return offset[0]*offset[0] + offset[1]*offset[1] + offset[2]*offset[2];
}

/* ----------------------------------------------------------------------
   allocate data structs that have to wait until sites exist
   so that nlocal,nghost,maxneigh are set
------------------------------------------------------------------------- */

void AppDiffusionCustom::allocate_data()
{
  // for no_energy or linear:
  //   make esites large enough for 2 sites and their 1,2 neighbors
  //   do not need psites
  // for nonlinear:
  //   make esites large enough for 2 sites and their 1,2,3 neighbors
  //   make psites large enough for 2 sites and their 1st neighbors
  // Schwoebel hops add one level of neighbor dependence to esites

  if ((engstyle == NO_ENERGY || engstyle == LINEAR) && hopstyle == NNHOP) {
    int emax = 1 + maxneigh + maxneigh*maxneigh;
    esites = new int[2*emax];
    psites = NULL;
  } else if ((engstyle == NO_ENERGY || engstyle == LINEAR) && 
	     hopstyle == SCHWOEBEL) {
    int emax = 1 + maxneigh + maxneigh*maxneigh + maxneigh*maxneigh*maxneigh;
    esites = new int[2*emax];
    psites = NULL;
  } else if (engstyle == NONLINEAR && hopstyle == NNHOP) {
    int emax = 1 + maxneigh + maxneigh*maxneigh + maxneigh*maxneigh*maxneigh;
    int pmax = 1 + maxneigh;
    esites = new int[2*emax];
    psites = new int[2*pmax];
  } else if (engstyle == NONLINEAR && hopstyle == SCHWOEBEL) {
    int emax = 1 + maxneigh + maxneigh*maxneigh + 
      maxneigh*maxneigh*maxneigh + maxneigh*maxneigh*maxneigh*maxneigh;
    int pmax = 1 + maxneigh;
    esites = new int[2*emax];
    psites = new int[2*pmax];
  }

  echeck = new int[nlocal+nghost];
  pcheck = new int[nlocal+nghost];

  memory->create(firstevent,nlocal,"app:firstevent");

  ecoord = new double[maxneigh+1];
  for (int i = 0; i <= maxneigh; i++) ecoord[i] = 0.0;

  memory->create(hbarrier,maxneigh+1,maxneigh+1,"app:hbarrier");
  memory->create(sbarrier,maxneigh+1,maxneigh+1,"app:sbarrier");

  for (int i = 0; i <= maxneigh; i++)
    for (int j = 0; j <= maxneigh; j++)
      hbarrier[i][j] = sbarrier[i][j] = 0.0;

  hopsite = new int[maxneigh*maxneigh + maxneigh];
  neigh_check = new int[nlocal +nghost];

  marklist = new int[maxneigh*maxneigh];

  mark = NULL;
  if (hopstyle == SCHWOEBEL) memory->create(mark,nlocal+nghost,"app:mark");
  if (mark)
    for (int i = 0; i < nlocal+nghost; i++) mark[i] = 0;

  box_dims = new double[3];
  box_dims[0] =  domain->boxxhi-domain->boxxlo;
  box_dims[1] = domain->boxyhi-domain->boxylo;
  box_dims[2] = domain->boxzhi-domain->boxzlo;
  half_box_dims = new double[3];
  half_box_dims[0] = box_dims[0] * 0.5;
  half_box_dims[1] = box_dims[1] * 0.5;
  half_box_dims[2] = box_dims[2] * 0.5;

}
