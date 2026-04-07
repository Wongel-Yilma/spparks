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

#ifdef APP_CLASS
AppStyle(diffusion/multiphase/gcn,AppDiffusionMultiphaseGCN)

#else

#ifndef SPK_APP_DIFFUSION_MULTIPHASE_GCN_H
#define SPK_APP_DIFFUSION_MULTIPHASE_GCN_H

#include "app_lattice.h"
#include <map>
#include <set>
#include <utility>
#include <vector>
#include <torch/torch.h>
#include <torch/script.h>

namespace SPPARKS_NS {

class AppDiffusionMultiphaseGCN : public AppLattice {

 public:
  AppDiffusionMultiphaseGCN(class SPPARKS *, int, char **);
  ~AppDiffusionMultiphaseGCN();
  void input_app(char *, int, char **);
  void grow_app();
  void init_app();
  void setup_app();

  double site_energy(int);
  void site_event_rejection(int, class RandomPark *);
  double site_propensity(int);
  void site_event(int, class RandomPark *);

 private:
  int engstyle;
  int allocated;
  int *esites;
  int *echeck;
  int dimension;
  int *lattice;
  int *type;
  double *x_md;
  double *y_md;
  double *z_md;
  std::vector<std::vector<double>> first_shell_coords;
  std::vector<std::vector<int64_t>> edge_index_vec;
  std::vector<std::vector<double>> destinations;
  // torch::Device device;
  torch::Tensor batch;
  torch::jit::script::Module gcn;
  std::vector<int64_t> edge_index_linearized;
  torch::Tensor edge_index;
  int64_t gnn_edge_size;
  int64_t gnn_node_size;
  
  ///// Debugging variables for statistics /////
  // int num_evaluations;
  // int num_events_removed;
  // std::vector<int> removed_events;
  // std::vector<int> active_events;
  // std::vector<int> all_events_generated;
  // std::vector<double>eb_values;
  //////////////////////

  double *box_dims;
  double *half_box_dims;
  double *box_dims_md;

  struct Event {           // one event for an owned site
    double propensity;     // propensity of this event
    int destination;       // local ID of destination site
    int next;              // index of next event for this site
  };
  
  Event *events;           // list of events for all owned sites
  int nevents;             // # of events for all owned sites
  int maxevent;            // max # of events list can hold
  int *firstevent;         // index of 1st event for each owned site
  int freeevent;           // index of 1st unused event in list
  int *neigh_check;         // list of neighbors 
  int *hopsite;         // list of possible hops for one site 
  
  // phases and pairwise weights used for site energy calculation
  std::set<int> phase_labels;
  std::map<int,bool> is_pinned;
  std::map<std::pair<int,int>,double> weights;
  void parse_diffmultiphase(int narg, char **arg);
  
  void site_event_linear(int, class RandomPark*);
  void clear_events(int);
  void add_event(int, int, double);

  /* ----------------------------------------------------------------------
  Function and Variable definitions for jump rate calculation
  ---------------------------------------------------------------------- */

  // jump rate calculation
  double site_propensity_linear(int);
  // Class variables for jump rate calculation
  std::vector<int> local_neigh_check;
  std::vector<int> neigh_types;
  std::vector<double> xi;
  std::vector<double> xmdi;
  std::vector<double> kmc_coords;
  std::vector<double> md_coords;
  // helper functions for jump rate calculation
  double calculate_barrier_energy(int, int, int );
  // variables definition for calculating the barrier energy
  std::vector<int64_t> first_shell;
  std::vector<double> new_x;
  std::vector<double> new_y;
  std::vector<double> new_z;
  std::vector<std::vector<double>> original_cs;
  std::vector<std::vector<double>> rotated_cs;
  std::vector<std::vector<double>> R;
  std::vector<double> rotated_coords;
  std::vector<double> selected_coords;
  std::vector<float> edge_attr_vec;
  // helper functions for calculating the barrier energy
  void rotate_vectors( int num_vectors);
  void reset_vectors();
  double vec_dot(const std::vector<double> &,const std::vector<double> &);
  double vec_norm(const std::vector<double> &);
  std::vector<double> vec_cross(const std::vector<double> &,const std::vector<double> &);
  double calculate_distance(int , int);
  std::vector<int64_t> linearize_int(const std::vector<std::vector<int64_t>> &);
  
  /* ----------------------------------------------------------------------
    Function and Variable definitions for Deposition of H or Vacancy
  ---------------------------------------------------------------------- */
  void app_update(double );
  void shuffle_indices(std::vector<int>&, class RandomPark*);
  int depmode;
  int ndeposit_total;
  int cummulative_ndeposit;
  int maxbatch;
  double deposit_start_time;
  double sim_progress;
  class RandomPark *rand_deposition;
  struct DepInfo {
    int proc;
    int site;
  };
  DepInfo *depinfo;
  DepInfo *depinfo_global;
  DepInfo *depinfo_selected;
  

  
  void allocate_data();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Cannot use %s command until sites exist

This command requires sites exist before using it in an input script.

E: Can only use ecoord command with app_style diffusion nonlinear

Self-explanatory.

E: Cannot define Schwoebel barrier without Schwoebel model

Self-explanatory.

E: Unrecognized command

The command is assumed to be application specific, but is not
known to SPPARKS.  Check the input script.

E: Cannot perform deposition in parallel

UNDOCUMENTED

E: Cannot perform deposition with multiple sectors

UNDOCUMENTED

E: One or more sites have invalid values

The application only allows sites to be initialized with specific
values.

E: Did not reach event propensity threshhold

UNDOCUMENTED

E: BAD DONE

UNDOCUMENTED

*/
