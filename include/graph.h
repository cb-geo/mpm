#ifndef MPM_GRAPH_H_
#define MPM_GRAPH_H_

#include <algorithm>
#include <array>
#include <limits>
#include <memory>
#include <numeric>
#include <vector>

// Eigen
#include "Eigen/Dense"
// MPI
#ifdef USE_MPI
#include "mpi.h"
#endif
// TBB
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>

#include <tsl/robin_map.h>

//#include "../libparmetis/defs.h"
//#include "../libparmetis/macros.h"
//#include "../libparmetis/proto.h"
//#include "../libparmetis/rename.h"
//#include "../libparmetis/struct.h"
//#include "../metis/GKlib/gk_mkmemory.h"
//#include "../metis/libmetis/gklib_defs.h"
//#include <GKlib.h>
#include <parmetis.h>

//#include "proto.h"

#define MAXNCON 1
#define PMV3_OPTION_DBGLVL  1
#define PMV3_OPTION_SEED  2


#include "cell.h"
#include "container.h"
#include "factory.h"
#include "geometry.h"
#include "hdf5.h"
#include "logger.h"
#include "material/material.h"
#include "node.h"
#include "particle.h"
#include "particle_base.h"

namespace mpm {
// using Index = unsigned long long;

//! Base class of graph
//! \brief Base class that stores the information about graph
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Graph {
 public:
  //! Constructor
  Graph();

  Graph(Container<Cell<Tdim>>* cells, int num_threads);

  Graph& operator=(const Graph& graph);

  //! Initialize the graph
  void initialize(Container<Cell<Tdim>>* cells, int num_threads);

  //! Get the xadj
  idx_t* get_xadj();

  //! Get the adjncy
  idx_t* get_adjncy();

  //! Get the vtxdist
  idx_t* get_vtxdist();

  //! Get the vwgt
  idx_t* get_vwgt();

  //! Change the option value
  void change_options_PMV3_OPTION_DBGLVL(idx_t a);

  //! Change the option value
  void change_options_PMV3_OPTION_SEED(idx_t a);

  //! Chane the option value
  void change_options_0(idx_t a);

  //! nparts
  void assign_nparts(idx_t nparts);

  //! Tdim
  void assign_ndims(idx_t a);

  //! optype
  void assign_optype(idx_t a);

  //! adptf
  void assign_adptf(idx_t a);

  //! ipc2redist
  void assign_ipc2redist(real_t a);
  //! ndims
  idx_t get_ndims() { return ndims; }

  idx_t* part;
  idx_t* sizes;

  idx_t get_nvtxs();
  idx_t get_nparts();

  idx_t numflag = 0;
  idx_t wgtflag = 2;

  idx_t ncon;
  idx_t nparts;
  real_t ubvec[MAXNCON];
  idx_t options[10];
  real_t* xyz = NULL;
  idx_t ndims;

  real_t* tpwgts = NULL;
  idx_t* adjwgt; /* Array that stores the weights of the adjacency lists */
  idx_t nvtxs;

 private:
  Container<Cell<Tdim>>* cells_;

  real_t ipc2resit;

  idx_t adptf;
  idx_t optype;
  idx_t edgecut;
  idx_t gnvtxs, nedges, nobj;
  idx_t* xadj;    /* Pointers to the locally stored vertices */
  idx_t* vwgt;    /* Vertex weights */
  real_t* nvwgt;  /* Vertex weights */
  idx_t* vsize;   /* Vertex size */
  idx_t* adjncy;  /* Array that stores the adjacency lists of nvtxs */
  idx_t* vtxdist; /* Distribution of vertices */
  idx_t* home;    /* The initial partition of the vertex */
};                // namespace graph
}  // namespace mpm

#include "graph.tcc"

#endif  // MPM_GRAPH_H_
