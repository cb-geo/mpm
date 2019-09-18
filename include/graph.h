#ifndef MPM_GRAPH_H_
#define MPM_GRAPH_H_

#include <memory>
#include <vector>

// MPI
#ifdef USE_MPI
#include "mpi.h"
#endif

#include "cell.h"
#include "container.h"
#include "particle.h"
#include <parmetis.h>

namespace mpm {
const int MAXNCON = 1;

//! Base class of graph
//! \brief Base class that stores the information about graph
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Graph {
 public:
  //! Constructor
  Graph();

  Graph& operator=(const Graph& graph);

  //! Initialize the graph
  void initialize(Container<Cell<Tdim>>* cells, int num_threads, int mype);

  //! Do the partition
  bool make_partition(MPI_Comm* comm);

  //! Do the collection
  void collect_partition(int ncells, int npes, int mype, MPI_Comm* comm);

  //! Get the xadj
  idx_t* get_xadj();

  //! Get the adjncy
  idx_t* get_adjncy();

  //! Get the vtxdist
  idx_t* get_vtxdist();

  //! Get the vwgt
  idx_t* get_vwgt();

  //! Tdim
  void assign_ndims(idx_t a);

  //! Get nparts
  idx_t get_nparts();

  //! Get partition
  idx_t* get_partition();

 private:
  Container<Cell<Tdim>>* cells_;

  real_t ipc2resit;
  idx_t numflag = 0;
  idx_t wgtflag = 2;

  idx_t ncon;
  idx_t nparts;
  real_t ubvec[MAXNCON];
  idx_t options[1];
  real_t* xyz = nullptr;
  idx_t ndims;
  idx_t edgecut = 0;

  real_t* tpwgts = nullptr;
  idx_t* adjwgt; /* Array that stores the weights of the adjacency lists */
  idx_t nvtxs;
  idx_t* part = nullptr;
  idx_t* partition = nullptr;

  idx_t adptf;
  idx_t optype;
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
