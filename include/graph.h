#ifndef MPM_GRAPH_H_
#define MPM_GRAPH_H_

#include <memory>
#include <vector>

// MPI
#ifdef USE_MPI
#include "mpi.h"
#endif

#include "parmetis.h"

#include "cell.h"
#include "container.h"
#include "particle.h"

namespace mpm {
const int MAXNCON = 1;

//! Base class of graph
//! \brief Base class that stores the information about graph
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Graph {
 public:
  //! Constructor with cells, n
  Graph(Container<Cell<Tdim>> cells, int size, int rank);

  //! Create graph partition
  bool create_partitions(MPI_Comm* comm);

  //! Collect partitions
  void collect_partitions(int ncells, int npes, int rank, MPI_Comm* comm);

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
  // Container of cells
  Container<Cell<Tdim>> cells_;

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
  // Array that stores the weights of the adjacency lists
  idx_t* adjwgt;
  idx_t nvtxs;
  idx_t* part = nullptr;
  idx_t* partition = nullptr;

  idx_t adptf;
  idx_t optype;
  idx_t gnvtxs, nedges, nobj;
  // Pointers to the locally stored vertices
  idx_t* xadj;
  // Vertex weights
  idx_t* vwgt;
  // Vertex weights
  real_t* nvwgt;
  // Vertex size
  idx_t* vsize;
  // Array that stores the adjacency lists of nvtxs
  idx_t* adjncy;
  // Distribution of vertices
  idx_t* vtxdist;
  // The initial partition of the vertex
  idx_t* home;
};  // namespace graph
}  // namespace mpm

#include "graph.tcc"

#endif  // MPM_GRAPH_H_
