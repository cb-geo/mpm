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

  //! Return xadj
  idx_t* xadj();

  //! Return adjncy
  idx_t* adjncy();

  //! Return vtxdist
  idx_t* vtxdist();

  //! Return vwgt
  idx_t* vwgt();

  //! Tdim
  void assign_ndims(idx_t a);

  //! Return nparts
  idx_t nparts();

  //! partition
  idx_t* partition();

 private:
  // Container of cells
  Container<Cell<Tdim>> cells_;

  idx_t numflag_ = 0;
  idx_t wgtflag_ = 2;

  idx_t ncon_;
  idx_t nparts_;
  real_t ubvec_[MAXNCON];
  idx_t options_[1];
  real_t* xyz_ = nullptr;
  idx_t ndims_;
  idx_t edgecut_ = 0;

  real_t* tpwgts_ = nullptr;
  // Array that stores the weights of the adjacency lists
  idx_t* adjwgt_;
  idx_t nvtxs_;
  idx_t* part_ = nullptr;
  idx_t* partition_ = nullptr;

  // Pointers to the locally stored vertices
  idx_t* xadj_;
  // Vertex weights
  idx_t* vwgt_;
  // Vertex weights
  real_t* nvwgt;
  // Vertex size
  idx_t* vsize;
  // Array that stores the adjacency lists of nvtxs
  idx_t* adjncy_;
  // Distribution of vertices
  idx_t* vtxdist_;
};  // namespace graph
}  // namespace mpm

#include "graph.tcc"

#endif  // MPM_GRAPH_H_
