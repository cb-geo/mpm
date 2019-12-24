#ifndef MPM_GRAPH_H_
#define MPM_GRAPH_H_

#include <memory>
#include <vector>

// MPI
#ifdef USE_MPI
#include "mpi.h"
#endif

#ifdef USE_PARMETIS
#include <parhip_interface.h>

#include "cell.h"
#include "container.h"
#include "particle.h"

namespace mpm {

//! Base class of graph
//! \brief Base class that stores the information about graph
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Graph {
 public:
  //! Constructor with cells, size and rank
  Graph(Container<Cell<Tdim>> cells, int size, int rank);

  //! Destructor
  ~Graph();

  //! Create graph partition
  bool create_partitions(MPI_Comm* comm);

  //! Collect partitions
  void collect_partitions(int ncells, int mpi_size, int rank, MPI_Comm* comm);

  //! Return xadj
  idxtype* xadj();

  //! Return adjncy
  idxtype* adjncy();

  //! Return vtxdist
  idxtype* vtxdist();

  //! Return vwgt
  idxtype* vwgt();

  //! Tdim
  void assign_ndims(idxtype a);

  //! Return nparts
  int nparts();

 private:
  // Container of cells
  Container<Cell<Tdim>> cells_;

  idxtype numflag_ = 0;
  idxtype wgtflag_ = 2;

  idxtype ncon_ = 0;
  int nparts_ = 0;
  idxtype options_[1];
  idxtype ndims_ = 0;
  int edgecut_ = 0;

  // Array that stores the weights of the adjacency lists
  idxtype* adjwgt_ = nullptr;
  idxtype nvtxs_ = 0;
  std::vector<mpm::Index> part_;
  idxtype* partition_ = nullptr;

  // Pointers to the locally stored vertices
  idxtype* xadj_ = nullptr;
  // Vertex weights
  idxtype* vwgt_ = nullptr;
  // Array that stores the adjacency lists of nvtxs
  idxtype* adjncy_ = nullptr;
  // Distribution of vertices
  idxtype* vtxdist_ = nullptr;
};  // namespace graph
}  // namespace mpm

#include "graph.tcc"
#endif

#endif  // MPM_GRAPH_H_
