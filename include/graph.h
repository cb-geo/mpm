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

  //! Create graph partition
  bool create_partitions(MPI_Comm* comm);

  //! Collect partitions
  void collect_partitions(int mpi_size, int rank, MPI_Comm* comm);

  //! Return xadj
  std::vector<idxtype> xadj() const;

  //! Return adjncy
  std::vector<idxtype> adjncy() const;

  //! Return vtxdist
  std::vector<idxtype> vtxdist() const;

  //! Return vwgt
  std::vector<idxtype> vwgt() const;

  //! Tdim
  void assign_ndims(idxtype a);

  //! Return nparts
  int nparts();

 private:
  // Container of cells
  Container<Cell<Tdim>> cells_;
  // Number of partitions
  int nparts_ = 0;
  // Number of dimensions
  idxtype ndims_ = 0;
  // Edge cut
  int edgecut_ = 0;

  // Partition ids
  std::vector<mpm::Index> part_;
  // Array that stores the weights of the adjacency lists
  std::vector<idxtype> adjwgt_;
  // Pointers to the locally stored vertices
  std::vector<idxtype> xadj_;
  // Vertex weights
  std::vector<idxtype> vwgt_;
  // Array that stores the adjacency lists of nvtxs
  std::vector<idxtype> adjncy_;
  // Distribution of vertices
  std::vector<idxtype> vtxdist_;
};  // namespace graph
}  // namespace mpm

#include "graph.tcc"
#endif

#endif  // MPM_GRAPH_H_
