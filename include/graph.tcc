//! Constructor with cells, size and rank
template <unsigned Tdim>
mpm::Graph<Tdim>::Graph(Vector<Cell<Tdim>> cells) {
  this->cells_ = cells;
}

//! Constructor with cells, size and rank
template <unsigned Tdim>
void mpm::Graph<Tdim>::construct_graph(int mpi_size, int mpi_rank) {
  // Clear all graph properties
  this->xadj_.clear();
  this->vwgt_.clear();
  this->adjncy_.clear();
  this->vtxdist_.clear();
  this->part_.clear();

  //! There is no adjacency weight (edge weight)
  this->adjwgt_.clear();

  idxtype sum = cells_.size();

  idxtype part = 0;
  part = sum / mpi_size;
  idxtype rest = sum % mpi_size;
  if (rest != 0) part = part + 1;

  idxtype start = 0;
  vtxdist_.emplace_back(start);
  start = start + part;
  //! Insert the local cells for each processor
  if (sum != 1 && part != 0) {
    while (start < sum) {
      vtxdist_.emplace_back(start);
      start = start + part;
    }
  }

  //! If the number of processors can not be evenly distributed, then the last
  //! processor will handle the rest of cells
  if (rest != 0) {
    // start = start - part;
    start = sum;
    vtxdist_.emplace_back(start);
  } else {
    vtxdist_.emplace_back(start);
  }

  this->xadj_.emplace_back(0);

  mpm::Index offset = 0;
  start = vtxdist_[mpi_rank];
  idxtype end = vtxdist_[mpi_rank + 1];

  for (auto citr = cells_.cbegin(); citr != cells_.cend(); ++citr) {

    if ((*citr)->id() >= start && (*citr)->id() < end) {
      if (offset == 0) this->ndims_ = (*citr)->centroid().rows();

      //! Insert the offset of the size of cell's neighbour
      offset += (*citr)->nneighbours();

      this->xadj_.emplace_back(offset);

      //! get the neighbours
      auto neighbours = (*citr)->neighbours();

      //! get the id of neighbours
      for (const auto& neighbour : neighbours) {
        adjncy_.emplace_back(neighbour);
        adjwgt_.emplace_back(1.);
      }
      vwgt_.emplace_back((*citr)->nglobal_particles());
    }
  }

  //! assign nparts
  //! nparts is different from mpi_size, but here we can set them equal
  nparts_ = mpi_size;

  //! allocate space for part
  part_.reserve(cells_.size());
}

//! Return xadj
template <unsigned Tdim>
std::vector<idxtype> mpm::Graph<Tdim>::xadj() const {
  return this->xadj_;
}

//! Return adjncy
template <unsigned Tdim>
std::vector<idxtype> mpm::Graph<Tdim>::adjncy() const {
  return this->adjncy_;
}
//! Return vtxdist
template <unsigned Tdim>
std::vector<idxtype> mpm::Graph<Tdim>::vtxdist() const {
  return this->vtxdist_;
}

//! Return vwgt
template <unsigned Tdim>
std::vector<idxtype> mpm::Graph<Tdim>::vwgt() const {
  return this->vwgt_;
}

template <unsigned Tdim>
void mpm::Graph<Tdim>::assign_ndims(idxtype n) {
  this->ndims_ = n;
}

//! Return nparts
template <unsigned Tdim>
int mpm::Graph<Tdim>::nparts() {
  return this->nparts_;
}

//! Create partition
template <unsigned Tdim>
bool mpm::Graph<Tdim>::create_partitions(MPI_Comm* comm, int mode) {

  // The amount of imbalance that is allowed. (3%)
  double imbalance = 0.03;

  int seed = 0;

  // Suppress output from the partitioning library.
  bool suppress_output = true;

  ParHIPPartitionKWay(
      this->vtxdist_.data(), this->xadj_.data(), this->adjncy_.data(),
      this->vwgt_.data(), this->adjwgt_.data(), &this->nparts_, &imbalance,
      suppress_output, seed, mode, &this->edgecut_, this->part_.data(), comm);
  return true;
}

//! Collect the partitions and store it in the graph
template <unsigned Tdim>
std::vector<mpm::Index> mpm::Graph<Tdim>::collect_partitions(int mpi_size,
                                                             int mpi_rank,
                                                             MPI_Comm* comm) {
  //! allocate space to partition
  mpm::Index ncells = this->cells_.size();
  std::vector<mpm::Index> partition(ncells, 0);
  // ID of cells, which should transfer particles
  std::vector<mpm::Index> exchange_cells;
  if (mpi_rank == 0) {
    int par = 0;
    for (int i = 0; i < this->vtxdist_[1]; ++i) {
      partition[par] = this->part_[i];
      par = par + 1;
    }

    for (int penum = 1; penum < mpi_size; ++penum) {
      idxtype rnvtxs = this->vtxdist_[penum + 1] - this->vtxdist_[penum];
      std::vector<idxtype> rpart(rnvtxs);
      //! penum is the source process
      MPI_Recv(rpart.data(), rnvtxs, MPI_UNSIGNED_LONG_LONG, penum, 1, *comm,
               MPI_STATUS_IGNORE);

      for (int i = 0; i < rnvtxs; ++i) {
        partition[par] = rpart[i];
        par = par + 1;
      }
    }

    for (int penum = 1; penum < mpi_size; ++penum)
      MPI_Send(partition.data(), ncells, MPI_UNSIGNED_LONG_LONG, penum, 1,
               *comm);

  } else {
    MPI_Send(this->part_.data(),
             this->vtxdist_[mpi_rank + 1] - this->vtxdist_[mpi_rank],
             MPI_UNSIGNED_LONG_LONG, 0, 1, *comm);
    MPI_Recv(partition.data(), ncells, MPI_UNSIGNED_LONG_LONG, 0, 1, *comm,
             MPI_STATUS_IGNORE);
  }

  // Assign partition to cells
  for (auto citr = this->cells_.cbegin(); citr != this->cells_.cend(); ++citr) {
    auto current_rank = partition[(*citr)->id()];
    auto previous_rank = (*citr)->rank();
    // If the current rank is different from cell rank
    if (current_rank != previous_rank) {
      // Assign current MPI rank
      (*citr)->rank(current_rank);
      // Add cell id to list of cells to transfer particles if there are
      // particles
      if ((*citr)->nglobal_particles() > 0)
        exchange_cells.emplace_back((*citr)->id());
    }
  }
  return exchange_cells;
}
