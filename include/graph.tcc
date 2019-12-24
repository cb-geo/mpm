//! Constructor with cells, size and rank
template <unsigned Tdim>
mpm::Graph<Tdim>::Graph(Container<Cell<Tdim>> cells, int mpi_size,
                        int mpi_rank) {

  this->cells_ = cells;
  //! Basic parameters used in ParMETIS
  std::vector<idxtype> vxadj;
  std::vector<idxtype> vadjncy;
  std::vector<idxtype> vvtxdist;
  std::vector<idxtype> vvwgt;

  //! There is no weight to adjwgt
  this->adjwgt_ = {};
  //! There is only one weight of one vertex
  this->ncon_ = 1;

  //! Use default value to fill the options[1]
  this->options_[0] = 0;

  idxtype sum = cells_.size();

  idxtype part = 0;
  part = sum / mpi_size;
  idxtype rest = sum % mpi_size;
  if (rest != 0) part = part + 1;

  idxtype start = 0;
  vvtxdist.push_back(start);
  start = start + part;
  //! Insert the local cells for each processor
  if (sum != 1 && part != 0) {
    while (start < sum) {
      vvtxdist.push_back(start);
      start = start + part;
    }
  }

  //! If the numbr of processor can not be evenly distributed, then the last
  //! processor handle the rest of cells
  if (rest != 0) {
    start = start - part;
    start = sum;
    vvtxdist.push_back(start);
  } else {
    vvtxdist.push_back(start);
  }

  vxadj.push_back(0);

  long long counter = 0;
  start = vvtxdist[mpi_rank];
  idxtype end = vvtxdist[mpi_rank + 1];

  for (auto stcl = cells_.cbegin(); stcl != cells_.cend(); ++stcl) {

    if ((*stcl)->id() >= start && (*stcl)->id() < end) {
      if (counter == 0) this->ndims_ = (*stcl)->centroid().rows();

      //! Insert the offset of the size of cell's neighbour
      counter += (*stcl)->nneighbours();

      vxadj.push_back(counter);

      //! get the neighbours
      auto neighbours = (*stcl)->neighbours();

      //! get the id of neighbours
      for (auto neighbour : neighbours) vadjncy.push_back(neighbour);

      vvwgt.push_back((*stcl)->nparticles());
    }
  }

  idxtype* final_xadj = (idxtype*)malloc(vxadj.size() * sizeof(idxtype));
  idxtype* final_adjncy = (idxtype*)malloc(vadjncy.size() * sizeof(idxtype));
  idxtype* final_vtxdist = (idxtype*)malloc(vvtxdist.size() * sizeof(idxtype));
  idxtype* final_vwgt = (idxtype*)malloc(vvwgt.size() * sizeof(idxtype));
  //! Assign the value
  for (long i = 0; i < vxadj.size(); ++i) final_xadj[i] = vxadj.at(i);
  for (long i = 0; i < vadjncy.size(); ++i) final_adjncy[i] = vadjncy.at(i);
  for (long i = 0; i < vvtxdist.size(); ++i) final_vtxdist[i] = vvtxdist.at(i);
  for (long i = 0; i < vvwgt.size(); ++i) final_vwgt[i] = vvwgt.at(i);

  //! Assign the pointer
  this->adjncy_ = final_adjncy;
  this->xadj_ = final_xadj;
  this->vtxdist_ = final_vtxdist;
  this->vwgt_ = final_vwgt;
  std::vector<idxtype>(vadjncy).swap(vadjncy);
  std::vector<idxtype>(vxadj).swap(vxadj);
  std::vector<idxtype>(vvtxdist).swap(vvtxdist);
  std::vector<idxtype>(vvwgt).swap(vvwgt);

  //! assign nparts
  //! nparts is different from mpi_size, but here we can set them equal
  nparts_ = mpi_size;

  //! nvtxs
  this->nvtxs_ = vtxdist_[mpi_rank + 1] - vtxdist_[mpi_rank];

  //! allocate space for part
  part_.reserve(cells_.size());

  //! assign edgecut
  this->edgecut_ = 0;
}

//! Destructor
template <unsigned Tdim>
mpm::Graph<Tdim>::~Graph() {
  part_.clear();
}

//! Get the xadj
template <unsigned Tdim>
idxtype* mpm::Graph<Tdim>::xadj() {
  return this->xadj_;
}

//! Get the adjncy
template <unsigned Tdim>
idxtype* mpm::Graph<Tdim>::adjncy() {
  return this->adjncy_;
}
//! Get the vtxdist
template <unsigned Tdim>
idxtype* mpm::Graph<Tdim>::vtxdist() {
  return this->vtxdist_;
}

//! Get the vwgt
template <unsigned Tdim>
idxtype* mpm::Graph<Tdim>::vwgt() {
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
bool mpm::Graph<Tdim>::create_partitions(MPI_Comm* comm) {

  // The amount of imbalance that is allowed. (3%)
  double imbalance = 0.03;

  int mode = 4;  // Fast Mode
  int seed = 0;

  // Suppress output from the partitioning library.
  bool suppress_output = true;

  ParHIPPartitionKWay(this->vtxdist(), this->xadj(), this->adjncy(),
                      this->vwgt(), this->adjwgt_, &this->nparts_, &imbalance,
                      suppress_output, seed, mode, &this->edgecut_,
                      this->part_.data(), comm);
  return true;
}

//! Collect the partitions and store it in the graph
template <unsigned Tdim>
void mpm::Graph<Tdim>::collect_partitions(int ncells, int mpi_size,
                                          int mpi_rank, MPI_Comm* comm) {
  //! allocate space to partition
  MPI_Status status;
  this->partition_ = (idxtype*)malloc(ncells * sizeof(idxtype));

  if (mpi_rank == 0) {
    int par = 0;
    for (int i = 0; i < this->vtxdist_[1]; ++i) {
      this->partition_[par] = this->part_[i];
      par = par + 1;
    }

    for (int penum = 1; penum < mpi_size; ++penum) {
      idxtype rnvtxs = this->vtxdist_[penum + 1] - this->vtxdist_[penum];
      idxtype* rpart = (idxtype*)malloc(rnvtxs * sizeof(idxtype));
      //! penum is the source process
      MPI_Recv((void*)rpart, rnvtxs, MPI_UNSIGNED_LONG_LONG, penum, 1, *comm,
               &status);

      for (int i = 0; i < rnvtxs; ++i) {
        this->partition_[par] = rpart[i];
        par = par + 1;
      }
      free(rpart);
    }

    for (int penum = 1; penum < mpi_size; ++penum)
      MPI_Send((void*)this->partition_, ncells, MPI_UNSIGNED_LONG_LONG, penum,
               1, *comm);

  } else {
    MPI_Send(this->part_.data(),
             this->vtxdist_[mpi_rank + 1] - this->vtxdist_[mpi_rank],
             MPI_UNSIGNED_LONG_LONG, 0, 1, *comm);
    MPI_Recv((void*)this->partition_, ncells, MPI_UNSIGNED_LONG_LONG, 0, 1,
             *comm, &status);
  }

  // Assign partition to cells
  for (auto citr = this->cells_.cbegin(); citr != this->cells_.cend(); ++citr)
    (*citr)->rank(this->partition_[(*citr)->id()]);

  free(this->partition_);
}
