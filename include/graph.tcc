  //! Constructor with cells, size and rank
template <unsigned Tdim>
mpm::Graph<Tdim>::Graph(Container<Cell<Tdim>> cells, int mpi_size,
                        int mpi_rank) {

  this->cells_ = cells;
  //! Basic parameters used in ParMETIS
  std::vector<idx_t> vxadj;
  std::vector<idx_t> vadjncy;
  std::vector<idx_t> vvtxdist;
  std::vector<idx_t> vvwgt;

  //! There is no weight to adjwgt
  this->adjwgt_ = {};
  //! There is only one weight of one vertex
  this->ncon_ = 1;

  //! Use default value to fill the options[1]
  this->options_[0] = 0;

  long sum = cells_.size();

  long part = 0;
  part = sum / mpi_size;
  long rest = sum % mpi_size;
  if (rest != 0) part = part + 1;

  long start = 0;
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
  idx_t end = vvtxdist[mpi_rank + 1];

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

  idx_t* final_xadj = (idx_t*)malloc(vxadj.size() * sizeof(idx_t));
  idx_t* final_adjncy = (idx_t*)malloc(vadjncy.size() * sizeof(idx_t));
  idx_t* final_vtxdist = (idx_t*)malloc(vvtxdist.size() * sizeof(idx_t));
  idx_t* final_vwgt = (idx_t*)malloc(vvwgt.size() * sizeof(idx_t));
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
  std::vector<idx_t>(vadjncy).swap(vadjncy);
  std::vector<idx_t>(vxadj).swap(vxadj);
  std::vector<idx_t>(vvtxdist).swap(vvtxdist);
  std::vector<idx_t>(vvwgt).swap(vvwgt);

  //! assign ubvec (ParMETIS suggests 1.05)
  for (int nncon = 0; nncon < MAXNCON; ++nncon) ubvec_[nncon] = 1.05;

  //! assign nparts
  //! nparts is different from mpi_size, but here we can set them equal
  nparts_ = mpi_size;

  //! assign tpwgts
  std::vector<real_t> ttpwgts;
  real_t sub_total = 0.0;
  for (int ntpwgts = 0; ntpwgts < ((nparts_) * this->ncon_); ++ntpwgts) {
    if (ntpwgts != (nparts_ * this->ncon_) - 1) {
      ttpwgts.push_back(1.0 / (real_t)nparts_);
      sub_total = sub_total + 1.0 / (real_t)nparts_;
    } else {
      ttpwgts.push_back(1.0 - sub_total);
    }
  }
  real_t* mtpwts = (real_t*)malloc(ttpwgts.size() * sizeof(real_t));
  for (int ntpwgts = 0; ntpwgts < ttpwgts.size(); ++ntpwgts)
    mtpwts[ntpwgts] = ttpwgts.at(ntpwgts);
  this->tpwgts_ = mtpwts;
  std::vector<real_t>(ttpwgts).swap(ttpwgts);

  //! nvtxs
  this->nvtxs_ = vtxdist_[mpi_rank + 1] - vtxdist_[mpi_rank];

  //! assign xyz
  std::vector<real_t> mxyz;

  for (auto stcl = cells_.cbegin(); stcl != cells_.cend(); ++stcl) {
    if ((*stcl)->id() >= start && (*stcl)->id() < end) {
      for (int dimension = 0; dimension < (*stcl)->centroid().rows();
           dimension++) {
        mxyz.push_back(((*stcl)->centroid())(dimension, 0));
      }
    }
  }

  real_t* txyz = (real_t*)malloc(mxyz.size() * sizeof(real_t));
  for (long i = 0; i < mxyz.size(); ++i) txyz[i] = mxyz.at(i);
  this->xyz_ = txyz;

  std::vector<real_t>(mxyz).swap(mxyz);

  //! allocate space for part
  this->part_ = (idx_t*)malloc(this->nvtxs_ * sizeof(idx_t));
  for (int mpart = 0; mpart < this->nvtxs_; ++mpart)
    this->part_[mpart] = mpi_rank % this->nparts_;

  //! assign edgecut
  this->edgecut_ = 0;
}

//! Get the xadj
template <unsigned Tdim>
idx_t* mpm::Graph<Tdim>::xadj() {
  return this->xadj_;
}

//! Get the adjncy
template <unsigned Tdim>
idx_t* mpm::Graph<Tdim>::adjncy() {
  return this->adjncy_;
}
//! Get the vtxdist
template <unsigned Tdim>
idx_t* mpm::Graph<Tdim>::vtxdist() {
  return this->vtxdist_;
}

//! Get the vwgt
template <unsigned Tdim>
idx_t* mpm::Graph<Tdim>::vwgt() {
  return this->vwgt_;
}

template <unsigned Tdim>
void mpm::Graph<Tdim>::assign_ndims(idx_t n) {
  this->ndims_ = n;
}

//! Return nparts
template <unsigned Tdim>
idx_t mpm::Graph<Tdim>::nparts() {
  return this->nparts_;
}

//! Return partition
template <unsigned Tdim>
idx_t mpm::Graph<Tdim>::partition(idx_t id) {
  return this->partition_[id];
}

//! Create partition
template <unsigned Tdim>
bool mpm::Graph<Tdim>::create_partitions(MPI_Comm* comm) {
  //! assign part
  ParMETIS_V3_PartGeomKway(
      this->vtxdist(), this->xadj(), this->adjncy(), this->vwgt(),
      this->adjwgt_, &this->wgtflag_, &this->numflag_, &this->ndims_,
      this->xyz_, &this->ncon_, &this->nparts_, this->tpwgts_, this->ubvec_,
      this->options_, &this->edgecut_, this->part_, comm);
  return true;
}

//! Collect the partitions and store it in the graph
template <unsigned Tdim>
void mpm::Graph<Tdim>::collect_partitions(int ncells, int mpi_size,
                                          int mpi_rank, MPI_Comm* comm) {
  //! allocate space to partition
  MPI_Status status;
  this->partition_ = (idx_t*)malloc(ncells * sizeof(idx_t));

  if (mpi_rank == 0) {
    int par = 0;
    for (int i = 0; i < this->vtxdist_[1]; ++i) {
      this->partition_[par] = this->part_[i];
      par = par + 1;
    }

    for (int penum = 1; penum < mpi_size; ++penum) {
      idx_t rnvtxs = this->vtxdist_[penum + 1] - this->vtxdist_[penum];
      idx_t* rpart = (idx_t*)malloc(rnvtxs * sizeof(idx_t));
      //! penum is the source process
      MPI_Recv((void*)rpart, rnvtxs, IDX_T, penum, 1, *comm, &status);

      for (int i = 0; i < rnvtxs; ++i) {
        this->partition_[par] = rpart[i];
        par = par + 1;
      }
      free(rpart);
    }

    for (int penum = 1; penum < mpi_size; ++penum)
      MPI_Send((void*)this->partition_, ncells, IDX_T, penum, 1, *comm);

  } else {
    MPI_Send((void*)this->part_,
             this->vtxdist_[mpi_rank + 1] - this->vtxdist_[mpi_rank], IDX_T, 0,
             1, *comm);
    //！ free space
    free(this->part_);
    MPI_Recv((void*)this->partition_, ncells, IDX_T, 0, 1, *comm, &status);
  }
}
