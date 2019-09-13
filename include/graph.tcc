//! Operator = definition
template <unsigned Tdim>
mpm::Graph<Tdim>& mpm::Graph<Tdim>::operator=(const Graph& graph) {
  adjncy = graph.adjncy;
  xadj = graph.xadj;
  vwgt = graph.vwgt;
  vtxdist = graph.vtxdist;
  numflag = graph.numflag;
  wgtflag = graph.wgtflag;
  ncon = graph.ncon;
  nparts = graph.nparts;

  //! assign ubvec
  int i;
  for (i = 0; i < MAXNCON; ++i) {
    ubvec[i] = graph.ubvec[i];
  }

  //! assign option
  for (i = 0; i < 10; ++i) {
    options[i] = graph.options[i];
  }
  // options = graph.options;

  xyz = graph.xyz;
  ndims = graph.ndims;

  tpwgts = graph.tpwgts;
  adjwgt =
      graph.adjwgt; /* Array that stores the weights of the adjacency lists */
  cells_ = graph.cells_;
  part = graph.part;
  partition = graph.partition;
  ipc2resit = graph.ipc2resit;
  edgecut = graph.edgecut;
  adptf = graph.adptf;
  optype = graph.optype;
  edgecut = graph.edgecut;
  gnvtxs = graph.gnvtxs;
  nvtxs = graph.nvtxs;
  nedges = graph.nedges;
  nobj = graph.nobj;
  nvwgt = graph.nvwgt; /* Vertex weights */
  vsize = graph.vsize; /* Vertex size */
  home = graph.home;   /* The initial partition of the vertex */

  return *this;
}

//! Initialize the graph
template <unsigned Tdim>
void mpm::Graph<Tdim>::initialize(Container<Cell<Tdim>>* cells, int num_threads,
                                  int mype) {

  this->cells_ = cells;
  //! Basic parameters used in ParMETIS
  std::vector<idx_t> vxadj;
  std::vector<idx_t> vadjncy;
  std::vector<idx_t> vvtxdist;
  std::vector<idx_t> vvwgt;

  //! There is no weight to adjwgt
  this->adjwgt = {};
  //! There is only one weight of one vertex
  this->ncon = 1;

  //! Use default value to fill the options[10]
  this->options[0] = 1;
  //! Can change the value here
  this->options[PMV3_OPTION_DBGLVL] = 6;
  this->options[PMV3_OPTION_SEED] = 1;

  long sum = cells_->size();

  long part = 0;
  part = sum / num_threads;
  long rest = sum % num_threads;
  if (rest != 0) {
    part = part + 1;
  }

  long start = 0;
  vvtxdist.push_back(start);
  start = start + part;
  //! Insert the loal cells for each processor
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
  start = vvtxdist[mype];
  idx_t end = vvtxdist[mype + 1];

  for (auto stcl = cells_->cbegin(); stcl != cells_->cend(); ++stcl) {

    if ((*stcl)->id() >= start && (*stcl)->id() < end) {
      if (counter == 0) {
        this->ndims = (*stcl)->centroid().rows();
      }
      //! Insert the offset of the size of cell's neighbour
      counter += (*stcl)->nneighbours();

      vxadj.push_back(counter);

      //! get the neighbours
      auto neighbours = (*stcl)->neighbours();

      //! get the id of neighbours
      for (auto neigh = neighbours.cbegin(); neigh != neighbours.cend();
           ++neigh) {
        vadjncy.push_back((*neigh));
      }

      vvwgt.push_back((*stcl)->nparticles());
    }
  }

  idx_t* final_xadj = (idx_t*)malloc(vxadj.size() * sizeof(idx_t));
  idx_t* final_adjncy = (idx_t*)malloc(vadjncy.size() * sizeof(idx_t));
  idx_t* final_vtxdist = (idx_t*)malloc(vvtxdist.size() * sizeof(idx_t));
  idx_t* final_vwgt = (idx_t*)malloc(vvwgt.size() * sizeof(idx_t));
  long i;
  //! Assign the value
  for (i = 0; i < vxadj.size(); ++i) {
    final_xadj[i] = vxadj.at(i);
  }
  for (i = 0; i < vadjncy.size(); ++i) {
    final_adjncy[i] = vadjncy.at(i);
  }
  for (i = 0; i < vvtxdist.size(); ++i) {
    final_vtxdist[i] = vvtxdist.at(i);
  }
  for (i = 0; i < vvwgt.size(); ++i) {
    final_vwgt[i] = vvwgt.at(i);
  }

  //! Assign the pointer
  this->adjncy = final_adjncy;
  this->xadj = final_xadj;
  this->vtxdist = final_vtxdist;
  this->vwgt = final_vwgt;
  std::vector<idx_t>(vadjncy).swap(vadjncy);
  std::vector<idx_t>(vxadj).swap(vxadj);
  std::vector<idx_t>(vvtxdist).swap(vvtxdist);
  std::vector<idx_t>(vvwgt).swap(vvwgt);

  //! assign ubvec
  int nncon = 0;

  //! The guide suggests 1.05
  for (nncon = 0; nncon < MAXNCON; ++nncon) ubvec[nncon] = 1.05;
  //! assign nparts
  //! nparts is different from num_threads, but here we can set them equal
  nparts = num_threads;

  //! assign tpwgts
  std::vector<real_t> ttpwgts;
  int ntpwgts;
  real_t sub_total = 0.0;
  for (ntpwgts = 0; ntpwgts < ((nparts) * this->ncon); ++ntpwgts) {
    if (ntpwgts != (nparts * this->ncon) - 1) {
      ttpwgts.push_back(1.0 / (real_t)nparts);
      sub_total = sub_total + 1.0 / (real_t)nparts;
    } else {
      ttpwgts.push_back(1.0 - sub_total);
    }
  }
  real_t* mtpwts = (real_t*)malloc(ttpwgts.size() * sizeof(real_t));
  for (ntpwgts = 0; ntpwgts < ttpwgts.size(); ++ntpwgts) {
    mtpwts[ntpwgts] = ttpwgts.at(ntpwgts);
  }
  this->tpwgts = mtpwts;
  std::vector<real_t>(ttpwgts).swap(ttpwgts);

  //! nvtxs
  this->nvtxs = vtxdist[mype + 1] - vtxdist[mype];

  //! assign xyz
  std::vector<real_t> mxyz;

  for (auto stcl = cells_->cbegin(); stcl != cells_->cend(); ++stcl) {
    if ((*stcl)->id() >= start && (*stcl)->id() < end) {
      int dimension = 0;
      for (dimension = 0; dimension < (*stcl)->centroid().rows(); dimension++) {
        mxyz.push_back(((*stcl)->centroid())(dimension, 0));
      }
    }
  }

  real_t* txyz = (real_t*)malloc(mxyz.size() * sizeof(real_t));
  for (i = 0; i < mxyz.size(); ++i) {
    txyz[i] = mxyz.at(i);
  }
  this->xyz = txyz;

  std::vector<real_t>(mxyz).swap(mxyz);

  //! allocate space for part
  this->part = (idx_t*)malloc(this->nvtxs * sizeof(idx_t));
  int mpart = 0;
  for (mpart = 0; mpart < this->nvtxs; ++mpart) {
    this->part[mpart] = mype % this->nparts;
  }

  //! assign edgecut
  this->edgecut = 0;
}

//! Get the xadj
template <unsigned Tdim>
idx_t* mpm::Graph<Tdim>::get_xadj() {
  return this->xadj;
}

//! Get the adjncy
template <unsigned Tdim>
idx_t* mpm::Graph<Tdim>::get_adjncy() {
  return this->adjncy;
}
//! Get the vtxdist
template <unsigned Tdim>
idx_t* mpm::Graph<Tdim>::get_vtxdist() {
  return this->vtxdist;
}

//! Get the vwgt
template <unsigned Tdim>
idx_t* mpm::Graph<Tdim>::get_vwgt() {
  return this->vwgt;
}

template <unsigned Tdim>
mpm::Graph<Tdim>::Graph() {
  //! Add something to ameliorate
  //! We can add some more codes here
  this->adjwgt = nullptr;
  this->xadj = nullptr;
  this->vwgt = nullptr;
  this->nvwgt = nullptr;
  this->vsize = nullptr;
  this->adjncy = nullptr;
  this->vtxdist = nullptr;
  this->home = nullptr;
  this->part = nullptr;
  this->partition = nullptr;
  this->edgecut = 0;
}

template <unsigned Tdim>
void mpm::Graph<Tdim>::assign_ndims(idx_t n) {
  this->ndims = n;
}

//! get nparts
template <unsigned Tdim>
idx_t mpm::Graph<Tdim>::get_nparts() {
  return this->nparts;
}

//! get partition
template <unsigned Tdim>
idx_t* mpm::Graph<Tdim>::get_partition() {
  return this->partition;
}

//! do the partition
template <unsigned Tdim>
bool mpm::Graph<Tdim>::make_partition(MPI_Comm* comm) {
  //! assign part
  ParMETIS_V3_PartGeomKway(
      this->get_vtxdist(), this->get_xadj(), this->get_adjncy(),
      this->get_vwgt(), this->adjwgt, &this->wgtflag, &this->numflag,
      &this->ndims, this->xyz, &this->ncon, &this->nparts, this->tpwgts,
      this->ubvec, this->options, &this->edgecut, this->part, comm);
  return true;
}

//! collect the partition and store it in the graph
template <unsigned Tdim>
void mpm::Graph<Tdim>::collect_partition(int ncells, int npes, int mype,
                                         MPI_Comm* comm) {
  //! allocate space to partition
  MPI_Status status;
  this->partition = (idx_t*)malloc(ncells * sizeof(idx_t));
  int penum;
  if (mype == 0) {
    int par = 0;
    int i;
    for (i = 0; i < this->vtxdist[1]; ++i) {
      this->partition[par] = this->part[i];
      par = par + 1;
    }
    for (penum = 1; penum < npes; ++penum) {
      idx_t rnvtxs = this->vtxdist[penum + 1] - this->vtxdist[penum];
      idx_t* rpart;
      rpart = (idx_t*)malloc(rnvtxs * sizeof(idx_t));
      //! penum is the source process
      MPI_Recv((void*)rpart, rnvtxs, IDX_T, penum, 1, *comm, &status);
      int i;
      for (i = 0; i < rnvtxs; ++i) {
        this->partition[par] = rpart[i];
        par = par + 1;
      }
      free(rpart);
    }
    for (penum = 1; penum < npes; ++penum) {
      MPI_Send((void*)this->partition, ncells, IDX_T, penum, 1, *comm);
    }

  } else {
    MPI_Send((void*)this->part, this->vtxdist[mype + 1] - this->vtxdist[mype],
             IDX_T, 0, 1, *comm);
    //！ free space
    free(this->part);
    MPI_Recv((void*)this->partition, ncells, IDX_T, 0, 1, *comm, &status);
  }
}