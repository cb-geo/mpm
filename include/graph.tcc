// Constructor with id
template <unsigned Tdim>
mpm::Graph<Tdim>::Graph(Container<Cell<Tdim>>* cells, int num_threads) {
  // initialize(cells,num_threads);
  this->cells_ = cells;
  // int num = num_threads;
  //! Basic parameters used in ParMETIS
  std::vector<idx_t> xadj;
  std::vector<idx_t> adjncy;
  std::vector<idx_t> vtxdist;
  std::vector<idx_t> vwgt;

  //! Insert the 0th position
  xadj.push_back(0);
  //! Iterate through all the cells
  long long counter = 0;
  for (auto stcl = cells_->cbegin(); stcl != cells_->cend(); ++stcl) {
    //! Insert the offset of the size of cell's neighbour
    counter += (*stcl)->nneighbours();
    xadj.push_back(counter);

    //! Insert the id of neighbours of cell
    //! Iterate through the neighbours

    auto neighbours = (*stcl)->get_neighbours_();
    for (auto neigh = neighbours->begin(); neigh != neighbours->end();
         ++neigh) {
      //! Here key is the id
      adjncy.push_back(neigh.key());
    }

    //! Insert the weights of each vertex(number of particles in each cell)
    vwgt.push_back((*stcl)->nparticles());
  }

  //! Assign the processor and assign value to vtxdist
  long sum = cells_->size();
  long part = sum / num_threads;
  long rest = sum % num_threads;
  long start = 0;
  vtxdist.push_back(start);
  start = start + part;
  //! Insert the loal cells for each processor
  while (start < sum) {
    vtxdist.push_back(start);
  }
  //! If the numbr of processor can not be evenly distributed, then the last
  //! processor handle the rest of cells
  if (rest != 0) {
    start = start - part;
    start = start + rest;
    vtxdist.push_back(start);
  }
  //! Initialize the non-dynamic array
  idx_t final_xadj[xadj.size()];
  idx_t final_adjncy[adjncy.size()];
  idx_t final_vtxdist[vtxdist.size()];
  idx_t final_vwgt[vwgt.size()];
  long i;
  //! Assign the value
  for (i = 0; i < xadj.size(); i++) {
    final_xadj[i] = xadj.at(i);
  }
  for (i = 0; i < adjncy.size(); i++) {
    final_adjncy[i] = adjncy.at(i);
  }
  for (i = 0; i < vtxdist.size(); i++) {
    final_vtxdist[i] = vtxdist.at(i);
  }
  for (i = 0; i < vwgt.size(); i++) {
    final_vwgt[i] = vwgt.at(i);
  }

  //! Assign the pointer
  this->adjncy = final_adjncy;
  this->xadj = final_xadj;
  this->vtxdist = final_vtxdist;
  this->vwgt = final_vwgt;
}

//! Operator = definition
template <unsigned Tdim>
mpm::Graph<Tdim>& mpm::Graph<Tdim>::operator=(const Graph& graph) {
  adjncy = graph.adjncy;
  xadj = graph.xadj;
  vwgt = graph.vwgt;
  vtxdist = graph.vtxdist;
  part = graph.part;
  sizes = graph.sizes;

  numflag = graph.numflag;
  wgtflag = graph.wgtflag;

  ncon = graph.ncon;

  nparts = graph.nparts;

  int i;
  for (i = 0; i < MAXNCON; i++) {
    ubvec[i] = graph.ubvec[i];
  }

  for (i = 0; i < 10; i++) {
    options[i] = graph.options[i];
  }
  // options = graph.options;

  xyz = graph.xyz;
  ndims = graph.ndims;

  tpwgts = graph.tpwgts;
  adjwgt =
      graph.adjwgt; /* Array that stores the weights of the adjacency lists */

  cells_ = graph.cells_;

  ipc2resit = graph.ipc2resit;

  adptf = graph.adptf;
  optype = graph.optype;
  edgecut = graph.edgecut;
  gnvtxs = graph.gnvtxs;
  nvtxs = graph.nvtxs;
  nedges = graph.nedges;
  nobj = graph.nobj;
  // idx_t *xadj;		        /* Pointers to the locally stored
  // vertices
  // */ idx_t *vwgt;		        /* Vertex weights */
  nvwgt = graph.nvwgt; /* Vertex weights */
  vsize = graph.vsize; /* Vertex size */
  // idx_t *adjncy;	        /* Array that stores the adjacency lists of
  // nvtxs
  // */ idx_t *vtxdist;	        /* Distribution of vertices */
  home = graph.home; /* The initial partition of the vertex */

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
      auto neighbours = (*stcl)->get_neighbours_();

      //! get the id of neighbours
      for (auto neigh = neighbours->cbegin(); neigh != neighbours->cend();
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
  for (i = 0; i < vxadj.size(); i++) {
    final_xadj[i] = vxadj.at(i);
  }
  for (i = 0; i < vadjncy.size(); i++) {
    final_adjncy[i] = vadjncy.at(i);
  }
  for (i = 0; i < vvtxdist.size(); i++) {
    final_vtxdist[i] = vvtxdist.at(i);
  }
  for (i = 0; i < vvwgt.size(); i++) {
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
  for (nncon = 0; nncon < MAXNCON; nncon++) {
    ubvec[nncon] = 1.05;
  }
  //! assign nparts
  //! nparts is different from num_threads, but here we can set them equal
  nparts = num_threads;

  //! assign tpwgts
  std::vector<real_t> ttpwgts;
  int ntpwgts;
  real_t sub_total = 0.0;
  for (ntpwgts = 0; ntpwgts < ((nparts) * this->ncon); ntpwgts++) {
    if (ntpwgts != (nparts * this->ncon) - 1) {
      ttpwgts.push_back(1.0 / (real_t)nparts);
      sub_total = sub_total + 1.0 / (real_t)nparts;
    } else {
      ttpwgts.push_back(1.0 - sub_total);
    }
  }
  real_t* mtpwts = (real_t*)malloc(ttpwgts.size() * sizeof(real_t));
  for (ntpwgts = 0; ntpwgts < ttpwgts.size(); ntpwgts++) {
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
  for (i = 0; i < mxyz.size(); i++) {
    txyz[i] = mxyz.at(i);
  }
  this->xyz = txyz;

  std::vector<real_t>(mxyz).swap(mxyz);
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
  //! Add somethng to ameliorate
  //! We can add some more codes here
}

template <unsigned Tdim>
void mpm::Graph<Tdim>::assign_ndims(idx_t n) {
  this->ndims = n;
}

//! get_nparts
template <unsigned Tdim>
idx_t mpm::Graph<Tdim>::get_nparts() {
  return this->nparts;
}