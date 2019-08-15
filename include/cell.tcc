//! Constructor with cell id, number of nodes and element
template <unsigned Tdim>
mpm::Cell<Tdim>::Cell(Index id, unsigned nnodes,
                      const std::shared_ptr<const Element<Tdim>>& elementptr,
                      bool isoparametric)
    : id_{id}, nnodes_{nnodes}, isoparametric_{isoparametric} {
  // Check if the dimension is between 1 & 3
  static_assert((Tdim >= 1 && Tdim <= 3), "Invalid global dimension");

  //! Logger
  std::string logger =
      "cell" + std::to_string(Tdim) + "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);

  try {
    if (elementptr->nfunctions() == this->nnodes_) {
      element_ = elementptr;
    } else {
      throw std::runtime_error(
          "Specified number of shape functions and nodes don't match");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
}

//! Return the initialisation status of cells
template <unsigned Tdim>
bool mpm::Cell<Tdim>::initialise() {
  bool status = false;
  try {
    // Check if node pointers are present and are equal to the expected number
    if (this->nnodes_ == this->nodes_.size()) {
      // Initialise cell properties (volume, centroid, length)
      this->compute_volume();
      this->compute_centroid();
      this->compute_mean_length();
      this->compute_nodal_coordinates();

      // Get centroid of a cell in natural coordinates which are zeros
      Eigen::Matrix<double, Tdim, 1> xi_centroid;
      xi_centroid.setZero();

      Eigen::Matrix<double, Tdim, 1> zero;
      zero.setZero();

      // Get B-Matrix at the centroid
      bmatrix_centroid_ =
          element_->bmatrix(xi_centroid, this->nodal_coordinates_, zero, zero);

      status = true;
    } else {
      throw std::runtime_error(
          "Specified number of nodes for a cell is not present");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return status;
}

//! Return the initialisation status of cells
//! \retval initialisation_status Cell has nodes, shape functions and volumes
template <unsigned Tdim>
bool mpm::Cell<Tdim>::is_initialised() const {
  // Check if node pointers are present and are equal to the # shape_fns
  return ((this->nnodes_ != 0 && this->nfunctions() == this->nodes_.size()) &&
          // Check if shape function is assigned
          this->nfunctions() != 0 &&
          // Check if volume of a cell is initialised
          (std::fabs(this->volume_ - std::numeric_limits<double>::max()) >
           1.0E-10) &&
          // Check if mean length of a cell is initialised
          (std::fabs(this->mean_length_ - std::numeric_limits<double>::max()) >
           1.0E-10));
}

//! Assign quadrature
template <unsigned Tdim>
void mpm::Cell<Tdim>::assign_quadrature(unsigned nquadratures) {
  this->quadrature_ = element_->quadrature(nquadratures);
}

//! Assign quadrature
template <unsigned Tdim>
std::vector<Eigen::Matrix<double, Tdim, 1>> mpm::Cell<Tdim>::generate_points() {
  // Assign a default quadrature of 1
  if (this->quadrature_ == nullptr) this->assign_quadrature(1);

  const auto quadratures = quadrature_->quadratures();

  // Vector of gauss points
  std::vector<Eigen::Matrix<double, Tdim, 1>> points;

  // Get indices of corner nodes
  Eigen::VectorXi indices = element_->corner_indices();
  // Matrix of nodal coordinates
  Eigen::MatrixXd nodal_coords;
  nodal_coords.resize(Tdim, indices.size());

  for (unsigned j = 0; j < indices.size(); ++j)
    nodal_coords.col(j) = nodes_[indices(j)]->coordinates();

  // Zeros
  Eigen::Matrix<double, Tdim, 1> zeros = Eigen::Matrix<double, Tdim, 1>::Zero();

  // Get local coordinates of gauss points and transform to global
  for (unsigned i = 0; i < quadratures.cols(); ++i) {
    const auto lpoint = quadratures.col(i);
    // Get shape functions
    const auto sf = element_->shapefn(lpoint, zeros, zeros);
    const auto point = nodal_coords * sf;
    const auto xi = this->transform_real_to_unit_cell(point);
    bool status = true;
    // Check if point is within the cell
    for (unsigned i = 0; i < xi.size(); ++i)
      if (xi(i) < -1. || xi(i) > 1. || std::isnan(xi(i))) status = false;

    if (status)
      points.emplace_back(point);
    else
      console_->warn("Cannot generate point: ({}, {}) in cell xi: ({}, {})",
                     point(0), point(1), xi(0), xi(1));
  }

  return points;
}

//! Add a node pointer and return the status of addition of a node
template <unsigned Tdim>
bool mpm::Cell<Tdim>::add_node(
    unsigned local_id, const std::shared_ptr<mpm::NodeBase<Tdim>>& node_ptr) {
  bool insertion_status = false;
  try {
    // If number of node ptrs in a cell is less than the maximum number of nodes
    // per cell
    // The local id should be between 0 and maximum number of nodes
    if (nodes_.size() < this->nnodes_ &&
        (local_id >= 0 && local_id < this->nnodes_)) {
      nodes_.emplace_back(node_ptr);
      insertion_status = true;
    } else {
      throw std::runtime_error(
          "Number nodes in a cell exceeds the maximum allowed per cell");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return insertion_status;
}

//! Activate nodes if particle is present
template <unsigned Tdim>
void mpm::Cell<Tdim>::activate_nodes() {
  // If number of particles are present, set all associated nodes as active
  if (particles_.size() > 0) {
    // std::lock_guard<std::mutex> guard(cell_mutex_);
    for (unsigned i = 0; i < nodes_.size(); ++i) nodes_[i]->assign_status(true);
  }
}

//! Return a vector of side node id pairs
template <unsigned Tdim>
std::vector<std::array<mpm::Index, 2>> mpm::Cell<Tdim>::side_node_pairs()
    const {
  // Create a vector of node_pairs
  std::vector<std::array<mpm::Index, 2>> node_pairs;
  // Get the indices of sides
  Eigen::MatrixXi indices = element_->sides_indices();
  // Iterate over indices and get node ids
  for (unsigned i = 0; i < indices.rows(); ++i)
    node_pairs.emplace_back(std::array<mpm::Index, 2>(
        {nodes_[indices(i, 0)]->id(), nodes_[indices(i, 1)]->id()}));

  return node_pairs;
}

//! Add a neighbour cell and return the status of addition of a node
template <unsigned Tdim>
bool mpm::Cell<Tdim>::add_neighbour(mpm::Index neighbour_id) {
  bool insertion_status = false;
  try {
    // If cell id is not the same as the current cell
    if (neighbour_id != this->id())
      insertion_status = (neighbours_.insert(neighbour_id)).second;
    else
      throw std::runtime_error("Invalid local id of a cell neighbour");

  } catch (std::exception& exception) {
    console_->error("{} {}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return insertion_status;
}

//! Add a particle id and return the status of addition of a particle id
template <unsigned Tdim>
bool mpm::Cell<Tdim>::add_particle_id(Index id) {
  bool status = false;
  std::lock_guard<std::mutex> guard(cell_mutex_);
  // Check if it is found in the container
  auto itr = std::find(particles_.begin(), particles_.end(), id);
  if (itr == particles_.end()) {
    particles_.emplace_back(id);
    status = true;
  }
  return status;
}

//! Remove a particle id
template <unsigned Tdim>
void mpm::Cell<Tdim>::remove_particle_id(Index id) {
  std::lock_guard<std::mutex> guard(cell_mutex_);
  particles_.erase(std::remove(particles_.begin(), particles_.end(), id),
                   particles_.end());
}

//! Compute volume of a 1D cell
//! Computes the length of cell
template <>
inline void mpm::Cell<1>::compute_volume() {
  this->volume_ =
      std::fabs(nodes_[0]->coordinates()[0] - nodes_[1]->coordinates()[1]);
}

//! Compute volume of a 2D cell
//! Computes the volume of a triangle and a quadrilateral
template <>
inline void mpm::Cell<2>::compute_volume() {
  try {
    Eigen::VectorXi indices = element_->corner_indices();
    // Triangle
    if (indices.size() == 3) {

      //   2 0
      //     |`\
      //     |  `\
      //     |    `\
      //     |      `\
      //     |        `\
      //   0 0----------0 1
      //
      auto node0 = nodes_[indices(0)]->coordinates();
      auto node1 = nodes_[indices(1)]->coordinates();
      auto node2 = nodes_[indices(2)]->coordinates();

      // 2 * Area = (x1 * y2 - x2 * y1) - (x0 * y2 - x2 * y0) + (x0 * y1 - x1 *
      // y0)
      volume_ = std::fabs(((node1(0) * node2(1)) - (node2(0) * node1(1))) -
                          ((node0(0) * node2(1)) - (node2(0) * node0(1))) +
                          ((node0(0) * node1(1)) - (node1(0) * node0(1)))) /
                2.0;
      // Quadrilateral
    } else if (indices.size() == 4) {

      //        b
      // 3 0---------0 2
      //   | \   q / |
      // a |   \  /  | c
      //   |   p \   |
      //   |  /    \ |
      // 0 0---------0 1
      //         d
      const double a = (nodes_[indices(0)]->coordinates() -
                        nodes_[indices(3)]->coordinates())
                           .norm();
      const double b = (nodes_[indices(2)]->coordinates() -
                        nodes_[indices(3)]->coordinates())
                           .norm();
      const double c = (nodes_[indices(1)]->coordinates() -
                        nodes_[indices(2)]->coordinates())
                           .norm();
      const double d = (nodes_[indices(0)]->coordinates() -
                        nodes_[indices(1)]->coordinates())
                           .norm();
      const double p = (nodes_[indices(0)]->coordinates() -
                        nodes_[indices(2)]->coordinates())
                           .norm();
      const double q = (nodes_[indices(1)]->coordinates() -
                        nodes_[indices(3)]->coordinates())
                           .norm();

      // K = 1/4 * sqrt ( 4p^2q^2 - (a^2 + c^2 - b^2 -d^2)^2)
      volume_ =
          0.25 * std::sqrt(4 * p * p * q * q -
                           std::pow((a * a + c * c - b * b - d * d), 2.0));
    } else {
      throw std::runtime_error(
          "Unable to compute volume, number of vertices is incorrect");
    }
    // Check negative volume
    if (this->volume_ <= 0)
      throw std::runtime_error(
          "Negative or zero volume cell, misconfigured cell!");
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
}

//! Compute volume of a 3D cell
//! Computes the volume of a hexahedron
template <>
inline void mpm::Cell<3>::compute_volume() {
  try {
    Eigen::VectorXi indices = element_->corner_indices();
    // Hexahedron
    if (indices.size() == 8) {
      // Node numbering as read in by mesh file
      //        d               c
      //          *_ _ _ _ _ _*
      //         /|           /|
      //        / |          / |
      //     a *_ |_ _ _ _ _* b|
      //       |  |         |  |
      //       |  |         |  |
      //       |  *_ _ _ _ _|_ *
      //       | / h        | / g
      //       |/           |/
      //       *_ _ _ _ _ _ *
      //     e               f
      //
      // Calculation of hexahedron volume from
      // https://arc.aiaa.org/doi/pdf/10.2514/3.9013

      const auto a = nodes_[indices(7)]->coordinates();
      const auto b = nodes_[indices(6)]->coordinates();
      const auto c = nodes_[indices(2)]->coordinates();
      const auto d = nodes_[indices(3)]->coordinates();
      const auto e = nodes_[indices(4)]->coordinates();
      const auto f = nodes_[indices(5)]->coordinates();
      const auto g = nodes_[indices(1)]->coordinates();
      const auto h = nodes_[indices(0)]->coordinates();

      volume_ =
          (1.0 / 12) *
              (a - g).dot(((b - d).cross(c - a)) + ((e - b).cross(f - a)) +
                          ((d - e).cross(h - a))) +
          (1.0 / 12) *
              (b - g).dot(((b - d).cross(c - a)) + ((c - g).cross(c - f))) +
          (1.0 / 12) *
              (e - g).dot(((e - b).cross(f - a)) + ((f - g).cross(h - f))) +
          (1.0 / 12) *
              (d - g).dot(((d - e).cross(h - a)) + ((h - g).cross(h - c)));
    } else {
      throw std::runtime_error(
          "Unable to compute volume, number of vertices is incorrect");
    }
    // Check negative volume
    if (this->volume_ <= 0)
      throw std::runtime_error(
          "Negative or zero volume cell, misconfigured cell!");
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
}

//! Compute the centroid of the cell
template <unsigned Tdim>
void mpm::Cell<Tdim>::compute_centroid() {
  // Get indices of corner nodes
  Eigen::VectorXi indices = element_->corner_indices();

  // Calculate the centroid of the cell
  centroid_.setZero();
  for (unsigned i = 0; i < indices.size(); ++i)
    centroid_ += nodes_[indices(i)]->coordinates();

  centroid_ /= indices.size();
}

//! Compute mean length of cell
template <unsigned Tdim>
void mpm::Cell<Tdim>::compute_mean_length() {
  // Get the indices of sub-triangles
  Eigen::MatrixXi indices = element_->sides_indices();
  this->mean_length_ = 0.;
  // Calculate the mean length
  for (unsigned i = 0; i < indices.rows(); ++i)
    this->mean_length_ += (nodes_[indices(i, 0)]->coordinates() -
                           nodes_[indices(i, 1)]->coordinates())
                              .norm();
  this->mean_length_ /= indices.rows();
}

//! Return nodal coordinates
template <unsigned Tdim>
void mpm::Cell<Tdim>::compute_nodal_coordinates() {
  nodal_coordinates_.resize(this->nnodes_, Tdim);
  try {
    // If cell is initialised
    if (this->is_initialised()) {
      for (unsigned i = 0; i < nodes_.size(); ++i)
        nodal_coordinates_.row(i) = nodes_[i]->coordinates().transpose();
    } else {
      throw std::runtime_error(
          "Cell is not initialised to return nodal coordinates!");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
}

//! Check if a point is in a 1D cell by checking bounding box range
template <>
inline bool mpm::Cell<1>::point_in_cartesian_cell(
    const Eigen::Matrix<double, 1, 1>& point) {
  if (nodes_[0]->coordinates()(0) <= point(0) &&
      nodes_[1]->coordinates()(0) >= point(0))
    return true;
  else
    return false;
}

//! Check if a point is in a 2D cell by checking bounding box range
template <>
inline bool mpm::Cell<2>::point_in_cartesian_cell(
    const Eigen::Matrix<double, 2, 1>& point) {
  // Check if a point is within x and y ranges
  if (nodes_[0]->coordinates()(0) <= point(0) &&
      nodes_[1]->coordinates()(0) >= point(0) &&
      nodes_[1]->coordinates()(1) <= point(1) &&
      nodes_[2]->coordinates()(1) >= point(1))
    return true;
  else
    return false;
}

//! Check if a point is in a 3D cell by checking bounding box range
template <>
inline bool mpm::Cell<3>::point_in_cartesian_cell(
    const Eigen::Matrix<double, 3, 1>& point) {
  // Check if a point is within x, y and z ranges
  if (nodes_[0]->coordinates()(0) <= point(0) &&
      nodes_[1]->coordinates()(0) >= point(0) &&
      nodes_[1]->coordinates()(1) <= point(1) &&
      nodes_[2]->coordinates()(1) >= point(1) &&
      nodes_[0]->coordinates()(2) <= point(2) &&
      nodes_[4]->coordinates()(2) >= point(2))
    return true;
  else
    return false;
}

//! Check approximately if a point is in a cell by using cell length
template <unsigned Tdim>
inline bool mpm::Cell<Tdim>::approx_point_in_cell(
    const Eigen::Matrix<double, Tdim, 1>& point) {
  const double length = (point - this->centroid_).norm();
  if (length < (this->mean_length_ * 2.))
    return true;
  else
    return false;
}

//! Check if a point is in a cell by affine transformation and newton-raphson
template <unsigned Tdim>
inline bool mpm::Cell<Tdim>::is_point_in_cell(
    const Eigen::Matrix<double, Tdim, 1>& point,
    Eigen::Matrix<double, Tdim, 1>* xi) {

  // Set an initial value of Xi
  (*xi).fill(std::numeric_limits<double>::max());

  // Check if point is approximately in the cell
  if (!this->approx_point_in_cell(point)) return false;

  bool status = true;

  // Check if cell is cartesian, if so use cartesian local coordinates
  if (!isoparametric_) (*xi) = this->local_coordinates_point(point);
  // Isoparametric element
  else
    (*xi) = this->transform_real_to_unit_cell(point);

  // Check if the transformed coordinate is within the unit cell:
  // between 0 and 1-xi(1-i) if the element is a triangle, and between
  // -1 and 1 if otherwise
  if (this->element_->corner_indices().size() == 3) {
    for (unsigned i = 0; i < (*xi).size(); ++i)
      if ((*xi)(i) < 0. || (*xi)(i) > 1. - (*xi)(1 - i) || std::isnan((*xi)(i)))
        status = false;
  } else {
    for (unsigned i = 0; i < (*xi).size(); ++i)
      if ((*xi)(i) < -1. || (*xi)(i) > 1. || std::isnan((*xi)(i)))
        status = false;
  }
  return status;
}

//! Return the local coordinates of a point in a 1D cell
template <>
inline Eigen::Matrix<double, 1, 1> mpm::Cell<1>::local_coordinates_point(
    const Eigen::Matrix<double, 1, 1>& point) {
  // Local point coordinates
  Eigen::Matrix<double, 1, 1> xi;

  try {
    // Indices of corner nodes
    Eigen::VectorXi indices = element_->corner_indices();

    // Linear
    if (indices.size() == 2) {
      //
      // 0 0---------0 1
      //        l
      const double length = (nodes_[indices(0)]->coordinates() -
                             nodes_[indices(1)]->coordinates())
                                .norm();

      const Eigen::Matrix<double, 1, 1> centre =
          (nodes_[indices(0)]->coordinates() +
           nodes_[indices(1)]->coordinates()) /
          2.0;

      xi(0) = 2. * (point(0) - centre(0)) / length;
    } else {
      throw std::runtime_error("Unable to compute local coordinates");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return xi;
}

//! Return the local coordinates of a point in a 2D cell
template <>
inline Eigen::Matrix<double, 2, 1> mpm::Cell<2>::local_coordinates_point(
    const Eigen::Matrix<double, 2, 1>& point) {
  // Local point coordinates
  Eigen::Matrix<double, 2, 1> xi;

  try {
    // Indices of corner nodes
    Eigen::VectorXi indices = element_->corner_indices();

    // Triangle
    if (indices.size() == 3) {
      //   2 0
      //     |\
      //     | \  
      //   c |  \ b
      //     |   \
      //     |    \
      //   0 0-----0 1
      //        a
      //

      auto node0 = nodes_[indices(0)]->coordinates();
      auto node1 = nodes_[indices(1)]->coordinates();
      auto node2 = nodes_[indices(2)]->coordinates();

      const double area = ((node1(0) - node0(0)) * (node2(1) - node0(1)) -
                           (node2(0) - node0(0)) * (node1(1) - node0(1))) /
                          2.0;

      xi(0) = 1 / (2 * area) *
              ((point(0) - node0(0)) * (node2(1) - node0(1)) -
               (node2(0) - node0(0)) * (point(1) - node0(1)));

      xi(1) = -1 / (2 * area) *
              ((point(0) - node0(0)) * (node1(1) - node0(1)) -
               (node1(0) - node0(0)) * (point(1) - node0(1)));
      // Quadrilateral
    } else if (indices.size() == 4) {
      //        b
      // 3 0--------0 2
      //   | \   / |
      // a |  \ /  | c
      //   |  / \  |
      //   | /   \ |
      // 0 0---------0 1
      //         d
      const double xlength = (nodes_[indices(0)]->coordinates() -
                              nodes_[indices(1)]->coordinates())
                                 .norm();
      const double ylength = (nodes_[indices(1)]->coordinates() -
                              nodes_[indices(2)]->coordinates())
                                 .norm();

      const Eigen::Matrix<double, 2, 1> centre =
          (nodes_[indices(0)]->coordinates() +
           nodes_[indices(1)]->coordinates() +
           nodes_[indices(2)]->coordinates() +
           nodes_[indices(3)]->coordinates()) /
          4.0;

      xi(0) = 2. * (point(0) - centre(0)) / xlength;
      xi(1) = 2. * (point(1) - centre(1)) / ylength;
    } else {
      throw std::runtime_error("Unable to compute local coordinates");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return xi;
}

//! Return the local coordinates of a point in a 3D cell
template <>
inline Eigen::Matrix<double, 3, 1> mpm::Cell<3>::local_coordinates_point(
    const Eigen::Matrix<double, 3, 1>& point) {

  // Local point coordinates
  Eigen::Matrix<double, 3, 1> xi;

  try {
    // Indices of corner nodes
    Eigen::VectorXi indices = element_->corner_indices();

    // Hexahedron
    if (indices.size() == 8) {
      // Node numbering as read in by mesh file
      //        3               2
      //          *_ _ _ _ _ _*
      //         /|           /|
      //        / |          / |
      //     7 *_ |_ _ _ _ _* 6|
      //       |  |         |  |
      //       |  |         |  |
      //       |  *_ _ _ _ _|_ *
      //       | / 0        | / 1
      //       |/           |/
      //       *_ _ _ _ _ _ *
      //     4               5
      //

      const double xlength = (nodes_[indices(0)]->coordinates() -
                              nodes_[indices(1)]->coordinates())
                                 .norm();
      const double ylength = (nodes_[indices(1)]->coordinates() -
                              nodes_[indices(2)]->coordinates())
                                 .norm();
      const double zlength = (nodes_[indices(1)]->coordinates() -
                              nodes_[indices(5)]->coordinates())
                                 .norm();

      // Compute centre
      Eigen::Matrix<double, 3, 1> centre;
      centre.setZero();
      for (unsigned i = 0; i < indices.size(); ++i)
        centre += nodes_[indices(i)]->coordinates();
      centre /= indices.size();

      xi(0) = 2. * (point(0) - centre(0)) / xlength;
      xi(1) = 2. * (point(1) - centre(1)) / ylength;
      xi(2) = 2. * (point(2) - centre(2)) / zlength;

    } else {
      throw std::runtime_error("Unable to compute local coordinates");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return xi;
}

//! Return the local coordinates of a point in a 2D cell
template <>
inline Eigen::Matrix<double, 1, 1> mpm::Cell<1>::local_coordinates_point_2d(
    const Eigen::Matrix<double, 1, 1>& point) {
  // Local point coordinates
  Eigen::Matrix<double, 1, 1> xi;
  xi << std::numeric_limits<double>::max();
  console_->error("{} #{}: The analytical solution is valid only for 2D\n",
                  __FILE__, __LINE__);
  return xi;
}

//! Return the local coordinates of a point in a 2D cell
//! Analytical solution based on A consistent point-searching algorithm for
//! solution interpolation in unstructured meshes consisting of 4-node bilinear
//! quadrilateral elements - Zhao et al., 1999

template <>
inline Eigen::Matrix<double, 2, 1> mpm::Cell<2>::local_coordinates_point_2d(
    const Eigen::Matrix<double, 2, 1>& point) {
  // Local point coordinates
  Eigen::Matrix<double, 2, 1> xi;
  xi.fill(std::numeric_limits<double>::max());

  const double xa = point(0);
  const double ya = point(1);
  const double x1 = nodes_[0]->coordinates()(0);
  const double y1 = nodes_[0]->coordinates()(1);
  const double x2 = nodes_[1]->coordinates()(0);
  const double y2 = nodes_[1]->coordinates()(1);
  const double x3 = nodes_[2]->coordinates()(0);
  const double y3 = nodes_[2]->coordinates()(1);
  const double x4 = nodes_[3]->coordinates()(0);
  const double y4 = nodes_[3]->coordinates()(1);

  const double a1 = x1 + x2 + x3 + x4;
  const double a2 = -x1 + x2 + x3 - x4;
  const double a3 = -x1 - x2 + x3 + x4;
  const double a4 = x1 - x2 + x3 - x4;

  const double b1 = y1 + y2 + y3 + y4;
  const double b2 = -y1 + y2 + y3 - y4;
  const double b3 = -y1 - y2 + y3 + y4;
  const double b4 = y1 - y2 + y3 - y4;

  const double c1 = 4. * xa - a1;
  const double c2 = 4. * ya - b1;

  // General solution of xi and eta based on solving with Sympy
  // a2 * xi(0) + a3 * xi(1) + a4 * xi(0) * xi(1) = 4 x_a - a1
  // b2 * xi(0) + b3 * xi(1) + b4 * xi(0) * xi(1) = 4 y_a - b1
  const double u1 =
      (-a1 * b4 + a4 * b1 - 4 * a4 * ya + 4 * b4 * xa +
       (-a3 * b4 + a4 * b3) *
           ((-a1 * b4 + a2 * b3 - a3 * b2 + a4 * b1 - 4 * a4 * ya +
             4 * b4 * xa) /
                (2 * (a3 * b4 - a4 * b3)) -
            (std::sqrt(a1 * a1 * b4 * b4 - 2 * a1 * a2 * b3 * b4 -
                       2 * a1 * a3 * b2 * b4 - 2 * a1 * a4 * b1 * b4 +
                       4 * a1 * a4 * b2 * b3 + 8 * a1 * a4 * b4 * ya -
                       8 * a1 * b4 * b4 * xa + a2 * a2 * b3 * b3 +
                       4 * a2 * a3 * b1 * b4 - 2 * a2 * a3 * b2 * b3 -
                       16 * a2 * a3 * b4 * ya - 2 * a2 * a4 * b1 * b3 +
                       8 * a2 * a4 * b3 * ya + 8 * a2 * b3 * b4 * xa +
                       a3 * a3 * b2 * b2 - 2 * a3 * a4 * b1 * b2 +
                       8 * a3 * a4 * b2 * ya + 8 * a3 * b2 * b4 * xa +
                       a4 * a4 * b1 * b1 - 8 * a4 * a4 * b1 * ya +
                       16 * a4 * a4 * ya * ya + 8 * a4 * b1 * b4 * xa -
                       16 * a4 * b2 * b3 * xa - 32 * a4 * b4 * xa * ya +
                       16 * b4 * b4 * xa * xa)) /
                (2 * (a3 * b4 - a4 * b3)))) /
      (a2 * b4 - a4 * b2);

  const double u2 =
      (-a1 * b4 + a2 * b3 - a3 * b2 + a4 * b1 - 4 * a4 * ya + 4 * b4 * xa) /
          (2 * (a3 * b4 - a4 * b3)) -
      (std::sqrt(
          a1 * a1 * b4 * b4 - 2 * a1 * a2 * b3 * b4 - 2 * a1 * a3 * b2 * b4 -
          2 * a1 * a4 * b1 * b4 + 4 * a1 * a4 * b2 * b3 +
          8 * a1 * a4 * b4 * ya - 8 * a1 * b4 * b4 * xa + a2 * a2 * b3 * b3 +
          4 * a2 * a3 * b1 * b4 - 2 * a2 * a3 * b2 * b3 -
          16 * a2 * a3 * b4 * ya - 2 * a2 * a4 * b1 * b3 +
          8 * a2 * a4 * b3 * ya + 8 * a2 * b3 * b4 * xa + a3 * a3 * b2 * b2 -
          2 * a3 * a4 * b1 * b2 + 8 * a3 * a4 * b2 * ya +
          8 * a3 * b2 * b4 * xa + a4 * a4 * b1 * b1 - 8 * a4 * a4 * b1 * ya +
          16 * a4 * a4 * ya * ya + 8 * a4 * b1 * b4 * xa -
          16 * a4 * b2 * b3 * xa - 32 * a4 * b4 * xa * ya +
          16 * b4 * b4 * xa * xa)) /
          (2 * (a3 * b4 - a4 * b3));

  // Second solution of a quadratic equation
  const double v1 =
      (-a1 * b4 + a4 * b1 - 4 * a4 * ya + 4 * b4 * xa +
       (-a3 * b4 + a4 * b3) *
           ((-a1 * b4 + a2 * b3 - a3 * b2 + a4 * b1 - 4 * a4 * ya +
             4 * b4 * xa) /
                (2 * (a3 * b4 - a4 * b3)) +
            (std::sqrt(a1 * a1 * b4 * b4 - 2 * a1 * a2 * b3 * b4 -
                       2 * a1 * a3 * b2 * b4 - 2 * a1 * a4 * b1 * b4 +
                       4 * a1 * a4 * b2 * b3 + 8 * a1 * a4 * b4 * ya -
                       8 * a1 * b4 * b4 * xa + a2 * a2 * b3 * b3 +
                       4 * a2 * a3 * b1 * b4 - 2 * a2 * a3 * b2 * b3 -
                       16 * a2 * a3 * b4 * ya - 2 * a2 * a4 * b1 * b3 +
                       8 * a2 * a4 * b3 * ya + 8 * a2 * b3 * b4 * xa +
                       a3 * a3 * b2 * b2 - 2 * a3 * a4 * b1 * b2 +
                       8 * a3 * a4 * b2 * ya + 8 * a3 * b2 * b4 * xa +
                       a4 * a4 * b1 * b1 - 8 * a4 * a4 * b1 * ya +
                       16 * a4 * a4 * ya * ya + 8 * a4 * b1 * b4 * xa -
                       16 * a4 * b2 * b3 * xa - 32 * a4 * b4 * xa * ya +
                       16 * b4 * b4 * xa * xa)) /
                (2 * (a3 * b4 - a4 * b3)))) /
      (a2 * b4 - a4 * b2);

  const double v2 =
      (-a1 * b4 + a2 * b3 - a3 * b2 + a4 * b1 - 4 * a4 * ya + 4 * b4 * xa) /
          (2 * (a3 * b4 - a4 * b3)) +
      (std::sqrt(
          a1 * a1 * b4 * b4 - 2 * a1 * a2 * b3 * b4 - 2 * a1 * a3 * b2 * b4 -
          2 * a1 * a4 * b1 * b4 + 4 * a1 * a4 * b2 * b3 +
          8 * a1 * a4 * b4 * ya - 8 * a1 * b4 * b4 * xa + a2 * a2 * b3 * b3 +
          4 * a2 * a3 * b1 * b4 - 2 * a2 * a3 * b2 * b3 -
          16 * a2 * a3 * b4 * ya - 2 * a2 * a4 * b1 * b3 +
          8 * a2 * a4 * b3 * ya + 8 * a2 * b3 * b4 * xa + a3 * a3 * b2 * b2 -
          2 * a3 * a4 * b1 * b2 + 8 * a3 * a4 * b2 * ya +
          8 * a3 * b2 * b4 * xa + a4 * a4 * b1 * b1 - 8 * a4 * a4 * b1 * ya +
          16 * a4 * a4 * ya * ya + 8 * a4 * b1 * b4 * xa -
          16 * a4 * b2 * b3 * xa - 32 * a4 * b4 * xa * ya +
          16 * b4 * b4 * xa * xa)) /
          (2 * (a3 * b4 - a4 * b3));

  // Choosing a quadratic solution
  if (u1 >= -1. && u1 <= 1. && u2 >= -1. && u2 <= 1.) {
    xi(0) = u1;
    xi(1) = u2;
    return xi;
  } else if (v1 >= -1. && v1 <= 1. && v2 >= -1. && v2 <= 1.) {
    xi(0) = v1;
    xi(1) = v2;
    return xi;
  }

  // Case1: a4 == 0 and b4 != 0: Eq 10
  if (a4 == 0 && b4 == 0) {
    xi(0) = (b3 * c1 - a3 * c2) / (a2 * b3 - a3 * b2);
    xi(1) = (-b2 * c1 + a2 * c2) / (a2 * b3 - a3 * b2);
  } else if (a4 == 0 && b4 != 0) {  // Case 2: Eq 11
    if (a2 == 0 && a3 != 0) {       // Case 2.1 Eq 12
      xi(1) = c1 / a3;
      xi(0) = (c2 - b3 * xi(1)) / (b2 + b4 * xi(1));
    } else if (a2 != 0 && a3 == 0) {  // Case 2.2 Eq 13
      xi(0) = c1 / a2;
      xi(1) = (c2 - b2 * xi(0)) / (b3 + b4 * xi(0));
    } else {  // a2 != 0 && a3 != 0 // Case 2.3 Eq 14
      const double aa = b4 * a3 / a2;
      const double bb = ((b2 * a3 - b4 * c1) / a2) - b3;
      const double cc = -(b2 * c1 / a2) + c2;
      // There are two possible solutions
      const double u2 = (-bb + std::sqrt(bb * bb - 4 * aa * cc)) / (2 * aa);
      const double u1 = (c1 - a3 * u2) / a2;
      // Second solution of a quadratic equation
      const double v2 = (-bb - std::sqrt(bb * bb - 4 * aa * cc)) / (2 * aa);
      const double v1 = (c1 - a3 * v2) / a2;
      if (u1 >= -1. && u1 <= 1. && u2 >= -1. && u2 <= 1.) {
        xi(0) = u1;
        xi(1) = u2;
      } else if (v1 >= -1. && v1 <= 1. && v2 >= -1. && v2 <= 1.) {
        xi(0) = v1;
        xi(1) = v2;
      }
    }
  } else if (a4 != 0 && b4 == 0) {  // Case 3 Eq 16
    if (b2 == 0 && b3 != 0) {       // Case 3.1 Eq 17
      xi(1) = c2 / b3;
      xi(0) = (c1 - a3 * xi(1)) / (a2 + a4 * xi(1));
    } else if (b2 != 0 && b3 == 0) {  // Case 3.2 Eq 18
      xi(0) = c2 / b2;
      xi(1) = (c1 - a2 * xi(0)) / (a3 + a4 * xi(0));
    } else {  // b2 != 0 && b3 != 0  // Case 3.3 Eq 19
      const double aa = a4 * b3 / b2;
      const double bb = ((a2 * b3 - a4 * c2) / b2) - b3;
      const double cc = -(a2 * c2 / b2) + c1;
      // There are two possible solutions
      const double u2 = (-bb + std::sqrt((bb * bb - 4 * aa * cc))) / (2 * aa);
      const double u1 = (c2 - b3 * u2) / b2;
      // Second solution of a quadratic equation
      const double v2 = (-bb - std::sqrt((bb * bb - 4 * aa * cc))) / (2 * aa);
      const double v1 = (c2 - b3 * v2) / b2;
      if (u1 >= -1. && u1 <= 1. && u2 >= -1. && u2 <= 1.) {
        xi(0) = u1;
        xi(1) = u2;
      } else if (v1 >= -1. && v1 <= 1. && v2 >= -1. && v2 <= 1.) {
        xi(0) = v1;
        xi(1) = v2;
      }
    }
  } else {  // a4 != 0 && b4 != 0   // Case 4 Eq 21
    const double a2s = a2 / a4;
    const double a3s = a3 / a4;

    const double b2s = b2 / b4;
    const double b3s = b3 / b4;

    const double c1s = c1 / a4;
    const double c2s = c2 / b4;

    if ((a2s - b2s) == 0) {  // Case 4.1 Eq 25
      xi(1) = (c1s - c2s) / (a3s - b3s);
      xi(0) = (c1s - a3s * xi(1)) / (a2s + xi(1));
    } else {
      const double alpha = (c1s - c2s) / (a2s - b2s);
      const double beta = (a3s - b3s) / (a2s - b2s);

      if (beta == 0) {  // Case 4.2a Eq 28
        xi(0) = alpha;
        xi(1) = (c1s - a2s * xi(0)) / (a3s + xi(0));
      } else {  // Case 4.2b Eq 29
        // There are two possible solutions
        const double u2 =
            (-(a2s * beta + a3s - alpha) +
             std::sqrt((a2s * beta + a3s - alpha) * (a2s * beta + a3s - alpha) -
                       (4 * beta * (c1s - a2s * alpha)))) /
            (2. * beta);
        const double u1 = alpha - beta * u2;
        // Second solution of a quadratic equation
        const double v2 =
            (-(a2s * beta + a3s - alpha) -
             std::sqrt((a2s * beta + a3s - alpha) * (a2s * beta + a3s - alpha) -
                       (4 * beta * (c1s - a2s * alpha)))) /
            (2. * beta);
        const double v1 = alpha - beta * v2;
        if (u1 >= -1. && u1 <= 1. && u2 >= -1. && u2 <= 1.) {
          xi(0) = u1;
          xi(1) = u2;
        } else if (v1 >= -1. && v1 <= 1. && v2 >= -1. && v2 <= 1.) {
          xi(0) = v1;
          xi(1) = v2;
        }
      }
    }
  }
  return xi;
}

//! Return the local coordinates of a point in a 2D cell
template <>
inline Eigen::Matrix<double, 3, 1> mpm::Cell<3>::local_coordinates_point_2d(
    const Eigen::Matrix<double, 3, 1>& point) {
  // Local point coordinates
  Eigen::Matrix<double, 3, 1> xi;
  xi.fill(std::numeric_limits<double>::max());
  console_->error("{} #{}: The analytical solution is valid only for 2D\n",
                  __FILE__, __LINE__);
  return xi;
}

//! Return the local coordinates of a point in a 1D cell
template <>
inline Eigen::Matrix<double, 1, 1> mpm::Cell<1>::transform_real_to_unit_cell(
    const Eigen::Matrix<double, 1, 1>& point) {
  return this->local_coordinates_point(point);
}

//! Return the local coordinates of a point in a 2D cell
template <>
inline Eigen::Matrix<double, 2, 1> mpm::Cell<2>::transform_real_to_unit_cell(
    const Eigen::Matrix<double, 2, 1>& point) {

  // If not isoparametric then use cartesian transformation
  if (!this->isoparametric_)
    return mpm::Cell<2>::local_coordinates_point(point);

  // Local coordinates of a point in an unit cell
  Eigen::Matrix<double, 2, 1> xi;
  xi.setZero();

  // Zeros
  const Eigen::Matrix<double, 2, 1> zero = Eigen::Matrix<double, 2, 1>::Zero();

  // Get indices of corner nodes
  Eigen::VectorXi indices = element_->corner_indices();

  // Matrix of nodal coordinates
  Eigen::MatrixXd nodal_coords;
  nodal_coords.resize(2, indices.size());
  for (unsigned j = 0; j < indices.size(); ++j) {
    Eigen::Matrix<double, 2, 1> node = nodes_[indices(j)]->coordinates();
    for (unsigned i = 0; i < 2; ++i) {
      nodal_coords(i, j) = node[i];
    }
  }

  // Try analytical solution
  xi = this->local_coordinates_point_2d(point);

  // Analytical tolerance
  const double tolerance = 1.0E-16 * mean_length_ * mean_length_;

  bool status = true;
  // Check if xi is within the cell
  for (unsigned i = 0; i < xi.size(); ++i)
    if (xi(i) < -1. || xi(i) > 1. || std::isnan(xi(i))) status = false;
  // Local shape function
  const auto sf = element_->shapefn_local(xi, zero, zero);
  // f(x) = p(x) - p, where p is the real point
  const auto analytical_residual = (nodal_coords * sf) - point;
  // Early exit
  if ((analytical_residual.squaredNorm() < tolerance) && status) return xi;
  // If tolerance is high
  const auto analytical_xi = xi;

  // Affine guess of xi
  Eigen::Matrix<double, 2, 1> affine_guess;
  // Boolean to check if affine is nan
  bool affine_nan = false;

  // Affine tolerance
  const double affine_tolerance = 1.0E-16 * mean_length_ * mean_length_;

  // Coordinates of a unit cell
  const auto unit_cell = element_->unit_cell_coordinates();

  // Affine residual
  Eigen::Matrix<double, 2, 1> affine_residual;

  // Affine transformation, using linear interpolation for the initial guess
  if (element_->degree() == mpm::ElementDegree::Linear) {
    // A = vertex * KA
    const Eigen::Matrix<double, 2, 2> A =
        nodal_coords * mpm::TransformR2UAffine<2, 4>::KA;

    // b = vertex * Kb
    const Eigen::Matrix<double, 2, 1> b =
        point - (nodal_coords * mpm::TransformR2UAffine<2, 4>::Kb);

    // Affine transform: A^-1 * b
    // const Eigen::Matrix<double, 2, 1>
    affine_guess = A.inverse() * b;

    // Check for nan
    for (unsigned i = 0; i < affine_guess.size(); ++i)
      if (std::isnan(affine_guess(i))) affine_nan = true;

    // Set xi to affine guess
    if (!affine_nan) xi = affine_guess;
    // If guess is nan set xi to zero
    else
      xi.setZero();

    // Local shape function
    const auto sf = element_->shapefn_local(xi, zero, zero);

    // f(x) = p(x) - p, where p is the real point
    affine_residual = (nodal_coords * sf) - point;

    // Early exit
    if ((affine_residual.squaredNorm() < affine_tolerance) && !affine_nan)
      return xi;
  }
  // Select the best of analytical and affine transformation
  Eigen::Matrix<double, 2, 1> geometry_xi;
  Eigen::Matrix<double, 2, 1> geometry_residual;
  geometry_residual.fill(std::numeric_limits<double>::max());

  if (analytical_residual.norm() < affine_residual.norm()) {
    geometry_residual = analytical_residual;
    geometry_xi = analytical_xi;
  } else {
    geometry_residual = affine_residual;
    geometry_xi = affine_guess;
  }

  // Maximum iterations of newton raphson
  const unsigned max_iterations = 100;

  // Tolerance for newton raphson
  const double Tolerance = 1.0E-10;

  // Newton Raphson iteration to solve for x
  // x_{n+1} = x_n - f(x)/f'(x)
  // f(x) = p(x) - p, where p is the real point
  // p(x) is the computed point.
  Eigen::Matrix<double, 2, 1> nr_residual;
  unsigned iter = 0;
  for (; iter < max_iterations; ++iter) {
    // Calculate local Jacobian
    const Eigen::Matrix<double, 2, 2> jacobian =
        element_->jacobian_local(xi, unit_cell, zero, zero);

    // Set guess xi to zero
    if (std::abs(jacobian.determinant()) < 1.0E-10) xi.setZero();

    // Local shape function
    const auto sf = element_->shapefn_local(xi, zero, zero);

    // Residual (f(x))
    // f(x) = p(x) - p, where p is the real point
    nr_residual = (nodal_coords * sf) - point;

    // f(x)/f'(x)
    const Eigen::Matrix<double, 2, 1> delta = jacobian.inverse() * nr_residual;

    // Line search
    double step_length = 1.;
    for (unsigned line_trials = 0; line_trials < 10; ++line_trials) {
      // Trial xi
      // x_{n+1} = x_n - f(x)/f'(x)
      const Eigen::Matrix<double, 2, 1> xi_trial = xi - (step_length * delta);

      // Trial shape function
      const auto sf_trial = element_->shapefn_local(xi_trial, zero, zero);

      // Trial residual: f(x) = p(x) - p, where p is the real point
      const Eigen::Matrix<double, 2, 1> nr_residual_trial =
          (nodal_coords * sf_trial) - point;

      if (nr_residual_trial.norm() < nr_residual.norm()) {
        xi = xi_trial;
        nr_residual = nr_residual_trial;
        break;
      } else if (step_length > 0.05)
        step_length /= 2.;
      else {
        // Line search failed
        break;
      }
    }
    // Convergence criteria
    if ((step_length * delta).norm() < Tolerance) break;

    // Check for nan and set to a trial xi
    if (std::isnan(xi(0)) || std::isnan(xi(1))) xi.setZero();
  }

  // At end of iteration return affine or xi based on lowest norm
  if ((iter == max_iterations) && !affine_nan &&
      (element_->degree() == mpm::ElementDegree::Linear))
    return geometry_residual.norm() < nr_residual.norm() ? geometry_xi : xi;

  return xi;
}

//! Return the local coordinates of a point in a 2D/3D cell
template <>
inline Eigen::Matrix<double, 3, 1> mpm::Cell<3>::transform_real_to_unit_cell(
    const Eigen::Matrix<double, 3, 1>& point) {

  // If not isoparametric then use cartesian transformation
  if (!this->isoparametric_)
    return mpm::Cell<3>::local_coordinates_point(point);

  // Local coordinates of a point in an unit cell
  Eigen::Matrix<double, 3, 1> xi;
  xi.setZero();

  // Get indices of corner nodes
  Eigen::VectorXi indices = element_->corner_indices();

  // Matrix of nodal coordinates
  Eigen::MatrixXd nodal_coords;
  nodal_coords.resize(3, indices.size());

  for (unsigned j = 0; j < indices.size(); ++j) {
    Eigen::Matrix<double, 3, 1> node = nodes_[indices(j)]->coordinates();
    for (unsigned i = 0; i < 3; ++i) {
      nodal_coords(i, j) = node[i];
    }
  }

  // Affine transformation, using linear interpolation for the initial guess
  // Affine guess of xi
  Eigen::Matrix<double, 3, 1> affine_guess;
  // Boolean to check if affine is nan
  bool affine_nan = false;
  // Zeros
  const Eigen::Matrix<double, 3, 1> zero = Eigen::Matrix<double, 3, 1>::Zero();

  // Affine tolerance
  const double affine_tolerance = 1.0E-16 * mean_length_ * mean_length_;

  // Coordinates of a unit cell
  const auto unit_cell = element_->unit_cell_coordinates();

  // Affine residual
  Eigen::Matrix<double, 3, 1> affine_residual;

  // Affine transformation, using linear interpolation for the initial guess
  if (element_->degree() == mpm::ElementDegree::Linear) {
    // A = vertex * KA
    const Eigen::Matrix<double, 3, 3> A =
        nodal_coords * mpm::TransformR2UAffine<3, 8>::KA;

    // b = vertex * Kb
    const Eigen::Matrix<double, 3, 1> b =
        point - (nodal_coords * mpm::TransformR2UAffine<3, 8>::Kb);

    // Affine transform: A^-1 * b
    // const Eigen::Matrix<double, 3, 1>
    affine_guess = A.inverse() * b;

    // Check for nan
    for (unsigned i = 0; i < affine_guess.size(); ++i)
      if (std::isnan(affine_guess(i))) affine_nan = true;

    // Set xi to affine guess
    if (!affine_nan) xi = affine_guess;
    // If guess is nan set xi to zero
    else
      xi.setZero();

    // Local shape function
    const auto sf = element_->shapefn_local(xi, zero, zero);

    // f(x) = p(x) - p, where p is the real point
    affine_residual = (nodal_coords * sf) - point;

    // Early exit
    if ((affine_residual.squaredNorm() < affine_tolerance) && !affine_nan)
      return xi;
  }

  // Maximum iterations of newton raphson
  const unsigned max_iterations = 100;

  // Tolerance for newton raphson
  const double Tolerance = 1.0E-10;

  // Newton Raphson iteration to solve for x
  // x_{n+1} = x_n - f(x)/f'(x)
  // f(x) = p(x) - p, where p is the real point
  // p(x) is the computed point.
  Eigen::Matrix<double, 3, 1> nr_residual;
  unsigned iter = 0;
  for (; iter < max_iterations; ++iter) {
    // Calculate local Jacobian
    const Eigen::Matrix<double, 3, 3> jacobian =
        element_->jacobian_local(xi, unit_cell, zero, zero);

    // Set guess xi to zero
    if (std::abs(jacobian.determinant()) < 1.0E-10) xi.setZero();

    // Local shape function
    const auto sf = element_->shapefn_local(xi, zero, zero);

    // Residual (f(x))
    // f(x) = p(x) - p, where p is the real point
    nr_residual = (nodal_coords * sf) - point;

    // f(x)/f'(x)
    const Eigen::Matrix<double, 3, 1> delta = jacobian.inverse() * nr_residual;

    // Line search
    double step_length = 1.;
    for (unsigned line_trials = 0; line_trials < 10; ++line_trials) {
      // Trial xi
      // x_{n+1} = x_n - f(x)/f'(x)
      const Eigen::Matrix<double, 3, 1> xi_trial = xi - (step_length * delta);

      // Trial shape function
      const auto sf_trial = element_->shapefn_local(xi_trial, zero, zero);

      // Trial residual: f(x) = p(x) - p, where p is the real point
      const Eigen::Matrix<double, 3, 1> nr_residual_trial =
          (nodal_coords * sf_trial) - point;

      if (nr_residual_trial.norm() < nr_residual.norm()) {
        xi = xi_trial;
        nr_residual = nr_residual_trial;
        break;
      } else if (step_length > 0.05)
        step_length /= 2.;
      else {
        // Line search failed
        break;
      }
    }
    // Convergence criteria
    if ((step_length * delta).norm() < Tolerance) break;

    // Check for nan and set to a trial xi
    if (std::isnan(xi(0)) || std::isnan(xi(1))) xi.setZero();
  }

  // At end of iteration return affine or xi based on lowest norm
  if ((iter == max_iterations) && !affine_nan &&
      (element_->degree() == mpm::ElementDegree::Linear))
    return affine_residual.norm() < nr_residual.norm() ? affine_guess : xi;

  return xi;
}

//! Map particle mass to nodes
template <unsigned Tdim>
void mpm::Cell<Tdim>::map_particle_mass_to_nodes(const Eigen::VectorXd& shapefn,
                                                 unsigned phase, double pmass) {
  for (unsigned i = 0; i < shapefn.size(); ++i) {
    nodes_[i]->update_mass(true, phase, shapefn(i) * pmass);
  }
}

//! Map particle volume to nodes
template <unsigned Tdim>
void mpm::Cell<Tdim>::map_particle_volume_to_nodes(
    const Eigen::VectorXd& shapefn, unsigned phase, double pvolume) {

  for (unsigned i = 0; i < shapefn.size(); ++i) {
    nodes_[i]->update_volume(true, phase, shapefn(i) * pvolume);
  }
}

//! Map particle mass and momentum to nodes for a given phase
template <unsigned Tdim>
void mpm::Cell<Tdim>::map_mass_momentum_to_nodes(
    const Eigen::VectorXd& shapefn, unsigned phase, double pmass,
    const Eigen::VectorXd& pvelocity) {
  for (unsigned i = 0; i < this->nfunctions(); ++i) {
    nodes_[i]->update_mass(true, phase, shapefn(i) * pmass);
    nodes_[i]->update_momentum(true, phase, shapefn(i) * pmass * pvelocity);
  }
}

//! Map particle pressure to nodes for a given phase
template <unsigned Tdim>
void mpm::Cell<Tdim>::map_pressure_to_nodes(const Eigen::VectorXd& shapefn,
                                            unsigned phase, double pmass,
                                            double ppressure) {

  for (unsigned i = 0; i < this->nfunctions(); ++i)
    nodes_[i]->update_mass_pressure(phase, shapefn(i) * pmass * ppressure);
}

//! Compute nodal momentum from particle mass and velocity for a given phase
template <unsigned Tdim>
void mpm::Cell<Tdim>::compute_nodal_momentum(const Eigen::VectorXd& shapefn,
                                             unsigned phase, double pmass,
                                             const Eigen::VectorXd& pvelocity) {

  for (unsigned i = 0; i < this->nfunctions(); ++i)
    nodes_[i]->update_momentum(true, phase, shapefn(i) * pmass * pvelocity);
}

//! Compute strain rate
template <unsigned Tdim>
Eigen::VectorXd mpm::Cell<Tdim>::compute_strain_rate(
    const std::vector<Eigen::MatrixXd>& bmatrix, unsigned phase) {
  // Define strain rate
  Eigen::Matrix<double, Tdof, 1> strain_rate =
      Eigen::Matrix<double, Tdof, 1>::Zero(bmatrix.at(0).rows());

  try {
    // Check if B-Matrix size and number of nodes match
    if (this->nfunctions() != bmatrix.size() ||
        this->nnodes() != bmatrix.size())
      throw std::runtime_error(
          "Number of nodes / shapefn doesn't match BMatrix");

    for (unsigned i = 0; i < this->nnodes(); ++i)
      strain_rate += bmatrix.at(i) * nodes_[i]->velocity(phase);
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return strain_rate;
}

//! Compute strain rate for reduced integration at the centroid of cell
template <unsigned Tdim>
Eigen::VectorXd mpm::Cell<Tdim>::compute_strain_rate_centroid(unsigned phase) {

  // Define strain rate at centroid
  Eigen::Matrix<double, Tdof, 1> strain_rate_centroid =
      Eigen::Matrix<double, Tdof, 1>::Zero(bmatrix_centroid_.at(0).rows());

  // Compute strain rate
  for (unsigned i = 0; i < this->nnodes(); ++i)
    strain_rate_centroid +=
        bmatrix_centroid_.at(i) * nodes_[i]->velocity(phase);
  return strain_rate_centroid;
}

//! Compute the nodal body force of a cell from particle mass and gravity
template <unsigned Tdim>
void mpm::Cell<Tdim>::compute_nodal_body_force(const Eigen::VectorXd& shapefn,
                                               unsigned phase, double pmass,
                                               const VectorDim& pgravity) {
  // Map external forces from particle to nodes
  for (unsigned i = 0; i < this->nfunctions(); ++i)
    nodes_[i]->update_external_force(true, phase,
                                     (shapefn(i) * pgravity * pmass));
}

//! Compute the nodal traction force of a cell from particle
template <unsigned Tdim>
void mpm::Cell<Tdim>::compute_nodal_traction_force(
    const Eigen::VectorXd& shapefn, unsigned phase, const VectorDim& traction) {
  // Map external forces from particle to nodes
  for (unsigned i = 0; i < this->nfunctions(); ++i)
    nodes_[i]->update_external_force(true, phase, (shapefn(i) * traction));
}

//! Compute the nodal internal force  of a cell from particle stress and
//! volume
template <unsigned Tdim>
inline void mpm::Cell<Tdim>::compute_nodal_internal_force(
    const std::vector<Eigen::MatrixXd>& bmatrix, unsigned phase, double pvolume,
    const Eigen::Matrix<double, 6, 1>& pstress) {
  // Define strain rate
  Eigen::VectorXd stress;

  switch (Tdim) {
    case (1): {
      stress.resize(1);
      stress.setZero();
      stress(0) = pstress(0);
      break;
    }
    case (2): {
      stress.resize(3);
      stress.setZero();
      stress(0) = pstress(0);
      stress(1) = pstress(1);
      stress(2) = pstress(3);
      break;
    }
    default: {
      stress.resize(6);
      stress = pstress;
      break;
    }
  }
  // Map internal forces from particle to nodes
  for (unsigned j = 0; j < this->nfunctions(); ++j)
    nodes_[j]->update_internal_force(
        true, phase, (pvolume * bmatrix.at(j).transpose() * stress));
}

//! Return velocity at a given point by interpolating from nodes
template <unsigned Tdim>
Eigen::Matrix<double, Tdim, 1> mpm::Cell<Tdim>::interpolate_nodal_velocity(
    const Eigen::VectorXd& shapefn, unsigned phase) {
  Eigen::Matrix<double, Tdim, 1> velocity;
  velocity.setZero();

  for (unsigned i = 0; i < this->nfunctions(); ++i)
    velocity += shapefn(i) * nodes_[i]->velocity(phase);

  return velocity;
}

//! Return acceleration at a point by interpolating from nodes
template <unsigned Tdim>
Eigen::Matrix<double, Tdim, 1> mpm::Cell<Tdim>::interpolate_nodal_acceleration(
    const Eigen::VectorXd& shapefn, unsigned phase) {
  Eigen::Matrix<double, Tdim, 1> acceleration;
  acceleration.setZero();
  for (unsigned i = 0; i < this->nfunctions(); ++i)
    acceleration += shapefn(i) * nodes_[i]->acceleration(phase);

  return acceleration;
}

//! Return pressure at a point by interpolating from nodes
template <unsigned Tdim>
double mpm::Cell<Tdim>::interpolate_nodal_pressure(
    const Eigen::VectorXd& shapefn, unsigned phase) {
  double pressure = 0;
  for (unsigned i = 0; i < this->nfunctions(); ++i)
    pressure += shapefn(i) * nodes_[i]->pressure(phase);

  return pressure;
}

//! Assign velocity constraint
template <unsigned Tdim>
bool mpm::Cell<Tdim>::assign_velocity_constraint(unsigned face_id, unsigned dir,
                                                 double velocity) {
  bool status = true;
  try {
    //! Constraint directions can take values between 0 and Dim * Nphases - 1
    if (dir >= 0 && face_id < element_->nfaces()) {
      this->velocity_constraints_[face_id].emplace_back(
          std::make_pair<unsigned, double>(static_cast<unsigned>(dir),
                                           static_cast<double>(velocity)));
    } else
      throw std::runtime_error("Constraint direction is out of bounds");
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Compute all face normals 2d
template <>
inline void mpm::Cell<2>::compute_normals() {

  //! Set number of faces from element
  for (unsigned face_id = 0; face_id < element_->nfaces(); ++face_id) {
    // Get the nodes of the face
    const Eigen::VectorXi indices = element_->face_indices(face_id);

    // Compute the vector to calculate normal (perpendicular)
    // a = node(0) - node(1)
    Eigen::Matrix<double, 2, 1> a = (this->nodes_[indices(0)])->coordinates() -
                                    (this->nodes_[indices(1)])->coordinates();

    // Compute normal and make unit vector
    // The normal vector n to vector a is defined such that the dot product
    // between a and n is always 0 In 2D, n(0) = -a(1), n(1) = a(0) Note that
    // the reverse does not work to produce normal that is positive pointing out
    // of the element
    Eigen::Matrix<double, 2, 1> normal_vector;
    normal_vector << -a(1), a(0);
    normal_vector = normal_vector.normalized();

    face_normals_.insert(std::make_pair<unsigned, Eigen::VectorXd>(
        static_cast<unsigned>(face_id), normal_vector));
  }
}

//! Compute all face normals 3d
template <>
inline void mpm::Cell<3>::compute_normals() {

  //! Set number of faces from element
  for (unsigned face_id = 0; face_id < element_->nfaces(); ++face_id) {
    // Get the nodes of the face
    const Eigen::VectorXi indices = element_->face_indices(face_id);

    // Compute two vectors to calculate normal
    // a = node(1) - node(0)
    // b = node(3) - node(0)
    Eigen::Matrix<double, 3, 1> a = (this->nodes_[indices(1)])->coordinates() -
                                    (this->nodes_[indices(0)])->coordinates();
    Eigen::Matrix<double, 3, 1> b = (this->nodes_[indices(3)])->coordinates() -
                                    (this->nodes_[indices(0)])->coordinates();

    // Compute normal and make unit vector
    // normal = a x b
    // Note that definition of a and b are such that normal is always out of
    // page
    Eigen::Matrix<double, 3, 1> normal_vector = a.cross(b);
    normal_vector = normal_vector.normalized();

    face_normals_.insert(std::make_pair<unsigned, Eigen::VectorXd>(
        static_cast<unsigned>(face_id), normal_vector));
  }
}

//! Return a sorted list of face node ids
template <unsigned Tdim>
inline std::vector<std::vector<mpm::Index>>
    mpm::Cell<Tdim>::sorted_face_node_ids() {
  std::vector<std::vector<mpm::Index>> set_face_nodes;
  //! Set number of faces from element
  for (unsigned face_id = 0; face_id < element_->nfaces(); ++face_id) {
    std::vector<mpm::Index> face_nodes;

    // Get the nodes of the face
    const Eigen::VectorXi indices = element_->face_indices(face_id);
    for (int id = 0; id < indices.size(); ++id)
      face_nodes.emplace_back(nodes_[indices(id)]->id());

    // Sort in ascending order
    std::sort(face_nodes.begin(), face_nodes.end());
    set_face_nodes.emplace_back(face_nodes);
  }
  return set_face_nodes;
}
