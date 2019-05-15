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
bool mpm::Cell<Tdim>::add_neighbour(
    unsigned local_id, const std::shared_ptr<mpm::Cell<Tdim>>& cell_ptr) {
  bool insertion_status = false;
  try {
    // If number of cell ptrs id is not the current cell id
    if (cell_ptr->id() != this->id()) {
      insertion_status = neighbour_cells_.insert(local_id, cell_ptr);
    } else {
      throw std::runtime_error("Invalid local id of a cell neighbour");
    }
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
//! Computes the volume of a quadrilateral
template <>
inline void mpm::Cell<2>::compute_volume() {
  try {
    Eigen::VectorXi indices = element_->corner_indices();
    // Quadrilateral
    if (indices.size() == 4) {

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
    const Eigen::Matrix<double, Tdim, 1>& point) {

  // Check if point is approximately in the cell
  if (!this->approx_point_in_cell(point)) return false;

  // Check if cell is cartesian, if so use cartesian checker
  if (!isoparametric_) return mpm::Cell<Tdim>::point_in_cartesian_cell(point);

  bool status = true;
  // Get local coordinates
  Eigen::Matrix<double, Tdim, 1> xi = this->transform_real_to_unit_cell(point);
  // Check if the transformed coordinate is within the unit cell (-1, 1)
  for (unsigned i = 0; i < xi.size(); ++i)
    if (xi(i) < -1. || xi(i) > 1. || std::isnan(xi(i))) status = false;
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

    // Quadrilateral
    if (indices.size() == 4) {
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

  // Check for an analytical solution
  const long double x = point(0);
  const long double y = point(1);

  const double x0 = nodal_coords(0, 0);
  const double x1 = nodal_coords(0, 1);
  const double x2 = nodal_coords(0, 2);
  const double x3 = nodal_coords(0, 3);

  const double y0 = nodal_coords(1, 0);
  const double y1 = nodal_coords(1, 1);
  const double y2 = nodal_coords(1, 2);
  const double y3 = nodal_coords(1, 3);

  const long double a = (x1 - x3) * (y0 - y2) - (x0 - x2) * (y1 - y3);
  const long double b = -(x0 - x1 - x2 + x3) * y + (x - 2 * x1 + x3) * y0 -
                        (x - 2 * x0 + x2) * y1 - (x - x1) * y2 + (x - x0) * y3;
  const long double c = (x0 - x1) * y - (x - x1) * y0 + (x - x0) * y1;

  const long double discriminant = b * b - 4 * a * c;

  // Discriminant is negative if the point is not in the cell
  if (discriminant > 0.0) {
    long double eta1, eta2;
    // Special case #1: if a is zero, then use the linear formula
    if (a == 0.0 && b != 0.0) {
      eta1 = -c / b;
      eta2 = -c / b;
    }
    // Special case #2: a is zero for parallelograms and very small for
    // near-parallelograms:
    else if (std::abs(a) < 1e-8 * std::abs(b)) {
      // if both a and c are very small then the root should be near
      // zero: this first case will capture that
      eta1 = 2 * c / (-b - std::sqrt(discriminant));
      eta2 = 2 * c / (-b + std::sqrt(discriminant));
    }
    // finally, use the plain version:
    else {
      eta1 = (-b - std::sqrt(discriminant)) / (2 * a);
      eta2 = (-b + std::sqrt(discriminant)) / (2 * a);
    }

    // pick the one closer to the center of the cell.
    // eta
    const long double eta =
        (std::abs(eta1 - 0.5) < std::abs(eta2 - 0.5)) ? eta1 : eta2;
    // Scale from 0 to 1 -> -1 to 1
    xi(1) = (2. * eta) - 1.;

    // There are two ways to compute xi from eta, but either one may have a
    // zero denominator.
    const long double subexpr0 = -eta * x2 + x0 * (eta - 1);
    const long double xi_denominator0 = eta * x3 - x1 * (eta - 1) + subexpr0;
    const double max_x = std::max(std::max(std::abs(x0), std::abs(x1)),
                                  std::max(std::abs(x2), std::abs(x3)));

    if (std::abs(xi_denominator0) > 1.0E-10 * max_x) {
      const long double xi0 = (x + subexpr0) / xi_denominator0;
      xi(0) = (2. * xi0) - 1.;
      bool status = true;
      // Check if xi is within the cell
      for (unsigned i = 0; i < xi.size(); ++i)
        if (xi(i) < -1. || xi(i) > 1. || std::isnan(xi(i))) status = false;

      if (status) return xi;
    } else {
      const double max_y = std::max(std::max(std::abs(y0), std::abs(y1)),
                                    std::max(std::abs(y2), std::abs(y3)));
      const long double subexpr1 = -eta * y2 + y0 * (eta - 1);
      const long double xi_denominator1 = eta * y3 - y1 * (eta - 1) + subexpr1;
      if (std::abs(xi_denominator1) > 1.0E-10 * max_y) {
        const long double xi0 = (subexpr1 + y) / xi_denominator1;
        xi(0) = (2. * xi0) - 1.;
        bool status = true;
        // Check if xi is within the cell
        for (unsigned i = 0; i < xi.size(); ++i)
          if (xi(i) < -1. || xi(i) > 1. || std::isnan(xi(i))) status = false;

        if (status) return xi;
      }
    }
  }

  // Affine guess of xi
  Eigen::Matrix<double, 2, 1> affine_guess;
  // Boolean to check if affine is nan
  bool affine_nan = false;
  // Zeros
  const Eigen::Matrix<double, 2, 1> zero = Eigen::Matrix<double, 2, 1>::Zero();

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
    return affine_residual.norm() < nr_residual.norm() ? affine_guess : xi;

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

//! Map mass and momentum of a particle in a subdomain to nodes for a given
//! phase
template <unsigned Tdim>
void mpm::Cell<Tdim>::map_mass_momentum_to_nodes_subdomain(
    const Eigen::VectorXd& shapefn, unsigned phase, double pmass,
    const Eigen::VectorXd& pvelocity, unsigned mid) {

  for (unsigned i = 0; i < this->nfunctions(); ++i) {
    nodes_[i]->update_mass_subdomain(true, phase, shapefn(i) * pmass, mid);
    nodes_[i]->update_momentum_subdomain(true, phase,
                                         shapefn(i) * pmass * pvelocity, mid);
  }
}

//! Map particle coordinates to nodes
template <unsigned Tdim>
void mpm::Cell<Tdim>::map_coordinates_to_nodes(
    const Eigen::VectorXd& shapefn, double pmass,
    const Eigen::VectorXd& coordinates) {
  for (unsigned i = 0; i < this->nfunctions(); ++i) {
    nodes_[i]->update_coordinates(true, pmass, shapefn(i) * coordinates);
  }
}

//! Map coordinates of a particle in a subdomain to nodes
template <unsigned Tdim>
void mpm::Cell<Tdim>::map_coordinates_to_nodes_subdomain(
    const Eigen::VectorXd& shapefn, double pmass,
    const Eigen::VectorXd& coordinates, unsigned mid) {
  for (unsigned i = 0; i < this->nfunctions(); ++i) {
    nodes_[i]->update_coordinates_subdomain(true, pmass,
                                            shapefn(i) * coordinates, mid);
  }
}

//! Map particle pressure to nodes for a given phase
template <unsigned Tdim>
void mpm::Cell<Tdim>::map_pressure_to_nodes(const Eigen::VectorXd& shapefn,
                                            unsigned phase, double pmass,
                                            double ppressure) {

  for (unsigned i = 0; i < this->nfunctions(); ++i)
    nodes_[i]->update_pressure(true, phase, shapefn(i) * pmass * ppressure);
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

//! Compute strain rate of a subdomain
template <unsigned Tdim>
Eigen::VectorXd mpm::Cell<Tdim>::compute_strain_rate_subdomain(
    const std::vector<Eigen::MatrixXd>& bmatrix, unsigned phase, unsigned mid) {
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
      strain_rate += bmatrix.at(i) * nodes_[i]->velocity_subdomain(phase, mid);
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

//! Compute strain rate of a subdomain for reduced integration at the centroid
//! of cell
template <unsigned Tdim>
Eigen::VectorXd mpm::Cell<Tdim>::compute_strain_rate_centroid_subdomain(
    unsigned phase, unsigned mid) {

  // Define strain rate at centroid
  Eigen::Matrix<double, Tdof, 1> strain_rate_centroid =
      Eigen::Matrix<double, Tdof, 1>::Zero(bmatrix_centroid_.at(0).rows());

  // Compute strain rate
  for (unsigned i = 0; i < this->nnodes(); ++i)
    strain_rate_centroid +=
        bmatrix_centroid_.at(i) * nodes_[i]->velocity_subdomain(phase, mid);
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

//! Compute the nodal body force of a subdomain of a cell from particle mass and
//! gravity
template <unsigned Tdim>
void mpm::Cell<Tdim>::compute_nodal_body_force_subdomain(
    const Eigen::VectorXd& shapefn, unsigned phase, double pmass,
    const VectorDim& pgravity, unsigned mid) {
  // Map subdomain external forces from particle to nodes
  for (unsigned i = 0; i < this->nfunctions(); ++i)
    nodes_[i]->update_external_force_subdomain(
        true, phase, (shapefn(i) * pgravity * pmass), mid);
}

//! Compute the nodal traction force of a cell from particle
template <unsigned Tdim>
void mpm::Cell<Tdim>::compute_nodal_traction_force(
    const Eigen::VectorXd& shapefn, unsigned phase, const VectorDim& traction) {
  // Map external forces from particle to nodes
  for (unsigned i = 0; i < this->nfunctions(); ++i)
    nodes_[i]->update_external_force(true, phase, (shapefn(i) * traction));
}

//! Compute the nodal traction force of a subdomain of a cell from the particle
template <unsigned Tdim>
void mpm::Cell<Tdim>::compute_nodal_traction_force_subdomain(
    const Eigen::VectorXd& shapefn, unsigned phase, const VectorDim& traction,
    unsigned mid) {
  // Map external forces from particle to nodes
  for (unsigned i = 0; i < this->nfunctions(); ++i)
    nodes_[i]->update_external_force_subdomain(true, phase,
                                               (shapefn(i) * traction), mid);
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

//! Compute the noal internal forcee of a subdomain of a cell from particle
//! stress and volume
template <unsigned Tdim>
inline void mpm::Cell<Tdim>::compute_nodal_internal_force_subdomain(
    const std::vector<Eigen::MatrixXd>& bmatrix, unsigned phase, double pvolume,
    const Eigen::Matrix<double, 6, 1>& pstress, unsigned mid) {
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
    nodes_[j]->update_internal_force_subdomain(
        true, phase, (pvolume * bmatrix.at(j).transpose() * stress), mid);
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

//! Return velocity of a subdomain at given location by interpolating from nodes
template <unsigned Tdim>
Eigen::Matrix<double, Tdim, 1>
    mpm::Cell<Tdim>::interpolate_nodal_velocity_subdomain(
        const Eigen::VectorXd& shapefn, unsigned phase, unsigned mid) {
  Eigen::Matrix<double, Tdim, 1> velocity;
  velocity.setZero();

  for (unsigned i = 0; i < this->nfunctions(); ++i)
    velocity += shapefn(i) * nodes_[i]->velocity_subdomain(phase, mid);

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

//! Return acceleration of a subdomain at given location by interpolating from
//! nodes
template <unsigned Tdim>
Eigen::Matrix<double, Tdim, 1>
    mpm::Cell<Tdim>::interpolate_nodal_acceleration_subdomain(
        const Eigen::VectorXd& shapefn, unsigned phase, unsigned mid) {
  Eigen::Matrix<double, Tdim, 1> acceleration;
  acceleration.setZero();

  for (unsigned i = 0; i < this->nfunctions(); ++i)
    acceleration += shapefn(i) * nodes_[i]->acceleration_subdomain(phase, mid);

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
