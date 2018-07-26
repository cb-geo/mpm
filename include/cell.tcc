//! Constructor with cell id, number of nodes and shapefn
template <unsigned Tdim>
mpm::Cell<Tdim>::Cell(Index id, unsigned nnodes,
                      const std::shared_ptr<ShapeFn<Tdim>>& shapefnptr)
    : id_{id}, nnodes_{nnodes} {
  // Check if the dimension is between 1 & 3
  static_assert((Tdim >= 1 && Tdim <= 3), "Invalid global dimension");

  try {
    if (shapefnptr->nfunctions() >= this->nnodes_) {
      shapefn_ = shapefnptr;
    } else {
      throw std::runtime_error(
          "Specified number of shape functions is not defined");
    }
    if (!(nnodes > Tdim)) {
      throw std::runtime_error(
          "Specified number of nodes for a cell is too low");
    }
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
    std::abort();
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
      status = true;
    } else {
      throw std::runtime_error(
          "Specified number of nodes for a cell is not present");
    }
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
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
          // Check if mean lenght of a cell is initialised
          (std::fabs(this->mean_length_ - std::numeric_limits<double>::max()) >
           1.0E-10));
}

//! Assign a shape function to cell
template <unsigned Tdim>
bool mpm::Cell<Tdim>::shapefn(
    const std::shared_ptr<ShapeFn<Tdim>>& shapefnptr) {
  bool status = false;
  try {
    if (shapefnptr->nfunctions() >= this->nnodes_) {
      shapefn_ = shapefnptr;
      status = true;
    } else {
      throw std::runtime_error(
          "Specified number of shape functions is not defined");
    }
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
  }
  return status;
}

//! Add a node pointer and return the status of addition of a node
template <unsigned Tdim>
bool mpm::Cell<Tdim>::add_node(unsigned local_id,
                               std::shared_ptr<mpm::NodeBase<Tdim>> node_ptr) {
  bool insertion_status = false;
  try {
    // If number of node ptrs in a cell is less than the maximum number of nodes
    // per cell
    // The local id should be between 0 and maximum number of nodes
    if (nodes_.size() < this->nnodes_ &&
        (local_id >= 0 && local_id < this->nnodes_)) {
      insertion_status = nodes_.insert(local_id, node_ptr);
    } else {
      throw std::runtime_error(
          "Number nodes in a cell exceeds the maximum allowed per cell");
    }
  } catch (std::exception& exception) {
    std::cerr << exception.what() << "\n";
  }
  return insertion_status;
}

//! Activate nodes if particle is present
template <unsigned Tdim>
bool mpm::Cell<Tdim>::activate_nodes() {
  bool status = true;
  try {
    // If number of particles are present, set node status to active
    if (particles_.size() > 0) {
      // Activate all nodes
      for (unsigned i = 0; i < nodes_.size(); ++i)
        nodes_[i]->assign_status(true);
    } else {
      throw std::runtime_error("No particles in cell, can't activate nodes");
    }
  } catch (std::exception& exception) {
    std::cerr << exception.what() << "\n";
    status = false;
  }
  return status;
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
    std::cerr << exception.what() << "\n";
  }
  return insertion_status;
}

//! Add a particle id and return the status of addition of a particle id
template <unsigned Tdim>
bool mpm::Cell<Tdim>::add_particle_id(Index id) {
  bool status = false;
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
  particles_.erase(std::remove(particles_.begin(), particles_.end(), id),
                   particles_.end());
}

//! Compute volume of a 2D cell
//! Computes the volume of a quadrilateral
template <>
inline void mpm::Cell<2>::compute_volume() {
  try {
    Eigen::VectorXi indices = shapefn_->corner_indices();
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
      const double a = (nodes_[indices[0]]->coordinates() -
                        nodes_[indices[3]]->coordinates())
                           .norm();
      const double b = (nodes_[indices[2]]->coordinates() -
                        nodes_[indices[3]]->coordinates())
                           .norm();
      const double c = (nodes_[indices[1]]->coordinates() -
                        nodes_[indices[2]]->coordinates())
                           .norm();
      const double d = (nodes_[indices[0]]->coordinates() -
                        nodes_[indices[1]]->coordinates())
                           .norm();
      const double p = (nodes_[indices[0]]->coordinates() -
                        nodes_[indices[2]]->coordinates())
                           .norm();
      const double q = (nodes_[indices[1]]->coordinates() -
                        nodes_[indices[3]]->coordinates())
                           .norm();

      // K = 1/4 * sqrt ( 4p^2q^2 - (a^2 + c^2 - b^2 -d^2)^2)
      volume_ =
          0.25 * std::sqrt(4 * p * p * q * q -
                           std::pow((a * a + c * c - b * b - d * d), 2.0));
    } else {
      throw std::runtime_error(
          "Unable to compute volume, number of vertices is incorrect");
    }
  } catch (std::exception& except) {
    std::cout << __FILE__ << __LINE__
              << "Compute volume of a cell: " << except.what() << '\n';
  }
}

//! Compute volume of a 3D cell
//! Computes the volume of a hexahedron
template <>
inline void mpm::Cell<3>::compute_volume() {
  try {
    Eigen::VectorXi indices = shapefn_->corner_indices();
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

      const auto a = nodes_[indices[7]]->coordinates();
      const auto b = nodes_[indices[6]]->coordinates();
      const auto c = nodes_[indices[2]]->coordinates();
      const auto d = nodes_[indices[3]]->coordinates();
      const auto e = nodes_[indices[4]]->coordinates();
      const auto f = nodes_[indices[5]]->coordinates();
      const auto g = nodes_[indices[1]]->coordinates();
      const auto h = nodes_[indices[0]]->coordinates();

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
  } catch (std::exception& except) {
    std::cout << __FILE__ << __LINE__
              << "Compute volume of a cell: " << except.what() << '\n';
  }
}

//! Compute the centroid of the cell
template <unsigned Tdim>
void mpm::Cell<Tdim>::compute_centroid() {
  // Get indices of corner nodes
  Eigen::VectorXi indices = shapefn_->corner_indices();

  // Calculate the centroid of the cell
  centroid_.setZero();
  for (unsigned i = 0; i < indices.size(); ++i)
    centroid_ += nodes_[indices[i]]->coordinates();

  centroid_ /= indices.size();
}

//! Compute mean length of cell
template <unsigned Tdim>
void mpm::Cell<Tdim>::compute_mean_length() {
  // Get the indices of sub-triangles
  Eigen::MatrixXi indices = shapefn_->sides_indices();
  this->mean_length_ = 0.;
  // Calculate the mean length
  for (unsigned i = 0; i < indices.rows(); ++i)
    this->mean_length_ += (nodes_[indices(i, 0)]->coordinates() -
                           nodes_[indices(i, 1)]->coordinates())
                              .norm();
  this->mean_length_ /= indices.rows();
}

//! Check if a point is in a 2D cell by breaking the cell into sub-volumes
template <unsigned Tdim>
inline bool mpm::Cell<Tdim>::point_in_cell(
    const Eigen::Matrix<double, 2, 1>& point) {

  // Tolerance for volume / area comparison
  const double tolerance = 1.0E-10;

  if (std::fabs(volume_ - std::numeric_limits<double>::max()) < tolerance)
    this->compute_volume();

  // Get the indices of sub-triangles
  Eigen::MatrixXi indices = shapefn_->inhedron_indices();

  // Initialise sub-triangle areas to 0
  double triareas = 0.;

  // Return status
  bool status = false;

  // Iterate over each sub-triangles to calculate the area
  // of each sub-triangle. If the sum of area of all
  // sub-triangle is equal to the area of hexahedron then
  // the point is inside the Hexahedron, if not the point is outside
  for (unsigned i = 0; i < indices.rows(); ++i) {

    // Get the 2 vertices of the sub-triangle
    const Eigen::Matrix<double, 2, 1> a = nodes_[indices(i, 0)]->coordinates();
    const Eigen::Matrix<double, 2, 1> b = nodes_[indices(i, 1)]->coordinates();

    // Compute the area of a triangle
    // Area = |1/2(x1(y2−y3)+x2(y3−y1)+x3(y1−y2))|
    Eigen::Matrix<double, 3, 3> area;
    // clang-format off
    area << 1.  , 1.  , 1.,
            a(0), b(0), point(0),
            a(1), b(1), point(1);
    // clang-format on
    const double triarea = 0.5 * std::fabs(area.determinant());

    triareas += triarea;

    // Optimisation check, if the sub-tetrahedra area exceeds the area of
    // hexahedron, abort and return false (point is outside).
    if ((triareas > volume_) && (std::fabs(triareas - volume_) > tolerance)) {
      return false;
    }
  }

  // Check if the point is inside the hedron
  if (std::fabs(triareas - volume_) < tolerance) {
    status = true;
  }
  return status;
}

//! Check if a point is in a 3D cell by breaking the cell into sub-volumes
template <unsigned Tdim>
inline bool mpm::Cell<Tdim>::point_in_cell(
    const Eigen::Matrix<double, 3, 1>& point) {

  // Tolerance for volume / area comparison
  const double tolerance = 1.0E-10;

  if (std::fabs(volume_ - std::numeric_limits<double>::max()) < tolerance)
    this->compute_volume();

  // Get the indices of sub-tetrahedron
  Eigen::MatrixXi indices = shapefn_->inhedron_indices();

  // Initialise sub-tetrahedra volumes to 0
  double tetvolumes = 0.;

  // Return status
  bool status = false;

  // Iterate over each sub-tetrahedrons to calculate the volume
  // of each sub-tetrahedron. If the sum of volume of all
  // sub-tetrahedron is equal to the volume of hexahedron then
  // the point is inside the Hexahedron, if not the point is outside
  for (unsigned i = 0; i < indices.rows(); ++i) {

    // Get the 3 vertices of the sub-tetrahedron
    const Eigen::Matrix<double, 3, 1> a = nodes_[indices(i, 0)]->coordinates();
    const Eigen::Matrix<double, 3, 1> b = nodes_[indices(i, 1)]->coordinates();
    const Eigen::Matrix<double, 3, 1> c = nodes_[indices(i, 2)]->coordinates();

    // Compute the volume of a tetrahedron
    // Volume = 1/6 | (a - d).((b - d) x (c - d)) |
    const double tetvolume =
        (1. / 6.) * std::fabs((a - point).dot((b - point).cross(c - point)));

    tetvolumes += tetvolume;
    // Optimisation check, if the sub-tetrahedra volume exceeds the volume of
    // hexahedron, abort and return false (point is outside).
    if ((tetvolumes > volume_) &&
        (std::fabs(tetvolumes - volume_) > tolerance)) {
      return false;
    }
  }
  // Check if the point is inside the hedron
  if (std::fabs(tetvolumes - volume_) < tolerance) {
    status = true;
  }
  return status;
}

//! Check if a point is in a 3D cell by affine transformation and newton-raphson
template <unsigned Tdim>
inline bool mpm::Cell<Tdim>::is_point_in_cell(
    const Eigen::Matrix<double, Tdim, 1>& point) {
  bool status = true;
  // Get local coordinates
  Eigen::Matrix<double, Tdim, 1> xi = this->transform_real_to_unit_cell(point);
  // Check if the transformed coordinate is within the unit cell (-1, 1)
  for (unsigned i = 0; i < xi.size(); ++i)
    if (xi(i) < -1. || xi(i) > 1.) status = false;

  return status;
}
//! Return the local coordinates of a point in a 1D cell
template <unsigned Tdim>
inline Eigen::Matrix<double, 1, 1> mpm::Cell<Tdim>::local_coordinates_point(
    const Eigen::Matrix<double, 1, 1>& point) {
  // Local point coordinates
  Eigen::Matrix<double, 1, 1> xi;

  try {
    // Indices of corner nodes
    Eigen::VectorXi indices = shapefn_->corner_indices();

    // Linear
    if (indices.size() == 2) {
      //
      // 0 0---------0 1
      //        l
      const double length = (nodes_[indices[0]]->coordinates() -
                             nodes_[indices[1]]->coordinates())
                                .norm();

      const Eigen::Matrix<double, 1, 1> centre =
          (nodes_[indices[0]]->coordinates() +
           nodes_[indices[1]]->coordinates()) /
          2.0;

      xi(0) = 2. * (point(0) - centre(0)) / length;
    } else {
      throw std::runtime_error("Unable to compute local coordinates");
    }
  } catch (std::exception& except) {
    std::cout << __FILE__ << __LINE__
              << "Compute local coordinate of a point in a cell: "
              << except.what() << '\n';
  }
  return xi;
}

//! Return the local coordinates of a point in a 2D cell
template <unsigned Tdim>
inline Eigen::Matrix<double, 2, 1> mpm::Cell<Tdim>::local_coordinates_point(
    const Eigen::Matrix<double, 2, 1>& point) {
  // Local point coordinates
  Eigen::Matrix<double, 2, 1> xi;

  try {
    // Indices of corner nodes
    Eigen::VectorXi indices = shapefn_->corner_indices();

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
      const double xlength = (nodes_[indices[0]]->coordinates() -
                              nodes_[indices[1]]->coordinates())
                                 .norm();
      const double ylength = (nodes_[indices[1]]->coordinates() -
                              nodes_[indices[2]]->coordinates())
                                 .norm();

      const Eigen::Matrix<double, 2, 1> centre =
          (nodes_[indices[0]]->coordinates() +
           nodes_[indices[1]]->coordinates() +
           nodes_[indices[2]]->coordinates() +
           nodes_[indices[3]]->coordinates()) /
          4.0;

      xi(0) = 2. * (point(0) - centre(0)) / xlength;
      xi(1) = 2. * (point(1) - centre(1)) / ylength;
    } else {
      throw std::runtime_error("Unable to compute local coordinates");
    }
  } catch (std::exception& except) {
    std::cout << __FILE__ << __LINE__
              << "Compute local coordinate of a point in a cell: "
              << except.what() << '\n';
  }
  return xi;
}

//! Return the local coordinates of a point in a 3D cell
template <unsigned Tdim>
inline Eigen::Matrix<double, 3, 1> mpm::Cell<Tdim>::local_coordinates_point(
    const Eigen::Matrix<double, 3, 1>& point) {

  // Local point coordinates
  Eigen::Matrix<double, 3, 1> xi;

  try {
    // Indices of corner nodes
    Eigen::VectorXi indices = shapefn_->corner_indices();

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

      const double xlength = (nodes_[indices[0]]->coordinates() -
                              nodes_[indices[1]]->coordinates())
                                 .norm();
      const double ylength = (nodes_[indices[1]]->coordinates() -
                              nodes_[indices[2]]->coordinates())
                                 .norm();
      const double zlength = (nodes_[indices[1]]->coordinates() -
                              nodes_[indices[5]]->coordinates())
                                 .norm();

      // Compute centre
      Eigen::Matrix<double, 3, 1> centre;
      centre.setZero();
      for (unsigned i = 0; i < indices.size(); ++i)
        centre += nodes_[indices[i]]->coordinates();
      centre /= indices.size();

      xi(0) = 2. * (point(0) - centre(0)) / xlength;
      xi(1) = 2. * (point(1) - centre(1)) / ylength;
      xi(2) = 2. * (point(2) - centre(2)) / ylength;

    } else {
      throw std::runtime_error("Unable to compute local coordinates");
    }
  } catch (std::exception& except) {
    std::cout << __FILE__ << __LINE__
              << "Compute local coordinate of a point in a cell: "
              << except.what() << '\n';
  }
  return xi;
}

//! Return the local coordinates of a point in a 1D cell
template <unsigned Tdim>
inline Eigen::Matrix<double, 1, 1> mpm::Cell<Tdim>::transform_real_to_unit_cell(
    const Eigen::Matrix<double, 1, 1>& point) {
  return this->local_coordinates_point(point);
}

//! Return the local coordinates of a point in a 2D/3D cell
template <unsigned Tdim>
inline Eigen::Matrix<double, 2, 1> mpm::Cell<Tdim>::transform_real_to_unit_cell(
    const Eigen::Matrix<double, 2, 1>& point) {

  // Local coordinates of a point in an unit cell
  Eigen::Matrix<double, Tdim, 1> xi;
  xi.setZero();

  // Maximum iterations of newton raphson
  const unsigned max_iterations = 500;
  // Tolerance for newton raphson
  const double tolerance = 1.e-10;

  // Matrix of nodal coordinates
  Eigen::MatrixXd nodal_coords;
  nodal_coords.resize(Tdim, this->nfunctions());

  for (unsigned j = 0; j < this->nfunctions(); ++j) {
    Eigen::Matrix<double, Tdim, 1> node = nodes_[j]->coordinates();
    for (unsigned i = 0; i < Tdim; ++i) {
      nodal_coords(i, j) = node[i];
    }
  }

  // Coordinates of a unit cell
  const auto unit_cell = shapefn_->unit_cell_coordinates();

  // Affine transformation, using linear interpolation for the initial guess
  if (shapefn_->degree() == mpm::ShapeFnDegree::Linear) {
    // A = vertex * KA
    Eigen::Matrix<double, 2, 2> A;
    A = nodal_coords * mpm::TransformR2UAffine<2, 4>::KA;

    // b = vertex * Kb
    Eigen::Matrix<double, 2, 1> b =
        point - (nodal_coords * mpm::TransformR2UAffine<2, 4>::Kb);

    // Affine transform: A^-1 * b
    Eigen::Matrix<double, 2, 1> affine_guess = A.inverse() * b;

    // Check for nan
    bool guess_nan = false;
    for (unsigned i = 0; i < affine_guess.size(); ++i)
      if (std::isnan(std::fabs(affine_guess(i)))) bool guess_nan = true;

    // Set xi to affine guess
    if (!guess_nan) xi = affine_guess;

    // Shape function
    const auto sf = shapefn_->shapefn(xi);

    // f(x) = p(x) - p, where p is the real point
    Eigen::Matrix<double, Tdim, 1> fx = (nodal_coords * sf) - point;

    // Early exit
    if (fx.squaredNorm() < (1e-24 * this->mean_length_ * this->mean_length_))
      return xi;
  }

  // Newton Raphson iteration to solve for x
  // x_{n+1} = x_n - f(x)/f'(x)
  // f(x) = p(x) - p, where p is the real point
  // p(x) is the computed point.
  for (unsigned iter = 0; iter < max_iterations; ++iter) {

    // Calculate Jacobian
    Eigen::Matrix<double, Tdim, Tdim> jacobian;
    const auto grad_sf = shapefn_->grad_shapefn(xi);
    jacobian = unit_cell.transpose() * grad_sf;

    // Shape function
    const auto sf = shapefn_->shapefn(xi);

    // Residual (f(x))
    // f(x) = p(x) - p, where p is the real point
    Eigen::Matrix<double, Tdim, 1> residual = (nodal_coords * sf) - point;

    // x_{n+1} = x_n - f(x)/f'(x)
    xi -= (jacobian.inverse() * residual);

    // Convergence criteria
    if (residual.norm() < tolerance) break;
  }
  return xi;
}

//! Return the local coordinates of a point in a 2D/3D cell
template <unsigned Tdim>
inline Eigen::Matrix<double, 3, 1> mpm::Cell<Tdim>::transform_real_to_unit_cell(
    const Eigen::Matrix<double, 3, 1>& point) {

  // Local coordinates of a point in an unit cell
  Eigen::Matrix<double, Tdim, 1> xi;
  xi.setZero();

  // Maximum iterations of newton raphson
  const unsigned max_iterations = 100;
  // Tolerance for newton raphson
  const double tolerance = 1.e-11;

  // Matrix of nodal coordinates
  Eigen::MatrixXd nodal_coords;
  nodal_coords.resize(Tdim, this->nfunctions());

  for (unsigned j = 0; j < this->nfunctions(); ++j) {
    Eigen::Matrix<double, Tdim, 1> node = nodes_[j]->coordinates();
    for (unsigned i = 0; i < Tdim; ++i) {
      nodal_coords(i, j) = node[i];
    }
  }

  // Coordinates of a unit cell
  const auto unit_cell = shapefn_->unit_cell_coordinates();

  // Affine transformation, using linear interpolation for the initial guess
  if (shapefn_->degree() == mpm::ShapeFnDegree::Linear) {
    // A = vertex * KA
    Eigen::Matrix<double, 3, 3> A;
    A = nodal_coords * mpm::TransformR2UAffine<3, 8>::KA;

    // b = vertex * Kb
    Eigen::Matrix<double, 3, 1> b =
        point - (nodal_coords * mpm::TransformR2UAffine<3, 8>::Kb);

    // Affine transform: A^-1 * b
    Eigen::Matrix<double, 3, 1> affine_guess = A.inverse() * b;

    // Check for nan
    bool guess_nan = false;
    for (unsigned i = 0; i < affine_guess.size(); ++i)
      if (std::isnan(std::fabs(affine_guess(i)))) bool guess_nan = true;

    // Set xi to affine guess
    if (!guess_nan) xi = affine_guess;

    // Shape function
    const auto sf = shapefn_->shapefn(xi);

    // f(x) = p(x) - p, where p is the real point
    Eigen::Matrix<double, Tdim, 1> fx = (nodal_coords * sf) - point;

    // Early exit
    if (fx.squaredNorm() < (1e-24 * this->mean_length_ * this->mean_length_))
      return xi;
  }

  // Newton Raphson iteration to solve for x
  // x_{n+1} = x_n - f(x)/f'(x)
  // f(x) = p(x) - p, where p is the real point
  // p(x) is the computed point.
  for (unsigned iter = 0; iter < max_iterations; ++iter) {
    // Calculate Jacobian
    Eigen::Matrix<double, Tdim, Tdim> jacobian;
    const auto grad_sf = shapefn_->grad_shapefn(xi);
    jacobian = unit_cell.transpose() * grad_sf;

    // Shape function
    const auto sf = shapefn_->shapefn(xi);

    // Residual f(x)
    Eigen::Matrix<double, Tdim, 1> residual;
    // f(x) = p(x) - p
    residual = (nodal_coords * sf) - point;

    // x_{n+1} = x_n - f(x)/f'(x)
    xi -= (jacobian.inverse() * residual);

    // Convergence criteria
    if (residual.norm() < tolerance) break;
  }
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
void mpm::Cell<Tdim>::map_particle_volume_to_nodes(const VectorDim& xi,
                                                   unsigned phase,
                                                   double pvolume) {
  const auto shapefns = shapefn_->shapefn(xi);
  for (unsigned i = 0; i < shapefns.size(); ++i) {
    nodes_[i]->update_volume(true, phase, shapefns(i) * pvolume);
  }
}

//! Map particle mass and momentum to nodes for a given phase
template <unsigned Tdim>
void mpm::Cell<Tdim>::map_mass_momentum_to_nodes(
    const Eigen::VectorXd& shapefn, unsigned phase, double pmass,
    const Eigen::VectorXd& pvelocity) {

  // if (pmass.size() == pvelocity.cols()) mass = pmass.asDiagonal();
  for (unsigned i = 0; i < this->nfunctions(); ++i) {
    nodes_[i]->update_mass(true, phase, shapefn(i) * pmass);
    nodes_[i]->update_momentum(true, phase, shapefn(i) * pmass * pvelocity);
  }
}

//! Compute nodal momentum from particle mass and velocity for a given phase
template <unsigned Tdim>
void mpm::Cell<Tdim>::compute_nodal_momentum(const Eigen::VectorXd& shapefn,
                                             unsigned phase, double pmass,
                                             const Eigen::VectorXd& pvelocity) {

  // if (pmass.size() == pvelocity.cols()) mass = pmass.asDiagonal();
  for (unsigned i = 0; i < this->nfunctions(); ++i) {
    nodes_[i]->update_momentum(true, phase, shapefn(i) * pmass * pvelocity);
  }
}

//! Compute strain rate
template <unsigned Tdim>
Eigen::VectorXd mpm::Cell<Tdim>::compute_strain_rate(
    const std::vector<Eigen::MatrixXd>& bmatrix, unsigned phase) {
  // Define strain rate
  Eigen::VectorXd strain_rate;

  switch (Tdim) {
    case (1): {
      strain_rate.resize(1);
      break;
    }
    case (2): {
      strain_rate.resize(3);
      break;
    }
    default: {
      strain_rate.resize(6);
      break;
    }
  }

  strain_rate.setZero();

  try {
    // Check if B-Matrix size and number of nodes match
    if (this->nfunctions() != bmatrix.size() ||
        this->nnodes() != bmatrix.size())
      throw std::runtime_error(
          "Number of nodes / shapefn doesn't match BMatrix");

    for (unsigned i = 0; i < this->nnodes(); ++i) {
      Eigen::Matrix<double, Tdim, 1> node_velocity = nodes_[i]->velocity(phase);
      strain_rate += bmatrix.at(i) * node_velocity;
    }
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
  }
  return strain_rate;
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
Eigen::VectorXd mpm::Cell<Tdim>::interpolate_nodal_velocity(const VectorDim& xi,
                                                            unsigned phase) {
  Eigen::Matrix<double, Tdim, 1> velocity =
      Eigen::Matrix<double, Tdim, 1>::Zero();
  const auto shapefns = shapefn_->shapefn(xi);
  for (unsigned i = 0; i < this->nfunctions(); ++i)
    velocity += shapefns(i) * nodes_[i]->velocity(phase);

  return velocity;
}

//! Return acceleration at a point by interpolating from nodes
template <unsigned Tdim>
Eigen::VectorXd mpm::Cell<Tdim>::interpolate_nodal_acceleration(
    const VectorDim& xi, unsigned phase) {
  Eigen::Matrix<double, Tdim, 1> acceleration =
      Eigen::Matrix<double, Tdim, 1>::Zero();
  const auto shapefns = shapefn_->shapefn(xi);
  for (unsigned i = 0; i < this->nfunctions(); ++i)
    acceleration += shapefns(i) * nodes_[i]->acceleration(phase);

  return acceleration;
}
