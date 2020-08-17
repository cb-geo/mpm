#ifndef MPM_DISCONTINUITY_H_
#define MPM_DISCONTINUITY_H_

#include "cell.h"
#include "data_types.h"
#include "io_mesh.h"
#include "logger.h"
#include "memory.h"
#include "node_base.h"
#include "vector.h"

namespace mpm {

template <unsigned Tdim>
struct discontinuity_point;

//! class for describe the discontinuous surface
//! \brief
//! \tparam Tdim Dimension
template <unsigned Tdim>
class DiscontinuityBase {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Constructor with id
  //! \param[in] discontinuity id
  //! \param[in] discontinuity properties json
  DiscontinuityBase(unsigned id, const Json& discontinuity_props);

  //! Destructor
  virtual ~DiscontinuityBase(){};

  //! Delete copy constructor
  DiscontinuityBase(const DiscontinuityBase<Tdim>&) = delete;

  //! Delete assignement operator
  DiscontinuityBase& operator=(const DiscontinuityBase<Tdim>&) = delete;

  //！ initialization
  //! \param[in] the coordinates of all points
  //! \param[in] the point index of each surface
  virtual bool initialize(
      const std::vector<VectorDim>& points,
      const std::vector<std::vector<mpm::Index>>& surfs) = 0;

  //! create points from file
  //! \param[in] points the coordinates list of points
  bool create_points(const std::vector<VectorDim>& points);

  //! create elements from file
  //! \param[in] surfs the point index list of each surface
  virtual bool create_surfaces(const std::vector<std::vector<mpm::Index>>& surfs) {
    return true;
  };

  // return the levelset values of each coordinates
  //! \param[in] coordinates coordinates
  //! \param[in] phi_list the reference of phi for all coordinates
  virtual void compute_levelset(const std::vector<VectorDim>& coordinates,
                                std::vector<double>& phi_list) = 0;

  //! compute the normal vectors of coordinates
  //! \param[in] coordinates The coordinates
  //! \param[in] normal vector the normal vector of the given coordinates
  virtual void compute_normal(const VectorDim& coordinates,
                              VectorDim& normal_vector) = 0;

  //! return self_contact
  bool self_contact() const { return self_contact_; };

  //! return the friction coefficient
  double friction_coef() const { return friction_coef_; };

  //! return the number of the points
  mpm::Index npoints() const { return points_.size(); };

  //! Locate points in a cell
  //! \param[in] cells vector of cells
  //! \param[in] map_cells map of cells
  void locate_discontinuity_mesh(const Vector<Cell<Tdim>>& cells,
                                 const Map<Cell<Tdim>>& map_cells) noexcept;

  //! Compute updated position
  //! \param[in] dt Time-step
  void compute_updated_position(const double dt) noexcept;

  //! Compute shape function
  void compute_shapefn() noexcept;

  //! Assign point friction coefficient
  virtual void assign_point_friction_coef() noexcept = 0;

 protected:
  //! Logger
  std::unique_ptr<spdlog::logger> console_;

  //! vector of points
  std::vector<mpm::discontinuity_point<Tdim>> points_;

  //! self-contact
  bool self_contact_{true};

  //! friction coefficient
  double friction_coef_;

};  // DiscontinuityBase class


//! struct of discontinuity point
template <unsigned Tdim>
struct discontinuity_point {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! construct with coordinate
  discontinuity_point(const VectorDim& coordinate) {

    friction_coef_ = 0;
    coordinates_ = coordinate;
    cell_ = nullptr;
    //! Logger
    console_ = spdlog::get("discontinuity_point");
  }

  //! Return coordinates
  //! \retval return coordinates
  VectorDim coordinates() const { return coordinates_; }

  //! Return cell_id
  //! \retval return cell_id_
  Index cell_id() const { return cell_id_; }

  //! Assign a cell to point
  //! \param[in] cellptr Pointer to a cell
  //! \param[in] xi Local coordinates of the point in reference cell
  bool assign_cell_xi(const std::shared_ptr<Cell<Tdim>>& cellptr,
                      const Eigen::Matrix<double, Tdim, 1>& xi);

  //! Return cell ptr status
  bool cell_ptr() const { return cell_ != nullptr; }

  //! Assign a cell to point
  //! \param[in] cellptr Pointer to a cell
  bool assign_cell(const std::shared_ptr<Cell<Tdim>>& cellptr);

  //! Compute reference coordinates in a cell
  bool compute_reference_location() noexcept;

  //! Locate particles in a cell
  //! \param[in] cells vector of cells
  //! \param[in] map_cells map of cells
  void locate_discontinuity_mesh(const Vector<Cell<Tdim>>& cells,
                                const Map<Cell<Tdim>>& map_cells) noexcept;

  //! Compute updated position
  void compute_updated_position(double dt) noexcept;

  //! Compute shape function
  void compute_shapefn() noexcept;

  //! Assign point friction coefficient
  //! \param[in] friction_coef
  void assign_friction_coef(double friction_coef) noexcept {friction_coef_ = friction_coef;};

 private:
  //! point coordinates
  VectorDim coordinates_;
  //! Cell id
  Index cell_id_{std::numeric_limits<Index>::max()};
  //! Shape functions
  Eigen::VectorXd shapefn_;
  //! Cell
  std::shared_ptr<Cell<Tdim>> cell_;
  //! Reference coordinates (in a cell)
  Eigen::Matrix<double, Tdim, 1> xi_;
  //! Vector of nodal pointers
  std::vector<std::shared_ptr<NodeBase<Tdim>>> nodes_;
  //! Logger
  std::shared_ptr<spdlog::logger> console_;
  //! friction coefficient
  double friction_coef_{0.};
};

//! struct of discontinuity line: for 2d, need to be done
struct discontinuity_line {
 public:
  //! Return points indices
  Eigen::Matrix<int, 2, 1> points() const { return points_; };

 private:
  //! points index of the line
  Eigen::Matrix<int, 2, 1> points_;
};

//! struct of discontinuity surface: triangle
template <unsigned Tdim>
struct discontinuity_surface {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! construct with points indices
  discontinuity_surface(const std::vector<mpm::Index>& points) {
    for (int i = 0; i < points.size(); ++i) points_[i] = points[i];
  }
  //! Return points indices
  Eigen::Matrix<mpm::Index, 3, 1> points() const { return points_; }

  //! assign the surface center
  //! \param[in] center coordinates of the surface center
  inline void assign_center(VectorDim& center) { center_ = center; }

  //! assign the surface normal vector
  //! \param[in] normal normal vector of the surface
  inline void assign_normal(VectorDim& normal) { normal_ = normal; }

  //! Reture normal of the elements
  VectorDim normal() const { return normal_; }

  //! Return the vertical distance to the surface
  //! \param[in]  coor coordinates
  double vertical_distance(const VectorDim& coor) const {
    return (coor[0] - center_[0]) * normal_[0] +
           (coor[1] - center_[1]) * normal_[1] +
           (coor[2] - center_[2]) * normal_[2];
  };

 private:
  //! points indices
  Eigen::Matrix<mpm::Index, 3, 1> points_;

  // the center coordinates
  VectorDim center_;

  // the normal vector
  VectorDim normal_;
};

}  // namespace mpm

#include "discontinuity_base.tcc"
#include "discontinuity_point.tcc"

#endif  // MPM_DiscontinuityBase_H_
