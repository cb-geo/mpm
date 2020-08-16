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
struct discontinuity_point {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  discontinuity_point(const VectorDim& coordinate) {
    coordinates_ = coordinate;
    cell_ = nullptr;
    //! Logger
    console_ = spdlog::get("discontinuity_point");
  }

  //! Return coordinates
  //! \retval coordinates_ return coordinates of the nodebase
  VectorDim coordinates() const { return coordinates_; }

  //! Return cell_id
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

  //! Locate points in a cell
  void locate_discontinuity_mesh(Vector<Cell<Tdim>>& cells,
                                 Map<Cell<Tdim>>& map_cells) noexcept;

  //! Compute updated position
  void compute_updated_position(double dt) noexcept;

  //! Compute shape function
  void compute_shapefn() noexcept;

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
};

//! class for to describe the discontinuous surface
//! \brief
//! \details nodes, lines and areas
//! \tparam Tdim Dimension
template <unsigned Tdim>
class DiscontinuityBase {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Constructor with id
  //! \param[in] discontinuity_properties discontinuity properties
  DiscontinuityBase(unsigned id, const Json& discontinuity_props);

  //! Destructor
  virtual ~DiscontinuityBase(){};

  //! Delete copy constructor
  DiscontinuityBase(const DiscontinuityBase<Tdim>&) = delete;

  //! Delete assignement operator
  DiscontinuityBase& operator=(const DiscontinuityBase<Tdim>&) = delete;

  // initialization
  virtual bool initialize(
      const std::vector<VectorDim>& coordinates,
      const std::vector<std::vector<mpm::Index>>& pointsets) = 0;

  //! create points from file
  bool create_points(const std::vector<VectorDim>& coordinates);

  //! create elements from file
  virtual bool create_elements(
      const std::vector<std::vector<mpm::Index>>& elements) {
    return true;
  };

  // return the levelset values of each doordinates
  //! \param[in] the vector of the coordinates
  virtual void compute_levelset(const std::vector<VectorDim>& coordinates,
                                std::vector<double>& phi_list) = 0;

  // return the normal vectors of given coordinates
//! \param[in] the coordinates
virtual void compute_normal(
  const  VectorDim& coordinates, VectorDim& normal_vector) = 0;

  // return self_contact
  bool self_contact() { return self_contact_; };

  // return the friction coefficient
  double friction_coef() { return friction_coef_; };

  // return the number of the points
  mpm::Index npoints() { return points_.size(); };

  void points_list(std::vector<mpm::discontinuity_point<Tdim>>& points) {
    points = points_;
  }
  //! Locate points in a cell
  void locate_discontinuity_mesh(Vector<Cell<Tdim>>& cells,
                                 Map<Cell<Tdim>>& map_cells) noexcept;

  //! Compute updated position
  void compute_updated_position(double dt) noexcept;

  //! Compute shape function
  void compute_shapefn() noexcept;

 protected:
  //! Logger
  std::unique_ptr<spdlog::logger> console_;

  std::vector<mpm::discontinuity_point<Tdim>> points_;

  // number of points
  mpm::Index numpoint_; //delete

  // self-contact
  bool self_contact_{true};

  double friction_coef_;

};  // DiscontinuityBase class

struct discontinuity_line {
 public:
  //! Return points indices
  Eigen::Matrix<int, 2, 1> points() const { return points_; };

 private:
  //! points index of the line
  Eigen::Matrix<int, 2, 1> points_;
};

template <unsigned Tdim>
struct discontinuity_element {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  discontinuity_element(const std::vector<mpm::Index>& points) {
    for (int i = 0; i < points.size(); ++i) points_[i] = points[i];
  }
  //! Return points indices
  Eigen::Matrix<mpm::Index, 3, 1> points() const { return points_; }

  inline void set_center(VectorDim& center) { center_ = center; }

  inline void set_normal(VectorDim& normal) { normal_ = normal; }
  
  //! Reture normal of the elements
  VectorDim normal() const {return normal_;}

  double Vertical_distance(const VectorDim& coor) const {
    return (coor[0] - center_[0]) * normal_[0] +
           (coor[1] - center_[1]) * normal_[1] +
           (coor[2] - center_[2]) * normal_[2];
  };

 private:
  //! points indices
  Eigen::Matrix<mpm::Index, 3, 1> points_;

  // the center of the triangular elements
  VectorDim center_;

  // the normal of the triangular elements
  VectorDim normal_;
};

}  // namespace mpm

#include "discontinuity_base.tcc"

#endif  // MPM_DiscontinuityBase_H_
