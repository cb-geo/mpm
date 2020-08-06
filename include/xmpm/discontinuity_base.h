#ifndef MPM_DISCONTINUITY_H_
#define MPM_DISCONTINUITY_H_


#include "logger.h"

#include "data_types.h"

namespace mpm {

template <unsigned Tdim>
struct discontinuous_point {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  discontinuous_point(const VectorDim& coordinate) {
    coordinates_ = coordinate;
  }

  //! Return coordinates
  //! \retval coordinates_ return coordinates of the nodebase
  VectorDim coordinates() const { return coordinates_; }

 private:
  //! point coordinates
  VectorDim coordinates_;
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

  // Constructor
  DiscontinuityBase();

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

  bool self_contact() { return self_contact_; };

  void set_frictional_coef(double coef) { frictional_coef_ = coef; };

  double frictional_coef() { return frictional_coef_; };

 protected:
  std::vector<mpm::discontinuous_point<Tdim>> points_;

  // number of points
  mpm::Index numpoint_;

  //! Logger
  std::unique_ptr<spdlog::logger> console_;

  // self-contact
  bool self_contact_{true};

  double frictional_coef_;

};  // DiscontinuityBase class

struct discontinuous_line {
 public:
  //! Return points indices
  Eigen::Matrix<int, 2, 1> points() const { return points_; };

 private:
  //! points index of the line
  Eigen::Matrix<int, 2, 1> points_;
};

struct discontinuous_element {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, 3, 1>;

  discontinuous_element(const std::vector<mpm::Index>& points) {
    for (int i = 0; i < 3; ++i) points_[i] = points[i];
  }
  //! Return points indices
  Eigen::Matrix<mpm::Index, 3, 1> points() const { return points_; }

  inline void set_center(VectorDim& center) { center_ = center; }

  inline void set_normal(VectorDim& normal) { normal_ = normal; }

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

template <unsigned Tdim>
mpm::DiscontinuityBase<Tdim>::DiscontinuityBase() {
  numpoint_ = 0;

  frictional_coef_ = -1;

  std::string logger = "discontinuitybase";
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

//! create points from file
template <unsigned Tdim>
bool mpm::DiscontinuityBase<Tdim>::create_points(
    const std::vector<VectorDim>& coordinates) {
  bool status = true;
  try {
    // Check if point coordinates is empty
    if (coordinates.empty())
      throw std::runtime_error("List of coordinates is empty");
    // Iterate over all coordinates
    for (const auto& point_coordinates : coordinates) {

      // Add point
      mpm::discontinuous_point<Tdim> point(point_coordinates);

      points_.emplace_back(point);  //
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

#endif  // MPM_DiscontinuityBase_H_
