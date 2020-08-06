#ifndef MPM_DISCONTINUITY_3D_H_
#define MPM_DISCONTINUITY_3D_H_

#include "discontinuity_base.h"

//! MPM namespace
namespace mpm {

template <unsigned Tdim>
class Discontinuity_3D : public DiscontinuityBase<Tdim> {

 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;
  // constructor
  Discontinuity_3D();

  // initialization
  virtual bool initialize(
      const std::vector<VectorDim>& coordinates,
      const std::vector<std::vector<mpm::Index>>& pointsets) {
    bool status = true;
    // Create points from file
    bool point_status = this->create_points(coordinates);
    if (!point_status) {
      status = false;
      throw std::runtime_error(
          "Addition of points in discontinuity to mesh failed");
    }
    // Create elements from file
    bool element_status = create_elements(pointsets);
    if (!element_status) {
      status = false;
      throw std::runtime_error(
          "Addition of elements in discontinuity to mesh failed");
    }

    bool normal_status = initialize_center_normal();
    if (!normal_status) {
      status = false;
      throw std::runtime_error(
          "initialized the center and normal of the discontunity failed");
    }
    return status;
  };

  //! create elements from file
  virtual bool create_elements(
      const std::vector<std::vector<mpm::Index>>& elements) override;

  // initialize the center and normal of the triangular elements
  bool initialize_center_normal();

  // return the cross product of ab and bc
  VectorDim ThreeCross(const VectorDim& a, const VectorDim& b,
                       const VectorDim& c);

  // return the levelset values of each doordinates
  //! \param[in] the vector of the coordinates
  virtual void compute_levelset(const std::vector<VectorDim>& coordinates,
                                std::vector<double>& phi_list) override;

 protected:
  using mpm::DiscontinuityBase<Tdim>::points_;

  using mpm::DiscontinuityBase<Tdim>::console_;

  using mpm::DiscontinuityBase<Tdim>::numpoint_;

 private:
  // vector of elements
  std::vector<discontinuous_element> elements_;

  // number of elements
  mpm::Index numelement_;
};

}  // namespace mpm
#include "discontinuity_3d.tcc"

#endif  // MPM_HEXAHEDRON_ELEMENT_H_
