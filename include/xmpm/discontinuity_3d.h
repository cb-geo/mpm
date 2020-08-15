#ifndef MPM_DISCONTINUITY_3D_H_
#define MPM_DISCONTINUITY_3D_H_

#include "discontinuity_base.h"

//! MPM namespace
namespace mpm {

template <unsigned Tdim>
class Discontinuity3D : public DiscontinuityBase<Tdim> {

 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Constructor with id
  //! \param[in] discontinuity_properties discontinuity properties
  Discontinuity3D(unsigned id, const Json& discontinuity_props);

  // initialization
  virtual bool initialize(
      const std::vector<VectorDim>& coordinates,
      const std::vector<std::vector<mpm::Index>>& pointsets);

  //! create elements from file
  virtual bool create_elements(
      const std::vector<std::vector<mpm::Index>>& elements) override;

  // initialize the center and normal of the triangular elements
  bool initialize_center_normal();

  // return the cross product of ab and bc
  VectorDim three_cross_product(const VectorDim& a, const VectorDim& b,
                       const VectorDim& c);

  // return the levelset values of each doordinates
  //! \param[in] the vector of the coordinates
 void compute_levelset(const std::vector<VectorDim>& coordinates,
                                std::vector<double>& phi_list) override;
  // return the normal vectors of given coordinates
  //! \param[in] the coordinates
  void compute_normal(
   const VectorDim& coordinates, VectorDim& normal_vector) override;

 protected:
  using mpm::DiscontinuityBase<Tdim>::points_;

  using mpm::DiscontinuityBase<Tdim>::console_;

  using mpm::DiscontinuityBase<Tdim>::numpoint_;

  using mpm::DiscontinuityBase<Tdim>::friction_coef_;

 private:
  // vector of elements
  std::vector<discontinuous_element<Tdim>> elements_;

  // number of elements //delete
  mpm::Index numelement_;
};

}  // namespace mpm
#include "discontinuity_3d.tcc"

#endif  // MPM_HEXAHEDRON_ELEMENT_H_
