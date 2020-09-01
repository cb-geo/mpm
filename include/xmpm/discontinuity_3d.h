#ifndef MPM_DISCONTINUITY_3D_H_
#define MPM_DISCONTINUITY_3D_H_

#include "discontinuity_base.h"

//! MPM namespace
namespace mpm {
//! Discontinuity3D class derived from DiscontinuityBase class for 3D
template <unsigned Tdim>
class Discontinuity3D : public DiscontinuityBase<Tdim> {

 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Constructor with id
  //! \param[in] discontinuity_props discontinuity properties
  Discontinuity3D(unsigned id, const Json& discontinuity_props);

  //! initialization
  //! \param[in] the coordinates of all points
  //! \param[in] the point index of each surface
  virtual bool initialize(const std::vector<VectorDim>& points,
                          const std::vector<std::vector<mpm::Index>>& surfs);
  //! create elements from file
  //! \param[in] surfs the point index list of each surface
  virtual bool create_surfaces(
      const std::vector<std::vector<mpm::Index>>& surfs) override;

  //! initialize the center and normal vector of each surface
  bool initialize_center_normal();

  //! return the cross product of ab and bc
  //! \param[in] a,b,c coordinates of three points
  VectorDim three_cross_product(const VectorDim& a, const VectorDim& b,
                                const VectorDim& c);

  // return the levelset values of each coordinates
  //! \param[in] coordinates coordinates
  //! \param[in] phi_list the reference of phi for all coordinates
  void compute_levelset(const VectorDim& coordinates,
                        double& phi_particle) override;

  //! compute the normal vectors of coordinates
  //! \param[in] coordinates The coordinates
  //! \param[in] normal vector the normal vector of the given coordinates
  void compute_normal(const VectorDim& coordinates,
                      VectorDim& normal_vector) override;

  //! Assign point friction coefficient
  void assign_point_friction_coef() noexcept override;

 protected:
  //! vector of points
  using mpm::DiscontinuityBase<Tdim>::points_;
  //! Logger
  using mpm::DiscontinuityBase<Tdim>::console_;
  //! friction coefficient
  using mpm::DiscontinuityBase<Tdim>::friction_coef_;

 private:
  // vector of surfaces
  std::vector<discontinuity_surface<Tdim>> surfaces_;
};

}  // namespace mpm
#include "discontinuity_3d.tcc"

#endif  // MPM_HEXAHEDRON_ELEMENT_H_
