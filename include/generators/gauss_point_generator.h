#ifndef MPM_GAUSS_POINT_GENERATOR_H_
#define MPM_GAUSS_POINT_GENERATOR_H_

#include "point_generator.h"

namespace mpm {

//! GaussPointGenerator class
//! \brief Base class that defines generation of material points
//! \tparam Tdim Dimension
template <unsigned Tdim>
class GaussPointGenerator : public PointGenerator<Tdim> {
 public:
  //! Constructor with mesh pointer and generator properties
  GaussPointGenerator(const std::shared_ptr<mpm::Mesh<Tdim>>& mesh,
                      const std::shared_ptr<mpm::IO>& io,
                      const Json& generator);

  //! Generate material points
  std::vector<Eigen::Matrix<double, Tdim, 1>> generate_points() override;

 private:
  // Mesh object
  using PointGenerator<Tdim>::mesh_;
  // IO object
  using PointGenerator<Tdim>::io_;
  // Generator
  using PointGenerator<Tdim>::generator_;

};  // GaussPointGenerator class
}  // namespace mpm

#include "gauss_point_generator.tcc"

#endif  // MPM_GAUSS_POINT_GENERATOR_H_
