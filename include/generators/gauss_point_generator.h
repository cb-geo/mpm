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
                 const Json& generator);

  //! Generate material points
  std::vector<Eigen::Matrix<double, Tdim, 1>> generate_points() override;

};  // GaussPointGenerator class
}  // namespace mpm

#include "gauss_point_generator.tcc"

#endif  // MPM_GAUSS_POINT_GENERATOR_H_
