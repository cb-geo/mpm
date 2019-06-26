#ifndef MPM_FILE_POINT_GENERATOR_H_
#define MPM_FILE_POINT_GENERATOR_H_


#include "point_generator.h"

namespace mpm {

//! FilePointGenerator class
//! \brief Base class that defines generation of material points
//! \tparam Tdim Dimension
template <unsigned Tdim>
class FilePointGenerator : public PointGenerator<Tdim> {
 public:
  //! Constructor with mesh pointer and generator properties
  FilePointGenerator(const std::shared_ptr<mpm::Mesh<Tdim>>& mesh,
                 const Json& generator);

  //! Generate material points
  std::vector<Eigen::Matrix<double, Tdim, 1>> generate_points() override;

};  // FilePointGenerator class
}  // namespace mpm

#include "file_point_generator.tcc"

#endif  // MPM_FILE_POINT_GENERATOR_H_
