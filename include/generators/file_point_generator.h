#ifndef MPM_FILE_POINT_GENERATOR_H_
#define MPM_FILE_POINT_GENERATOR_H_

#include "point_generator.h"

// JSON
using Json = nlohmann::json;

namespace mpm {

//! FilePointGenerator class
//! \brief Base class that defines generation of material points
//! \tparam Tdim Dimension
template <unsigned Tdim>
class FilePointGenerator : public PointGenerator<Tdim> {
 public:
  //! Constructor with mesh pointer and generator properties
  FilePointGenerator(const std::shared_ptr<mpm::Mesh<Tdim>>& mesh,
                     const std::shared_ptr<IO>& io, const Json& generator);

  //! Generate material points
  std::vector<Eigen::Matrix<double, Tdim, 1>> generate_points() override;

 private:
  // Mesh object
  using PointGenerator<Tdim>::mesh_;
  // IO object
  using PointGenerator<Tdim>::io_;
  // Generator
  using PointGenerator<Tdim>::generator_;
  // Logger
  using PointGenerator<Tdim>::console_;

};  // namespace mpm
}  // namespace mpm

#include "file_point_generator.tcc"
#endif  // MPM_FILE_POINT_GENERATOR_H_
