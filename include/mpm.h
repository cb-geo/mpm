#ifndef MPM_MPM_H_
#define MPM_MPM_H_

#include <memory>
#include <vector>

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

#include "io.h"
#include "mesh.h"
#include "read_mesh.h"
#include "read_mesh_ascii.h"

namespace mpm {
//! MPM class
//! \brief MPM class calls solver and algorithm
//! \details MPM class: implicit and explicit MPM
//! \tparam Tdim Dimension
template <unsigned Tdim>
class MPM {
 public:
  //! Constructor
  MPM(std::unique_ptr<IO>&& io) : io_(std::move(io)) {
    // Unique id
    uuid_ =
        boost::lexical_cast<std::string>(boost::uuids::random_generator()());
    meshes_.clear();
  }

  //! Initialise mesh
  virtual bool initialise() = 0;

 protected:
  //! A unique id for the analysis
  std::string uuid_;
  //! Time step size
  double dt_{std::numeric_limits<double>::max()};
  //! Number of steps
  mpm::Index nsteps_{std::numeric_limits<mpm::Index>::max()};
  //! Gravity
  Eigen::Matrix<double, Tdim, 1> gravity_;
  //! Mesh object
  std::vector<std::unique_ptr<mpm::Mesh<Tdim>>> meshes_;
  //! A unique ptr to IO object
  std::unique_ptr<mpm::IO> io_;
  //! JSON analysis object
  Json analysis_;
  //! Logger
  std::shared_ptr<spdlog::logger> console_;
};
}  // namespace mpm

#endif  // MPM_MPM_H_
