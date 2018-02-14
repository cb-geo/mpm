#ifndef MPM_MPM_H_
#define MPM_MPM_H_

#include <memory>
#include <vector>

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

#include "mesh.h"
#include "container.h"
#include "particle.h"

namespace mpm {

//! MPM class
//! \brief MPM class calls solver and algorithm
//! \details MPM class: implicit and explicit MPM
//! \tparam Tdim Dimension
template <unsigned Tdim>
class MPM {
 public:
  //! Constructor
  MPM() {
    // Unique id
    analysis_uuid_ =
        boost::lexical_cast<std::string>(boost::uuids::random_generator()());
    meshes_.clear();
  };

  //! Map particle mass to nodes
// particles.for_each(std::bind(

 private:
  // A unique id for the analysis
  std::string analysis_uuid_;

  //! Mesh object
  std::vector<std::unique_ptr<mpm::Mesh<Tdim>>> meshes_;
  //! Particle container
  mpm::Container<mpm::ParticleBase<Tdim>> particles_;
  //! Nodes container
  mpm::Container<mpm::NodeBase<Tdim>> nodes_;
  //! Cells container
  mpm::Container<mpm::Cell<Tdim>> cells_;
};
}  // namespace mpm

#endif  // MPM_MPM_H_
