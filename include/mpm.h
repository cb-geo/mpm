#ifndef MPM_MPM_H_
#define MPM_MPM_H_

#include <memory>
#include <vector>

#include "mesh.h"

namespace mpm {

//! MPM class
//! \brief MPM class calls solver and algorithm
//! \details MPM class: implicit and explicit MPM
//! \tparam Tdim Dimension
template <unsigned Tdim>
class MPM {
  //! Default constructor
  MPM(){};

 private:
  //! Mesh object
  std::vector<std::unique_ptr<lem::Mesh<Tdim>>> meshes_;
};
}

#endif  // MPM_MPM_H_
