#ifndef MPM_DATA_STORE_H_
#define MPM_DATA_STORE_H_

#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <vector>

#include "Eigen/Dense"

#include "data_types.h"

#include <iostream>

namespace mpm {

//! DataStore class
//! \brief Base class that stores information
template <unsigned Tdim>
class DataStore {
 public:
  //! Constructor with id, number of nodes
  DataStore(unsigned nnodes, unsigned nphases)
      : nnodes_{nnodes}, nphases_{nphases} {
    this->initialise();
  }

  void initialise() {
    internal_forces.clear();
    external_forces.clear();
    // Initialise vector quantities
    for (unsigned i = 0; i < nphases_; ++i) {
      std::vector<Eigen::Matrix<double, Tdim, 1>> phase_internal_forces;
      std::vector<Eigen::Matrix<double, Tdim, 1>> phase_external_forces;
      // Reserve
      phase_internal_forces.reserve(nnodes_);
      phase_external_forces.reserve(nnodes_);
      for (unsigned j = 0; j < nnodes_; ++j) {
        phase_internal_forces.emplace_back(
            Eigen::Matrix<double, Tdim, 1>::Zero());
        phase_external_forces.emplace_back(
            Eigen::Matrix<double, Tdim, 1>::Zero());
      }
      internal_forces.emplace_back(phase_internal_forces);
      external_forces.emplace_back(phase_external_forces);
    }
  }

  //! Internal forces
  std::vector<std::vector<Eigen::Matrix<double, Tdim, 1>>> internal_forces;
  //! External forces
  std::vector<std::vector<Eigen::Matrix<double, Tdim, 1>>> external_forces;

 private:
  // Number of nodes
  Index nnodes_;
  // Number of phases
  unsigned nphases_{1};

};  // Datastore class
}  // namespace mpm
#endif  // MPM_DATA_STORE_H_
