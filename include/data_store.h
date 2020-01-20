#ifndef MPM_DATA_STORE_H_
#define MPM_DATA_STORE_H_

#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <vector>

#include "Eigen/Dense"

namespace mpm {

//! DataStore class
//! \brief Base class that stores information
template <unsigned Tdim>
class DataStore {
 public:
  //! Constructor with id, number of nodes
  DataStore(unsigned nnodes) {
    internal_forces.reserve(nnodes);
    external_forces.reserve(nnodes);
    for (unsigned i = 0; i < nnodes; ++i) {
      internal_forces.emplace_back(Eigen::Matrix<double, Tdim, 1>::Zero());
      external_forces.emplace_back(Eigen::Matrix<double, Tdim, 1>::Zero());
    }
  }

  void initialise() {
    for (unsigned i = 0; i < internal_forces.size(); ++i) {
      internal_forces[i].setZero();
      external_forces[i].setZero();
    }
  }
  
  //! Internal forces
  std::vector<Eigen::Matrix<double, Tdim, 1>> internal_forces;
  //! External forces
  std::vector<Eigen::Matrix<double, Tdim, 1>> external_forces;

};  // Datastore class
}  // namespace mpm
#endif  // MPM_DATA_STORE_H_
