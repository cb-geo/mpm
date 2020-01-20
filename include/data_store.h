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
class DataStore {
 public:
  //! Constructor with id, number of nodes
  DataStore(unsigned nnodes, unsigned dim, unsigned nphases) {
    internal_forces.resize(nnodes * dim);
    external_forces.resize(nnodes * dim);
  }

  void initialise() {
    internal_forces.setZero();
    external_forces.setZero();
  }
  
  //! Internal forces
  Eigen::VectorXd internal_forces;
  //! External forces
  Eigen::VectorXd external_forces;

};  // Datastore class
}  // namespace mpm
#endif  // MPM_DATA_STORE_H_
