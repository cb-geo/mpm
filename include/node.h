#ifndef MPM_NODE_H_
#define MPM_NODE_H_

#include <array>
#include <iostream>
#include <limits>
#include <vector>

#include "node_base.h"
#include "serialize.h"

namespace mpm {

//! Global index type for the node
using Index = unsigned long long;

// Node base class
//! \brief Base class that stores the information about nodes
//! \details Node class: id_ and coordinates.
//! \tparam Tdim Dimension
//! \tparam Tdof Degrees of Freedom
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
class Node : public NodeBase<Tdim> {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Constructor with id, coordinates and dof
  Node(Index id, const VectorDim& coord);

  //! Destructor
  virtual ~Node(){};

  //! Delete copy constructor
  Node(const Node<Tdim, Tdof, Tnphases>&) = delete;

  //! Delete assignement operator
  Node& operator=(const Node<Tdim, Tdof, Tnphases>&) = delete;

  //! Initialise properties
  void initialise();

  //! Return id of the nodebase
  Index id() const { return id_; }

  //! Assign coordinates
  //! \param[in] coord Assign coord as coordinates of the nodebase
  void coordinates(const VectorDim& coord) { coordinates_ = coord; }

  //! Return coordinates
  //! \retval coordinates_ return coordinates of the nodebase
  VectorDim coordinates() const { return coordinates_; }

  //! Return degrees of freedom
  unsigned dof() const { return dof_; }

  //! Update nodal mass
  void update_mass(bool update, unsigned nphase, double mass);

  //! Return mass
  double mass(unsigned nphase) const { return mass_(0, nphase); }

  //! Update external force (body force / traction force)
  void update_external_force(bool update, unsigned nphase, const Eigen::VectorXd& force);

  //! Return external force
  Eigen::VectorXd external_force(unsigned nphase) const { return external_force_.col(nphase); }

  //! Update internal force (body force / traction force)
  void update_internal_force(bool update, unsigned nphase,
                             const Eigen::VectorXd& force);

  //! Return internal force
  Eigen::VectorXd internal_force(unsigned nphase) const { return internal_force_.col(nphase); }

  //! Update momentum
  void update_momentum(bool update, unsigned nphase, const Eigen::VectorXd& momentum);

  //! Return momentum
  Eigen::VectorXd momentum(unsigned nphase) const {
    return momentum_.col(nphase);
  }

  //! Assign acceleration
  void assign_acceleration(unsigned nphase,
                           const Eigen::VectorXd& acceleration);

  //! Return acceleration
  Eigen::VectorXd acceleration(unsigned nphase) const {
    return acceleration_.col(nphase);
  }

  //! Compute velocity from the momentum
  void compute_velocity();

  //! Assign velocity
  void assign_velocity(unsigned nphase, const Eigen::VectorXd& velocity);

  //! Return velocity
  Eigen::VectorXd velocity(unsigned nphase) const {
    return velocity_.col(nphase);
  }

 private:
  //! nodebase id
  Index id_{std::numeric_limits<Index>::max()};
  //! nodal coordinates
  VectorDim coordinates_;
  //! Degrees of freedom
  unsigned dof_{std::numeric_limits<unsigned>::max()};
  //! Mass solid
  Eigen::Matrix<double, 1, Tnphases> mass_;
  //! External force
  Eigen::Matrix<double, Tdim, Tnphases> external_force_;
  //! Internal force
  Eigen::Matrix<double, Tdim, Tnphases> internal_force_;
  //! Velocity
  Eigen::Matrix<double, Tdim, Tnphases> velocity_;
  //! Momentum
  Eigen::Matrix<double, Tdim, Tnphases> momentum_;
  //! Acceleration
  Eigen::Matrix<double, Tdim, Tnphases> acceleration_;
};  // Node class
}  // namespace mpm

#include "node.tcc"

#endif  // MPM_NODE_H_
