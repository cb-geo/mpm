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

  //! Assign nodal mass
  void assign_mass(unsigned nphase, double mass) { mass_(0, nphase) = mass; }

  //! Update nodal mass
  void update_mass(const Eigen::VectorXd& mass) { mass_ += mass; }

  //! Return mass
  double mass(unsigned nphase) const { return mass_(0, nphase); }

  //! Assign force
  void assign_force(unsigned nphase, const Eigen::VectorXd& force);

  //! Return force
  Eigen::VectorXd force(unsigned nphase) const { return force_.col(nphase); }

  //! Assign velocity
  void assign_velocity(unsigned nphase, const Eigen::VectorXd& velocity);

  //! Return velocity
  Eigen::VectorXd velocity(unsigned nphase) const {
    return velocity_.col(nphase);
  }

  //! Assign momentum
  void assign_momentum(unsigned nphase, const Eigen::VectorXd& momentum);

  //! Update nodal momentum
  void update_momentum(const Eigen::MatrixXd& momentum) {
    momentum_ += momentum;
  }

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

  //! Update external force (body force / traction force)
  void update_external_force(const Eigen::MatrixXd& external_force) {
    ext_force_ += external_force;
  }

  //! Update internal force
  void update_internal_force(const Eigen::MatrixXd& internal_force) {
    int_force_ += internal_force;
  }

  //! Compute velocity from the momentum
  void compute_velocity();

 private:
  //! nodebase id
  Index id_{std::numeric_limits<Index>::max()};
  //! nodal coordinates
  VectorDim coordinates_;
  //! Degrees of freedom
  unsigned dof_{std::numeric_limits<unsigned>::max()};
  //! Mass solid
  Eigen::Matrix<double, 1, Tnphases> mass_;
  //! Force
  Eigen::Matrix<double, Tdof, Tnphases> force_;
  //! External force
  Eigen::Matrix<double, Tdim, Tnphases> ext_force_;
  //! Internal force
  Eigen::Matrix<double, Tdim, Tnphases> int_force_;
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
