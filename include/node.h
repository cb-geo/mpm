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

// Node class
//! \brief Base class that stores the information about nodes
//! \details Node class: id_ and coordinates.
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Node : public NodeBase<Tdim> {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Constructor with id, coordinates and dof
  Node(Index id, const VectorDim& coord, unsigned dof);

  //! Destructor
  virtual ~Node(){};

  //! Delete copy constructor
  Node(const Node<Tdim>&) = delete;

  //! Delete assignement operator
  Node& operator=(const Node<Tdim>&) = delete;

  //! Initialise properties
  void initialise();

  //! Return degrees of freedom
  unsigned dof() const { return dof_; }

  //! Assign nodal mass
  void assign_mass(double mass) { mass_ = mass; }

  //! Return mass
  double mass() const { return mass_; }

  //! Assign force
  void assign_force(const Eigen::VectorXd& force);

  //! Return force
  Eigen::VectorXd force() { return force_; }

  //! Assign velocity
  void assign_velocity(const Eigen::VectorXd& velocity);

  //! Return velocity
  Eigen::VectorXd velocity() { return velocity_; }

  //! Assign momentum
  void assign_momentum(const Eigen::VectorXd& momentum);

  //! Return momentum
  Eigen::VectorXd momentum() { return momentum_; }

  //! Assign acceleration
  void assign_acceleration(const Eigen::VectorXd& acceleration);

  //! Return acceleration
  Eigen::VectorXd acceleration() { return acceleration_; }

 protected:
  //! node id
  using NodeBase<Tdim>::id_;

  //! nodal coordinates
  using NodeBase<Tdim>::coordinates_;

 private:
  //! Degrees of freedom
  unsigned dof_{std::numeric_limits<unsigned>::max()};
  //! Mass solid
  double mass_{std::numeric_limits<double>::max()};
  //! Force
  Eigen::VectorXd force_;
  //! Velocity
  Eigen::VectorXd velocity_;
  //! Momentum
  Eigen::VectorXd momentum_;
  //! Acceleration
  Eigen::VectorXd acceleration_;
};  // Node class
}  // mpm namespace

#include "node.tcc"

#endif  // MPM_NODE_H_
