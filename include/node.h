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

  //! Return degrees of freedom
  unsigned dof() const { return dof_; }

  //! Assign nodal mass
  void assign_mass(unsigned nphase, double mass) { mass_(0, nphase) = mass; }

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

 protected:
  //! node id
  using NodeBase<Tdim>::id_;
  //! nodal coordinates
  using NodeBase<Tdim>::coordinates_;

 private:
  //! Degrees of freedom
  unsigned dof_{std::numeric_limits<unsigned>::max()};
  //! Mass solid
  Eigen::Matrix<double, 1, Tnphases> mass_;
  //! Force
  Eigen::Matrix<double, Tdof, Tnphases> force_;
  //! Velocity
  Eigen::Matrix<double, Tdof, Tnphases> velocity_;
  //! Momentum
  Eigen::Matrix<double, Tdof, Tnphases> momentum_;
  //! Acceleration
  Eigen::Matrix<double, Tdof, Tnphases> acceleration_;
};  // Node class
}  // namespace mpm

#include "node.tcc"

#endif  // MPM_NODE_H_
