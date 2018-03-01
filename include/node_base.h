#ifndef MPM_NODE_BASE_H_
#define MPM_NODE_BASE_H_

#include <array>
#include <iostream>
#include <limits>
#include <vector>

#include "serialize.h"

namespace mpm {

//! Global index type for the node_base
using Index = unsigned long long;

//! NodeBase base class for nodes
//! \brief Base class that stores the information about node_bases
//! \details NodeBase class: id_ and coordinates.
//! \tparam Tdim Dimension
template <unsigned Tdim>
class NodeBase {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Default constructor
  NodeBase() = default;

  //! Destructor
  virtual ~NodeBase(){};

  //! Delete copy constructor
  NodeBase(const NodeBase<Tdim>&) = delete;

  //! Delete assignement operator
  NodeBase& operator=(const NodeBase<Tdim>&) = delete;

  //! Return id of the nodebase
  virtual Index id() const = 0;

  //! Assign coordinates
  virtual void coordinates(const VectorDim& coord) = 0;

  //! Return coordinates
  //! \retval coordinates_ return coordinates of the nodebase
  virtual VectorDim coordinates() const = 0;

  //! Initialise properties
  virtual void initialise() = 0;

  //! Return degrees of freedom
  virtual unsigned dof() const = 0;

  //! Assign nodal mass
  virtual void update_mass(bool update, unsigned nphase, double mass) = 0;

  //! Return mass
  virtual double mass(unsigned nphase) const = 0;

  //! Update external force (body force / traction force)
  virtual void update_external_force(bool update, unsigned nphase,
                                     const Eigen::VectorXd& force) = 0;

  //! Return external force
  virtual Eigen::VectorXd external_force(unsigned nphase) const = 0;

  //! Update internal force (body force / traction force)
  virtual void update_internal_force(bool update, unsigned nphase,
                                     const Eigen::VectorXd& force) = 0;

  //! Return internal force
  virtual Eigen::VectorXd internal_force(unsigned nphase) const = 0;

  //! Update momentum
  virtual void update_momentum(bool update, unsigned nphase,
                               const Eigen::VectorXd& momentum) = 0;

  //! Return momentum
  virtual Eigen::VectorXd momentum(unsigned nphase) const = 0;

  //! Assign acceleration
  virtual void assign_acceleration(unsigned nphase,
                                   const Eigen::VectorXd& acceleration) = 0;

  //! Return acceleration
  virtual Eigen::VectorXd acceleration(unsigned nphase) const = 0;

  //! Compute velocity from the momentum
  virtual void compute_velocity() = 0;

  //! Assign velocity
  virtual void assign_velocity(unsigned nphase,
                               const Eigen::VectorXd& velocity) = 0;

  //! Return velocity
  virtual Eigen::VectorXd velocity(unsigned nphase) const = 0;

};  // NodeBase class
}  // namespace mpm

#endif  // MPM_NODE_BASE_H_
