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

  //! Update mass at the nodes from particle
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] nphase Index corresponding to the phase
  //! \param[in] mass Mass from the particles in a cell
  virtual void update_mass(bool update, unsigned nphase, double mass) = 0;

  //! Return mass at a given node for a given phase
  virtual double mass(unsigned nphase) const = 0;

  //! Update volume at the nodes from particle
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] nphase Index corresponding to the phase
  //! \param[in] mass Mass from the particles in a cell
  virtual void update_volume(bool update, unsigned nphase, double volume) = 0;

  //! Return volume at a given node for a given phase
  virtual double volume(unsigned nphase) const = 0;

  //! Update external force (body force / traction force)
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] nphase Index corresponding to the phase
  //! \param[in] force External force from the particles in a cell
  virtual void update_external_force(bool update, unsigned nphase,
                                     const Eigen::VectorXd& force) = 0;

  //! Return external force
  //! \param[in] nphase Index corresponding to the phase
  virtual Eigen::VectorXd external_force(unsigned nphase) const = 0;

  //! Update internal force (body force / traction force)
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] nphase Index corresponding to the phase
  //! \param[in] force Internal force from the particles in a cell
  virtual void update_internal_force(bool update, unsigned nphase,
                                     const Eigen::VectorXd& force) = 0;

  //! Return internal force
  //! \param[in] nphase Index corresponding to the phase
  virtual Eigen::VectorXd internal_force(unsigned nphase) const = 0;

  //! Update nodal momentum
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] nphase Index corresponding to the phase
  //! \param[in] momentum Momentum from the particles in a cell
  virtual void update_momentum(bool update, unsigned nphase,
                               const Eigen::VectorXd& momentum) = 0;

  //! Return momentum
  //! \param[in] nphase Index corresponding to the phase
  virtual Eigen::VectorXd momentum(unsigned nphase) const = 0;

  //! Compute velocity from the momentum
  virtual void compute_velocity() = 0;

  //! Return velocity
  //! \param[in] nphase Index corresponding to the phase
  virtual Eigen::VectorXd velocity(unsigned nphase) const = 0;

  //! Update nodal acceleration
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] nphase Index corresponding to the phase
  //! \param[in] acceleration Acceleration from the particles in a cell
  virtual void update_acceleration(bool update, unsigned nphase,
                                   const Eigen::VectorXd& acceleration) = 0;

  //! Return acceleration
  //! \param[in] nphase Index corresponding to the phase
  virtual Eigen::VectorXd acceleration(unsigned nphase) const = 0;

};  // NodeBase class
}  // namespace mpm

#endif  // MPM_NODE_BASE_H_
