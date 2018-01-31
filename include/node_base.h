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

// NodeBase base class for nodes
//! \brief Base class that stores the information about node_bases
//! \details NodeBase class: id_ and coordinates.
//! \tparam Tdim Dimension
template <unsigned Tdim>
class NodeBase {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Constructor with id and coordinates
  //! \param[in] id Node id
  //! \param[in] coord coordinates of the nodebase
  NodeBase(Index id, const VectorDim& coord) : id_{id} {
    // Check if the dimension is between 1 & 3
    static_assert((Tdim >= 1 && Tdim <= 3), "Invalid global dimension");
    coordinates_ = coord;
  };

  //! Destructor
  virtual ~NodeBase(){};

  //! Delete copy constructor
  NodeBase(const NodeBase<Tdim>&) = delete;

  //! Delete assignement operator
  NodeBase& operator=(const NodeBase<Tdim>&) = delete;

  //! Return id of the nodebase
  Index id() const { return id_; }

  //! Assign coordinates
  //! \param[in] coord Assign coord as coordinates of the nodebase
  void coordinates(const VectorDim& coord) { coordinates_ = coord; }

  //! Return coordinates
  //! \retval coordinates_ return coordinates of the nodebase
  VectorDim coordinates() const { return coordinates_; }

  //! Initialise properties
  virtual void initialise() = 0;

  //! Return degrees of freedom
  virtual unsigned dof() const = 0;

  //! Assign nodal mass
  virtual void assign_mass(double mass) = 0;

  //! Return mass
  virtual double mass() const = 0;

  //! Assign force
  virtual void assign_force(unsigned nphase, const Eigen::VectorXd& force) = 0;

  //! Return force
  virtual Eigen::VectorXd force(unsigned nphase) const = 0;

  //! Assign velocity
  virtual void assign_velocity(unsigned nphase, const Eigen::VectorXd& velocity) = 0;

  //! Return velocity
  virtual Eigen::VectorXd velocity(unsigned nphase) const = 0;

  //! Assign momentum
  virtual void assign_momentum(unsigned nphase, const Eigen::VectorXd& momentum) = 0;

  //! Return momentum
  virtual Eigen::VectorXd momentum(unsigned nphase) const = 0;

  //! Assign acceleration
  virtual void assign_acceleration(unsigned nphase, const Eigen::VectorXd& acceleration) = 0;

  //! Return acceleration
  virtual Eigen::VectorXd acceleration(unsigned nphase) const = 0;

 protected:
  //! nodebase id
  Index id_{std::numeric_limits<Index>::max()};
  //! nodal coordinates
  VectorDim coordinates_;

};  // NodeBase class
}  // namespace mpm

#endif  // MPM_NODE_BASE_H_
