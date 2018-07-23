#ifndef MPM_NODE_H_
#define MPM_NODE_H_

#include <array>
#include <iostream>
#include <limits>
#include <mutex>
#include <vector>

#include "node_base.h"

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
  //! \param[in] id Node id
  //! \param[in] coord coordinates of the node
  Node(Index id, const VectorDim& coord);

  //! Virtual destructor
  virtual ~Node(){};

  //! Delete copy constructor
  Node(const Node<Tdim, Tdof, Tnphases>&) = delete;

  //! Delete assignement operator
  Node& operator=(const Node<Tdim, Tdof, Tnphases>&) = delete;

  //! Initialise nodal properties
  void initialise();

  //! Return id of the nodebase
  Index id() const { return id_; }

  //! Assign coordinates
  //! \param[in] coord Assign coord as coordinates of the nodebase
  void assign_coordinates(const VectorDim& coord) { coordinates_ = coord; }

  //! Return coordinates
  //! \retval coordinates_ return coordinates of the nodebase
  VectorDim coordinates() const { return coordinates_; }

  //! Return degrees of freedom
  unsigned dof() const { return dof_; }

  //! Update mass at the nodes from particle
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] nphase Index corresponding to the phase
  //! \param[in] mass Mass from the particles in a cell
  void update_mass(bool update, unsigned nphase, double mass);

  //! Return mass at a given node for a given phase
  //! \param[in] nphase Index corresponding to the phase
  double mass(unsigned nphase) const { return mass_(nphase); }

  //! Update volume at the nodes from particle
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] nphase Index corresponding to the phase
  //! \param[in] volume Volume from the particles in a cell
  void update_volume(bool update, unsigned nphase, double volume);

  //! Return volume at a given node for a given phase
  //! \param[in] nphase Index corresponding to the phase
  double volume(unsigned nphase) const { return volume_(nphase); }

  //! Update external force (body force / traction force)
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] nphase Index corresponding to the phase
  //! \param[in] force External force from the particles in a cell
  //! \retval status Update status
  bool update_external_force(bool update, unsigned nphase,
                             const Eigen::VectorXd& force);

  //! Return external force at a given node for a given phase
  //! \param[in] nphase Index corresponding to the phase
  Eigen::VectorXd external_force(unsigned nphase) const {
    return external_force_.col(nphase);
  }

  //! Update internal force (body force / traction force)
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] nphase Index corresponding to the phase
  //! \param[in] force Internal force from the particles in a cell
  //! \retval status Update status
  bool update_internal_force(bool update, unsigned nphase,
                             const Eigen::VectorXd& force);

  //! Return internal force at a given node for a given phase
  //! \param[in] nphase Index corresponding to the phase
  Eigen::VectorXd internal_force(unsigned nphase) const {
    return internal_force_.col(nphase);
  }

  //! Update momentum at the nodes
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] nphase Index corresponding to the phase
  //! \param[in] momentum Momentum from the particles in a cell
  //! \retval status Update status
  bool update_momentum(bool update, unsigned nphase,
                       const Eigen::VectorXd& momentum);

  //! Return momentum at a given node for a given phase
  //! \param[in] nphase Index corresponding to the phase
  Eigen::VectorXd momentum(unsigned nphase) const {
    return momentum_.col(nphase);
  }

  //! Compute velocity from the momentum
  void compute_velocity();

  //! Return velocity at a given node for a given phase
  //! \param[in] nphase Index corresponding to the phase
  Eigen::VectorXd velocity(unsigned nphase) const {
    return velocity_.col(nphase);
  }

  //! Update nodal acceleration
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] nphase Index corresponding to the phase
  //! \param[in] acceleration Acceleration from the particles in a cell
  //! \retval status Update status
  bool update_acceleration(bool update, unsigned nphase,
                           const Eigen::VectorXd& acceleration);

  //! Return acceleration at a given node for a given phase
  //! \param[in] nphase Index corresponding to the phase
  Eigen::VectorXd acceleration(unsigned nphase) const {
    return acceleration_.col(nphase);
  }

  // TODO: REMOVE
  void stats() {
    std::string out = "Node: " + std::to_string(id_) +
                      "mass: " + std::to_string(mass_(0)) + "\n";
    std::cout << out;
  }

 private:
  //! Mutex
  std::mutex node_mutex_;
  //! nodebase id
  Index id_{std::numeric_limits<Index>::max()};
  //! nodal coordinates
  VectorDim coordinates_;
  //! Degrees of freedom
  unsigned dof_{std::numeric_limits<unsigned>::max()};
  //! Mass
  Eigen::Matrix<double, 1, Tnphases> mass_;
  //! Volume
  Eigen::Matrix<double, 1, Tnphases> volume_;
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
