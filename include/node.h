#ifndef MPM_NODE_H_
#define MPM_NODE_H_

#include <array>
#include <limits>
#include <mutex>
#include <tuple>
#include <vector>

#include "logger.h"
#include "node_base.h"

namespace mpm {

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
  ~Node() override{};

  //! Delete copy constructor
  Node(const Node<Tdim, Tdof, Tnphases>&) = delete;

  //! Delete assignement operator
  Node& operator=(const Node<Tdim, Tdof, Tnphases>&) = delete;

  //! Initialise nodal properties
  void initialise() override;

  //! Return id of the nodebase
  Index id() const override { return id_; }

  //! Assign coordinates
  //! \param[in] coord Assign coord as coordinates of the nodebase
  void assign_coordinates(const VectorDim& coord) override {
    coordinates_ = coord;
  }

  //! Return coordinates
  //! \retval coordinates_ return coordinates of the nodebase
  VectorDim coordinates() const override { return coordinates_; }

  //! Return degrees of freedom
  unsigned dof() const override { return dof_; }

  //! Assign status
  void assign_status(bool status) override { status_ = status; }

  //! Return status
  bool status() const override { return status_; }

  //! Update mass at the nodes from particle
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] mass Mass from the particles in a cell
  void update_mass(bool update, unsigned phase, double mass) override;

  //! Return mass at a given node for a given phase
  //! \param[in] phase Index corresponding to the phase
  double mass(unsigned phase) const override { return mass_(phase); }

  //! Update volume at the nodes from particle
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] volume Volume from the particles in a cell
  void update_volume(bool update, unsigned phase, double volume) override;

  //! Return volume at a given node for a given phase
  //! \param[in] phase Index corresponding to the phase
  double volume(unsigned phase) const override { return volume_(phase); }

  //! Assign traction force to the node
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] direction Index corresponding to the direction of traction
  //! \param[in] traction Nodal traction in specified direction
  //! \retval status Assignment status
  bool assign_traction_force(unsigned phase, unsigned direction,
                             double traction) override;

  //! Update external force (body force / traction force)
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] force External force from the particles in a cell
  //! \retval status Update status
  bool update_external_force(bool update, unsigned phase,
                             const VectorDim& force) override;

  //! Return external force at a given node for a given phase
  //! \param[in] phase Index corresponding to the phase
  VectorDim external_force(unsigned phase) const override {
    return external_force_.col(phase);
  }

  //! Update internal force (body force / traction force)
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] force Internal force from the particles in a cell
  //! \retval status Update status
  bool update_internal_force(bool update, unsigned phase,
                             const VectorDim& force) override;

  //! Return internal force at a given node for a given phase
  //! \param[in] phase Index corresponding to the phase
  VectorDim internal_force(unsigned phase) const override {
    return internal_force_.col(phase);
  }

  //! Update pressure at the nodes from particle
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] mass_pressure Product of mass x pressure of a particle
  void update_mass_pressure(unsigned phase, double mass_pressure) override;

  //! Assign pressure at the nodes from particle
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] mass_pressure Product of mass x pressure of a particle
  void assign_pressure(unsigned phase, double mass_pressure) override;

  //! Return pressure at a given node for a given phase
  //! \param[in] phase Index corresponding to the phase
  double pressure(unsigned phase) const override { return pressure_(phase); }

  //! Update momentum at the nodes
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] momentum Momentum from the particles in a cell
  //! \retval status Update status
  bool update_momentum(bool update, unsigned phase,
                       const VectorDim& momentum) override;

  //! Return momentum at a given node for a given phase
  //! \param[in] phase Index corresponding to the phase
  VectorDim momentum(unsigned phase) const override {
    return momentum_.col(phase);
  }

  //! Compute velocity from the momentum
  void compute_velocity() override;

  //! Return velocity at a given node for a given phase
  //! \param[in] phase Index corresponding to the phase
  VectorDim velocity(unsigned phase) const override {
    return velocity_.col(phase);
  }

  //! Update nodal acceleration
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] acceleration Acceleration from the particles in a cell
  //! \retval status Update status
  bool update_acceleration(bool update, unsigned phase,
                           const VectorDim& acceleration) override;

  //! Return acceleration at a given node for a given phase
  //! \param[in] phase Index corresponding to the phase
  VectorDim acceleration(unsigned phase) const override {
    return acceleration_.col(phase);
  }

  //! Compute acceleration and velocity
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] dt Timestep in analysis
  bool compute_acceleration_velocity(unsigned phase, double dt) override;

  //! Assign velocity constraint
  //! Directions can take values between 0 and Dim * Nphases
  //! \param[in] dir Direction of velocity constraint
  //! \param[in] velocity Applied velocity constraint
  bool assign_velocity_constraint(unsigned dir, double velocity) override;

  //! Apply velocity constraints
  void apply_velocity_constraints() override;

  //! Assign friction constraint
  //! Directions can take values between 0 and Dim * Nphases
  //! \param[in] dir Direction of friction constraint
  //! \param[in] sign Sign of normal wrt coordinate system for friction
  //! \param[in] friction Applied friction constraint
  bool assign_friction_constraint(unsigned dir, int sign,
                                  double friction) override;

  //! Apply friction constraints
  //! \param[in] dt Time-step
  void apply_friction_constraints(double dt) override;

  //! Assign rotation matrix
  //! \param[in] rotation_matrix Rotation matrix of the node
  void assign_rotation_matrix(
      const Eigen::Matrix<double, Tdim, Tdim>& rotation_matrix) override {
    rotation_matrix_ = rotation_matrix;
    generic_boundary_constraints_ = true;
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
  //! Status
  bool status_{false};
  //! Mass
  Eigen::Matrix<double, 1, Tnphases> mass_;
  //! Volume
  Eigen::Matrix<double, 1, Tnphases> volume_;
  //! External force
  Eigen::Matrix<double, Tdim, Tnphases> external_force_;
  //! Internal force
  Eigen::Matrix<double, Tdim, Tnphases> internal_force_;
  //! Pressure
  Eigen::Matrix<double, 1, Tnphases> pressure_;
  //! Velocity
  Eigen::Matrix<double, Tdim, Tnphases> velocity_;
  //! Momentum
  Eigen::Matrix<double, Tdim, Tnphases> momentum_;
  //! Acceleration
  Eigen::Matrix<double, Tdim, Tnphases> acceleration_;
  //! Velocity constraints
  std::map<unsigned, double> velocity_constraints_;
  //! Rotation matrix for general velocity constraints
  Eigen::Matrix<double, Tdim, Tdim> rotation_matrix_;
  //! A general velocity (non-Cartesian/inclined) constraint is specified at the
  //! node
  bool generic_boundary_constraints_{false};
  //! Frictional constraints
  bool friction_{false};
  std::tuple<unsigned, int, double> friction_constraint_;
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
};  // Node class
}  // namespace mpm

#include "node.tcc"

#endif  // MPM_NODE_H_
