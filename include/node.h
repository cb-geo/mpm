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

//! Global index type for the node
using Index = unsigned long long;

// Node base class
//! \brief Base class that stores the information about nodes
//! \details Node class: id_ and coordinates.
//! \tparam Tdim Dimension
//! \tparam Tdof Degrees of Freedom
template <unsigned Tdim, unsigned Tdof>
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
  Node(const Node<Tdim, Tdof>&) = delete;

  //! Delete assignement operator
  Node& operator=(const Node<Tdim, Tdof>&) = delete;

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
  //! \param[in] mass Mass from the particles in a cell
  void update_mass(bool update, double mass) override;

  //! Return mass at a given node
  double mass() const override { return mass_; }

  //! Update volume at the nodes from particle
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] volume Volume from the particles in a cell
  void update_volume(bool update, double volume) override;

  //! Return volume at a given node
  double volume() const override { return volume_; }

  //! Assign traction force to the node
  //! \param[in] direction Index corresponding to the direction of traction
  //! \param[in] traction Nodal traction in specified direction
  //! \retval status Assignment status
  bool assign_traction_force(unsigned direction, double traction) override;

  //! Update external force (body force / traction force)
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] force External force from the particles in a cell
  //! \retval status Update status
  bool update_external_force(bool update, const VectorDim& force) override;

  //! Return external force at a given node
  VectorDim external_force() const override {
    return external_force_;
  }

  //! Update internal force (body force / traction force)
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] force Internal force from the particles in a cell
  //! \retval status Update status
  bool update_internal_force(bool update, const VectorDim& force) override;

  //! Return internal force at a given node
  VectorDim internal_force() const override { return internal_force_; }

  //! Update pressure at the nodes from particle
  //! \param[in] mass_pressure Product of mass x pressure of a particle
  void update_mass_pressure(double mass_pressure) override;

  //! Assign pressure at the nodes from particle
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] mass_pressure Product of mass x pressure of a particle
  void assign_pressure(double mass_pressure) override;

  //! Return pressure at a given node 
  double pressure() const override { return pressure_; }

  //! Update momentum at the nodes
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] momentum Momentum from the particles in a cell
  //! \retval status Update status
  bool update_momentum(bool update, const VectorDim& momentum) override;

  //! Return momentum at a given node
  VectorDim momentum() const override { return momentum_; }

  //! Compute velocity from the momentum
  void compute_velocity() override;

  //! Return velocity at a given node
  VectorDim velocity() const override { return velocity_; }

  //! Update nodal acceleration
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] acceleration Acceleration from the particles in a cell
  //! \retval status Update status
  bool update_acceleration(bool update, const VectorDim& acceleration) override;

  //! Return acceleration at a given node
  VectorDim acceleration() const override { return acceleration_; }

  //! Compute acceleration and velocity
  //! \param[in] dt Timestep in analysis
  bool compute_acceleration_velocity(double dt) override;

  //! Assign velocity constraint
  //! Directions can take values between 0 and Dim
  //! \param[in] dir Direction of velocity constraint
  //! \param[in] velocity Applied velocity constraint
  bool assign_velocity_constraint(unsigned dir, double velocity) override;

  //! Apply velocity constraints
  void apply_velocity_constraints() override;

  //! Assign friction constraint
  //! Directions can take values between 0 and Dim
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
  double mass_;
  //! Volume
  double volume_;
  //! External force
  Eigen::Matrix<double, Tdim, 1> external_force_;
  //! Internal force
  Eigen::Matrix<double, Tdim, 1> internal_force_;
  //! Pressure
  double pressure_;
  //! Velocity
  Eigen::Matrix<double, Tdim, 1> velocity_;
  //! Momentum
  Eigen::Matrix<double, Tdim, 1> momentum_;
  //! Acceleration
  Eigen::Matrix<double, Tdim, 1> acceleration_;
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
