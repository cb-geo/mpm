#ifndef MPM_NODE_H_
#define MPM_NODE_H_

#include "logger.h"
#include "mutex.h"
#include "nodal_properties.h"
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
  void initialise() noexcept override;

  //! Return id of the nodebase
  Index id() const override { return id_; }

  //! Initialise shared pointer to nodal properties pool
  //! \param[in] prop_id Property id in the nodal property pool
  //! \param[in] nodal_properties Shared pointer to nodal properties pool
  void initialise_property_handle(
      unsigned prop_id,
      std::shared_ptr<mpm::NodalProperties> property_handle) noexcept override;

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
  void update_mass(bool update, unsigned phase, double mass) noexcept override;

  //! Return mass at a given node for a given phase
  //! \param[in] phase Index corresponding to the phase
  double mass(unsigned phase) const override { return mass_(phase); }

  //! Update volume at the nodes from particle
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] volume Volume from the particles in a cell
  void update_volume(bool update, unsigned phase,
                     double volume) noexcept override;

  //! Return volume at a given node for a given phase
  //! \param[in] phase Index corresponding to the phase
  double volume(unsigned phase) const override { return volume_(phase); }

  //! Assign concentrated force to the node
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] direction Index corresponding to the direction of traction
  //! \param[in] force Nodal concentrated force in specified direction
  //! \param[in] function math function
  //! \retval status Assignment status
  bool assign_concentrated_force(
      unsigned phase, unsigned direction, double force,
      const std::shared_ptr<FunctionBase>& function) override;

  //! Apply concentrated force to external force
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] current time
  void apply_concentrated_force(unsigned phase, double current_time) override;

  //! Update external force (body force / traction force)
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] force External force from the particles in a cell
  void update_external_force(bool update, unsigned phase,
                             const VectorDim& force) noexcept override;

  //! Return external force at a given node for a given phase
  //! \param[in] phase Index corresponding to the phase
  VectorDim external_force(unsigned phase) const override {
    return external_force_.col(phase);
  }

  //! Update internal force (body force / traction force)
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] force Internal force from the particles in a cell
  void update_internal_force(bool update, unsigned phase,
                             const VectorDim& force) noexcept override;

  //! Return internal force at a given node for a given phase
  //! \param[in] phase Index corresponding to the phase
  VectorDim internal_force(unsigned phase) const override {
    return internal_force_.col(phase);
  }

  //! Update pressure at the nodes from particle
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] mass_pressure Product of mass x pressure of a particle
  void update_mass_pressure(unsigned phase,
                            double mass_pressure) noexcept override;

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
  void update_momentum(bool update, unsigned phase,
                       const VectorDim& momentum) noexcept override;

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
  void update_acceleration(bool update, unsigned phase,
                           const VectorDim& acceleration) noexcept override;

  //! Return acceleration at a given node for a given phase
  //! \param[in] phase Index corresponding to the phase
  VectorDim acceleration(unsigned phase) const override {
    return acceleration_.col(phase);
  }

  //! Compute acceleration and velocity
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] dt Timestep in analysis
  bool compute_acceleration_velocity(unsigned phase,
                                     double dt) noexcept override;

  //! Compute acceleration and velocity with cundall damping factor
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] dt Timestep in analysis
  //! \param[in] damping_factor Damping factor
  bool compute_acceleration_velocity_cundall(
      unsigned phase, double dt, double damping_factor) noexcept override;

  //! Assign velocity constraint
  //! Directions can take values between 0 and Dim * Nphases
  //! \param[in] dir Direction of velocity constraint
  //! \param[in] velocity Applied velocity constraint
  bool assign_velocity_constraint(unsigned dir, double velocity) override;

  //! Apply velocity constraints
  void apply_velocity_constraints() override;

  //! Assign absorbing constraint
  //! \param[in] pwave_v P-wave velocity
  //! \param[in] swave_v S-wave velocity
  bool assign_absorbing_constraint(unsigned dir, double pwave_v,
                                   double swave_v) override;

  //! Apply absorbing constraint
  void apply_absorbing_constraint() override;

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

  //! Add material id from material points to list of materials in materials_
  //! \param[in] id Material id to be stored at the node
  void append_material_id(unsigned id) override;

  //! Return material ids in node
  std::set<unsigned> material_ids() const override { return material_ids_; }

  //! Assign MPI rank to node
  //! \param[in] rank MPI Rank of the node
  bool mpi_rank(unsigned rank) override;

  //! Assign MPI rank to node
  //! \param[in] rank MPI Rank of the node
  std::set<unsigned> mpi_ranks() const override { return mpi_ranks_; }

  //! Clear MPI rank
  void clear_mpi_ranks() override { mpi_ranks_.clear(); }

  //! Return ghost id
  Index ghost_id() const override { return ghost_id_; }

  //! Set ghost id
  void ghost_id(Index gid) override { ghost_id_ = gid; }

  //! Update nodal property at the nodes from particle
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] property Property name
  //! \param[in] property_value Property quantity from the particles in the cell
  //! \param[in] mat_id Id of the material within the property data
  //! \param[in] nprops Dimension of property (1 if scalar, Tdim if vector)
  void update_property(bool update, const std::string& property,
                       const Eigen::MatrixXd& property_value, unsigned mat_id,
                       unsigned nprops) noexcept override;

  //! Compute multimaterial change in momentum
  void compute_multimaterial_change_in_momentum() override;

  //! Compute multimaterial separation vector
  void compute_multimaterial_separation_vector() override;

  //! Compute multimaterial normal unit vector
  void compute_multimaterial_normal_unit_vector() override;

 private:
  //! Mutex
  SpinMutex node_mutex_;
  //! nodebase id
  Index id_{std::numeric_limits<Index>::max()};
  //! nodal property id
  unsigned prop_id_{std::numeric_limits<unsigned>::max()};
  //! shared ghost id
  Index ghost_id_{std::numeric_limits<Index>::max()};
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
  //! Displacement
  Eigen::Matrix<double, Tdim, 1> contact_displacement_;
  //! Velocity
  Eigen::Matrix<double, Tdim, Tnphases> velocity_;
  //! Momentum
  Eigen::Matrix<double, Tdim, Tnphases> momentum_;
  //! Acceleration
  Eigen::Matrix<double, Tdim, Tnphases> acceleration_;
  //! Velocity constraints
  std::map<unsigned, double> velocity_constraints_;
  //! Absorbing Constraints
  std::tuple<unsigned, double, double> absorbing_constraints_;
  //! Absorbing Traction
  Eigen::Matrix<double, Tdim, Tnphases> absorbing_traction_;
  //! Rotation matrix for general velocity constraints
  Eigen::Matrix<double, Tdim, Tdim> rotation_matrix_;
  //! Material ids whose information was passed to this node
  std::set<unsigned> material_ids_;
  //! A general velocity (non-Cartesian/inclined) constraint is specified at the
  //! node
  bool generic_boundary_constraints_{false};
  //! Frictional constraints
  bool friction_{false};
  std::tuple<unsigned, int, double> friction_constraint_;
  //! Concentrated force
  Eigen::Matrix<double, Tdim, Tnphases> concentrated_force_;
  //! Mathematical function for force
  std::shared_ptr<FunctionBase> force_function_{nullptr};
  //! Nodal property pool
  std::shared_ptr<mpm::NodalProperties> property_handle_{nullptr};
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
  //! MPI ranks
  std::set<unsigned> mpi_ranks_;
};  // Node class
}  // namespace mpm

#include "node.tcc"

#endif  // MPM_NODE_H_
