#ifndef MPM_NODE_BASE_H_
#define MPM_NODE_BASE_H_

#include <array>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <tuple>
#include <vector>

#include <Eigen/Dense>
#include <flat/flat_map.hpp>

#include "data_types.h"
#include "function_base.h"
#include "mpm_properties.h"
#include "nodal_properties.h"

namespace mpm {

//! NodeBase base class for nodes
//! \brief Base class that stores the information about node_bases
//! \details NodeBase class: id_ and coordinates.
//! \tparam Tdim Dimension
template <unsigned Tdim>
class NodeBase {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  // Constructor with id and coordinates
  //! \param[in] id assign as the id_ of the node
  //! \param[in] coords coordinates of the node
  NodeBase(mpm::Index id, const VectorDim& coords){};

  //! Destructor
  virtual ~NodeBase(){};

  //! Delete copy constructor
  NodeBase(const NodeBase<Tdim>&) = delete;

  //! Delete assignement operator
  NodeBase& operator=(const NodeBase<Tdim>&) = delete;

  //! Return id of the nodebase
  virtual Index id() const = 0;

  //! Initialise shared pointer to nodal properties pool
  //! \param[in] prop_id Property id in the nodal property pool
  //! \param[in] nodal_properties Shared pointer to nodal properties pool
  virtual void initialise_property_handle(
      unsigned prop_id,
      std::shared_ptr<mpm::NodalProperties> property_handle) noexcept = 0;

  //! Assign coordinates
  virtual void assign_coordinates(const VectorDim& coord) = 0;

  //! Return coordinates
  //! \retval coordinates_ return coordinates of the nodebase
  virtual VectorDim coordinates() const = 0;

  //! Initialise properties
  virtual void initialise() noexcept = 0;

  //! Return degrees of freedom
  virtual unsigned dof() const = 0;

  //! Assign status
  virtual void assign_status(bool status) = 0;

  //! Return status
  virtual bool status() const = 0;

  //! Assign boolean property at the nodes
  //! \param[in] property Name of the property to assign
  //! \param[in] boolean Property boolean (true/false) of the node
  virtual void assign_boolean_property(mpm::properties::Boolean property,
                                       bool boolean) noexcept = 0;

  //! Return boolean property
  //! \param[in] property Name of the property to update
  //! \retval boolean property at node
  virtual bool boolean_property(mpm::properties::Boolean property) const = 0;

  //! Update scalar property at the nodes
  //! \param[in] property Name of the property to update
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] value Property value from the particles in a cell
  virtual void update_scalar_property(mpm::properties::Scalar property,
                                      bool update, unsigned phase,
                                      double value) noexcept = 0;

  //! Return property at a given node for a given phase
  //! \param[in] property Name of the property to return
  //! \param[in] phase Index corresponding to the phase
  //! \retval scalar property at the designated phase
  virtual double scalar_property(mpm::properties::Scalar property,
                                 unsigned phase) const = 0;

  //! Update vector property at the nodes
  //! \param[in] property Name of the property to update
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] value Property value from the particles in a cell
  virtual void update_vector_property(
      mpm::properties::Vector property, bool update, unsigned phase,
      const Eigen::Matrix<double, Tdim, 1>& value) noexcept = 0;

  //! Return property at a given node for a given phase
  //! \param[in] property Name of the property to return
  //! \param[in] phase Index corresponding to the phase
  //! \retval vector property at the designated phase
  virtual Eigen::Matrix<double, Tdim, 1> vector_property(
      mpm::properties::Vector property, unsigned phase) const = 0;

  //! Update mass at the nodes from particle
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] mass Mass from the particles in a cell
  virtual void update_mass(bool update, unsigned phase,
                           double mass) noexcept = 0;

  //! Return mass at a given node for a given phase
  virtual double mass(unsigned phase) const = 0;

  //! Update volume at the nodes from particle
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] volume Volume from the particles in a cell
  virtual void update_volume(bool update, unsigned phase,
                             double volume) noexcept = 0;

  //! Return volume at a given node for a given phase
  virtual double volume(unsigned phase) const = 0;

  //! Assign concentrated force to the node
  //! \param[in] direction Index corresponding to the direction of traction
  //! \param[in] traction Nodal concentrated force in specified direction
  //! \param[in] function math function
  //! \retval status Assignment status
  virtual bool assign_concentrated_force(
      unsigned phase, unsigned direction, double traction,
      const std::shared_ptr<FunctionBase>& function) = 0;

  //! Apply concentrated force to external force
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] current time
  virtual void apply_concentrated_force(unsigned phase,
                                        double current_time) = 0;

  //! Update external force (body force / traction force)
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] force External force from the particles in a cell
  virtual void update_external_force(bool update, unsigned phase,
                                     const VectorDim& force) noexcept = 0;

  //! Return external force
  //! \param[in] phase Index corresponding to the phase
  virtual VectorDim external_force(unsigned phase) const = 0;

  //! Update internal force (body force / traction force)
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] force Internal force from the particles in a cell
  virtual void update_internal_force(bool update, unsigned phase,
                                     const VectorDim& force) noexcept = 0;

  //! Return internal force
  //! \param[in] phase Index corresponding to the phase
  virtual VectorDim internal_force(unsigned phase) const = 0;

  //! Update pressure at the nodes from particle
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] mass_pressure Product of mass x pressure of a particle
  virtual void update_mass_pressure(bool update, unsigned phase,
                                    double mass_pressure) noexcept = 0;

  //! Compute pressure from the mass pressure
  virtual void compute_pressure() = 0;

  //! Assign pressure at the nodes from particle
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] pressure Pressure of a particle
  virtual void update_pressure(bool update, unsigned phase,
                               double pressure) = 0;

  //! Return pressure at a given node for a given phase
  //! \param[in] phase Index corresponding to the phase
  virtual double pressure(unsigned phase) const = 0;

  //! Update nodal momentum
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] momentum Momentum from the particles in a cell
  virtual void update_momentum(bool update, unsigned phase,
                               const VectorDim& momentum) noexcept = 0;

  //! Return momentum
  //! \param[in] phase Index corresponding to the phase
  virtual VectorDim momentum(unsigned phase) const = 0;

  //! Compute velocity from the momentum
  virtual void compute_velocity() = 0;

  //! Return velocity
  //! \param[in] phase Index corresponding to the phase
  virtual VectorDim velocity(unsigned phase) const = 0;

  //! Update nodal acceleration
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] acceleration Acceleration from the particles in a cell
  virtual void update_acceleration(bool update, unsigned phase,
                                   const VectorDim& acceleration) = 0;

  //! Return acceleration
  //! \param[in] phase Index corresponding to the phase
  virtual VectorDim acceleration(unsigned phase) const = 0;

  //! Compute acceleration
  //! \param[in] dt Time-step
  virtual bool compute_acceleration_velocity(unsigned phase,
                                             double dt) noexcept = 0;

  //! Compute acceleration and velocity with cundall damping factor
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] dt Timestep in analysis
  //! \param[in] damping_factor Damping factor
  virtual bool compute_acceleration_velocity_cundall(
      unsigned phase, double dt, double damping_factor) noexcept = 0;

  //! Assign velocity constraint
  //! Directions can take values between 0 and Dim * Nphases
  //! \param[in] dir Direction of velocity constraint
  //! \param[in] velocity Applied velocity constraint
  virtual bool assign_velocity_constraint(unsigned dir, double velocity) = 0;

  //! Apply velocity constraints
  virtual void apply_velocity_constraints() = 0;

  //! Assign friction constraint
  //! Directions can take values between 0 and Dim * Nphases
  //! \param[in] dir Direction of friction constraint
  //! \param[in] sign Sign of normal wrt coordinate system for friction
  //! \param[in] friction Applied friction constraint
  virtual bool assign_friction_constraint(unsigned dir, int sign,
                                          double friction) = 0;

  //! Apply friction constraints
  //! \param[in] dt Time-step
  virtual void apply_friction_constraints(double dt) = 0;

  //! Assign rotation matrix
  //! \param[in] rotation_matrix Rotation matrix of the node
  virtual void assign_rotation_matrix(
      const Eigen::Matrix<double, Tdim, Tdim>& rotation_matrix) = 0;

  //! Add material id from material points to list of materials in materials_
  //! \param[in] id Material id to be stored at the node
  virtual void append_material_id(unsigned id) = 0;

  //! Return material ids in node
  virtual std::set<unsigned> material_ids() const = 0;

  //! Assign MPI rank to node
  //! \param[in] rank MPI Rank of the node
  virtual bool mpi_rank(unsigned rank) = 0;

  //! Assign MPI rank to node
  //! \param[in] rank MPI Rank of the node
  virtual std::set<unsigned> mpi_ranks() const = 0;

  //! Clear MPI ranks on node
  virtual void clear_mpi_ranks() = 0;

  //! Return ghost id
  virtual Index ghost_id() const = 0;

  //! Set ghost id
  virtual void ghost_id(Index gid) = 0;

  //! Update nodal property at the nodes from particle
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] property Property name
  //! \param[in] property_value Property quantity from the particles in the cell
  //! \param[in] mat_id Id of the material within the property data
  //! \param[in] nprops Dimension of property (1 if scalar, Tdim if vector)
  virtual void update_property(bool update, const std::string& property,
                               const Eigen::MatrixXd& property_value,
                               unsigned mat_id, unsigned nprops) noexcept = 0;

  //! Compute multimaterial change in momentum
  virtual void compute_multimaterial_change_in_momentum() = 0;

  //! Compute multimaterial separation vector
  virtual void compute_multimaterial_separation_vector() = 0;

  //! Compute multimaterial normal unit vector
  virtual void compute_multimaterial_normal_unit_vector() = 0;

};  // NodeBase class
}  // namespace mpm

#endif  // MPM_NODE_BASE_H_
