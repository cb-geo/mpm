#ifndef MPM_NODE_BASE_H_
#define MPM_NODE_BASE_H_

#include <array>
#include <limits>
#include <map>
#include <vector>

#include <Eigen/Dense>

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

  //! Assign coordinates
  virtual void assign_coordinates(const VectorDim& coord) = 0;

  //! Return coordinates
  //! \retval coordinates_ return coordinates of the nodebase
  virtual VectorDim coordinates() const = 0;

  //! Initialise properties
  virtual void initialise() = 0;

  //! Initialise nodal properties of a subdomain
  virtual void initialise_subdomain(unsigned mid) = 0;

  //! Return degrees of freedom
  virtual unsigned dof() const = 0;

  //! Assign status
  virtual void assign_status(bool status) = 0;

  //! Return status
  virtual bool status() const = 0;

  //! Update mass at the nodes from particle
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] mass Mass from the particles in a cell
  virtual void update_mass(bool update, unsigned phase, double mass) = 0;

  //! Update mass of a subdomain at the nodes from particle
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] mass Mass from the particles in a cell
  //! \param[in] mid Subdomain id (material id)
  virtual void update_mass_subdomain(bool update, unsigned phase, double mass,
                                     unsigned mid) = 0;

  //! Return mass at a given node for a given phase
  virtual double mass(unsigned phase) const = 0;

  //! Update volume at the nodes from particle
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] volume Volume from the particles in a cell
  virtual void update_volume(bool update, unsigned phase, double volume) = 0;

  //! Return volume at a given node for a given phase
  virtual double volume(unsigned phase) const = 0;

  //! Assign traction force to the node
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] direction Index corresponding to the direction of traction
  //! \param[in] traction Nodal traction in specified direction
  //! \retval status Assignment status
  virtual bool assign_traction_force(unsigned phase, unsigned direction,
                                     double traction) = 0;

  //! Assign traction force of a subdomain to the node
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] direction Index corresponding to the direction of traction
  //! \param[in] traction Nodal traction in specified direction
  //! \param[in] mid Subdomain id (material id)
  //! \retval status Assignment status
  virtual bool assign_traction_force_subdomain(unsigned phase,
                                               unsigned direction,
                                               double traction,
                                               unsigned mid) = 0;

  //! Update external force (body force / traction force)
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] force External force from the particles in a cell
  //! \retval status Update status
  virtual bool update_external_force(bool update, unsigned phase,
                                     const VectorDim& force) = 0;

  //! Update external force of a subdomain (body force / traction force)
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] force External force from the particles in a cell
  //! \param[in] mid Subdomain id (material id)
  //! \retval status Update status
  virtual bool update_external_force_subdomain(bool update, unsigned phase,
                                               const VectorDim& force,
                                               unsigned mid) = 0;

  //! Return external force
  //! \param[in] phase Index corresponding to the phase
  virtual VectorDim external_force(unsigned phase) const = 0;

  //! Update internal force (body force / traction force)
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] force Internal force from the particles in a cell
  //! \retval status Update status
  virtual bool update_internal_force(bool update, unsigned phase,
                                     const VectorDim& force) = 0;

  //! Update internal force of a subdomain (body force / traction force)
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] force Internal force from the particles in a cell
  //! \param[in] mid Subdomain id (material id)
  //! \retval status Update status
  virtual bool update_internal_force_subdomain(bool update, unsigned phase,
                                               const VectorDim& force,
                                               unsigned mid) = 0;

  //! Return internal force
  //! \param[in] phase Index corresponding to the phase
  virtual VectorDim internal_force(unsigned phase) const = 0;

  //! Update pressure at the nodes from particle
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] mass_pressure Product of mass x pressure of a particle
  virtual void update_pressure(bool update, unsigned phase,
                               double mass_pressure) = 0;

  //! Return pressure at a given node for a given phase
  //! \param[in] phase Index corresponding to the phase
  virtual double pressure(unsigned phase) const = 0;

  //! Update nodal momentum
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] momentum Momentum from the particles in a cell
  //! \retval status Update status
  virtual bool update_momentum(bool update, unsigned phase,
                               const VectorDim& momentum) = 0;

  //! Update momentum of a subdomain at the nodes
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] momentum Momentum from the particles in a cell
  //! \param[in] mid subdomain id (material id)
  //! \retval status Update status
  virtual bool update_momentum_subdomain(bool update, unsigned phase,
                                         const VectorDim& momentum,
                                         unsigned mid) = 0;

  //! Update coordinates at the nodes
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] shapefn_coordinates Shapefn * coordinates from the particles in
  //! a cell \retval status Update status
  virtual bool update_coordinates(
      bool update, double pmass,
      const Eigen::Matrix<double, Tdim, 1>& shapefn_coordinates) = 0;

  //! Update coordinates of a subdomain at the nodes
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] shapefn_coordinates Shapefn * coordinates from the particles in
  //! a cell \param[in] mid subdomain id (material id) \retval status Update
  //! status
  virtual bool update_coordinates_subdomain(
      bool update, double pmass,
      const Eigen::Matrix<double, Tdim, 1>& shapefn_coordinates,
      unsigned mid) = 0;

  //! Return momentum
  //! \param[in] phase Index corresponding to the phase
  virtual VectorDim momentum(unsigned phase) const = 0;

  //! Compute velocity from the momentum
  virtual void compute_velocity() = 0;

  //! Compute velocity of a subdomain from the momentum
  virtual void compute_velocity_subdomain(unsigned mid) = 0;

  //! Return velocity
  //! \param[in] phase Index corresponding to the phase
  virtual VectorDim velocity(unsigned phase) const = 0;

  //! Return velocity of a subdomain at a given node for a given phase
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] mid Subdomain id (material id)
  virtual VectorDim velocity_subdomain(unsigned phase, unsigned mid) const = 0;

  //! Update nodal acceleration
  //! \param[in] update A boolean to update (true) or assign (false)
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] acceleration Acceleration from the particles in a cell
  //! \retval status Update status
  virtual bool update_acceleration(bool update, unsigned phase,
                                   const VectorDim& acceleration) = 0;

  //! Return acceleration
  //! \param[in] phase Index corresponding to the phase
  virtual VectorDim acceleration(unsigned phase) const = 0;

  //! Return acceleration of a subdomain at a given node for a given phase
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] mid Subdomain id (material id)
  virtual VectorDim acceleration_subdomain(unsigned phase,
                                           unsigned mid) const = 0;

  //! Compute acceleration
  //! \param[in] dt Time-step
  virtual bool compute_acceleration_velocity(unsigned phase, double dt) = 0;

  //! Compute acceleration and velocity of a subdomain
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] dt Timestep in analysis
  //! \param[in] mid Subdomain id (material id)
  virtual bool compute_acceleration_velocity_subdomain(unsigned phase,
                                                       double dt,
                                                       unsigned mid) = 0;

  //! Assign velocity constraint
  //! Directions can take values between 0 and Dim * Nphases
  //! \param[in] dir Direction of velocity constraint
  //! \param[in] velocity Applied velocity constraint
  virtual bool assign_velocity_constraint(unsigned dir, double velocity) = 0;

  //! Apply velocity constraints
  virtual void apply_velocity_constraints() = 0;

  //! Apply velocity constraints of a subdomain
  virtual void apply_velocity_constraints_subdomain(unsigned mid) = 0;

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

  //! Apply friction constraints of a subdomain
  //! \param[in] dt Time-step
  //! \param[in] mid Subdomain id (material id)
  virtual void apply_friction_constraints_subdomain(double dt,
                                                    unsigned mid) = 0;

  //! Assign rotation matrix
  //! \param[in] rotation_matrix Rotation matrix of the node
  virtual void assign_rotation_matrix(
      const Eigen::Matrix<double, Tdim, Tdim>& rotation_matrix) = 0;

  //! Compute normal vector of contact interface of a subdomain
  //! \param[in] normal_vector_type Computation method
  //! \param[in] normal_vector Default normal vector
  //! \param[in] mid Material id of subdomain
  virtual void compute_normal_vector(
      const std::string normal_vector_type,
      const Eigen::Matrix<double, Tdim, 1> normal_vector, unsigned mid) = 0;

  //! Compute contact components
  //! \param[in] separation_vector Separation vector
  //! \param[in] separation_normal_length Length of separation vector in normal
  //! direction \param[in] separation_tangent_length Length of separation vector
  //! in tangent direction
  //! \param[in] mid Material id of subdomain
  virtual void compute_contact_components(
      const Eigen::Matrix<double, Tdim, 1> separation_vector,
      double* separation_normal_length, double* separation_tangent_length,
      unsigned mid) = 0;

  //! Compute correction momentum
  //! \param[in] friction_type Iterface friction type (rough or frictional)
  //! \param[in] friction_coefficient Friction coefficient
  //! \param[in] momentum_difference Difference between total momentum and
  //! subdomain momentum \param[in] separation_normal_length Length of
  //! separation vector in normal direction
  //! \param[in] mid Material id of subdomain
  virtual void compute_correction_momentum(
      const std::string friction_type, const double friction_coefficient,
      const Eigen::Matrix<double, Tdim, 1> momentum_difference,
      const double separation_normal_length, unsigned mid) = 0;

  //! Compute contact force
  //! \param[in] separation_normal_length Length of separation vector in normal
  //! direction \param[in] separation_normal_length Length of separation vector
  //! \param[in] separation_cut_off Cut off value of separation check
  //! \param[in] dc_n Normal contact stiffness
  //! \param[in] dc_t Tangent contact stiffness
  //! in tangent direction \param[in] mid Subdomain id (material id)
  virtual void compute_contact_force(const double separation_normal_length,
                                     const double separation_tangent_length,
                                     const double separation_cut_off,
                                     const double dc_n, const double dc_t,
                                     unsigned mid) = 0;

  //! Compute contact interface on node
  //! \param[in] friction_type Iterface friction type (rough or frictional)
  //! \param[in] friction_coefficient Friction coefficient
  //! \param[in] separation_cut_off Cut off value of separation check
  //! \param[in] dc_n Normal contact stiffness
  //! \param[in] dc_t Tangent contact stiffness
  //! \param[in] mid Material id of subdomain
  virtual void compute_contact_interface(bool contact_force,
                                         const std::string friction_type,
                                         const double friction_coefficient,
                                         const double separation_cut_off,
                                         const double dc_n, const double dc_t,
                                         unsigned mid) = 0;
};  // NodeBase class
}  // namespace mpm

#endif  // MPM_NODE_BASE_H_
