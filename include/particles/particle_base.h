#ifndef MPM_PARTICLEBASE_H_
#define MPM_PARTICLEBASE_H_

// MPI
#ifdef USE_MPI
#include "mpi.h"
#endif

#include <array>
#include <limits>
#include <memory>
#include <vector>

#include "cell.h"
#include "data_types.h"
#include "function_base.h"
#include "material.h"
#include "pod_particle.h"
#include "pod_particle_twophase.h"

namespace mpm {

// Forward declaration of Material
template <unsigned Tdim>
class Material;

//! Particle phases
enum ParticlePhase : unsigned int {
  SinglePhase = 0,
  Mixture = 0,
  Solid = 0,
  Liquid = 1,
  Gas = 2
};

//! Particle type
extern std::map<std::string, int> ParticleType;
extern std::map<int, std::string> ParticleTypeName;
extern std::map<std::string, std::string> ParticleHDF5TypeName;

//! ParticleBase class
//! \brief Base class that stores the information about particleBases
//! \details ParticleBase class: id_ and coordinates.
//! \tparam Tdim Dimension
template <unsigned Tdim>
class ParticleBase {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Constructor with id and coordinates
  //! \param[in] id Particle id
  //! \param[in] coord coordinates of the particle
  ParticleBase(Index id, const VectorDim& coord);

  //! Constructor with id, coordinates and status
  //! \param[in] id Particle id
  //! \param[in] coord coordinates of the particle
  //! \param[in] status Particle status (active / inactive)
  ParticleBase(Index id, const VectorDim& coord, bool status);

  //! Destructor
  virtual ~ParticleBase(){};

  //! Delete copy constructor
  ParticleBase(const ParticleBase<Tdim>&) = delete;

  //! Delete assignement operator
  ParticleBase& operator=(const ParticleBase<Tdim>&) = delete;

  //! Initialise particle POD data
  //! \param[in] particle POD data of particle
  //! \retval status Status of reading POD particle
  virtual bool initialise_particle(PODParticle& particle) = 0;

  //! Initialise particle POD data and material
  //! \param[in] particle POD data of particle
  //! \param[in] materials Material associated with the particle arranged in a
  //! vector
  //! \retval status Status of reading POD particle
  virtual bool initialise_particle(
      PODParticle& particle,
      const std::vector<std::shared_ptr<Material<Tdim>>>& materials) = 0;

  //! Return particle data as POD
  //! \retval particle POD of the particle
  virtual std::shared_ptr<void> pod() = 0;

  //! Return id of the particleBase
  Index id() const { return id_; }

  //! Assign coordinates
  //! \param[in] coord Assign coord as coordinates of the particleBase
  void assign_coordinates(const VectorDim& coord) { coordinates_ = coord; }

  //! Return coordinates
  //! \retval coordinates_ return coordinates of the particleBase
  VectorDim coordinates() const { return coordinates_; }

  //! Compute reference coordinates in a cell
  virtual bool compute_reference_location() = 0;

  //! Return reference location
  virtual VectorDim reference_location() const = 0;

  //! Assign cell
  virtual bool assign_cell(const std::shared_ptr<Cell<Tdim>>& cellptr) = 0;

  //! Assign cell and xi
  virtual bool assign_cell_xi(const std::shared_ptr<Cell<Tdim>>& cellptr,
                              const Eigen::Matrix<double, Tdim, 1>& xi) = 0;

  //! Assign cell id
  virtual bool assign_cell_id(Index id) = 0;

  //! Return cell id
  virtual Index cell_id() const = 0;

  //! Return cell ptr status
  virtual bool cell_ptr() const = 0;

  //! Remove cell
  virtual void remove_cell() = 0;

  //! Compute shape functions
  virtual void compute_shapefn() noexcept = 0;

  //! Assign volume
  virtual bool assign_volume(double volume) = 0;

  //! Return volume
  virtual double volume() const = 0;

  //! Return the approximate particle diameter
  virtual double diameter() const = 0;

  //! Return size of particle in natural coordinates
  virtual VectorDim natural_size() const = 0;

  //! Compute volume of particle
  virtual void compute_volume() noexcept = 0;

  //! Update volume based on centre volumetric strain rate
  virtual void update_volume() noexcept = 0;

  //! Return mass density
  virtual double mass_density() const = 0;

  //! Compute mass of particle
  virtual void compute_mass() noexcept = 0;

  //! Map particle mass and momentum to nodes
  virtual void map_mass_momentum_to_nodes() noexcept = 0;

  //! Map multimaterial properties to nodes
  virtual void map_multimaterial_mass_momentum_to_nodes() noexcept = 0;

  //! Map multimaterial displacements to nodes
  virtual void map_multimaterial_displacements_to_nodes() noexcept = 0;

  //! Map multimaterial domain gradients to nodes
  virtual void map_multimaterial_domain_gradients_to_nodes() noexcept = 0;

  //! Assign material
  virtual bool assign_material(const std::shared_ptr<Material<Tdim>>& material,
                               unsigned phase = mpm::ParticlePhase::Solid) = 0;

  //! Return material of particle
  //! \param[in] phase Index to indicate material phase
  std::shared_ptr<Material<Tdim>> material(
      unsigned phase = mpm::ParticlePhase::Solid) const {
    return material_[phase];
  }

  //! Return material id
  //! \param[in] phase Index to indicate material phase
  unsigned material_id(unsigned phase = mpm::ParticlePhase::Solid) const {
    return material_id_[phase];
  }

  //! Assign material state variables
  virtual bool assign_material_state_vars(
      const mpm::dense_map& state_vars,
      const std::shared_ptr<mpm::Material<Tdim>>& material,
      unsigned phase = mpm::ParticlePhase::Solid) = 0;

  //! Return state variables
  //! \param[in] phase Index to indicate material phase
  mpm::dense_map state_variables(
      unsigned phase = mpm::ParticlePhase::Solid) const {
    return state_variables_[phase];
  }

  //! Assign a state variable
  virtual void assign_state_variable(
      const std::string& var, double value,
      unsigned phase = mpm::ParticlePhase::Solid) = 0;

  //! Return a state variable
  virtual double state_variable(
      const std::string& var,
      unsigned phase = mpm::ParticlePhase::Solid) const = 0;

  //! Assign status
  void assign_status(bool status) { status_ = status; }

  //! Status
  bool status() const { return status_; }

  //! Initialise properties
  virtual void initialise() = 0;

  //! Assign mass
  virtual void assign_mass(double mass) = 0;

  //! Return mass
  virtual double mass() const = 0;

  //! Assign pressure
  virtual void assign_pressure(double pressure,
                               unsigned phase = mpm::ParticlePhase::Solid) = 0;

  //! Return pressure
  virtual double pressure(unsigned phase = mpm::ParticlePhase::Solid) const = 0;

  //! Compute strain
  virtual void compute_strain(double dt) noexcept = 0;

  //! Strain
  virtual Eigen::Matrix<double, 6, 1> strain() const = 0;

  //! Strain rate
  virtual Eigen::Matrix<double, 6, 1> strain_rate() const = 0;

  //! Volumetric strain of centroid
  virtual double volumetric_strain_centroid() const = 0;

  //! dvolumetric strain
  virtual double dvolumetric_strain() const = 0;

  //! Initial stress
  virtual void initial_stress(const Eigen::Matrix<double, 6, 1>& stress) = 0;

  //! Compute stress
  virtual void compute_stress() noexcept = 0;

  //! Return stress
  virtual Eigen::Matrix<double, 6, 1> stress() const = 0;

  //! Map body force
  virtual void map_body_force(const VectorDim& pgravity) noexcept = 0;

  //! Map internal force
  virtual void map_internal_force() noexcept = 0;

  //! Map particle pressure to nodes
  virtual bool map_pressure_to_nodes(
      unsigned phase = mpm::ParticlePhase::Solid) noexcept = 0;

  //! Compute pressure smoothing of the particle based on nodal pressure
  virtual bool compute_pressure_smoothing(
      unsigned phase = mpm::ParticlePhase::Solid) noexcept = 0;

  //! Assign velocity
  virtual bool assign_velocity(const VectorDim& velocity) = 0;

  //! Return velocity
  virtual VectorDim velocity() const = 0;

  //! Return displacement of the particle
  virtual VectorDim displacement() const = 0;

  //! Assign traction
  virtual bool assign_traction(unsigned direction, double traction) = 0;

  //! Return traction
  virtual VectorDim traction() const = 0;

  //! Map traction force
  virtual void map_traction_force() noexcept = 0;

  //! Compute updated position
  virtual void compute_updated_position(
      double dt, bool velocity_update = false) noexcept = 0;

  //! Return scalar data of particles
  //! \param[in] property Property string
  //! \retval data Scalar data of particle property
  virtual double scalar_data(const std::string& property) const = 0;

  //! Return vector data of particles
  //! \param[in] property Property string
  //! \retval data Vector data of particle property
  virtual VectorDim vector_data(const std::string& property) const = 0;

  //! Return tensor data of particles
  //! \param[in] property Property string
  //! \retval data Tensor data of particle property
  virtual Eigen::VectorXd tensor_data(const std::string& property) const = 0;

  //! Apply particle velocity constraints
  //! \param[in] dir Direction of particle velocity constraint
  //! \param[in] velocity Applied particle velocity constraint
  virtual void apply_particle_velocity_constraints(unsigned dir,
                                                   double velocity) = 0;

  //! Assign material id of this particle to nodes
  virtual void append_material_id_to_nodes() const = 0;

  //! Assign particle free surface
  virtual void assign_free_surface(bool free_surface) = 0;

  //! Assign particle free surface
  virtual bool free_surface() = 0;

  //! Compute free surface in particle level by density ratio comparison
  virtual bool compute_free_surface_by_density(
      double density_ratio_tolerance = 0.65) = 0;

  //! Assign normal vector
  virtual void assign_normal(const VectorDim& normal) = 0;

  //! Return normal vector
  virtual VectorDim normal() = 0;

  //! Return the number of neighbour particles
  virtual unsigned nneighbours() const = 0;

  //! Assign neighbour particles
  //! \param[in] neighbours set of id of the neighbouring particles
  //! \retval insertion_status Return the successful addition of a node
  virtual void assign_neighbours(const std::vector<mpm::Index>& neighbours) = 0;

  //! Return neighbour ids
  virtual std::vector<mpm::Index> neighbours() const = 0;

  //! Type of particle
  virtual std::string type() const = 0;

  //! Serialize
  //! \retval buffer Serialized buffer data
  virtual std::vector<uint8_t> serialize() = 0;

  //! Deserialize
  //! \param[in] buffer Serialized buffer data
  //! \param[in] material Particle material pointers
  virtual void deserialize(
      const std::vector<uint8_t>& buffer,
      std::vector<std::shared_ptr<mpm::Material<Tdim>>>& materials) = 0;

  //! Navier-Stokes functions----------------------------------
  //! Assigning beta parameter to particle
  //! \param[in] pressure parameter determining type of projection
  virtual void assign_projection_parameter(double parameter) {
    throw std::runtime_error(
        "Calling the base class function (assign_projection_parameter) in "
        "ParticleBase:: illegal operation!");
  };

  //! Map laplacian element matrix to cell (used in poisson equation LHS)
  virtual bool map_laplacian_to_cell() {
    throw std::runtime_error(
        "Calling the base class function (map_laplacian_to_cell) in "
        "ParticleBase:: "
        "illegal operation!");
    return 0;
  };

  //! Map poisson rhs element matrix to cell (used in poisson equation RHS)
  virtual bool map_poisson_right_to_cell() {
    throw std::runtime_error(
        "Calling the base class function (map_poisson_right_to_cell) in "
        "ParticleBase:: "
        "illegal operation!");
    return 0;
  };

  //! Map correction matrix element matrix to cell (used to correct velocity)
  virtual bool map_correction_matrix_to_cell() {
    throw std::runtime_error(
        "Calling the base class function (map_correction_matrix_to_cell) in "
        "ParticleBase:: "
        "illegal operation!");
    return 0;
  };

  //! Update pressure after solving poisson equation
  virtual bool compute_updated_pressure() {
    throw std::runtime_error(
        "Calling the base class function (compute_updated_pressure) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  //! TwoPhase functions--------------------------------------------------------
  //! Update porosity
  //! \param[in] dt Analysis time step
  virtual void update_porosity(double dt) {
    throw std::runtime_error(
        "Calling the base class function (update_porosity) in "
        "ParticleBase:: illegal operation!");
  };

  //! Assign saturation degree
  virtual bool assign_saturation_degree() {
    throw std::runtime_error(
        "Calling the base class function (assign_saturation_degree) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  //! Assign velocity to the particle liquid phase
  //! \param[in] velocity A vector of particle liquid phase velocity
  //! \retval status Assignment status
  virtual bool assign_liquid_velocity(const VectorDim& velocity) {
    throw std::runtime_error(
        "Calling the base class function (assign_liquid_velocity) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  //! Compute pore pressure
  //! \param[in] dt Time step size
  virtual void compute_pore_pressure(double dt) {
    throw std::runtime_error(
        "Calling the base class function (compute_pore_pressure) in "
        "ParticleBase:: illegal operation!");
  };

  //! Map drag force coefficient
  virtual bool map_drag_force_coefficient() {
    throw std::runtime_error(
        "Calling the base class function (map_drag_force_coefficient) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  //! Initialise particle pore pressure by watertable
  virtual bool initialise_pore_pressure_watertable(
      const unsigned dir_v, const unsigned dir_h, const VectorDim& gravity,
      std::map<double, double>& reference_points) {
    throw std::runtime_error(
        "Calling the base class function "
        "(initial_pore_pressure_watertable) in "
        "ParticleBase:: illegal operation!");
    return false;
  };

  //! Initialise particle pore pressure by watertable
  virtual bool assign_porosity() {
    throw std::runtime_error(
        "Calling the base class function "
        "(assign_porosity) in "
        "ParticleBase:: illegal operation!");
    return false;
  };

  //! Initialise particle pore pressure by watertable
  virtual bool assign_permeability() {
    throw std::runtime_error(
        "Calling the base class function "
        "(assign_permeability) in "
        "ParticleBase:: illegal operation!");
    return false;
  };

  //! Return liquid mass
  //! \retval liquid mass Liquid phase mass
  virtual double liquid_mass() const {
    throw std::runtime_error(
        "Calling the base class function (liquid_mass) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  //! Return velocity of the particle liquid phase
  //! \retval liquid velocity Liquid phase velocity
  virtual VectorDim liquid_velocity() const {
    auto error = VectorDim::Zero();
    throw std::runtime_error(
        "Calling the base class function (liquid_velocity) in "
        "ParticleBase:: illegal operation!");
    return error;
  };

  //! Return porosity
  //! \retval porosity Porosity
  virtual double porosity() const {
    throw std::runtime_error(
        "Calling the base class function (porosity) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };
  //----------------------------------------------------------------------------

 protected:
  //! particleBase id
  Index id_{std::numeric_limits<Index>::max()};
  //! coordinates
  VectorDim coordinates_;
  //! Cell id
  Index cell_id_{std::numeric_limits<Index>::max()};
  //! Status
  bool status_{true};
  //! Reference coordinates (in a cell)
  Eigen::Matrix<double, Tdim, 1> xi_;
  //! Cell
  std::shared_ptr<Cell<Tdim>> cell_;
  //! Vector of nodal pointers
  std::vector<std::shared_ptr<NodeBase<Tdim>>> nodes_;
  //! Material
  std::vector<std::shared_ptr<Material<Tdim>>> material_;
  //! Unsigned material id
  std::vector<unsigned> material_id_;
  //! Material state history variables
  std::vector<mpm::dense_map> state_variables_;
  //! Vector of particle neighbour ids
  std::vector<mpm::Index> neighbours_;
};  // ParticleBase class
}  // namespace mpm

#include "particle_base.tcc"

#endif  // MPM_PARTICLEBASE_H__
