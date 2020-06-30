#ifndef MPM_PARTICLEBASE_H_
#define MPM_PARTICLEBASE_H_

#include <array>
#include <limits>
#include <memory>
#include <vector>

#include "cell.h"
#include "data_types.h"
#include "function_base.h"
#include "hdf5_particle.h"
#include "material.h"

namespace mpm {

// Forward declaration of Material
template <unsigned Tdim>
class Material;

//! Particle phases
enum ParticlePhase : unsigned int {
  Mixture = 0,
  Solid = 0,
  Liquid = 1,
  Gas = 2
};

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

  //! Initialise particle HDF5 data
  //! \param[in] particle HDF5 data of particle
  //! \retval status Status of reading HDF5 particle
  virtual bool initialise_particle(const HDF5Particle& particle) = 0;

  //! Initialise particle HDF5 data and material
  //! \param[in] particle HDF5 data of particle
  //! \param[in] material Material associated with the particle
  //! \retval status Status of reading HDF5 particle
  virtual bool initialise_particle(
      const HDF5Particle& particle,
      const std::shared_ptr<Material<Tdim>>& material) = 0;

  //! Assign material history variables
  //! \param[in] state_vars State variables
  //! \param[in] material Material associated with the particle
  //! \retval status Status of cloning HDF5 particle
  virtual bool assign_material_state_vars(
      const mpm::dense_map& state_vars,
      const std::shared_ptr<mpm::Material<Tdim>>& material) = 0;

  //! Retrun particle data as HDF5
  //! \retval particle HDF5 data of the particle
  virtual HDF5Particle hdf5() const = 0;

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

  //! Assign material
  virtual bool assign_material(
      const std::shared_ptr<Material<Tdim>>& material) = 0;

  //! Return material id
  unsigned material_id() const { return material_id_; }

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

  //! Return pressure
  virtual double pressure() const = 0;

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
  virtual void initial_stress(const Eigen::Matrix<double, 6, 1>&) = 0;

  //! Compute stress
  virtual void compute_stress() noexcept = 0;

  //! Return stress
  virtual Eigen::Matrix<double, 6, 1> stress() const = 0;

  //! Map body force
  virtual void map_body_force(const VectorDim& pgravity) noexcept = 0;

  //! Map internal force
  virtual void map_internal_force() noexcept = 0;

  //! Map particle pressure to nodes
  virtual bool map_pressure_to_nodes() noexcept = 0;

  //! Compute pressure smoothing of the particle based on nodal pressure
  virtual bool compute_pressure_smoothing() noexcept = 0;

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

  //! Return a state variable
  virtual double state_variable(const std::string& var) const = 0;

  //! Return tensor data of particles
  //! \param[in] property Property string
  //! \retval vecdata Tensor data of particle property
  virtual Eigen::VectorXd tensor_data(const std::string& property) = 0;

  //! Apply particle velocity constraints
  //! \param[in] dir Direction of particle velocity constraint
  //! \param[in] velocity Applied particle velocity constraint
  virtual void apply_particle_velocity_constraints(unsigned dir,
                                                   double velocity) = 0;

  //! Assign material id of this particle to nodes
  virtual void append_material_id_to_nodes() const = 0;

  //! Return the number of neighbour particles
  virtual unsigned nneighbours() const = 0;

  //! Assign neighbour particles
  //! \param[in] neighbours set of id of the neighbouring particles
  //! \retval insertion_status Return the successful addition of a node
  virtual void assign_neighbours(const std::vector<mpm::Index>& neighbours) = 0;

  //! Return neighbour ids
  virtual std::vector<mpm::Index> neighbours() const = 0;

  // Twophase particle functions------------------------------------------------

  //! Assign porosity
  virtual bool assign_porosity() {
    throw std::runtime_error(
        "Calling the base class function "
        "(assign_porosity) in "
        "ParticleBase:: illegal operation!");
    return false;
  };

  //! Assign particle free surface
  virtual bool free_surface() {
    throw std::runtime_error(
        "Calling the base class function "
        "(free_surface) in "
        "ParticleBase:: illegal operation!");
    return false;
  };

  //! Assign particle pressure constraints
  virtual bool assign_particle_pore_pressure_constraint(double pressure) {
    throw std::runtime_error(
        "Calling the base class function "
        "(assign_particle_pore_pressure_constraint) in "
        "ParticleBase:: illegal operation!");
    return false;
  };

  //! Initialise liquid phase
  virtual void initialise_liquid_phase() {
    throw std::runtime_error(
        "Calling the base class function (initialise_liquid_phase) in "
        "ParticleBase:: illegal operation!");
  };

  //! Assign material
  //! \param[in] material Pointer to a material
  virtual bool assign_liquid_material(
      const std::shared_ptr<Material<Tdim>>& material) {
    throw std::runtime_error(
        "Calling the base class function (assign_liquid_material) in "
        "ParticleBase:: illegal operation!");
    return false;
  };

  //! Assign pore pressure
  //! \param[in] pressure Pore liquid pressure
  virtual void assign_pore_pressure(double pressure) {
    throw std::runtime_error(
        "Calling the base class function (assign_pore_pressure) in "
        "ParticleBase:: illegal operation!");
  };

  //! Assign liquid traction
  //! \param[in] direction Index corresponding to the direction of traction
  //! \param[in] traction Particle traction in specified direction
  //! \retval status Assignment status
  virtual bool assign_liquid_traction(unsigned direction, double traction) {
    throw std::runtime_error(
        "Calling the base class function (assign_liquid_traction) in "
        "ParticleBase:: illegal operation!");
    return false;
  };

  //! Return liquid phase traction
  virtual VectorDim liquid_traction() const {
    auto error = VectorDim::Zero();
    throw std::runtime_error(
        "Calling the base class function (liquid_traction) in "
        "ParticleBase:: illegal operation!");
    return error;
  };

  //! Return liquid mass
  //! \retval liquid mass Liquid phase mass
  virtual double liquid_mass() const {
    throw std::runtime_error(
        "Calling the base class function (liquid_mass) in "
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
    return false;
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

  //! Return strain of the particle liquid phase
  //! \retval liquid strain Liquid phase strain
  virtual Eigen::Matrix<double, 6, 1> liquid_strain() const {
    auto error = Eigen::Matrix<double, 6, 1>::Zero();
    throw std::runtime_error(
        "Calling the base class function (liquid_strain) in "
        "ParticleBase:: illegal operation!");
    return error;
  };

  //! Assign pore pressure to nodes
  virtual void map_pore_pressure_to_nodes() {
    throw std::runtime_error(
        "Calling the base class function (map_pore_pressure_to_nodes) in "
        "ParticleBase:: illegal operation!");
  };

  //! Compute pore pressure somoothening by interpolating nodal pressure
  virtual bool compute_pore_pressure_smoothing() {
    throw std::runtime_error(
        "Calling the base class function (compute_pore_pressure_smoothing) in "
        "ParticleBase:: illegal operation!");
    return false;
  };

  //! Compute pore pressure
  //! \param[in] dt Time step size
  virtual void compute_pore_pressure(double dt) {
    throw std::runtime_error(
        "Calling the base class function (compute_pore_pressure) in "
        "ParticleBase:: illegal operation!");
  };

  //! Return pore pressure
  //! \retval pore pressure Pore liquid pressure
  virtual double pore_pressure() const {
    throw std::runtime_error(
        "Calling the base class function (pore_pressure) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  //! Return free surface status
  //! \retval free surface Free surface status
  virtual bool free_surface() const {
    throw std::runtime_error(
        "Calling the base class function (free_surface) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  //! Return excessive pore pressure
  //! \retval excessive pore pressure Excessive pore pressure
  virtual double excessive_pore_pressure() const {
    throw std::runtime_error(
        "Calling the base class function (excessive_pore_pressure) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };

  //! Update particle permeability
  virtual VectorDim update_permeability() {
    auto error = VectorDim::Zero();
    throw std::runtime_error(
        "Calling the base class function (update_permeability) in "
        "ParticleBase:: illegal operation!");
    return error;
  };

  //! Update porosity
  virtual bool update_porosity(double dt) {
    throw std::runtime_error(
        "Calling the base class function (update_porosity) in "
        "ParticleBase:: illegal operation!");
    return false;
  };

  //! Map drag force coefficient
  virtual bool map_drag_force_coefficient() {
    throw std::runtime_error(
        "Calling the base class function (map_drag_force_coefficient) in "
        "ParticleBase:: illegal operation!");
    return false;
  };

  //! Assign particle liquid phase velocity constraints
  virtual bool assign_particle_liquid_velocity_constraint(unsigned dir,
                                                          double velocity) {
    throw std::runtime_error(
        "Calling the base class function "
        "(assign_particle_liquid_velocity_constraint) in "
        "ParticleBase:: illegal operation!");
    return false;
  };

  //! Apply particle liquid phase velocity constraints
  virtual void apply_particle_liquid_velocity_constraints() {
    throw std::runtime_error(
        "Calling the base class function "
        "(apply_particle_liquid_velocity_constraints) in "
        "ParticleBase:: illegal operation!");
  };

  //! Initialise particle pore pressure by watertable
  virtual bool initialise_pore_pressure_watertable(
      const unsigned dir_v, const unsigned dir_h,
      std::map<double, double>& refernece_points) {
    throw std::runtime_error(
        "Calling the base class function "
        "(initialise_pore_pressure_watertable) in "
        "ParticleBase:: illegal operation!");
    return false;
  };

  // ---------------------------------------------------------------------------

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
  std::shared_ptr<Material<Tdim>> material_;
  //! Unsigned material id
  unsigned material_id_{std::numeric_limits<unsigned>::max()};
  //! Material state history variables
  mpm::dense_map state_variables_;
  //! Vector of particle neighbour ids
  std::vector<mpm::Index> neighbours_;
};  // ParticleBase class
}  // namespace mpm

#include "particle_base.tcc"

#endif  // MPM_PARTICLEBASE_H__
