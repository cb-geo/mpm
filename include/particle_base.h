#ifndef MPM_PARTICLEBASE_H_
#define MPM_PARTICLEBASE_H_

#include <array>
#include <limits>
#include <memory>
#include <vector>

#include "cell.h"
#include "hdf5.h"
#include "material/material.h"

namespace mpm {

// Forward declaration of Material
template <unsigned Tdim>
class Material;

//! Global index type for the particleBase
using Index = unsigned long long;

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

  //! Assign cell id
  virtual bool assign_cell_id(Index id) = 0;

  //! Return cell id
  virtual Index cell_id() const = 0;

  //! Return cell ptr status
  virtual bool cell_ptr() const = 0;

  //! Remove cell
  virtual void remove_cell() = 0;

  //! Compute shape functions
  virtual bool compute_shapefn() = 0;

  //! Assign volume
  virtual bool assign_volume(unsigned phase, double volume) = 0;

  //! Return volume
  virtual double volume(unsigned phase) const = 0;

  //! Return size of particle in natural coordinates
  virtual VectorDim natural_size() const = 0;

  //! Compute volume of particle
  virtual bool compute_volume(unsigned phase) = 0;

  //! Update volume based on centre volumetric strain rate
  virtual bool update_volume_strainrate(unsigned phase, double dt) = 0;

  //! Compute mass of particle
  virtual bool compute_mass(unsigned phase) = 0;

  //! Map particle mass and momentum to nodes
  virtual bool map_mass_momentum_to_nodes(unsigned phase) = 0;

  // Assign material
  virtual bool assign_material(
      const std::shared_ptr<Material<Tdim>>& material) = 0;

  //! Assign status
  void assign_status(bool status) { status_ = status; }

  //! Status
  bool status() const { return status_; }

  //! Initialise properties
  virtual void initialise() = 0;

  //! Assign mass
  virtual void assign_mass(unsigned phase, double mass) = 0;

  //! Return mass
  virtual double mass(unsigned phase) const = 0;

  //! Return pressure
  virtual double pressure(unsigned phase) const = 0;

  //! Compute strain
  virtual void compute_strain(unsigned phase, double dt) = 0;

  //! Strain
  virtual Eigen::Matrix<double, 6, 1> strain(unsigned phase) const = 0;

  //! Strain rate
  virtual Eigen::Matrix<double, 6, 1> strain_rate(unsigned phase) const = 0;

  //! Volumetric strain of centroid
  virtual double volumetric_strain_centroid(unsigned phase) const = 0;

  //! Initial stress
  virtual void initial_stress(unsigned phase,
                              const Eigen::Matrix<double, 6, 1>&) = 0;

  //! Compute stress
  virtual bool compute_stress(unsigned phase) = 0;

  //! Return stress
  virtual Eigen::Matrix<double, 6, 1> stress(unsigned phase) const = 0;

  //! Map body force
  virtual void map_body_force(unsigned phase, const VectorDim& pgravity) = 0;

  //! Map internal force
  virtual bool map_internal_force(unsigned phase) = 0;

  //! Update pressure of the particles
  virtual bool update_pressure(unsigned phase, double dvolumetric_strain) = 0;

  //! Map particle pressure to nodes
  virtual bool map_pressure_to_nodes(unsigned phase) = 0;

  //! Compute pressure smoothing of the particle based on nodal pressure
  virtual bool compute_pressure_smoothing(unsigned phase) = 0;

  //! Assign velocity
  virtual bool assign_velocity(unsigned phase, const VectorDim& velocity) = 0;

  //! Return velocity
  virtual VectorDim velocity(unsigned phase) const = 0;

  //! Assign traction
  virtual bool assign_traction(unsigned phase, unsigned direction,
                               double traction) = 0;

  //! Return traction
  virtual VectorDim traction(unsigned phase) const = 0;

  //! Map traction force
  virtual void map_traction_force(unsigned phase) = 0;

  //! Compute updated position
  virtual bool compute_updated_position(unsigned phase, double dt) = 0;

  //! Compute updated position based on nodal velocity
  virtual bool compute_updated_position_velocity(unsigned phase, double dt) = 0;

  //! Return a state variable
  virtual double state_variable(const std::string& var) const = 0;

  //! Assign particle velocity constraint
  //! Directions can take values between 0 and Dim * Nphases
  //! \param[in] dir Direction of particle velocity constraint
  //! \param[in] velocity Applied particle velocity constraint
  virtual bool assign_particle_velocity_constraint(unsigned dir,
                                                   double velocity) = 0;

  //! Apply particle velocity constraints
  virtual void apply_particle_velocity_constraints() = 0;

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
  //! Material
  std::shared_ptr<Material<Tdim>> material_;
  //! Material state history variables
  mpm::dense_map state_variables_;
};  // ParticleBase class
}  // namespace mpm

#include "particle_base.tcc"

#endif  // MPM_PARTICLEBASE_H__
