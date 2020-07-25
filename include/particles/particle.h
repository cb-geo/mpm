#ifndef MPM_PARTICLE_H_
#define MPM_PARTICLE_H_

#include <array>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "cell.h"
#include "logger.h"
#include "particle_base.h"

namespace mpm {

//! Particle class
//! \brief Base class that stores the information about particles
//! \details Particle class: id_ and coordinates.
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Particle : public ParticleBase<Tdim> {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Define DOFs
  static const unsigned Tdof = (Tdim == 1) ? 1 : 3 * (Tdim - 1);

  //! Construct a particle with id and coordinates
  //! \param[in] id Particle id
  //! \param[in] coord coordinates of the particle
  Particle(Index id, const VectorDim& coord);

  //! Construct a particle with id, coordinates and status
  //! \param[in] id Particle id
  //! \param[in] coord coordinates of the particle
  //! \param[in] status Particle status (active / inactive)
  Particle(Index id, const VectorDim& coord, bool status);

  //! Destructor
  ~Particle() override{};

  //! Delete copy constructor
  Particle(const Particle<Tdim>&) = delete;

  //! Delete assignment operator
  Particle& operator=(const Particle<Tdim>&) = delete;

  //! Initialise particle from HDF5 data
  //! \param[in] particle HDF5 data of particle
  //! \retval status Status of reading HDF5 particle
  bool initialise_particle(const HDF5Particle& particle) override;

  //! Initialise particle HDF5 data and material
  //! \param[in] particle HDF5 data of particle
  //! \param[in] material Material associated with the particle
  //! \retval status Status of reading HDF5 particle
  bool initialise_particle(
      const HDF5Particle& particle,
      const std::shared_ptr<Material<Tdim>>& material) override;

  //! Assign material history variables
  //! \param[in] state_vars State variables
  //! \param[in] material Material associated with the particle
  //! \retval status Status of cloning HDF5 particle
  bool assign_material_state_vars(
      const mpm::dense_map& state_vars,
      const std::shared_ptr<mpm::Material<Tdim>>& material) override;

  //! Retrun particle data as HDF5
  //! \retval particle HDF5 data of the particle
  HDF5Particle hdf5() const override;

  //! Initialise properties
  void initialise() override;

  //! Compute reference coordinates in a cell
  bool compute_reference_location() noexcept override;

  //! Return reference location
  VectorDim reference_location() const override { return xi_; }

  //! Assign a cell to particle
  //! If point is in new cell, assign new cell and remove particle id from old
  //! cell. If point can't be found in the new cell, check if particle is still
  //! valid in the old cell, if it is leave it as is. If not, set cell as null
  //! \param[in] cellptr Pointer to a cell
  bool assign_cell(const std::shared_ptr<Cell<Tdim>>& cellptr) override;

  //! Assign a cell to particle
  //! If point is in new cell, assign new cell and remove particle id from old
  //! cell. If point can't be found in the new cell, check if particle is still
  //! valid in the old cell, if it is leave it as is. If not, set cell as null
  //! \param[in] cellptr Pointer to a cell
  //! \param[in] xi Local coordinates of the point in reference cell
  bool assign_cell_xi(const std::shared_ptr<Cell<Tdim>>& cellptr,
                      const Eigen::Matrix<double, Tdim, 1>& xi) override;

  //! Assign cell id
  //! \param[in] id Cell id
  bool assign_cell_id(Index id) override;

  //! Return cell id
  Index cell_id() const override { return cell_id_; }

  //! Return cell ptr
  std::shared_ptr<Cell<Tdim>> cell() const override { return cell_; }

  //! Return cell ptr status
  bool cell_ptr() const override { return cell_ != nullptr; }

  //! Remove cell associated with the particle
  void remove_cell() override;

  //! Compute shape functions of a particle, based on local coordinates
  void compute_shapefn() noexcept override;

  //! Assign volume
  //! \param[in] volume Volume of particle
  bool assign_volume(double volume) override;

  //! Return volume
  double volume() const override {
    return this->scalar_property(mpm::properties::Scalar::Volume);
  }

  //! Return size of particle in natural coordinates
  VectorDim natural_size() const override { return natural_size_; }

  //! Return mass density
  double mass_density() const override {
    return this->scalar_property(mpm::properties::Scalar::MassDensity);
  }

  //! Map multimaterial properties to nodes
  void map_multimaterial_mass_momentum_to_nodes() noexcept override;

  //! Map multimaterial displacements to nodes
  void map_multimaterial_displacements_to_nodes() noexcept override;

  //! Map multimaterial domain gradients to nodes
  void map_multimaterial_domain_gradients_to_nodes() noexcept override;

  //! Assign nodal mass to particles
  //! \param[in] mass Mass from the particles in a cell
  void assign_mass(double mass) override {
    scalar_properties_.at(mpm::properties::Scalar::Mass) = mass;
  }

  //! Return mass of the particles
  double mass() const override {
    return this->scalar_property(mpm::properties::Scalar::Mass);
  }

  //! Assign material
  //! \param[in] material Pointer to a material
  bool assign_material(
      const std::shared_ptr<Material<Tdim>>& material) override;

  //! Return material
  //! \retval material Pointer to a material
  std::shared_ptr<Material<Tdim>> material() const override;

  //! Compute strain
  //! \param[in] dt Analysis time step
  void compute_strain(double dt) noexcept override;

  //! Return strain of the particle
  Eigen::Matrix<double, 6, 1> strain() const override { return strain_; }

  //! Return strain rate of the particle
  Eigen::Matrix<double, 6, 1> strain_rate() const override {
    return strain_rate_;
  };

  //! Return dvolumetric strain of centroid
  //! \retval dvolumetric strain at centroid
  double dvolumetric_strain() const override { return dvolumetric_strain_; }

  //! Return volumetric strain of centroid
  //! \retval volumetric strain at centroid
  double volumetric_strain_centroid() const override {
    return volumetric_strain_centroid_;
  }

  //! Initial stress
  //! \param[in] stress Initial sress
  void initial_stress(const Eigen::Matrix<double, 6, 1>& stress) override {
    this->stress_ = stress;
  }

  //! Compute stress
  void compute_stress() noexcept override;

  //! Return stress of the particle
  Eigen::Matrix<double, 6, 1> stress() const override { return stress_; }

  //! Map internal force
  inline void map_internal_force() noexcept override;

  //! Assign velocity to the particle
  //! \param[in] velocity A vector of particle velocity
  void assign_velocity(const VectorDim& velocity) override {
    vector_properties_.at(mpm::properties::Vector::Velocity) = velocity;
  };

  //! Return velocity of the particle
  VectorDim velocity() const override {
    return this->vector_property(mpm::properties::Vector::Velocity);
  }

  //! Return displacement of the particle
  VectorDim displacement() const override {
    return this->vector_property(mpm::properties::Vector::Displacement);
  }

  //! Assign traction to the particle
  //! \param[in] direction Index corresponding to the direction of traction
  //! \param[in] traction Particle traction in specified direction
  //! \retval status Assignment status
  bool assign_traction(unsigned direction, double traction) override;

  //! Return traction of the particle
  //! \param[in] phase Index corresponding to the phase
  VectorDim traction() const override { return traction_; }

  //! Compute updated position of the particle
  //! \param[in] dt Analysis time step
  //! \param[in] velocity_update Update particle velocity from nodal vel
  void compute_updated_position(double dt,
                                bool velocity_update = false) noexcept override;

  //! Assign a state variable
  //! \param[in] var State variable
  //! \param[in] value State variable to be assigned
  void assign_state_variable(const std::string& var, double value) override {
    if (state_variables_.find(var) != state_variables_.end())
      state_variables_.at(var) = value;
  }

  //! Return a state variable
  //! \param[in] var State variable
  //! \retval Quantity of the state history variable
  double state_variable(const std::string& var) const override {
    return (state_variables_.find(var) != state_variables_.end())
               ? state_variables_.at(var)
               : std::numeric_limits<double>::quiet_NaN();
  }

  //! Assign a state variable
  //! \param[in] value Particle pressure to be assigned
  void assign_pressure(double pressure) override {
    this->assign_state_variable("pressure", pressure);
  }

  //! Return pressure of the particles
  double pressure() const override {
    return (state_variables_.find("pressure") != state_variables_.end())
               ? state_variables_.at("pressure")
               : std::numeric_limits<double>::quiet_NaN();
  }

  //! Return tensor data of particles
  //! \param[in] property Property string
  //! \retval vecdata Tensor data of particle property
  Eigen::VectorXd tensor_data(const std::string& property) override;

  //! Apply particle velocity constraints
  //! \param[in] dir Direction of particle velocity constraint
  //! \param[in] velocity Applied particle velocity constraint
  void apply_particle_velocity_constraints(unsigned dir,
                                           double velocity) override;

  //! Assign material id of this particle to nodes
  void append_material_id_to_nodes() const override;

  //! Return the number of neighbour particles
  unsigned nneighbours() const override { return neighbours_.size(); };

  //! Assign neighbour particles
  //! \param[in] neighbours set of id of the neighbouring particles
  //! \retval insertion_status Return the successful addition of a node
  void assign_neighbours(const std::vector<mpm::Index>& neighbours) override;

  //! Return neighbour ids
  std::vector<mpm::Index> neighbours() const override { return neighbours_; };

 private:
  //! Compute strain rate
  inline Eigen::Matrix<double, 6, 1> compute_strain_rate(
      const Eigen::MatrixXd& dn_dx, unsigned phase) noexcept;

 private:
  //! particle id
  using ParticleBase<Tdim>::id_;
  //! coordinates
  using ParticleBase<Tdim>::coordinates_;
  //! Reference coordinates (in a cell)
  using ParticleBase<Tdim>::xi_;
  //! Cell
  using ParticleBase<Tdim>::cell_;
  //! Cell id
  using ParticleBase<Tdim>::cell_id_;
  //! Nodes
  using ParticleBase<Tdim>::nodes_;
  //! Status
  using ParticleBase<Tdim>::status_;
  //! Material
  using ParticleBase<Tdim>::material_;
  //! Material id
  using ParticleBase<Tdim>::material_id_;
  //! State variables
  using ParticleBase<Tdim>::state_variables_;
  //! Neighbour particles
  using ParticleBase<Tdim>::neighbours_;
  //! Scalar properties
  using ParticleBase<Tdim>::boolean_properties_;
  //! Scalar properties
  using ParticleBase<Tdim>::scalar_properties_;
  //! Vector properties
  using ParticleBase<Tdim>::vector_properties_;
  //! Shape functions
  using ParticleBase<Tdim>::shapefn_;
  //! Size of particle
  Eigen::Matrix<double, 1, Tdim> size_;
  //! Size of particle in natural coordinates
  Eigen::Matrix<double, 1, Tdim> natural_size_;
  //! Stresses
  Eigen::Matrix<double, 6, 1> stress_;
  //! Strains
  Eigen::Matrix<double, 6, 1> strain_;
  //! dvolumetric strain
  double dvolumetric_strain_{0.};
  //! Volumetric strain at centroid
  double volumetric_strain_centroid_{0.};
  //! Strain rate
  Eigen::Matrix<double, 6, 1> strain_rate_;
  //! dstrains
  Eigen::Matrix<double, 6, 1> dstrain_;
  //! Particle velocity constraints
  std::map<unsigned, double> particle_velocity_constraints_;
  //! Surface Traction (given as a stress; force/area)
  Eigen::Matrix<double, Tdim, 1> traction_;
  //! dN/dX
  Eigen::MatrixXd dn_dx_;
  //! dN/dX at cell centroid
  Eigen::MatrixXd dn_dx_centroid_;
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
  //! Map of vector properties
  std::map<std::string, std::function<Eigen::VectorXd()>> properties_;

};  // Particle class
}  // namespace mpm

#include "particle.tcc"
#include "particle_functions.tcc"

#endif  // MPM_PARTICLE_H__
