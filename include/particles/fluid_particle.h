#ifndef MPM_FLUID_PARTICLE_H_
#define MPM_FLUID_PARTICLE_H_

#include <array>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "logger.h"
#include "particle.h"

namespace mpm {

//! Fluid Particle class
//! \brief Class with function specific to fluid particles
//! \tparam Tdim Dimension
template <unsigned Tdim>
class FluidParticle : public mpm::Particle<Tdim> {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Construct a particle with id and coordinates
  //! \param[in] id Particle id
  //! \param[in] coord Coordinates of the particles
  FluidParticle(Index id, const VectorDim& coord);

  //! Destructor
  ~FluidParticle() override{};

  //! Delete copy constructor
  FluidParticle(const FluidParticle<Tdim>&) = delete;

  //! Delete assignment operator
  FluidParticle& operator=(const FluidParticle<Tdim>&) = delete;

  //! Compute both solid and liquid mass
  void compute_mass() noexcept override;

  //! Compute stress
  void compute_stress() noexcept override;

  //! Map internal force
  inline void map_internal_force() noexcept override;

  // FIXME: NAMING ERROR
  //! Initial pressure
  //! \param[in] pressure Initial pressure
  void initial_pore_pressure(double pressure) override {
    state_variables_.at("pressure") = pressure;
  }

  // //! Map drag force coefficient
  // bool map_drag_force_coefficient() override;

  // //! Update particle permeability
  // //! \retval status Update status
  // bool update_permeability() override;

  // //! Assign particle pressure constraints
  // //! \retval status Assignment status
  // bool assign_particle_pressure_constraint(double pressure) override;

  //------------------------------------------------------------
  //! Semi-implict mpm
  //! Assigning beta parameter to particle
  void assign_semi_implicit_param(double parameter) override {
    this->beta_ = parameter;
  };

  //! Map laplacian element matrix to cell
  bool map_L_to_cell() override;

  //! Map F element matrix to cell (used in poisson equation RHS)
  bool map_F_to_cell() override;

  //! Map K_cor element matrix to cell
  bool map_K_cor_to_cell() override;

  bool compute_updated_pressure() override;

 private:
  // //! Assign particle permeability
  // //! \retval status Assignment status
  // virtual bool assign_permeability();

  //! Compute stress
  virtual Eigen::Matrix<double, 6, 1> compute_turbulent_stress();

 protected:
  //! Cell
  using ParticleBase<Tdim>::cell_;
  //! Nodes
  using ParticleBase<Tdim>::nodes_;
  //! State variables
  using ParticleBase<Tdim>::state_variables_;
  //! Shape functions
  using Particle<Tdim>::shapefn_;
  //! dN/dX
  using Particle<Tdim>::dn_dx_;
  //! Fluid strain rate
  using Particle<Tdim>::strain_rate_;
  //! Fluid material
  using Particle<Tdim>::material_;
  //! Effective stress of soil skeleton
  using Particle<Tdim>::stress_;
  //! Particle mass density
  using Particle<Tdim>::mass_density_;
  //! Particle mass density
  using Particle<Tdim>::mass_;
  //! Particle total volume
  using Particle<Tdim>::volume_;
  //! Porosity
  using Particle<Tdim>::porosity_;
  //! Beta parameter for semi-implicit update
  double beta_{0.0};
  //! Pressure constraint
  double pressure_constraint_{std::numeric_limits<unsigned>::max()};
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
  //! ---------------------------------------------------------------------
  //! semi-implicit
  //! Degree of saturation
  double liquid_saturation_{1.0};
  //! Permeability
  VectorDim permeability_;
  //! Permeability parameter
  VectorDim c1_;
};  // FluidParticle class
}  // namespace mpm

#include "fluid_particle.tcc"

#endif  // MPM_FLUID_PARTICLE_H__
