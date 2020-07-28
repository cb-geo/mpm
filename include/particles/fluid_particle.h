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

  //! Compute stress
  void compute_stress() noexcept override;

  //! Map internal force
  inline void map_internal_force() noexcept override;

  //! ----------------------------------------------------------------
  //! Semi-Implicit integration functions based on Chorin's Projection
  //! ----------------------------------------------------------------

  //! Assigning beta parameter to particle
  //! \param[in] pressure parameter determining type of projection
  void assign_projection_parameter(double parameter) override {
    this->projection_param_ = parameter;
  };

  //! Map laplacian element matrix to cell (used in poisson equation LHS)
  bool map_laplacian_to_cell() override;

  //! Map poisson rhs element matrix to cell (used in poisson equation RHS)
  bool map_poisson_right_to_cell() override;

  //! Map correction matrix element matrix to cell (used to correct velocity)
  bool map_correction_matrix_to_cell() override;

  //! Update pressure after solving poisson equation
  bool compute_updated_pressure() override;

 private:
  //! Compute turbulent stress
  virtual Eigen::Matrix<double, 6, 1> compute_turbulent_stress();

 private:
  //! Cell
  using ParticleBase<Tdim>::cell_;
  //! Nodes
  using ParticleBase<Tdim>::nodes_;
  //! Fluid material
  using ParticleBase<Tdim>::material_;
  //! State variables
  using ParticleBase<Tdim>::state_variables_;
  //! Shape functions
  using Particle<Tdim>::shapefn_;
  //! dN/dX
  using Particle<Tdim>::dn_dx_;
  //! Fluid strain rate
  using Particle<Tdim>::strain_rate_;
  //! Effective stress of soil skeleton
  using Particle<Tdim>::stress_;
  //! Particle mass density
  using Particle<Tdim>::mass_density_;
  //! Particle mass density
  using Particle<Tdim>::mass_;
  //! Particle total volume
  using Particle<Tdim>::volume_;
  //! Projection parameter for semi-implicit update
  double projection_param_{0.0};
  //! Pressure constraint
  double pressure_constraint_{std::numeric_limits<unsigned>::max()};
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
};  // FluidParticle class
}  // namespace mpm

#include "fluid_particle.tcc"

#endif  // MPM_FLUID_PARTICLE_H__
