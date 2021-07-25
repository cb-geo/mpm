#ifndef MPM_MPM_SCHEME_MUSL_H_
#define MPM_MPM_SCHEME_MUSL_H_

#ifdef USE_GRAPH_PARTITIONING
#include "graph.h"
#endif

#include "mpm_scheme.h"

namespace mpm {

//! MPMSchemeMUSL class
//! \brief MPMSchemeUSL Derived class for MUSL stress update scheme
//! \tparam Tdim Dimension
template <unsigned Tdim>
class MPMSchemeMUSL : public MPMScheme<Tdim> {
 public:
  //! Default constructor with mesh class
  MPMSchemeMUSL(const std::shared_ptr<mpm::Mesh<Tdim>>& mesh, double dt);

  //! Precompute stress
  //! \param[in] phase Phase to smooth postssure
  //! \param[in] postssure_smoothing Enable or disable postssure smoothing
  virtual inline void precompute_stress_strain(
      unsigned phase, bool pressure_smoothing) override;
  //! Postcompute stress
  //! \param[in] phase Phase to smooth postssure
  //! \param[in] postssure_smoothing Enable or disable postssure smoothing
  virtual inline void postcompute_stress_strain(
      unsigned phase, bool pressure_smoothing) override;

  //! Compute acceleration velocity position
  //! \param[in] velocity_update Velocity or acceleration update flag
  //! \param[in] phase Phase of particle
  //! \param[in] damping_type Type of damping
  //! \param[in] damping_factor Value of critical damping
  virtual inline void compute_particle_kinematics(
      bool velocity_update, unsigned phase, const std::string& damping_type,
      double damping_factor) override;

  //! Compute position
  //! \param[in] phase Phase of particle
  virtual inline void compute_particle_updated_position(bool velocity_update, unsigned phase) override;

  //! Stress update scheme
  //! \retval scheme Stress update scheme
  virtual inline std::string scheme() const override;

 protected:
  //! Mesh object
  using mpm::MPMScheme<Tdim>::mesh_;
  //! MPI Size
  using mpm::MPMScheme<Tdim>::mpi_size_;
  //! MPI rank
  using mpm::MPMScheme<Tdim>::mpi_rank_;
  //! Time increment
  using mpm::MPMScheme<Tdim>::dt_;

};  // MPMSchemeMUSL class
}  // namespace mpm

#include "mpm_scheme_musl.tcc"

#endif  // MPM_MPM_SCHEME_MUSL_H_
