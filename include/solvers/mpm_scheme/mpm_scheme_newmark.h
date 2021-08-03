#ifndef MPM_MPM_SCHEME_NEWMARK_H_
#define MPM_MPM_SCHEME_NEWMARK_H_

#ifdef USE_GRAPH_PARTITIONING
#include "graph.h"
#endif

#include "mpm_scheme.h"

namespace mpm {

//! MPMSchemeNewmark class
//! \brief MPMSchemeNewmark Derived class for Newmark predictor-corrector scheme
//! \tparam Tdim Dimension
template <unsigned Tdim>
class MPMSchemeNewmark : public MPMScheme<Tdim> {
 public:
  //! Default constructor with mesh class
  MPMSchemeNewmark(const std::shared_ptr<mpm::Mesh<Tdim>>& mesh, double dt);

  //! Compute nodal kinematics - map mass, momentum and inertia to nodes
  //! \param[in] phase Phase to smooth pressure
  virtual inline void compute_nodal_kinematics(unsigned phase) override;

  //! Predict nodal kinematics - Predictor step of Newmark scheme
  //! \param[in] newmark_beta Parameter beta of Newmark scheme
  //! \param[in] newmark_gamma Parameter gamma of Newmark scheme
  //! \param[in] phase Phase to smooth pressure
  virtual inline void predict_nodal_kinematics(unsigned phase, double newmark_beta,
                                               double newmark_gamma) override;

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

  //! Compute forces
  //! \param[in] gravity Acceleration due to gravity
  //! \param[in] step Number of step in solver
  //! \param[in] concentrated_nodal_forces Boolean for if a concentrated force
  //! is applied or not
  virtual inline void compute_forces(
      const Eigen::Matrix<double, Tdim, 1>& gravity, unsigned phase,
      unsigned step, bool concentrated_nodal_forces) override;

  //! Postcompute nodal kinematics - map mass and momentum to nodes
  //! \param[in] phase Phase to smooth pressure
  virtual inline void postcompute_nodal_kinematics(unsigned phase) override;

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

};  // MPMSchemeNewmark class
}  // namespace mpm

#include "mpm_scheme_newmark.tcc"

#endif  // MPM_MPM_SCHEME_NEWMARK_H_
