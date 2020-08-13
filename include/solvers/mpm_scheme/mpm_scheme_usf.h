#ifndef MPM_MPM_SCHEME_USF_H_
#define MPM_MPM_SCHEME_USF_H_

#ifdef USE_GRAPH_PARTITIONING
#include "graph.h"
#endif

#include "mpm_scheme.h"

namespace mpm {

//! MPMSchemeUSF class
//! \brief MPMSchemeUSF Derived class for USF stress update scheme
//! \tparam Tdim Dimension
template <unsigned Tdim>
class MPMSchemeUSF : public MPMScheme<Tdim> {
 public:
  //! Default constructor with mesh class
  MPMSchemeUSF(const std::shared_ptr<mpm::Mesh<Tdim>>& mesh, double dt);

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

};  // MPMSchemeUSF class
}  // namespace mpm

#include "mpm_scheme_usf.tcc"

#endif  // MPM_MPM_SCHEME_H_
