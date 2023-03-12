#ifndef MPM_MPM_SCHEME_USL_H_
#define MPM_MPM_SCHEME_USL_H_

#ifdef USE_GRAPH_PARTITIONING
#include "graph.h"
#endif

#include "mpm_scheme.h"

namespace mpm {

//! MPMSchemeUSL class
//! \brief MPMSchemeUSL Derived class for USL stress update scheme
//! \tparam Tdim Dimension
template <unsigned Tdim>
class MPMSchemeUSL : public MPMScheme<Tdim> {
 public:
  //! Default constructor with mesh class
  MPMSchemeUSL(const std::shared_ptr<mpm::Mesh<Tdim>>& mesh, double dt);

  //! Precompute stress
  //! \param[in] phase Phase to smooth postssure
  //! \param[in] postssure_smoothing Enable or disable postssure smoothing
  //! \param[in] interface Enable or disable interface
  virtual inline void precompute_stress_strain(unsigned phase,
                                               bool pressure_smoothing,
                                               bool interface) override;
  //! Postcompute stress
  //! \param[in] phase Phase to smooth postssure
  //! \param[in] postssure_smoothing Enable or disable postssure smoothing
  //! \param[in] interface Enable or disable interface
  virtual inline void postcompute_stress_strain(unsigned phase,
                                                bool pressure_smoothing,
                                                bool interface) override;

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

};  // MPMSchemeUSL class
}  // namespace mpm

#include "mpm_scheme_usl.tcc"

#endif  // MPM_MPM_SCHEME_USL_H_
