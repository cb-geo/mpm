#ifndef _MPM_STRESS_UPDATE_USL_
#define _MPM_STRESS_UPDATE_USL_

#ifdef USE_GRAPH_PARTITIONING
#include "graph.h"
#endif

#include "stress_update.h"

namespace mpm {

//! StressUpdateUSL class
//! \brief StressUpdateUSL Derived class for USL stress update scheme
//! \tparam Tdim Dimension
template <unsigned Tdim>
class StressUpdateUSL : public StressUpdate<Tdim> {
 public:
  //! Default constructor with mesh class
  StressUpdateUSL(const std::shared_ptr<mpm::Mesh<Tdim>>& mesh, double dt);

  //! Precompute stress
  //! \param[in] phase Phase to smooth postssure
  //! \param[in] postssure_smoothing Enable or disable postssure smoothing
  virtual void precompute_stress_strain(unsigned phase,
                                        bool pressure_smoothing) override;
  //! Postcompute stress
  //! \param[in] phase Phase to smooth postssure
  //! \param[in] postssure_smoothing Enable or disable postssure smoothing
  virtual void postcompute_stress_strain(unsigned phase,
                                         bool pressure_smoothing) override;

  //! Stress update scheme
  //! \retval scheme Stress update scheme
  virtual std::string scheme() const override;

 protected:
  //! Mesh object
  using mpm::StressUpdate<Tdim>::mesh_;
  //! MPI Size
  using mpm::StressUpdate<Tdim>::mpi_size_;
  //! MPI rank
  using mpm::StressUpdate<Tdim>::mpi_rank_;
  //! Time increment
  using mpm::StressUpdate<Tdim>::dt_;

};  // StressUpdateUSL class
}  // namespace mpm

#include "stress_update_usl.tcc"

#endif  // _MPM_STRESS_UPDATE_USL_
