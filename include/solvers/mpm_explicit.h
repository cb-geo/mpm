#ifndef MPM_MPM_EXPLICIT_H_
#define MPM_MPM_EXPLICIT_H_

#ifdef USE_GRAPH_PARTITIONING
#include "graph.h"
#endif

#include "mpm_base.h"

namespace mpm {

//! MPMExplicit class
//! \brief A class that implements the fully explicit one phase mpm
//! \details A single-phase explicit MPM
//! \tparam Tdim Dimension
template <unsigned Tdim>
class MPMExplicit : public MPMBase<Tdim> {
 public:
  //! Default constructor
  MPMExplicit(std::unique_ptr<IO>&& io);

  //! Domain decomposition
  void mpi_domain_decompose();

  //! Solve
  bool solve() override;

  //! Pressure smoothing
  //! \param[in] phase Phase to smooth pressure
  void pressure_smoothing(unsigned phase);

  //! Compute stress strain
  //! \param[in] phase Phase to smooth pressure
  void compute_stress_strain(unsigned phase);

 protected:
  // Generate a unique id for the analysis
  using mpm::MPMBase<Tdim>::uuid_;
  //! Time step size
  using mpm::MPMBase<Tdim>::dt_;
  //! Current step
  using mpm::MPMBase<Tdim>::step_;
  //! Number of steps
  using mpm::MPMBase<Tdim>::nsteps_;
  //! Output steps
  using mpm::MPMBase<Tdim>::output_steps_;
  //! A unique ptr to IO object
  using mpm::MPMBase<Tdim>::io_;
  //! JSON analysis object
  using mpm::MPMBase<Tdim>::analysis_;
  //! JSON post-process object
  using mpm::MPMBase<Tdim>::post_process_;
  //! Logger
  using mpm::MPMBase<Tdim>::console_;

#ifdef USE_GRAPH_PARTITIONING
  //! Graph
  using mpm::MPMBase<Tdim>::graph_;
#endif

  //! velocity update
  using mpm::MPMBase<Tdim>::velocity_update_;
  //! Gravity
  using mpm::MPMBase<Tdim>::gravity_;
  //! Mesh object
  using mpm::MPMBase<Tdim>::mesh_;
  //! Materials
  using mpm::MPMBase<Tdim>::materials_;
  //! VTK attributes
  using mpm::MPMBase<Tdim>::vtk_attributes_;
  //! Node concentrated force
  using mpm::MPMBase<Tdim>::set_node_concentrated_force_;

 private:
  //! Boolean to switch between USL and USF
  bool usl_{false};
  //! Pressure smoothing
  bool pressure_smoothing_{false};
  //! Interface
  bool interface_{false};

};  // MPMExplicit class
}  // namespace mpm

#include "mpm_explicit.tcc"

#endif  // MPM_MPM_EXPLICIT_H_
