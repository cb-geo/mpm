#ifndef MPM_MPM_IMPLICIT_LINEAR_H_
#define MPM_MPM_IMPLICIT_LINEAR_H_

#ifdef USE_GRAPH_PARTITIONING
#include "graph.h"
#endif

#include "mpm_base.h"

#include "assembler_base.h"
#include "solver_base.h"

namespace mpm {

//! MPMImplicitLinear class
//! \brief A class that implements the fully implicit linear one phase mpm
//! \details A single-phase implicit linear MPM
//! \tparam Tdim Dimension
template <unsigned Tdim>
class MPMImplicitLinear : public MPMBase<Tdim> {
 public:
  //! Default constructor
  MPMImplicitLinear(const std::shared_ptr<IO>& io);

  //! Solve
  bool solve() override;

  //! Compute stress strain
  //! \param[in] phase Phase to smooth pressure
  void compute_stress_strain(unsigned phase);

  //! Class private functions
 private:
  /**
   * \defgroup Implicit Functions dealing with implicit MPM
   */
  /**@{*/
  //! Initialise matrix
  //! \ingroup Implicit
  bool initialise_matrix();

  //! Initialise matrix
  //! \ingroup Implicit
  bool reinitialise_matrix();

  //! Compute equilibrium equation
  //! \ingroup Implicit
  bool compute_equilibrium_equation();
  /**@}*/

  //! Class private variables
 private:
  // Generate a unique id for the analysis
  using mpm::MPMBase<Tdim>::uuid_;
  //! Time step size
  using mpm::MPMBase<Tdim>::dt_;
  //! Current step
  using mpm::MPMBase<Tdim>::step_;
  //! Number of steps
  using mpm::MPMBase<Tdim>::nsteps_;
  //! Number of steps
  using mpm::MPMBase<Tdim>::nload_balance_steps_;
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
  //! MPM Scheme
  using mpm::MPMBase<Tdim>::mpm_scheme_;
  //! Stress update method
  using mpm::MPMBase<Tdim>::stress_update_;

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
  //! Node concentrated force
  using mpm::MPMBase<Tdim>::set_node_concentrated_force_;
  //! Damping type
  using mpm::MPMBase<Tdim>::damping_type_;
  //! Damping factor
  using mpm::MPMBase<Tdim>::damping_factor_;
  //! Locate particles
  using mpm::MPMBase<Tdim>::locate_particles_;
  //! Nonlocal neighbourhood
  using mpm::MPMBase<Tdim>::node_neighbourhood_;
  //! Pressure smoothing
  bool pressure_smoothing_{false};

  /**
   * \defgroup ImplicitVariables Variables dealing with implicit MPM
   */
  /**@{*/
  //! Parameter beta of Newmark scheme
  using mpm::MPMBase<Tdim>::newmark_beta_;
  //! Parameter gamma of Newmark scheme
  using mpm::MPMBase<Tdim>::newmark_gamma_;
  //! Assembler object
  std::shared_ptr<mpm::AssemblerBase<Tdim>> assembler_;
  //! Linear solver object
  tsl::robin_map<std::string,
                 std::shared_ptr<mpm::SolverBase<Eigen::SparseMatrix<double>>>>
      linear_solver_;
  /**@}*/

};  // MPMImplicitLinear class
}  // namespace mpm

#include "mpm_implicit_linear.tcc"

#endif  // MPM_MPM_IMPLICIT_LINEAR_H_
