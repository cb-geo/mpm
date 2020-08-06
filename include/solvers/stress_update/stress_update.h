#ifndef _MPM_STRESS_UPDATE_
#define _MPM_STRESS_UPDATE_

#ifdef USE_GRAPH_PARTITIONING
#include "graph.h"
#endif

#include "mesh.h"

namespace mpm {

//! StressUpdate class
//! \brief StressUpdate base class to support different stress update schemes
//! \tparam Tdim Dimension
template <unsigned Tdim>
class StressUpdate {
 public:
  //! Default constructor with mesh class
  StressUpdate(const std::shared_ptr<mpm::Mesh<Tdim>>& mesh, double dt);

  //! Intialize
  void initialise();

  //! Map mass and momentum to nodes
  void momentum_nodes(unsigned phase);

  //! Compute stress and strain
  //! \param[in] phase Phase to smooth pressure
  //! \param[in] pressure_smoothing Enable or disable pressure smoothing
  void compute_stress_strain(unsigned phase, bool pressure_smoothing);

  //! Precompute stress and strain (empty call)
  //! \param[in] phase Phase to smooth pressure
  //! \param[in] pressure_smoothing Enable or disable pressure smoothing
  virtual void precompute_stress_strain(unsigned phase,
                                        bool pressure_smoothing) = 0;

  //! Postcompute stress and strain (empty call)
  //! \param[in] phase Phase to smooth postssure
  //! \param[in] postssure_smoothing Enable or disable postssure smoothing
  virtual void postcompute_stress_strain(unsigned phase,
                                         bool pressure_smoothing) = 0;

  //! Pressure smoothing
  //! \param[in] phase Phase to smooth pressure
  //! \param[in] pressure_smoothing Enable or disable pressure smoothing
  void pressure_smoothing(unsigned phase);

  //! Compute forces
  //! \param[in] gravity Acceleration due to gravity
  //! \param[in] step Number of step in solver
  //! \param[in] concentrated_nodal_forces Boolean for if a concentrated force
  //! is applied or not
  void compute_forces(const Eigen::Matrix<double, Tdim, 1>& gravity,
                      unsigned phase, unsigned step,
                      bool concentrated_nodal_forces);

  //! Compute acceleration velocity position
  //! \param[in] velocity_update Velocity or acceleration update flag
  //! \param[in] phase Phase of particle
  //! \param[in] damping_type Type of damping
  //! \param[in] damping_factor Value of critical damping
  void compute_particle_kinematics(bool velocity_update, unsigned phase,
                                   const std::string& damping_type,
                                   double damping_factor);

  //! Compute particle location
  //! \param[in] locate_particles Flag to enable locate particles, if set to
  //! false, unlocated particles will be removed
  void locate_particles(bool locate_particles);

  //! Stress update scheme
  //! \retval scheme Stress update scheme
  virtual std::string scheme() const = 0;

 protected:
  //! Mesh object
  std::shared_ptr<mpm::Mesh<Tdim>> mesh_;
  //! Time increment
  double dt_;
  //! MPI Size
  int mpi_size_ = 1;
  //! MPI rank
  int mpi_rank_ = 0;
};  // StressUpdate class
}  // namespace mpm

#include "stress_update.tcc"

#endif  // _MPM_STRESS_UPDATE_
