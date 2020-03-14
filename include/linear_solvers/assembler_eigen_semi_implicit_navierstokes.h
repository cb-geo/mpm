#ifndef MPM_ASSEMBLER_EIGEN_SEMI_IMPLICIT_NAVIERSTOKES_H_
#define MPM_ASSEMBLER_EIGEN_SEMI_IMPLICIT_NAVIERSTOKES_H_

#include <Eigen/Sparse>
#include <string>

// Speed log
#include "assembler_base.h"
#include "spdlog/spdlog.h"

#include "cg_eigen.h"
#include "mesh.h"

namespace mpm {
template <unsigned Tdim>
class AssemblerEigenSemiImplicitNavierStokes : public AssemblerBase<Tdim> {
 public:
  //! Constructor
  AssemblerEigenSemiImplicitNavierStokes();

  //! Create a pair between nodes and index in Matrix / Vector
  bool assign_global_node_indices(unsigned active_dof) override;

  //! Return laplacian matrix
  Eigen::SparseMatrix<double>& laplacian_matrix() override {
    return laplacian_matrix_;
  }

  //! Assemble laplacian matrix
  bool assemble_laplacian_matrix(double dt) override;

  //! Return poisson RHS vector
  Eigen::VectorXd& poisson_rhs_vector() override { return poisson_rhs_vector_; }

  //! Assemble poisson RHS vector
  bool assemble_poisson_right(std::shared_ptr<mpm::Mesh<Tdim>>& mesh_,
                              double dt) override;

  //! Assign free surface node id
  void assign_free_surface(
      const std::set<mpm::Index>& free_surface_id) override {
    free_surface_ = free_surface_id;
  }

  //! Return free surface node id
  std::set<mpm::Index> free_surface() override { return free_surface_; }

  //! Apply pressure constraints to poisson equation
  void apply_pressure_constraints();

  //! Assign pressure increment
  void assign_pressure_increment(Eigen::VectorXd pressure_increment) override {
    pressure_increment_ = pressure_increment;
  }

  //! Return pressure increment
  Eigen::VectorXd& pressure_increment() override { return pressure_increment_; }

  Eigen::SparseMatrix<double>& K_cor_matrix() override { return K_cor_matrix_; }

  bool assign_pressure_constraints(double beta,
                                   const double current_time) override;

  //! Assemble K_cor matrix (used in correcting nodal velocity)
  bool assemble_K_cor_matrix(std::shared_ptr<mpm::Mesh<Tdim>>& mesh_,
                             double dt) override;

 protected:
  //! Logger
  std::shared_ptr<spdlog::logger> console_;

 private:
  //! number of nodes
  using AssemblerBase<Tdim>::active_dof_;
  //! Mesh object
  using AssemblerBase<Tdim>::mesh_;
  //! Laplacian matrix
  Eigen::SparseMatrix<double> laplacian_matrix_;
  //! Poisson RHS vector
  Eigen::VectorXd poisson_rhs_vector_;
  //! Free surface
  std::set<mpm::Index> free_surface_;
  //! Pressure constraints
  Eigen::SparseVector<double> pressure_constraints_;
  //! \delta p^(t+1) = p^(t+1) - beta * p^(t)
  Eigen::VectorXd pressure_increment_;

  //! K_cor_matrix
  Eigen::SparseMatrix<double> K_cor_matrix_;
  //! Laplacian coefficient
  Eigen::VectorXd poisson_right_vector_;
  //! Solver base
  std::shared_ptr<mpm::SolverBase<Tdim>> solver_;
  //! Global node indices
  std::vector<Eigen::VectorXi> global_node_indices_;
  //! Velocity constraints
  Eigen::SparseMatrix<double> velocity_constraints_;
};
}  // namespace mpm

#include "assembler_eigen_semi_implicit_navierstokes.tcc"
#endif  // MPM_ASSEMBLER_EIGEN_SEMI_IMPLICIT_NAVIERSTOKES_H_
