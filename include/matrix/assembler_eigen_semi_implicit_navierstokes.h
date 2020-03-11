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

  //! Resize containers of matrix
  bool resize_semi_implicit_matrix();

  //! Assemble K_cor matrix (used in correcting nodal velocity)
  bool assemble_K_cor_matrix(std::shared_ptr<mpm::Mesh<Tdim>>& mesh_,
                             double dt) override;

  //! Assemble_laplacian_matrix
  bool assemble_laplacian_matrix(double dt) override;

  //! Assemble_poisson_right
  bool assemble_poisson_right(std::shared_ptr<mpm::Mesh<Tdim>>& mesh_,
                              double dt) override;

  void assign_pressure_increment(Eigen::VectorXd pressure_increment) override {
    pressure_increment_ = pressure_increment;
  }

  Eigen::SparseMatrix<double>& K_cor_matrix() override { return K_cor_matrix_; }

  Eigen::SparseMatrix<double>& laplacian_matrix() override {
    return laplacian_matrix_;
  }

  Eigen::VectorXd& force_laplacian_matrix() override {
    return force_laplacian_matrix_;
  }

  Eigen::VectorXd& pressure_increment() override { return pressure_increment_; }

  std::set<mpm::Index> free_surface() override { return free_surface_; }

  void assign_free_surface(
      const std::set<mpm::Index>& free_surface_id) override {
    free_surface_ = free_surface_id;
  }

  bool assign_pressure_constraints(double beta,
                                   const double current_time) override;

  void apply_pressure_constraints();

 protected:
  //! Logger
  std::shared_ptr<spdlog::logger> console_;

 private:
  //! number of nodes
  using AssemblerBase<Tdim>::active_dof_;
  //! Mesh object
  using AssemblerBase<Tdim>::mesh_;
  //! poisson equation RHS (F31 and F32)
  Eigen::VectorXd force_laplacian_matrix_;
  //! Laplacian matrix
  Eigen::SparseMatrix<double> laplacian_matrix_;
  //! K_cor_matrix
  Eigen::SparseMatrix<double> K_cor_matrix_;
  //! p^(t+1) - beta * p^(t)
  Eigen::VectorXd pressure_increment_;
  //! Laplacian coefficient
  Eigen::VectorXd poisson_right_vector_;
  //! Solver base
  std::shared_ptr<mpm::SolverBase<Tdim>> solver_;
  //! Global node indices
  std::vector<Eigen::VectorXi> global_node_indices_;
  //! Velocity constraints
  Eigen::SparseMatrix<double> velocity_constraints_;
  //! Pressure constraints
  Eigen::SparseVector<double> pressure_constraints_;
  //! Free surface
  std::set<mpm::Index> free_surface_;
};
}  // namespace mpm

#include "assembler_eigen_semi_implicit_navierstokes.tcc"
#endif  // MPM_ASSEMBLER_EIGEN_SEMI_IMPLICIT_NAVIERSTOKES_H_
