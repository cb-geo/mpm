#ifndef MPM_ASSEMBLER_PARALLEL_SEMI_IMPLICIT_NAVIERSTOKES_H_
#define MPM_ASSEMBLER_PARALLEL_SEMI_IMPLICIT_NAVIERSTOKES_H_

#include <Eigen/Sparse>
#include <string>

// Speed log
#include "assembler_base.h"
#include "spdlog/spdlog.h"

#include "cg_eigen.h"
#include "mesh.h"

namespace mpm {
template <unsigned Tdim>
class AssemblerParallelSemiImplicitNavierStokes : public AssemblerBase<Tdim> {
 public:
  //! Constructor
  AssemblerParallelSemiImplicitNavierStokes();

  //! Create a pair between nodes and index in Matrix / Vector
  bool assign_global_node_indices(unsigned nactive_node, unsigned nglobal_active_node) override;

  //! Return laplacian matrix
  Eigen::SparseMatrix<double>& laplacian_matrix() override {
    return laplacian_matrix_;
  }

  //! Assemble laplacian matrix
  bool assemble_laplacian_matrix(double dt) override;

  //! Return poisson RHS vector
  Eigen::VectorXd& poisson_rhs_vector() override { return poisson_rhs_vector_; }

  //! Assemble poisson RHS vector
  bool assemble_poisson_right(double dt) override;

  //! Assign free surface node id
  void assign_free_surface(
      const std::set<mpm::Index>& free_surface_id) override {
    free_surface_ = free_surface_id;
  }

  //! Assign pressure constraints
  bool assign_pressure_constraints(double beta,
                                   const double current_time) override;

  //! Apply pressure constraints to poisson equation
  void apply_pressure_constraints();

  //! Return pressure increment
  Eigen::VectorXd& pressure_increment() override { return pressure_increment_; }

  //! Assign pressure increment
  void assign_pressure_increment(Eigen::VectorXd pressure_increment) override {
    pressure_increment_ = pressure_increment;
  }

  //! Return correction matrix
  Eigen::SparseMatrix<double>& correction_matrix() override {
    return correction_matrix_;
  }

  //! Assemble corrector RHS
  bool assemble_corrector_right(double dt) override;

  unsigned global_active_dof() override {return global_active_dof_;};

  std::vector<mpm::Index> rank_global_mapper() override {return rank_global_mapper_;};

 protected:
  //! Logger
  std::shared_ptr<spdlog::logger> console_;

 private:
  //! number of nodes
  using AssemblerBase<Tdim>::active_dof_;
  //! Mesh object
  using AssemblerBase<Tdim>::mesh_;
  //! Number of total active_dof in all rank
  unsigned global_active_dof_;
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
  //! correction_matrix
  Eigen::SparseMatrix<double> correction_matrix_;
  //! Global node indices
  std::vector<Eigen::VectorXi> global_node_indices_;
  //! Rank to Global mapper
  std::vector<mpm::Index> rank_global_mapper_;
};
}  // namespace mpm

#include "assembler_parallel_semi_implicit_navierstokes.tcc"
#endif  // MPM_ASSEMBLER_PARALLEL_SEMI_IMPLICIT_NAVIERSTOKES_H_
