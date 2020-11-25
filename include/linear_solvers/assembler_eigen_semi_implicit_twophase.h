#ifndef MPM_ASSEMBLER_EIGEN_SEMI_IMPLICIT_TWOPHASE_H_
#define MPM_ASSEMBLER_EIGEN_SEMI_IMPLICIT_TWOPHASE_H_

#include <Eigen/Sparse>
#include <string>

// Speed log
#include "assembler_eigen_semi_implicit_navierstokes.h"
#include "spdlog/spdlog.h"

#include "cg_eigen.h"
#include "mesh.h"

namespace mpm {
template <unsigned Tdim>
class AssemblerEigenSemiImplicitTwoPhase
    : public AssemblerEigenSemiImplicitNavierStokes<Tdim> {
 public:
  //! Constructor
  AssemblerEigenSemiImplicitTwoPhase();

  //! Assemble coefficient matrix for two-phase predictor
  bool assemble_predictor_left(double dt) override;

  //! Return predictor coefficient LHS matrix
  Eigen::SparseMatrix<double>& predictor_lhs_matrix(unsigned dir) override {
    return predictor_lhs_matrix_.at(dir);
  }

  //! Assemble predictor RHS force vector
  bool assemble_predictor_right(double dt) override;

  //! Return predictor RHS vector
  Eigen::MatrixXd& predictor_rhs_vector() override {
    return predictor_rhs_vector_;
  }

  //! Assign intermediate acceleration
  void assign_intermediate_acceleration(
      unsigned dir, Eigen::VectorXd acceleration_inter) override {
    intermediate_acceleration_.col(dir) = acceleration_inter;
  }

  //! Return intermediate acceleration
  Eigen::MatrixXd& intermediate_acceleration() override {
    return intermediate_acceleration_;
  }

  //! Assemble poisson RHS vector
  bool assemble_poisson_right(double dt) override;

  //! Assemble corrector RHS
  bool assemble_corrector_right(double dt) override;

  //! Assign velocity constraints
  bool assign_velocity_constraints() override;

  //! Apply velocity constraints for matrix and vector
  bool apply_velocity_constraints() override;

  //! Assign pressure constraints
  bool assign_pressure_constraints(double beta, double current_time) override;

 protected:
  //! number of nodes
  using AssemblerBase<Tdim>::active_dof_;
  //! Mesh object
  using AssemblerBase<Tdim>::mesh_;
  //! Logger
  using AssemblerEigenSemiImplicitNavierStokes<Tdim>::console_;
  //! Global node indices
  using AssemblerEigenSemiImplicitNavierStokes<Tdim>::global_node_indices_;
  //! Poisson RHS vector
  using AssemblerEigenSemiImplicitNavierStokes<Tdim>::poisson_rhs_vector_;
  //! Pressure constraints
  using AssemblerEigenSemiImplicitNavierStokes<Tdim>::pressure_constraints_;
  //! Correction_matrix
  using AssemblerEigenSemiImplicitNavierStokes<Tdim>::correction_matrix_;
  //! Coefficient matrix for two-phase predictor
  std::map<unsigned, Eigen::SparseMatrix<double>> predictor_lhs_matrix_;
  //! RHS vector for two-phase predictor
  Eigen::MatrixXd predictor_rhs_vector_;
  //! Intermediate acceleration vector (each column represent one direction)
  Eigen::MatrixXd intermediate_acceleration_;
  //! Velocity constraints
  Eigen::SparseMatrix<double> velocity_constraints_;
};
}  // namespace mpm

#include "assembler_eigen_semi_implicit_twophase.tcc"
#endif  // MPM_ASSEMBLER_EIGEN_SEMI_IMPLICIT_TWOPHASE_H_
