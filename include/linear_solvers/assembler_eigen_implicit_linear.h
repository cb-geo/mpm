#ifndef MPM_ASSEMBLER_EIGEN_IMPLICIT_LINEAR_H_
#define MPM_ASSEMBLER_EIGEN_IMPLICIT_LINEAR_H_

#include <Eigen/Sparse>
#include <string>

// Speed log
#include "assembler_base.h"
#include "spdlog/spdlog.h"

#include "mesh.h"

namespace mpm {
template <unsigned Tdim>
class AssemblerEigenImplicitLinear : public AssemblerBase<Tdim> {
 public:
  //! Constructor
  //! \param[in] node_neighbourhood Number of node neighbourhood considered
  AssemblerEigenImplicitLinear(unsigned node_neighbourhood);

  /**
   * \defgroup Implicit Functions dealing with implicit MPM
   */
  /**@{*/
  //! Return stiffness matrix
  //! \ingroup Implicit
  Eigen::SparseMatrix<double>& stiffness_matrix() override {
    return stiffness_matrix_;
  }

  //! Assemble stiffness matrix
  //! \ingroup Implicit
  bool assemble_stiffness_matrix() override;

  //! Return residual force RHS vector
  //! \ingroup Implicit
  Eigen::VectorXd& residual_force_rhs_vector() override {
    return residual_force_rhs_vector_;
  }

  //! Assemble residual force RHS vector
  //! \ingroup Implicit
  bool assemble_residual_force_right() override;

  //! Assign displacement constraints
  //! \ingroup Implicit
  bool assign_displacement_constraints(double current_time) override;

  //! Apply displacement constraints to equilibrium equation
  //! \ingroup Implicit
  void apply_displacement_constraints() override;

  //! Return displacement increment
  //! \ingroup Implicit
  Eigen::VectorXd& displacement_increment() override {
    return displacement_increment_;
  }

  //! Assign displacement increment
  //! \ingroup Implicit
  void assign_displacement_increment(
      const Eigen::VectorXd& displacement_increment) override {
    displacement_increment_ = displacement_increment;
  }
  /**@{*/

 protected:
  //! number of nodes
  using AssemblerBase<Tdim>::active_dof_;
  //! Mesh object
  using AssemblerBase<Tdim>::mesh_;
  //! Number of sparse matrix container size
  using AssemblerBase<Tdim>::sparse_row_size_;
  //! Global node indices
  using AssemblerBase<Tdim>::global_node_indices_;
  //! Number of total active_dof in all rank
  using AssemblerBase<Tdim>::global_active_dof_;
  //! Rank to Global mapper
  using AssemblerBase<Tdim>::rank_global_mapper_;
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
  //! Stiffness matrix
  Eigen::SparseMatrix<double> stiffness_matrix_;
  //! Residual force RHS vector
  Eigen::VectorXd residual_force_rhs_vector_;
  //! Displacement constraints
  Eigen::SparseVector<double> displacement_constraints_;
  //! Displacement increment
  Eigen::VectorXd displacement_increment_;
  /**@{*/
};
}  // namespace mpm

#include "assembler_eigen_implicit_linear.tcc"
#endif  // MPM_ASSEMBLER_EIGEN_IMPLICIT_LINEAR_H_
