#ifndef MPM_ASSEMBLER_BASE_H_
#define MPM_ASSEMBLER_BASE_H_

#include <algorithm>
#include <array>
#include <memory>
#include <utility>
#include <vector>

#include <Eigen/Sparse>

#include "mesh.h"
#include "node_base.h"

namespace mpm {

// Linear solver assembler base class
//! \brief Perform matrix assembly procedures
//! \details Build global matrices considering local element matrices
//! \tparam Tdim Dimension
template <unsigned Tdim>
class AssemblerBase {
 public:
  //! Constructor
  //! \param[in] node_neighbourhood Number of node neighbourhood considered
  AssemblerBase(unsigned node_neighbourhood) {
    //! Global degrees of freedom
    active_dof_ = 0;
    //! Assign sparse row size
    switch (node_neighbourhood) {
      case 0:
        sparse_row_size_ = (Tdim == 2) ? 9 : 27;
        break;
      case 1:
        sparse_row_size_ = (Tdim == 2) ? 25 : 125;
        break;
      case 2:
        sparse_row_size_ = (Tdim == 2) ? 49 : 343;
        break;
      default:
        sparse_row_size_ = (Tdim == 2) ? 9 : 27;
    }
  }

  // Virtual destructor
  virtual ~AssemblerBase() = default;

  //! Copy constructor
  AssemblerBase(const AssemblerBase<Tdim>&) = default;

  //! Assignment operator
  AssemblerBase& operator=(const AssemblerBase<Tdim>&) = default;

  //! Move constructor
  AssemblerBase(AssemblerBase<Tdim>&&) = default;

  //! Assign mesh pointer
  //! \param[in] mesh mesh pointer
  void assign_mesh_pointer(const std::shared_ptr<mpm::Mesh<Tdim>>& mesh) {
    mesh_ = mesh;
  }

  //! Create a pair between nodes and index in Matrix / Vector
  virtual bool assign_global_node_indices(unsigned nactive_node,
                                          unsigned nglobal_active_node) = 0;

  /**
   * \defgroup Implicit Functions dealing with implicit MPM
   */
  /**@{*/
  //! Return stiffness matrix
  //! \ingroup Implicit
  virtual Eigen::SparseMatrix<double>& stiffness_matrix() {
    throw std::runtime_error(
        "Calling the base class function (stiffness_matrix) in "
        "AssemblerBase:: illegal operation!");
  };

  //! Assemble stiffness matrix
  //! \ingroup Implicit
  virtual bool assemble_stiffness_matrix() {
    throw std::runtime_error(
        "Calling the base class function (assemble_stiffness_matrix) in "
        "AssemblerBase:: illegal operation!");
    return 0;
  };

  //! Return residual force RHS vector
  //! \ingroup Implicit
  virtual Eigen::VectorXd& residual_force_rhs_vector() {
    throw std::runtime_error(
        "Calling the base class function (residual_force_rhs_vector) in "
        "AssemblerBase:: illegal operation!");
  };

  //! Assemble residual force RHS vector
  //! \ingroup Implicit
  virtual bool assemble_residual_force_right() {
    throw std::runtime_error(
        "Calling the base class function (assemble_residual_force_right) in "
        "AssemblerBase:: illegal operation!");
    return 0;
  };

  //! Assign displacement constraints
  //! \ingroup Implicit
  virtual bool assign_displacement_constraints(double current_time) {
    throw std::runtime_error(
        "Calling the base class function (assign_displacement_constraints) in "
        "AssemblerBase:: illegal operation!");
    return 0;
  };

  //! Apply displacement constraints to equilibrium equation
  //! \ingroup Implicit
  virtual void apply_displacement_constraints() {
    throw std::runtime_error(
        "Calling the base class function (apply_displacement_constraints) in "
        "AssemblerBase:: illegal operation!");
  };

  //! Return displacement increment
  //! \ingroup Implicit
  virtual Eigen::VectorXd& displacement_increment() {
    throw std::runtime_error(
        "Calling the base class function (displacement_increment) in "
        "AssemblerBase:: illegal operation!");
  };

  //! Assign displacement increment
  //! \ingroup Implicit
  virtual void assign_displacement_increment(
      const Eigen::VectorXd& displacement_increment) {
    throw std::runtime_error(
        "Calling the base class function (assign_displacement_increment) in "
        "AssemblerBase:: illegal operation!");
  };
  /**@{*/

  //! Navier-Stokes
  //! functions-------------------------------------------------------- Return
  //! laplacian matrix
  virtual Eigen::SparseMatrix<double>& laplacian_matrix() {
    throw std::runtime_error(
        "Calling the base class function (laplacian_matrix) in "
        "AssemblerBase:: illegal operation!");
  };

  //! Assemble laplacian matrix
  virtual bool assemble_laplacian_matrix(double dt) {
    throw std::runtime_error(
        "Calling the base class function (assemble_laplacian_matrix) in "
        "AssemblerBase:: illegal operation!");
    return 0;
  };

  //! Return poisson RHS vector
  virtual Eigen::VectorXd& poisson_rhs_vector() {
    throw std::runtime_error(
        "Calling the base class function (poisson_rhs_vector) in "
        "AssemblerBase:: illegal operation!");
  };

  //! Assemble poisson RHS vector
  virtual bool assemble_poisson_right(double dt) {
    throw std::runtime_error(
        "Calling the base class function (assemble_poisson_right) in "
        "AssemblerBase:: illegal operation!");
    return 0;
  };

  //! Assign free surface node id
  virtual void assign_free_surface(
      const std::set<mpm::Index>& free_surface_id) {
    throw std::runtime_error(
        "Calling the base class function (assign_free_surface) in "
        "AssemblerBase:: illegal operation!");
  };

  //! Assign pressure constraints
  virtual bool assign_pressure_constraints(double beta, double current_time) {
    throw std::runtime_error(
        "Calling the base class function (assign_pressure_constraints) in "
        "AssemblerBase:: illegal operation!");
    return 0;
  };

  //! Apply pressure constraints to poisson equation
  virtual void apply_pressure_constraints() {
    throw std::runtime_error(
        "Calling the base class function (apply_pressure_constraints) in "
        "AssemblerBase:: illegal operation!");
  };

  //! Return pressure increment
  virtual Eigen::VectorXd& pressure_increment() {
    throw std::runtime_error(
        "Calling the base class function (pressure_increment) in "
        "AssemblerBase:: illegal operation!");
  };

  //! Assign pressure increment
  virtual void assign_pressure_increment(
      const Eigen::VectorXd& pressure_increment) {
    throw std::runtime_error(
        "Calling the base class function (assign_pressure_increment) in "
        "AssemblerBase:: illegal operation!");
  };

  //! Return correction matrix
  virtual Eigen::SparseMatrix<double>& correction_matrix() {
    throw std::runtime_error(
        "Calling the base class function (correction_matrix) in "
        "AssemblerBase:: illegal operation!");
  };

  //! Assemble corrector RHS
  virtual bool assemble_corrector_right(double dt) {
    throw std::runtime_error(
        "Calling the base class function (assemble_corrector_right) in "
        "AssemblerBase:: illegal operation!");
  };

  //! Return the total size of global dof in all rank
  virtual unsigned global_active_dof() = 0;

  //! Return number of total active_dof
  virtual unsigned active_dof() { return active_dof_; };

  //! Return a vector to map local (rank) index to global index
  virtual std::vector<int> rank_global_mapper() = 0;

  //! TwoPhase functions--------------------------------------------------------
  //! Assemble coefficient matrix for two-phase predictor
  virtual bool assemble_predictor_left(double dt) {
    throw std::runtime_error(
        "Calling the base class function (assemble_predictor_left) in "
        "AssemblerBase:: illegal operation!");
    return 0;
  };

  //! Return predictor coefficient LHS matrix
  virtual Eigen::SparseMatrix<double>& predictor_lhs_matrix(unsigned dir) {
    throw std::runtime_error(
        "Calling the base class function (predictor_lhs_matrix) in "
        "AssemblerBase:: illegal operation!");
  };

  //! Assemble RHS force vector for two-phase predictor
  virtual bool assemble_predictor_right(double dt) {
    throw std::runtime_error(
        "Calling the base class function (assemble_predictor_right) in "
        "AssemblerBase:: illegal operation!");
    return 0;
  };

  //! Return predictor RHS force vector
  virtual Eigen::MatrixXd& predictor_rhs_vector() {
    throw std::runtime_error(
        "Calling the base class function (predictor_rhs_vector) in "
        "AssemblerBase:: illegal operation!");
  };

  //! Assign velocity constraints for matrix and vector
  virtual bool assign_velocity_constraints() {
    throw std::runtime_error(
        "Calling the base class function (assign_velocity_constraints) in "
        "AssemblerBase:: illegal operation!");
    return 0;
  };

  //! Apply velocity constraints for matrix and vector
  virtual bool apply_velocity_constraints() {
    throw std::runtime_error(
        "Calling the base class function (apply_velocity_constraints) in "
        "AssemblerBase:: illegal operation!");
    return 0;
  };

  //! Assign nodal intermediate acceleration
  virtual void assign_intermediate_acceleration(
      unsigned dim, Eigen::VectorXd acceleration_inter) {
    throw std::runtime_error(
        "Calling the base class function (assign_intermediate_acceleration) in "
        "AssemblerBase:: illegal operation!");
  };

  //! Return intermediate acceleration
  virtual Eigen::MatrixXd& intermediate_acceleration() {
    throw std::runtime_error(
        "Calling the base class function (intermediate_acceleration) in "
        "AssemblerBase:: illegal operation!");
  };

 protected:
  //! Number of total active_dof
  unsigned active_dof_;
  //! Mesh object
  std::shared_ptr<mpm::Mesh<Tdim>> mesh_;
  //! Number of sparse matrix container size
  unsigned sparse_row_size_;
};
}  // namespace mpm

#endif  // MPM_ASSEMBLER_BASE_H_