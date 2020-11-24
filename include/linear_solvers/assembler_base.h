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
  AssemblerBase() {
    //! Global degrees of freedom
    active_dof_ = 0;
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

  //! Assemble laplacian matrix
  virtual Eigen::SparseMatrix<double>& laplacian_matrix() = 0;

  //! Assemble laplacian matrix
  virtual bool assemble_laplacian_matrix(double dt) = 0;

  //! Assemble poisson RHS vector
  virtual Eigen::VectorXd& poisson_rhs_vector() = 0;

  //! Assemble poisson RHS vector
  virtual bool assemble_poisson_right(double dt) = 0;

  //! Assign free surface node id
  virtual void assign_free_surface(
      const std::set<mpm::Index>& free_surface_id) = 0;

  //! Assign pressure constraints
  virtual bool assign_pressure_constraints(double beta,
                                           double current_time) = 0;

  //! Apply pressure constraints to poisson equation
  virtual void apply_pressure_constraints() = 0;

  //! Return pressure increment
  virtual Eigen::VectorXd& pressure_increment() = 0;

  //! Assign pressure increment
  virtual void assign_pressure_increment(
      const Eigen::VectorXd& pressure_increment) = 0;

  //! Return correction matrix
  virtual Eigen::SparseMatrix<double>& correction_matrix() = 0;

  //! Assemble corrector RHS
  virtual bool assemble_corrector_right(double dt) = 0;

  //! Return the total size of global dof in all rank
  virtual unsigned global_active_dof() = 0;

  //! Return number of total active_dof
  virtual unsigned active_dof() { return active_dof_; };

  //! Return a vector to map local (rank) index to global index
  virtual std::vector<int> rank_global_mapper() = 0;

  //! TwoPhase functions--------------------------------------------------------
  //! Assemble coefficient matrix for two-phase predictor
  virtual bool assemble_predictor_left(unsigned dir, double dt) {
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
};
}  // namespace mpm

#endif  // MPM_ASSEMBLER_BASE_H_