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
  void assign_mesh_pointer(std::shared_ptr<mpm::Mesh<Tdim>>& mesh) {
    mesh_ = mesh;
  }

  //! Create a pair between nodes and index in Matrix / Vector
  virtual bool assign_global_node_indices(unsigned active_dof) = 0;

  //! Assemble laplacian matrix
  virtual Eigen::SparseMatrix<double>& laplacian_matrix() = 0;

  //! Assemble laplacian matrix
  virtual bool assemble_laplacian_matrix(double dt) = 0;

  //! Assemble poisson RHS vector
  virtual Eigen::VectorXd& poisson_rhs_vector() = 0;

  //! Assemble poisson RHS vector
  virtual bool assemble_poisson_right(std::shared_ptr<mpm::Mesh<Tdim>>& mesh_,
                                      double dt) = 0;

  //! Assign free surface node id
  virtual void assign_free_surface(
      const std::set<mpm::Index>& free_surface_id) = 0;

  //! Assign pressure constraints
  virtual bool assign_pressure_constraints(double beta,
                                           const double current_time) = 0;

  //! Apply pressure constraints to poisson equation
  virtual void apply_pressure_constraints() = 0;

  //! Return pressure increment
  virtual Eigen::VectorXd& pressure_increment() = 0;

  //! Assign pressure increment
  virtual void assign_pressure_increment(
      Eigen::VectorXd pressure_increment) = 0;

  //! Return correction matrix
  virtual Eigen::SparseMatrix<double>& correction_matrix() = 0;

  //! Assemble corrector RHS
  virtual bool assemble_corrector_right(std::shared_ptr<mpm::Mesh<Tdim>>& mesh_,
                                        double dt) = 0;

  //! Assemble displacement vector
  // virtual void assemble_disp_vector() = 0;

  //! Apply displacement to nodes
  // virtual void apply_displacement_nodes() = 0;

  //! Apply forces to nodes
  // virtual void apply_forces_nodes() = 0;

  //! Apply restraints
  // virtual Eigen::VectorXd global_restraints() = 0;

  //! Initialise force vector to zero
  // virtual void initialise_force_zero() = 0;

  virtual std::set<mpm::Index> free_surface() = 0;

  virtual unsigned active_dof() { return active_dof_; };

  //! Assemble stiffness matrix (semi-implicit)
  virtual bool assemble_stiffness_matrix(unsigned dir, double dt) {
    throw std::runtime_error(
        "Calling the base class function "
        "(assemble_stiffness_matrix) in "
        "AssemblerBase:: illegal operation!");
    return false;
  };

  //! Assemble force vector (semi-implicit)
  virtual bool assemble_force_vector(double dt) {
    throw std::runtime_error(
        "Calling the base class function "
        "(assemble_force_vector) in "
        "AssemblerBase:: illegal operation!");
    return false;
  };

  //! Return stiffness matrix
  virtual Eigen::SparseMatrix<double>& stiffness_matrix(unsigned dir) {
    static Eigen::SparseMatrix<double> error;
    throw std::runtime_error(
        "Calling the base class function "
        "(stiffness_matrix) in "
        "AssemblerBase:: illegal operation!");
    return error;
  };

  //! Return intermediate force
  virtual Eigen::MatrixXd& force_inter() {
    static Eigen::MatrixXd error;
    throw std::runtime_error(
        "Calling the base class function "
        "(force_inter) in "
        "AssemblerBase:: illegal operation!");
    return error;
  };

  //! Assign intermediate acceleration
  virtual void assign_intermediate_acceleration(
      unsigned dim, Eigen::VectorXd acceleration_inter) {
    throw std::runtime_error(
        "Calling the base class function "
        "(assign_intermediate_acceleration) in "
        "AssemblerBase:: illegal operation!");
  };

  //! Return intermediate velocity
  virtual Eigen::MatrixXd& acceleration_inter() {
    static Eigen::MatrixXd error;
    throw std::runtime_error(
        "Calling the base class function "
        "(acceleration_inter) in "
        "AssemblerBase:: illegal operation!");
    return error;
  };

  virtual bool assign_velocity_constraints() {
    throw std::runtime_error(
        "Calling the base class function "
        "(assign_velocity_constraints) in "
        "AssemblerBase:: illegal operation!");
    return false;
  };

  //! Apply velocity constraints
  virtual bool apply_velocity_constraints() {
    throw std::runtime_error(
        "Calling the base class function "
        "(apply_velocity_constraints) in "
        "AssemblerBase:: illegal operation!");
    return false;
  };

 protected:
  //! Number of total active_dof
  unsigned active_dof_;
  //! Mesh object
  std::shared_ptr<mpm::Mesh<Tdim>> mesh_;
  //! Sparse Stiffness Matrix
  std::shared_ptr<Eigen::SparseMatrix<double>> stiffness_matrix_;
  //! Force vector
  std::shared_ptr<Eigen::MatrixXd> force_inter_;
  //! Intermediate velocity vector
  std::shared_ptr<Eigen::VectorXd> velocity_inter_vector_;
  //! Displacement vector
  std::shared_ptr<Eigen::MatrixXd> displacement_vector_;
};
}  // namespace mpm

#endif  // MPM_ASSEMBLER_BASE_H_