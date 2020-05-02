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
  virtual bool assemble_corrector_right(double dt) = 0;

  virtual unsigned global_active_dof(){};

  virtual std::vector<int> rank_global_mapper(){};

 protected:
  //! Number of total active_dof
  unsigned active_dof_;
  //! Mesh object
  std::shared_ptr<mpm::Mesh<Tdim>> mesh_;
};
}  // namespace mpm

#endif  // MPM_ASSEMBLER_BASE_H_