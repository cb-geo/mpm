#ifndef MPM_ITERATIVE_EIGEN_H_
#define MPM_ITERATIVE_EIGEN_H_

#include <cmath>

#include "factory.h"
#include "solver_base.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

namespace mpm {

//! MPM Iterative Eigen solver class
//! \brief Iterative linear sparse matrix solver class using Eigen library
template <typename Traits>
class IterativeEigen : public SolverBase<Traits> {
 public:
  //! Constructor
  //! \param[in] max_iter Maximum number of iterations
  //! \param[in] tolerance Tolerance for solver to achieve convergence
  IterativeEigen(unsigned max_iter, double tolerance)
      : mpm::SolverBase<Traits>(max_iter, tolerance) {
    //! Logger
    std::string logger = "EigenIterativeSolver::";
    console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
  };

  //! Destructor
  ~IterativeEigen(){};

  //! Matrix solver with default initial guess
  Eigen::VectorXd solve(const Eigen::SparseMatrix<double>& A,
                        const Eigen::VectorXd& b) override;

  //! Return the type of solver
  std::string solver_type() const { return "Eigen"; }

  //! Assign global active dof
  void assign_global_active_dof(unsigned global_active_dof) override {}

  //! Assign rank to global mapper
  void assign_rank_global_mapper(std::vector<int> rank_global_mapper) override {
  }

 protected:
  //! Solver type
  using SolverBase<Traits>::sub_solver_type_;
  //! Preconditioner type
  using SolverBase<Traits>::preconditioner_type_;
  //! Maximum number of iterations
  using SolverBase<Traits>::max_iter_;
  //! Tolerance
  using SolverBase<Traits>::tolerance_;
  //! Verbosity
  using SolverBase<Traits>::verbosity_;
  //! Logger
  using SolverBase<Traits>::console_;
};
}  // namespace mpm

#include "iterative_eigen.tcc"

#endif  // MPM_ITERATIVE_EIGEN_H_
