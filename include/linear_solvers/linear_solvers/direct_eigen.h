#ifndef MPM_DIRECT_EIGEN_H_
#define MPM_DIRECT_EIGEN_H_

#include <cmath>

#include "factory.h"
#include "solver_base.h"
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>

namespace mpm {

//! MPM Direct Eigen solver class
//! \brief Direct linear sparse matrix solver class using Eigen library
template <typename Traits>
class DirectEigen : public SolverBase<Traits> {
 public:
  //! Constructor
  DirectEigen(unsigned max_iter, double tolerance)
      : mpm::SolverBase<Traits>(max_iter, tolerance) {
    //! Logger
    std::string logger = "EigenDirectSolver::";
    console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
    //! Default sub solver type
    sub_solver_type_ = "lu";
  };

  //! Destructor
  ~DirectEigen(){};

  //! Matrix solver
  Eigen::VectorXd solve(const Eigen::SparseMatrix<double>& A,
                        const Eigen::VectorXd& b) override;

  //! Return the type of solver
  std::string solver_type() const { return "EigenDirectSolver"; }

  //! Assign global active dof
  void assign_global_active_dof(unsigned global_active_dof) override {}

  //! Assign rank to global mapper
  void assign_rank_global_mapper(std::vector<int> rank_global_mapper) override {
  }

 protected:
  //! Solver type
  using SolverBase<Traits>::sub_solver_type_;
  //! Verbosity
  using SolverBase<Traits>::verbosity_;
  //! Logger
  using SolverBase<Traits>::console_;
};
}  // namespace mpm

#include "direct_eigen.tcc"

#endif  // MPM_DIRECT_EIGEN_H_