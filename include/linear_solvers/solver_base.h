#ifndef MPM_SOLVER_BASE_H_
#define MPM_SOLVER_BASE_H_

#include "logger.h"
namespace mpm {
template <typename Traits>
class SolverBase {
 public:
  //! Constructor with min and max iterations and tolerance
  SolverBase(unsigned max_iter, double tolerance) {
    max_iter_ = max_iter;
    tolerance_ = tolerance;
    //! Logger
    console_ = std::make_unique<spdlog::logger>("SolverBase", mpm::stdout_sink);
  };

  //! Matrix solver with default initial guess
  virtual Eigen::VectorXd solve(const Eigen::SparseMatrix<double>& A,
                                const Eigen::VectorXd& b,
                                std::string solver_type) = 0;

 protected:
  //! Maximum number of iterations
  unsigned max_iter_;
  //! Tolerance
  double tolerance_;
  //! Logger
  std::shared_ptr<spdlog::logger> console_;
};
}  // namespace mpm

#endif