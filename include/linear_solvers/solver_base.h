#ifndef MPM_SOLVER_BASE_H_
#define MPM_SOLVER_BASE_H_

#include "data_types.h"
#include "logger.h"

namespace mpm {
template <typename Traits>
class SolverBase {
 public:
  //! Constructor with min and max iterations and tolerance
  //! \param[in] max_iter Maximum number of iterations
  //! \param[in] tolerance Tolerance for solver to achieve convergence
  SolverBase(unsigned max_iter, double tolerance) {
    max_iter_ = max_iter;
    tolerance_ = tolerance;
    //! Logger
    std::string logger = "SolverBase::";
    console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
  };

  //! Destructor
  virtual ~SolverBase(){};

  //! Matrix solver with default initial guess
  virtual Eigen::VectorXd solve(const Eigen::SparseMatrix<double>& A,
                                const Eigen::VectorXd& b,
                                std::string solver_type) = 0;

  //! Assign global active dof
  virtual void assign_global_active_dof(unsigned global_active_dof) = 0;

  //! Assign rank to global mapper
  virtual void assign_rank_global_mapper(
      std::vector<int> rank_global_mapper) = 0;

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