#ifndef MPM_SOLVER_BASE_H_
#define MPM_SOLVER_BASE_H_

#include "data_types.h"
#include "logger.h"
#include <chrono>

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
                                const Eigen::VectorXd& b) = 0;

  //! Assign global active dof
  virtual void assign_global_active_dof(unsigned global_active_dof) = 0;

  //! Assign rank to global mapper
  virtual void assign_rank_global_mapper(
      std::vector<int> rank_global_mapper) = 0;

  //! Set sub solver type
  void set_sub_solver_type(const std::string& type) noexcept {
    sub_solver_type_ = type;
  }

  //! Set preconditioner type
  void set_preconditioner_type(const std::string& type) noexcept {
    preconditioner_type_ = type;
  }

  //! Set maximum number of iterations
  void set_max_iteration(unsigned max_iter) noexcept { max_iter_ = max_iter; }

  //! Set relative iteration tolerance
  void set_tolerance(double tol) noexcept { tolerance_ = tol; }

  //! Set absolute iteration tolerance
  void set_abs_tolerance(double tol) noexcept { abs_tolerance_ = tol; }

  //! Set divergence iteration tolerance
  void set_div_tolerance(double tol) noexcept { div_tolerance_ = tol; }

  //! Set verbosity
  void set_verbosity(unsigned v) noexcept { verbosity_ = v; }

 protected:
  //! Solver type
  std::string sub_solver_type_{"cg"};
  //! Preconditioner type
  std::string preconditioner_type_{"none"};
  //! Maximum number of iterations
  unsigned max_iter_;
  //! Relative tolerance
  double tolerance_;
  //! Absolute tolerance
  double abs_tolerance_;
  //! Divergence tolerance
  double div_tolerance_;
  //! Verbosity
  unsigned verbosity_{0};
  //! Logger
  std::shared_ptr<spdlog::logger> console_;
};
}  // namespace mpm

#endif