#ifndef MPM_CG_EIGEN_H_
#define MPM_CG_EIGEN_H_

#include <cmath>

#include "factory.h"
#include "solver_base.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <iostream>

namespace mpm {

//! MPM Eigen CG class
//! \brief Conjugate Gradient solver class using Eigen
template <unsigned Tdim>
class CGEigen : public SolverBase<Tdim> {
 public:
  //! Constructor
  //! \param[in] max_iter Maximum number of iterations

  //! \param[in] tolerance Tolerance for solver to achieve convergence
  CGEigen(unsigned max_iter, double tolerance)
      : mpm::SolverBase<Tdim>(max_iter, tolerance) {
    //! Logger
    console_ = spdlog::stdout_color_mt("EigenSolver");
  };

  //! Matrix solver with default initial guess
  Eigen::VectorXd solve(const Eigen::SparseMatrix<double>& A,
                        const Eigen::VectorXd& b,
                        std::string solver_type) override;

  //! Matrix solver with defined initial guess
  Eigen::VectorXd solve(const Eigen::SparseMatrix<double>& A,
                        const Eigen::VectorXd& b, std::string solver_type,
                        const Eigen::VectorXd& initial_guess) override;

  //! Return the type of solver
  std::string solver_type() const { return "Eigen"; }

 protected:
  //! Maximum number of iterations
  using SolverBase<Tdim>::max_iter_;
  //! Tolerance
  using SolverBase<Tdim>::tolerance_;
  //! Logger
  using SolverBase<Tdim>::console_;
  //! cg_type_ (leastSquaresConjugateGradient or ConjugateGradient)
  std::string cg_type_;
};
}  // namespace mpm

#include "cg_eigen.tcc"

#endif  // MPM_CG_EIGEN_H_
