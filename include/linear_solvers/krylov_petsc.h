#ifndef MPM_KRYLOV_PETSC_H_
#define MPM_KRYLOV_PETSC_H_

#include <cmath>

#include "factory.h"
#include "solver_base.h"
#include <iostream>
#include <petscksp.h>

namespace mpm {

//! MPM Eigen CG class
//! \brief Conjugate Gradient solver class using Eigen
template <typename Traits>
class KrylovPETSC : public SolverBase<Traits> {
 public:
  //! Constructor
  //! \param[in] max_iter Maximum number of iterations

  //! \param[in] tolerance Tolerance for solver to achieve convergence
  KrylovPETSC(unsigned max_iter, double tolerance)
      : mpm::SolverBase<Traits>(max_iter, tolerance) {
    //! Logger
    console_ = spdlog::stdout_color_mt("PETSCSolver");
  };

  //! Matrix solver with default initial guess
  Eigen::VectorXd solve(const Eigen::SparseMatrix<double>& A,
                        const Eigen::VectorXd& b,
                        std::string solver_type) override;

  //! Return the type of solver
  std::string solver_type() const { return "PETSC"; }

 protected:
  //! Maximum number of iterations
  using SolverBase<Traits>::max_iter_;
  //! Tolerance
  using SolverBase<Traits>::tolerance_;
  //! Logger
  using SolverBase<Traits>::console_;
  //! cg_type_ (leastSquaresConjugateGradient or ConjugateGradient)
  std::string cg_type_;
};
}  // namespace mpm

#include "krylov_petsc.tcc"

#endif  // MPM_KRYLOV_PETSC_H_
