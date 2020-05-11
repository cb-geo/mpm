#ifndef MPM_KRYLOV_PETSC_H_
#define MPM_KRYLOV_PETSC_H_

#include <cmath>

#include "factory.h"
#include "mpi.h"
#include "solver_base.h"
#include <iostream>

#ifdef USE_PETSC
#include <petscksp.h>
#endif

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
    std::string logger = "PETSCKrylovSolver::";
    console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);

#ifdef USE_PETSC
    // Initialize PETSc
    int petsc_argc = 0;
    char* petsc_arg0 = NULL;
    char** petsc_argv = &petsc_arg0;
    PetscErrorCode ierr;
    PetscInitialize(&petsc_argc, &petsc_argv, 0, 0);
    if (ierr) {
      throw std::runtime_error("Error when initializaing Petsc.");
    }
#endif
  };

  //! Destructor
  ~KrylovPETSC() {
#ifdef USE_PETSC
    // Finalize PETSc
    PetscFinalize();
#endif
  }

  //! Matrix solver with default initial guess
  Eigen::VectorXd solve(const Eigen::SparseMatrix<double>& A,
                        const Eigen::VectorXd& b,
                        std::string solver_type) override;

  //! Return the type of solver
  std::string solver_type() const { return "PETSC"; }

  void assign_global_active_dof(unsigned global_active_dof) override {
    global_active_dof_ = global_active_dof;
  };

  void assign_rank_global_mapper(std::vector<int> rank_global_mapper) override {
    rank_global_mapper_ = rank_global_mapper;
  };

 protected:
  //! Maximum number of iterations
  using SolverBase<Traits>::max_iter_;
  //! Tolerance
  using SolverBase<Traits>::tolerance_;
  //! Logger
  using SolverBase<Traits>::console_;
  //! cg_type_ (leastSquaresConjugateGradient or ConjugateGradient)
  std::string cg_type_;
  //! Global active dof
  unsigned global_active_dof_;
  //! Rank global Mapper
  std::vector<int> rank_global_mapper_;
};
}  // namespace mpm

#include "krylov_petsc.tcc"
#endif  // MPM_KRYLOV_PETSC_H_
