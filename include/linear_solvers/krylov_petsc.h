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
      throw std::runtime_error("Error when initialiazing Petsc.");
    }
    abs_tolerance_ = PETSC_DEFAULT;
    div_tolerance_ = PETSC_DEFAULT;
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
                        const Eigen::VectorXd& b) override;

  //! Return the type of solver
  std::string solver_type() const { return "PETSC"; }

  //! Assign global active dof
  void assign_global_active_dof(unsigned global_active_dof) override {
    global_active_dof_ = global_active_dof;
  };

  //! Assign rank to global mapper
  void assign_rank_global_mapper(std::vector<int> rank_global_mapper) override {
    rank_global_mapper_ = rank_global_mapper;
  };

 protected:
  //! Solver type
  using SolverBase<Traits>::sub_solver_type_;
  //! Preconditioner type
  using SolverBase<Traits>::preconditioner_type_;
  //! Maximum number of iterations
  using SolverBase<Traits>::max_iter_;
  //! Relative tolerance
  using SolverBase<Traits>::tolerance_;
  //! Absolute tolerance
  using SolverBase<Traits>::abs_tolerance_;
  //! Divergence tolerance
  using SolverBase<Traits>::div_tolerance_;
  //! Verbosity
  using SolverBase<Traits>::verbosity_;
  //! Logger
  using SolverBase<Traits>::console_;
  //! Global active dof
  unsigned global_active_dof_;
  //! Rank global Mapper
  std::vector<int> rank_global_mapper_;
};
}  // namespace mpm

#include "krylov_petsc.tcc"
#endif  // MPM_KRYLOV_PETSC_H_
