//! Conjugate Gradient with default initial guess
template <typename Traits>
Eigen::VectorXd mpm::KrylovPETSC<Traits>::solve(
    const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b) {
  //! Initialize solution vector x
  Eigen::VectorXd x(b.size());
  try {
#if USE_PETSC

    // Initialise MPI mpi_rank and size
    int mpi_rank = 0;
    int mpi_size = 1;

    // Get MPI mpi_rank
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    // Get number of MPI mpi_ranks
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    // Solver start
    auto solver_begin = std::chrono::steady_clock::now();
    if (verbosity_ > 0 && mpi_rank == 0)
      console_->info("Type: \"{}\", Preconditioner: \"{}\", Begin!",
                     sub_solver_type_, preconditioner_type_);

    // Initialize PETSC matrix, vectors, and parameters
    KSP solver;
    Mat petsc_A;
    Vec petsc_b, petsc_x;
    KSPConvergedReason reason;

    // Initialize vector b across the ranks
    VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE, global_active_dof_, &petsc_b);
    VecDuplicate(petsc_b, &petsc_x);

    // Initialize Matrix A across the ranks
    MatCreateAIJ(MPI_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, global_active_dof_,
                 global_active_dof_, 0, NULL, 0, NULL, &petsc_A);
    MatSetOption(petsc_A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

    // Copying Eigen vector b to petsc b vector
    VecSetValues(petsc_b, rank_global_mapper_.size(),
                 rank_global_mapper_.data(), b.data(), ADD_VALUES);
    VecAssemblyBegin(petsc_b);
    VecAssemblyEnd(petsc_b);

    // Output assembled right-hand-side vector in each rank
    PetscViewer viewer;
    if (verbosity_ == 3) {
      PetscViewerASCIIGetStdout(MPI_COMM_WORLD, &viewer);
      VecView(petsc_b, viewer);
    }

    // Copying Eigen matrix A to petsc A matrix
    for (int k = 0; k < A.outerSize(); ++k) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
        MatSetValue(petsc_A, rank_global_mapper_[it.row()],
                    rank_global_mapper_[k], it.value(), ADD_VALUES);
      }
    }
    MatAssemblyBegin(petsc_A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(petsc_A, MAT_FINAL_ASSEMBLY);

    // Output coefficient matrix A
    if (verbosity_ == 3) {
      PetscViewerASCIIGetStdout(MPI_COMM_WORLD, &viewer);
      MatView(petsc_A, viewer);
    }

    // Initiate PETSC solver
    KSPCreate(MPI_COMM_WORLD, &solver);
    KSPSetOperators(solver, petsc_A, petsc_A);

    // Set solver_type
    if (sub_solver_type_ == "cg") {
      KSPSetType(solver, KSPCG);
      KSPCGSetType(solver, KSP_CG_SYMMETRIC);
    } else if (sub_solver_type_ == "gmres")
      KSPSetType(solver, KSPGMRES);
    else if (sub_solver_type_ == "bicgstab")
      KSPSetType(solver, KSPBCGS);
    else if (sub_solver_type_ == "lsqr")
      KSPSetType(solver, KSPLSQR);
    else {
      if (verbosity_ > 0 && mpi_rank == 0)
        console_->warn(
            "Sub solver type is not available! Using \"gmres\" as default "
            "type. Available sub solver type implemented in KrylovPETSC "
            "class are: \"cg\", \"gmres\", \"lsqr\", and "
            "\"bicgstab\".");
    }

    // Set tolerance
    KSPSetTolerances(solver, tolerance_, abs_tolerance_, div_tolerance_,
                     max_iter_);

    // Set preconditioner
    if (preconditioner_type_ != "none") {
      PC pc;
      KSPSetInitialGuessNonzero(solver, PETSC_TRUE);
      KSPGetPC(solver, &pc);
      PCFactorSetShiftType(pc, MAT_SHIFT_POSITIVE_DEFINITE);
      if (preconditioner_type_ == "jacobi") PCSetType(pc, PCJACOBI);
      if (preconditioner_type_ == "bjacobi") PCSetType(pc, PCBJACOBI);
      if (preconditioner_type_ == "pbjacobi") PCSetType(pc, PCPBJACOBI);
      if (preconditioner_type_ == "asm") PCSetType(pc, PCASM);
      if (preconditioner_type_ == "eisenstat") PCSetType(pc, PCEISENSTAT);
      if (preconditioner_type_ == "icc") PCSetType(pc, PCICC);
    }

    // Solve linear system of equation x = A^(-1) b
    KSPSolve(solver, petsc_b, petsc_x);
    KSPGetConvergedReason(solver, &reason);

    // Print residual in each iteration
    if (verbosity_ >= 1) {
      PetscViewerAndFormat* vf;
      PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,
                                 PETSC_VIEWER_DEFAULT, &vf);
      PetscInt its;
      KSPGetIterationNumber(solver, &its);
      PetscPrintf(PETSC_COMM_WORLD, "\nConvergence in %d iterations.\n",
                  (int)its);
      PetscReal rnorm;
      if (verbosity_ >= 2) {
        for (int i = 0; i < its; i++) {
          KSPMonitorTrueResidualNorm(solver, i, rnorm, vf);
        }
      }
    }

    // Warn if solver does not converge
    if (reason < 0) {
      PetscPrintf(MPI_COMM_WORLD,
                  "\nKrylov PETSC solver \"%s\" with \"%s\" preconditioner "
                  "DIVERGED, try to modify the preconditioner, set tolerance "
                  "and maximum iteration.\n",
                  sub_solver_type_.c_str(), preconditioner_type_.c_str());
    }

    // Scatter and gather for cloning procedure
    if (mpi_size > 1) {
      // Initiate scatter arrays
      VecScatter ctx;
      Vec x_seq;
      PetscScalar* x_data;
      VecScatterCreateToAll(petsc_x, &ctx, &x_seq);
      VecScatterBegin(ctx, petsc_x, x_seq, INSERT_VALUES, SCATTER_FORWARD);
      VecScatterEnd(ctx, petsc_x, x_seq, INSERT_VALUES, SCATTER_FORWARD);
      VecGetArray(x_seq, &x_data);

      // Copy petsc x to Eigen x
      for (unsigned i = 0; i < x.size(); i++) {
        const int global_index = rank_global_mapper_[i];
        x(i) = x_data[global_index];
      }

      // Destroy scatter arrays
      VecRestoreArray(x_seq, &x_data);
      VecScatterDestroy(&ctx);
      VecDestroy(&x_seq);
    } else {
      PetscScalar value;
      // Copy petsc x to Eigen x
      for (unsigned i = 0; i < x.size(); i++) {
        const int global_index = rank_global_mapper_[i];
        VecGetValues(petsc_x, 1, &global_index, &value);
        x(i) = value;
      }
    }

    // Free work space
    VecDestroy(&petsc_x);
    VecDestroy(&petsc_b);
    MatDestroy(&petsc_A);
    KSPDestroy(&solver);

    // Solver End
    auto solver_end = std::chrono::steady_clock::now();
    if (verbosity_ > 0 && mpi_rank == 0)
      console_->info(
          "Type: \"{}\", Preconditioner: \"{}\", End! Duration: {} ms.",
          sub_solver_type_, preconditioner_type_,
          std::chrono::duration_cast<std::chrono::milliseconds>(solver_end -
                                                                solver_begin)
              .count());

#endif
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }

  return x;
}