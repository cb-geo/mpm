//! Conjugate Gradient with default initial guess
template <typename Traits>
Eigen::VectorXd mpm::KrylovPETSC<Traits>::solve(
    const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b,
    std::string solver_type) {
  Eigen::VectorXd x;

  try {

    // Task
    // 1. Convert eigen to PETSC type
    // 2. Call PETSC solver
    // 3. Convert back from PETSC type to Eigen for output
    KSP solver;
    Mat petsc_A;
    Vec petsc_b, petsc_x;
    PetscErrorCode ierr;
    MPI_Comm comm;
    PetscInt dim = b.size();
    PetscInt rank, vi, mi, mj;
    PetscScalar v, m;
    KSPConvergedReason reason;

    PetscInt low, high, rlow, rhigh;

    int petsc_argc = 1;
    char* petsc_arg = "petsc_argv";
    char** petsc_argv = &petsc_arg;
    PetscInitialize(&petsc_argc, &petsc_argv, 0, 0);

    comm = PETSC_COMM_WORLD;
    MPI_Bcast(&dim, 1, MPI_INT, 0, comm);

    VecCreate(comm, &petsc_b);
    VecSetSizes(petsc_b, PETSC_DECIDE, dim);
    VecSetType(petsc_b, VECMPI);
    VecDuplicate(petsc_b, &petsc_x);

    MatCreate(PETSC_COMM_WORLD, &petsc_A);
    MatSetSizes(petsc_A, PETSC_DECIDE, PETSC_DECIDE, dim, dim);
    MatSetType(petsc_A, MATSEQAIJ);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank == 0) {
      for (vi = 0; vi < dim; vi++) {
        v = b(vi);
        VecSetValues(petsc_b, 1, &vi, &v, INSERT_VALUES);
      }
      // This for loop should be optimized for sparse A
      for (mi = 0; mi < dim; mi++) {
        for (mj = 0; mj < dim; mj++) {
          m = A.coeff(mi, mj);
          MatSetValue(petsc_A, mi, mj, m, INSERT_VALUES);
        }
      }
    }

    // If b and A are broadcasted, parallel seting value
    /*VecGetOwnershipRange(petsc_b, &low, &high);
    for (vi = low; vi < high; vi++){
      v = b(vi);
      VecSetValues(petsc_b, 1, &vi, &v, INSERT_VALUES);
    }
    MatGetOwnershipRange(petsc_A, &rlow, &rhigh);
    for (mi = rlow; mi < rhigh; mi++){
      for(mj = 0; mj < dim; mj++){
            m = A.coeffRef(mi, mj);
            MatSetValue(petsc_A, 1, &mi, 1, &mj, &m, INSERT_VALUES);
          }
    }*/

    MatAssemblyBegin(petsc_A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(petsc_A, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(petsc_b);
    VecAssemblyEnd(petsc_b);

    if (solver_type == "cg") {
      //   Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;

      //   solver.setMaxIterations(max_iter_);
      //   solver.setTolerance(tolerance_);
      //   solver.compute(A);

      //   x = solver.solve(b);

      //   if (solver.info() != Eigen::Success) {
      //     throw std::runtime_error("Fail to solve linear systems!\n");
      //   }
      KSPCreate(comm, &solver);
      KSPSetOperators(solver, petsc_A, petsc_A);
      KSPSetType(solver, KSPCG);
      KSPSolve(solver, petsc_b, petsc_x);
      KSPGetConvergedReason(solver, &reason);
      if (reason < 0) {
        PetscPrintf(PETSC_COMM_WORLD, "\nKSPCG solver Diverged;\n");
      }

      if (rank == 0) {
        for (vi = 0; vi < dim; vi++) {
          VecGetValues(petsc_x, 1, &vi, &v);
          x(vi) = v;
        }
      }

    } else if (solver_type == "lscg") {

      //   // Another option is LDLT, but not accurate as lscg
      //   // Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
      //   Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>>
      //   solver;

      //   solver.setMaxIterations(max_iter_);
      //   solver.setTolerance(tolerance_);
      //   solver.compute(A);

      //   x = solver.solve(b);

      //   if (solver.info() != Eigen::Success) {
      //     throw std::runtime_error("Fail to solve linear systems!\n");
      //   }
    }

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return x;
}