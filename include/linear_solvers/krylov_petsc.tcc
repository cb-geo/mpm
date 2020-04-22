//! Conjugate Gradient with default initial guess
template <typename Traits>
Eigen::VectorXd mpm::KrylovPETSC<Traits>::solve(
    const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b,
    std::string solver_type) {
  Eigen::VectorXd x;

  try {

#if USE_PETSC


    std::cout << "TEST OUTPUT: 0: " << std::endl;

    // Initialize PETSC parameters
    KSP solver;
    Mat petsc_A;
    Vec petsc_b, petsc_x;
    PetscErrorCode ierr;
    MPI_Comm comm;
    PetscInt dim = b.size();
    PetscInt rank, vi, mi, mj;
    PetscScalar v, m;
    KSPConvergedReason reason;
    PetscViewer viewer;
    PetscInt low, high, rlow, rhigh;

    std::cout << "TEST OUTPUT: 1: " << std::endl;

    int petsc_argc = 1;
    char* petsc_arg = "petsc_argv";
    char** petsc_argv = &petsc_arg;
    PetscInitialize(&petsc_argc, &petsc_argv, 0, 0);

    std::cout << "TEST OUTPUT: 2: " << std::endl;

    comm = PETSC_COMM_WORLD;
    // Broadcast DIM to all ranks
    MPI_Bcast(&dim, 1, MPI_INT, 0, comm);

    std::cout << "TEST OUTPUT: 3: " << std::endl;

    // TODO: Now do all rank create matrix and vector? Or only rank 0 should run
    // these lines Create PETSC_B and PETSC_X with global dim = DIM,
    // PETSC_DECIDE local dim
    VecCreate(comm, &petsc_b);
    VecSetSizes(petsc_b, PETSC_DECIDE, dim);
    VecSetType(petsc_b, VECMPI);
    VecDuplicate(petsc_b, &petsc_x);

    std::cout << "TEST OUTPUT: 4: " << std::endl;

    // Create PETSC_A with global dim = DIM by DIM, PETSC_DECIDE local dim
    MatCreate(PETSC_COMM_WORLD, &petsc_A);
    MatSetSizes(petsc_A, PETSC_DECIDE, PETSC_DECIDE, dim, dim);
    MatSetType(petsc_A, MATMPIAIJ);

    std::cout << "TEST OUTPUT: 5: " << std::endl;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    // TODO: Optimize send from value-to-value to array-to-array
    if (rank == 0) {
      // Copy b to petsc_b
      for (vi = 0; vi < dim; vi++) {
        v = b(vi);
        VecSetValues(petsc_b, 1, &vi, &v, INSERT_VALUES);
      }
      // Copy A to petsc_A. Can be optimized for sparse A
      for (mi = 0; mi < dim; mi++) {
        for (mj = 0; mj < dim; mj++) {
          m = A.coeff(mi, mj);
          MatSetValue(petsc_A, mi, mj, m, INSERT_VALUES);
        }
      }
    }

    std::cout << "TEST OUTPUT: 6: " << rank << std::endl;

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

    MPI_Barrier(comm);

    // Assemble PETSC_A, PETSC_b
    MatAssemblyBegin(petsc_A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(petsc_A, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(petsc_b);
    VecAssemblyEnd(petsc_b);

    std::cout << "TEST OUTPUT: 7: " << std::endl;

    if (solver_type == "cg") {
      KSPCreate(comm, &solver);
      KSPSetOperators(solver, petsc_A, petsc_A);
      KSPSetType(solver, KSPCG);
      KSPSolve(solver, petsc_b, petsc_x);
      KSPGetConvergedReason(solver, &reason);
      if (reason < 0) {
        PetscPrintf(PETSC_COMM_WORLD, "\nKSPCG solver Diverged;\n");
      }

      // Copy petsc_x to x
      if (rank == 0) {
        for (vi = 0; vi < dim; vi++) {
          VecGetValues(petsc_x, 1, &vi, &v);
          x(vi) = v;
        }
      }
      PetscFinalize();
    }

    std::cout << "TEST OUTPUT: 8: " << std::endl;

#endif

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return x;
}