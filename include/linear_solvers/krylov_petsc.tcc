//! Conjugate Gradient with default initial guess
template <typename Traits>
Eigen::VectorXd mpm::KrylovPETSC<Traits>::solve(
    const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b,
    std::string solver_type) {
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

    // Initialize PETSC matrix, vectors, and parameters
    KSP solver;
    Mat petsc_A;
    Vec petsc_b, petsc_x;
    PetscInt dim = global_active_dof_;
    KSPConvergedReason reason;
    PetscViewer viewer;

    PetscInt vi, mi, mj;
    PetscScalar v, m;
    PetscInt low, high, row_low, row_high;

    //! Initialize Petsc
    int petsc_argc = 1;
    char* petsc_arg = "p";
    char** petsc_argv = &petsc_arg;
    PetscInitialize(&petsc_argc, &petsc_argv, 0, 0);

    //! Initialize vector b across the ranks
    VecCreateMPI(MPI_COMM_WORLD, b.size(), global_active_dof_, &petsc_b);
    VecDuplicate(petsc_b, &petsc_x);

    // Initialize Matrix A across the ranks
    MatCreateAIJ(MPI_COMM_WORLD, A.rows(), A.cols(), global_active_dof_,
                 global_active_dof_, 0, NULL, 0, NULL, &petsc_A);
    MatSetOption(petsc_A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

    //! Copying Eigen matrix b to petsc b vector
    VecSetValues(petsc_b, rank_global_mapper_.size(),
                 rank_global_mapper_.data(), b.data(), INSERT_VALUES);
    VecAssemblyBegin(petsc_b);
    VecAssemblyEnd(petsc_b);

    // PetscViewerASCIIGetStdout(MPI_COMM_WORLD, &viewer);
    // VecView(petsc_b, viewer);
    // MPI_Barrier(MPI_COMM_WORLD);

    // for (unsigned i = 0; i < b.size(); i++) {
    //   auto value = b[i];
    //   unsigned global_index = rank_global_mapper_[i];
    //   VecSetValues(petsc_b, 1, &global_index, &value, INSERT_VALUES);
    // }

    std::cout << "TEST OUTPUT: 5.1: " << std::endl;
    MatGetOwnershipRange(petsc_A, &row_low, &row_high);
    for (mi = row_low; mi < row_high; mi++) {
      for (mj = 0; mj < dim; mj++) {
        m = A.coeff(mi, mj);
        MatSetValue(petsc_A, mi, mj, m, INSERT_VALUES);
      }
    }

    std::cout << "TEST OUTPUT: 6: " << mpi_rank << std::endl;

    // If b and A are broadcasted, parallel seting value
    /*VecGetOwnershipRange(petsc_b, &low, &high);
    for (vi = low; vi < high; vi++){
      v = b(vi);
      VecSetValues(petsc_b, 1, &vi, &v, INSERT_VALUES);
    }
    MatGetOwnershipRange(petsc_A, &row_low, &row_high);
    for (mi = row_low; mi < row_high; mi++){
      for(mj = 0; mj < dim; mj++){
            m = A.coeffRef(mi, mj);
            MatSetValue(petsc_A, 1, &mi, 1, &mj, &m, INSERT_VALUES);
          }
    }*/

    MPI_Barrier(MPI_COMM_WORLD);

    // Assemble PETSC_A, PETSC_b
    MatAssemblyBegin(petsc_A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(petsc_A, MAT_FINAL_ASSEMBLY);
    if (mpi_rank == 0) {
      for (vi = 0; vi < dim; vi++) {
        VecGetValues(petsc_b, 1, &vi, &v);
        std::cout << v << std::endl;
      }
    }
    std::cout << "TEST OUTPUT: 7: " << std::endl;

    if (solver_type == "cg") {
      KSPCreate(MPI_COMM_WORLD, &solver);
      std::cout << "TEST OUTPUT: 7.1: " << std::endl;
      KSPSetOperators(solver, petsc_A, petsc_A);
      std::cout << "TEST OUTPUT: 7.2: " << std::endl;
      KSPSetType(solver, KSPCG);
      std::cout << "TEST OUTPUT: 7.3: " << std::endl;
      KSPSolve(solver, petsc_b, petsc_x);
      std::cout << "TEST OUTPUT: 7.4: " << std::endl;
      KSPGetConvergedReason(solver, &reason);
      std::cout << "TEST OUTPUT: 7.5: " << std::endl;
      if (reason < 0) {
        PetscPrintf(MPI_COMM_WORLD, "\nKSPCG solver Diverged;\n");
      }

      // Copy petsc_x to x
      if (mpi_rank == 0) {
        for (vi = 0; vi < dim; vi++) {
          VecGetValues(petsc_x, 1, &vi, &v);
          x(vi) = v;
        }
      }
      PetscFinalize();
    }

    std::cout << "TEST OUTPUT: 8: solution is: " << x << std::endl;

#endif

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return x;
}