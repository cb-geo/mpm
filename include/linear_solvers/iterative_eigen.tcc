//! Conjugate Gradient with default initial guess
template <typename Traits>
Eigen::VectorXd mpm::IterativeEigen<Traits>::solve(
    const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b) {
  Eigen::VectorXd x;
  try {

    // Solver start
    auto solver_begin = std::chrono::steady_clock::now();
    if (verbosity_ > 0)
      console_->info("Type: \"{}\", Preconditioner: \"{}\", Begin!",
                     sub_solver_type_, preconditioner_type_);

    if (verbosity_ == 3) {
      std::cout << "Coefficient Matrix A: " << A << std::endl;
      std::cout << "RHS Vector b: " << b << std::endl;
    }

    if (sub_solver_type_ == "cg") {
      if (preconditioner_type_ == "none") {
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;

        solver.setMaxIterations(max_iter_);
        solver.setTolerance(tolerance_);
        solver.compute(A);

        x = solver.solve(b);

        if (verbosity_ >= 1) {
          std::cout << "#iterations:     " << solver.iterations() << std::endl;
          std::cout << "estimated error: " << solver.error() << std::endl;
        }

        if (solver.info() != Eigen::Success) {
          throw std::runtime_error("Fail to solve linear systems!\n");
        }
      } else if (preconditioner_type_ == "icc") {
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                                 Eigen::Lower | Eigen::Upper,
                                 Eigen::IncompleteCholesky<double>>
            solver;

        solver.setMaxIterations(max_iter_);
        solver.setTolerance(tolerance_);
        solver.compute(A);

        x = solver.solve(b);

        if (verbosity_ >= 1) {
          std::cout << "#iterations:     " << solver.iterations() << std::endl;
          std::cout << "estimated error: " << solver.error() << std::endl;
        }

        if (solver.info() != Eigen::Success) {
          throw std::runtime_error("Fail to solve linear systems!\n");
        }
      } else {
        throw std::runtime_error("Preconditioner type is unavailable!\n");
      }

    } else if (sub_solver_type_ == "lscg") {
      Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solver;

      solver.setMaxIterations(max_iter_);
      solver.setTolerance(tolerance_);
      solver.compute(A);

      x = solver.solve(b);

      if (verbosity_ >= 1) {
        std::cout << "#iterations:     " << solver.iterations() << std::endl;
        std::cout << "estimated error: " << solver.error() << std::endl;
      }

      if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Fail to solve linear systems!\n");
      }

    } else if (sub_solver_type_ == "bicgstab") {
      Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;

      solver.setMaxIterations(max_iter_);
      solver.setTolerance(tolerance_);
      solver.compute(A);

      x = solver.solve(b);

      if (verbosity_ >= 1) {
        std::cout << "#iterations:     " << solver.iterations() << std::endl;
        std::cout << "estimated error: " << solver.error() << std::endl;
      }

      if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Fail to solve linear systems!\n");
      }
    }

    // Solver End
    auto solver_end = std::chrono::steady_clock::now();
    if (verbosity_ > 0)
      console_->info(
          "Type: \"{}\", Preconditioner: \"{}\", End! Duration: {} ms.",
          sub_solver_type_, preconditioner_type_,
          std::chrono::duration_cast<std::chrono::milliseconds>(solver_end -
                                                                solver_begin)
              .count());

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return x;
}