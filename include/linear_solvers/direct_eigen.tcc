//! Conjugate Gradient with default initial guess
template <typename Traits>
Eigen::VectorXd mpm::DirectEigen<Traits>::solve(
    const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b) {
  Eigen::VectorXd x;
  try {

    // Solver start
    auto solver_begin = std::chrono::steady_clock::now();
    if (verbosity_ > 0)
      console_->info("Type: \"{}\", Begin!", sub_solver_type_);

    if (verbosity_ == 3) {
      std::cout << "Coefficient Matrix A: " << A << std::endl;
      std::cout << "RHS Vector b: " << b << std::endl;
    }

    if (sub_solver_type_ == "lu") {
      Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

      solver.analyzePattern(A);
      solver.factorize(A);

      x = solver.solve(b);

      if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Fail to solve linear systems!\n");
      }
    } else if (sub_solver_type_ == "ldlt") {
      Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

      solver.analyzePattern(A);
      solver.factorize(A);

      x = solver.solve(b);

      if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Fail to solve linear systems!\n");
      }
    } else {
      throw std::runtime_error(
          "Sub solver type is not available! Available sub solver type "
          "implemented in DirectEigen class are: \"lu\" and \"ldlt\".\n");
    }

    // Solver End
    auto solver_end = std::chrono::steady_clock::now();
    if (verbosity_ > 0)
      console_->info("Type: \"{}\", End! Duration: {} ms.", sub_solver_type_,
                     std::chrono::duration_cast<std::chrono::milliseconds>(
                         solver_end - solver_begin)
                         .count());

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return x;
}