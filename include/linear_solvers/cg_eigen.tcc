//! Conjugate Gradient with default initial guess
template <typename Traits>
Eigen::VectorXd mpm::CGEigen<Traits>::solve(
    const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b,
    std::string solver_type) {
  Eigen::VectorXd x;
  try {

    if (solver_type == "cg") {
      Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;

      solver.setMaxIterations(max_iter_);
      solver.setTolerance(tolerance_);
      solver.compute(A);

      x = solver.solve(b);

      if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Fail to solve linear systems!\n");
      }

    } else if (solver_type == "lscg") {

      // Another option is LDLT, but not accurate as lscg
      // Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
      Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solver;

      solver.setMaxIterations(max_iter_);
      solver.setTolerance(tolerance_);
      solver.compute(A);

      x = solver.solve(b);

      if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Fail to solve linear systems!\n");
      }
    }

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return x;
}