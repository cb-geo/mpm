//! Conjugate Gradient Solver
template <typename Traits>
bool lem::CGMKL<Traits>::solve() {
  bool convergence = false;
  double init_delta, new_delta;

  const MKL_INT n = this->vec_b_->size(), incr = 1;
  ptr_t<double> vs(n), vd(n), vq(n), vm(n);

  double* vec_x = this->vec_x_->data();
  double* vec_b = this->vec_b_->data();
  const double* vrestraints = this->vrestraints_.data();

  // Residual
  // (*vec_b_) -= mat_a * (*vec_x_);
  for (const auto& mat_a : *this->mat_a_) {
    // y := alpha*A*x + beta*y
    const double a = -1.0, b = 1.0;
    mkl_dcsrmv("N", &n, &n, &a, "GiiC", mat_a.valuePtr(), mat_a.innerIndexPtr(),
               mat_a.outerIndexPtr(), mat_a.outerIndexPtr() + 1, vec_x, &b,
               vec_b);
  }

  Eigen::Map<Eigen::VectorXd>(vm, n) = this->precondition_jacobian();

  // vd = vm * vec_b
  vdMul(n, vm, vec_b, vd);

  init_delta = new_delta = ddot(&n, vd, &incr, vec_b, &incr);

  unsigned i = 0;
  double alpha = 0., old_delta = 0.;
  for (i = 0; i < this->max_iter_; ++i) {
    // vq = mat_a_ * vd;
    memset(vq, 0, n * sizeof(double));
    for (const auto& mat_a : *this->mat_a_) {
      const double a = 1.0, b = 1.0;
      // y := alpha*A*x + beta*y
      mkl_dcsrmv("N", &n, &n, &a, "GiiC", mat_a.valuePtr(),
                 mat_a.innerIndexPtr(), mat_a.outerIndexPtr(),
                 mat_a.outerIndexPtr() + 1, vd, &b, vq);
    }
    // Apply restraints
    // vq = vq * vrestraints
    vdMul(n, vq, vrestraints, vq);

    // alpha = new_delta / vd.dot(vq)
    alpha = new_delta / ddot(&n, vd, &incr, vq, &incr);

    {  // vec_x += alpha * vd
      double a = alpha;
      daxpy(&n, &a, vd, &incr, vec_x, &incr);
    }
    {  // vec_b -= alpha * vq
      double a = -alpha;
      daxpy(&n, &a, vq, &incr, vec_b, &incr);
    }

    // vs = vm * vec_b
    vdMul(n, vm, vec_b, vs);

    old_delta = new_delta;
    new_delta = ddot(&n, vec_b, &incr, vs, &incr);

    {  // vd = vs + (new_delta / old_delta) * vd;
      double a = 1.0;
      double b = new_delta / old_delta;
      daxpby(&n, &a, vs, &incr, &b, vd, &incr);
    }

    // Break if convergence criterion is achieved
    if (new_delta <
        std::fabs(this->tolerance_ * this->tolerance_ * init_delta)) {
      convergence = true;
      break;
    }
  }

  // Update external nodal force in the system
  (this->vec_b_)->setZero();
  // (*vec_b_) = mat_a * (*vec_x_);
  for (const auto& mat_a : *this->mat_a_) {
    // y := alpha*A*x + beta*y
    const double a = 1., b = 1.;
    mkl_dcsrmv("N", &n, &n, &a, "GiiC", mat_a.valuePtr(), mat_a.innerIndexPtr(),
               mat_a.outerIndexPtr(), mat_a.outerIndexPtr() + 1, vec_x, &b,
               vec_b);
  }

  // Update delta and iterations
  this->delta_ = new_delta;
  this->niterations_ = i;

  console_->info("Iteration: {} of {}; delta: {}", i, this->max_iter_,
                 new_delta);

  return convergence;
}