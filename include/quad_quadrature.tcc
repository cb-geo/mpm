// Default constructor
template <unsigned Tdim, unsigned Tnquadratures>
mpm::QuadrilateralQuadrature<Tdim, Tnquadratures>::QuadrilateralQuadrature()
    : mpm::QuadratureBase<Tdim, Tnquadratures>() {

  static_assert(Tdim == 2, "Invalid dimension for a quadrilateral element");
  static_assert(
      ((Tnquadratures == 1) || (Tnquadratures == 4) || (Tnquadratures == 9)),
      "Invalid number of quadratures");

  if (Tnquadratures == 1) {
    // Define quadratures
    qpoints_(0, 0) = 0.;
    qpoints_(0, 1) = 0.;

    weights_.at(0) = 8.;

  } else if (Tnquadratures == 4) {
    // Define quadratures
    const double Val_1_by_sqrt3 = 1. / std::sqrt(3.);

    qpoints_(0, 0) = -Val_1_by_sqrt3;
    qpoints_(0, 1) = -Val_1_by_sqrt3;
    qpoints_(1, 0) = Val_1_by_sqrt3;
    qpoints_(1, 1) = -Val_1_by_sqrt3;
    qpoints_(2, 0) = Val_1_by_sqrt3;
    qpoints_(2, 1) = Val_1_by_sqrt3;
    qpoints_(3, 0) = -Val_1_by_sqrt3;
    qpoints_(3, 1) = Val_1_by_sqrt3;

    weights_.at(0) = 1.;
    weights_.at(1) = 1.;
    weights_.at(2) = 1.;
    weights_.at(3) = 1.;

  } else if (Tnquadratures == 9) {

    const double Val_sqrt_3by5 = std::sqrt(3. / 5.);
    const double Val_25_by_81 = 25. / 81.;
    const double Val_40_by_81 = 40. / 81.;
    const double Val_64_by_81 = 64. / 81.;

    qpoints_(0, 0) = -Val_sqrt_3by5;
    qpoints_(0, 1) = -Val_sqrt_3by5;

    qpoints_(1, 0) = Val_sqrt_3by5;
    qpoints_(1, 1) = -Val_sqrt_3by5;

    qpoints_(2, 0) = Val_sqrt_3by5;
    qpoints_(2, 1) = Val_sqrt_3by5;

    qpoints_(3, 0) = -Val_sqrt_3by5;
    qpoints_(3, 1) = Val_sqrt_3by5;

    qpoints_(4, 0) = 0;
    qpoints_(4, 1) = -Val_sqrt_3by5;

    qpoints_(5, 0) = Val_sqrt_3by5;
    qpoints_(5, 1) = 0;

    qpoints_(6, 0) = 0;
    qpoints_(6, 1) = Val_sqrt_3by5;

    qpoints_(7, 0) = -Val_sqrt_3by5;
    qpoints_(7, 1) = 0;

    qpoints_(8, 0) = 0;
    qpoints_(8, 1) = 0;

    weights_.at(0) = Val_25_by_81;
    weights_.at(1) = Val_25_by_81;
    weights_.at(2) = Val_25_by_81;
    weights_.at(3) = Val_25_by_81;
    weights_.at(4) = Val_40_by_81;
    weights_.at(5) = Val_40_by_81;
    weights_.at(6) = Val_40_by_81;
    weights_.at(7) = Val_40_by_81;
    weights_.at(8) = Val_64_by_81;
  }
}