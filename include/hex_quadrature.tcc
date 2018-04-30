// Default constructor
template <unsigned Tdim, unsigned Tnquadratures>
mpm::HexahedronQuadrature<Tdim, Tnquadratures>::HexahedronQuadrature()
    : mpm::QuadratureBase<Tdim, Tnquadratures>() {

  static_assert(Tdim == 3, "Invalid dimension for a 3D hexahedron element");
  static_assert(
      (Tnquadratures == 1) || (Tnquadratures == 8) || (Tnquadratures == 27),
      "Invalid number of quadratures");

  if (Tnquadratures == 1) {
    // Define quadratures
    qpoints_(0, 0) = 0.;
    qpoints_(0, 1) = 0.;
    qpoints_(0, 2) = 0.;

    weights_.at(0) = 8.;

  } else if (Tnquadratures == 8) {
    qpoints_(0, 0) = -1. / std::sqrt(3.);
    qpoints_(0, 1) = -1. / std::sqrt(3.);
    qpoints_(0, 2) = -1. / std::sqrt(3.);

    qpoints_(1, 0) = 1. / std::sqrt(3.);
    qpoints_(1, 1) = -1. / std::sqrt(3.);
    qpoints_(1, 2) = -1. / std::sqrt(3.);

    qpoints_(2, 0) = 1. / std::sqrt(3.);
    qpoints_(2, 1) = 1. / std::sqrt(3.);
    qpoints_(2, 2) = -1. / std::sqrt(3.);

    qpoints_(3, 0) = -1. / std::sqrt(3.);
    qpoints_(3, 1) = 1. / std::sqrt(3.);
    qpoints_(3, 2) = -1. / std::sqrt(3.);

    qpoints_(4, 0) = -1. / std::sqrt(3.);
    qpoints_(4, 1) = -1. / std::sqrt(3.);
    qpoints_(4, 2) = 1. / std::sqrt(3.);

    qpoints_(5, 0) = 1. / std::sqrt(3.);
    qpoints_(5, 1) = -1. / std::sqrt(3.);
    qpoints_(5, 2) = 1. / std::sqrt(3.);

    qpoints_(6, 0) = 1. / std::sqrt(3.);
    qpoints_(6, 1) = 1. / std::sqrt(3.);
    qpoints_(6, 2) = 1. / std::sqrt(3.);

    qpoints_(7, 0) = -1. / std::sqrt(3.);
    qpoints_(7, 1) = 1. / std::sqrt(3.);
    qpoints_(7, 2) = 1. / std::sqrt(3.);

    weights_.at(0) = 1.;
    weights_.at(1) = 1.;
    weights_.at(2) = 1.;
    weights_.at(3) = 1.;
    weights_.at(4) = 1.;
    weights_.at(5) = 1.;
    weights_.at(6) = 1.;
    weights_.at(7) = 1.;

  } else if (Tnquadratures == 27) {
    const double Qpoint_a = std::sqrt(3. / 5.);
    const double Weight_b = 5. / 9.;
    const double Weight_c = 8. / 9.;

    // Along xi(0) = -1 plane
    // Line xi(1) = -1
    qpoints_(0, 0) = -Qpoint_a;
    qpoints_(0, 1) = -Qpoint_a;
    qpoints_(0, 2) = -Qpoint_a;

    qpoints_(1, 0) = -Qpoint_a;
    qpoints_(1, 1) = -Qpoint_a;
    qpoints_(1, 2) = 0.;

    qpoints_(2, 0) = -Qpoint_a;
    qpoints_(2, 1) = -Qpoint_a;
    qpoints_(2, 2) = Qpoint_a;

    // Line xi(1) = 0
    qpoints_(3, 0) = -Qpoint_a;
    qpoints_(3, 1) = 0.;
    qpoints_(3, 2) = -Qpoint_a;

    qpoints_(4, 0) = -Qpoint_a;
    qpoints_(4, 1) = 0.;
    qpoints_(4, 2) = 0.;

    qpoints_(5, 0) = -Qpoint_a;
    qpoints_(5, 1) = 0.;
    qpoints_(5, 2) = Qpoint_a;

    // Line xi(1) = +1
    qpoints_(6, 0) = -Qpoint_a;
    qpoints_(6, 1) = Qpoint_a;
    qpoints_(6, 2) = -Qpoint_a;

    qpoints_(7, 0) = -Qpoint_a;
    qpoints_(7, 1) = Qpoint_a;
    qpoints_(7, 2) = 0.;

    qpoints_(8, 0) = -Qpoint_a;
    qpoints_(8, 1) = Qpoint_a;
    qpoints_(8, 2) = Qpoint_a;

    // Along xi(0) = 0 plane
    // Line xi(1) = -1
    qpoints_(9, 0) = 0.;
    qpoints_(9, 1) = -Qpoint_a;
    qpoints_(9, 2) = -Qpoint_a;

    qpoints_(10, 0) = 0.;
    qpoints_(10, 1) = -Qpoint_a;
    qpoints_(10, 2) = 0.;

    qpoints_(11, 0) = 0.;
    qpoints_(11, 1) = -Qpoint_a;
    qpoints_(11, 2) = Qpoint_a;

    // Line xi(1) = 0
    qpoints_(12, 0) = 0.;
    qpoints_(12, 1) = 0.;
    qpoints_(12, 2) = -Qpoint_a;

    qpoints_(13, 0) = 0.;
    qpoints_(13, 1) = 0.;
    qpoints_(13, 2) = 0.;

    qpoints_(14, 0) = 0.;
    qpoints_(14, 1) = 0.;
    qpoints_(14, 2) = Qpoint_a;

    // Line xi(1) = +1
    qpoints_(15, 0) = 0.;
    qpoints_(15, 1) = Qpoint_a;
    qpoints_(15, 2) = -Qpoint_a;

    qpoints_(16, 0) = 0.;
    qpoints_(16, 1) = Qpoint_a;
    qpoints_(16, 2) = 0.;

    qpoints_(17, 0) = 0.;
    qpoints_(17, 1) = Qpoint_a;
    qpoints_(17, 2) = Qpoint_a;

    // Along xi(0) = +1 plane
    // Line xi(1) = -1
    qpoints_(18, 0) = Qpoint_a;
    qpoints_(18, 1) = -Qpoint_a;
    qpoints_(18, 2) = -Qpoint_a;

    qpoints_(19, 0) = Qpoint_a;
    qpoints_(19, 1) = -Qpoint_a;
    qpoints_(19, 2) = 0.;

    qpoints_(20, 0) = Qpoint_a;
    qpoints_(20, 1) = -Qpoint_a;
    qpoints_(20, 2) = Qpoint_a;

    // Line xi(1) = 0
    qpoints_(21, 0) = Qpoint_a;
    qpoints_(21, 1) = 0.;
    qpoints_(21, 2) = -Qpoint_a;

    qpoints_(22, 0) = Qpoint_a;
    qpoints_(22, 1) = 0.;
    qpoints_(22, 2) = 0.;

    qpoints_(23, 0) = Qpoint_a;
    qpoints_(23, 1) = 0.;
    qpoints_(23, 2) = Qpoint_a;

    // Line xi(1) = +1
    qpoints_(24, 0) = Qpoint_a;
    qpoints_(24, 1) = Qpoint_a;
    qpoints_(24, 2) = -Qpoint_a;

    qpoints_(25, 0) = Qpoint_a;
    qpoints_(25, 1) = Qpoint_a;
    qpoints_(25, 2) = 0.;

    qpoints_(26, 0) = Qpoint_a;
    qpoints_(26, 1) = Qpoint_a;
    qpoints_(26, 2) = Qpoint_a;

    // Weights
    weights_.at(0) = Weight_b * Weight_b * Weight_b;
    weights_.at(1) = Weight_b * Weight_b * Weight_c;
    weights_.at(2) = Weight_b * Weight_b * Weight_b;
    weights_.at(3) = Weight_b * Weight_c * Weight_b;
    weights_.at(4) = Weight_b * Weight_c * Weight_c;
    weights_.at(5) = Weight_b * Weight_c * Weight_b;
    weights_.at(6) = Weight_b * Weight_b * Weight_b;
    weights_.at(7) = Weight_b * Weight_b * Weight_c;
    weights_.at(8) = Weight_b * Weight_b * Weight_b;

    weights_.at(9) = Weight_c * Weight_b * Weight_b;
    weights_.at(10) = Weight_c * Weight_b * Weight_c;
    weights_.at(11) = Weight_c * Weight_b * Weight_b;
    weights_.at(12) = Weight_c * Weight_c * Weight_b;
    weights_.at(13) = Weight_c * Weight_c * Weight_c;
    weights_.at(14) = Weight_c * Weight_c * Weight_b;
    weights_.at(15) = Weight_c * Weight_b * Weight_b;
    weights_.at(16) = Weight_c * Weight_b * Weight_c;
    weights_.at(17) = Weight_c * Weight_b * Weight_b;

    weights_.at(18) = Weight_b * Weight_b * Weight_b;
    weights_.at(19) = Weight_b * Weight_b * Weight_c;
    weights_.at(20) = Weight_b * Weight_b * Weight_b;
    weights_.at(21) = Weight_b * Weight_c * Weight_b;
    weights_.at(22) = Weight_b * Weight_c * Weight_c;
    weights_.at(23) = Weight_b * Weight_c * Weight_b;
    weights_.at(24) = Weight_b * Weight_b * Weight_b;
    weights_.at(25) = Weight_b * Weight_b * Weight_c;
    weights_.at(26) = Weight_b * Weight_b * Weight_b;
  }
}