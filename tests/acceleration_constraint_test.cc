#include <memory>

#include "catch.hpp"

#include "acceleration_constraint.h"
#include "function_base.h"
#include "linear_function.h"

//! Check acceleration constraint
TEST_CASE("Acceleration constraint is checked", "[acceleration][constraint]") {
  // Tolerance
  const double Tolerance = 1.E-9;

  SECTION("Check acceleration constraint") {

    // Json property
    Json jfunctionproperties;
    jfunctionproperties["id"] = 0;
    std::vector<double> x_values{{0.0, 0.5, 1.0, 1.5}};
    std::vector<double> fx_values{{0.0, 1.0, 1.0, 0.0}};
    jfunctionproperties["xvalues"] = x_values;
    jfunctionproperties["fxvalues"] = fx_values;

    // math function
    std::shared_ptr<mpm::FunctionBase> mfunction =
        std::make_shared<mpm::LinearFunction>(0, jfunctionproperties);

    int set_id = -1;
    unsigned dir = 1;
    double acceleration = 10;
    auto nacceleration = std::make_shared<mpm::AccelerationConstraint>(
        set_id, mfunction, dir, acceleration);

    // Check particle set id
    REQUIRE(nacceleration->setid() == set_id);
    // Check direction
    REQUIRE(nacceleration->dir() == dir);
    // Check acceleration at x = 0
    REQUIRE(nacceleration->acceleration(0.00) == Approx(0.).epsilon(Tolerance));
    // Check acceleration at x = 0.25
    REQUIRE(nacceleration->acceleration(0.25) == Approx(5.).epsilon(Tolerance));
    // Check acceleration at x = 0.5
    REQUIRE(nacceleration->acceleration(0.50) ==
            Approx(10.).epsilon(Tolerance));
    // Check acceleration at x = 0.75
    REQUIRE(nacceleration->acceleration(0.75) ==
            Approx(10.).epsilon(Tolerance));
    // Check acceleration at x = 1.0
    REQUIRE(nacceleration->acceleration(1.00) ==
            Approx(10.).epsilon(Tolerance));
    // Check acceleration at x = 1.25
    REQUIRE(nacceleration->acceleration(1.25) == Approx(5.).epsilon(Tolerance));
    // Check acceleration at x = 5
    REQUIRE(nacceleration->acceleration(5.00) == Approx(0.).epsilon(Tolerance));
  }
}
