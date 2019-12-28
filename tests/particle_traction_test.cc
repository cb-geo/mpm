#include <memory>

#include "catch.hpp"

#include "function_base.h"
#include "linear_function.h"
#include "particle_traction.h"

//! Check particle traction
TEST_CASE("Particle traction is checked", "[load][traction]") {
  // Tolerance
  const double Tolerance = 1.E-9;

  SECTION("Check particle traction") {

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
    double traction = 10;
    auto ptraction = std::make_shared<mpm::ParticleTraction>(set_id, mfunction,
                                                             dir, traction);

    // Check particle set id
    REQUIRE(ptraction->setid() == set_id);
    // Check direction
    REQUIRE(ptraction->dir() == dir);
    // Check traction at x = 0
    REQUIRE(ptraction->traction(0.00) == Approx(0.).epsilon(Tolerance));
    // Check traction at x = 0.25
    REQUIRE(ptraction->traction(0.25) == Approx(5.).epsilon(Tolerance));
    // Check traction at x = 0.5
    REQUIRE(ptraction->traction(0.50) == Approx(10.).epsilon(Tolerance));
    // Check traction at x = 0.75
    REQUIRE(ptraction->traction(0.75) == Approx(10.).epsilon(Tolerance));
    // Check traction at x = 1.0
    REQUIRE(ptraction->traction(1.00) == Approx(10.).epsilon(Tolerance));
    // Check traction at x = 1.25
    REQUIRE(ptraction->traction(1.25) == Approx(5.).epsilon(Tolerance));
    // Check traction at x = 5
    REQUIRE(ptraction->traction(5.00) == Approx(0.).epsilon(Tolerance));
  }
}
