#include <cmath>
#include <limits>
#include <memory>
#include <vector>

#include "catch.hpp"
#include "json.hpp"

#include "function_base.h"
#include "sin_function.h"

//! \brief Check Functions class
TEST_CASE("Sin function is checked", "[sinfn]") {

  // Tolerance
  const double Tolerance = 1.E-9;

  // Sin function properties
  unsigned id = 0;

  // Json property
  Json jfunctionproperties;
  jfunctionproperties["id"] = 0;

  SECTION("Check incorrect sin function initialisation") {
    bool status = true;
    try {
      std::shared_ptr<mpm::FunctionBase> sinfn =
          std::make_shared<mpm::SinFunction>(id, jfunctionproperties);
    } catch (std::exception& exception) {
      status = false;
    }
    REQUIRE(status == false);
  }

  SECTION("Check correct sin function initialisation") {
    bool status = true;
    try {
      // Initialize
      jfunctionproperties["x0"] = 0.0;
      jfunctionproperties["a"] = 2.0;
      std::shared_ptr<mpm::FunctionBase> sinfn =
          std::make_shared<mpm::SinFunction>(id, jfunctionproperties);
    } catch (std::exception& exception) {
      status = false;
    }
    REQUIRE(status == true);
  }

  SECTION("Check sin function for x0 = 0") {
    // Initialize
    jfunctionproperties["x0"] = 0.0;
    jfunctionproperties["a"] = 2.0;
    std::shared_ptr<mpm::FunctionBase> sinfn =
        std::make_shared<mpm::SinFunction>(id, jfunctionproperties);
    // check id
    REQUIRE(sinfn->id() == id);

    // check values for different x values
    double x = 0.0;
    REQUIRE(sinfn->value(x) == Approx(0.0).epsilon(Tolerance));
    x = 0.5;
    REQUIRE(sinfn->value(x) == Approx(0.8414709848).epsilon(Tolerance));
    x = 1.0;
    REQUIRE(sinfn->value(x) == Approx(0.9092974268).epsilon(Tolerance));
    x = 2.5;
    REQUIRE(sinfn->value(x) == Approx(-0.9589242747).epsilon(Tolerance));
  }

  SECTION("Check sin function for x0 = 0.5") {
    // Initialize
    jfunctionproperties["x0"] = 0.5;
    jfunctionproperties["a"] = 2.0;
    std::shared_ptr<mpm::FunctionBase> sinfn =
        std::make_shared<mpm::SinFunction>(id, jfunctionproperties);
    // check id
    REQUIRE(sinfn->id() == id);

    // check values for different x values
    double x = 0.0;
    REQUIRE(sinfn->value(x) == Approx(-0.8414709848).epsilon(Tolerance));
    x = 0.5;
    REQUIRE(sinfn->value(x) == Approx(0.0).epsilon(Tolerance));
    x = 1.0;
    REQUIRE(sinfn->value(x) == Approx(0.8414709848).epsilon(Tolerance));
    x = 2.5;
    REQUIRE(sinfn->value(x) == Approx(-0.7568024953).epsilon(Tolerance));
  }
}
