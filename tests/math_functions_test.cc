#include <cmath>
#include <limits>
#include <memory>
#include <vector>

#include "catch.hpp"
#include "json.hpp"

#include "functions/function_base.h"
#include "functions/step_function.h"

//! \brief Check Functions class
TEST_CASE("Linear step function is checked", "[stepfunction]") {

  // Tolerance
  const double Tolerance = 1.E-9;

  // Step function properties
  unsigned id = 0;
  std::vector<double> x_values{{0.0, 0.5, 1.0, 1.5}};
  std::vector<double> fx_values{{0.0, 1.0, 1.0, 0.0}};

  // Json property
  Json jfunctionproperties;
  jfunctionproperties["id"] = id;
  jfunctionproperties["xvalues"] = x_values;
  jfunctionproperties["fxvalues"] = fx_values;

  SECTION("Check incorrect step function initialisation") {
    bool status = true;
    fx_values.emplace_back(0.0);
    jfunctionproperties["fxvalues"] = fx_values;
    try {
      std::shared_ptr<mpm::FunctionBase> stepfunctionptr =
          std::make_shared<mpm::StepFunction>(id, jfunctionproperties);
    } catch (std::exception& exception) {
      status = false;
    }
    REQUIRE(status == false);
  }

  SECTION("Check correct step function initialisation") {
    bool status = true;
    try {
      std::shared_ptr<mpm::FunctionBase> stepfunctionptr =
          std::make_shared<mpm::StepFunction>(id, jfunctionproperties);
    } catch (std::exception& exception) {
      status = false;
    }
    REQUIRE(status == true);
  }

  SECTION("Check step function") {
    std::shared_ptr<mpm::FunctionBase> stepfunctionptr =
        std::make_shared<mpm::StepFunction>(id, jfunctionproperties);
    // check id
    REQUIRE(stepfunctionptr->id() == id);

    // check values for different x values
    double x = 0.0;
    REQUIRE(stepfunctionptr->value(x) == Approx(0.0).epsilon(Tolerance));
    x = 0.0000001;
    REQUIRE(stepfunctionptr->value(x) == Approx(0.0000002).epsilon(Tolerance));
    x = 0.4999999;
    REQUIRE(stepfunctionptr->value(x) == Approx(0.9999998).epsilon(Tolerance));
    x = 0.5;
    REQUIRE(stepfunctionptr->value(x) == Approx(1.0).epsilon(Tolerance));
    x = 0.9999999;
    REQUIRE(stepfunctionptr->value(x) == Approx(1.0).epsilon(Tolerance));
    x = 1;
    REQUIRE(stepfunctionptr->value(x) == Approx(1.0).epsilon(Tolerance));
    x = 1.2;
    REQUIRE(stepfunctionptr->value(x) == Approx(0.6).epsilon(Tolerance));
    x = 1.5;
    REQUIRE(stepfunctionptr->value(x) == Approx(0.0).epsilon(Tolerance));
    x = 1.50000001;
    REQUIRE(stepfunctionptr->value(x) == Approx(0.0).epsilon(Tolerance));
    x = 10.0;
    REQUIRE(stepfunctionptr->value(x) == Approx(0.0).epsilon(Tolerance));
  }
}
