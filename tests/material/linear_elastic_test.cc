#include <limits>

#include "Eigen/Dense"
#include "catch.hpp"

#include "material/linear_elastic.h"

//! \brief Check linearelastic class
TEST_CASE("LinearElastic is checked", "[material][linear_elastic]") {

  //! Check for id = 0
  SECTION("LinearElastic id is zero") {
    const unsigned id = 0;
    auto material = std::make_shared<mpm::LinearElastic>(id);
    REQUIRE(material->id() == 0);
  }

  SECTION("LinearElastic id is positive") {
    //! Check for id is a positive value
    const unsigned id = std::numeric_limits<unsigned>::max();
    auto material = std::make_shared<mpm::LinearElastic>(id);
    REQUIRE(material->id() == std::numeric_limits<unsigned>::max());
  }
}
