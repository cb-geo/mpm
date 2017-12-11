#include <limits>

#include "Eigen/Dense"
#include "catch.hpp"

#include "material/material.h"

//! \brief Check material class
TEST_CASE("Material is checked", "[material]") {

  //! Check for id = 0
  SECTION("Material id is zero") {
    const unsigned id = 0;
    auto material = std::make_shared<mpm::Material>(id);
    REQUIRE(material->id() == 0);
  }

  SECTION("Material id is positive") {
    //! Check for id is a positive value
    const unsigned id = std::numeric_limits<unsigned>::max();
    auto material = std::make_shared<mpm::Material>(id);
    REQUIRE(material->id() == std::numeric_limits<unsigned>::max());
  }
}
