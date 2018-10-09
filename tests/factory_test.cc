#include "catch.hpp"
#include "spdlog/spdlog.h"

#include "factory.h"
#include "logger.h"

class ElementBase {
 public:
  // Default constructor
  ElementBase() = default;
  // Return class name
  virtual std::string name() const = 0;
};

class Quadrilateral : public ElementBase {
 public:
  // Default constructor
  Quadrilateral() = default;
  // Return class name
  std::string name() const override { return "Quadrilateral"; }
};

class Hexahedron : public ElementBase {
 public:
  // Default constructor
  Hexahedron() = default;
  // Return class name
  std::string name() const override { return "Hexahedron"; }
};

//! \brief Check factory class
TEST_CASE("Check factory class", "[factory]") {
  // Quadrilateral element
  static Register<ElementBase, Quadrilateral> quad("Quad");
  // Hexahedron element
  static Register<ElementBase, Hexahedron> hex("Hex");

  // Initialise logger
  auto console = spdlog::stdout_color_mt("factory");

  // Create a quadrilateral element
  REQUIRE(Factory<ElementBase>::instance()->check("Quad") == true);
  auto quad_element = Factory<ElementBase>::instance()->create("Quad");
  REQUIRE(quad_element->name() == "Quadrilateral");

  // Create a hexahedron element
  REQUIRE(Factory<ElementBase>::instance()->check("Hex") == true);
  auto hex_element = Factory<ElementBase>::instance()->create("Hex");
  REQUIRE(hex_element->name() == "Hexahedron");

  // Create a non-existant element
  REQUIRE(Factory<ElementBase>::instance()->check("New") == false);
  try {
    auto element = Factory<ElementBase>::instance()->create("New");
  } catch (std::exception& exception) {
    console->error("Factory error: {}", exception.what());
  }
}
