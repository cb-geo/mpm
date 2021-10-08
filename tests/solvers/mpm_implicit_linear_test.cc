#include "catch.hpp"

//! Alias for JSON
#include "json.hpp"
using Json = nlohmann::json;

#include "mpm_implicit_linear.h"
#include "write_mesh_particles.h"

// Check MPM Implicit Linear
TEST_CASE("MPM 2D Implicit Linear implementation is checked",
          "[MPM][2D][Implicit][1Phase]") {
  // Dimension
  const unsigned Dim = 2;

  // Write JSON file
  const std::string fname = "mpm-implicit-linear";
  const std::string analysis = "MPMImplicitLinear2D";
  const std::string mpm_scheme = "newmark";
  const std::string lin_solver_type = "IterativeEigen";
  bool resume = false;
  REQUIRE(mpm_test::write_json_implicit_linear(2, resume, analysis, mpm_scheme,
                                               fname, lin_solver_type) == true);

  // Write JSON Entity Sets file
  REQUIRE(mpm_test::write_entity_set() == true);

  // Write Mesh
  REQUIRE(mpm_test::write_mesh_2d() == true);

  // Write Particles
  REQUIRE(mpm_test::write_particles_2d() == true);

  // Assign argc and argv to input arguments of MPM
  int argc = 5;
  // clang-format off
  char* argv[] = {(char*)"./mpm",
                  (char*)"-f",  (char*)"./",
                  (char*)"-i",  (char*)"mpm-implicit-linear-2d.json"};
  // clang-format on

  SECTION("Check initialisation") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run Implicit Linear MPM
    auto mpm = std::make_unique<mpm::MPMImplicitLinear<Dim>>(std::move(io));

    // Initialise materials
    REQUIRE_NOTHROW(mpm->initialise_materials());
    // Initialise mesh
    REQUIRE_NOTHROW(mpm->initialise_mesh());
    // Initialise particles
    REQUIRE_NOTHROW(mpm->initialise_particles());
    // Initialise external loading
    REQUIRE_NOTHROW(mpm->initialise_loads());

    // Renitialise materials
    REQUIRE_THROWS(mpm->initialise_materials());
  }

  SECTION("Check solver") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run Implicit Linear MPM
    auto mpm = std::make_unique<mpm::MPMImplicitLinear<Dim>>(std::move(io));
    // Solve
    REQUIRE(mpm->solve() == true);
    // Test check point restart
    REQUIRE(mpm->checkpoint_resume() == false);
  }

  SECTION("Check resume") {
    // Write JSON file
    const std::string fname = "mpm-implicit-linear";
    const std::string analysis = "MPMImplicitLinear2D";
    const std::string mpm_scheme = "newmark";
    const std::string lin_solver_type = "IterativeEigen";
    bool resume = true;
    REQUIRE(mpm_test::write_json_implicit_linear(2, resume, analysis,
                                                 mpm_scheme, fname,
                                                 lin_solver_type) == true);

    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run Implicit Linear MPM
    auto mpm = std::make_unique<mpm::MPMImplicitLinear<Dim>>(std::move(io));

    // Initialise materials
    REQUIRE_NOTHROW(mpm->initialise_materials());
    // Initialise mesh
    REQUIRE_NOTHROW(mpm->initialise_mesh());

    // Test check point restart
    REQUIRE(mpm->checkpoint_resume() == true);
    {
      // Create an IO object
      auto io = std::make_unique<mpm::IO>(argc, argv);
      // Run Implicit Linear MPM
      auto mpm_resume =
          std::make_unique<mpm::MPMImplicitLinear<Dim>>(std::move(io));
      REQUIRE(mpm_resume->solve() == true);
    }
  }
}

// Check MPM Implicit Linear
TEST_CASE("MPM 3D Implicit Linear implementation is checked",
          "[MPM][3D][Implicit][1Phase]") {
  // Dimension
  const unsigned Dim = 3;

  // Write JSON file
  const std::string fname = "mpm-implicit-linear";
  const std::string analysis = "MPMImplicitLinear3D";
  const std::string mpm_scheme = "newmark";
  const std::string lin_solver_type = "IterativeEigen";
  const bool resume = false;
  REQUIRE(mpm_test::write_json_implicit_linear(3, resume, analysis, mpm_scheme,
                                               fname, lin_solver_type) == true);

  // Write JSON Entity Sets file
  REQUIRE(mpm_test::write_entity_set() == true);

  // Write Mesh
  REQUIRE(mpm_test::write_mesh_3d() == true);

  // Write Particles
  REQUIRE(mpm_test::write_particles_3d() == true);

  // Assign argc and argv to input arguments of MPM
  int argc = 5;
  // clang-format off
  char* argv[] = {(char*)"./mpm",
                  (char*)"-f",  (char*)"./",
                  (char*)"-i",  (char*)"mpm-implicit-linear-3d.json"};
  // clang-format on

  SECTION("Check initialisation") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run Implicit Linear MPM
    auto mpm = std::make_unique<mpm::MPMImplicitLinear<Dim>>(std::move(io));

    // Initialise materials
    REQUIRE_NOTHROW(mpm->initialise_materials());
    // Initialise mesh
    REQUIRE_NOTHROW(mpm->initialise_mesh());
    // Initialise particles
    REQUIRE_NOTHROW(mpm->initialise_particles());

    // Renitialise materials
    REQUIRE_THROWS(mpm->initialise_materials());
  }

  SECTION("Check solver") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run Implicit Linear MPM
    auto mpm = std::make_unique<mpm::MPMImplicitLinear<Dim>>(std::move(io));
    // Solve
    REQUIRE(mpm->solve() == true);
    // Test check point restart
    REQUIRE(mpm->checkpoint_resume() == false);
  }

  SECTION("Check resume") {
    // Write JSON file
    const std::string fname = "mpm-implicit-linear";
    const std::string analysis = "MPMImplicitLinear3D";
    const std::string mpm_scheme = "newmark";
    const std::string lin_solver_type = "IterativeEigen";
    bool resume = true;
    REQUIRE(mpm_test::write_json_implicit_linear(3, resume, analysis,
                                                 mpm_scheme, fname,
                                                 lin_solver_type) == true);

    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run Implicit Linear MPM
    auto mpm = std::make_unique<mpm::MPMImplicitLinear<Dim>>(std::move(io));

    // Initialise materials
    REQUIRE_NOTHROW(mpm->initialise_materials());
    // Initialise mesh
    REQUIRE_NOTHROW(mpm->initialise_mesh());

    // Test check point restart
    REQUIRE(mpm->checkpoint_resume() == true);
    {
      // Solve
      auto io = std::make_unique<mpm::IO>(argc, argv);
      // Run Implicit Linear MPM
      auto mpm_resume =
          std::make_unique<mpm::MPMImplicitLinear<Dim>>(std::move(io));
      REQUIRE(mpm_resume->solve() == true);
    }
  }
}
