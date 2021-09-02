#include "catch.hpp"

//! Alias for JSON
#include "json.hpp"
using Json = nlohmann::json;

#include "mpm_semi_implicit_twophase.h"
#include "write_mesh_particles.h"

// Check MPM Semi-implicit TwoPhase
TEST_CASE("MPM 2D Semi-implicit TwoPhase implementation is checked",
          "[MPM][2D][Semi-implicit][2Phase]") {
  // Dimension
  const unsigned Dim = 2;

  // Write JSON file
  const std::string fname = "mpm-semi-implicit-twophase";
  const std::string analysis = "MPMSemiImplicitTwoPhase2D";
  const std::string mpm_scheme = "usf";
  const std::string fsd_type = "density";
  const std::string lin_solver_type = "IterativeEigen";
  bool resume = false;
  REQUIRE(mpm_test::write_json_twophase(2, resume, analysis, mpm_scheme, fname,
                                        fsd_type, lin_solver_type) == true);

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
                  (char*)"-i",  (char*)"mpm-semi-implicit-twophase-2d.json"};
  // clang-format on

  SECTION("Check initialisation") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run Semi Implicit MPM
    auto mpm =
        std::make_unique<mpm::MPMSemiImplicitTwoPhase<Dim>>(std::move(io));

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
    // Run Semi Implicit MPM
    auto mpm =
        std::make_unique<mpm::MPMSemiImplicitTwoPhase<Dim>>(std::move(io));
    // Solve
    REQUIRE(mpm->solve() == true);
    // Test check point restart
    REQUIRE(mpm->checkpoint_resume() == false);
  }

  SECTION("Check resume") {
    // Write JSON file
    const std::string fname = "mpm-semi-implicit-twophase";
    const std::string analysis = "MPMSemiImplicitTwoPhase2D";
    const std::string mpm_scheme = "usf";
    const std::string fsd_type = "density";
    const std::string lin_solver_type = "IterativeEigen";
    bool resume = true;
    REQUIRE(mpm_test::write_json_twophase(2, resume, analysis, mpm_scheme,
                                          fname, fsd_type,
                                          lin_solver_type) == true);

    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run Semi Implicit MPM
    auto mpm =
        std::make_unique<mpm::MPMSemiImplicitTwoPhase<Dim>>(std::move(io));

    // Initialise materials
    REQUIRE_NOTHROW(mpm->initialise_materials());
    // Initialise mesh
    REQUIRE_NOTHROW(mpm->initialise_mesh());

    // Test check point restart
    REQUIRE(mpm->checkpoint_resume() == true);
    {
      // Solve
      auto io = std::make_unique<mpm::IO>(argc, argv);
      // Run explicit MPM
      auto mpm_resume =
          std::make_unique<mpm::MPMSemiImplicitTwoPhase<Dim>>(std::move(io));
      REQUIRE(mpm_resume->solve() == true);
    }
  }

  SECTION("Check pressure smoothing") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run Semi Implicit MPM
    auto mpm =
        std::make_unique<mpm::MPMSemiImplicitTwoPhase<Dim>>(std::move(io));
    // Pressure smoothing
    REQUIRE_NOTHROW(mpm->pressure_smoothing(mpm::ParticlePhase::Solid));
    REQUIRE_NOTHROW(mpm->pressure_smoothing(mpm::ParticlePhase::Liquid));
  }
}

// Check MPM Semi Implicit
TEST_CASE("MPM 3D Semi-implicit TwoPhase implementation is checked",
          "[MPM][3D][Semi-implicit][2Phase]") {
  // Dimension
  const unsigned Dim = 3;

  // Write JSON file
  const std::string fname = "mpm-semi-implicit-twophase";
  const std::string analysis = "MPMSemiImplicitTwoPhase3D";
  const std::string mpm_scheme = "usf";
  const std::string fsd_type = "density";
  const std::string lin_solver_type = "IterativeEigen";
  const bool resume = false;
  REQUIRE(mpm_test::write_json_twophase(3, resume, analysis, mpm_scheme, fname,
                                        fsd_type, lin_solver_type) == true);

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
                  (char*)"-i",  (char*)"mpm-semi-implicit-twophase-3d.json"};
  // clang-format on

  SECTION("Check initialisation") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run Semi Implicit MPM
    auto mpm =
        std::make_unique<mpm::MPMSemiImplicitTwoPhase<Dim>>(std::move(io));

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
    // Run Semi Implicit MPM
    auto mpm =
        std::make_unique<mpm::MPMSemiImplicitTwoPhase<Dim>>(std::move(io));
    // Solve
    REQUIRE(mpm->solve() == true);
    // Test check point restart
    REQUIRE(mpm->checkpoint_resume() == false);
  }

  SECTION("Check resume") {
    // Write JSON file
    const std::string fname = "mpm-semi-implicit-twophase";
    const std::string analysis = "MPMSemiImplicitTwoPhase3D";
    const std::string mpm_scheme = "usf";
    const std::string fsd_type = "density";
    const std::string lin_solver_type = "IterativeEigen";
    bool resume = true;
    REQUIRE(mpm_test::write_json_twophase(3, resume, analysis, mpm_scheme,
                                          fname, fsd_type,
                                          lin_solver_type) == true);

    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run Semi Implicit MPM
    auto mpm =
        std::make_unique<mpm::MPMSemiImplicitTwoPhase<Dim>>(std::move(io));

    // Initialise materials
    REQUIRE_NOTHROW(mpm->initialise_materials());
    // Initialise mesh
    REQUIRE_NOTHROW(mpm->initialise_mesh());

    // Test check point restart
    REQUIRE(mpm->checkpoint_resume() == true);
    {
      // Solve
      auto io = std::make_unique<mpm::IO>(argc, argv);
      // Run explicit MPM
      auto mpm_resume =
          std::make_unique<mpm::MPMSemiImplicitTwoPhase<Dim>>(std::move(io));
      REQUIRE(mpm_resume->solve() == true);
    }
  }

  SECTION("Check pressure smoothing") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run Semi Implicit MPM
    auto mpm =
        std::make_unique<mpm::MPMSemiImplicitTwoPhase<Dim>>(std::move(io));
    // Pressure smoothing
    REQUIRE_NOTHROW(mpm->pressure_smoothing(mpm::ParticlePhase::Solid));
    REQUIRE_NOTHROW(mpm->pressure_smoothing(mpm::ParticlePhase::Liquid));
  }
}
