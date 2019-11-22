#include "catch.hpp"

//! Alias for JSON
#include "json.hpp"
using Json = nlohmann::json;

#include "io.h"

#include "map.h"

// Check IO for input string
TEST_CASE("IO is checked for input parsing", "[IO][JSON]") {

  // Check input JSON
  SECTION("Check input JSON object") {
    // Make json object with input files
    Json json_file = {
        {"title", "Example JSON Input for MPM"},
        {"input_files",
         {{"config", "mpm.json"},
          {"mesh", "mesh-3d.txt"},
          {"constraints", "mesh_constraints.txt"},
          {"particles", "particles-3d.txt"},
          {"initial_stresses", "initial_soil_stress.txt"},
          {"materials", "materials.txt"},
          {"traction", "traction.txt"},
          {"entity_sets", "entity_sets.json"}}},
        {"mesh", {{"mesh_reader", "Ascii3D"}, {"cell_type", "ED3H8"}}},
        {"particle",
         {{"material_id", "0"},
          {"particle_type", "P3D"},
          {"particle_sets",
           {{{"set_id", {0}}, {"material_id", 0}},
            {{"set_id", {1}}, {"material_id", 1}}}}}},
        {"materials",
         {{{"id", 0},
           {"type", "LinearElastic"},
           {"density", 1000.},
           {"youngs_modulus", 1.0E+8},
           {"poisson_ratio", 0.495}},
          {{"id", 1},
           {"type", "LinearElastic"},
           {"density", 2300.},
           {"youngs_modulus", 1.5E+6},
           {"poisson_ratio", 0.25}}}},
        {"analysis",
         {{"type", "MPMExplicit3D"},
          {"stress_update", "usf"},
          {"dt", 0.001},
          {"nsteps", 1000},
          {"gravity", true},
          {"boundary_friction", 0.5},
          {"damping", {{"damping", true}, {"damping_ratio", 0.02}}},
          {"newmark", {{"newmark", true}, {"gamma", 0.5}, {"beta", 0.25}}}}},
        {"post_processing", {{"path", "results/"}, {"output_steps", 10}}}};

    // Make json object with entity sets files
    Json json_set_file = {
        {"node_sets",
         {{{"id", 0}, {"set", {0, 1}}}, {{"id", 1}, {"set", {2, 3}}}}},
        {"particle_sets",
         {{{"id", 0}, {"set", {0, 1}}}, {{"id", 1}, {"set", {2, 3}}}}}};
    // Dump JSON as an input file to be read
    std::ofstream file;
    file.open("mpm.json");
    file << json_file.dump(2);
    file.close();

    // Dump JSON as an entity sets file to be read
    std::ofstream set_file;
    set_file.open("entity_sets.json");
    set_file << json_set_file.dump(2);
    set_file.close();

    // Assign argc and argv to nput arguments of MPM
    int argc = 7;
    // clang-format off
    char* argv[] = {(char*)"./mpm",
                    (char*)"-f",  (char*)"./",
                    (char*)"-p",  (char*)"8",
                    (char*)"-i",  (char*)"mpm.json"};
    // clang-format on

    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);

    // Check number of threads
    REQUIRE(io->nthreads() == 8);

    // Check analysis type
    REQUIRE(io->analysis_type() == "MPMExplicit3D");

    // Check cmake JSON object
    REQUIRE(io->file_name("config") == "./mpm.json");

    // Material file should return an empty string, as file is missing
    std::string material_file = io->file_name("material");
    REQUIRE(material_file.empty() == true);

    // Check if mpm.json exists
    REQUIRE(io->check_file("./mpm.json") == true);

    // Check if entity_sets.json exists
    REQUIRE(io->check_file("./entity_sets.json") == true);

    // Check if a non-existant file is present
    REQUIRE(io->check_file("../fail.txt") == false);

    // Get analysis object
    Json analysis = io->analysis();
    // Check analysis dt
    REQUIRE(analysis["dt"] == json_file["analysis"]["dt"]);
    // Check analysis number_steps
    REQUIRE(analysis["nsteps"] == json_file["analysis"]["nsteps"]);
    // Check analysis gravity
    REQUIRE(analysis["gravity"] == json_file["analysis"]["gravity"]);
    // Check analysis soil_particle_spacing
    REQUIRE(analysis["soil_particle_spacing"] ==
            json_file["analysis"]["soil_particle_spacing"]);
    // Check analysis boundary_friction
    REQUIRE(analysis["boundary_friction"] ==
            json_file["analysis"]["boundary_friction"]);
    // Check analysis gravity
    REQUIRE(analysis["damping"] == json_file["analysis"]["damping"]);
    // Check analysis gravity
    REQUIRE(analysis["newmark"] == json_file["analysis"]["newmark"]);

    // Check return of a named JSON object
    const std::string obj_name = "mesh";
    REQUIRE(io->json_object(obj_name) == json_file["mesh"]);

    // Check return of materials JSON object
    const std::string materials = "materials";
    REQUIRE(io->json_object(materials) == json_file["materials"]);

    // Get post processing object
    Json post_processing = io->post_processing();
    // Check post processing data
    REQUIRE(post_processing == json_file["post_processing"]);

    // Check input file name
    const std::string attribute = "geometry";
    const std::string extension = ".vtp";
    const std::string uuid_ = "MPM";
    unsigned step = 57;
    unsigned max_steps = 100;
    auto meshfile =
        io->output_file(attribute, extension, uuid_, step, max_steps).string();
    REQUIRE(meshfile == "./results/MPM/geometry057.vtp");
    // Check output folder
    REQUIRE(io->output_folder() == "results/");

    // Check entity sets
    tsl::robin_map<mpm::Index, std::vector<mpm::Index>> node_sets;
    node_sets.insert(std::pair<mpm::Index, std::vector<mpm::Index>>(0, {0, 1}));
    node_sets.insert(std::pair<mpm::Index, std::vector<mpm::Index>>(1, {2, 3}));
    tsl::robin_map<mpm::Index, std::vector<mpm::Index>> check =
        io->entity_sets(io->file_name("entity_sets"), "node_sets");
    REQUIRE(std::equal(check.begin(), check.end(), node_sets.begin()) == true);
    REQUIRE(check.size() == node_sets.size());
  }
}
