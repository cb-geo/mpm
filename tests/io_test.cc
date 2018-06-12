#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <utility>

#include "Eigen/Dense"
#include "catch.hpp"

//! Alias for JSON
#include "json.hpp"
using json = nlohmann::json;

#include "io.h"

//! \brief Check IO for input string
TEST_CASE("IO for input string, boolean and double",
          "[IO][string][boolean][double]") {

  // Check add particle
  SECTION("Check for string input") {
    //! Make json object with input files
    json json_file =  {
                        {"inputMeshFileName", "/inputFiles/mesh.smf"},  
                        {"inputSubMeshFileName", "/inputFiles/submesh.dat"},
                        {"constraintsFileName", "/inputFiles/mesh.constraints"},  
                        {"inputSoilParticleFileName", "/inputFiles/soilParticles.dat"},  
                        {"initStressSoilPFileName", "/inputFiles/initStressSoilP.dat"},  
                        {"materialFileName", "/inputFiles/material.dat"},  
                        {"tractionsSoilPFileName", "/inputFiles/Traction.dat"},  
                        {"boundaryFrictionMiu", 0}, 
                        {"soilParticleSpacing", 0},
                        {"timeInterval", 0},
                        {"numberOfSteps", 0},
                        {"numberOfSubStepsOS", 0},
                        {"gamma", 0},
                        {"beta", 0},
                        {"dampingRatio", 0},
                        {"gravityFlag", 0},
                        {"newmarkMethod", 0},
                        {"dampingFlag", 0}
                      };
                      
    //! Make pointers to io
    auto io = std::unique_ptr<IO>(new IO(json_file));

    //! Check for only string inputs
    // REQUIRE(io->inputMeshFileName() == "folder/inputFiles/mesh.smf");
  }
}


