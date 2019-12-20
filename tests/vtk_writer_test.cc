#include <boost/filesystem.hpp>

#include "catch.hpp"

#include "vtk_writer.h"

// Check ReadMeshAscii
TEST_CASE("VTK Writer is checked", "[vtk][writer]") {

  SECTION("Check parallel vtk file ") {
    std::vector<Eigen::Matrix<double, 3, 1>> coordinates;
    // VTK PolyData writer
    auto vtk_writer = std::make_unique<VtkWriter>(coordinates);

    const std::string parallel_vtk_file = "parallel_vtk.pvtp";
    const std::string attribute = "stress";
    int mpi_size = 2;
    unsigned step = 1000;
    unsigned max_steps = 10000;
    vtk_writer->write_parallel_vtk(parallel_vtk_file, attribute, mpi_size, step,
                                   max_steps);

    REQUIRE(boost::filesystem::exists(parallel_vtk_file) == true);
  }
}
