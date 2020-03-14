#include <boost/filesystem.hpp>

#include "catch.hpp"

#include "vtk_writer.h"

// Check Write VTK
TEST_CASE("VTK Writer is checked", "[vtk][writer]") {

  SECTION("Check scalar/vector/tensor file ") {
    std::vector<Eigen::Matrix<double, 3, 1>> coordinates;
    // VTK PolyData writer
    auto vtk_writer = std::make_unique<VtkWriter>(coordinates);

    SECTION("Write scalar data") {
      const std::string vtk_file = "scalar_vtk.vtp";
      const std::string attribute = "mass";
      std::vector<double> data;

      REQUIRE_NOTHROW(
          vtk_writer->write_scalar_point_data(vtk_file, data, attribute));
    }

    SECTION("Write vector data") {
      const std::string vtk_file = "vector_vtk.vtp";
      const std::string attribute = "displacements";
      std::vector<Eigen::Vector3d> data;

      REQUIRE_NOTHROW(
          vtk_writer->write_vector_point_data(vtk_file, data, attribute));
    }

    SECTION("Write tensor data") {
      const std::string vtk_file = "tensor_vtk.vtp";
      const std::string attribute = "stresses";
      std::vector<Eigen::Matrix<double, 6, 1>> data;

      REQUIRE_NOTHROW(
          vtk_writer->write_tensor_point_data(vtk_file, data, attribute));
    }
  }

  SECTION("Check parallel vtk vector file ") {
    std::vector<Eigen::Matrix<double, 3, 1>> coordinates;
    // VTK PolyData writer
    auto vtk_writer = std::make_unique<VtkWriter>(coordinates);

    const std::string parallel_vtk_file = "parallel_vector_vtk.pvtp";
    const std::string attribute = "displacements";
    int mpi_size = 2;
    unsigned step = 1000;
    unsigned max_steps = 10000;
    vtk_writer->write_parallel_vtk(parallel_vtk_file, attribute, mpi_size, step,
                                   max_steps);

    // Check file data
    std::string ppolydata =
        "<?xml version=\"1.0\"?>\n<VTKFile type=\"PPolyData\" version=\"0.1\" "
        "byte_order=\"LittleEndian\" "
        "compressor=\"vtkZLibDataCompressor\">\n<PPolyData "
        "GhostLevel=\"0\">\n\t<PPointData Vectors=\"" +
        attribute +
        "\">\n\t\t<PDataArray "
        "type=\"Float64\" Name=\"" +
        attribute +
        "\" "
        "NumberOfComponents=\"3\"/>\n\t</"
        "PPointData>\n\n\t<PPoints>\n\t\t<PDataArray "
        "type=\"Float32\" Name=\"Points\" "
        "NumberOfComponents=\"3\"/>\n\t</PPoints>\n";

    for (unsigned i = 0; i < mpi_size; ++i)
      ppolydata += "\n\t<Piece Source=\"displacements-" + std::to_string(i) +
                   "_" + std::to_string(mpi_size) + "-01000.vtp\"/>";

    ppolydata += "\n</PPolyData>\n\n</VTKFile>";

    // Check if file exists
    REQUIRE(boost::filesystem::exists(parallel_vtk_file) == true);

    std::ifstream ifs(parallel_vtk_file);
    std::string content((std::istreambuf_iterator<char>(ifs)),
                        (std::istreambuf_iterator<char>()));

    // Check file content
    REQUIRE(content.length() == ppolydata.length());
    REQUIRE(content == ppolydata);
  }

  SECTION("Check parallel vtk scalar file ") {
    std::vector<Eigen::Matrix<double, 3, 1>> coordinates;
    // VTK PolyData writer
    auto vtk_writer = std::make_unique<VtkWriter>(coordinates);

    const std::string parallel_vtk_file = "parallel_scalar_vtk.pvtp";
    const std::string attribute = "pdstrain";
    int mpi_size = 2;
    unsigned step = 1000;
    unsigned max_steps = 10000;
    vtk_writer->write_parallel_vtk(parallel_vtk_file, attribute, mpi_size, step,
                                   max_steps, 1);

    // Check file data
    std::string ppolydata =
        "<?xml version=\"1.0\"?>\n<VTKFile type=\"PPolyData\" version=\"0.1\" "
        "byte_order=\"LittleEndian\" "
        "compressor=\"vtkZLibDataCompressor\">\n<PPolyData "
        "GhostLevel=\"0\">\n\t<PPointData Scalars=\"" +
        attribute +
        "\">\n\t\t<PDataArray "
        "type=\"Float64\" Name=\"" +
        attribute +
        "\"/>\n\t</"
        "PPointData>\n\n\t<PPoints>\n\t\t<PDataArray "
        "type=\"Float32\" Name=\"Points\" "
        "NumberOfComponents=\"3\"/>\n\t</PPoints>\n";

    for (unsigned i = 0; i < mpi_size; ++i)
      ppolydata += "\n\t<Piece Source=\"" + attribute + "-" +
                   std::to_string(i) + "_" + std::to_string(mpi_size) +
                   "-01000.vtp\"/>";

    ppolydata += "\n</PPolyData>\n\n</VTKFile>";

    // Check if file exists
    REQUIRE(boost::filesystem::exists(parallel_vtk_file) == true);

    std::ifstream ifs(parallel_vtk_file);
    std::string content((std::istreambuf_iterator<char>(ifs)),
                        (std::istreambuf_iterator<char>()));

    // Check file content
    REQUIRE(content.length() == ppolydata.length());
    REQUIRE(content == ppolydata);
  }
}
