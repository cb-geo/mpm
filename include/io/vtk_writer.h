#ifndef VTK_WRITER_H_
#define VTK_WRITER_H_

#include <fstream>

#include <Eigen/Dense>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkLine.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkVersion.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkZLibDataCompressor.h>

#include <fstream>
#include <string>
#include <vector>

#include "data_types.h"

//! VTK Writer class
//! \brief VTK writer class
class VtkWriter {
 public:
  // Constructor with coordinates
  VtkWriter(const std::vector<Eigen::Matrix<double, 3, 1>>& coordinates);

  //! Write coordinates
  void write_geometry(const std::string& filename);

  //! Write tensor data
  //! \param[in] filename Output file to write geometry
  //! \param[in] data Tensor data
  //! \param[in] data_field Field name ("Displacement", "Forces")
  void write_tensor_point_data(
      const std::string& filename,
      const std::vector<Eigen::Matrix<double, 6, 1>>& data,
      const std::string& data_fields);

  //! Write vector data
  //! \param[in] filename Output file to write geometry
  //! \param[in] data Vector data
  //! \param[in] data_field Field name ("Displacement", "Forces")
  void write_vector_point_data(const std::string& filename,
                               const std::vector<Eigen::Vector3d>& data,
                               const std::string& data_fields);

  //! \brief Write scalar data
  //! \param[in] filename Output file to write geometry
  //! \param[in] data Scalar field data
  //! \param[in] data_field Field name ("pdstrain")
  void write_scalar_point_data(const std::string& filename,
                               const std::vector<double>& data,
                               const std::string& data_field);

  //! Write mesh
  //! \param[in] filename Mesh VTP file
  //! \param[in] coordinates Nodal coordinates
  //! \param[in] node_pairs Node pairs ids
  void write_mesh(const std::string& filename,
                  const std::vector<Eigen::Vector3d>& coordinates,
                  const std::vector<std::array<mpm::Index, 2>>& node_pairs);

  //! Write Parallel VTK file
  //! \param[in] filename Mesh PVTP file name
  //! \param[in] attribute VTK data attribute to be written
  //! \param[in] mpi_size Number of MPI tasks
  //! \param[in] step Current time step
  //! \param[in] max_steps Maximum number of steps in the simulation
  //! \param[in] ncomponents Number of components to write
  void write_parallel_vtk(const std::string& filename,
                          const std::string& attribute, int mpi_size,
                          unsigned step, unsigned max_steps,
                          unsigned ncomponents = 3);

 private:
  //! Vector of nodal coordinates
  vtkSmartPointer<vtkPoints> points_;
};

#endif  // VTK_WRITER_H_
