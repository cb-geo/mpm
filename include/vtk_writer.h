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

#include <string>
#include <vector>

namespace mpm {

//! Global index type for the cell
using Index = unsigned long long;

}  // namespace mpm

//! VTK Writer class
//! \brief VTK writer class
class VtkWriter {
 public:
  // Constructor with coordinates
  VtkWriter(const std::vector<Eigen::Matrix<double, 3, 1>>& coordinates);

  //! Write coordinates
  void write_geometry(const std::string& filename);

  //! Write vector data
  void write_vector_point_data(const std::string& filename,
                               const std::vector<Eigen::Vector3d>& data,
                               const std::string& data_fields);

  //! Write mesh
  //! \param[in] filename Mesh VTP file
  //! \param[in] coordinates Nodal coordinates
  //! \param[in] node_pairs Node pairs ids
  void write_mesh(const std::string& filename,
                  const std::vector<Eigen::Vector3d>& coordinates,
                  const std::vector<std::array<mpm::Index, 2>>& node_pairs);

 private:
  //! Vector of nodal coordinates
  vtkSmartPointer<vtkPoints> points_;
};

#endif  // VTK_WRITER_H_
