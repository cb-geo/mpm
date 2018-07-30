#ifndef VTK_WRITER_H_
#define VTK_WRITER_H_

#include <fstream>

#include <Eigen/Dense>

#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkVersion.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkZLibDataCompressor.h>

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

 private:
  //! Vector of nodal coordinates
  vtkSmartPointer<vtkPoints> points_;
};

#endif  // VTK_WRITER_H_
