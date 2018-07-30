#include "vtk_writer.h"

//! VTK Writer class Constructor with coordniates
//! \param[in] coordinate Point coordinates
//! \param[in] node_pairs Node ID pairs to form elements
VtkWriter::VtkWriter(
    const std::vector<Eigen::Matrix<double, 3, 1>>& coordinates) {
  // Assign points
  points_ = vtkSmartPointer<vtkPoints>::New();
  unsigned id = 0;
  for (const auto& coordinate : coordinates) {
    const double* point = coordinate.data();
    points_->InsertPoint(id, point);
    ++id;
  }
}

//! Write coordinates
//! \param[in] filename Output file to write geometry
void VtkWriter::write_geometry(const std::string& filename) {

  // Create a polydata to store everything in it
  auto pdata = vtkSmartPointer<vtkPolyData>::New();

  // Add the points to the dataset
  pdata->SetPoints(points_);

  // Write file
  auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

  writer->SetFileName(filename.c_str());

#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(pdata);
#else
  writer->SetInputData(pdata);
#endif
  writer->Write();
}

//! \brief Write vector data
//! \param[in] filename Output file to write geometry
//! \param[in] data Vector data
//! \param[in] data_field Field name ("Displacement", "Forces")
void VtkWriter::write_vector_point_data(
    const std::string& filename, const std::vector<Eigen::Vector3d>& data,
    const std::string& data_field) {

  // Create a polydata to store everything in it
  auto pdata = vtkSmartPointer<vtkPolyData>::New();

  // Add the points to the dataset
  pdata->SetPoints(points_);

  // Create an array to hold distance information
  auto vectordata = vtkSmartPointer<vtkDoubleArray>::New();
  vectordata->SetNumberOfComponents(3);
  vectordata->SetNumberOfTuples(pdata->GetNumberOfPoints());
  vectordata->SetName(data_field.c_str());

  // Evaluate the signed distance function at all of the grid points
  for (vtkIdType id = 0; id < pdata->GetNumberOfPoints(); ++id) {
    const double* vdata = data.at(id).data();
    vectordata->SetTuple(id, vdata);
  }

  // Add the SignedDistances to the grid
  pdata->GetPointData()->SetVectors(vectordata);

  // Write file
  auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

  writer->SetFileName(filename.c_str());

#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(pdata);
#else
  writer->SetInputData(pdata);
#endif

  writer->SetCompressor(vtkZLibDataCompressor::New());
  writer->Write();
}
