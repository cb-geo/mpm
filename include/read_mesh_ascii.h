#ifndef MPM_READ_MESH_ASCII_H_
#define MPM_READ_MESH_ASCII_H_

#include <iostream>
#include <vector>

#include "Eigen/Dense"

#include "read_mesh.h"

//! MPM namespace
namespace mpm {

//! ReadMeshAscii class
//! \brief Derived class that returns mesh and particles locataions from ascii
//! file \tparam Tdim Dimension
template <unsigned Tdim>
class ReadMeshAscii : public ReadMesh<Tdim> {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Constructor
  ReadMeshAscii() : mpm::ReadMesh<Tdim>() {
    //! Logger
    console_ = spdlog::stdout_color_st("ReadMeshAscii");
  }

  //! Destructor
  ~ReadMeshAscii() = default;

  //! Read mesh nodes file
  //! \param[in] mesh_file file name with nodes and cells
  //! \retval coordinates Vector of nodal coordinates
  std::vector<VectorDim> read_mesh_nodes(const std::string& mesh_file);

  //! Read mesh cells file
  //! \param[in] mesh_file file name with nodes and cells
  //! \retval cells Vector of nodal indices of cells
  std::vector<std::vector<unsigned>> read_mesh_cells(
      const std::string& mesh_file);

 private:
  //! Logger
  std::shared_ptr<spdlog::logger> console_;
};  // ReadAscii class
}  // namespace mpm

#include "read_mesh_ascii.tcc"

#endif  // MPM_READ_MESH_ASCII_H_
