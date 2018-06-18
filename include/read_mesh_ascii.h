#ifndef MPM_READ_ASCII_H_
#define MPM_READ_ASCII_H_

#include <iostream>
#include <vector>

#include "Eigen/Dense"

#include "read_mesh.h"

namespace mpm {

//! ReadAscii class
//! \brief Derived class that returns mesh and particles locataions from ascii
//! file \tparam Tdim Dimension
template <unsigned Tdim>
class ReadAscii : public ReadMesh {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Constructor
  ReadMeshAscii() {
    //! Logger
    console_ = spdlog::stdout_color_st("ReadMeshAscii");
  }

  //! Destructor
  ~Mesh() = default;

  //! Read mesh nodes file
  //! \param[in] read_ascii file name with nodes and cells
  //! \retval coordinates Vector of nodal coordinates
  std::vector<VectorDim> read_ascii_nodes(const std::string& mesh_file) = 0;

  //! Read mesh cells file
  //! \param[in] read_ascii file name with nodes and cells
  //! \retval cells Vector of nodal indices of cells
  std::vector<std::vector<unsigned>> read_ascii_cells(
      const std::string& mesh_file) = 0;

 private:
  //! Logger
  std::shared_ptr<spdlog::logger> console_;
};  // ReadAscii class
}  // namespace mpm

#endif  // MPM_READ_ASCII_H_
