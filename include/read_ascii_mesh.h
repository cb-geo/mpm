#ifndef MPM_READ_MESH_ASCII_H_
#define MPM_READ_MESH_ASCII_H_

#include "read_mesh.h"

namespace mpm {

//! ReadMeshAscii class
//! \brief Derived class that returns mesh and particles locataions from ascii
//! file \tparam Tdim Dimension
template <unsigned Tdim>
class ReadMeshAscii : public ReadMeshMesh {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Read mesh nodes file
  //! \param[in] read_ascii file name with nodes and cells
  //! \retval nodes Vector of nodal coordinates
  std::vector<VectorDim> read_ascii_nodes(const std::string& mesh_file);

  //! Read mesh cells file
  //! \param[in] read_ascii file name with nodes and cells
  //! \retval cells Vector of nodal indices of cells
  std::vector<std::vector<unsigned>> read_ascii_cells(
      const std::string& mesh_file);

};  // ReadMeshAscii class
}  // namespace mpm

#include "read_mesh_ascii.tcc"

#endif  // MPM_READ_MESH_ASCII_H_
