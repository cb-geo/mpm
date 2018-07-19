#ifndef MPM_READ_MESH_H_
#define MPM_READ_MESH_H_

#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

#include "Eigen/Dense"
// Speed log
#include "spdlog/spdlog.h"

#include "logger.h"

namespace mpm {

//! ReadMesh class
//! \brief Abstract Base class that returns mesh and particles locataions
//! \tparam Tdim Dimension
template <unsigned Tdim>
class ReadMesh {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  // Constructor
  ReadMesh() = default;

  //! Default destructor
  virtual ~ReadMesh() = default;

  //! Delete copy constructor
  ReadMesh(const ReadMesh<Tdim>&) = delete;

  //! Delete assignement operator
  ReadMesh& operator=(const ReadMesh<Tdim>&) = delete;

  //! Read mesh nodes file
  //! \param[in] mesh file name with nodes and cells
  //! \retval coordinates Vector of nodal coordinates
  virtual std::vector<VectorDim> read_mesh_nodes(const std::string& mesh) = 0;

  //! Read mesh cells file
  //! \param[in] mesh file name with nodes and cells
  //! \retval cells Vector of nodal indices of cells
  virtual std::vector<std::vector<unsigned>> read_mesh_cells(
      const std::string& mesh) = 0;

  //! Read particles file
  //! \param[in] particles_files file name with particle coordinates
  //! \retval coordinates Vector of particle coordinates
  virtual std::vector<VectorDim> read_particles(
      const std::string& particles_file) = 0;
};  // ReadMesh class
}  // namespace mpm

#endif  // MPM_READ_MESH_H_
