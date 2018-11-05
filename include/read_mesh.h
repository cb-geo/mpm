#ifndef MPM_READ_MESH_H_
#define MPM_READ_MESH_H_

#include <array>
#include <exception>
#include <fstream>
#include <memory>
#include <tuple>
#include <vector>

// Boost string algorithm
#include "Eigen/Dense"
#include <boost/algorithm/string.hpp>
// Speed log
#include "spdlog/spdlog.h"

#include "logger.h"

namespace mpm {

//! Global index type for the cell
using Index = unsigned long long;

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
  virtual std::vector<std::vector<mpm::Index>> read_mesh_cells(
      const std::string& mesh) = 0;

  //! Read particles file
  //! \param[in] particles_files file name with particle coordinates
  //! \retval coordinates Vector of particle coordinates
  virtual std::vector<VectorDim> read_particles(
      const std::string& particles_file) = 0;

  //! Read particle stresses
  //! \param[in] particles_stresses file name with particle stresses
  //! \retval stresses Vector of particle stresses
  virtual std::vector<Eigen::Matrix<double, 6, 1>> read_particles_stresses(
      const std::string& particles_stresses) = 0;

  //! Read velocity constraints file
  //! \param[in] velocity_constraints_files file name with constraints
  virtual std::vector<std::tuple<mpm::Index, unsigned, double>>
      read_velocity_constraints(
          const std::string& velocity_constraints_file) = 0;

  //! Read particles traction file
  //! \param[in] traction_files file name with particle tractions
  virtual std::vector<std::tuple<mpm::Index, unsigned, double>>
      read_particles_tractions(const std::string& traction_file) = 0;

  //! Read particles cells file
  //! \param[in] particles_cells_file file name with particle cell ids
  virtual std::vector<std::array<mpm::Index, 2>> read_particles_cells(
      const std::string& particles_cells_file) = 0;

  //! Write particles cells file
  //! \param[in] particle_cells List of particles and cells
  //! \param[in] particles_cells_file file name with particle cell ids
  virtual void write_particles_cells(
      const std::string& particles_cells_file,
      const std::vector<std::array<mpm::Index, 2>>& particles_cells) = 0;

};  // ReadMesh class
}  // namespace mpm

#endif  // MPM_READ_MESH_H_
