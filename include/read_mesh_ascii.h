#ifndef MPM_READ_MESH_ASCII_H_
#define MPM_READ_MESH_ASCII_H_

#include <vector>

#include "Eigen/Dense"

#include "read_mesh.h"

//! MPM namespace
namespace mpm {

//! Global index type for the cell
using Index = unsigned long long;

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
    console_ = spdlog::get("ReadMeshAscii");
  }

  //! Destructor
  ~ReadMeshAscii() override = default;

  //! Read mesh nodes file
  //! \param[in] mesh file name with nodes and cells
  //! \retval coordinates Vector of nodal coordinates
  std::vector<VectorDim> read_mesh_nodes(const std::string& mesh) override;

  //! Read mesh cells file
  //! \param[in] mesh file name with nodes and cells
  //! \retval cells Vector of nodal indices of cells
  std::vector<std::vector<mpm::Index>> read_mesh_cells(
      const std::string& mesh) override;

  //! Read particles file
  //! \param[in] particles_files file name with particle coordinates
  //! \retval coordinates Vector of particle coordinates
  std::vector<VectorDim> read_particles(
      const std::string& particles_file) override;

  //! Read particle stresses
  //! \param[in] particles_stresses file name with particle stresses
  //! \retval stresses Vector of particle stresses
  std::vector<Eigen::Matrix<double, 6, 1>> read_particles_stresses(
      const std::string& particles_stresses) override;

  //! Read constraints file
  //! \param[in] velocity_constraints_files file name with constraints
  std::vector<std::tuple<mpm::Index, unsigned, double>>
      read_velocity_constraints(
          const std::string& velocity_constraints_file) override;

  //! Read volume file
  //! \param[in] volume_files file name with particle volumes
  std::vector<std::tuple<mpm::Index, double>> read_particles_volumes(
      const std::string& volume_file) override;

  //! Read traction file
  //! \param[in] traction_files file name with particle tractions
  std::vector<std::tuple<mpm::Index, unsigned, double>>
      read_particles_tractions(const std::string& traction_file) override;

  //! Read particles cells file
  //! \param[in] particles_cells_file file name with particle cell ids
  std::vector<std::array<mpm::Index, 2>> read_particles_cells(
      const std::string& particles_cells_file) override;

  //! Write particles cells file
  //! \param[in] particle_cells List of particles and cells
  //! \param[in] particles_cells_file file name with particle cell ids
  void write_particles_cells(
      const std::string& particles_cells_file,
      const std::vector<std::array<mpm::Index, 2>>& particles_cells) override;

 private:
  //! Logger
  std::shared_ptr<spdlog::logger> console_;
};  // ReadAscii class
}  // namespace mpm

#include "read_mesh_ascii.tcc"

#endif  // MPM_READ_MESH_ASCII_H_
