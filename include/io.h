#ifndef MPM_IO_H_
#define MPM_IO_H_

#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <utility>

#include <boost/filesystem.hpp>

#include "tclap/CmdLine.h"
//! Alias for JSON
#include "json.hpp"
using json = nlohmann::json;

//! \brief Input/Output handler
class IO {
 public:
  //! Constructor
  IO(int argc, char** argv);

  //! Return input mesh file name
  std::string mesh_file_name() const;

  //! Return input sub mesh file name
  std::string submesh_file_name() const;

  //! Return input constraint file name
  std::string constraints_file_name() const;

  //! Return input soil particle file name
  std::string soil_particle_file_name() const;

  //! Return initial stress soil particle file name
  std::string initial_stress_file_name() const;

  //! Return material file name
  std::string material_file_name() const;

  //! Return traction soil particle file name
  std::string traction_file_name() const;

  //! Create output file names
  boost::filesystem::path output_file(const std::string& attribute,
                                      const std::string& file_extension);

  //! Return gravity flag
  bool gravity_flag() const;

  //! Return boundary friction miu
  double boundary_friction() const;

  //! Return soil particle spacing
  double soil_particle_spacing() const;

  //! Return time step or interval
  double dt() const;

  //! Return number of steps of the problem
  unsigned number_steps() const;

  //! Return number of substeps of the problem
  unsigned number_substeps() const;

  //! Return Cundall Damping flag
  bool damping_flag() const;

  //! Return damping ratio for Cundall Damping
  double damping_ratio() const;

  //! Return newmark flag
  bool newmark_flag() const;

  //! Return gamma - newmark integration coefficient
  double newmark_gamma() const;

  //! Return beta - newmark integration coefficient
  double newmark_beta() const;

 private:
  //! Input directory
  std::string working_dir_;

  //! Input json object
  json json_;
};

#endif  // MPM_IO_H_
