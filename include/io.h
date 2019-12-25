#ifndef MPM_IO_H_
#define MPM_IO_H_

#include <fstream>
#include <memory>
#include <string>

#include <boost/filesystem.hpp>

#include "tclap/CmdLine.h"

#include "tsl/robin_map.h"
//! Alias for JSON
#include "json.hpp"
using Json = nlohmann::json;
// Speed log
#include "spdlog/spdlog.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

#include "data_types.h"

namespace mpm {
//! \brief Input/Output handler

class IO {
 public:
  //! Constructor with argc and argv
  //! \param[in] argc Number of input arguments
  //! \param[in] argv Input arguments
  IO(int argc, char** argv);

  //! Return number of tbb threads
  unsigned nthreads() const { return nthreads_; }

  //! Return working directory
  //! \retval string of working directory
  std::string working_directory() { return working_dir_; }

  //! Return input file name of mesh/submesh/soil particles
  //! or an empty string if specified file for the key is not found
  //! \param[in] key Input key in JSON for the filename of
  //! mesh/submesh/soilparticles
  std::string file_name(const std::string& key);

  //! Check if a file is present and valid
  //! \param[in] file_name Name of the file to check if it is present
  bool check_file(const std::string& file_name);

  //! Return analysis
  std::string analysis_type() const {
    return json_["analysis"]["type"].template get<std::string>();
  }

  //! Return json analysis object
  Json analysis() const { return json_["analysis"]; }

  //! Return json object
  Json json_object(const std::string& name) const { return json_[name]; }

  //! Return the entity sets from the input set JSON file
  //! \param[in] filename File name
  //! \param[in] sets_type type of sets
  //! \retval entity_sets map of entity sets
  tsl::robin_map<mpm::Index, std::vector<mpm::Index>> entity_sets(
      const std::string& filename, const std::string& sets_type);

  //! Return post processing object
  Json post_processing() const { return json_["post_processing"]; }

  //! Return the output folder for the analysis
  std::string output_folder() const;

  //! Create output VTK file names (eg. velocity0000*.vtk)
  //! Generates a file based on attribute, current step and maxsteps
  //! \param[in] attribute Attribute being written (eg., velocity / stress)
  //! \param[in] file_extension File Extension (*.vtk or *.vtp)
  //! \param[in] step Current step
  //! \param[in] max_steps Total number of steps to be solved
  //! \param[in] parallel Write output as parallel file system
  //! \return file_name File name with the correct attribute and a VTK extension
  boost::filesystem::path output_file(const std::string& attribute,
                                      const std::string& file_extension,
                                      const std::string& analysis_id,
                                      unsigned step, unsigned max_steps,
                                      bool parallel = true);

 private:
  //! Number of parallel threads
  unsigned nthreads_{0};
  //! Working directory
  std::string working_dir_;
  //! Input file name
  std::string input_file_{"mpm.json"};
  //! Input JSON object
  Json json_;
  //! Logger
  std::shared_ptr<spdlog::logger> console_;
};

}  // namespace mpm
#endif  // MPM_IO_H_
