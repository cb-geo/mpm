#include "io.h"

//! Constructor with argc and argv
mpm::IO::IO(int argc, char** argv) {
  //! Logger
  console_ = spdlog::get("IO");
  try {
    // Set title
    TCLAP::CmdLine cmd("Material Point Method (CB-Geo)", ' ', "Alpha V1.0");

    // Define working directory
    TCLAP::ValueArg<std::string> cwd_arg(
        "f", "working_dir", "Current working folder", true, "", "working_dir");
    cmd.add(cwd_arg);

    // Define input file
    TCLAP::ValueArg<std::string> input_arg("i", "input_file",
                                           "Input JSON file [mpm.json]", false,
                                           "mpm.json", "input_file");
    cmd.add(input_arg);

    // Define # TBB parallel threads
    TCLAP::ValueArg<unsigned int> tbb_arg("p", "tbb_parallel",
                                          "Number of parallel TBB threads",
                                          false, 0, "tbb_parallel");
    cmd.add(tbb_arg);

    // Parse arguments
    cmd.parse(argc, argv);

    // Set working directory
    working_dir_ = cwd_arg.getValue();

    // Set input file if the optional argument is not empty
    input_file_ = input_arg.getValue();

    // Set number of threads
    nthreads_ = tbb_arg.getValue();

  } catch (TCLAP::ArgException& except) {  // catch any exceptions
    console_->error("error: {}  for arg {}", except.error(), except.argId());
  }

  // Get input JSON file
  std::string file = working_dir_ + input_file_;
  std::ifstream ifs(file);

  try {
    if (!ifs.is_open())
      throw std::runtime_error(
          std::string("Input file not found in specified location: ") + file);
  } catch (const std::runtime_error& except) {
    console_->error("{}", except.what());
    std::terminate();
  }

  json_ = Json::parse(ifs);
}

//! Return input file name of mesh/submesh/soil particles
//! or an empty string if specified file for the key is not found
std::string mpm::IO::file_name(const std::string& filename) {

  std::string file_name;
  // Read input file name from the JSON object
  try {
    file_name = working_dir_ +
                json_["input_files"][filename].template get<std::string>();
  } catch (const std::exception& except) {
    console_->warn("Invalid JSON argument: {}; error: {}", filename,
                   except.what());
    file_name.clear();
    return file_name;
  }

  // Check if a file is present, if not set file_name to empty
  if (!this->check_file(file_name)) file_name.clear();

  return file_name;
}

//! Check if a file is present
bool mpm::IO::check_file(const std::string& filename) {
  bool status = false;

  // Check if file is present
  std::ifstream file;
  file.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  try {
    file.open(filename);
    status = true;
    file.close();
  } catch (std::exception& exception) {
    status = false;
    console_->error("Failed to find file {}: {}", filename, exception.what());
  }
  return status;
}

//! Create output VTK file names (eg. Velocity0000*.vtk)
boost::filesystem::path mpm::IO::output_file(const std::string& attribute,
                                             const std::string& file_extension,
                                             const std::string& analysis_id,
                                             unsigned step,
                                             unsigned max_steps) {
  std::stringstream file_name;
  std::string path = this->output_folder();

  file_name.str(std::string());
  file_name << attribute;
  file_name.fill('0');
  int digits = log10(max_steps) + 1;
  file_name.width(digits);
  file_name << step;

#ifdef USE_MPI
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  // Get number of MPI ranks
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  if (mpi_size > 1) {
    const std::string rank_size =
        "-" + std::to_string(mpi_rank) + "_" + std::to_string(mpi_size);
    file_name << rank_size;
  }
#endif

  file_name << file_extension;

  // Include path
  if (!path.empty()) path = working_dir_ + path;

  // Create results folder if not present
  boost::filesystem::path dir(path);
  if (!boost::filesystem::exists(dir)) boost::filesystem::create_directory(dir);

  // Create analysis folder
  path += analysis_id + "/";
  dir = path;
  if (!boost::filesystem::exists(dir)) boost::filesystem::create_directory(dir);

  boost::filesystem::path file_path(path + file_name.str().c_str());
  return file_path;
}

//! Return output folder
std::string mpm::IO::output_folder() const {
  std::string path{"results/"};

  Json json_postprocess = this->post_processing();

  try {
    auto results = json_postprocess.at("path");
    if (!results.empty()) path = results;

  } catch (std::exception& except) {
    console_->error("Output file creation: {}", except.what());
    console_->warn("Using default path: {}", path);
  }
  return path;
}
