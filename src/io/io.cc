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

    // Define # parallel threads
    TCLAP::ValueArg<unsigned int> parallel_arg(
        "p", "parallel", "Number of parallel threads", false, 0, "parallel");
    cmd.add(parallel_arg);

    // Parse arguments
    cmd.parse(argc, argv);

    // Set working directory
    working_dir_ = cwd_arg.getValue();

    // Set input file if the optional argument is not empty
    input_file_ = input_arg.getValue();

    // Set number of threads
    nthreads_ = parallel_arg.getValue();

  } catch (TCLAP::ArgException& except) {  // catch any exceptions
    console_->error("error: {}  for arg {}", except.error(), except.argId());
  }

  // Get input JSON file
  std::string file = working_dir_ + input_file_;
  std::ifstream ifs(file);

  if (!ifs.is_open())
    throw std::runtime_error(
        std::string("Input file not found in specified location: ") + file);

  json_ = Json::parse(ifs);
}

//! Return input file name of mesh/submesh/soil particles
//! or an empty string if specified file for the key is not found
std::string mpm::IO::file_name(const std::string& filename) {

  std::string file_name;
  // Read input file name from the JSON object
  try {
    file_name = working_dir_ + filename;
    // Check if a file is present, if not set file_name to empty
    if (!this->check_file(file_name))
      throw std::runtime_error("no file found!");

  } catch (const std::exception& except) {
    console_->warn("Fetching file: {}; failed with: {}", filename,
                   except.what());
    file_name.clear();
    return file_name;
  }
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
                                             unsigned step, unsigned max_steps,
                                             bool parallel) {
  std::stringstream file_name;
  std::string path = this->output_folder();

  file_name.str(std::string());
  file_name << attribute;

#ifdef USE_MPI
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  // Get number of MPI ranks
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  if (mpi_size > 1 && parallel) {
    const std::string rank_size =
        "-" + std::to_string(mpi_rank) + "_" + std::to_string(mpi_size) + "-";
    file_name << rank_size;
  }
#endif

  file_name.fill('0');
  int digits = log10(max_steps) + 1;
  file_name.width(digits);
  file_name << step;
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

//! Return map of entity sets from the JSON file
tsl::robin_map<mpm::Index, std::vector<mpm::Index>> mpm::IO::entity_sets(
    const std::string& filename, const std::string& sets_type) {

  // Input file stream for sets JSON file
  std::ifstream sets_file;
  sets_file.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  // Map of the entity sets
  tsl::robin_map<mpm::Index, std::vector<mpm::Index>> entity_sets;

  try {
    sets_file.open(filename);
    // Entity sets JSON object
    Json json_sets_ = Json::parse(sets_file);
    // type of sets (e.g: node_sets or particle_sets)
    Json sets = json_sets_[sets_type];

    if (sets.size() > 0) {
      for (Json::iterator itr = sets.begin(); itr != sets.end(); ++itr) {
        // Get the entity set ids
        mpm::Index id = (*itr)["id"].template get<mpm::Index>();
        // Get the entity ids
        std::vector<mpm::Index> entity_ids = (*itr).at("set");
        // Add the entity set to the list
        entity_sets.insert(
            std::pair<mpm::Index, std::vector<mpm::Index>>(id, entity_ids));
      }
    }
    sets_file.close();
  } catch (const std::out_of_range& range_error) {
    console_->warn("{} {} reading {}: {}", __FILE__, __LINE__, sets_type,
                   filename, range_error.what());
  } catch (const std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }

  return entity_sets;
}

//! Return analysis
std::string mpm::IO::analysis_type() const {
  return json_["analysis"]["type"].template get<std::string>();
}

//! Return json analysis object
Json mpm::IO::analysis() const { return json_["analysis"]; }

//! Return json object
Json mpm::IO::json_object(const std::string& key) const {
  Json empty;
  if (json_.find(key) != json_.end()) {
    return json_.at(key);
  } else {
    throw std::runtime_error("No object found, returning an empty object");
    return empty;
  }
}

//! Return post processing object
Json mpm::IO::post_processing() const { return json_["post_processing"]; }

//! Return number of threads
unsigned mpm::IO::nthreads() const { return nthreads_; }
