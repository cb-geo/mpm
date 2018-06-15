#include "io.h"

//! Constructor with argc and argv
IO::IO(int argc, char** argv) {
  //! Logger
  console_ = spdlog::stdout_color_mt("IO");
  try {

    // Set title
    TCLAP::CmdLine cmd("Material Point Method (CB-Geo)", ' ', "Alpha");

    // Define working directory
    TCLAP::ValueArg<std::string> cwd_arg("f", "working_dir",
                                         "Current working folder", true, "",
                                         "Working_folder");
    cmd.add(cwd_arg);

    // Define input file
    TCLAP::ValueArg<std::string> input_arg("i", "input_file",
                                           "Input JSON file [mpm.json]", false,
                                           "mpm.json", "input_file");
    cmd.add(input_arg);

    // Parse arguments
    cmd.parse(argc, argv);

    // Set working directory
    working_dir_ = cwd_arg.getValue();

    // Set input file if the optional argument is not empty
    input_file_ = input_arg.getValue();

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
std::string IO::file_name(const std::string& filename) {

  std::string file_name;
  // Read input file name from the JSON object
  try {
    file_name = working_dir_ +
                json_["input_files"][filename].template get<std::string>();
  } catch (const std::exception& except) {
    console_->warn("Invalid JSON argument: {}", except.what());
  }

  // Check if a file is present, if not set file_name to empty
  if (!this->check_file(file_name)) file_name.clear();

  return file_name;
}

//! Check if a file is present
bool IO::check_file(const std::string& filename) {
  bool status = false;

  // Check if file is present
  std::ifstream file;
  file.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  try {
    file.open(filename);
    status = true;
    file.close();
  } catch (std::ifstream::failure& exception) {
    status = false;
    console_->error("Failed to find file: {}", exception.what());
  }
  return status;
}
