#include "io.h"

//! Constructor
IO::IO(int argc, char** argv) {
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

//! \brief Return user-specified mesh file name
std::string IO::file_name(const std::string& file) {

  std::string mesh_file_name;

  //! Check if file is present
  mesh_file_name =
      working_dir_ + json_["input_directory"].template get<std::string>() +
      json_["input_files"]["mesh_filename"].template get<std::string>();
  std::ifstream meshfile;
  meshfile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  try {
    meshfile.open(mesh_file_name);
  } catch (const std::ifstream::failure& except) {
    std::cerr << "Exception opening/reading mesh file";
  }
  meshfile.close();

  return mesh_file_name;
}
