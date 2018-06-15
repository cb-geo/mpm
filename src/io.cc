#include "io.h"

//! Constructor
//! \param[in] argc Number of input arguments
//! \param[in] argv Array of input arguments
IO::IO(int argc, char** argv) {

  // Set title
  TCLAP::CmdLine cmd("Material Point Method (CB-Geo)", ' ', "0.0.1");

  // Define working directory
  TCLAP::ValueArg<std::string> cwd_arg(
      "f", "working_dir", "Current working folder", true, "", "Working_folder");
  cmd.add(cwd_arg);

  // Define input file
  TCLAP::ValueArg<std::string> input_arg("i", "input_file",
                                         "Input JSON file [test.json]", false,
                                         "test.json", "input_file");
  cmd.add(input_arg);

  // Parse arguments
  cmd.parse(argc, argv);

  //! Set working directory
  working_dir_ = cwd_arg.getValue();

  //! Set input file if the optional argument is not empty
  const auto input_file = input_arg.getValue();

  const auto json_filename = working_dir_ + input_file;

  //! Check if json file is present
  std::ifstream inputFile(json_filename);

  try {
    if (!inputFile.is_open())
      throw std::runtime_error(
          std::string("Input file not found in specified location: ") +
          json_filename);
  } catch (const std::runtime_error& except) {
    std::cerr << "Exception opening/reading json file";
  }

  //! Store json object as private variable
  //! Read json file and store to private variables
  json_ = json::parse(inputFile);
}

//! \brief Return user-specified mesh file name
std::string IO::mesh_file_name() const {

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

//! \brief Return user-specified submesh file name
std::string IO::submesh_file_name() const {

  std::string submesh_file_name;

  //! Check if file is present
  submesh_file_name =
      working_dir_ + json_["input_directory"].template get<std::string>() +
      json_["input_files"]["submesh_filename"].template get<std::string>();
  std::ifstream meshfile;
  meshfile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  try {
    meshfile.open(submesh_file_name);
  } catch (const std::ifstream::failure& except) {
    std::cerr << "Exception opening/reading submesh file";
  }
  meshfile.close();

  return submesh_file_name;
}

//! \brief Return user-specified constraints file name
std::string IO::constraints_file_name() const {

  std::string constraints_file_name;

  //! Check if file is present
  constraints_file_name =
      working_dir_ + json_["input_directory"].template get<std::string>() +
      json_["input_files"]["constraints_filename"].template get<std::string>();
  std::ifstream meshfile;
  meshfile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  try {
    meshfile.open(constraints_file_name);
  } catch (const std::ifstream::failure& except) {
    std::cerr << "Exception opening/reading constraints file";
  }
  meshfile.close();

  return constraints_file_name;
}

//! \brief Return user-specified soil particle file name
std::string IO::soil_particle_file_name() const {

  std::string soil_particle_file_name;

  //! Check if file is present
  soil_particle_file_name =
      working_dir_ + json_["input_directory"].template get<std::string>() +
      json_["input_files"]["soil_particle_filename"]
          .template get<std::string>();
  std::ifstream meshfile;
  meshfile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  try {
    meshfile.open(soil_particle_file_name);
  } catch (const std::ifstream::failure& except) {
    std::cerr << "Exception opening/reading soil particle file";
  }
  meshfile.close();

  return soil_particle_file_name;
}

//! \brief Return user-specified initial stress of soil particle file name
std::string IO::initial_stress_file_name() const {

  std::string initial_stress_file_name;

  //! Check if file is present
  initial_stress_file_name =
      working_dir_ + json_["input_directory"].template get<std::string>() +
      json_["input_files"]["initial_stress_filename"]
          .template get<std::string>();
  std::ifstream meshfile;
  meshfile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  try {
    meshfile.open(initial_stress_file_name);
  } catch (const std::ifstream::failure& except) {
    std::cerr << "Exception opening/reading initial stress file";
  }
  meshfile.close();

  return initial_stress_file_name;
}

//! \brief Return user-specified material file name
std::string IO::material_file_name() const {

  std::string material_file_name;

  //! Check if file is present
  material_file_name =
      working_dir_ + json_["input_directory"].template get<std::string>() +
      json_["input_files"]["material_filename"].template get<std::string>();
  std::ifstream meshfile;
  meshfile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  try {
    meshfile.open(material_file_name);
  } catch (const std::ifstream::failure& except) {
    std::cerr << "Exception opening/reading material file";
  }
  meshfile.close();

  return material_file_name;
}

//! \brief Return user-specified traction file name
std::string IO::traction_file_name() const {

  std::string traction_file_name;

  //! Check if file is present
  traction_file_name =
      working_dir_ + json_["input_directory"].template get<std::string>() +
      json_["input_files"]["traction_filename"].template get<std::string>();
  std::ifstream meshfile;
  meshfile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  try {
    meshfile.open(traction_file_name);
  } catch (const std::ifstream::failure& except) {
    std::cerr << "Exception opening/reading traction file";
  }
  meshfile.close();

  return traction_file_name;
}

//! \brief Write output file names and store them in private member
//! \param[in] attribute Attribute being written (eg., material_points / stress)
//! \param[in] file_extension File Extension (*.txt)
boost::filesystem::path IO::output_file(const std::string& attribute,
                                        const std::string& file_extension) {

  std::stringstream file_name;
  std::string path{"results/"};

  try {
    std::string results = json_["post_processing"]["output_directory"]
                              .template get<std::string>();
    if (!results.empty()) path = results;

  } catch (std::exception& except) {
    std::cerr << except.what() << '\n';
    // console_->error("Output file creation: {}", except.what());
    // console_->warn("Using default path: {}", path);
  }

  //! Make the file_name
  file_name.str(std::string());
  file_name << attribute << file_extension;

  //! Include path
  if (!path.empty()) path = working_dir_ + path;

  //! Create results folder if not present
  boost::filesystem::path dir(path);
  if (!boost::filesystem::exists(dir)) boost::filesystem::create_directory(dir);

  //! Create full path with working directory path and file name
  boost::filesystem::path file_path(path + file_name.str().c_str());
  return file_path;
}

//! \brief Return user-specified gravity flag
bool IO::gravity_flag() const {

  bool gravity_flag;

  //! If not specified, set default value of 1
  try {
    if (json_.at("analysis").at("gravity_glag").size())
      gravity_flag = json_["analysis"]["gravity_flag"].template get<bool>();
  } catch (json::out_of_range& out_of_range) {
    std::cerr << out_of_range.what() << '\n';
    std::cout << "gravity_flag not specified. Using a default value of 1\n";
    gravity_flag = 1;
  }

  return gravity_flag;
}

//! \brief Return user-specified boundary friction miu
double IO::boundary_friction() const {

  double boundary_friction;

  //! If not specified, set default value of 0
  try {
    if (json_.at("analysis").at("boundary_friction").size())
      boundary_friction =
          json_["analysis"]["boundary_friction"].template get<double>();
  } catch (json::out_of_range& out_of_range) {
    std::cerr << out_of_range.what() << '\n';
    std::cout
        << "boundary_friction not specified. Using a default value of 0\n";
    boundary_friction = 0;
  }

  return boundary_friction;
}

//! \brief Return user-specified soil particle spacing
double IO::soil_particle_spacing() const {

  double soil_particle_spacing;

  //! If not specified, set default value of 1
  try {
    if (json_.at("analysis").at("soil_particle_spacing").size())
      soil_particle_spacing =
          json_["analysis"]["soil_particle_spacing"].template get<double>();
  } catch (json::out_of_range& out_of_range) {
    std::cerr << out_of_range.what() << '\n';
    std::cout
        << "soil_particle_spacing not specified. Using a default value of 1\n";
    soil_particle_spacing = 1;
  }

  return soil_particle_spacing;
}

//! \brief Return user-specified time interval dt
double IO::dt() const {

  double dt;

  //! If not specified, set default value of 0.001
  try {
    if (json_.at("analysis").at("dt").size())
      dt = json_["analysis"]["dt"].template get<double>();
  } catch (json::out_of_range& out_of_range) {
    std::cerr << out_of_range.what() << '\n';
    std::cout << "dt not specified. Using a default value of 0.001\n";
    dt = 0.001;
  }

  return dt;
}

//! \brief Return user-specified number of steps
unsigned IO::number_steps() const {

  unsigned number_steps;

  //! If not specified, set default value of 1000
  try {
    if (json_.at("analysis").at("number_steps").size())
      number_steps = json_["analysis"]["number_steps"].template get<unsigned>();
  } catch (json::out_of_range& out_of_range) {
    std::cerr << out_of_range.what() << '\n';
    std::cout << "number_steps not specified. Using a default value of 1000\n";
    number_steps = 1000;
  }

  return number_steps;
}

//! \brief Return user-specified number of sub steps for post-processing
unsigned IO::number_substeps() const {

  unsigned number_substeps;

  //! If not specified, set default value of 10
  try {
    if (json_.at("post_processing").at("output_steps").size())
      number_substeps =
          json_["post_processing"]["output_steps"].template get<unsigned>();
  } catch (json::out_of_range& out_of_range) {
    std::cerr << out_of_range.what() << '\n';
    std::cout
        << "number_substeps not specified. Using a default value of 1000\n";
    number_substeps = 10;
  }

  return number_substeps;
}

//! \brief Return user-specified damping flag
bool IO::damping_flag() const {

  bool damping_flag;

  //! If not specified, set default value of 0
  try {
    if (json_.at("analysis").at("damping").at("damping_flag").size())
      damping_flag =
          json_["analysis"]["damping"]["damping_flag"].template get<bool>();
  } catch (json::out_of_range& out_of_range) {
    std::cerr << out_of_range.what() << '\n';
    std::cout << "damping_flag not specified. Using a default value of 0\n";
    damping_flag = 0;
  }

  return damping_flag;
}

//! \brief Return user-specified damping ratio
double IO::damping_ratio() const {

  double damping_ratio;

  //! If not specified, set default value of 0.05
  try {
    if (json_.at("analysis").at("damping").at("damping_ratio").size())
      damping_ratio =
          json_["analysis"]["damping"]["damping_ratio"].template get<double>();
  } catch (json::out_of_range& out_of_range) {
    std::cerr << out_of_range.what() << '\n';
    std::cout << "damping_ratio not specified. Using a default value of 0.05\n";
    damping_ratio = 0.05;
  }

  return damping_ratio;
}

//! \brief Return user-specified newmark method flag
bool IO::newmark_flag() const {

  bool newmark_flag;

  //! If not specified, set default value of 0
  try {
    if (json_.at("analysis").at("newmark").at("newmark_flag").size())
      newmark_flag =
          json_["analysis"]["newmark"]["newmark_flag"].template get<bool>();
  } catch (json::out_of_range& out_of_range) {
    std::cerr << out_of_range.what() << '\n';
    std::cout << "newmark_flag not specified. Using a default value of 0\n";
    newmark_flag = 0;
  }

  return newmark_flag;
}

//! \brief Return user-specified newmark method gamma coefficient
double IO::newmark_gamma() const {

  double newmark_gamma;

  //! If not specified, set default value of 0.5
  try {
    if (json_.at("analysis").at("newmark").at("newmark_gamma").size())
      newmark_gamma =
          json_["analysis"]["newmark"]["newmark_gamma"].template get<double>();
  } catch (json::out_of_range& out_of_range) {
    std::cerr << out_of_range.what() << '\n';
    std::cout << "newmark_gamma not specified. Using a default value of 0.5\n";
    newmark_gamma = 0.5;
  }

  return newmark_gamma;
}

//! \brief Return user-specified newmark method beta coefficient
double IO::newmark_beta() const {

  double newmark_beta;

  //! If not specified, set default value of 0.25
  try {
    if (json_.at("analysis").at("newmark").at("newmark_beta").size())
      newmark_beta =
          json_["analysis"]["newmark"]["newmark_beta"].template get<double>();
  } catch (json::out_of_range& out_of_range) {
    std::cerr << out_of_range.what() << '\n';
    std::cout << "newmark_beta not specified. Using a default value of 0.25\n";
    newmark_beta = 0.25;
  }

  return newmark_beta;
}