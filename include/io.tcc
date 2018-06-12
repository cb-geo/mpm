#include "io.h"

//! Constructor
//! \param[in] json object input 
IO::IO(json json_filename) {

  /*
  // Set title
  TCLAP::CmdLine cmd("Material Point Generator (CB-Geo)", ' ', "0.0.1");

  // Define working directory
  TCLAP::ValueArg<std::string> cwd_arg(
      "f", "working_dir", "Current working folder", true, "", "Working_folder");
  cmd.add(cwd_arg);

  // Define input file
  TCLAP::ValueArg<std::string> input_arg("i", "input_file",
                                         "Input JSON file [cube.json]", false,
                                         "cube.json", "input_file");
  cmd.add(input_arg);

  // Parse arguments
  cmd.parse(argc, argv);

  //! Set working directory
  working_dir_ = cwd_arg.getValue();

  //! Set input file if the optional argument is not empty
  const auto input_file = input_arg.getValue();

  const auto json_filename = working_dir_ + input_file;
  */

  working_dir_ = "folder";

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

  //! Store json object directory to private variables
  //! Check if mesh file is present
  inputcheckfileName_ =
      working_dir_ + json_["inputcheckfileName"].template get<std::string>();
  std::ifstream checkfile;
  checkfile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  try {
    checkfile.open(inputcheckfileName_);
  } catch (const std::ifstream::failure& except) {
    std::cerr << "Exception opening/reading mesh file";
  }
  checkfile.close();

  //! Check if submesh file is present
  inputSubcheckfileName_ =
      working_dir_ + json_["inputSubcheckfileName"].template get<std::string>();
  std::ifstream checkfile;
  checkfile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  try {
    checkfile.open(inputSubcheckfileName_);
  } catch (const std::ifstream::failure& except) {
    std::cerr << "Exception opening/reading submesh file";
  }
  checkfile.close();

  //! Check if constraints file is present
  constraintsFileName_ =
      working_dir_ + json_["constraintsFileName"].template get<std::string>();
  std::ifstream checkfile;
  checkfile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  try {
    checkfile.open(constraintsFileName_);
  } catch (const std::ifstream::failure& except) {
    std::cerr << "Exception opening/reading constraints file";
  }
  checkfile.close();

  //! Check if input particle soil file is present
  inputSoilParticleFileName_ =
      working_dir_ + json_["inputSoilParticleFileName"].template get<std::string>();
  std::ifstream checkfile;
  checkfile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  try {
    checkfile.open(inputSoilParticleFileName_);
  } catch (const std::ifstream::failure& except) {
    std::cerr << "Exception opening/reading particle soil file";
  }
  checkfile.close();

  //! Check if input initial stress for soil particles file is present
  initStressSoilPFileName_ =
      working_dir_ + json_["initStressSoilPFileName"].template get<std::string>();
  std::ifstream checkfile;
  checkfile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  try {
    checkfile.open(initStressSoilPFileName_);
  } catch (const std::ifstream::failure& except) {
    std::cerr << "Exception opening/reading initial stress file";
  }
  checkfile.close();

  //! Check if material file is present
  materialFileName_ =
      working_dir_ + json_["materialFileName"].template get<std::string>();
  std::ifstream checkfile;
  checkfile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  try {
    checkfile.open(materialFileName_);
  } catch (const std::ifstream::failure& except) {
    std::cerr << "Exception opening/reading material file";
  }
  checkfile.close();

  //! Check if traction file is present
  tractionsSoilPFileName_ =
      working_dir_ + json_["tractionsSoilPFileName"].template get<std::string>();
  std::ifstream checkfile;
  checkfile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  try {
    checkfile.open(tractionsSoilPFileName_);
  } catch (const std::ifstream::failure& except) {
    std::cerr << "Exception opening/reading traction file";
  }
  checkfile.close();


  //! Store json object constants to private variables 

  //! Read and store gravity flag
  //! If not specified, set default value of 1
  try {
    if (json_.at("gravityFlag").size())
      gravityFlag_ = json_["gravityFlag"].template get<unsigned>();
  } catch (json::out_of_range& out_of_range) {
    std::cerr << out_of_range.what() << '\n';
    std::cout << "gravityFlag not specified. Using a default value of 1\n";
    gravityFlag_ = 1;
  }

  //! Read and store boundary friction miu
  //! If not specified, set default value of 0
  try {
    if (json_.at("boundaryFrictionMiu").size())
      boundaryFrictionMiu_ = json_["boundaryFrictionMiu"].template get<unsigned>();
  } catch (json::out_of_range& out_of_range) {
    std::cerr << out_of_range.what() << '\n';
    std::cout << "boundaryFrictionMiu not specified. Using a default value of 0\n";
    boundaryFrictionMiu_ = 0;
  }


  //! Read and store particle spacing
  //! If not specified, set default value of 1
  try {
    if (json_.at("soilParticleSpacing").size())
      soilParticleSpacing_ = json_["soilParticleSpacing"].template get<unsigned>();
  } catch (json::out_of_range& out_of_range) {
    std::cerr << out_of_range.what() << '\n';
    std::cout << "soilParticleSpacing not specified. Using a default value of 1\n";
    soilParticleSpacing_ = 1;
  }

  //! Read and store time interval dt
  //! If not specified, set default value of 0.001
  try {
    if (json_.at("timeInterval").size())
      dt_ = json_["timeInterval"].template get<unsigned>();
  } catch (json::out_of_range& out_of_range) {
    std::cerr << out_of_range.what() << '\n';
    std::cout << "timeInterval not specified. Using a default value of 0.001\n";
    dt_ = 0.001;
  }

  //! Read and store number of steps
  //! If not specified, set default value of 1000
  try {
    if (json_.at("numberOfSteps").size())
      numberOfSteps_ = json_["numberOfSteps"].template get<unsigned>();
  } catch (json::out_of_range& out_of_range) {
    std::cerr << out_of_range.what() << '\n';
    std::cout << "numberOfSteps not specified. Using a default value of 1000\n";
    numberOfSteps_ = 1000;
  } 

  //! Read and store number of sub steps
  //! If not specified, set default value of 1000
  try {
    if (json_.at("numberOfSubStepsOS").size())
      numberOfSubStepsOS_ = json_["numberOfSubStepsOS"].template get<unsigned>();
  } catch (json::out_of_range& out_of_range) {
    std::cerr << out_of_range.what() << '\n';
    std::cout << "numberOfSubStepsOS not specified. Using a default value of 1\n";
    numberOfSubStepsOS_ = 1;
  }  

  //! Read and store newmarkMethod flag
  //! If not specified, set default value of 0
  try {
    if (json_.at("newmarkMethod").size())
      newmarkMethod_ = json_["newmarkMethod"].template get<unsigned>();
  } catch (json::out_of_range& out_of_range) {
    std::cerr << out_of_range.what() << '\n';
    std::cout << "newmarkMethod not specified. Using a default value of 0\n";
    newmarkMethod_ = 0;
  }  

  //! Read and store gamma coefficient of newmark
  //! If not specified, set default value of 0.5
  try {
    if (json_.at("gamma").size())
      gamma_ = json_["gamma"].template get<unsigned>();
  } catch (json::out_of_range& out_of_range) {
    std::cerr << out_of_range.what() << '\n';
    std::cout << "gamma not specified. Using a default value of 0.5\n";
    gamma_ = 0.5;
  }  

  //! Read and store beta coefficient of newmark
  //! If not specified, set default value of 0.25
  try {
    if (json_.at("beta").size())
      beta_ = json_["beta"].template get<unsigned>();
  } catch (json::out_of_range& out_of_range) {
    std::cerr << out_of_range.what() << '\n';
    std::cout << "beta not specified. Using a default value of 0.25\n";
    beta_ = 0.25;
  }  

  //! Read and store damping flag
  //! If not specified, set default value of 0
  try {
    if (json_.at("dampingFlag").size())
      dampingFlag_ = json_["dampingFlag"].template get<unsigned>();
  } catch (json::out_of_range& out_of_range) {
    std::cerr << out_of_range.what() << '\n';
    std::cout << "dampingFlag not specified. Using a default value of 0\n";
    dampingFlag_ = 0;
  }  

  //! Read and store damping ratio
  //! If not specified, set default value of 0.05
  try {
    if (json_.at("dampingRatio").size())
      dampingRatio_ = json_["dampingRatio"].template get<unsigned>();
  } catch (json::out_of_range& out_of_range) {
    std::cerr << out_of_range.what() << '\n';
    std::cout << "dampingRatio not specified. Using a default value of 0.05\n";
    dampingRatio_ = 0.05;
  }  
}