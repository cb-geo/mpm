//! Return coordinates of nodes in a mesh from input file
template <unsigned Tdim>
std::vector<Eigen::Matrix<double, Tdim, 1>>
    mpm::ReadMeshAscii<Tdim>::read_mesh_nodes(const std::string& mesh) {
  // Nodal coordinates
  std::vector<VectorDim> coordinates;
  coordinates.clear();

  // input file stream
  std::fstream file;
  file.open(mesh.c_str(), std::ios::in);

  try {
    if (file.is_open() && file.good()) {
      // Line
      std::string line;
      // bool to check firstline
      bool read_first_line = false;
      // Number of coordinate lines
      unsigned nlines = 0;
      // Coordinates
      Eigen::Matrix<double, Tdim, 1> coords;
      // # of nodes and cells
      unsigned nnodes = 0, ncells = 0;
      // ignore stream
      double ignore;

      while (std::getline(file, line)) {
        boost::algorithm::trim(line);
        std::istringstream istream(line);
        // ignore comment lines (# or !) or blank lines
        if ((line.find('#') == std::string::npos) &&
            (line.find('!') == std::string::npos) && (line != "")) {
          while (istream.good()) {
            if (!read_first_line) {
              // Read number of nodes and cells
              istream >> nnodes >> ncells;
              coordinates.reserve(nnodes);
              read_first_line = true;
              break;
            }
            // Read until nodal information is present
            if (nlines <= nnodes) {
              // Read to coordinates
              for (unsigned i = 0; i < Tdim; ++i) istream >> coords[i];
              coordinates.emplace_back(coords);
              break;

            } else {
              // Ignore stream
              istream >> ignore;
            }
          }
          ++nlines;
        }
      }
    }
    file.close();
  } catch (std::exception& exception) {
    console_->error("Read mesh nodes: {}", exception.what());
    file.close();
  }

  return coordinates;
}

//! Return indices of nodes of cells in a mesh from input file
template <unsigned Tdim>
std::vector<std::vector<mpm::Index>> mpm::ReadMeshAscii<Tdim>::read_mesh_cells(
    const std::string& mesh) {
  // Indices of nodes
  std::vector<std::vector<mpm::Index>> cells;
  cells.clear();

  std::fstream file;
  file.open(mesh.c_str(), std::ios::in);

  try {
    if (file.is_open() && file.good()) {
      std::string line;
      // bool to check firstline
      bool read_first_line = false;
      // Number of coordinate lines
      unsigned nlines = 0;
      // Coordinates
      Eigen::Matrix<double, Tdim, 1> coords;
      // # of nodes and cells
      mpm::Index nnodes = 0, ncells = 0;
      // ignore stream
      double ignore;

      while (std::getline(file, line)) {
        boost::algorithm::trim(line);
        std::istringstream istream(line);
        // Vector of node ids for a cell
        std::vector<mpm::Index> nodes;
        nodes.clear();
        // ignore comment lines (# or !) or blank lines
        if ((line.find('#') == std::string::npos) &&
            (line.find('!') == std::string::npos) && (line != "")) {
          while (istream.good()) {
            if (!read_first_line) {
              // Read number of nodes and cells
              istream >> nnodes >> ncells;
              cells.reserve(ncells);
              read_first_line = true;
              break;
            }
            // Ignore nodal coordinates
            if (nlines > nnodes) {
              // Read node ids of each cell
              mpm::Index nid;
              istream >> nid;
              nodes.emplace_back(nid);
            } else {
              // Ignore stream not related to node ids of cells
              istream >> ignore;
            }
          }
          ++nlines;
          // Check if nodes is not empty, before adding to cell
          if (!nodes.empty()) {
            cells.emplace_back(nodes);
          }
        }
      }
    }
    file.close();
  } catch (std::exception& exception) {
    console_->error("Read mesh cells: {}", exception.what());
    file.close();
  }

  return cells;
}

//! Return coordinates of particles
template <unsigned Tdim>
std::vector<Eigen::Matrix<double, Tdim, 1>>
    mpm::ReadMeshAscii<Tdim>::read_particles(
        const std::string& particles_file) {

  // Nodal coordinates
  std::vector<VectorDim> coordinates;
  coordinates.clear();

  // Expected number of particles
  mpm::Index nparticles;

  // bool to check firstline
  bool read_first_line = false;

  // input file stream
  std::fstream file;
  file.open(particles_file.c_str(), std::ios::in);

  try {
    if (file.is_open() && file.good()) {
      // Line
      std::string line;
      while (std::getline(file, line)) {
        boost::algorithm::trim(line);
        std::istringstream istream(line);
        // ignore comment lines (# or !) or blank lines
        if ((line.find('#') == std::string::npos) &&
            (line.find('!') == std::string::npos) && (line != "")) {
          while (istream.good()) {
            if (!read_first_line) {
              // Read number of nodes and cells
              istream >> nparticles;
              coordinates.reserve(nparticles);
              read_first_line = true;
              break;
            }
            // Coordinates
            Eigen::Matrix<double, Tdim, 1> coords;
            // Read to coordinates
            for (unsigned i = 0; i < Tdim; ++i) istream >> coords[i];
            coordinates.emplace_back(coords);
            break;
          }
        }
      }
    }
    file.close();
  } catch (std::exception& exception) {
    console_->error("Read particle coordinates: {}", exception.what());
    file.close();
  }

  return coordinates;
}

//! Return stresses of particles
template <unsigned Tdim>
std::vector<Eigen::Matrix<double, 6, 1>>
    mpm::ReadMeshAscii<Tdim>::read_particles_stresses(
        const std::string& particles_stresses) {

  // Nodal stresses
  std::vector<Eigen::Matrix<double, 6, 1>> stresses;
  stresses.clear();

  // Expected number of particles
  mpm::Index nparticles;

  // bool to check firstline
  bool read_first_line = false;

  // input file stream
  std::fstream file;
  file.open(particles_stresses.c_str(), std::ios::in);

  try {
    if (file.is_open() && file.good()) {
      // Line
      std::string line;
      while (std::getline(file, line)) {
        boost::algorithm::trim(line);
        std::istringstream istream(line);
        // ignore comment lines (# or !) or blank lines
        if ((line.find('#') == std::string::npos) &&
            (line.find('!') == std::string::npos) && (line != "")) {
          while (istream.good()) {
            if (!read_first_line) {
              // Read number of nodes and cells
              istream >> nparticles;
              stresses.reserve(nparticles);
              read_first_line = true;
              break;
            }
            // Stresses
            Eigen::Matrix<double, 6, 1> stress;
            // Read to stress
            for (unsigned i = 0; i < stress.size(); ++i) istream >> stress[i];
            stresses.emplace_back(stress);
            break;
          }
        }
      }
    }
    file.close();
  } catch (std::exception& exception) {
    console_->error("Read particle stresses: {}", exception.what());
    file.close();
  }
  return stresses;
}

//! Return velocity constraints of particles
template <unsigned Tdim>
std::vector<std::tuple<mpm::Index, unsigned, double>>
    mpm::ReadMeshAscii<Tdim>::read_velocity_constraints(
        const std::string& velocity_constraints_file) {

  // Nodal velocity constraints
  std::vector<std::tuple<mpm::Index, unsigned, double>> constraints;
  constraints.clear();

  // input file stream
  std::fstream file;
  file.open(velocity_constraints_file.c_str(), std::ios::in);

  try {
    if (file.is_open() && file.good()) {
      // Line
      std::string line;
      while (std::getline(file, line)) {
        boost::algorithm::trim(line);
        std::istringstream istream(line);
        // ignore comment lines (# or !) or blank lines
        if ((line.find('#') == std::string::npos) &&
            (line.find('!') == std::string::npos) && (line != "")) {
          while (istream.good()) {
            // ID
            mpm::Index id;
            // Direction
            unsigned dir;
            // Velocity
            double velocity;
            // Read stream
            istream >> id >> dir >> velocity;
            constraints.emplace_back(std::make_tuple(id, dir, velocity));
          }
        }
      }
    }
    file.close();
  } catch (std::exception& exception) {
    console_->error("Read velocity constraints: {}", exception.what());
    file.close();
  }
  return constraints;
}

//! Return particles traction
template <unsigned Tdim>
std::vector<std::tuple<mpm::Index, unsigned, double>>
    mpm::ReadMeshAscii<Tdim>::read_particles_tractions(
        const std::string& traction_file) {

  // particle tractions
  std::vector<std::tuple<mpm::Index, unsigned, double>> tractions;
  tractions.clear();

  // input file stream
  std::fstream file;
  file.open(traction_file.c_str(), std::ios::in);

  try {
    if (file.is_open() && file.good()) {
      // Line
      std::string line;
      while (std::getline(file, line)) {
        boost::algorithm::trim(line);
        std::istringstream istream(line);
        // ignore comment lines (# or !) or blank lines
        if ((line.find('#') == std::string::npos) &&
            (line.find('!') == std::string::npos) && (line != "")) {
          while (istream.good()) {
            // ID
            mpm::Index id;
            // Direction
            unsigned dir;
            // Traction
            double traction;
            // Read stream
            istream >> id >> dir >> traction;
            tractions.emplace_back(std::make_tuple(id, dir, traction));
          }
        }
      }
    }
    file.close();
  } catch (std::exception& exception) {
    console_->error("Read traction : {}", exception.what());
    file.close();
  }
  return tractions;
}
