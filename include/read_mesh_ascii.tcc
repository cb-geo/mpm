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
        std::istringstream istream(line);
        // ignore comment lines (# or !) or blank lines
        if ((line.find('#') == std::string::npos) &&
            (line.find('!') == std::string::npos) && (line != "")) {
          while (istream.good()) {
            if (!read_first_line) {
              // Read number of nodes and cells
              istream >> nnodes >> ncells;
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
  } catch (std::exception& exception) {
    console_->error("Read mesh nodes: {}", exception.what());
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
              read_first_line = true;
              break;
            }
            // Ignore nodal coordinates
            if (nlines > nnodes) {
              // Read node ids of each cell
              mpm::Index nid;
              if (!istream.str().empty()) {
                istream >> nid;
                nodes.emplace_back(nid);
              }
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
  } catch (std::exception& exception) {
    console_->error("Read mesh cells: {}", exception.what());
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

  // input file stream
  std::fstream file;
  file.open(particles_file.c_str(), std::ios::in);

  try {
    if (file.is_open() && file.good()) {
      // Line
      std::string line;
      while (std::getline(file, line)) {
        std::istringstream istream(line);
        // ignore comment lines (# or !) or blank lines
        if ((line.find('#') == std::string::npos) &&
            (line.find('!') == std::string::npos) && (line != "")) {
          // Coordinates
          Eigen::Matrix<double, Tdim, 1> coords;
          while (istream.good()) {
            // Read to coordinates
            for (unsigned i = 0; i < Tdim; ++i) istream >> coords[i];
            break;
          }
          coordinates.emplace_back(coords);
        }
      }
    }
  } catch (std::exception& exception) {
    console_->error("Read mesh nodes: {}", exception.what());
  }

  return coordinates;
}
