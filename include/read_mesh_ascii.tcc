//! Return coordinates of nodes in a mesh from input file
template <unsigned Tdim>
std::vector<Eigen::Matrix<double, Tdim, 1>>
    mpm::ReadMeshAscii<Tdim>::read_mesh_nodes(const std::string& mesh) {
  // Nodal coordinates
  std::vector<std::array<double, Tdim>> coordinates;
  coordinates.clear();

  std::fstream file;
  file.open(mesh.c_str(), std::ios::in);

  try {
    if (file.is_open() && file.good()) {
      std::string line;
      Eigen::Matrix<double, Tdim, 1> coords;

      while (std::getline(file, line)) {
        std::istringstream istream(line);
        // ignore comment lines (# or !) or blank lines
        if ((line.find('#') == std::string::npos) &&
            (line.find('!') == std::string::npos) && (line != "")) {
          while (istream.good()) {
            // Read to coordinates
            for (unsigned i = 0; i < Tdim; ++i) istream >> coords[i];
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

//! Return indices of nodes of cells in a mesh from input file
template <unsigned Tdim>
std::vector<std::vector<unsigned>> mpm::ReadMeshAscii<Tdim>::read_mesh_cells(
    const std::string& mesh) {
  // Indices of nodes
  std::vector<std::vector<unsigned>> cells;
  cells.clear();

  std::fstream file;
  file.open(mesh.c_str(), std::ios::in);

  try {
    if (file.is_open() && file.good()) {
      std::string line;

      while (std::getline(file, line)) {
        // Create a vector of node ids of each cell
        std::vector<unsigned> nodes;
        nodes.clear();
        std::istringstream istream(line);
        // ignore comment lines (# or !) or blank lines
        if ((line.find('#') == std::string::npos) &&
            (line.find('!') == std::string::npos) && (line != "")) {
          while (istream.good()) {
            // Read node id
            unsigned nid;
            istream >> nid;
            nodes.emplace_back(nid);
          }
          cells.emplace_back(nodes);
        }
      }
    }
  } catch (std::exception& exception) {
    console_->error("Read mesh cells: {}", exception.what());
  }

  return cells;
}
