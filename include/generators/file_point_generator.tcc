#include <fstream>
//! Constructor with mesh pointer and generator properties
template <unsigned Tdim>
mpm::FilePointGenerator<Tdim>::FilePointGenerator(
    const std::shared_ptr<mpm::Mesh<Tdim>>& mesh,
    const std::shared_ptr<mpm::IO>& io, const Json& generator)
    : PointGenerator<Tdim>(mesh, io, generator) {
  console_->error("{} {}", __FILE__, __LINE__);
}

//! Generate material points
template <unsigned Tdim>
std::vector<Eigen::Matrix<double, Tdim, 1>>
    mpm::FilePointGenerator<Tdim>::generate_points() {
  console_->error("{} {}", __FILE__, __LINE__);

// Dump JSON as an input file to be read
  std::ofstream file;
  file.open("test.json");
  file << generator_.dump(2);
  file.close();


  console_->info("Reader {}", generator_["generator"]["reader"].template get<std::string>());
  // Get particle reader from JSON object
  const std::string reader = generator_["generator"]["reader"].template get<std::string>();
  console_->error("{} {}", __FILE__, __LINE__);

  console_->info("Reader {}", reader);
  
  // Create a particle reader
  auto particle_reader =
      Factory<mpm::ReadMesh<Tdim>>::instance()->create(reader);

  console_->error("{} {}", __FILE__, __LINE__);
  // Read particles from file : this needs modification in IO class
  const auto pfile =
      io_->working_dir() + generator_["generator"]["location"].template get<std::string>();
  console_->error("{} {}", __FILE__, __LINE__);
  return particle_reader->read_particles(pfile);
}
