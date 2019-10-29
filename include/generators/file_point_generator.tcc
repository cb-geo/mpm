//! Constructor with mesh pointer and generator properties
template <unsigned Tdim>
mpm::FilePointGenerator<Tdim>::FilePointGenerator(
    const std::shared_ptr<mpm::Mesh<Tdim>>& mesh,
    const std::shared_ptr<mpm::IO>& io, const Json& generator)
    : PointGenerator<Tdim>(mesh, io, generator) {}

//! Generate material points
template <unsigned Tdim>
std::vector<Eigen::Matrix<double, Tdim, 1>>
    mpm::FilePointGenerator<Tdim>::generate_points() {

  // Get particle reader from JSON object
  const std::string reader = generator_["reader"].template get<std::string>();

  // Create a particle reader
  auto particle_reader =
      Factory<mpm::ReadMesh<Tdim>>::instance()->create(reader);

  // Read particles from file : this needs modification in IO class
  const auto pfile =
      io_->working_dir() + generator_["location"].template get<std::string>();
  return particle_reader->read_particles(pfile);
}
