//! Constructor with mesh pointer and generator properties
template <unsigned Tdim>
mpm::GaussPointGenerator<Tdim>::GaussPointGenerator(
    const std::shared_ptr<mpm::Mesh<Tdim>>& mesh,
    const std::shared_ptr<mpm::IO>& io, const Json& generator)
    : PointGenerator<Tdim>(mesh, io, generator) {}

//! Generate material points
template <unsigned Tdim>
std::vector<Eigen::Matrix<double, Tdim, 1>>
    mpm::GaussPointGenerator<Tdim>::generate_points() {
  std::vector<Eigen::Matrix<double, Tdim, 1>> coordinates;
  // Get number of particles per cell
  unsigned nparticles_cell =
      generator_["nparticles_cells"].template get<unsigned>();

  if (nparticles_cell > 0)
    coordinates = mesh_->generate_material_points(nparticles_cell);
  return coordinates;
}
