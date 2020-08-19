//! Constructor of contact with mesh
template <unsigned Tdim>
mpm::Contact<Tdim>::Contact(const std::shared_ptr<mpm::Mesh<Tdim>>& mesh) {
  // Assign mesh
  mesh_ = mesh;
}
