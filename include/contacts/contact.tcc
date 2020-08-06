//! Constructor of contact with mesh
template <unsigned Tdim>
mpm::Contact<Tdim>::Contact(const std::shared_ptr<mpm::Mesh<Tdim>>& mesh) {
  // Assign mesh
  mesh_ = mesh;
  // Initialise MPI rank and size
  int mpi_rank_ = 0;
  int mpi_size_ = 1;

#ifdef USE_MPI
  // Get MPI rank
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
  // Get number of MPI ranks
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);
#endif
}
