// Initialise
template <unsigned Tdim>
bool mpm::MPMExplicit<Tdim>::initialise() {
  // Unique id
  uuid_ = boost::lexical_cast<std::string>(boost::uuids::random_generator()());
  return true;
}
