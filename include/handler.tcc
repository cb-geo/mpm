//! Insert a pointer
//! \param[in] ptr A shared pointer 
//! \tparam T A class with a template argument Tdim
template<class T>
bool mpm::Handler<T>::insert(const std::shared_ptr<T>& ptr) {
  bool insertion_status = elements_.insert(std::make_pair(ptr->id(), ptr)).second;
  return insertion_status;
}
