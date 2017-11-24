//! Insert a pointer
//! \param[in] ptr A shared pointer 
//! \tparam T A class with a template argument Tdim
template<class T>
bool mpm::Handler<T>::insert(const std::shared_ptr<T>& ptr) {
  bool insertion_status = elements_.insert(std::make_pair(ptr->id(), ptr)).second;
  return insertion_status;
}

//! Insert a pointer at a given id
//! \param[in] id Global/local index of the pointer
//! \param[in] ptr A shared pointer 
//! \tparam T A class with a template argument Tdim
template<class T>
bool mpm::Handler<T>::insert(mpm::Index id, const std::shared_ptr<T>& ptr) {
  bool insertion_status = elements_.insert(std::make_pair(id, ptr)).second;
  return insertion_status;
}
