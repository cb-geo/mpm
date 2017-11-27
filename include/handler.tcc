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

//! Iterate over elements in the container
//! \tparam T A class with a template argument Tdim
//! \tparam Tfn A functor
template <class T>
template <class Titr, class Tfn>
Tfn mpm::Handler<T>::for_each(Titr first, Titr last, Tfn fn) {
  for (; first != last; ++first) {
    fn((*first).second);
  }
  return fn;
}
