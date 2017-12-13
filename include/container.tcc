//! Add an element pointer
//! \param[in] ptr A shared pointer
//! \tparam T A class with a template argument Tdim
template <class T>
bool mpm::Container<T>::add(const std::shared_ptr<T>& ptr) {
  bool insertion_status = false;
  // Check if it is found in the container
  auto itr = std::find_if(this->begin(), this->end(),
                          [ptr](std::shared_ptr<T> const& element) {
                            return element->id() == ptr->id();
                          });

  if (itr == this->end()) {
    elements_.push_back(ptr);
    insertion_status = true;
  }

  return insertion_status;
}

//! Remove a pointer
//! \param[in] ptr A shared pointer
//! \tparam T A class with a template argument Tdim
template <class T>
bool mpm::Container<T>::add(const std::shared_ptr<T>& ptr) {
  bool insertion_status = false;
  // Check if it is found in the container
  auto itr = std::find_if(this->begin(), this->end(),
                          [ptr](std::shared_ptr<T> const& element) {
                            return element->id() == ptr->id();
                          });

  if (itr == this->end()) {
    elements_.push_back(ptr);
    insertion_status = true;
  }

  return insertion_status;
}

//! Iterate over elements in the container
//! \tparam T A class with a template argument Tdim
//! \tparam Tunaryfn A unary function
template <class T>
template <class Tunaryfn>
Tunaryfn mpm::Container<T>::for_each(Tunaryfn fn) {
  return std::for_each(this->begin(), this->end(), fn);
}
