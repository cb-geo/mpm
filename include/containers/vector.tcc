//! Add an element pointer
template <class T>
bool mpm::Vector<T>::add(const std::shared_ptr<T>& ptr, bool check_duplicates) {
  bool insertion_status = false;
  if (check_duplicates) {
    // Check if it is found in the Vector
    auto itr = std::find_if(this->cbegin(), this->cend(),
                            [ptr](std::shared_ptr<T> const& element) {
                              return element->id() == ptr->id();
                            });

    if (itr == this->cend()) {
      elements_.push_back(ptr);
      insertion_status = true;
    }
  } else {
    elements_.push_back(ptr);
    insertion_status = true;
  }
  return insertion_status;
}

//! Remove a pointer
template <class T>
bool mpm::Vector<T>::remove(const std::shared_ptr<T>& ptr) {
  auto size = elements_.size();
  // Check if it is found in the Vector
  elements_.erase(std::remove(elements_.begin(), elements_.end(), ptr),
                  elements_.end());
  return !(size == elements_.size());
}

//! Iterate over elements in the Vector
template <class T>
template <class Tunaryfn>
Tunaryfn mpm::Vector<T>::for_each(Tunaryfn fn) {
  return std::for_each(elements_.begin(), elements_.end(), fn);
}
