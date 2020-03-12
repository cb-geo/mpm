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
  bool removal_status = false;

  // Check if it is found in the Vector
  auto itr = std::find_if(this->cbegin(), this->cend(),
                          [ptr](std::shared_ptr<T> const& element) {
                            return element->id() == ptr->id();
                          });

  // If Itr is present create a new set of elements
  if (itr != this->cend()) {
    tbb::concurrent_vector<std::shared_ptr<T>> new_elements;
    new_elements.reserve(elements_.size() - 1);
    auto it = std::copy_if(elements_.begin(), elements_.end(),
                           std::back_inserter(new_elements),
                           [ptr](std::shared_ptr<T> const& element) {
                             return element->id() != ptr->id();
                           });

    elements_ = new_elements;
    removal_status = true;
  }
  return removal_status;
}

//! Iterate over elements in the Vector
template <class T>
template <class Tunaryfn>
Tunaryfn mpm::Vector<T>::for_each(Tunaryfn fn) {
  return std::for_each(elements_.begin(), elements_.end(), fn);
}
