#ifndef MPM_CONTAINER_H_
#define MPM_CONTAINER_H_

#include <algorithm>
#include <vector>

// TBB
#include <tbb/concurrent_vector.h>

#include "data_types.h"

namespace mpm {

// container class
//! \brief A class that offers a container and iterators
//! \tparam T A class with a template argument Tdim
template <class T>
class Container {
 public:
  //! Default constructor
  Container<T>() = default;

  //! Add a pointer to an element
  //! \param[in] ptr A shared pointer
  //! \param[in] check_duplicates Parameter to check duplicates
  bool add(const std::shared_ptr<T>&, bool check_duplicates = true);

  //! Remove an element pointer
  //! \param[in] ptr A shared pointer
  bool remove(const std::shared_ptr<T>&);

  //! Return number of elements in the container
  std::size_t size() const { return elements_.size(); }

  //! Reserve the size of container
  void reserve(const mpm::Index size) { elements_.reserve(size); }

  //! Clear
  void clear() { elements_.clear(); }

  //! Return begin iterator of nodes
  typename tbb::concurrent_vector<std::shared_ptr<T>>::const_iterator cbegin()
      const {
    return elements_.cbegin();
  }

  //! Return end iterator of nodes
  typename tbb::concurrent_vector<std::shared_ptr<T>>::const_iterator cend()
      const {
    return elements_.cend();
  }

  //! Return value at a given index
  std::shared_ptr<T> operator[](Index id) const { return elements_.at(id); }

  //! Iterate over elements in the container
  //! \tparam T A class with a template argument Tdim
  //! \tparam Tunaryfn A unary function

  template <class Tunaryfn>
  Tunaryfn for_each(Tunaryfn fn);

 private:
  // Unordered map of index and pointer
  tbb::concurrent_vector<std::shared_ptr<T>> elements_;
};  // Container class

#include "container.tcc"

}  // namespace mpm
#endif  // MPM_CONTAINER_H_
