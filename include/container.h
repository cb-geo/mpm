#ifndef MPM_CONTAINER_H_
#define MPM_CONTAINER_H_

#include <algorithm>
#include <vector>

// TBB
#include <tbb/concurrent_vector.h>

namespace mpm {

// Global index type for the node
using Index = unsigned long long;

// container class
//! \brief A class that offers a container and iterators
//! \tparam T A class with a template argument Tdim
template <class T>
class Container {
 public:
  //! Default constructor
  Container<T>() = default;

  //! Add a pointer to an element
  bool add(const std::shared_ptr<T>&);

  //! Remove an element pointer
  bool remove(const std::shared_ptr<T>&);

  //! Return number of elements in the container
  std::size_t size() const { return elements_.size(); }

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

  //! Iterate over elements in the container
  template <class Tunaryfn>
  Tunaryfn for_each(Tunaryfn fn);

 private:
  // Unordered map of index and pointer
  tbb::concurrent_vector<std::shared_ptr<T>> elements_;
};  // Container class

#include "container.tcc"

}  // namespace mpm
#endif  // MPM_CONTAINER_H_
