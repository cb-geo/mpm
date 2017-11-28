#ifndef MPM_CONTAINER_H_
#define MPM_CONTAINER_H_

#include <algorithm>
#include <vector>

namespace mpm {

// Global index type for the node
using Index = unsigned long long;

// container class
//! \brief A class that offers a container and iterators
//! \tparam T A class with a template argument Tdim
//! \tparam Tdim Dimension
template <class T>
class Container {
 public:
  //! Default constructor
  Container<T>() = default;
  
  //! Insert a pointer
  bool insert(const std::shared_ptr<T>&);

  //! Return number of elements in the container
  std::size_t size() const { return elements_.size(); }
  
  //! Return begin iterator of nodes
  typename std::vector<std::shared_ptr<T>>::const_iterator begin() const {
    return elements_.cbegin();
  }

  //! Return end iterator of nodes
  typename std::vector<std::shared_ptr<T>>::const_iterator end() const {
    return elements_.cend();
  }

  //! Iterate over elements in the container
  template <class Tunaryfn>
  Tunaryfn for_each(Tunaryfn fn);

 private:
  // Unordered map of index and pointer
  std::vector<std::shared_ptr<T>> elements_;
};  // Container class

#include "container.tcc"
  
}  // mpm namespace
#endif  // MPM_CONTAINER_H_
