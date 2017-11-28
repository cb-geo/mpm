#ifndef MPM_HANDLER_H_
#define MPM_HANDLER_H_

#include <algorithm>
#include <unordered_map>

namespace mpm {

// Global index type for the node
using Index = unsigned long long;

// handler class
//! \brief A class that offers a container and iterators
//! \tparam T A class with a template argument Tdim
//! \tparam Tdim Dimension
template <class T>
class Handler {
 public:
  //! Default constructor
  Handler<T>() = default;

  //! Insert a pointer
  bool insert(const std::shared_ptr<T>&);

  //! Insert an id and a pointer
  bool insert(Index, const std::shared_ptr<T>&);

  //! Return number of elements in the container
  std::size_t size() const { return elements_.size(); }

  //! Return begin iterator of nodes
  typename std::unordered_map<Index, std::shared_ptr<T>>::const_iterator begin()
      const {
    return elements_.cbegin();
  }

  //! Return end iterator of nodes
  typename std::unordered_map<Index, std::shared_ptr<T>>::const_iterator end()
      const {
    return elements_.cend();
  }

  //! Iterate over elements in the container
  template <class Tunaryfn>
  Tunaryfn for_each(Tunaryfn fn);

 private:
  // Unordered map of index and pointer
  std::unordered_map<Index, std::shared_ptr<T>> elements_;
};  // Handler class

#include "handler.tcc"

}  // mpm namespace
#endif  // MPM_HANDLER_H_
