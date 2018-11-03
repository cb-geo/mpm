#ifndef MPM_MAP_H_
#define MPM_MAP_H_

#include <algorithm>
#include <unordered_map>
#include <tsl/robin_map.h>

namespace mpm {

// Global index type for the node
using Index = unsigned long long;

// Map class
//! \brief A class that offers a container and iterators
//! \tparam T A class with a template argument Tdim
template <class T>
class Map {
 public:
  //! Default constructor
  Map<T>() = default;

  //! Insert a pointer
  //! \param[in] ptr A shared pointer
  bool insert(const std::shared_ptr<T>& ptr);

  //! Insert a pointer at a given id
  //! \param[in] id Global/local index of the pointer
  //! \param[in] ptr A shared pointer
  bool insert(Index id, const std::shared_ptr<T>& ptr);

  //! Remove a pointer at a given id
  //! \param[in] id Global/local index of the pointer
  bool remove(Index id);

  //! Return number of elements in the container
  std::size_t size() const { return elements_.size(); }

  //! Return value at a given index
  std::shared_ptr<T> operator[](Index id) const { return elements_.at(id); }

  //! Return begin iterator of nodes
  typename tsl::robin_map<Index, std::shared_ptr<T>>::const_iterator begin()
      const {
    return elements_.cbegin();
  }

  //! Return end iterator of nodes
  typename tsl::robin_map<Index, std::shared_ptr<T>>::const_iterator end()
      const {
    return elements_.cend();
  }

  //! Iterate over elements in the container
  //! \tparam Tunaryfn A unary function
  template <class Tunaryfn>
  Tunaryfn for_each(Tunaryfn fn);

 private:
  // Unordered map of index and pointer
  tsl::robin_map<Index, std::shared_ptr<T>> elements_;
};  // Map class

#include "map.tcc"

}  // namespace mpm
#endif  // MPM_MAP_H_
