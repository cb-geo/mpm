#ifndef MPM_VECTOR_H_
#define MPM_VECTOR_H_

#include <algorithm>
#include <vector>

// TBB
#include <tbb/concurrent_vector.h>

#include "data_types.h"

namespace mpm {

// vector class
//! \brief A class that offers a vector and iterators
//! \tparam T A class with a template argument Tdim
template <class T>
class Vector {
 public:
  //! Default constructor
  Vector<T>() = default;

  //! Add a pointer to an element
  //! \param[in] ptr A shared pointer
  //! \param[in] check_duplicates Parameter to check duplicates
  bool add(const std::shared_ptr<T>&, bool check_duplicates = true);

  //! Remove an element pointer
  //! \param[in] ptr A shared pointer
  bool remove(const std::shared_ptr<T>&);

  //! Return number of elements in the vector
  std::size_t size() const { return elements_.size(); }

  //! Reserve the size of vector
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

  //! Iterate over elements in the vector
  //! \tparam T A class with a template argument Tdim
  //! \tparam Tunaryfn A unary function

  template <class Tunaryfn>
  Tunaryfn for_each(Tunaryfn fn);

 private:
  // Unordered map of index and pointer
  tbb::concurrent_vector<std::shared_ptr<T>> elements_;
};  // Vector class

#include "vector.tcc"

}  // namespace mpm
#endif  // MPM_VECTOR_H_
