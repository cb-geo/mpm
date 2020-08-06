#ifndef _MPM_INTERFACE_
#define _MPM_INTERFACE_

#include "mesh.h"

namespace mpm {

//! Interface class
//! \brief Interface base class
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Interface {
 public:
  //! Default constructor with mesh class
  Interface(const std::shared_ptr<mpm::Mesh<Tdim>>& mesh);

  //! Intialize
  virtual inline void initialise(){};

  //! Compute contact forces
  virtual inline void compute_contact_forces(){};

 protected:
  //! Mesh object
  std::shared_ptr<mpm::Mesh<Tdim>> mesh_;
  //! MPI Size
  int mpi_size_ = 1;
  //! MPI rank
  int mpi_rank_ = 0;
};  // Interface class
}  // namespace mpm

#include "interface.tcc"

#endif  // _MPM_INTERFACE_
