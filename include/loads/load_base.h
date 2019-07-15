#ifndef MPM_LOAD_BASE_H_
#define MPM_LOAD_BASE_H_

#include "logger.h"

namespace mpm {

//! Load Base class
//! \brief Base class that handles external loading
//! \details load base class to apply external loads on nodes and particles
class LoadBase {
 public:
  // Construct a Load Base with a global unique id
  //! \param[in] id Global id
  LoadBase(unsigned id){};

  //! Default destructor
  virtual ~LoadBase(){};

  //! Delete copy constructor
  LoadBase(const LoadBase&) = delete;

  //! Delete assignement operator
  LoadBase& operator=(const LoadBase&) = delete;

  virtual double value(const double current_time,
                       const double magnitude) const = 0;

protected:
  //! Index
  unsigned id_;
};  // LoadBase class
}  // namespace mpm


#endif  // MPM_LOAD_BASE_H_
