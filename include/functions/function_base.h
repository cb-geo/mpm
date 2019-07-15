#ifndef MPM_FUNCTION_BASE_H_
#define MPM_FUNCTION_BASE_H_

#include <vector>
#include <limits>
#include <stdexcept>

#include "logger.h"

namespace mpm {

//! Load Base class
//! \brief Base class that implements a mathematical function
//! \details function base class to compute the value of a mathematical function
//! at a given input
class FunctionBase {
 public:
  // Construct a Function Base with a unique id
  //! \param[in] id Global id
  FunctionBase(unsigned id){};

  //! Default destructor
  virtual ~FunctionBase(){};

  //! Delete copy constructor
  FunctionBase(const FunctionBase&) = delete;

  //! Delete assignement operator
  FunctionBase& operator=(const FunctionBase&) = delete;

  //! Return the value of the function at given input
  //! \param[in] input x
  //! \retval f(x)
  virtual double value(const double x_input) = 0;
};  // FunctionBase class
}  // namespace mpm


#endif  // MPM_FUNCTION_BASE_H_
