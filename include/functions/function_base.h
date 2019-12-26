#ifndef MPM_FUNCTION_BASE_H_
#define MPM_FUNCTION_BASE_H_

#include <limits>
#include <stdexcept>
#include <vector>

#include "json.hpp"

// JSON
using Json = nlohmann::json;

namespace mpm {

//! Mathematical function base class
//! \brief Base class that implements a mathematical function
//! \details function base class to compute the value of a mathematical function
//! at a given input
class FunctionBase {
 public:
  //! Construct a Function Base with a unique id
  //! \param[in] id Global id
  //! \param[in] json object of function properties
  FunctionBase(unsigned id, const Json& function_properties) : id_{id} {};

  //! Default destructor
  virtual ~FunctionBase(){};

  //! Delete copy constructor
  FunctionBase(const FunctionBase&) = delete;

  //! Delete assignement operator
  FunctionBase& operator=(const FunctionBase&) = delete;

  //! Return id of the function
  unsigned id() const { return id_; }

  //! Return the value of the function at given input
  //! \param[in] input x
  //! \retval f(x)
  virtual double value(double x_input) = 0;

 protected:
  //! function id
  unsigned id_{std::numeric_limits<unsigned>::max()};
  // Function properties
  Json properties_;
};  // FunctionBase class
}  // namespace mpm

#endif  // MPM_FUNCTION_BASE_H_
