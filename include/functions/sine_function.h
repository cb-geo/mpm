#ifndef MPM_SINE_FUNCTION_H_
#define MPM_SINE_FUNCTION_H_

#include <map>

#include "function_base.h"

namespace mpm {

//! Sine Function class
//! \brief class that computes a sine function
//! \details Sine Function class computes the values of a sine function at a
//! given input
class SineFunction : public FunctionBase {
 public:
  // Construct a Sine function with a unique id
  //! \param[in] id Global id
  //! \param[in] json object of function properties
  SineFunction(unsigned id, const Json& function_properties)
    : FunctionBase(id, function_properties) {
    properties_ = function_properties;
  }

  //! Default destructor
  ~SineFunction() override{};

  //! Delete copy constructor
  SineFunction(const SineFunction&) = delete;

  //! Delete assignement operator
  SineFunction& operator=(const SineFunction&) = delete;

  double value(const double x_input) override {
    return 0.0;
  };

 private:
  //! function id
  using FunctionBase::id_;
  //! function properties
  using FunctionBase::properties_;
};  // SineFunction class
}  // namespace mpm


#endif  // MPM_SINE_FUNCTION_H_
