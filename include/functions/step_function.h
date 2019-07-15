#ifndef MPM_STEP_FUNCTION_H_
#define MPM_STEP_FUNCTION_H_

#include <map>

#include "function_base.h"

namespace mpm {

//! Step Function class
//! \brief class that computes a mathematical step function
//! \details Step Function class computes the value of a step function at given
//! input function
class StepFunction : public FunctionBase {
 public:
  // Construct a Step function with a unique id
  //! \param[in] id Global id
  //! \param[in] json object of function properties
  StepFunction(unsigned id, const Json& function_properties);

  //! Default destructor
  ~StepFunction() override{};

  //! Delete copy constructor
  StepFunction(const StepFunction&) = delete;

  //! Delete assignement operator
  StepFunction& operator=(const StepFunction&) = delete;

  //! Return the value of the step function at given input
  //! \param[in] input x
  //! \retval f(x)
  double value(const double x_input) override;

 private:
  //! function id
  using FunctionBase::id_;
  //! function properties
  using FunctionBase::properties_;
  //! Tabular data of x
  //! first->key (unique table index)
  //! second->x
  std::map<unsigned, double> xvalues_;
  //! Tabular data f(x)
  //! first->key (unique table index)
  //! second->f(x)
  std::map<unsigned, double> fxvalues_;
};  // StepFunction class
}  // namespace mpm

#include "step_function.tcc"

#endif  // MPM_STEP_TRACTION_H_
