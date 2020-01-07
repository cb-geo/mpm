#ifndef MPM_SIN_FUNCTION_H_
#define MPM_SIN_FUNCTION_H_

#include <cmath>

#include "function_base.h"

namespace mpm {

//! Sin Function class
//! \brief class that computes a mathematical sin function
//! \details Sin Function computes the value of a sin function
class SinFunction : public FunctionBase {
 public:
  // Construct a Sin function with a unique id
  //! \param[in] id Global id
  //! \param[in] json object of function properties
  SinFunction(unsigned id, const Json& properties);

  //! Default destructor
  ~SinFunction() override{};

  //! Delete copy constructor
  SinFunction(const SinFunction&) = delete;

  //! Delete assignement operator
  SinFunction& operator=(const SinFunction&) = delete;

  //! Return the value of the sin function at given input
  //! \param[in] input x
  //! \retval f(x) = sin (a * (x - x0))
  double value(double x) const override;

 private:
  //! function id
  using FunctionBase::id_;
  //! x0
  double x0_{0.};
  //! a_
  double a_{0.};
};  // SinFunction class
}  // namespace mpm

#endif  // MPM_SIN_FUNCTION_H_
