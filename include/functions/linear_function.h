#ifndef MPM_LINEAR_FUNCTION_H_
#define MPM_LINEAR_FUNCTION_H_

#include <map>

#include "function_base.h"

namespace mpm {

//! Linear Function class
//! \brief class that computes a mathematical linear function
//! \details Linear Function class computes the value of a linear function at
//! given input function
class LinearFunction : public FunctionBase {
 public:
  // Construct a Linear function with a unique id
  //! \param[in] id Global id
  //! \param[in] json object of function properties
  LinearFunction(unsigned id, const Json& function_properties);

  //! Default destructor
  ~LinearFunction() override = default;

  //! Return the value of the linear function at given input
  //! \param[in] input x
  //! \retval f(x)
  double value(double x) const override;

 private:
  //! function id
  using FunctionBase::id_;
  //! Tabular data of x and f(x)
  std::map<unsigned, std::pair<double, double>> x_fx_;
};  // LinearFunction class
}  // namespace mpm

#endif  // MPM_LINEAR_TRACTION_H_
