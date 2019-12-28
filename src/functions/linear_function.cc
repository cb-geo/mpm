#include "linear_function.h"
// Constructor
mpm::LinearFunction::LinearFunction(unsigned id,
                                    const Json& function_properties)
    : mpm::FunctionBase(id, function_properties) {
  properties_ = function_properties;
  // get tabular values for the linear linear relationship
  if (properties_.at("xvalues").is_array() &&
      properties_.at("fxvalues").is_array() &&
      properties_.at("xvalues").size() == properties_.at("fxvalues").size()) {
    std::vector<double> xvalue = properties_.at("xvalues");
    std::vector<double> fx = properties_.at("fxvalues");
    for (unsigned index = 0; index < xvalue.size(); ++index)
      x_fx_.insert({index, std::make_pair(xvalue.at(index), fx.at(index))});
  } else
    throw std::runtime_error(
        "Cannot initialise linear function; x and f(x) are invalid");
}

// Return f(x) for a given x
double mpm::LinearFunction::value(double x_input) {
  // Check if the linear relationship for the linear function is defined, if not
  // return an error
  double value = std::numeric_limits<double>::quiet_NaN();

  if (x_fx_.empty())
    throw std::runtime_error(
        "Cannot find the f(x); no linear function is defined");
  // If the given 'x' is less than x-begin, return f(x_begin)
  if (x_input < (x_fx_.cbegin()->second.first))
    return x_fx_.cbegin()->second.second;
  // If the given 'x' is greater than x-end, return f(x_end)
  if (x_input > (x_fx_.rbegin()->second.first))
    return x_fx_.rbegin()->second.second;
  // If the given 'x' is within the range, compute the relevant f(x)
  for (unsigned i = 0; i < (x_fx_.size() - 1); ++i) {
    double x_factor = (x_input - x_fx_.at(i).first) /
                      (x_fx_.at(i + 1).first - x_fx_.at(i).first);
    if (x_factor > (-std::numeric_limits<double>::epsilon()) &&
        x_factor < (1.0 + std::numeric_limits<double>::epsilon()))
      return (x_fx_.at(i).second +
              x_factor * (x_fx_.at(i + 1).second - x_fx_.at(i).second));
  }
  return value;
}
