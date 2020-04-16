#include "sin_function.h"
// Constructor
mpm::SinFunction::SinFunction(unsigned id, const Json& properties)
    : mpm::FunctionBase(id, properties) {
  x0_ = properties.at("x0");
  a_ = properties.at("a");
  if (properties.contains("xrange"))
    if (properties.at("xrange").is_array() &&
        properties.at("xrange").size() == 2) {
      for (unsigned index = 0; index < 2; ++index)
        xrange_.at(index) = properties.at("xrange").at(index);
    } else
      throw std::runtime_error(
          "Cannot initialise sine function; x range is invalid");
}

// Return f(x) for a given x
double mpm::SinFunction::value(double x_input) const {
  if ((x_input < xrange_[0]) || (x_input > xrange_[1])) return 0.0;
  return std::sin(a_ * (x_input - x0_));
}
