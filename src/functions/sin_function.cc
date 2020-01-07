#include "sin_function.h"
// Constructor
mpm::SinFunction::SinFunction(unsigned id, const Json& properties)
    : mpm::FunctionBase(id, properties) {
  x0_ = properties.at("x0");
  a_ = properties.at("a");
}

// Return f(x) for a given x
double mpm::SinFunction::value(double x_input) const {
  return std::sin(a_ * (x_input - x0_));
}
