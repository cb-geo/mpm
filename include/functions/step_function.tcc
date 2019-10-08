// Constructor
mpm::StepFunction::StepFunction(unsigned id, const Json& function_properties)
    : mpm::FunctionBase(id, function_properties) {
  properties_ = function_properties;
  // get tabular values for the linear step relationship
  if (properties_.at("xvalues").is_array() &&
      properties_.at("fxvalues").is_array() &&
      properties_.at("xvalues").size() == properties_.at("fxvalues").size()) {
    std::vector<double> xvalue = properties_.at("xvalues");
    std::vector<double> fx = properties_.at("fxvalues");
    for (unsigned index = 0; index < xvalue.size(); ++index) {
      xvalues_.insert({index, xvalue.at(index)});
      fxvalues_.insert({index, fx.at(index)});
    }
  } else
    throw std::runtime_error(
        "Cannot initialise step function; x and f(x) are invalid");
}

// Return f(x) for a given x
double mpm::StepFunction::value(const double x_input) {
  // Check if the linear relationship for the step function is defined, if not
  // return an error
  if (xvalues_.empty() || fxvalues_.empty())
    throw std::runtime_error(
        "Cannot find the f(x); no linear function is defined");
  // If the given 'x' is less than x-begin, return f(x-begin)
  if (x_input < (xvalues_.cbegin()->second)) return fxvalues_.cbegin()->second;
  // If the given 'x' is greater than x-end, return f(x-end)
  if (x_input > (xvalues_.rbegin()->second)) return fxvalues_.rbegin()->second;
  // If the given 'x' is within the range, compute the relevant f(x)
  for (unsigned i = 0; i < (xvalues_.size() - 1); ++i) {
    double x_factor =
        (x_input - xvalues_.at(i)) / (xvalues_.at(i + 1) - xvalues_.at(i));
    if (x_factor > (-std::numeric_limits<double>::epsilon()) &&
        x_factor < (1.0 + std::numeric_limits<double>::epsilon()))
      return (fxvalues_.at(i) +
              x_factor * (fxvalues_.at(i + 1) - fxvalues_.at(i)));
  }
}
