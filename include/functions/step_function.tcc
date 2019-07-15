// Constructor
mpm::StepFunction::StepFunction(unsigned id, const std::vector<double>& xvalue,
                            const std::vector<double>& fx)
    : mpm::FunctionBase(id) {
  // unique id
  id_ = id;
  // get tabular values for the linear step relationship
  if (xvalue.size() != fx.size())
    throw std::runtime_error(
        "Cannot initialise step function; x and f(x) are invalid");
  for (unsigned index = 0; index < xvalue.size(); ++index) {
    xvalues_.insert({index, xvalue.at(index)});
    fxvalues_.insert({index, fx.at(index)});
  }
}

// Return f(x) for a given x
double mpm::StepFunction::value(const double x_input) {
  // Check if the linear relationship for the step function is defined, if not
  // return an error
  if (xvalues_.empty() || fxvalues_.empty())
    throw std::runtime_error(
        "Cannot find the f(x); no linear function is defined");
  // If the given 'x' is less than x-begin, return f(x-begin)
  if(x_input < (xvalues_.cbegin()->first))
    return fxvalues_.cbegin()->first;
  // If the given 'x' is greater than x-end, return f(x-end)
  if(x_input > (xvalues_.rbegin()->first))
    return fxvalues_.rbegin()->first;
  // If the given 'x' is within the range, compute the relevant f(x)
  for (unsigned i = 0; i < (xvalues_.size() - 1); ++i) {
    double x_factor =
        (x_input - xvalues_.at(i)) / (xvalues_.at(i + 1) - xvalues_.at(i));
    if (x_factor > (-std::numeric_limits<double>::epsilon()) && x_factor < 1.0)
      return (fxvalues_.at(i) +
              x_factor * (fxvalues_.at(i + 1) - fxvalues_.at(i)));
  }
}
