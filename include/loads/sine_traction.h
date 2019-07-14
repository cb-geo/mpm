#ifndef MPM_SINE_TRACTION_H_
#define MPM_SINE_TRACTION_H_

#include <map>

#include "load_base.h"

namespace mpm {

//! Sine Traction class
//! \brief class that computes traction from a mathematical function
//! \details Sine Traction class computes the traction force from a sinusoidal
//! function
//! \tparam Tdim Dimension
template <unsigned Tdim>
class SineTraction : public LoadBase<Tdim> {
 public:
  // Construct a Sine traction with a global unique id
  //! \param[in] id Global id
  SineTraction(unsigned id);

  //! Default destructor
  ~SineTraction() override{};

  //! Delete copy constructor
  SineTraction(const SineTraction<Tdim>&) = delete;

  //! Delete assignement operator
  SineTraction& operator=(const SineTraction<Tdim>&) = delete;

  double value(const double current_time,
               const double magnitude) const override{return 0.0;};

 private:
  unsigned id_;
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
};  // StepTraction class
}  // namespace mpm


#endif  // MPM_SINE_TRACTION_H_
