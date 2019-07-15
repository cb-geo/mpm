#ifndef MPM_SINE_LOAD_H_
#define MPM_SINE_LOAD_H_

#include <map>

#include "load_base.h"

namespace mpm {

//! Sine Traction class
//! \brief class that computes traction from a mathematical function
//! \details Sine Traction class computes the traction force from a sinusoidal
//! function
class SineLoad : public LoadBase {
 public:
  // Construct a Sine traction with a global unique id
  //! \param[in] id Global id
  SineLoad(unsigned id) : LoadBase(id){}

  //! Default destructor
  ~SineLoad() override{};

  //! Delete copy constructor
  SineLoad(const SineLoad&) = delete;

  //! Delete assignement operator
  SineLoad& operator=(const SineLoad&) = delete;

  double value(const double current_time,
               const double magnitude) const override{return 0.0;};

 private:
  //! index
  using LoadBase::id_;
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
};  // StepTraction class
}  // namespace mpm


#endif  // MPM_SINE_TRACTION_H_
