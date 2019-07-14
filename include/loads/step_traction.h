#ifndef MPM_STEP_TRACTION_H_
#define MPM_STEP_TRACTION_H_

#include <map>

#include "load_base.h"

namespace mpm {

//! Step Traction class
//! \brief class that computes traction from a mathematical function
//! \details Step Traction class computes the traction force from a step
//! function
//! \tparam Tdim Dimension
template <unsigned Tdim>
class StepTraction : public LoadBase<Tdim> {
 public:
  // Construct a Step traction with a global unique id
  //! \param[in] id Global id
  StepTraction(unsigned id, const std::map<double, double>& time_vs_load);

  //! Default destructor
  ~StepTraction() override{};

  //! Delete copy constructor
  StepTraction(const StepTraction<Tdim>&) = delete;

  //! Delete assignement operator
  StepTraction& operator=(const StepTraction<Tdim>&) = delete;

  double value(const double current_time,
               const double magnitude) const override{return 0.0;};

 private:
  unsigned id_;
  //! Tabular data of time and relative loads
  //! first->time (this is total time)
  //! second->relative load
  std::map<double, double> time_vs_load_;

  //! Logger
  std::unique_ptr<spdlog::logger> console_;
};  // StepTraction class
}  // namespace mpm


#endif  // MPM_STEP_TRACTION_H_
