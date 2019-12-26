#ifndef MPM_PARTICLE_TRACTION_H_
#define MPM_PARTICLE_TRACTION_H_

//! Alias for JSON
#include "json.hpp"
using Json = nlohmann::json;

#include "function_base.h"

namespace mpm {

class ParticleTraction {
 public:
  // Constructor
  //! \param[setid] setid Particle set id
  //! \param[in] mfunction Math function if defined
  //! \param[dir] dir Direction of traction load
  //! \param[traction] traction Particle traction
  ParticleTraction(int setid,
                   const std::shared_ptr<mpm::FunctionBase>& traction_fn,
                   unsigned dir, double traction)
      : setid_{setid},
        traction_fn_{traction_fn},
        dir_{dir},
        traction_{traction} {};

  // Set id
  int setid() const { return setid_; }

  // Direction
  unsigned dir() const { return dir_; }

  // Return traction
  double traction(double current_time) const {
    return traction_ * (this->traction_fn_)->value(current_time);
  }

 private:
  // ID
  int setid_;
  // Math function
  std::shared_ptr<mpm::FunctionBase> traction_fn_;
  // Direction
  unsigned dir_;
  // Traction
  double traction_;
};
}  // namespace mpm
#endif  // MPM_PARTICLE_TRACTION_H_
