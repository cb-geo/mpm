#ifndef MPM_INTERFACE_INTERFACE_H_
#define MPM_INTERFACE_INTERFACE_H_

#include <limits>
#include <set>

#include "Eigen/Dense"
#include "json.hpp"

#include "logger.h"

// JSON
using Json = nlohmann::json;

namespace mpm {

//! Interface base class
//! \brief Interface class that stores the information about interface law
//! \details Interface class with interface law
class Interface {
 public:
  // Constructor with id
  //! \param[in] id Interface id
  Interface(unsigned id, const Json& interface_properties);

  //! Destructor
  virtual ~Interface(){};

  //! Delete copy constructor
  Interface(const Interface&) = delete;

  //! Delete assignement operator
  Interface& operator=(const Interface&) = delete;

  //! Return id of the interface
  unsigned id() const { return id_; }

  //! Return friction coefficient of the interface
  double friction() const { return friction_; }

 private:
  //! interface id
  unsigned id_{std::numeric_limits<unsigned>::max()};
  //! Friction coefficient
  double friction_{0};
  //! Set of materials that share this interface law
  std::set<unsigned> materials_;
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
};  // Interface class
}  // namespace mpm

#endif  // MPM_INTERFACE_INTERFACE_H_
