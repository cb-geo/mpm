#include "interface.h"

//! Read interface properties
mpm::Interface::Interface(unsigned id, const Json& interface_properties)
    : id_{id} {
  //! Logger
  std::string logger = "interface::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);

  //! Read interface properties
  try {
    // Read Friction
    friction_ = interface_properties.at("friction").template get<double>();

    // Read set of nodes
    std::set<unsigned> materials_ = interface_properties.at("material_ids");

  } catch (Json::exception& except) {
    console_->error("Interface parameter not set: {} {}\n", except.what(),
                    except.id);
  }
}
