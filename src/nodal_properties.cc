#include "nodal_properties.h"

// Function to create new property with given name and size (rows x cols)
bool mpm::NodalProperties::create_property(const std::string& property,
                                           unsigned rows, unsigned columns) {
  // Initialize a matrix with size of rows times columns and insert it to the
  // property database map
  Eigen::MatrixXd property_data = Eigen::MatrixXd::Zero(rows, columns);
  properties_.insert(
      std::pair<std::string, Eigen::MatrixXd>(property, property_data));
  return true;
}

// Return data in the nodal properties map at a specific index
Eigen::MatrixXd mpm::NodalProperties::property(const std::string& property,
                                               unsigned node_id,
                                               unsigned mat_id,
                                               unsigned nprops = 1) const {
  return properties_.at(property).block(node_id * nprops, mat_id, nprops, 1);
}

// Assign property value to a pair of node and material
bool mpm::NodalProperties::assign_property(const std::string& property,
                                           unsigned node_id, unsigned mat_id,
                                           Eigen::MatrixXd property_value,
                                           unsigned nprops = 1) {
  properties_.at(property).block(node_id * nprops, mat_id, nprops, 1) =
      property_value;
  return true;
}
