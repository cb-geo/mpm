#include "nodal_properties.h"

// Function to create new property with given name and size (rows x cols)
bool mpm::NodalProperties::create_property(const std::string& property,
                                           unsigned rows, unsigned columns) {
  // Initialize a matrix with size of rows times columns and insert it to the
  // property database map
  Eigen::MatrixXd property_data = Eigen::MatrixXd::Zero(rows, columns);
  std::pair<std::map<std::string, Eigen::MatrixXd>::iterator, bool> status =
      properties_.insert(
          std::pair<std::string, Eigen::MatrixXd>(property, property_data));
  return status.second;
}

// Return data in the nodal properties map at a specific index
Eigen::MatrixXd mpm::NodalProperties::property(const std::string& property,
                                               unsigned node_id,
                                               unsigned mat_id,
                                               unsigned nprops) const {
  // Const pointer to location of property: node_id * nprops x mat_id
  const double* position = &properties_.at(property)(node_id * nprops, mat_id);
  mpm::MapProperty property_map(position, nprops);
  return property_map;
}

// Assign property value to a pair of node and material
void mpm::NodalProperties::assign_property(const std::string& property,
                                           unsigned node_id, unsigned mat_id,
                                           Eigen::MatrixXd property_value,
                                           unsigned nprops) {
  // Assign a property value matrix with dimensions nprops x 1 to its proper
  // location in the properties_ matrix that stores all nodal properties
  properties_.at(property).block(node_id * nprops, mat_id, nprops, 1) =
      property_value;
}
