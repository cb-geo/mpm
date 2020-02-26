#include "nodal_properties.h"

// Constructor
mpm::NodalProperties::NodalProperties() {
  // Initialize map of properties as an empty map
  properties_ = std::map<std::string, Eigen::MatrixXd>();
}

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

// Fetch data in the nodal properties map
Eigen::MatrixXd mpm::NodalProperties::property_data(const std::string& property,
                                                    unsigned id_node,
                                                    unsigned id_mat,
                                                    unsigned nprops = 1) {
  return properties_.at(property).block(id_node * nprops, id_mat, nprops, 1);
}

// Assign property value to a pair of node and material
bool mpm::NodalProperties::assign_property(const std::string& property,
                                           Eigen::MatrixXd property_value,
                                           unsigned id_node, unsigned id_mat,
                                           unsigned nprops = 1) {
  properties_.at(property).block(id_node * nprops, id_mat, nprops, 1) =
      property_value;
  return true;
}
