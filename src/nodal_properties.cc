#include "nodal_properties.h"

// Function to create new property with given name and size (rows x cols)
bool mpm::NodalProperties::create_property(const std::string& property,
                                           unsigned rows, unsigned columns) {
  // Create a matrix with size of rows times columns and insert it to the
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
  mpm::MapProperty property_handle(position, nprops);
  return property_handle;
}

// Assign property value to a pair of node and material
void mpm::NodalProperties::assign_property(
    const std::string& property, unsigned node_id, unsigned mat_id,
    const Eigen::MatrixXd& property_value, unsigned nprops) {
  // Assign a property value matrix to its proper location in the properties_
  // matrix that stores all nodal properties
  properties_.at(property).block(node_id * nprops, mat_id,
                                 property_value.rows(), property_value.cols()) =
      property_value;
}

// Update property value according to a pair of node and material
void mpm::NodalProperties::update_property(
    const std::string& property, unsigned node_id, unsigned mat_id,
    const Eigen::MatrixXd& property_value, unsigned nprops) {
  // Update a property value matrix with dimensions nprops x 1 considering its
  // proper location in the properties_ matrix that stores all nodal properties
  properties_.at(property).block(node_id * nprops, mat_id, nprops, 1) =
      property_value + this->property(property, node_id, mat_id, nprops);
}

// Initialise all the nodal values for all properties in the property pool
void mpm::NodalProperties::initialise_nodal_properties() {
  // Iterate over all properties in the property map
  for (auto prop_itr = properties_.begin(); prop_itr != properties_.end();
       ++prop_itr) {
    // Create Matrix with zero values that has same size of the current property
    // in the iteration. The referred size is equal to rows * cols, where:
    // rows = number of nodes * size of property (1 if property is scalar, Tdim
    // if property is vector)
    // cols = number of materials
    Eigen::MatrixXd zeroed_property =
        Eigen::MatrixXd::Zero(prop_itr->second.rows(), prop_itr->second.cols());
    if (prop_itr->first != "current_velocities")
      this->assign_property(prop_itr->first, 0, 0, zeroed_property);
  }
}
