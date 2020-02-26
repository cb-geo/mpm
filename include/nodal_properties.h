#ifndef MPM_NODAL_PROPERTIES_H_
#define MPM_NODAL_PROPERTIES_H_

#include <Eigen/Dense>
#include <map>

namespace mpm {
// \brief Multimaterial parameters on each node
struct NodalProperties {

  //! Constructor
  NodalProperties();

  //! Function to create new property with given name and size (rows x cols)
  //! \param[in] property Property name
  //! \param[in] rows Number of nodes times the number of the dimension of the
  //! property (1 if scalar, Tdim if vector)
  //! \param[in] columns Number of materials
  bool create_property(const std::string& property, unsigned rows,
                       unsigned columns);

  // Fetch data in the nodal properties map
  // \param[in] property Property name
  // \param[in] id_node Id of the node within the property data
  // \param[in] id_mat Id of the material within the property data
  // \param[in] nprops Dimension of property (1 if scalar, Tdim if vector)
  Eigen::MatrixXd property_data(const std::string& property, unsigned id_node,
                                unsigned id_mat, unsigned nprops = 1);

  // Assign property value to a pair of node and material
  // \param[in] property_value Property value to be assigned
  // \param[in] property Property name
  // \param[in] id_node Id of the node within the property data
  // \param[in] id_mat Id of the material within the property data
  // \param[in] nprops Dimension of property (1 if scalar, Tdim if vector)
  bool assign_property(const std::string& property,
                       Eigen::MatrixXd property_value, unsigned id_node,
                       unsigned id_mat, unsigned nprops = 1);

  std::map<std::string, Eigen::MatrixXd> properties_;
};  // NodalProperties struct
}  // namespace mpm

#endif  // MPM_NODAL_PROPERTIES_H_
