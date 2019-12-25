//! Compute Jacobian with particle size and deformation gradient
template <unsigned Tdim>
inline Eigen::Matrix<double, Tdim, Tdim> mpm::Element<Tdim>::jacobian(
    const VectorDim& xi, const Eigen::MatrixXd& nodal_coordinates,
    const VectorDim& particle_size,
    const VectorDim& deformation_gradient) const {

  // Get gradient shape functions
  const Eigen::MatrixXd grad_shapefn =
      this->grad_shapefn(xi, particle_size, deformation_gradient);

  // Jacobian dx_i/dxi_j
  return (grad_shapefn.transpose() * nodal_coordinates);
}

//! Compute Jacobian local with particle size and deformation gradient
template <unsigned Tdim>
inline Eigen::Matrix<double, Tdim, Tdim> mpm::Element<Tdim>::jacobian_local(
    const VectorDim& xi, const Eigen::MatrixXd& nodal_coordinates,
    const VectorDim& particle_size,
    const VectorDim& deformation_gradient) const {
  // Jacobian dx_i/dxi_j
  return this->jacobian(xi, nodal_coordinates, particle_size,
                        deformation_gradient);
}

template <unsigned Tdim>
inline Eigen::MatrixXd mpm::Element<Tdim>::dn_dx(
    const VectorDim& xi, const Eigen::MatrixXd& nodal_coordinates,
    const VectorDim& particle_size,
    const VectorDim& deformation_gradient) const {
  // Get gradient shape functions
  Eigen::MatrixXd grad_sf =
      this->grad_shapefn(xi, particle_size, deformation_gradient);

  // Jacobian dx_i/dxi_j
  Eigen::Matrix<double, Tdim, Tdim> jacobian =
      (grad_sf.transpose() * nodal_coordinates);

  // Gradient shapefn of the cell
  // dN/dx = [J]^-1 * dN/dxi
  return grad_sf * (jacobian.inverse()).transpose();
}
