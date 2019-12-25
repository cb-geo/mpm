#ifndef MPM_GIMP_HEX_ELEMENT_H_
#define MPM_GIMP_HEX_ELEMENT_H_

#include "hexahedron_element.h"

namespace mpm {

//! Hexahedron GIMP element class derived from Hexahedron
//! \brief Hexahedron GIMP element
//! \details 64-noded Hexahedron GIMP element, see online document for
//! details \n

//! \tparam Tdim Dimension
//! \tparam Tnfunctions Number of functions
template <unsigned Tdim, unsigned Tnfunctions>
class HexahedronGIMPElement : public HexahedronElement<3, 8> {

 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! constructor with number of shape functions
  HexahedronGIMPElement() : HexahedronElement<3, 8>() {
    static_assert(Tdim == 3, "Invalid dimension for a GIMP element");
    static_assert((Tnfunctions == 64),
                  "Specified number of shape functions is not defined");

    //! Logger
    std::string logger = "hex_gimp::<" + std::to_string(Tdim) + ", " +
                         std::to_string(Tnfunctions) + ">";
    console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
  }

  //! Evaluate shape functions at given local coordinates
  //! \param[in] xi given local coordinates
  //! \param[in] particle_size Particle size
  //! \param[in] deformation_gradient Deformation gradient
  //! \retval shapefn Shape function of a given cell
  Eigen::VectorXd shapefn(const VectorDim& xi, const VectorDim& particle_size,
                          const VectorDim& deformation_gradient) const override;

  //! Evaluate gradient of shape functions
  //! \param[in] xi given local coordinates
  //! \param[in] particle_size Particle size
  //! \param[in] deformation_gradient Deformation gradient
  //! \retval grad_shapefn Gradient of shape function of a given cell
  Eigen::MatrixXd grad_shapefn(
      const VectorDim& xi, const VectorDim& particle_size,
      const VectorDim& deformation_gradient) const override;

  //! Evaluate local shape functions at given local coordinates
  //! \param[in] xi given local coordinates
  //! \param[in] particle_size Particle size
  //! \param[in] deformation_gradient Deformation gradient
  //! \retval shapefn Shape function of a given cell
  Eigen::VectorXd shapefn_local(
      const VectorDim& xi, const VectorDim& particle_size,
      const VectorDim& deformation_gradient) const override;
  //! Compute Jacobian
  //! \param[in] xi given local coordinates
  //! \param[in] nodal_coordinates Coordinates of nodes forming the cell
  //! \param[in] particle_size Particle size
  //! \param[in] deformation_gradient Deformation gradient
  //! \retval jacobian Jacobian matrix
  Eigen::Matrix<double, Tdim, Tdim> jacobian(
      const Eigen::Matrix<double, 3, 1>& xi,
      const Eigen::MatrixXd& nodal_coordinates,
      const Eigen::Matrix<double, 3, 1>& particle_size,
      const Eigen::Matrix<double, 3, 1>& deformation_gradient) const override;

  //! Compute Jacobian local
  //! \param[in] xi given local coordinates
  //! \param[in] nodal_coordinates Coordinates of nodes forming the cell
  //! \param[in] particle_size Particle size
  //! \param[in] deformation_gradient Deformation gradient
  //! \retval jacobian Jacobian matrix
  Eigen::Matrix<double, Tdim, Tdim> jacobian_local(
      const Eigen::Matrix<double, 3, 1>& xi,
      const Eigen::MatrixXd& nodal_coordinates,
      const Eigen::Matrix<double, 3, 1>& particle_size,
      const Eigen::Matrix<double, 3, 1>& deformation_gradient) const override;
  //! Evaluate the B matrix at given local coordinates for a real cell
  //! \param[in] xi given local coordinates
  //! \param[in] nodal_coordinates Coordinates of nodes forming the cell
  //! \param[in] particle_size Particle size
  //! \param[in] deformation_gradient Deformation gradient
  //! \retval bmatrix B matrix
  std::vector<Eigen::MatrixXd> bmatrix(
      const VectorDim& xi, const Eigen::MatrixXd& nodal_coordinates,
      const VectorDim& particle_size,
      const VectorDim& deformation_gradient) const override;

  //! Return the type of shape function
  mpm::ShapefnType shapefn_type() const override {
    return mpm::ShapefnType::GIMP;
  }

  //! Return number of shape functions
  unsigned nfunctions() const override { return Tnfunctions; }

 private:
  //! Return natural nodal coordinates
  Eigen::MatrixXd natural_nodal_coordinates() const;

  //! Logger
  std::unique_ptr<spdlog::logger> console_;
};

}  // namespace mpm
#include "hexahedron_gimp_element.tcc"

#endif  // MPM_GIMP_HEX_ELEMENT_H_
