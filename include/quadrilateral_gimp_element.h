#ifndef MPM_GIMP_ELEMENT_H_
#define MPM_GIMP_ELEMENT_H_

#include <exception>

#include <Eigen/Dense>

#include "logger.h"
#include "quadrilateral_element.h"

namespace mpm {

//!   13----------12----------11----------10
//!   |           |           |           |
//!   |           |           |           |
//!   |           |           |           |
//!   |           |           |           |
//!   |           |           |           |
//! 		        (-1, 1)	    (1,1)
//!   14----------3-----------2-----------9
//!   |           |           |           |
//!   |           |           |           |
//!   |           |   Point   |           |
//!   |           | location  |           |
//!   |           |           |           |
//!   15----------0-----------1-----------8
//!		         (-1,-1)	    (1,-1)
//!   |           |           |           |
//!   |           |           |           |
//!   |           |           |           |
//!   |           |           |           |
//!   |           |           |           |
//!   4-----------5-----------6-----------7
//!
//! </pre>
//! \tparam Tdim Dimension
//! \tparam Tnfunctions Number of functions
template <unsigned Tdim, unsigned Tnfunctions>
class QuadrilateralGIMPElement : public QuadrilateralElement<2, 4> {

 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! constructor with number of shape functions
  QuadrilateralGIMPElement() : QuadrilateralElement<2, 4>() {
    static_assert(Tdim == 2, "Invalid dimension for a GIMP element");
    static_assert((Tnfunctions == 16),
                  "Specified number of shape functions is not defined");

    //! Logger
    std::string logger = "Gimp::<" + std::to_string(Tdim) + ", " +
                         std::to_string(Tnfunctions) + ">";
    console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
  }

  //! Evaluate shape functions at given local coordinates
  //! \param[in] xi given local coordinates
  //! \param[in] particle_size Particle size
  //! \param[in] deformation_gradient Deformation gradient
  //! \retval shapefn Shape function of a given cell
  Eigen::VectorXd shapefn(const VectorDim& xi,
                          const unsigned& number_of_particles,
                          const VectorDim& deformation_gradient) const override;
  //! Evaluate gradient of shape functions
  //! \param[in] xi given local coordinates
  //! \param[in] particle_size Particle size
  //! \param[in] deformation_gradient Deformation gradient
  //! \retval grad_shapefn Gradient of shape function of a given cell
  Eigen::MatrixXd grad_shapefn(
      const VectorDim& xi, const unsigned& number_of_particles,
      const VectorDim& deformation_gradient) const override;

  //! Return the type of shape function
  mpm::ShapefnType shapefn_type() const { return mpm::ShapefnType::GIMP; }

 private:
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
};

}  // namespace mpm
#include "quadrilateral_gimp_element.tcc"

#endif  // MPM_QUADRILATERAL_ELEMENT_H_
