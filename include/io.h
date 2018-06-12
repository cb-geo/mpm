#ifndef MPM_IO_H_
#define MPM_IO_H_

#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <utility>

#include <boost/filesystem.hpp>
#include <eigen3/Eigen/Dense>

//! Alias for JSON
#include "json.hpp"
using json = nlohmann::json;

//! \brief Input/Output handler
class IO {
 public:
  //! Constructor
  IO(json json_filename);

  //! Return input mesh file name
  std::string inputMeshFileName() const { return inputMeshFileName_; }

  //! Return input sub mesh file name
  std::string inputSubMeshFileName() const { return inputSubMeshFileName_; }

  //! Return input constraint file name
  std::string constraintsFileName() const { return constraintsFileName_; }

  //! Return input soil particle file name
  std::string inputSoilParticleFileName() const { return inputSoilParticleFileName_; }

  //! Return initial stress soil particle file name
  std::string initStressSoilPFileName() const { return initStressSoilPFileName_; }

  //! Return material file name
  std::string materialFileName() const { return materialFileName_; }

  //! Return traction soil particle file name
  std::string tractionsSoilPFileName() const { return tractionsSoilPFileName_; }

  //! Return output file name
  std::string outputFileName() const { return outputFileName_; }  

  //! Return gravity flag
  bool gravityFlag() const { return gravityFlag_; }

  //! Return boundary friction miu
  double boundaryMiu() const { return boundaryMiu_; }  

  //! Return soil particle spacing
  double soilParticleSpacing() const { return soilParticleSpacing_; }

  //! Return time step or interval
  double dt() const { return dt_; }

  //! Return number of steps of the problem
  unsigned numberOfSteps() const { return numberOfSteps_; }

  //! Return number of substeps of the problem
  unsigned numberOfSubStepsOS() const { return numberOfSubStepsOS_; }

  //! Return newmark flag
  bool newmarkMethod() const { return newmarkMethod_; }

  //! Return gamma - newmark integration coefficient
  double gamma() const { return gamma_; }

  //! Return beta - newmark integration coefficient
  double beta() const { return beta_; }  

  //! Return Cundall Damping flag
  bool dampingFlag() const { return dampingFlag_; }

  //! Return damping ratio for Cundall Damping
  double dampingRatio() const { return dampingRatio_; }

 private:
  //! Input directory
  std::string working_dir_;

  //! Input json object
  json json_;

  //! Material properties json object
  json json_material_properties_;

  //! Input mesh file name
  std::string inputMeshFileName_;

  //! Input submesh file name
  std::string inputSubMeshFileName_;

  //! Input mesh constraints file name
  std::string constraintsFileName_;

  //! Input soil particle file name
  std::string inputSoilParticleFileName_;

  //! Input initial stress soil particle file name
  std::string initStressSoilPFileName_;

  //! Input material file name
  std::string materialFileName_;

  //! Input traction soil particle file name
  std::string tractionsSoilPFileName_;

  //! Input output file name
  std::string outputFileName_;

  //! Input gravity flag
  bool gravityFlag_;

  //! Input boundary friction miu
  double boundaryMiu_;

  //! Input soil particle spacing
  double soilParticleSpacing_;

  //! Input time step or interval
  double dt_;

  //! Input number of steps of the problem
  unsigned numberOfSteps_;

  //! Input number of substeps where it saves
  unsigned numberOfSubStepsOS_;

  //! Input newmark flag
  bool newmarkMethod_;

  //! Input gamma - newmark integration coefficient
  double gamma_;

  //! Input beta - newmark integration coefficient
  double beta_; 

  //! Input Cundall Damping flag
  bool dampingFlag_;

  //! Input damping ratio for Cundall Damping
  double dampingRatio_;

};

#endif  // MPM_IO_H_
