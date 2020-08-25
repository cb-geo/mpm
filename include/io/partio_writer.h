#ifndef PARTIO_WRITER_H_
#define PARTIO_WRITER_H_

#ifdef USE_PARTIO
#include <vector>

#include <Partio.h>

#include "data_types.h"
#include "hdf5_particle.h"

namespace mpm::partio {

//! Write mesh
//! \param[in] filename Mesh VTP file
//! \param[in] particles HDF5 particles
bool write_particles(const std::string& filename,
                     const std::vector<mpm::PODParticle>& particles);

}  // namespace mpm::partio

#endif  // USE_PARTIO
#endif  // PARTIO_WRITER_H_
