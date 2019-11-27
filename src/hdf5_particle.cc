#include "hdf5_particle.h"
namespace mpm {
namespace hdf5::particle {
const size_t dst_offset[NFIELDS] = {
    HOFFSET(HDF5Particle, id),
    HOFFSET(HDF5Particle, mass),
    HOFFSET(HDF5Particle, volume),
    HOFFSET(HDF5Particle, pressure),
    HOFFSET(HDF5Particle, coord_x),
    HOFFSET(HDF5Particle, coord_y),
    HOFFSET(HDF5Particle, coord_z),
    HOFFSET(HDF5Particle, displacement_x),
    HOFFSET(HDF5Particle, displacement_y),
    HOFFSET(HDF5Particle, displacement_z),
    HOFFSET(HDF5Particle, nsize_x),
    HOFFSET(HDF5Particle, nsize_y),
    HOFFSET(HDF5Particle, nsize_z),
    HOFFSET(HDF5Particle, velocity_x),
    HOFFSET(HDF5Particle, velocity_y),
    HOFFSET(HDF5Particle, velocity_z),
    HOFFSET(HDF5Particle, stress_xx),
    HOFFSET(HDF5Particle, stress_yy),
    HOFFSET(HDF5Particle, stress_zz),
    HOFFSET(HDF5Particle, tau_xy),
    HOFFSET(HDF5Particle, tau_yz),
    HOFFSET(HDF5Particle, tau_xz),
    HOFFSET(HDF5Particle, strain_xx),
    HOFFSET(HDF5Particle, strain_yy),
    HOFFSET(HDF5Particle, strain_zz),
    HOFFSET(HDF5Particle, gamma_xy),
    HOFFSET(HDF5Particle, gamma_yz),
    HOFFSET(HDF5Particle, gamma_xz),
    HOFFSET(HDF5Particle, epsilon_v),
    HOFFSET(HDF5Particle, status),
    HOFFSET(HDF5Particle, cell_id),
};

// Get size of particle
HDF5Particle particle;
const size_t dst_sizes[NFIELDS] = {
    sizeof(unsigned long long),  // id
    sizeof(double),              // mass
    sizeof(double),              // volume
    sizeof(double),              // pressure
    sizeof(double),              // coord_x
    sizeof(double),              // coord_y
    sizeof(double),              // coord_z
    sizeof(double),              // disp_x
    sizeof(double),              // disp_y
    sizeof(double),              // disp_z
    sizeof(double),              // nsize_x
    sizeof(double),              // nsize_y
    sizeof(double),              // nsize_z
    sizeof(double),              // vel_x
    sizeof(double),              // vel_y
    sizeof(double),              // vel_z
    sizeof(double),              // stress_xx
    sizeof(double),              // stress_yy
    sizeof(double),              // stress_zz
    sizeof(double),              // tau_xx
    sizeof(double),              // tau_yy
    sizeof(double),              // tau_zz
    sizeof(double),              // strain_xx
    sizeof(double),              // strain_yy
    sizeof(double),              // strain_zz
    sizeof(double),              // gamma_xy
    sizeof(double),              // gamma_yz
    sizeof(double),              // gamma_zx
    sizeof(double),              // epsv
    sizeof(bool),                // status
    sizeof(unsigned long long)   // cellid
};

// Define particle field information
const char* field_names[NFIELDS] = {"id",
                                    "mass",
                                    "volume",
                                    "pressure",
                                    "coord_x",
                                    "coord_y",
                                    "coord_z",
                                    "displacement_x",
                                    "displacement_y",
                                    "displacement_z",
                                    "nsize_x",
                                    "nsize_y",
                                    "nsize_z",
                                    "velocity_x",
                                    "velocity_y",
                                    "velocity_z",
                                    "stress_xx",
                                    "stress_yy",
                                    "stress_zz",
                                    "tau_xy",
                                    "tau_yz",
                                    "tau_xz",
                                    "strain_xx",
                                    "strain_yy",
                                    "strain_zz",
                                    "gamma_xy",
                                    "gamma_yz",
                                    "gamma_xz",
                                    "epsilon_v",
                                    "status",
                                    "cell_id"};

// Initialize field types
const hid_t field_type[NFIELDS] = {
    H5T_NATIVE_LLONG,  H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
    H5T_NATIVE_DOUBLE, H5T_NATIVE_HBOOL,  H5T_NATIVE_LLONG};
}  // namespace hdf5::particle
}  // namespace mpm
