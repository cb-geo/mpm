#ifndef MPM_INJECTION_H_
#define MPM_INJECTION_H_
namespace mpm {
struct Injection {
  // Number of particles in each direction
  unsigned nparticles_dir;
  // Particle type
  std::string particle_type;
  // Material id
  unsigned material_id;
  // Cell id
  int cell_set_id;
  // Particle velocity
  std::vector<double> particle_velocity;
};
}  // namespace mpm

#endif  // MPM_INJECTION_H_
