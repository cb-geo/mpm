#ifndef MPM_INJECTION_H_
#define MPM_INJECTION_H_
namespace mpm {
struct Injection {
  // Number of particles in each direction
  unsigned nparticles_dir;
  // Particle type
  std::string particle_type;
  // Material id
  std::vector<unsigned> material_ids;
  // Cell id
  int cell_set_id;
  // Start
  double start_time{0.};
  // End
  double end_time{std::numeric_limits<double>::max()};
  // Particle velocity
  std::vector<double> velocity;
};
}  // namespace mpm

#endif  // MPM_INJECTION_H_
