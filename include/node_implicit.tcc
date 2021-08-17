//! Update displacement increment at the node
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_displacement_increment(
    const Eigen::VectorXd& displacement_increment, unsigned phase, const unsigned nactive_node) {

  for(unsigned i = 0; i < Tdim; ++i){
      displacement_(i) += displacement_increment(nactive_node * i + active_id_);
  }
}