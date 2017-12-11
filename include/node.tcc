// Constructor with id, coordinates and dof
//! \param[in] id Node id
//! \param[in] coord coordinates of the node
//! \param[in] dof Degrees of freedom
template <unsigned Tdim>
mpm::Node<Tdim>::Node(Index id, const VectorDim& coord, unsigned dof)
    : NodeBase<Tdim>(id, coord), dof_{dof} {}
