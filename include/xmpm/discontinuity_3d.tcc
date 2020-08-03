template <unsigned Tdim>
mpm::Discontinuity_3D<Tdim>::Discontinuity_3D()
{
  numelement_ = 0;
  
  std::string logger =
      "discontinuity" + std::to_string(Tdim) + "d";
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

//! create elements from file
template <unsigned Tdim>
bool mpm::Discontinuity_3D<Tdim>::create_elements(const std::vector<std::vector<mpm::Index>>& elements)
  {

    bool status = true;
    try {
      // Check if elements is empty
      if (elements.empty())
        throw std::runtime_error("List of elements is empty");
      // Iterate over all elements
      for (const auto& points : elements) {
        
        mpm::discontinuous_element element(points);

        elements_.emplace_back(element); // 
      }
    } catch (std::exception& exception) {
      console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
      status = false;
    }
    return status;
  }

  //initialize the center and normal of the elements
  template <unsigned Tdim>
  bool mpm::Discontinuity_3D<Tdim>::initialize_center_normal()
  {
    bool status = true;
    try {
      VectorDim center;
      VectorDim normal;
      Eigen::Matrix<mpm::Index, 3, 1> points;

      for(auto& element : elements_)
      {
        points = element.points();

        //the center of the element
        for(int i=0; i<3; i++)
          center[i] = 1.0/3*(points_[points[0]].coordinates()[i] + points_[points[1]].coordinates()[i] +points_[points[2]].coordinates()[i]);

        element.set_center(center);
       
       //the normal of the element
        normal =  ThreeCross(points_[points[0]].coordinates(),points_[points[1]].coordinates(),points_[points[2]].coordinates());
        double det = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
        normal = 1.0/det * normal;

        element.set_normal(normal);
      }
    } catch (std::exception& exception) {
      console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
      status = false;
    }
    return status;    
  }

//return the cross product of ab X bc
template <unsigned Tdim>
Eigen::Matrix<double, Tdim, 1> mpm::Discontinuity_3D<Tdim>::ThreeCross(const VectorDim& a,const VectorDim& b,const VectorDim& c){
    
    VectorDim threecross;
    threecross[0]=(b[1]-a[1])*(c[2]-b[2])-(b[2]-a[2])*(c[1]-b[1]);
    threecross[1]=(b[2]-a[2])*(c[0]-b[0])-(b[0]-a[0])*(c[2]-b[2]);
    threecross[2]=(b[0]-a[0])*(c[1]-b[1])-(b[1]-a[1])*(c[0]-b[0]);
    return threecross;
}

//return the levelset values of each doordinates
//! \param[in] the vector of the coordinates
template <unsigned Tdim>
void  mpm::Discontinuity_3D<Tdim>::compute_levelset(const std::vector<VectorDim>& coordinates,std::vector<double>& phi_list){

  mpm::Index i = 0;
  for(const auto& coor:coordinates)
  {
    //find the nearest distance from particle to cell: need to do by global searching and local searching
    double distance = std::numeric_limits<double>::max();
    for(const auto& element:elements_)
    {
      double Vertical_distance_ = element.Vertical_distance(coor);// Vertical_distance(coor);
      distance = std::abs(distance) < std::abs(Vertical_distance_)? distance:Vertical_distance_;
      if(!distance) distance = 1e-16;
    }
    
    phi_list[i] = distance;
    ++i;
  }
}