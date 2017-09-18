/*
  Copyright (C) 2017 by the authors of the ASPECT and CB-Geo MPM code.

  This file is part of ASPECT and CB-Geo MPM.

  MPM is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  MPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with MPM; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include <mpm/particle/particle_iterator.h>

namespace mpm
{
  namespace Particle
  {
    template <int dim, int spacedim>
    ParticleIterator<dim,spacedim>::ParticleIterator (const std::multimap<types::LevelInd, Particle<dim,spacedim> > &map,
                                                      const typename std::multimap<types::LevelInd, Particle<dim,spacedim> >::iterator &particle)
      :
      accessor (map, particle)
    {}



    template <int dim, int spacedim>
    ParticleAccessor<dim,spacedim> &
    ParticleIterator<dim,spacedim>::operator *()
    {
      return accessor;
    }



    template <int dim, int spacedim>
    ParticleAccessor<dim,spacedim> *
    ParticleIterator<dim,spacedim>::operator ->()
    {
      return &(this->operator* ());
    }



    template <int dim, int spacedim>
    const ParticleAccessor<dim,spacedim> &
    ParticleIterator<dim,spacedim>::operator *() const
    {
      return accessor;
    }



    template <int dim, int spacedim>
    const ParticleAccessor<dim,spacedim> *
    ParticleIterator<dim,spacedim>::operator ->() const
    {
      return &(this->operator* ());
    }



    template <int dim, int spacedim>
    bool
    ParticleIterator<dim,spacedim>::operator != (const ParticleIterator<dim,spacedim> &other) const
    {
      return accessor != other.accessor;
    }



    template <int dim, int spacedim>
    bool
    ParticleIterator<dim,spacedim>::operator == (const ParticleIterator<dim,spacedim> &other) const
    {
      return accessor == other.accessor;
    }



    template <int dim, int spacedim>
    ParticleIterator<dim,spacedim> &
    ParticleIterator<dim,spacedim>::operator++()
    {
      accessor.next();
      return *this;
    }



    template <int dim, int spacedim>
    ParticleIterator<dim,spacedim>
    ParticleIterator<dim,spacedim>::operator++(int)
    {
      ParticleIterator tmp(*this);
      operator++ ();

      return tmp;
    }



    template <int dim, int spacedim>
    ParticleIterator<dim,spacedim> &
    ParticleIterator<dim,spacedim>::operator--()
    {
      accessor.prev();
      return *this;
    }



    template <int dim, int spacedim>
    ParticleIterator<dim,spacedim>
    ParticleIterator<dim,spacedim>::operator--(int)
    {
      ParticleIterator tmp(*this);
      operator-- ();

      return tmp;
    }
  }
}

// explicit instantiation of the functions we implement in this file
namespace mpm
{
  namespace Particle
  {
#define INSTANTIATE(dim) \
  template class ParticleIterator<dim>;

    MPM_INSTANTIATE(INSTANTIATE)
  }
}
