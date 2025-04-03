// Copyright (C) 2009-2010 Matthias Messner, Michael Messner, Franz
// Rammerstorfer, Peter Urthaler
//
// This file is part of HyENA - a C++ boundary element methods library.
//
// HyENA is free software: you can redistribute it and/or modify it under the
// terms of the GNU Lesser Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// HyENA is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser Public License for more
// details.
//
// You should have received a copy of the GNU Lesser Public License along with
// HyENA. If not, see <http://www.gnu.org/licenses/>.

/**
 * @file   extractor.hpp
 * @author TT, Bernhard
 * @date   created:     01.06.11
 *         last change: 01.06.11
 */
#ifndef extractor_hpp
#define extractor_hpp

#include "hyena/core/dofs/dofhandleraccessor.hpp"

namespace hyena
{

  template<ELEMENT_SHAPE E,
           APPROXIMATION G,
           APPROXIMATION F,
           PROBLEM P,
           SPACE_TYPE S>
  class Extractor
    : public DofHandlerAccessor<E,G,F,P,S>
  {
    typedef DofTraits<E,G,F,P,S>                                    dt_type;
    typedef typename dt_type::dof_handler_type             dof_handler_type;
    typedef typename dt_type::ldof_array_type               ldof_array_type;
    typedef typename dt_type::const_ldof_array_type   const_ldof_array_type;

  public:
    // ctor
    Extractor(const dof_handler_type& dof_handler)
      : DofHandlerAccessor<E,G,F,P,S>(dof_handler)
    {}
    
    // fill array with const ldofs with type ldof_type
    void operator() (const_ldof_array_type& ldofs, 
                     const LDOF_TYPE ldof_type) const
    {
      typename ldof_array_type::const_iterator ldof_it =
        this->getLDofs().begin();

      for(; ldof_it != this->getLDofs().end(); ++ldof_it){
		if ((*ldof_it)->getType() == ldof_type or ldof_type==ANY_LDOF_TYPE)
			ldofs.push_back(*ldof_it);
      }
    }

  };



  template<ELEMENT_SHAPE E,
           APPROXIMATION G,
           APPROXIMATION F,
           PROBLEM P>
  class ContourExtractor<E,G,F,P,CONTINUOUS>
    : public DofHandlerAccessor<E,G,F,P,CONTINUOUS>
  {
    typedef DofTraits<E,G,F,P,CONTINUOUS>                     dof_traits_type;
    typedef typename dof_traits_type::dof_handler_type       dof_handler_type;
    typedef typename dof_traits_type::const_ldof_array_type 
    const_ldof_array_type;
  
    static const unsigned int ngdpe = ApproxTraits<E,F>::num_gdofs_per_element;
 
  public:
  
    // ctor
    ContourExtractor(dof_handler_type& dummy)
      : DofHandlerAccessor<E,G,F,P,CONTINUOUS>(dummy)
    {}
  
    void operator()(const const_ldof_array_type& ldofs, 
                    const_ldof_array_type& contour_ldofs, 
                    const LDOF_TYPE ldof_type) const
    {
      typename const_ldof_array_type::const_iterator ldof_it = ldofs.begin();
    
      for(; ldof_it !=ldofs.end(); ++ldof_it){
        bool check = false;
        const unsigned nsupp=(*ldof_it)->getGDof()->getNumberSupportElements();
        for(unsigned int supp = 0; supp < nsupp; ++supp){
          if(check) continue;
          for(unsigned int ngd = 0; ngd < ngdpe; ++ngd){
            if(check) continue;
            else if( (*ldof_it)->getGDof()->getSuperElement(supp)->getGDof(ngd)
                     ->getLDof((*ldof_it)->getLIDX())->getType() != ldof_type ) 
              check = true;
          }
        }
        if(check) contour_ldofs.push_back(*ldof_it);
      }
    }
  };





  template<ELEMENT_SHAPE E,
           APPROXIMATION G,
           APPROXIMATION F,
           PROBLEM P>
  class ContourExtractor<E,G,F,P,DISCONTINUOUS>
    : public DofHandlerAccessor<E,G,F,P,DISCONTINUOUS>
  {
    enum {nld = ProblemTraits<P>::num_ldofs_per_gdof,
          ngd = ApproxTraits<E, F>::num_gdofs_per_element};
  
    typedef DofTraits<E,G,F,P,DISCONTINUOUS>                   dof_traits_type;
    typedef typename dof_traits_type::dof_handler_type        dof_handler_type;
    typedef typename dof_traits_type::const_ldof_array_type 
    const_ldof_array_type;
    typedef typename dof_traits_type::superelement_type                se_type;
    typedef typename dof_traits_type::superelement_array_type    se_array_type;


    typedef typename std::vector<LDOF_TYPE>           ltype_array_type;
    typedef typename std::vector<ltype_array_type>    gtype_array_type;


    std::map<unsigned int, gtype_array_type > ldof_type_map_;

  public:
  
    ContourExtractor(dof_handler_type& dof_handler)
      : DofHandlerAccessor<E,G,F,P,DISCONTINUOUS>(dof_handler)
    {
      typename se_array_type::const_iterator iter = 
        this->getSuperElements().begin();
      for(;iter!=this->getSuperElements().end(); ++iter) {
        gtype_array_type g_array;
        for (unsigned int gd=0; gd<ngd; ++gd) {
          ltype_array_type l_array;
          for (unsigned int ld=0; ld<nld; ++ld)
            l_array.push_back((*iter)->getGDof(gd)
                              ->getLDof(ld)->getType());
          g_array.push_back(l_array);
        }
        ldof_type_map_
          .insert(make_pair((*iter)->getElement()->getId(),g_array));
      }
    }


    /**
     * extracts @p contour_ldofs from @p ldofs if one of the neighboring ldofs
     * is of @p ldof_type
     */
    void operator()(const const_ldof_array_type& ldofs, 
                    const_ldof_array_type& contour_ldofs,
                    const LDOF_TYPE ldof_type) const
    {
      const unsigned int ntn = ShapeTraits<E>::num_topo_nodes;
      unsigned int other_id, ld;
 
      typename const_ldof_array_type::const_iterator ldof_it = ldofs.begin();
      for(; ldof_it !=ldofs.end(); ++ldof_it){
        bool check=false;
        const se_type *const self = (*ldof_it)->getGDof()->getSuperElement(0);
        for(unsigned int n=0; n<ntn; ++n){
          if(check) continue;
          const Node<E,G> *const node = self->getElement()->getNode(n);
          const unsigned int num_elem = node->getNumElements();
          for(unsigned int e=0; e<num_elem; ++e){
            if(check) continue;
            other_id = node->getElement(e)->getId();
            ld = (*ldof_it)->getLIDX();
            for(unsigned int gd=0; gd<ngd;++gd){
              if(check) continue;
              if(ldof_type_map_.find(other_id)->second[gd][ld] == ldof_type) 
                check = true;
            }
          }
        }
        if(check) contour_ldofs.push_back(*ldof_it);
      }
    }
  };





  template <ELEMENT_SHAPE E,
            APPROXIMATION G,
            APPROXIMATION F,
            PROBLEM P,
            SPACE_TYPE S>
  class SuperelementExtractor
  {
    typedef DofTraits<E,G,F,P,S>                               dof_traits_type;
    typedef typename dof_traits_type::gdof_type                      gdof_type;
    typedef typename dof_traits_type::ldof_type                      ldof_type;
    typedef typename dof_traits_type::superelement_type      superelement_type;
    typedef typename dof_traits_type::const_ldof_array_type 
    const_ldof_array_type;
    typedef typename dof_traits_type::const_superelement_array_type
    const_superelement_array_type;
  
  public:
    void operator() (const const_ldof_array_type& ldofs,
                     const_superelement_array_type& elements) const
    {
      boost::unordered_set<unsigned int> eids;
      eids.clear();
      elements.clear();
      elements.reserve(ldofs.size());

      typename const_ldof_array_type::const_iterator ldof_it = ldofs.begin();
      for (; ldof_it!=ldofs.end(); ++ldof_it) {
        const gdof_type* gd = (*ldof_it)->getGDof();
        for (unsigned int e=0; e<gd->getNumberSupportElements(); ++e) {
          const superelement_type* el = gd->getSuperElement(e);
          if ((eids.insert(el->getID())).second) elements.push_back(el);
        }
      }
    }

  };





  template <ELEMENT_SHAPE E,
            APPROXIMATION G,
            APPROXIMATION F,
            PROBLEM P,
            SPACE_TYPE S>
  class GDofExtractor
  {
    typedef DofTraits<E,G,F,P,S>                               dof_traits_type;
    typedef typename dof_traits_type::gdof_type                      gdof_type;
    typedef typename dof_traits_type::ldof_type                      ldof_type;
    typedef typename dof_traits_type::const_ldof_array_type 
    const_ldof_array_type;
    typedef typename dof_traits_type::const_gdof_array_type 
    const_gdof_array_type;


  public:

    void operator() (const const_ldof_array_type& ldofs,
                     const_gdof_array_type& gdofs) const
    {
      boost::unordered_set<unsigned int> gdof_ids;
      gdof_ids.clear();
      gdofs.clear();
      gdofs.reserve(ldofs.size());
    
      typename const_ldof_array_type::const_iterator ldof_it = ldofs.begin();
      for (; ldof_it!=ldofs.end(); ++ldof_it) {
        const gdof_type* gd = (*ldof_it)->getGDof();
        if ((gdof_ids.insert(gd->getID())).second)
          gdofs.push_back(gd);
      }
    }

  };

} // end namespace hyena

#endif

