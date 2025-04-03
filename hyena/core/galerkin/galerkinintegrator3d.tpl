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
 * @file    galerkintegrator3d.hpp
 * @author  MaM, TT, Rf
 * @date    created:     16.08.11
 *          last change: 16.08.11
 */
#ifndef galerkinintegrator3d_tpl
#define galerkinintegrator3d_tpl



template<typename KERNEL>
template<typename SUPERELEMENT_X,
         typename SUPERELEMENT_Y,
         typename KERNEL_ADAPTOR,
         typename RESULT> 
void GalerkinIntegrator3d<KERNEL>::
integrate (const SUPERELEMENT_X *const se_x,
           const SUPERELEMENT_Y *const se_y, 
           const KERNEL_ADAPTOR& call_kernel, 
           RESULT& result) const
{
  typedef SUPERELEMENT_X superelement_type_x;
  typedef SUPERELEMENT_Y superelement_type_y;

  enum{num_nodes_x = superelement_type_x::num_nodes,
       num_nodes_y = superelement_type_y::num_nodes,
       num_field_nodes_x = superelement_type_x::field_type::num_nodes,
       num_field_nodes_y = superelement_type_y::field_type::num_nodes};
	
  const ELEMENT_SHAPE shape_x = superelement_type_x::element_type::shape;
  const ELEMENT_SHAPE shape_y = superelement_type_y::element_type::shape;

  const APPROXIMATION geom_x = superelement_type_x::element_type::order;
  const APPROXIMATION geom_y = superelement_type_y::element_type::order;

  const APPROXIMATION field_x = superelement_type_x::field_type::order;
  const APPROXIMATION field_y = superelement_type_y::field_type::order;

  unsigned int mapped_pos_x[num_nodes_x];
  unsigned int mapped_pos_y[num_nodes_y];

  unsigned int reordered_mapped_pos_x[num_field_nodes_x];
  unsigned int reordered_mapped_pos_y[num_field_nodes_y];
  
  // returns the singularity of the element constellation and sets the right
  // node order of both elements with respect to the present element
  // constellation.
  unsigned int first_x  = 0; 
  unsigned int first_y  = 0;
  unsigned int second_x = -1;
  unsigned int second_y = -1;
  unsigned int num_coincident_nodes = 0;


  for (unsigned int a=0; a<ShapeTraits<KERNEL::shape>::num_topo_nodes; ++a)
    for (unsigned int b=0; b<ShapeTraits<KERNEL::shape>::num_topo_nodes; ++b)
      if ( se_x->getNode(a) == se_y->getNode(b) ) {
				
        // increment number of coincident nodes
        num_coincident_nodes++;
				
        // set number of coincident nodes
        if (second_x == -1) 
          first_x = second_x = a;
        else			  second_x = a;
				
        if (second_y == -1)
          first_y = second_y = b;
        else			  second_y = b;
      }


  // integrate depending on singularity  
  switch (num_coincident_nodes)
    {
    case 0: { // regular case
      // node ordering must not be changed
      integrateRegular(se_x, se_y, 
                       ApproxTraits<shape_x, geom_x>::mapped_pos, 
                       ApproxTraits<shape_y, geom_y>::mapped_pos,
                       call_kernel, result);
    }
      break;
    case 1: { // vertex adjacent case
      // node odering has to be performed in order to get a regular element
      // constellation after Sauter.
      ElementOrientation<shape_x, geom_x>::order(first_x, mapped_pos_x);
      ElementOrientation<shape_y, geom_y>::order(first_y, mapped_pos_y);
      integrateSingular(duffy_va, se_x, se_y, mapped_pos_x, mapped_pos_y, 
                         call_kernel, result);
      
      // reordering the block
      ElementOrientation<shape_x, field_x>::
        reOrder(first_x, reordered_mapped_pos_x);
      ElementOrientation<shape_y, field_y>::
        reOrder(first_y, reordered_mapped_pos_y);
      reOrderBlock(reordered_mapped_pos_x, reordered_mapped_pos_y, result);
    }
      break;
    case 2: { // edge adjacent case
      // Choose the first node for the mapping, keeping in mind that the
      // common edge on both elements has the SAME(!) parametrization
      // (i.e. se_x has positive orientation whereas se_y has negative
      // orientation).
      if(second_x - first_x != 1){
        first_x = second_x;
        first_y = second_y;
      }
      // node odering has to be performed in order to get a regular element
      // constellation after Sauter.
      ElementOrientation<shape_x, geom_x>::order(first_x, mapped_pos_x);
      ElementOrientation<shape_y, geom_y>::reverseOrder(first_y, mapped_pos_y);
      integrateSingular(duffy_ea, se_x, se_y, mapped_pos_x, mapped_pos_y, 
                         call_kernel, result);

      // reordering the block
      ElementOrientation<shape_x, field_x>::
        reOrder(first_x, reordered_mapped_pos_x);
      ElementOrientation<shape_y, field_y>::
        reOrderReverse(first_y, reordered_mapped_pos_y);
      reOrderBlock(reordered_mapped_pos_x, reordered_mapped_pos_y, result);
    }
      break;
    case 3:   // coincident case for triangles
    case 4: { // coincident case for quadrangles
      // node ordering must not be changed
      integrateSingular(duffy_co, se_x, se_y, 
                        ApproxTraits<shape_x, geom_x>::mapped_pos, 
                        ApproxTraits<shape_y, geom_y>::mapped_pos,
                        call_kernel, result);
    }
      break;
    default:
      HYENA_ERROR_MSG("Galerkinintegrator3d: number of coincident not ok");
    }
}






template<typename KERNEL>
template<typename SUPERELEMENT_X,
         typename SUPERELEMENT_Y,
         typename KERNEL_ADAPTOR,
         typename RESULT> 
void GalerkinIntegrator3d<KERNEL>::
integrateRegular (const SUPERELEMENT_X *const se_x,
                  const SUPERELEMENT_Y *const se_y,
                  const unsigned int *const mapped_pos_x,
                  const unsigned int *const mapped_pos_y,
                  const KERNEL_ADAPTOR& call_kernel, 
                  RESULT& result) const
{
  // estimate the distance between the two superelements in order to minimize
  // quadrature points used for the numerical integration scheme
  double h_x = se_x->getMeshSize();
  double h_y = se_y->getMeshSize();
  double distance_x2y = fabs( ( se_x->getMidPoint() - 
                                se_y->getMidPoint() ).norm()
                              - .5*( h_x + h_y ) );
  
  double mesh_size = ( h_x > h_y ? h_x : h_y );
  
  distance_x2y /= mesh_size;
  
  // integration order chosen based on the proximity of the two superelements
  unsigned int order;

   if(      distance_x2y >  5.                       )	order = 2;
   else if( distance_x2y >  1. && distance_x2y <= 5. )	order = 3;
   else                                               	order = 4;
  //if(      distance_x2y >  5.                       )	order = 3;
  //else if( distance_x2y >  1. && distance_x2y <= 5. )	order = 6;
  //else                                               	order = 9;

  //int order = 5;

  for (unsigned int i =0; i < reg_quad.getNumPoints(order); ++i)
    for (unsigned int k=0; k < reg_quad.getNumPoints(order); ++k)
      call_kernel(reg_quad.getPoint(order,i),
                  reg_quad.getPoint(order,k),
                  se_x,
                  se_y,
                  mapped_pos_x,
                  mapped_pos_y,
                  reg_quad.getWeight(order,k) * reg_quad.getWeight(order,i),
                  result);
}






template<typename KERNEL>
template<typename DUFFY,
         typename SUPERELEMENT_X,
         typename SUPERELEMENT_Y, 
         typename KERNEL_ADAPTOR,
         typename RESULT> 
void GalerkinIntegrator3d<KERNEL>::
integrateSingular (const DUFFY& duffy,
                   const SUPERELEMENT_X *const se_x,
                   const SUPERELEMENT_Y *const se_y,
                   const unsigned int *const mapped_pos_x,
                   const unsigned int *const mapped_pos_y,
                   const KERNEL_ADAPTOR& call_kernel, 
                   RESULT& result) const
{
  const unsigned int num_duffy_regions = duffy.getNumRegions();
  
  for (unsigned int i=0; i<sing_quad.getNumPoints(); ++i)	{
    for (unsigned int j=0; j<sing_quad.getNumPoints(); ++j) {
      const double weight = sing_quad.getWeight(i) * sing_quad.getWeight(j);
      for (unsigned int k=0; k<num_duffy_regions; ++k)
        call_kernel(duffy.getPointX(i, j, k),
                    duffy.getPointY(i, j, k),
                    se_x,
                    se_y,
                    mapped_pos_x,
                    mapped_pos_y,
                    duffy.getWeight(i, j, k) * weight,
                    result);
    } // loop j
  } // loop i
}



// ---------------------------- static functions -------------------------------

/**
 * responsible for the reordering of of the @p block due to the regular
 * element constellation in the @p VRTX_ADJACENT and @p EDGE_ADJACENT case
 */
template<typename BT>
static void reOrderBlock( const unsigned int* reordered_mapped_pos_x,
                          const unsigned int* reordered_mapped_pos_y,
                          BT& block)
{ 
  BT orig_block = block;
  for(unsigned int i = 0; i < BT::num_rows; i++)
    for(unsigned int j = 0; j < BT::num_cols; j++)
      block(i,j) = orig_block(reordered_mapped_pos_x[i],
                              reordered_mapped_pos_y[j]);
}


template<typename BT>
static void reOrderBlock( const unsigned int* reordered_mapped_pos_x,
                          const unsigned int* reordered_mapped_pos_y,
                          std::vector<BT>& blocks)
{ 
  const unsigned int size = blocks.size();
  for (unsigned int i=0; i<size; ++i) 
    reOrderBlock(reordered_mapped_pos_x, reordered_mapped_pos_y, blocks[i]);
}


#endif

