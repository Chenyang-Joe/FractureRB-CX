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

#ifndef collointegrator3d_tpl
#define collointegrator3d_tpl



// Main integration routine for all 'Collocation-type' matrices delegates the
// job to the specialized integration routines.
template<typename KERNEL>
template<typename GDOF,
         typename SUPERELEMENT,
         typename KERNEL_ADAPTOR,
         typename RESULT> 
void ColloIntegrator3d<KERNEL>::
integrate(const GDOF *const  gdof_x,
          const SUPERELEMENT *const se_y,
          const KERNEL_ADAPTOR& call_kernel, 
          RESULT& result) const
{
  HYENA_STATIC_ASSERT( KERNEL::approx_x == CONSTANT ||
                       KERNEL::approx_x == LINEAR ||
                       KERNEL::approx_x == QUADRATIC,
                       NOT_YET_IMPLEMENTED );

  enum{num_nodes_y       = SUPERELEMENT::element_type::num_nodes,
       num_field_nodes_y = SUPERELEMENT::field_type::num_nodes,
       num_topo_nodes    = SUPERELEMENT::element_type::num_topo_nodes};
	
  const ELEMENT_SHAPE shape_y = SUPERELEMENT::element_type::shape;
  const APPROXIMATION geom_y  = SUPERELEMENT::element_type::order;

  const APPROXIMATION field_y = SUPERELEMENT::field_type::order;

  // allocate memory, this has to happen here and temporarily, due to
  // thread-safety
  unsigned int mapped_pos_y[num_nodes_y];
  unsigned int reordered_mapped_pos_y[num_field_nodes_y];
 	
  unsigned int first_y = 0; 
  
  // default singularity type is REGULAR, this is changed below if neccessary
  SING_INT singularity_type = REGULAR;
 
  switch(KERNEL::space_x)
  {
  case DISCONTINUOUS:
    //////////// DISCONTINUOUS distribution of CollocationPoint 's ////////
    // the singularity check depends a great deal on the space type, if the
    // CollocationPoint 's are distributed in a DISCONTINUOUS pattern, the
    // only possible singular type is COINCIDENT
    if(gdof_x->getSuperElement(0)->getID() == se_y->getID())
      singularity_type = COINCIDENT;
    break;
  case CONTINUOUS:
    ////////////// CONTINUOUS distribution of CollocationPoint 's ///////
    // on the other side, for CollocationPoint 's distributed in a
    // CONTINUOUS pattern the only possible singularity types are VRTX -
    // and EDGE_ADJACENT
    const unsigned int num_supp = gdof_x->getNumberSupportElements();
				
    switch(num_supp)
    {
    case 1: 
      if( gdof_x->getSuperElement(0)->getID()
          == se_y->getID() ){
        // if the ReferenceElementIdx is smaller than the number of topo
        // nodes the collocation point has to be on a vertex 
        if( gdof_x->getReferenceElementIdx(0) < num_topo_nodes ){
          singularity_type = VRTX_ADJACENT;
          first_y = gdof_x->getReferenceElementIdx(0);
        }
        // if the ReferenceElementIdx is geater or equal to number of topo
        // nodes and smaller than twice the number of topo nodes, the
        // collocation point has to be on a vertex
        else if( gdof_x->getReferenceElementIdx(0) >= 
                 num_topo_nodes &&
                 gdof_x->getReferenceElementIdx(0) < 
                 2*num_topo_nodes){
          singularity_type = EDGE_ADJACENT;
          first_y = 
            gdof_x->getReferenceElementIdx(0) - num_topo_nodes;
        }
        // if the ReferenceElementIdx is equal to twice the number of topo
        // nodes the CollocationPoint lies in the inside (QUADRATIC
        // QUADRANGLE)
        else if( gdof_x->getReferenceElementIdx(0) ==
                 2*num_topo_nodes)
          singularity_type = COINCIDENT;			
        // the ReferenceElementIdx should not be greater than twice the
        // number of topo nodes
        else HYENA_ERROR_MSG("singularity_check failed");
      }
      break;
			
    case 2: 
      for( unsigned int i=0; i<2; ++i ){
        if( gdof_x->getSuperElement(i)->getID()
            == se_y->getID()){
          if( gdof_x->getReferenceElementIdx(i) < num_topo_nodes ){
            singularity_type = VRTX_ADJACENT;
            first_y = gdof_x->getReferenceElementIdx(i);
          }
          else{
            singularity_type = EDGE_ADJACENT;
            first_y = 
              gdof_x->getReferenceElementIdx(i)-num_topo_nodes;
          }
        }
      }
      break;
		
    default: 
      // In theory there could be any number of supporting elements, but
      // everything more than 20 probably does not make sense. In that
      // case an error message is thrown.
//      if( num_supp < 20 ){
      if( true ){ // DH 06.11.2013: for suspended meshes, this number can be really high, try your best
        for( unsigned int i=0; i<num_supp; ++i )
          if( gdof_x->getSuperElement(i)->getID()
              ==se_y->getID()){
            singularity_type = VRTX_ADJACENT;
            first_y = gdof_x->getReferenceElementIdx(i);
          }
      }
      else HYENA_ERROR_MSG("singularity_check failed: check your mesh!!");
      break;
    } // switch
    break;
  }

  switch(singularity_type)
  {
  case COINCIDENT:  
    // check if indentation of CollocationPoint 's in ColloDuffy is the
    // same as in the CoPt_Iterator
    HYENA_ASSERT( duffy_co.getIndent() == gdof_x->getIndent() );
    integrateSingular(duffy_co, 
                      gdof_x, 
                      se_y, 
                      ApproxTraits<shape_y, geom_y>::mapped_pos,
                      call_kernel,
                      result);
    break;
	
  case EDGE_ADJACENT: 
    // check if indentation of CollocationPoint 's in ColloDuffy is the
    // same as in the CoPt_Iterator
    HYENA_ASSERT( duffy_ea.getIndent() == gdof_x->getIndent() );
    ElementOrientation<shape_y, geom_y>::order(first_y, mapped_pos_y);
    ElementOrientation<shape_y, field_y>::reOrder(first_y, 
                                                    reordered_mapped_pos_y);
    integrateSingular(duffy_ea, 
                      gdof_x, 
                      se_y, 
                      mapped_pos_y,
                      call_kernel,
                      result);	
    reOrderBlock(reordered_mapped_pos_y, result);
    break;
	
  case VRTX_ADJACENT:  
    // check if indentation of CollocationPoint 's in ColloDuffy is the
    // same as in the CoPt_Iterator
    HYENA_ASSERT( duffy_va.getIndent() == gdof_x->getIndent() );
    setOrientation<shape_y, geom_y>(first_y, mapped_pos_y);
    ElementOrientation<shape_y, geom_y>::order(first_y, mapped_pos_y);
    ElementOrientation<shape_y, field_y>::reOrder(first_y, 
                                                  reordered_mapped_pos_y);
    integrateSingular(duffy_va, 
                      gdof_x, 
                      se_y, 
                      mapped_pos_y, 
                      call_kernel,
                      result);	
    reOrderBlock(reordered_mapped_pos_y, result);
    break;

  case REGULAR:
    integrateRegular(gdof_x, 
                     se_y,
                     ApproxTraits<shape_y, geom_y>::mapped_pos, 
                     call_kernel,
                     result);	
    break;
			
  default: HYENA_ERROR_MSG("unknown SINGULARITY TYPE");
  }
}	







// Integration over two distant SuperElement s, that don't share any common
// point. The accuracy (@p _order) needed for the numerical integration is
// estimated trough the distance between the two SuperElement s.
template<typename KERNEL>
template<typename GDOF,
         typename SUPERELEMENT,
         typename KERNEL_ADAPTOR,
         typename RESULT> 
void ColloIntegrator3d<KERNEL>::
integrateRegular(const GDOF *const gdof_x,
                 const SUPERELEMENT *const se_y,
                 const unsigned int *mapped_pos_y,
                 const KERNEL_ADAPTOR& call_kernel, 
                 RESULT& result) const
{
  const double mesh_size = se_y->getMeshSize();

  double distance_x2y = fabs( (gdof_x->getColloPoint() - 
                               se_y->getMidPoint() ).norm() -
                              (0.5*mesh_size) );
	
  distance_x2y /= mesh_size;
	
  unsigned int order;
	
  // order chosen based on the proximity of the two superelements
  // if     (distance_x2y >= 0. && distance_x2y <= 1.) order = 6;
  // else if(distance_x2y >  1. && distance_x2y <= 5.) order = 4;
  // else                                              order = 2;

  if     (distance_x2y >= 0. && distance_x2y <= 3.) order = 14;
  else if(distance_x2y >  3. && distance_x2y <= 5.) order =  9;
  else                                              order =  4;


  // inner loop goes over all quadrature points
  for(unsigned int k=0; k<reg_quad.getNumPoints(order); ++k)
    call_kernel(reg_quad.getPoint(order,k),
                gdof_x->getColloPoint(),
                se_y,
                mapped_pos_y,
                reg_quad.getWeight(order,k),
                result);
}









// Singular integration over two SuperElement s sharing a common
// edge. ColloDuffy coordinate transformation is used to perform this
// geometrically weakly singular integration by a numerical scheme. 
template<typename KERNEL>
template<typename DUFFY,
         typename GDOF,
         typename SUPERELEMENT,
         typename KERNEL_ADAPTOR,
         typename RESULT> 
void ColloIntegrator3d<KERNEL>::
integrateSingular(const DUFFY& duffy,
                  const GDOF *const gdof_x,
                  const SUPERELEMENT *const se_y,
                  const unsigned int *mapped_pos_y,
                  const KERNEL_ADAPTOR& call_kernel, 
                  RESULT& result) const
{
  // if CollocationPoint 's are distributed discontinuously edge/vrtx_adjacent
  // do not exist, hence this routine should not be called ever
  HYENA_ASSERT(not(KERNEL::space_x==DISCONTINUOUS && 
                   DUFFY::sing==EDGE_ADJACENT));
  HYENA_ASSERT(not(KERNEL::space_x==DISCONTINUOUS && 
                   DUFFY::sing==VRTX_ADJACENT));
  
  const unsigned int num_duffy_regions = duffy.getNumRegions();

  // find out the which support of gdof_x the actual se_y
  // is
  unsigned int actual_support = 0;
  for(unsigned int s=0; s<gdof_x->getNumberSupportElements(); ++s)
    if(gdof_x->getSuperElement(s)->getID() == se_y->getID())
      actual_support = s;

  // find out the position where to get duffy coordinates wrt to the actual
  // CollocationPoint
  const unsigned int pos = 
    ColloDuffyTraits<KERNEL::shape,DUFFY::singularity,KERNEL::approx_x>::
    getColloPos(gdof_x->getReferenceElementIdx(actual_support));

  // map pos to the reference position defined by ColloDuffy
  const unsigned int duffy_pos = 
    ApproxTraits<KERNEL::shape, KERNEL::approx_x>::duffy_position[pos];

  // inner loop goes over all quadrature points
  for(unsigned int i=0; i<duffy.getNumQuadPoints(); ++i)
    for(unsigned int k=0; k<num_duffy_regions; ++k){
      call_kernel(duffy.getPointY(i,duffy_pos,k),
                  gdof_x->getColloPoint(),
                  se_y,
                  mapped_pos_y,
                  duffy.getWeight(i,duffy_pos,k)*duffy.getQuadWeight(i),
                  result);
    }
}



// ---------------------------- static functions -------------------------------
/**
 * responsible for the reordering of of the @p block due to the regular
 * element constellation in the @p VRTX_ADJACENT and @p EDGE_ADJACENT case
 */

template<typename BT>
static void reOrderBlock(const unsigned int* reordered_mapped_pos_y, 
                         BT& block)
{ 
  BT orig_block = block;
  for(unsigned int i = 0; i < BT::num_cols; i++)
    block(0,i) = orig_block(0, reordered_mapped_pos_y[i]);
}

template<typename BT>
static void reOrderBlock(const unsigned int* reordered_mapped_pos_y, 
                         std::vector<BT>& blocks)
{
  const unsigned int size = blocks.size();
  for (unsigned int i=0; i<size; ++i) 
    reOrderBlock(reordered_mapped_pos_y, blocks[i]);
}


#endif // collointegrator3d_tpl
