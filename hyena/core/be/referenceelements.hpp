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
 * @file    referenceelements.hpp
 * @author  Mathias, Michael, Rf
 * @date    created:     04.08.09
 *          last change: 10.12.09
 */
#ifndef referenceelements_hpp
#define referenceelements_hpp

#ifdef _MSC_VER	
//#define or ||
#include <iso646.h>
#endif

// system includes
#include <vector>


// own includes
#include "hyena/core/tria/element.hpp"
#include "hyena/core/common/enumerators.H"
#include "hyena/core/common/definitions.hpp"
#include "hyena/core/common/mat.hpp"
#include "hyena/core/common/macros.H"
#include "hyena/core/traits/traits.H"
#
namespace hyena
{

  /**
   * @ingroup be
   *
   * This is @p ReferenceElement class. It describes the ansatz space of its
   * data. It contains the shape functions and its gradient, a geometry mapping
   * function and the location of the collocation points are defined here.
   * Herein the outward normal vector is definde by the numbering, with the
   * right-hand-rule.
   *
   * No objects of @p ReferenceElement are needed, because all members are static.
   *
   * Following Elements are implemented:
   *
   *  LINE in 2-dimensional real space.
   *  numbering:
   *  @code
   *  linear line:
   *   o-----------o -> xi
   *   0           1
   *
   *  quadratic line:
   *   o-----o-----o -> x
   *   0     2     1
   *  @endcode
   *
   *  TRIANGLE in 3-dimensional real space.
   *  numbering:
   *  @code
   *  linear triangle:
   *             2
   *             o
   *  xi2       /|
   *   ^      /  |
   *   |    /    |
   *      /      |
   *    /        |
   *   o---------o -> xi1
   *   0         1
   *
   *  quadratic triangle:
   *             2
   *             o
   *  xi2       /|
   *   ^    5 /  |
   *   |    o    o 4
   *      /      |
   *    /        |
   *   o----o----o -> xi1
   *   0    3    1
   *  @endcode
   *
   *  QUADRANGLE in 3-dimensional real space.
   *  numbering:
   *  @code
   *  linear quadrangle:
   *  xi2
   *   ^
   *   |           2
   * 3 o-----------o
   *   |           |
   *   |           |
   *   |           |
   *   |           |
   *   |           |
   *   o-----------o -> xi1
   *   0           1
   *
   *  quadratic quadrangle:
   *  xi2
   *   ^
   *   |     6     2
   * 3 o-----o-----o
   *   |           |
   *   |           |
   * 7 o     o8    o 5
   *   |           |
   *   |           |
   *   o-----o-----o -> xi1
   *   0     4     1
   *  @endcode
   *
   * @tparam SHAPE shape of element
   * @tparam ORDER approximation order of element
   */
  template<ELEMENT_SHAPE SHAPE, APPROXIMATION ORDER>
    class ReferenceElement
  {
  public:

    //! @name compile time constants
    //@{
    static const ELEMENT_SHAPE shape = SHAPE;
    static const APPROXIMATION order = ORDER;

    enum{
      dim                   = ShapeTraits<SHAPE>::dim,
      shape_dim             = ShapeTraits<SHAPE>::shape_dim,
      num_nodes             = ApproxTraits<SHAPE,ORDER>::num_gdofs_per_element,
      num_gdofs_per_element = ApproxTraits<SHAPE,ORDER>::num_gdofs_per_element,
      //gdofs sollen nodes heissen!
      num_collo_pts         = ApproxTraits<SHAPE,ORDER>::num_collo_pts,
      num_edges             = ShapeTraits<SHAPE>::num_edges,
      num_edge_nodes        = ApproxTraits<SHAPE,ORDER>::num_edge_nodes,
      num_inner_nodes       = ApproxTraits<SHAPE,ORDER>::num_inner_nodes,
      num_corner_nodes  = ApproxTraits<SHAPE,ORDER>::num_corner_nodes
    };

    typedef typename PointTraits<shape_dim>::point_type  local_point_type;
    typedef typename PointTraits<      dim>::point_type global_point_type;
    //@} 


    static const unsigned int inner_node_ids[num_inner_nodes];
    static const unsigned int corner_node_ids[num_corner_nodes];
#ifdef _MSC_VER	
//#pragma message ( "Using untested modifications by DH!" )
	static const unsigned int edge_node_ids[num_edges][num_edge_nodes+1];
#else
	static const unsigned int edge_node_ids[num_edges][num_edge_nodes];
#endif
    static const local_point_type local_points_[num_gdofs_per_element];


    /**
     * Evaluate shape functions at a given local point.
     * @param[in] p local point
     * @param[out] shape_fun num_nodes-by-1 vector of evaluated shape functions
     */
    static void evaluateShapeFun(const local_point_type& p,
                                 Mat<double, num_nodes, 1>& shape_fun);



    /**
     * Get local coordinates of the collocation points.  @param[in ] indent
     * value for indenting ColloPoint 's on DISCONTINUOUS boundaries, default
     * value set to zero for CollocationPoint 's on CONTINUOUS boundaries
     */
    static const Mat<local_point_type, num_collo_pts, 1>
    getColloPoints(const double indent = 0.0);



    /**
     * Evaluate gradient of shape functions.
     *\f$ \nabla := [ \frac{\partial}{\partial \xi_1},
     *                \frac{\partial}{\partial \xi_2} ]\f$
     * @param[in] p local point
     * @param[out] grad num_nodes-by-shape_dim matrix
     */
    static void evaluateGradient(const local_point_type& p,
                                 Mat<double, num_nodes, shape_dim>& grad);


    static const unsigned int getInnerNodeId(const unsigned int n) {
      HYENA_ASSERT(n < num_inner_nodes);
      return inner_node_ids[n];
    }

    static const unsigned int getCornerNodeId(const unsigned int n) {
      HYENA_ASSERT(n < num_corner_nodes);
      return corner_node_ids[n];
    }

    static const unsigned int getEdgeNodeId(const unsigned int edge_n,
                                            const unsigned int n) {
      HYENA_ASSERT(n < num_edge_nodes);
      HYENA_ASSERT(edge_n < num_edges);
      return edge_node_ids[edge_n][n];
    }


    static const local_point_type& getLocalPoint(const unsigned int id) {
      HYENA_ASSERT(id < num_gdofs_per_element);
      return local_points_[id];
    }


    /**
     * Get corresponding global point. This is needed in DofHandler.
     * @tparam element_type geometry approximation is needed.
     * @param[in] element pointer to element
     * @param[in] n local id of point
     * @return corresponding global point
     */
    template<typename element_type>
    static const global_point_type
    getGlobalPoint(const element_type *const element,
                   const unsigned int n);

  

    /**
     * Geometrie mapping from local to global point.
     * @param[in] element pointer to element
     * @param[in] lp local point
     * @return corresponding global point 
     */
    static const global_point_type
    mapLocalToGlobal( const Element<shape, order> *const element,
                      const local_point_type& lp )
    {
      Mat<double, num_nodes, 1> shape_fun;
      global_point_type gp(0.);
      evaluateShapeFun(lp,shape_fun);
      for(unsigned int i=0; i<num_nodes; ++i)
        gp += element->getNode(i)->getPoint() * shape_fun[i];
      return gp;
    }



  private:

    //! @name inhibited 
    //@{
    ReferenceElement( );
    ReferenceElement( const ReferenceElement& );
    ReferenceElement& operator = ( const ReferenceElement& );
    ~ReferenceElement();
    //@}
  };



  // At this point all member function definitions are included. They are
  // located in following file
#include "referenceelements.tpl"

} // end namespace hyena

#endif //include guard
 



// /**
//  * Get corresponding local point.
//  * @param[in] n local id of point
//  * @return local point
//  */
// static const local_point_type getLocalPoint(const unsigned int n);
