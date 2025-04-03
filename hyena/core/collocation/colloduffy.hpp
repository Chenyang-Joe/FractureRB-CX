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
 * @file    colloduffy.hpp
 * @author  Michael
 * @date    created:     12.10.09
 *          last change: 29.01.10
 */
#ifndef colloduffy_h
#define colloduffy_h

// own includes
#include "hyena/core/common/enumerators.H"
#include "hyena/core/common/mat.hpp"
#include "hyena/core/common/macros.H"


namespace hyena
{

  // definitions of static arrays in ApproxTraits
  const unsigned int ApproxTraits<TRIANGLE,CONSTANT>::inner_collo_pos[ 1]={0};
#ifdef _MSC_VER	
  const unsigned int ApproxTraits<TRIANGLE,CONSTANT>::edge_collo_pos[  1]={-1};
#else
  const unsigned int ApproxTraits<TRIANGLE,CONSTANT>::edge_collo_pos[  0]={ };
#endif
  const unsigned int ApproxTraits<TRIANGLE,CONSTANT>::duffy_position[  1]={0};

#ifdef _MSC_VER	
  const unsigned int ApproxTraits<TRIANGLE,LINEAR>::inner_collo_pos[   1]={-1};
  const unsigned int ApproxTraits<TRIANGLE,LINEAR>::edge_collo_pos[    1]={-1};
#else
  const unsigned int ApproxTraits<TRIANGLE,LINEAR>::inner_collo_pos[   0]={ };
  const unsigned int ApproxTraits<TRIANGLE,LINEAR>::edge_collo_pos[    0]={ };
#endif
  const unsigned int ApproxTraits<TRIANGLE,LINEAR>::duffy_position[    3]=
    {0, 0, 0};

#ifdef _MSC_VER	
  const unsigned int ApproxTraits<TRIANGLE,QUADRATIC>::inner_collo_pos[1]={-1};
#else
  const unsigned int ApproxTraits<TRIANGLE,QUADRATIC>::inner_collo_pos[0]={ };
#endif
  const unsigned int ApproxTraits<TRIANGLE,QUADRATIC>::edge_collo_pos[ 1]={3};
  const unsigned int ApproxTraits<TRIANGLE,QUADRATIC>::duffy_position[ 6]=
    {0, 0, 0, 1, 1, 1};

  const unsigned int ApproxTraits<QUADRANGLE,CONSTANT>::inner_collo_pos[ 1]={0};
#ifdef _MSC_VER	
  const unsigned int ApproxTraits<QUADRANGLE,CONSTANT>::edge_collo_pos[  1]={-1};
#else
  const unsigned int ApproxTraits<QUADRANGLE,CONSTANT>::edge_collo_pos[  0]={ };
#endif
  const unsigned int ApproxTraits<QUADRANGLE,CONSTANT>::duffy_position[  1]={0};

#ifdef _MSC_VER	
  const unsigned int ApproxTraits<QUADRANGLE,LINEAR>::inner_collo_pos[   1]={-1};
  const unsigned int ApproxTraits<QUADRANGLE,LINEAR>::edge_collo_pos[    1]={-1};
#else
  const unsigned int ApproxTraits<QUADRANGLE,LINEAR>::inner_collo_pos[    ]={ };
  const unsigned int ApproxTraits<QUADRANGLE,LINEAR>::edge_collo_pos[    0]={ };
#endif
  const unsigned int ApproxTraits<QUADRANGLE,LINEAR>::duffy_position[    4]=
    {0, 0, 0, 0};

#ifdef _MSC_VER	
  const unsigned int ApproxTraits<QUADRANGLE,QUADRATIC>::inner_collo_pos[1]={-1};
#else
  const unsigned int ApproxTraits<QUADRANGLE,QUADRATIC>::inner_collo_pos[0]={ };
#endif
  const unsigned int ApproxTraits<QUADRANGLE,QUADRATIC>::edge_collo_pos[ 1]={4};
  const unsigned int ApproxTraits<QUADRANGLE,QUADRATIC>::duffy_position[ 9]=
    {0, 0, 0, 0, 1, 1, 1, 1, 8};







  // ColloDuffyTraits
  template <ELEMENT_SHAPE SHAPE,
            SING_INT SINGULARITY,
            APPROXIMATION ORDER>
  struct ColloDuffyTraits;

  template<ELEMENT_SHAPE SHAPE, APPROXIMATION ORDER>
  struct ColloDuffyTraits<SHAPE, VRTX_ADJACENT, ORDER>
  { enum {num_duffy_collo_pts = 1 };
    static unsigned int getColloPos(unsigned int i)
    { return 0; }
  };

  template<ELEMENT_SHAPE SHAPE, APPROXIMATION ORDER>
  struct ColloDuffyTraits<SHAPE, EDGE_ADJACENT, ORDER>
  {
    enum {num_duffy_collo_pts = ApproxTraits<SHAPE, ORDER>::num_edge_collo_pts};
    static unsigned int getColloPos(unsigned int i)
    { return ApproxTraits<SHAPE, ORDER>::edge_collo_pos[i]; }
  };

  template<ELEMENT_SHAPE SHAPE, APPROXIMATION ORDER>
  struct  ColloDuffyTraits<SHAPE, COINCIDENT, ORDER>
  { 
    enum {num_duffy_collo_pts=ApproxTraits<SHAPE, ORDER>::num_inner_collo_pts};
    static unsigned int getColloPos(unsigned int i)
    { return ApproxTraits<SHAPE, ORDER>::inner_collo_pos[i]; }
  };





  /**
   * @ingroup collocation
   *
   * The BaseColloDuffy object depends on the template parameter @p
   * ELEMENT_SHAPE of the two Elements it treats and @p SING_INT containing the
   * information about their singular (regular, coincident, vertex adjacent,
   * edge adjacent) position to each other. This class is intendet as base class
   * to ColloDuffy, since the transformation rules only have to be specialized
   * for these two template parameters, but the finally required ColloDuffy
   * object depends on more parameters. How the tranformation rules work, is
   * explained best in two steps.
   *
   * -# STEP:
   * The actual "duffy transformation" is used to map QuadraturePoints(*) from
   * a Quadrangle(xi',eta') to a Triangle(xi,eta):
   * \verbatim

   ^ eta'                     ^	eta
   |                          |
   o--------------o           |               o (xi2,eta2)
   |              |           |             / |
   |   *      *   |           |           /   |
   |              |           |         /  *  |
   |              |    =>     |       /       |
   |              |           |     /         |
   |   *      *   |           |   / *      *  |
   |              |           | /   *         |
   o--------------o-> xi'     o---------------o-> xi
   (xi1,eta1)
   \endverbatim
   * -# STEP:
   * Afterwards these QuadraturePoints are mapped onto all subtriangles the
   * ReferenceElement y (corresponding to the original SuperElement y) is
   * divided into.This means that singular elements end up having more
   * quadrature points condensed around the singular point ( lim r->0 ).
   *
   * The mapping for quadrature Points from above to one subtriangle is given
   * by:
   *
   * y1 = x[1] + ( eta1 - x[1]) * y[0]	+ (eta2 - eta1) * y[1];
   *
   * det = fabs( (xi1 - x[0])*(eta2 - eta1) - (eta1 - x[1])*(xi2 - xi1) );
   *
   * - (xi1,eta1) first vertex after the singular point per subtriangle
   * - (xi2,eta2) second vertex after the singular point per subtriangle
   * - x local coordinates of the collocation point on superelement_x
   * - y local coordinates of the quadrature point on superelement_y
   * - y1 local coordinates of the quadrature point on superelement_y
   * - det jacobi determinant lifts the singularity
   *
   * @tparam SHAPE shape of the Element to treat
   * @tparam SINGULARITY relative position of CollocationPoint to the Element
   * which has to be integrated over
   *
   * @sa ColloDuffy
   */
  template<ELEMENT_SHAPE SHAPE, SING_INT SINGULARITY>
  class BaseColloDuffy
  {
    // number of subregions a ReferenceElement is divided into for the duffy
    // transformation of the actual SING_INT
    enum{ num_regions = SingularityTraits<SHAPE,SINGULARITY>::num_col_regions };



    // define the own point_type/weight_type depending on the number of regions
    // a ReferenceElement is divided into for the duffy transformation of the
    // actual SING_INT
    typedef Point2  point_type[num_regions];
    typedef double weight_type[num_regions];


    // the quadrature rule required for singular integration is allways known at
    // compile time, further it is only required by the ColloDuffy object. Since
    // both, ColloDuffy and QuadratureRule, are defined at compile time all
    // possible Quadrature Points for the singular quadrature are defined at
    // compile time too therfore it is sufficient to just compute them
    // once. This explains, why the singular QuadratureRule is only needed here.
    typedef QuadratureRule<QUADRANGLE, GAUSS>                sing_quad_type;











  public:
    static const SING_INT singularity = SINGULARITY;



    //! ctor
    BaseColloDuffy(const unsigned int sing_quad_order)
      : sing_quad(sing_quad_order),
        num_quad_points(sing_quad.getNumPoints()),
        reference_points_y(NULL),
        weights(NULL)
    {
      // allocate memory for all duffy transformed quadrature points
      reference_points_y = new  point_type*[num_quad_points];
      weights            = new weight_type*[num_quad_points];

      for(unsigned int q=0; q<num_quad_points; ++q) {
        reference_points_y[q] = NULL;
        weights[           q] = NULL;
      }
    }






    //! dtor
    ~BaseColloDuffy( )
    {
      for(unsigned int q=0; q<num_quad_points; ++q) {
        delete[] reference_points_y[q];
        delete[] weights[           q];
      }

      delete[] reference_points_y;
      delete[] weights;

    };




    //! @return transformed coordinates
    const Point2& getPointY(const unsigned int i,
                            const unsigned int j,
                            const unsigned int k) const
    {
      HYENA_ASSERT(i<num_quad_points);
      HYENA_ASSERT(j<num_duffy_collo_points);
      HYENA_ASSERT(k<num_regions);
      return reference_points_y[i][j][k];
    }




    //! @return modification of quadrature weights due to duffy
    const double getWeight(const unsigned int i,
                           const unsigned int j,
                           const unsigned int k) const
    {
      HYENA_ASSERT(i<num_quad_points);
      HYENA_ASSERT(j<num_duffy_collo_points);
      HYENA_ASSERT(k<num_regions);
      return weights[i][j][k];
    }




    //! @return number of duffy regions
    const unsigned int getNumRegions() const
    {
      return num_regions;
    }





    //! @return number quadrature points for sing_quad
    const unsigned int getNumQuadPoints() const
    {
      return num_quad_points;
    }




    //! @return number of quadrature points for sing_quad
    const double getQuadWeight(const unsigned int k) const
    {
      return sing_quad.getWeight(k);
    }



  protected:
    /**
     * A ColloDuffy object is not to be copied! Thus the copy constructor is
     * private and not defined.
     */
    BaseColloDuffy(const BaseColloDuffy&);
    BaseColloDuffy();



    /**
     * A ColloDuffy object is not to be assigned! Thus the operator= is
     * private and not defined.
     */
    const BaseColloDuffy& operator= (const BaseColloDuffy&);



    /**
     * Spezialized transformation rule for the COINCIDENT case, i.e. the
     * Collocation point lies on one Element.
     * @param[in] x local coordinate of the CollocationPoint
     * @param[in] y local coordinate of the gauss point on superelement_y
     * @param[out] reference_points_y transformed gauss points
     * @param[out] weights modified gauss weights points
     */
    void transformDuffy(const Point2& x,
                        const Point2& y,
                        point_type& reference_points_y,
                        weight_type& weights);


    /**
     * Only one quadrature rule for singular integration is required due to the
     * fact, that duffy transformation is used for the weakly singular
     * integration. No matter what elements (Triangle or Quadrangle) are used
     * this only needs quadrature points on the reference quadrangle.
     */
    const sing_quad_type sing_quad;

    const unsigned int num_quad_points;

    unsigned int num_duffy_collo_points;

    point_type** reference_points_y;

    weight_type** weights;
  };







  /**
   * @ingroup collocation
   *
   * ColloDuffy contains the final duffy transformation requiered for
   * Collocation boundary element methods. It is derived from BaseColloDuffy and
   * depends on the ssame template parameters @p ELEMENT_SHAPE of the
   * ansatz-support element and @p SING_INT but additionally also on @p
   * APPROXIMATION of the test-support (number of collocation points on the
   * 'test-element') and @p SPACE_TYPE of the test-support (whether the
   * CollocationPoint 's are distributed dis/continuously).
   *
   * @tparam SHAPE shape of the ansatz-support element tau_y 
   * @tparam SINGULARITY relative position of CollocationPoint to tau_y
   * @tparam APPROXIMATION number of CollocationPoint 's on the test-support
   * element tau_x
   * @tparam SPACE_TYPE of test-support element tau_x
   *
   * @sa BaseColloDuffy
   */
  template<ELEMENT_SHAPE SHAPE,
           SING_INT SINGULARITY,
           APPROXIMATION ORDER,
           SPACE_TYPE SPACE>
  class ColloDuffy : public BaseColloDuffy<SHAPE,SINGULARITY>
  {
    enum{ num_regions = SingularityTraits<SHAPE,SINGULARITY>::num_col_regions };



    // define the own point_type/weight_type depending on the number of regions
    // a ReferenceElement is divided into for the duffy transformation of the
    // actual SING_INT
    typedef Point2  point_type[num_regions];
    typedef double weight_type[num_regions];

  public:
    static const SING_INT sing = SINGULARITY;


    ColloDuffy(const unsigned int sing_quad_order)
      : BaseColloDuffy<SHAPE,SINGULARITY>(sing_quad_order),
        indent(SPACE==DISCONTINUOUS ? ShapeTraits<SHAPE>::indent : 0.)
    {
      Mat<Point2,ReferenceElement<SHAPE,ORDER>::num_collo_pts,1> collo_points_x;

      collo_points_x = ReferenceElement<SHAPE,ORDER>::getColloPoints(indent);

      // compute local coordinates of all CollocationPoint 's
      switch(SPACE)
      {
      case CONTINUOUS:
        // On CONTINUOUS boundaries the number of CollocationPoint 's for the
        // current ColloDuffy object is defined by the ColloDuffyTraits
        num_duffy_collo_points =
          ColloDuffyTraits<SHAPE, SINGULARITY, ORDER>::num_duffy_collo_pts;
        break;
      case DISCONTINUOUS:
        // On DISCONTINUOUS boundaries all CollocationPoint 's ly inside a
        // boundary element, hence the number of CollocationPoints is the amount
        // of actual CollocationPoints defined on the ReferenceElement
        num_duffy_collo_points = ApproxTraits<SHAPE, ORDER>::num_collo_pts;
        break;
      default:
        HYENA_ERROR();
      }

      for (unsigned int q=0; q<num_quad_points; ++q) {
        reference_points_y[q] = new  point_type[num_duffy_collo_points];
        weights[           q] = new weight_type[num_duffy_collo_points];
      }



      // perform the duffy transformation for the given SING_INT
      for (unsigned int i=0; i<num_quad_points; ++i)
        for (unsigned int _j=0; _j<num_duffy_collo_points; ++_j) {

          // get the local position of the collocation point on the
          // ReferenceElement
          const unsigned int j =
            ColloDuffyTraits<SHAPE, SINGULARITY, ORDER>::getColloPos(_j);

          this -> transformDuffy( collo_points_x[j],
                                  sing_quad.getPoint(i),
                                  reference_points_y[i][_j],
                                  weights[i][_j] );
        }
    }



    //! @return indent
    const double getIndent() const
    {
      return indent;
    }


  private:
    ColloDuffy();
    ColloDuffy(const ColloDuffy&);
    const ColloDuffy& operator= (const ColloDuffy&);

    using BaseColloDuffy<SHAPE,SINGULARITY>::sing_quad;
    using BaseColloDuffy<SHAPE,SINGULARITY>::num_quad_points;
    using BaseColloDuffy<SHAPE,SINGULARITY>::num_duffy_collo_points;
    using BaseColloDuffy<SHAPE,SINGULARITY>::reference_points_y;
    using BaseColloDuffy<SHAPE,SINGULARITY>::weights;

    const double indent;
  };









  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  ///////////////////////////                    //////////////////////////
  ///////////////////////////   implementation   //////////////////////////
  ///////////////////////////                    //////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  /**
   * Subdivision of the ReferenceElement into three subtriangles for the
   * COINCIDENT transformation (* marks the collocation point).
   *
   * \verbatim
   ^ eta
   |
   |                       o
   |                     //|
   |                   / / |
   |                 /  /  |
   |               /   /   |
   |             /    /    |
   |           /  3  /     |
   |         /      *   2  |
   |       /     _ / \     |
   |     /   _ /       \   |
   |   / _ /      1     \  |
   | /_/                  \|
   o-----------------------o--> xi
   \endverbatim
  */
  template<> inline
  void BaseColloDuffy<TRIANGLE,COINCIDENT>::
  transformDuffy(const Point2& x,
                 const Point2& y,
                 point_type& reference_points_y,
                 weight_type& weights)
  {
    double y0  = y[0];
    double y01 = y[0]*y[1];
    double det = fabs(y0);

    // Triangle 1
    reference_points_y[0][0] = x[0] - x[0]*y0 +  y01;
    reference_points_y[0][1] = x[1] - x[1]*y0;
    weights[0] = det*fabs( x[1] );

    // Triangle 2
    reference_points_y[1][0] = x[0] + ( 1.0 - x[0])*y0;
    reference_points_y[1][1] = x[1] -         x[1] *y0	+  y01;
    weights[1] = det*fabs( 1.0 - x[0] );

    // Triangle 3
    reference_points_y[2][0] = x[0] + ( 1.0 - x[0])*y0 - y01;
    reference_points_y[2][1] = x[1] + ( 1.0 - x[1])*y0 - y01;
    weights[2] = det*fabs( (1.0 - x[1]) - (1.0 - x[0]) );

  }









  /**
   * Subdivision of the ReferenceElement into two subtriangles for the
   * EDGE_ADJACENT transformation (* marks the collocation point).
   *
   * \verbatim
   ^ eta
   |
   |                     o
   |                   //|
   |                 / / |
   |               /  /  |
   |             /   /   |
   |           /    /    |
   |         /     /     |
   |       /   2  /      |
   |     /       /       |
   |   /        /   1    |
   | /         /         |
   o----------*----------o-->	xi
   \endverbatim
  */
  template<> inline
  void BaseColloDuffy<TRIANGLE,EDGE_ADJACENT>::
  transformDuffy(const Point2& x,
                 const Point2& y,
                 point_type& reference_points_y,
                 weight_type& weights)
  {
    double y0  = y[0];
    double y01 = y[0]*y[1];
    double det = fabs(y0);

    // projection of the collocation point on superelement_x on to the bottom
    // edge of superelement_y (only happens for quadratic approximation)
    double x0 = 0.5;
    double x1 = 0.0;

    // Triangle 1
    reference_points_y[0][0] = x0 + (1.0 - x0)*y0;
    reference_points_y[0][1] = x1 -        x1 *y0 + y01;
    weights[0] = det*fabs(1.0 - x0);


    // Triangle 2
    reference_points_y[1][0] = x0 + ( 1.0 - x0)*y0 - y01;
    reference_points_y[1][1] = x1 + ( 1.0 - x1)*y0 - y01;
    weights[1] = det* fabs( (1.0 - x1)-(1.0 - x0) );
  }








  /**
   * Standard Duffy transformation in VRTX_ADJACENT case (* marks the
   * collocation point).
   *
   * \verbatim
   ^ eta
   |
   |                    o
   |                  / |
   |                /   |
   |              /     |
   |            /       |
   |          /         |
   |        /           |
   |      /       1     |
   |    /               |
   |  /                 |
   |/                   |
   *--------------------o--> xi
   \endverbatim
  */
  template<> inline
  void BaseColloDuffy<TRIANGLE,VRTX_ADJACENT>::
  transformDuffy(const Point2& x,
                 const Point2& y,
                 point_type& reference_points_y,
                 weight_type& weights)
  {
    double y0  = y[0];
    double y01 = y[0]*y[1];
    double det = fabs(y0);

    // Triangle 1
    reference_points_y[0][0] = y0;
    reference_points_y[0][1] = y01;
    weights[0]   = det;
  }







  /**
   * Subdivision of the ReferenceElement into four subtriangles for the
   * COINCIDENT transformation (* marks the collocation point).
   *
   * \verbatim
   ^ eta
   |
   o---------------------o
   |\                   /|
   |  \       3       /  |
   |    \           /    |
   |      \       /      |
   |        \   /        |
   |   4      *      2   |
   |        /   \        |
   |      /       \      |
   |    /           \    |
   |  /       1       \  |
   |/                   \|
   o---------------------o--> xi
   \endverbatim
  */
  template<> inline
  void BaseColloDuffy<QUADRANGLE,COINCIDENT>::
  transformDuffy(const Point2& x,
                 const Point2& y,
                 point_type& reference_points_y,
                 weight_type& weights)
  {

    double y0  = y[0];
    double y01 = y[0]*y[1];
    double det = fabs(y0);

    // Triangle 1
    reference_points_y[0][0] = x[0] - x[0]*y0 + y01;
    reference_points_y[0][1] = x[1] - x[1]*y0;
    weights[0] = det*fabs( x[1] );

    // Triangle 2
    reference_points_y[1][0] = x[0] + ( 1.0 - x[0])*y0;
    reference_points_y[1][1] = x[1] - (       x[1])*y0 + y01;
    weights[1]= det*fabs(1.0 - x[0]);

    // Triangle 3
    reference_points_y[2][0] = x[0] + (1.0 - x[0])*y0 - y01;
    reference_points_y[2][1] = x[1] + (1.0 - x[1])*y0;
    weights[2]= det*fabs(1.0 - x[1]);

    // Triangle 4
    reference_points_y[3][0] = x[0] - (      x[0])*y0;
    reference_points_y[3][1] = x[1] + (1.0 - x[1])*y0 - y01;
    weights[3]= det*fabs( x[0] );
  }










  /**
   * Subdivision of the ReferenceElement into four subtriangles for the
   * EDGE_ADJACENT transformation (* marks the collocation point).
   *
   * \verbatim
   ^ eta
   |
   o-------------------o
   |\        |        /|
   | \       |       / |
   |  \   4  |  3   /  |
   |   \     |     /   |
   |    \    |    /    |
   |     \   |   /     |
   |  1   \  |  /   2  |
   |       \ | /       |
   |        \|/        |
   o---------*---------o-->	xi
   \endverbatim
  */
  template<> inline
  void BaseColloDuffy<QUADRANGLE,EDGE_ADJACENT>::
  transformDuffy(const Point2& x,
                 const Point2& y,
                 point_type& reference_points_y,
                 weight_type& weights)
  {
    double y0  = y[0];
    double y01 = y[0]*y[1];
    double det = fabs(y0);

    // projection of the collocation point on superelement_x on to the bottom
    // edge of superelement_y (only happens for quadratic approximation)
    double x0 = 0.5;
    double x1 = 0.0;

    // Triangle 1
    reference_points_y[0][0] = x0 - (      x0)*y0;
    reference_points_y[0][1] = x1 + (1.0 - x1)*y0 - y01;
    weights[0] = det*fabs( x0 );

    // Triangle 2
    reference_points_y[1][0] = x0 + (1.0 - x0)*y0;
    reference_points_y[1][1] = x1 - (      x1)*y0 + y01;
    weights[1] = det*fabs( 1.0 - x0 );

    // Triangle 3
    reference_points_y[2][0] = x0 + (1.0 - x0)*y0 + (x0 - x1)*y01;
    reference_points_y[2][1] = x1 + (1.0 - x1)*y0;
    weights[2] = det*fabs( (1.0- x1)*(x0 - 1.0) );

    // Triangle 4
    reference_points_y[3][0] = x0 -  x0*y01;
    reference_points_y[3][1] = x1 + (1.0 - x1)*y0;
    weights[3] = det*fabs( (1.0- x1)*x0 );
  }








  /**
   * Subdivision of the ReferenceElement into two subtriangles for the
   * VRTX_ADJACENT transformation (* marks the collocation point).
   *
   * \verbatim
   ^ eta
   |
   o-------------------o
   |                 / |
   |               /   |
   |     2       /     |
   |           /       |
   |         /         |
   |       /           |
   |     /       1     |
   |   /               |
   | /                 |
   *-------------------o--> xi
   \endverbatim
  */
  template<> inline
  void BaseColloDuffy<QUADRANGLE,VRTX_ADJACENT>::
  transformDuffy(const Point2& x,
                 const Point2& y,
                 point_type& reference_points_y,
                 weight_type& weights)
  {
    double y0  = y[0];
    double y01 = y[0]*y[1];
    double det = fabs(y0);

    // Triangle 1
    reference_points_y[0][0] = y0;
    reference_points_y[0][1] = y01;
    weights[0] = det;

    // Triangle 2
    reference_points_y[1][0] = y0 - y01;
    reference_points_y[1][1] = y0;
    weights[1] = det;
  }

} // end namespace hyena

#endif //colloduffy_h
