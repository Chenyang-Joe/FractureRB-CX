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
 * @file    referenceelements.tpl
 * @author  Mathias, Michael, Rf, TT, BK
 * @date    created:     04.08.09
 *          last change: 25.08.11
 */
#ifndef referenceelement_tpl
#define referenceelement_tpl


//------------------------------- specialization -------------------------------

//------------------------------- line, constant -------------------------------
template<>
const unsigned int ReferenceElement<LINE, CONSTANT>::
inner_node_ids[num_inner_nodes] = {0};

#ifdef _MSC_VER	
template<>
const unsigned int ReferenceElement<LINE, CONSTANT>::
corner_node_ids[]={-1};
#else
template<>
const unsigned int ReferenceElement<LINE, CONSTANT>::
corner_node_ids[num_corner_nodes]={};
#endif
template<>
const ReferenceElement<LINE, CONSTANT>::local_point_type
ReferenceElement<LINE, CONSTANT>::
local_points_[] = {local_point_type(.0)};

template<> inline
void ReferenceElement<LINE, CONSTANT>::
evaluateShapeFun(const Point1& p, Vec1d& shape_fun)
{
  shape_fun[0] = 1.;
}



//------------------------------- line, linear ---------------------------------
#ifdef _MSC_VER	
template<>
const unsigned int ReferenceElement<LINE, LINEAR>::
inner_node_ids[] = {-1};
#else
template<>
const unsigned int ReferenceElement<LINE, LINEAR>::
inner_node_ids[num_inner_nodes] = {};
#endif

template<>
const unsigned int ReferenceElement<LINE, LINEAR>::
corner_node_ids[num_corner_nodes]={0,1};

template<>
const ReferenceElement<LINE, LINEAR>::local_point_type
ReferenceElement<LINE, LINEAR>::
local_points_[] = {local_point_type(-1.0),
                   local_point_type(1.0)};
template<> inline
void ReferenceElement<LINE, LINEAR>::
evaluateShapeFun(const Point1& p, Vec2d& shape_fun)
{
  shape_fun[0] = .5*(1.-p[0]);
  shape_fun[1] = .5*(1.+p[0]);
}



//---------------------------- line, quadratic ---------------------------------
template<>
const unsigned int ReferenceElement<LINE, QUADRATIC>::
inner_node_ids[num_inner_nodes] = {1};

template<>
const unsigned int ReferenceElement<LINE, QUADRATIC>::
corner_node_ids[num_corner_nodes]={0,2};

template<>
const ReferenceElement<LINE, QUADRATIC>::local_point_type
ReferenceElement<LINE, QUADRATIC>::
local_points_[] = {local_point_type(-1.0),
                   local_point_type(.0),
                   local_point_type(1.0)};

template<> inline
void ReferenceElement<LINE, QUADRATIC>::
evaluateShapeFun(const Point1& p, Vec3d& shape_fun)
{
  shape_fun[0] = .5 * p[0] * (p[0]-1.);
  shape_fun[1] = .5 * p[0] * (p[0]+1.);
  shape_fun[2] = (1.- p[0] *  p[0]);
}



//----------------------- --- triangle, constant -------------------------------
template<>
const unsigned int ReferenceElement<TRIANGLE, CONSTANT>::
inner_node_ids[num_inner_nodes] = {0};

#ifdef _MSC_VER	
template<>
const unsigned int ReferenceElement<TRIANGLE, CONSTANT>::
corner_node_ids[]={-1};
#else
template<>
const unsigned int ReferenceElement<TRIANGLE, CONSTANT>::
corner_node_ids[num_corner_nodes]={};
#endif

#ifdef _MSC_VER	
template<>
const unsigned int ReferenceElement<TRIANGLE, CONSTANT>::
edge_node_ids[num_edges][num_edge_nodes+1] = {{-1},
                                            {-1},
                                            {-1}};
#else
template<>
const unsigned int ReferenceElement<TRIANGLE, CONSTANT>::
edge_node_ids[num_edges][num_edge_nodes] = {{},
                                            {},
                                            {}};
#endif

template<>
const ReferenceElement<TRIANGLE, CONSTANT>::local_point_type
ReferenceElement<TRIANGLE, CONSTANT>::
local_points_[] = {local_point_type(2.0/3.0,1.0/3.0)};

template<> inline
void ReferenceElement<TRIANGLE, CONSTANT>::
evaluateShapeFun(const Point2& p, Vec1d& shape_fun)
{
  shape_fun[0] = 1.;
}

//template<> inline
//void ReferenceElement<TRIANGLE, CONSTANT>::
//getColloPoints(Mat<Point2,1,1>& collo_pts, const double indent)
//{
//	collo_pts[0][0] = 2./3.;   	collo_pts[0][1] = 1./3.;
//}

template<> inline
const Mat<Point2,1,1> ReferenceElement<TRIANGLE, CONSTANT>::
getColloPoints(const double indent)
{
  Mat<Point2,1,1> collo_pts;
  collo_pts[0][0] = 2./3.;
  collo_pts[0][1] = 1./3.;
  return collo_pts;
}



template<> inline
void ReferenceElement<TRIANGLE, CONSTANT>::
evaluateGradient(const Point2& p, Mat<double,1,2>& gradient)
{
  gradient(0,0) = 0; 	gradient(0,1) = 0;
}

//----------------- triangle, linear, constant ansatz --------------------------
template<> template<> inline
const Point3 ReferenceElement<TRIANGLE, CONSTANT>::
getGlobalPoint(const Element<TRIANGLE, LINEAR> *const element,
               const unsigned int n)
{
  const Point3 p0 = element -> getNode(0) -> getPoint();
  const Point3 p1 = element -> getNode(1) -> getPoint();
  const Point3 p2 = element -> getNode(2) -> getPoint();

  return (p0 + p1 + p2) / 3;
}



//---------------- triangle quadratic, constant ansatz -------------------------
// TODO:FIXME: stub, only valid for flat triangles
template<> template<> inline
const Point3 ReferenceElement<TRIANGLE, CONSTANT>::
getGlobalPoint(const Element<TRIANGLE, QUADRATIC> *const element,
               const unsigned int n)
{
  const Point3 p0 = element -> getNode(0) -> getPoint();
  const Point3 p1 = element -> getNode(1) -> getPoint();
  const Point3 p2 = element -> getNode(2) -> getPoint();

  return (p0 + p1 + p2) / 3;
}



//--------------------------- triangle, linear ---------------------------------
#ifdef _MSC_VER	
template<>
const unsigned int ReferenceElement<TRIANGLE, LINEAR>::
inner_node_ids[] = {-1};
#else
template<>
const unsigned int ReferenceElement<TRIANGLE, LINEAR>::
inner_node_ids[num_inner_nodes] = {};
#endif

template<>
const unsigned int ReferenceElement<TRIANGLE, LINEAR>::
corner_node_ids[num_corner_nodes]={0,1,2};

#ifdef _MSC_VER	
template<>
const unsigned int ReferenceElement<TRIANGLE, LINEAR>::
edge_node_ids[num_edges][num_edge_nodes+1] = {{-1},
                                            {-1},
                                            {-1}};
#else
template<>
const unsigned int ReferenceElement<TRIANGLE, LINEAR>::
edge_node_ids[num_edges][num_edge_nodes] = {{},
                                            {},
                                            {}};
#endif

template<> inline
void ReferenceElement<TRIANGLE, LINEAR>::
evaluateShapeFun(const Point2& p, Vec3d& shape_fun)
{
  shape_fun[0] = 1. - p[0]       ;
  shape_fun[1] =      p[0] - p[1];
  shape_fun[2] =             p[1];
}

//template<> inline
//void ReferenceElement<TRIANGLE, LINEAR>::
//getColloPoints(Mat<Point2,3,1>& collo_pts, const double indent) 
//{
//	collo_pts[0][0] = 2. * indent;   collo_pts[0][1] =         indent;
//	collo_pts[1][0] = 1. - indent;   collo_pts[1][1] =         indent;
//	collo_pts[2][0] = 1. - indent;   collo_pts[2][1] = (1 - 2.*indent);
//}

template<> inline
const Mat<Point2,3,1> ReferenceElement<TRIANGLE, LINEAR>::
getColloPoints(const double indent) 
{
  Mat<Point2,3,1> collo_pts;
  collo_pts[0][0] = 2. * indent;   collo_pts[0][1] =         indent;
  collo_pts[1][0] = 1. - indent;   collo_pts[1][1] =         indent;
  collo_pts[2][0] = 1. - indent;   collo_pts[2][1] = (1 - 2.*indent);
  return collo_pts;
}

template<> inline
void ReferenceElement<TRIANGLE, LINEAR>::
evaluateGradient(const Point2& p, Mat<double,3,2>& gradient)
{
  gradient(0,0) = -1;   gradient(0,1) =  0;
  gradient(1,0) =  1; 	gradient(1,1) = -1;
  gradient(2,0) =  0; 	gradient(2,1) =  1;
}



//----------------- triangle, linear, linear ansatz ----------------------------
template<>
const ReferenceElement<TRIANGLE, LINEAR>::local_point_type
ReferenceElement<TRIANGLE, LINEAR>::
local_points_[] = {local_point_type(0,0),
                   local_point_type(1,0),
                   local_point_type(1,1)};

template<> template<> inline
const Point3 ReferenceElement<TRIANGLE, LINEAR>::
getGlobalPoint(const Element<TRIANGLE, LINEAR> *const element,
               const unsigned int n)
{
  return ( element -> getNode(n)->getPoint() );
}

// template<> inline
// const Point2 ReferenceElement<TRIANGLE, LINEAR>::
// getLocalPoint(const unsigned int n)
// {
//   Mat<Point2, 3, 1> local_points;
//   local_points[0] = Point2(0.,0.);
//   local_points[1] = Point2(1.,0.);
//   local_points[2] = Point2(1.,1.);

// 	return local_points[n];
// }



//--------------------------- triangle, quadratic ------------------------------
#ifdef _MSC_VER	
template<>
const unsigned int ReferenceElement<TRIANGLE, QUADRATIC>::
inner_node_ids[] = {-1};
#else
template<>
const unsigned int ReferenceElement<TRIANGLE, QUADRATIC>::
inner_node_ids[num_inner_nodes] = {};
#endif
template<>
const unsigned int ReferenceElement<TRIANGLE, QUADRATIC>::
corner_node_ids[num_corner_nodes]={0,1,2};

#ifdef _MSC_VER	
template<>
const unsigned int ReferenceElement<TRIANGLE, QUADRATIC>::
edge_node_ids[num_edges][num_edge_nodes+1] = {{3},
                                            {4},
                                            {5}};
#else
template<>
const unsigned int ReferenceElement<TRIANGLE, QUADRATIC>::
edge_node_ids[num_edges][num_edge_nodes] = {{3},
                                            {4},
                                            {5}};
#endif

template<>
const ReferenceElement<TRIANGLE, QUADRATIC>::local_point_type
ReferenceElement<TRIANGLE, QUADRATIC>::
local_points_[] = {local_point_type(0.0,0.0),
                   local_point_type(1.0,0.0),
                   local_point_type(1.0,1.0),
                   local_point_type(0.5,0.0),
                   local_point_type(1.0,0.5),
                   local_point_type(0.5,0.5)};

template<> inline
void ReferenceElement<TRIANGLE, QUADRATIC>::
evaluateShapeFun(const Point2& p, Vec6d& shape_fun)
{
  double lam1 = 1. - p[0]       ;
  double lam2 =      p[0] - p[1];
  double lam3 =             p[1];

  shape_fun[0] = lam1*(2.*lam1- 1);
  shape_fun[1] = lam2*(2.*lam2- 1);
  shape_fun[2] = lam3*(2.*lam3- 1);
  shape_fun[3] = 4.*lam2*lam1;
  shape_fun[4] = 4.*lam2*lam3;
  shape_fun[5] = 4.*lam3*lam1;
}

//template<> inline
//void ReferenceElement<TRIANGLE, QUADRATIC>::
//getColloPoints(Mat<Point2,6,1>& collo_pts, const double indent) 
//{
//	collo_pts[0][0]= 2.*indent;    	    collo_pts[0][1]= indent;
//	collo_pts[1][0]= 1. - indent;     	collo_pts[1][1]= indent;
//	collo_pts[2][0]= 1. - indent;     	collo_pts[2][1]= 1. - 2.*indent;
//	collo_pts[3][0]= 0.5*(1.+ indent);	collo_pts[3][1]= indent;
//	collo_pts[4][0]= 1. - indent;     	collo_pts[4][1]= 0.5*(1.+indent);
//	collo_pts[5][0]= 0.5*(1.+ indent);	collo_pts[5][1]= 0.5*(1.+indent);
//}

template<> inline
const Mat<Point2,6,1> ReferenceElement<TRIANGLE, QUADRATIC>::
getColloPoints(const double indent) 
{
  Mat<Point2,6,1> collo_pts;
  collo_pts[0][0]= 2.*indent;    	    collo_pts[0][1]= indent;
  collo_pts[1][0]= 1. - indent;     	collo_pts[1][1]= indent;
  collo_pts[2][0]= 1. - indent;     	collo_pts[2][1]= 1. - 2.*indent;
  collo_pts[3][0]= 0.5*(1.+ indent);	collo_pts[3][1]= indent;
  collo_pts[4][0]= 1. - indent;     	collo_pts[4][1]= 0.5*(1.+indent);
  collo_pts[5][0]= 0.5*(1.+ indent);	collo_pts[5][1]= 0.5*(1.+indent);
  return collo_pts;
}

template<> inline
void ReferenceElement<TRIANGLE, QUADRATIC>::
evaluateGradient(const Point2& p, Mat<double,6,2>& gradient)
{
  gradient(0,0) =  4. * p[0]             - 3.;   
  gradient(1,0) =  4. * p[0] - 4. * p[1] - 1.; 	
  gradient(2,0) =                          0.;  
  gradient(3,0) = -8. * p[0] + 4. * p[1] + 4.;   
  gradient(4,0) =              4. * p[1]     ;   
  gradient(5,0) =             -4. * p[1]     ;
	
  gradient(0,1) =                          0.;
  gradient(1,1) = -4. * p[0] + 4. * p[1] + 1.;
  gradient(2,1) =              4. * p[1] - 1.;
  gradient(3,1) =  4. * p[0]             - 4.;
  gradient(4,1) =  4. * p[0] - 8. * p[1]     ;
  gradient(5,1) = -4. * p[0]             + 4.;
		
}



//--------------------------- quadrangle, constant -----------------------------
template<>
const unsigned int ReferenceElement<QUADRANGLE, CONSTANT>::
inner_node_ids[num_inner_nodes] = {0};

#ifdef _MSC_VER	
template<>
const unsigned int ReferenceElement<QUADRANGLE, CONSTANT>::
corner_node_ids[]={-1};
#else
template<>
const unsigned int ReferenceElement<QUADRANGLE, CONSTANT>::
corner_node_ids[num_corner_nodes]={};
#endif

#ifdef _MSC_VER	
template<>
const unsigned int ReferenceElement<QUADRANGLE, CONSTANT>::
edge_node_ids[num_edges][num_edge_nodes+1] = {{-1},
                                            {-1},
                                            {-1},
                                            {-1}};
#else
template<>
const unsigned int ReferenceElement<QUADRANGLE, CONSTANT>::
edge_node_ids[num_edges][num_edge_nodes] = {{},
                                            {},
                                            {},
                                            {}};
#endif

template<>
const ReferenceElement<QUADRANGLE, CONSTANT>::local_point_type
ReferenceElement<QUADRANGLE, CONSTANT>::
local_points_[] = {local_point_type(.5,.5)};

template<> inline
void ReferenceElement<QUADRANGLE, CONSTANT>::
evaluateShapeFun(const Point2& p, Vec1d& shape_fun)
{
  shape_fun[0] = 1.;
}

// template<> inline
// void ReferenceElement<QUADRANGLE, CONSTANT>::
// getColloPoints(Mat<Point2,1,1>& collo_pts, const double indent) 
// {
// 	collo_pts[0][0] = 0.5;    	collo_pts[0][1] = 0.5;
// }

template<> inline
const Mat<Point2,1,1> ReferenceElement<QUADRANGLE, CONSTANT>::
getColloPoints(const double indent) 
{
  Mat<Point2,1,1> collo_pts;
  collo_pts[0][0] = 0.5;
  collo_pts[0][1] = 0.5;
  return collo_pts;
}

template<> inline
void ReferenceElement<QUADRANGLE, CONSTANT>::
evaluateGradient(const Point2& p, Mat<double,1,2>& gradient)
{
  gradient(0,0) = 0; 	gradient(0,1) = 0;
}

// specialisation for linear geometry, constant ansatz
template<> template<> inline
const Point3 ReferenceElement<QUADRANGLE, CONSTANT>::
getGlobalPoint(const Element<QUADRANGLE, LINEAR> *const element,
               const unsigned int n)
{
  const Point3 p0 = element -> getNode(0)->getPoint();
  const Point3 p1 = element -> getNode(1)->getPoint();
  const Point3 p2 = element -> getNode(2)->getPoint();
  const Point3 p3 = element -> getNode(3)->getPoint();

  return (p0 + p1 + p2 + p3) / 4;
}



//-------------------------- quadrangle, linear --------------------------------
#ifdef _MSC_VER	
template<>
const unsigned int ReferenceElement<QUADRANGLE, LINEAR>::
inner_node_ids[] = {-1};
#else
template<>
const unsigned int ReferenceElement<QUADRANGLE, LINEAR>::
inner_node_ids[num_inner_nodes] = {};
#endif

template<>
const unsigned int ReferenceElement<QUADRANGLE, LINEAR>::
corner_node_ids[num_corner_nodes]={0,1,2,3};

#ifdef _MSC_VER	
template<>
const unsigned int ReferenceElement<QUADRANGLE, LINEAR>::
edge_node_ids[num_edges][num_edge_nodes+1] = {{-1},
                                            {-1},
                                            {-1},
                                            {-1}};
#else
template<>
const unsigned int ReferenceElement<QUADRANGLE, LINEAR>::
edge_node_ids[num_edges][num_edge_nodes] = {{},
                                            {},
                                            {},
                                            {}};
#endif

template<>
const ReferenceElement<QUADRANGLE, LINEAR>::local_point_type
ReferenceElement<QUADRANGLE, LINEAR>::
local_points_[] = {local_point_type(0,0),
                   local_point_type(1,0),
                   local_point_type(1,1),
                   local_point_type(0,1)};

template<> inline
void ReferenceElement<QUADRANGLE, LINEAR>::
evaluateShapeFun(const Point2& p, Vec4d& shape_fun)
{
  shape_fun[0] = ( 1. - p[0] ) * ( 1. - p[1] );
  shape_fun[1] =        p[0]   * ( 1. - p[1] );
  shape_fun[2] =        p[0]   *        p[1];
  shape_fun[3] = ( 1. - p[0] ) *        p[1];
}

template<> inline
const Mat<Point2,4,1> ReferenceElement<QUADRANGLE, LINEAR>::
getColloPoints(const double indent) 
{
  Mat<Point2,4,1> collo_pts;
  collo_pts[0][0] =      indent;  	collo_pts[0][1] =      indent;
  collo_pts[1][0] = 1. - indent;  	collo_pts[1][1] =      indent;
  collo_pts[2][0] = 1. - indent;  	collo_pts[2][1] = 1. - indent;
  collo_pts[3][0] =      indent;  	collo_pts[3][1] = 1. - indent;
  return collo_pts;
}

template<> inline
void ReferenceElement<QUADRANGLE, LINEAR>::
evaluateGradient( const Point2& p, Mat<double,4,2>& grad )
{
  grad(0,0)=-1.+p[1];  grad(0,1)=-1.+p[0]; // phi0
  grad(1,0)= 1.-p[1];  grad(1,1)=   -p[0]; // phi1
  grad(2,0)=    p[1];  grad(2,1)=    p[0]; // phi2
  grad(3,0)=   -p[1];  grad(3,1)= 1.-p[0]; // phi3
}



//------------------ quadrangle, linear, linear ansatz -------------------------
template<> template<> inline
const Point3 ReferenceElement<QUADRANGLE, LINEAR>::
getGlobalPoint(const Element<QUADRANGLE, LINEAR> *const element,
               const unsigned int n)
{
  return ( element -> getNode(n)->getPoint() );
}


// template<> inline
// const Point2 ReferenceElement<QUADRANGLE, LINEAR>::
// getLocalPoint(const unsigned int n)
// {
//   Mat<Point2, 4, 1> local_points;
//   local_points[0] = Point2(0.,0.);
//   local_points[1] = Point2(1.,0.);
//   local_points[2] = Point2(1.,1.);
//   local_points[3] = Point2(0.,1.);

// 	return local_points[n];
// }



//-------------------------- quadrangle, quadratic -----------------------------
template<>
const unsigned int ReferenceElement<QUADRANGLE, QUADRATIC>::
inner_node_ids[num_inner_nodes] = {8};

template<>
const unsigned int ReferenceElement<QUADRANGLE, QUADRATIC>::
corner_node_ids[num_corner_nodes]={0,1,2,3};

#ifdef _MSC_VER
const unsigned int ReferenceElement<QUADRANGLE, QUADRATIC>::
edge_node_ids[num_edges][num_edge_nodes+1] = {{4},
                                            {5},
                                            {6},
                                            {7}};
#else
template<>
const unsigned int ReferenceElement<QUADRANGLE, QUADRATIC>::
edge_node_ids[num_edges][num_edge_nodes] = {{4},
                                            {5},
                                            {6},
                                            {7}};
#endif
template<>
const ReferenceElement<QUADRANGLE, QUADRATIC>::local_point_type
ReferenceElement<QUADRANGLE, QUADRATIC>::
local_points_[] = {local_point_type(0,0),
                   local_point_type(1,0),
                   local_point_type(1,1),
                   local_point_type(0,1),
                   local_point_type(0.5,0),
                   local_point_type(1,0.5),
                   local_point_type(0.5,1),
                   local_point_type(0,0.5),
                   local_point_type(0.5,0.5)};

template<> inline
const Mat<Point2,9,1> ReferenceElement<QUADRANGLE, QUADRATIC>::
getColloPoints(const double indent) 
{
  Mat<Point2,9,1> collo_pts;
  collo_pts[0][0] =      indent;   collo_pts[0][1] =      indent;
  collo_pts[1][0] = 1. - indent;   collo_pts[1][1] =      indent;
  collo_pts[2][0] = 1. - indent;   collo_pts[2][1] = 1. - indent;
  collo_pts[3][0] =      indent;   collo_pts[3][1] = 1. - indent;
  collo_pts[4][0] = 0.5        ;   collo_pts[4][1] = indent 		;
  collo_pts[5][0] = 1. - indent;   collo_pts[5][1] = 0.5 				;
  collo_pts[6][0] = 0.5        ;   collo_pts[6][1] = 1. - indent;
  collo_pts[7][0] = indent     ;   collo_pts[7][1] = 0.5 				;
  collo_pts[8][0] = 0.5        ;   collo_pts[8][1] = 0.5        ;   
  return collo_pts;
}

template<> inline
void ReferenceElement<QUADRANGLE, QUADRATIC>::
evaluateShapeFun(const Point2& p, Vec9d& shape_fun)
{ 
  shape_fun[0] = (2*p[0] - 1)*(2*p[1] - 1)*(p[0] - 1)*(p[1] - 1);
  shape_fun[1] = p[0]*(2*p[0] - 1)*(2*p[1] - 1)*(p[1] - 1);
  shape_fun[2] = p[0]*p[1]*(2*p[0] - 1)*(2*p[1] - 1);
  shape_fun[3] = p[1]*(2*p[0] - 1)*(2*p[1] - 1)*(p[0] - 1);
  shape_fun[4] = (-p[0])*(4*p[0] - 4)*(2*p[1] - 1)*(p[1] - 1);
  shape_fun[5] = (-p[0])*p[1]*(2*p[0] - 1)*(4*p[1] - 4);
  shape_fun[6] = (-p[0])*p[1]*(4*p[0] - 4)*(2*p[1] - 1);
  shape_fun[7] = (-p[1])*(2*p[0] - 1)*(4*p[1] - 4)*(p[0] - 1);
  shape_fun[8] = p[0]*p[1]*(4*p[0] - 4)*(4*p[1] - 4);
}

template<> inline
void ReferenceElement<QUADRANGLE, QUADRATIC>::
evaluateGradient( const Point2& p, Mat<double,9,2>& grad )
{
  grad(0,0) = 2*(2*p[1] - 1)*(p[0] - 1)*(p[1] - 1) 
    + (2*p[0] - 1)*(2*p[1] - 1)*(p[1] - 1);
  grad(1,0) =(2*p[0] - 1)*(2*p[1] - 1)*(p[1] - 1) 
    + 2*p[0]*(2*p[1] - 1)*(p[1] - 1);
  grad(2,0) =p[1]*(2*p[0] - 1)*(2*p[1] - 1) + 2*p[0]*p[1]*(2*p[1] - 1);
  grad(3,0) =p[1]*(2*p[0] - 1)*(2*p[1] - 1) + 2*p[1]*(2*p[1] - 1)*(p[0] - 1);
  grad(4,0) =- (4*p[0] - 4)*(2*p[1] - 1)*(p[1] - 1) 
    - 4*p[0]*(2*p[1] - 1)*(p[1] - 1);
  grad(5,0) =(-p[1])*(2*p[0] - 1)*(4*p[1] - 4) - 2*p[0]*p[1]*(4*p[1] - 4);
  grad(6,0) =(-p[1])*(4*p[0] - 4)*(2*p[1] - 1) - 4*p[0]*p[1]*(2*p[1] - 1);
  grad(7,0) =(-p[1])*(2*p[0] - 1)*(4*p[1] - 4) - 2*p[1]*(4*p[1] - 4)*(p[0] - 1);
  grad(8,0) =p[1]*(4*p[0] - 4)*(4*p[1] - 4) + 4*p[0]*p[1]*(4*p[1] - 4);

  grad(0,1) =2*(2*p[0] - 1)*(p[0] - 1)*(p[1] - 1) 
    + (2*p[0] - 1)*(2*p[1] - 1)*(p[0] - 1);
  grad(1,1) =p[0]*(2*p[0] - 1)*(2*p[1] - 1) + 2*p[0]*(2*p[0] - 1)*(p[1] - 1);
  grad(2,1) =p[0]*(2*p[0] - 1)*(2*p[1] - 1) + 2*p[0]*p[1]*(2*p[0] - 1);
  grad(3,1) =(2*p[0] - 1)*(2*p[1] - 1)*(p[0] - 1) 
    + 2*p[1]*(2*p[0] - 1)*(p[0] - 1);
  grad(4,1) =(-p[0])*(4*p[0] - 4)*(2*p[1] - 1) - 2*p[0]*(4*p[0] - 4)*(p[1] - 1);
  grad(5,1) =(-p[0])*(2*p[0] - 1)*(4*p[1] - 4) - 4*p[0]*p[1]*(2*p[0] - 1);
  grad(6,1) =(-p[0])*(4*p[0] - 4)*(2*p[1] - 1) - 2*p[0]*p[1]*(4*p[0] - 4);
  grad(7,1) =-(2*p[0] - 1)*(4*p[1] - 4)*(p[0] - 1) 
    - 4*p[1]*(2*p[0] - 1)*(p[0] - 1);
  grad(8,1) =p[0]*(4*p[0] - 4)*(4*p[1] - 4) + 4*p[0]*p[1]*(4*p[0] - 4);
}
#endif // include guard
