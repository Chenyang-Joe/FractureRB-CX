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

//file: enumerators.H
//--------------------------------------------------------------------------
#ifndef enumerators_h
#define enumerators_h
//--------------------------------------------------------------------------
/** \file
 *  Contains Enumerators for the BEM code
 */
//--------------------------------------------------------------------------
#include <iostream>
#include <string>

namespace hyena
{

  /** Enumerator for element shapes */
  //enum ELEM_SHAPE{
  enum ELEMENT_SHAPE{
    /* 0 */ NO_SHAPE,       //!< element shape is not assigned
    /* 1 */ LINE,           //!< line element
    /* 2 */ TRIANGLE,       //!< triangle element
    /* 3 */ QUADRANGLE,     //!< quadrangle element
    /* 4 */ HYPERCUBE,      //!< hypercube element
    /* 5 */ INFINITE_LINE,  //!< infinite 2-d line element
    /* 6 */ INFINITE_QUAD   //!< infinite 3-d surface element
  };

  /** Enumerator for approximation order */
  enum APPROXIMATION{
    /* 0 */ NO_APPROX_ORDER,  //!< no approximation order assigned
    /* 1 */ CONSTANT,         //!< constant approximation
    /* 2 */ LINEAR,           //!< linear approximation
    /* 3 */ QUADRATIC         //!< quadratic approximation
  };

  /** Enumerator for type of Boundary Conditions */
  enum BC_TYPE {
    /* 0 */ NO_BC_TYPE,           //!< Type of BC is not initialized
    /* 1 */ INTERFACE,            //!< Interface boundary
    /* 2 */ MIXED,                //!< BC differs per dof
    /* 3 */ ROBIN,                //!< Robin boundary condition
    /* 4 */ ANY_BC,               //!< any BC type
    /* 5 */ FIXED,                //!< fixed edge (zero Dirichlet conditions)
    /* 6 */ PRESSURE,             //!< (normal) pressure
    /* 7*/ SHEAR,                //!< tangential forces
    /* 8*/ MOMENT,               //!< bending moment  
    /* 9*/ NEUMANN,              //!< Neumann bcs
    /* 10*/ DIRICHLET,             //!< Dirichlet bcs
            CRACK
  };

  /** Enumerator for type of Boundary Conditions */
  enum LDOF_TYPE {
    /* 0 */ NO_LDOF_TYPE,         //!< Type of ldof is not initialized
    /* 1 */ HOMOGEN,              //!< homogeneous boundary condition
    /* 2 */ KNOWN,                //!< inhomogeneous boundary condition
    /* 3 */ UNKNOWN,              //!< Boundary data is unknown
    /* 4 */ RUNKNOWN,             //!< Boundary data is robin and unknwon
    /* 5 */ IUNKNOWN,             //!< Boundary is interface, data is unknwon
    /* 6 */ ANY_LDOF_TYPE,
            CUNKNOWN, // unkown crack opening displacement
			CHOMOGEN
  };

  /** Enumerator for type of shape functions*/
  enum SPACE_TYPE {
    /* 0 */  NO_SPACE_SET,        //!< type of shape function not yet set
    /* 1 */  CONTINUOUS,          //!< continuous shape function
    /* 2 */  DISCONTINUOUS        //!< discontinuous shape function
  };

  /** Enumerator describing the Type of Singular Integration */
  enum SING_INT {
    /* 0 */ NO_SING_INT,    //!< singularity type is not set yet
    /* 1 */ REGULAR,        //!< Regular integration
    /* 2 */ COINCIDENT,     //!< coincident
    /* 3 */ EDGE_ADJACENT,  //!< edge adjacent
    /* 4 */ VRTX_ADJACENT  //!< vertex adjacent
  };

  /** Enumerator for quadrature rule */
  enum QUADRATURE_RULE {
    /* 0 */ NO_RULE,          //!< Quadrature Rule not (yet) set
    /* 1 */ GAUSS,            //!<
  };

  //TT: pde_fsol_refactoring
//  /** Enumerator for partial differential equation */
//  enum PDE{
//    /* 0 */ NO_PDE,              //!< PDE is not (yet) set
//    /* 1 */ LAPLACE_3D,          //!<
//    /* 2 */ HELMHOLTZ_3D,        //!<
//    /* 3 */ WAVEEQUATION_3D,     //!<
//    /* 4 */ HEATEQUATION_3D,     //!<
//    /* 5 */ ELASTOSTATIC_3D,     //!<
//    /* 6 */ ELASTODYNAMIC_3D,    //!<
//    /* 7 */ SATPOROLAPLACE_3D,    //!<
//  };

  enum PROBLEM{
    /* 0 */ NO_PROBLEM_SET,       //!< PDE is not (yet) set
    /* 1 */ REAL1_3D,          //!<
    /* 2 */ REAL3_3D,        //!<
    /* 3 */ REAL4_3D,        //!<
    /* 4 */ COMPLEX1_3D,     //!<
    /* 5 */ COMPLEX3_3D,     //!<
    /* 6 */ COMPLEX4_3D     //!<
  }; //TODO: add problem for satporolaplace

  ///** Enumerator for material data */
  //enum MATERIAL {
  //	/* 0 */ NO_MATERIAL,          //!< Material is not (yet) set
  //	/* 1 */ FLUID,                //!< Fluid, e.g., scalar problem (Laplace eq.)
  //	/* 2 */ ELASTIC,              //!< Elastic material
  //	/* 3 */ POROELASTIC           //!< Poroelastodynamic material
  //};

  /** Enumerator for type of Analysis */
  enum ANALYSIS {
    /* 0 */ NO_ANALYSIS,          //!< Type of Analysis is not (yet) set
    /* 1 */ STATIC,               //!< Static analysis
    /* 2 */ LAPLACE,              //!< Analysis in Laplace domain
    /* 3 */ FREQUENCY,            //!< Analysis in Frequency domain
    /* 4 */ CQ_METHOD,            //!< Convolution Quadrature Method (Lubich)
    /* 5 */ TIME_DIRECT           //!< Analysis directly in time domain (Mansur)
  };

  /** Enumerator for projection method */
  enum PROJECTION {
    /* 0 */ NO_PROJECTION,        //!< Projection type is not set
    /* 1 */ COLLOCATION,          //!< collocate integral equation
    /* 2 */ GALERKIN              //!< weight the residual of integral equation
  };

  ///** Enumerator for the system structure */
  //enum STRUCTURE {
  //	/* 0 */ NO_STRUCTURE,         //!< structure is not set
  //	/* 1 */ UNSYMMETRIC,          //!< use only first integral equation
  //	/* 2 */ SYMMETRIC,            //!< use both integral equations
  //	/* 3 */ INTERIOR              //!< for the interior point evaluation
  //};

  ///** Enumerator for Approximation of unknowns, geometry and bcs */
  //enum FUNC {
  //	/* 0 */ NO_FUNC,              //!< Type of Approximation is not (yet) set
  //	/* 1 */ CONSTANT,             //!< Constant approximation
  //	/* 2 */ LINEAR,               //!< Linear approximation
  //	/* 3 */ QUADRATIC,            //!< Quadratic approximation
  //	/* 4 */ SEMILINEAR,           //!< Linear in just one coordinate direction
  //	/* 5 */ BILINEAR,             //!< Bilinear approximation
  //	/* 6 */ NONLINEAR             //!< Nonlinear approximation
  //};

  /** Enumerator for the types of matrix */
  enum POTENTIAL {
    /* 0 */ NO_POTENTIAL,   //!< no Potential
    /* 0 */ SINGLELAYER,    //!< Single layer potential
    /* 1 */ DOUBLELAYER,    //!< Double layer potential
    /* 2 */ ADJOINT_DLP,    //!< Adjoint double layer potential
    /* 3 */ HYPEROPERATOR,  //!< Hyper singular integral operator
    /* 4 */ NEWTONPOTENTIAL,//!< Newton Potential
    /* 5 */ IDENTITY        //!< Identity
  };

  ///** Enumerator for the types of H - matrix */
  //enum H_MAT_TYPE {
  //	/* 0 */ H_STA,    //!< standart unsymmetric H - matrix
  //	/* 1 */ H_SYM,    //!< symmetric H - matrix
  //	/* 2 */ H_LTR,    //!< lower triangle H - matrix
  //	/* 3 */ H_UTR,    //!< upper triangle H - matrix
  //	/* 4 */ H_TST,    //!< transposed standart H - matrix
  //	/* 5 */ H_TSY,    //!< transposed symmetric H - matrix
  //	/* 6 */ H_MAS,    //!< mass H - matrix
  //	/* 7 */ NO_H      //!< no H - matrix
  //};

  ///** Enumerator for time-dependent boundary conditions */
  //enum TIME_FUN {
  //	/* 0 */ HEAVISIDE,   //!< standard Heaviside function
  //	/* 1 */ SINE,        //!< sine wave
  //	/* 2 */ IMPULSE,     //!< impulse of certain length
  //	/* 3 */ SAWTOOTH,    //!< inconceivable creation
  //	/* 4 */ RELAX,       //!< function creeping to zero
  //	/* 5 */ SIGMOIDAL,   //!< sigmoidal function
  //	/* 6 */ NO_TIME_FUN  //!< not set!
  //};

  ///** Enumerator to identify fundamental solutions */
  //enum FUNDSOL {
  //	/* 0 */ NO_FUNDSOL,         //!< fundamental solution not set
  //	/* 1 */ FUND_LAPLACE,       //!< laplace equation
  //	/* 2 */ FUND_HELMHOLTZ,     //!< helmholtz equation
  //	/* 3 */ FUND_WAVE,          //!< wave equation
  //	/* 4 */ FUND_ELASTOSTATIC,  //!< system of elastostatic
  //	/* 5 */ FUND_ELASTODYNAMIC, //!< elastodynamics
  //	/* 6 */ FUND_VISCO,         //!< viscous media
  //	/* 7 */ FUND_PORO           //!< poroelasto-something
  //};

  enum FIGURE {
    /*  0 */ NO_FIGURE,
    /*  1 */ QUAD,
    /*  2 */ CUBOID,
    /*  3 */ CIRCLE,
    /*  4 */ SPHERE,
    /*  5 */ HALFSPACE,
    /*  6 */ QUAD2,
    /*  7 */ HALFPLANE,
    /*  8 */ CYLINDER,
    /*  9 */ CUBOID2,
    /* 10 */ ABAQUS,
    /* 11 */ EXODUS
  };

  /** Enumerator to set geometry type { OPEN | CLOSED } */
  enum TOPOLOGY {
    /* 0 */ NO_TOPOLOGY,        //!< geometry type not set
    /* 1 */ OPEN,               //!< open surface
    /* 2 */ CLOSED              //!< closed surface
  };

  ///** Enumerator to set fast method type { NONE | ACA } */
  //enum FASTMETHOD {
  //	/* 0 */ NO_FASTMETHOD,      //!< no fast method
  //	/* 1 */ ACA                 //!< adaptive cross approximation
  //};

  //** Enumerator to set element type */
  // enum ELEM_TYPE {
  //  /* 0 */ NO_ELEM_TYPE,       //!< element type not set
  //  /* 1 */ LINE,               //!< line element
  //  /* 2 */ TRIANGLE,           //!< generic triangular element
  //  /* 3 */ TRIA3,              //!< 3 point triangular element
  //  /* 4 */ TRIA6,              //!< 6 point triangular element
  //  /* 5 */ QUADRANGLE,         //!< generic quadratic element
  //  /* 6 */ QUAD4,              //!< 4 point quadrilateral element
  //  /* 7 */ QUAD9,              //!< 9 point quadrilateral element
  //  /* 8 */ INFINITE_LINE,      //!< infinite 2-d line element
  //  /* 9 */ INFINITE_QUAD       //!< infinite 3-d surface element
  // };

  ////** Enumerator to set if vertex of element is support or not */
  //enum SUPPORT {
  //	/* 0 */ NO_SUPPORT,         //!< vertex is no support
  //	/* 1 */ IS_SUPPORT,         //!< vertex is support
  //	/* 2 */ PRIMAL,             //!< primal support
  //	/* 3 */ DUAL                //!< dual support
  //};

  ////** Enumerator to set the sign of system matrix block */
  //enum SIGN {
  //	/*  0 */ NO_SIGN =  0,      //!< no sign
  //	/* -1 */ MINUS   = -1,      //!< minus
  //	/* +1 */ PLUS    =  1       //!< plus
  //};

  ///** Enumerator to define Newton potentials */
  //enum NEWTON_POTENTIAL {
  //	/* 0 */ NO_NEWTON_POTENTIAL,  //!< newton potential not set
  //	/* 1 */ POINT_LOAD,           //!< a point load
  //	/* 2 */ DOMAIN_LOAD           //!< a constant domain load
  //};



  //----------------------------------------------------------------------------
  /* overloading operator<< for enums */
  //----------------------------------------------------------------------------

  inline std::ostream& operator<<(std::ostream& out, const ELEMENT_SHAPE& shape)
  {
    switch(shape) {
    case NO_SHAPE:      out << "element has no shape";         break;
    case LINE:          out << "Line";                         break;
    case TRIANGLE:      out << "Triangle";                     break;
    case QUADRANGLE:    out << "Quadrangle";                   break;
    case HYPERCUBE:     out << "Hypercube element";            break;
    case INFINITE_LINE: out << "Infinite 2-d line element";    break;
    case INFINITE_QUAD: out << "Infinite 3-d surface element"; break;

    default:         out << "SPACE_TYPE not initialized"; break;
    }
    return out;
  }

  inline std::ostream& operator<<(std::ostream& out, const APPROXIMATION& approx)
  {
    switch(approx) {
    case NO_APPROX_ORDER: out << "no approximation order assigned";  break;
    case CONSTANT:        out << "CONSTANT";                         break;
    case LINEAR:          out << "LINEAR";                           break;
    case QUADRATIC:       out << "QUADRATIC";                        break;
    default:              out << "APPROX_ORDER not initialized";     break;
    }
    return out;
  }

  inline std::ostream& operator<<(std::ostream& out, const BC_TYPE& bt)
  {
    switch(bt) {
    case NO_BC_TYPE: out << "NO_BC_TYPE"; break;
    case INTERFACE:  out << "INTERFACE";  break;
    case MIXED:      out << "MIXED";      break;
    case ROBIN:      out << "ROBIN";      break;
    case ANY_BC:     out << "ANY_BC";     break;
    case NEUMANN:    out << "NEUMANN";    break;
    case DIRICHLET:  out << "DIRICHLET";  break;
	case CRACK:		 out << "CRACK";	  break;
    default:         out << "BC_TYPE not initialized";
    }
    return out;
  }

  inline std::ostream& operator<<(std::ostream& out, const LDOF_TYPE& bt)
  {
    switch(bt) {
    case NO_LDOF_TYPE: out << "NO_LDOF_TYPE"; break;
    case HOMOGEN:      out << "HOMOGEN";      break;
    case KNOWN:        out << "KNOWN";        break;
    case UNKNOWN:      out << "UNKNOWN";      break;
    case RUNKNOWN:     out << "RUNKNOWN";     break;
    case IUNKNOWN:     out << "IUNKNOWN";     break;
    case CUNKNOWN:      out << "CUNKNOWN";  break;
    default:           out << "LDOF_TYPE not initialized";
    }
    return out;
  }

  inline std::ostream& operator<<(std::ostream& out, const SPACE_TYPE& sf)
  {
    switch(sf) {
    case NO_SPACE_SET:  out << "NO_SPACE_SET"; break;
    case CONTINUOUS:    out << "CONTINUOUS";   break;
    case DISCONTINUOUS: out << "DISCONTINUOUS"; break;
    default:            out << "SPACE_TYPE not initialized";
    }
    return out;
  }

  inline std :: ostream & operator<<( std :: ostream & out,
                                      const SING_INT & si )
  {
    switch( si ) {
    case NO_SING_INT:    out << "NO_SING_INT"; break;
    case REGULAR:        out << "Regular integration"; break;
    case COINCIDENT:     out << "Coincident elements"; break;
    case EDGE_ADJACENT:  out << "EA"; break;
    case VRTX_ADJACENT:  out << "VA"; break;
    default:             out << "SING_INT not initialized";
    };
    return out;
  }


  inline std :: ostream & operator<<( std :: ostream & out,
                                      const  QUADRATURE_RULE& rule )
  {
    switch( rule ) {
    case NO_RULE:        out << "NO_QUADRATURE_RULE"; break;
    case GAUSS:          out << "Gauss quadrature"; break;
    default:             out << "QUADRATURE_RULE not initialized";
    };
    return out;
  }



  inline std :: ostream & operator<<( std :: ostream & out,
                                      const  PROBLEM& problem )
  {
    switch( problem ) {
    case NO_PROBLEM_SET : out << "NO_PROBLEM_SET"; break;
    case REAL1_3D:        out << "REAL1_3D"; break;
    case REAL3_3D:        out << "REAL3_3D"; break;
    case REAL4_3D:        out << "REAL4_3D"; break;
    case COMPLEX1_3D:        out << "COMPLEX1_3D"; break;
    case COMPLEX3_3D:        out << "COMPLEX3_3D"; break;
    case COMPLEX4_3D:        out << "COMPLEX4_3D"; break;
    default:                 out << "PROBLEM not initialized";
    };
    return out;
  }


  inline std::ostream& operator<<( std::ostream &out, const FIGURE & f )
  {
    switch( f ) {
    case NO_FIGURE: out << "NO_FIGURE"; break;
    case QUAD:      out << "Quad";      break;
    case QUAD2:     out << "Quad2";     break;
    case CUBOID:    out << "Cuboid";    break;
    case CUBOID2:   out << "Cuboid2";   break;
    case CIRCLE:    out << "Circle";    break;
    case SPHERE:    out << "Sphere";    break;
    case HALFSPACE: out << "Halfspace"; break;
    case HALFPLANE: out << "Halfplane"; break;
    case CYLINDER:  out << "Cylinder" ; break;
    case ABAQUS:    out << "Abaqus";    break;
    case EXODUS:    out << "Exodus";    break;
    default:        out << "FIGURE not initialized";
    }
    return out;
  }

  //inline std::ostream& operator<<(std::ostream &out, const MATERIAL &m) {
  //	switch( m ) {
  //	case NO_MATERIAL:  out << "NO_MATERIAL";  break;
  //	case FLUID:        out << "FLUID";        break;
  //	case ELASTIC:      out << "ELASTIC";      break;
  //	case POROELASTIC:  out << "POROELASTIC";  break;
  //	default:           out << "MATERIAL not initialized";
  //	}
  //	return out;
  //}

  inline std::ostream & operator<<( std :: ostream & out, const ANALYSIS & a )
  {
    switch( a ) {
    case NO_ANALYSIS: out << "NO_ANALYSIS"; break;
    case STATIC:      out << "STATIC";      break;
    case LAPLACE:     out << "LAPLACE";     break;
    case FREQUENCY:   out << "FREQUENCY";   break;
    case CQ_METHOD:   out << "CQ_METHOD";   break;
    case TIME_DIRECT: out << "TIME_DIRECT"; break;
    default:          out << "ANALYSIS not initialized";
    }
    return out;
  }

  inline std :: ostream & operator<<( std :: ostream & out,
                                      const PROJECTION & p )
  {
    switch( p ) {
    case NO_PROJECTION: out << "NO_PROJECTION"; break;
    case COLLOCATION:   out << "COLLOCATION";   break;
    case GALERKIN:      out << "GALERKIN";      break;
    default:            out << "PROJECTION not initialized";
    }
    return out;
  }

  //inline std :: ostream & operator<<( std :: ostream & out,
  //																		const STRUCTURE & s )
  //{
  //	switch( s ) {
  //	case NO_STRUCTURE:  out << "NO_STRUCTURE";  break;
  //	case UNSYMMETRIC:   out << "UNSYMMETRIC";   break;
  //	case SYMMETRIC:     out << "SYMMETRIC";     break;
  //	default:            out << "STRUCTURE not initialized";
  //	};
  //	return out;
  //}


  //inline std :: ostream & operator<<( std :: ostream & out, const FUNC & f )
  //{
  //	switch( f ) {
  //	case NO_FUNC:    out << "NO_FUNC";    break;
  //	case CONSTANT:   out << "CONSTANT";   break;
  //	case LINEAR:     out << "LINEAR";     break;
  //	case QUADRATIC:  out << "QUADRATIC";  break;
  //	case SEMILINEAR: out << "SEMILINEAR"; break;
  //	case BILINEAR:   out << "BILINEAR";   break;
  //	case NONLINEAR:  out << "NONLINEAR";  break;
  //	default:         out << "FUNC not initialized";
  //	}
  //	return out;
  //}


  inline std :: ostream & operator<<( std :: ostream & out,
                                      const POTENTIAL & mat )
  {
    switch( mat ) {
    case NO_POTENTIAL:    out << "Potential Type not specified yet"; break;
    case SINGLELAYER:     out << "Single Layer Potential"; break;
    case DOUBLELAYER:     out << "Double Layer Potential"; break;
    case HYPEROPERATOR:   out << "Hyper singular operator"; break;
    case NEWTONPOTENTIAL: out << "Newton potenital"; break;
    case IDENTITY:        out << "Identity"; break;
    default:              out << "Potential type not defined"; break;
    }
    return out;
  }


  //inline std :: ostream & operator<<( std :: ostream & out,
  //																		const H_MAT_TYPE & h_mat )
  //{
  //	switch( h_mat ) {
  //	case H_STA:   out << "unsymmetric H - matrix"; break;
  //	case H_SYM:   out << "symmetric H - matrix"; break;
  //	case H_LTR:   out << "lower triangle H - matrix"; break;
  //	case H_UTR:   out << "upper triangle H - matrix"; break;
  //	case H_TST:   out << "transposed standart  H - matrix"; break;
  //	case H_TSY:   out << "transposed symmetric H - matrix"; break;
  //	case H_MAS:   out << "mass H - matrix"; break;
  //	case NO_H:    out << "no H - matrix"; break;
  //	default:      out << "H - matrix type not defined"; break;
  //	}
  //	return out;
  //}


  //inline std :: ostream & operator<<( std :: ostream & out,
  //																		const TIME_FUN & tf )
  //{
  //	switch( tf ) {
  //	case HEAVISIDE: out << "Heaviside function";    break;
  //	case SINE:      out << "Sine wave";             break;
  //	case IMPULSE:   out << "Impulse";               break;
  //	case SAWTOOTH:  out << "Auf und ab";            break;
  //	case RELAX:     out << "Exponential Decay";     break;
  //	case SIGMOIDAL: out << "Sigmoidal function";    break;
  //	default:        out << "Time Function not set"; break;
  //	}
  //	return out;
  //}

  //inline std :: ostream & operator<<( std::ostream &out, const FUNDSOL &f ) {
  //	switch( f ) {
  //	case NO_FUNDSOL:         out << "NO_FUNDSOL";              break;
  //	case FUND_LAPLACE:       out << "LAPLACE";                 break;
  //	case FUND_HELMHOLTZ:     out << "HELMHOLTZ";               break;
  //	case FUND_WAVE:          out << "WAVE";                    break;
  //	case FUND_ELASTOSTATIC:  out << "ELASTOSTATIC";            break;
  //	case FUND_ELASTODYNAMIC: out << "ELASTODYNAMIC";           break;
  //	case FUND_VISCO:         out << "VISCO";                   break;
  //	case FUND_PORO:          out << "PORO";                    break;
  //	default:                 out << "FUNDSOL not initialized"; break;
  //	}
  //	return out;
  //}

  inline std :: ostream & operator<<( std::ostream &out, const TOPOLOGY &t ) {
    switch( t ) {
    case NO_TOPOLOGY: out << "NO_TOPOLOGY";              break;
    case OPEN:        out << "OPEN";                     break;
    case CLOSED:      out << "CLOSED";                   break;
    default:          out << "TOPOLOGY not initialized"; break;
    }
    return out;
  }

  //inline std :: ostream & operator<<( std::ostream &out,const FASTMETHOD &t )
  //{
  //	switch( t ) {
  //	case NO_FASTMETHOD: out << "NO_FASTMETHOD  ";                break;
  //	case ACA:           out << "ACA";                            break;
  //	default:            out << "FASTMETHOD not initialized"; break;
  //	}
  //	return out;
  //}

  //inline std :: ostream & operator<<( std::ostream &out, const ELEM_TYPE &e ) {
  //	switch( e ) {
  //	case NO_ELEM_TYPE:   out << "NO_ELEM_TYPE";              break;
  //	case LINE:           out << "LINE";                      break;
  //	case TRIANGLE:       out << "TRIANGLE";                  break;
  //	case TRIA3:          out << "TRIA3";                     break;
  //	case TRIA6:          out << "TRIA6";                     break;
  //	case QUADRANGLE:     out << "QUADRANGLE";                break;
  //	case QUAD4:          out << "QUAD4";                     break;
  //	case QUAD9:          out << "QUAD9";                     break;
  //	case INFINITE_LINE:  out << "INFINITE_LINE";             break;
  //	case INFINITE_QUAD:  out << "INFINITE_QUAD";             break;
  //	default:             out << "ELEM_TYPE not initialized"; break;
  //	}
  //	return out;
  //}

  //inline std :: ostream & operator<<( std::ostream &out, const SUPPORT &s ) {
  //	switch( s ) {
  //	case NO_SUPPORT:   out << "NO_SUPPORT";              break;
  //	case IS_SUPPORT:   out << "IS_SUPPORT";              break;
  //	case PRIMAL:       out << "PRIMAL";                  break;
  //	case DUAL:         out << "DUAL";                    break;
  //	default:           out << "SUPPORT not initialized"; break;
  //	}
  //	return out;
  //}

  //inline std :: ostream & operator<<( std::ostream &out, const SIGN &s ) {
  //	switch( s ) {
  //	case NO_SIGN:      out << "NO_SIGN";               break;
  //	case MINUS:        out << "MINUS";                 break;
  //	case PLUS:         out << "PLUS";                   break;
  //	default:           out << "SIGN not initialized";  break;
  //	}
  //	return out;
  //}

  //inline std :: ostream & operator<<( std::ostream &out,
  //																		const NEWTON_POTENTIAL &np ) {
  //	switch( np ) {
  //	case NO_NEWTON_POTENTIAL: out << "NO_NEWTON_POTENTIAL";              break;
  //	case POINT_LOAD:          out << "POINT_LOAD";                       break;
  //	case DOMAIN_LOAD:         out << "DOMAIN_LOAD";                      break;
  //	default:                  out << "NEWTON_POTENTIAL not initialized"; break;
  //	}
  //	return out;
  //}

  //! input operator for the approximation order enumerators
  inline std :: istream & operator>>( std :: istream & is, APPROXIMATION & f )
  {
    std :: string buffer;
    is >> buffer;
    if      ( buffer.find( "CONST"   ) != std::string :: npos ) f = CONSTANT;
    else if ( buffer.find( "LIN"     ) != std::string :: npos ) f = LINEAR;
    else if ( buffer.find( "QUA"     ) != std::string :: npos ) f = QUADRATIC;
    else f = NO_APPROX_ORDER;
    return is;
  }

} // end namespace hyena

#endif // enumerators_h
