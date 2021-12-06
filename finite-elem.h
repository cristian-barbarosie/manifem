
// finite-elem.h 2021.12.06

//   This file is part of maniFEM, a C++ library for meshes and finite elements on manifolds.

//   Copyright 2019, 2020, 2021 Cristian Barbarosie cristian.barbarosie@gmail.com
//   https://github.com/cristian-barbarosie/manifem

//   ManiFEM is free software: you can redistribute it and/or modify it
//   under the terms of the GNU Lesser General Public License as published
//   by the Free Software Foundation, either version 3 of the License
//   or (at your option) any later version.

//   ManiFEM is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//   See the GNU Lesser General Public License for more details.

//   You should have received a copy of the GNU Lesser General Public License
//   along with maniFEM.  If not, see <https://www.gnu.org/licenses/>.

#ifndef MANIFEM_FINITE_ELEM_H
#define MANIFEM_FINITE_ELEM_H

#include <iostream>
#include <vector>
#include <map>
#include <forward_list>
#include "assert.h"

#include "mesh.h"
#include "field.h"
#include "function.h"
#include "manifold.h"

namespace maniFEM {

namespace tag {
	struct On { };  static const On on;
	struct gauss { };  static const gauss Gauss;
	struct WithMaster { };  static const WithMaster with_master;
	enum gauss_quadrature { seg_2, seg_3, seg_4, seg_5, seg_6,
                          tri_3, tri_3_Oden, tri_4, tri_4_Oden, tri_6,
                          quad_4, quad_9                              };
	struct FromFiniteElementWithMaster { };
	static const FromFiniteElementWithMaster from_finite_element_with_master;
	struct FromFiniteElement { };  static const FromFiniteElement from_finite_element;
	struct ThroughDockedFiniteElement { };
	static const ThroughDockedFiniteElement through_docked_finite_element;
	struct EnumerateCells { };  static const EnumerateCells enumerate_cells;
	struct ForGiven { };  static const ForGiven for_given;
	struct ForAGiven { };  static const ForAGiven for_a_given;
	struct BasisFunctions { };  static const BasisFunctions basis_functions;
	struct IntegralOf { };  static const IntegralOf integral_of;
	struct UFL_FFC { };  static const UFL_FFC ufl_ffc;
}

class FiniteElement;

//-----------------------------------------------------------------------------------------//

	
class Integrator

// a wrapper for Gauss quadrature
// and perhaps other integration methods (e.g. symbolic integration)

{	public :

	class Core;  class Gauss;  class UFL_FFC;
	
	Integrator::Core * core;

	// constructor

	inline Integrator ( const tag::NonExistent & )	: core { nullptr }  { }
	inline Integrator ( const tag::gauss &, const tag::gauss_quadrature & q,
	                    const tag::FromFiniteElementWithMaster &, FiniteElement & fe );

	// operator() integrates a given function on a given mesh or cell
	
	inline double operator() ( const Function & f, const tag::On &, const Cell cll );
	inline double operator() ( const Function & f, const tag::On &, const Mesh msh );

	// there is also an operator() without a domain,
	// for integrators attached to a finite element,
	// in the event that the finite element is already docked on a cell

	// we assume the function f is already expressed in terms of master coordinates
	// (it is composed with fe.core->transf)
	// (this is necessary for Gauss quadrature, not for symbolic integration)

	inline double operator()
	( const Function & f, const tag::ThroughDockedFiniteElement &, const FiniteElement & fe );

	// UFL_FFC integrators require a 'pre_compute' step
	inline void pre_compute
	( const tag::ForAGiven &, const tag::BasisFunction &, Function bf,
	  tag::IntegralOf &, std::vector < Function > &                   );
	inline void pre_compute
	( const tag::ForGiven &, const tag::BasisFunctions &, Function bf1, Function bf2,
	  tag::IntegralOf &, std::vector < Function > &                                  );
														
	
};  // end of  class Integrator


//-----------------------------------------------------------------------------------------//

class Integrator::Core

// a base class for Gauss quadrature
// and perhaps other integration methods (e.g. symbolic integration)

{	public :

	// constructor

	inline Core ()  { };

	// destructor

	virtual ~Core()  { };

	// operator()
	
	virtual double action ( Function f, const FiniteElement & fe ) = 0;

	// 'pre_compute' is only meaningful for UFL_FFC integrators
	virtual void pre_compute ( Function bf, std::vector < Function > & v ) = 0;
	virtual void pre_compute ( Function bf1, Function bf2, std::vector < Function > & v ) = 0;
	
};  // end of  class Integrator::Core

//-----------------------------------------------------------------------------------------//

inline double Integrator::operator()
( const Function & f, const tag::ThroughDockedFiniteElement &, const FiniteElement & fe )

{	return this->core->action ( f, fe );  }


inline void Integrator::pre_compute  // only meaningful for UFL_FFC integrators
( const tag::ForAGiven &, const tag::BasisFunction &, Function bf,
  tag::IntegralOf &, std::vector < Function > & v                 )

// bf is an arbitrary basis function (Function::MereSymbol) in the finite element
// we prepare computations for fast evaluation of integrals of expressions listed in 'v'

{	this->core->pre_compute ( bf, v );  }
	

inline void Integrator::pre_compute  // only meaningful for UFL_FFC integrators
( const tag::ForGiven &, const tag::BasisFunctions &, Function bf1, Function bf2,
  tag::IntegralOf &, std::vector < Function > & v                                )

// bf1 and bf2 are arbitrary basis functions (Function::MereSymbol) in the finite element
// we prepare computations for fast evaluation of integrals of expressions listed in 'v'

{	this->core->pre_compute ( bf1, bf2, v );  }

//-----------------------------------------------------------------------------------------//

	
class Integrator::Gauss : public Integrator::Core

{	public :

	std::vector < Cell > points;
	std::vector < double > weights;
	
	// constructor

	inline Gauss ( const tag::gauss_quadrature & q ) : Integrator::Core ()
	{ switch ( q )
		{	case tag::tri_3 : std::cout << "tri 3" << std::endl;  break;
			case tag::tri_4 : std::cout << "tri 4" << std::endl;  break;
			case tag::tri_6 : std::cout << "tri 6" << std::endl;  break;
			default: std::cout << "default" << std::endl;                  }
	}

	Gauss ( const tag::gauss_quadrature & q,
	        const tag::FromFiniteElementWithMaster &, FiniteElement & fe );

	// operator()
	
	double action ( Function f, const FiniteElement & fe );
	// virtual from Integrator::Core
	
	//  pre_compute  is virtual from Integrator::Core, here execution forbidden
	void pre_compute ( Function bf, std::vector < Function > & v );
	void pre_compute ( Function bf1, Function bf2, std::vector < Function > & v );

};  // end of  class Integrator::Gauss


//-----------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------//

	
// Again ...
// What is a finite element ?

class FiniteElement

// wrapper for different types of finite elements	
	
{	public :

	class Core;  class WithMaster;  class StandAlone;
	
	FiniteElement::Core * core;

	std::vector < Cell::Numbering * > numbers;

	// constructor

	inline FiniteElement ( const tag::Triangle &,  // no master, stand-alone
                         const tag::lagrange &, const tag::OfDegree &, size_t deg );
	inline FiniteElement ( const tag::WithMaster &, const tag::Segment &,
                         const tag::lagrange &, const tag::OfDegree &, size_t deg );
	inline FiniteElement ( const tag::WithMaster &, const tag::Segment &,
	                       const tag::lagrange &, const tag::OfDegree &, size_t deg,
	                       const tag::EnumerateCells &                               );
	inline FiniteElement ( const tag::WithMaster &, const tag::Triangle &,
                         const tag::lagrange &, const tag::OfDegree &, size_t deg );
	inline FiniteElement ( const tag::WithMaster &, const tag::Triangle &,
                         const tag::lagrange &, const tag::OfDegree &, size_t deg,
	                       const tag::EnumerateCells &                               );
	inline FiniteElement ( const tag::WithMaster &, const tag::Quadrangle &,
                         const tag::lagrange &, const tag::OfDegree &, size_t deg );
	inline FiniteElement ( const tag::WithMaster &, const tag::Quadrangle &,
                         const tag::lagrange &, const tag::OfDegree &, size_t deg,
	                       const tag::EnumerateCells &                               );

	// destructor
	
	inline ~FiniteElement ();
	
	// forbid copying

	FiniteElement ( const FiniteElement & ) = delete;
	FiniteElement ( const FiniteElement && ) = delete;
	FiniteElement & operator= ( const FiniteElement & ) = delete;
	FiniteElement & operator= ( const FiniteElement && ) = delete;

	// get basis functions associated to vertices, to segments, etc
	inline Function basis_function ( const Cell cll );
	inline Function basis_function ( const Cell c1, const Cell c2 );

	inline Integrator set_integrator ( const tag::gauss &, const tag::gauss_quadrature & );
	inline Integrator set_integrator ( const tag::UFL_FFC & );

	inline void dock_on ( const Cell & cll );
	inline void dock_on ( const Cell & cll, const tag::Winding & );

	inline double integrate ( const Function & );

	inline void pre_compute
	( const tag::ForAGiven &, const tag::BasisFunction &, Function bf, 
	  const tag::IntegralOf &, std::vector < Function > v             );
	inline void pre_compute
	( const tag::ForGiven &, const tag::BasisFunctions &, Function bf1, Function bf2, 
	  const tag::IntegralOf &, std::vector < Function > v                            );

	inline Cell::Numbering & numbering ( const tag::Vertices & );
	inline Cell::Numbering & numbering ( const tag::Segments & );
	inline Cell::Numbering & numbering ( const tag::CellsOfDim &, const size_t d );
	
};  // end of  class FiniteElement
	
//-----------------------------------------------------------------------------------------//


class Integrator::UFL_FFC : public Integrator::Core

// sort of symbolic integrator, obtained using UFL and FFC
// https://fenics.readthedocs.io/projects/ufl/en/latest/
// https://fenics.readthedocs.io/projects/ffc/en/latest/
// resulting code manipulated by hand for each case
// TODO : further optimize the code using symbolic computation and genetic/stochastic optimization

{	public :

	// this type of integrator is tightly linked to a finite element,
	// so we include it as an attribute
	FiniteElement fe;

	//constructor
	
	UFL_FFC ( const tag::FromFiniteElement &, FiniteElement & fe );

	// operator()
	
	double action ( Function f, const FiniteElement & fe );
	// virtual from Integrator::Core
	
	// this type of integrator benefits from an early declaration of
	// the integrals we intend to compute later (after docking on a cell)
	//  pre_compute  is virtual from Integrator::Core
	void pre_compute ( Function bf, std::vector < Function > & v );
	void pre_compute ( Function bf1, Function bf2, std::vector < Function > & v );

};  // end of  class Integrator::UFL_FFC

//-----------------------------------------------------------------------------------------//

inline Integrator::Integrator ( const tag::gauss &, const tag::gauss_quadrature & q,
                                const tag::FromFiniteElementWithMaster &, FiniteElement & fe )
:	core { new Integrator::Gauss ( q, tag::from_finite_element_with_master, fe ) }  { }

//-----------------------------------------------------------------------------------------//


class FiniteElement::Core

// a base class for FiniteElements

{	public :

	Integrator integr;

	std::map < Cell::Core *, Function > base_fun_1;
	std::map < Cell::Core *, std::map < Cell::Core *, Function > > base_fun_2;

	// at the beginning, 'docked_on' is a non-existent cell
	// later, we shall dock the finite element on a physical cell
	Cell docked_on;

	// constructor

	inline Core () : integr ( tag::non_existent ), docked_on ( tag::non_existent )  { };

	// destructor

	virtual ~Core()  { };

	virtual void dock_on ( const Cell & cll ) = 0;
	virtual void dock_on ( const Cell & cll, const tag::Winding & ) = 0;
	
	virtual void pre_compute ( const Function bf, std::vector < Function > result ) = 0;
	virtual void pre_compute ( const Function bf1, const Function bf2,
                             std::vector < Function > result        ) = 0;

};  // end of  class FiniteElement::Core

//-----------------------------------------------------------------------------------------//


inline FiniteElement::~FiniteElement ()  {  delete this->core;  }


inline void FiniteElement::dock_on ( const Cell & cll )
{	this->core->dock_on ( cll );  }

inline void FiniteElement::dock_on ( const Cell & cll, const tag::Winding & )
{	this->core->dock_on ( cll, tag::winding );  }


inline Cell::Numbering & FiniteElement::numbering ( const tag::Vertices & )
{	assert ( this->numbers .size() );
	return * this->numbers [0];       }

inline Cell::Numbering & FiniteElement::numbering ( const tag::Segments & )
{	assert ( this->numbers .size() > 1 );
	return * this->numbers [1];           }

inline Cell::Numbering & FiniteElement::numbering ( const tag::CellsOfDim &, const size_t d )
{	assert ( this->numbers .size() > d );
	return * this->numbers [d];           }

//-----------------------------------------------------------------------------------------//


class FiniteElement::WithMaster : public FiniteElement::Core

// finite elements which use a master element

{	public :

	// attributes inherited from FiniteElement::Core :
	// Integrator integr
	// std::map < Cell::Core *, Function > base_fun_1
	// std::map < Cell::Core *, std::map < Cell::Core *, Function > > base_fun_2
	// Cell docked_on

	Manifold master_manif;
	Function transf;

	// constructor

	inline WithMaster ( Manifold m )
	: FiniteElement::Core (), master_manif (m), transf ( 0. )
	{ }

	// get basis functions associated to vertices, to segments, etc
	inline Function basis_function ( Cell::Core * cll );
	inline Function basis_function ( Cell::Core * c1, Cell::Core * c2 );

	// dock_on  virtual from FiniteElement::Core
	void dock_on ( const Cell & cll ) = 0;
	void dock_on ( const Cell & cll, const tag::Winding & ) = 0;

	//  pre_compute  virtual from FiniteElement::Core, here execution forbidden
	void pre_compute ( const Function bf, std::vector < Function > result );
	void pre_compute ( const Function bf1, const Function bf2,
                     std::vector < Function > result        );

	class Segment;  class Triangle;  class Quadrangle;
	
};  // end of  class FiniteElement::withMaster

//-----------------------------------------------------------------------------------------//

inline Function FiniteElement::WithMaster::basis_function ( Cell::Core * cll )
{	std::map < Cell::Core *, Function > :: iterator it =
		this->base_fun_1 .find ( cll );
	assert ( it != this->base_fun_1 .end() );
  return it->second;                                   }
	
inline Function FiniteElement::WithMaster::basis_function ( Cell::Core * c1, Cell::Core * c2 )
{	std::map < Cell::Core *, std::map < Cell::Core *, Function > >
		:: iterator it = this->base_fun_2 .find ( c1 );
	assert ( it != this->base_fun_2 .end() );
	std::map < Cell::Core *, Function > :: iterator itt =
		it->second .find ( c2 );
	assert ( itt != it->second .end() );
  return itt->second;                                             }

inline Function FiniteElement::basis_function ( const Cell cll )
{	FiniteElement::WithMaster * fe_core =
		dynamic_cast < FiniteElement::WithMaster * > ( this->core );
	assert ( fe_core );
	return fe_core->basis_function ( cll .core );                        }

inline Function FiniteElement::basis_function ( const Cell c1, const Cell c2 )
{	FiniteElement::WithMaster * fe_core =
		dynamic_cast < FiniteElement::WithMaster * > ( this->core );
	assert ( fe_core );
	return fe_core->basis_function ( c1 .core, c2 .core );          }

inline double FiniteElement::integrate ( const Function & f )
{	// assert that 'this' is already docked :
	FiniteElement::WithMaster * fe_core =
		dynamic_cast < FiniteElement::WithMaster * > ( this->core );
	assert ( fe_core );
	assert ( fe_core->docked_on .exists() );
	return this->core->integr ( f, tag::through_docked_finite_element, *this );  }

inline void FiniteElement::pre_compute
( const tag::ForAGiven &, const tag::BasisFunction &, Function bf,
  const tag::IntegralOf &, std::vector < Function > v             )
{	this->core->pre_compute ( bf, v );  }

inline void FiniteElement::pre_compute
( const tag::ForGiven &, const tag::BasisFunctions &, Function bf1, Function bf2,
  const tag::IntegralOf &, std::vector < Function > v                            )
{	this->core->pre_compute ( bf1, bf2, v );  }

//-----------------------------------------------------------------------------------------//


class FiniteElement::WithMaster::Segment : public FiniteElement::WithMaster

// triangular finite elements which use a master element

{	public :

	// attributes inherited from FiniteElement::Core :
	// Integrator integr
	// std::map < Cell::Core *, Function > base_fun_1
	// std::map < Cell::Core *, std::map < Cell::Core *, Function > > base_fun_2
	// Cell docked_on

	// attributes inherited from FiniteElement::WithMaster :
	// Manifold master_manif
	// Function transf

	// constructor

	inline Segment ( Manifold m ) : FiniteElement::WithMaster ( m )  { }
	
	// dock_on  virtual from FiniteElement::Core
	void dock_on ( const Cell & cll );
	void dock_on ( const Cell & cll, const tag::Winding & );

	// two versions of pre_compute virtual from FiniteElement::Core,
	// defined by FiniteElement::WithMaster, execution forbidden

};  // end of  class FiniteElement::WithMaster::Segment

//-----------------------------------------------------------------------------------------//

class FiniteElement::WithMaster::Triangle : public FiniteElement::WithMaster

// triangular finite elements which use a master element

{	public :

	// attributes inherited from FiniteElement::Core :
	// Integrator integr
	// std::map < Cell::Core *, Function > base_fun_1
	// std::map < Cell::Core *, std::map < Cell::Core *, Function > > base_fun_2
	// Cell docked_on

	// attributes inherited from FiniteElement::WithMaster :
	// Manifold master_manif
	// Function transf

	// constructor

	inline Triangle ( Manifold m ) : FiniteElement::WithMaster ( m )  { }
	
	// dock_on  virtual from FiniteElement::Core
	void dock_on ( const Cell & cll );
	void dock_on ( const Cell & cll, const tag::Winding & );

	// two versions of pre_compute virtual from FiniteElement::Core,
	// defined by FiniteElement::WithMaster, execution forbidden

};  // end of  class FiniteElement::withMaster::Triangle

//-----------------------------------------------------------------------------------------//

class FiniteElement::WithMaster::Quadrangle : public FiniteElement::WithMaster

// quadrangular finite elements which use a master element

{	public :

	// attributes inherited from FiniteElement::Core :
	// Integrator integr
	// std::map < Cell::Core *, Function > base_fun_1
	// std::map < Cell::Core *, std::map < Cell::Core *, Function > > base_fun_2
	// Cell docked_on

	// attributes inherited from FiniteElement::WithMaster :
	// Manifold master_manif
	// Function transf

	// constructor

	inline Quadrangle ( Manifold m ) : FiniteElement::WithMaster ( m )  { }
	
	// dock_on  virtual from FiniteElement::Core
	void dock_on ( const Cell & cll );
	void dock_on ( const Cell & cll, const tag::Winding & );

	// two versions of pre_compute virtual from FiniteElement::Core,
	// defined by FiniteElement::WithMaster, execution forbidden

};  // end of  class FiniteElement::withMaster::Quadrangle

//-----------------------------------------------------------------------------------------//


inline FiniteElement::FiniteElement
( const tag::WithMaster &, const tag::Segment &,
  const tag::lagrange &, const tag::OfDegree &, size_t deg )
	
:	core { nullptr }

{	assert ( deg == 1 );

	// we keep the working manifold and restore it at the end
	Manifold work_manif = Manifold::working;
	
	Manifold RR_master ( tag::Euclid, tag::of_dim, 1 );
	Function t = RR_master.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	// we should take advantage of the memory space already reserved for x and y

	this->core = new FiniteElement::WithMaster::Segment ( RR_master );
	
	work_manif.set_as_working_manifold();                                                   }


inline FiniteElement::FiniteElement
( const tag::WithMaster &, const tag::Segment &,
  const tag::lagrange &, const tag::OfDegree &, size_t deg,
  const tag::EnumerateCells &                               )
	
:	core { nullptr }

{	assert ( deg == 1 );

	// we keep the working manifold and restore it at the end
	Manifold work_manif = Manifold::working;
	
	Manifold RR_master ( tag::Euclid, tag::of_dim, 1 );
	Function t = RR_master .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	// we should take advantage of the memory space already reserved for x and y

	// std::vector < Cell::Numbering::Field * > numbers;
	// this->numbers .emplace ( this->numbers .end(), tag::vertices );
	this->numbers .push_back ( new Cell::Numbering::Field ( tag::vertices ) );

	this->core = new FiniteElement::WithMaster::Segment ( RR_master );
	
	work_manif .set_as_working_manifold();                                                 }


inline FiniteElement::FiniteElement
(	const tag::WithMaster &, const tag::Triangle &,
	const tag::lagrange &, const tag::OfDegree &, size_t deg )
	
:	core { nullptr }

{	assert ( deg == 1 );

	// we keep the working manifold and restaure it at the end
	Manifold work_manif = Manifold::working;
	
	Manifold RR2_master ( tag::Euclid, tag::of_dim, 2 );
	Function xi_eta = RR2_master .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	// we should take advantage of the memory space already reserved for x and y
	// Function xi = xi_eta [0], eta = xi_eta [1];

	this->core = new FiniteElement::WithMaster::Triangle ( RR2_master );
	
	work_manif .set_as_working_manifold();                                                      }


inline FiniteElement::FiniteElement
( const tag::WithMaster &, const tag::Triangle &,
  const tag::lagrange &, const tag::OfDegree &, size_t deg,
  const tag::EnumerateCells &                               )
	
:	core { nullptr }

{	assert ( deg == 1 );

	// we keep the working manifold and restaure it at the end
	Manifold work_manif = Manifold::working;
	
	Manifold RR2_master ( tag::Euclid, tag::of_dim, 2 );
	Function xi_eta = RR2_master .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	// we should take advantage of the memory space already reserved for x and y
	// Function xi = xi_eta [0], eta = xi_eta [1];

	// std::vector < Cell::Numbering::Field * > numbers;
	// this->numbers .emplace ( this->numbers .end(), tag::vertices );
	this->numbers .push_back ( new Cell::Numbering::Field ( tag::vertices ) );

	this->core = new FiniteElement::WithMaster::Triangle ( RR2_master );
	
	work_manif .set_as_working_manifold();                                                      }


inline FiniteElement::FiniteElement
(	const tag::WithMaster &, const tag::Quadrangle &,
	const tag::lagrange &, const tag::OfDegree &, size_t deg )
	
:	core { nullptr }

{	assert ( deg == 1 );

	// we keep the working manifold and restaure it at the end
	Manifold work_manif = Manifold::working;
	
	Manifold RR2_master ( tag::Euclid, tag::of_dim, 2 );
	Function xi_eta = RR2_master .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	// we should take advantage of the memory space already reserved for x and y
	// Function xi = xi_eta [0], eta = xi_eta [1];

	this->core = new FiniteElement::WithMaster::Quadrangle ( RR2_master );
	
	work_manif .set_as_working_manifold();                                                      }


inline FiniteElement::FiniteElement
( const tag::WithMaster &, const tag::Quadrangle &,
  const tag::lagrange &, const tag::OfDegree &, size_t deg,
  const tag::EnumerateCells &                               )
	
:	core { nullptr }

{	assert ( deg == 1 );

	// we keep the working manifold and restaure it at the end
	Manifold work_manif = Manifold::working;
	
	Manifold RR2_master ( tag::Euclid, tag::of_dim, 2 );
	Function xi_eta = RR2_master .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	// we should take advantage of the memory space already reserved for x and y
	// Function xi = xi_eta [0], eta = xi_eta [1];

	// std::vector < Cell::Numbering::Field * > numbers;
	// this->numbers .emplace ( this->numbers .end(), tag::vertices );
	this->numbers .push_back ( new Cell::Numbering::Field ( tag::vertices ) );

	this->core = new FiniteElement::WithMaster::Quadrangle ( RR2_master );
	
	work_manif .set_as_working_manifold();                                                      }


inline Integrator FiniteElement::set_integrator
( const tag::gauss &, const tag::gauss_quadrature & q )
	
{ FiniteElement::WithMaster * this_core = tag::Util::assert_cast
		< FiniteElement::Core *, FiniteElement::WithMaster * > ( this->core );
	this_core->integr .core =
		new Integrator::Gauss ( q, tag::from_finite_element_with_master, *this );
	return this_core->integr;                                                    }

//-----------------------------------------------------------------------------------------//


class FiniteElement::StandAlone : public FiniteElement::Core

// finite elements with no master

{	public :

	// attributes inherited from FiniteElement::Core :
	// std::map < Cell::Core *, Function > base_fun_1
	// std::map < Cell::Core *, std::map < Cell::Core *, Function > > base_fun_2

	std::map < Cell::Core *, size_t > local_numbering_1;
	std::map < Function::Core *, size_t > basis_numbering_1;

	// at doking, this finite element will perform all computations
	// (previously declared in 'pre_compute') and store the results in 'result_of_integr'
	// the first index identifies a basis function (with the aid of 'basis_numbering')
	// the second index simply identifies the required computation
	std::vector < std::vector < double > > result_of_integr;
	
	// constructor

	inline StandAlone ( ) : FiniteElement::Core ()  { }

	// get basis functions associated to vertices, to segments, etc
	inline Function basis_function ( Cell::Core * cll );
	inline Function basis_function ( Cell::Core * c1, Cell::Core * c2 );

	//  dock_on  virtual from FiniteElement::Core
	void dock_on ( const Cell & cll ) = 0;
	void dock_on ( const Cell & cll, const tag::Winding & ) = 0;

	// two versions of  pre_compute  stay pure virtual from FiniteElement::Core

	class TypeOne;
	
};  // end of  class FiniteElement::StandAlone

//-----------------------------------------------------------------------------------------//


class FiniteElement::StandAlone::TypeOne : public FiniteElement::StandAlone

// finite elements with no master
// searches for the basis function provided by 'integrate'
// uses an std::map < Function, size_t > index_of_basis_function

{	public :

	// attributes inherited from FiniteElement::Core :
	// std::map < Cell::Core *, Function > base_fun_1
	// std::map < Cell::Core *, std::map < Cell::Core *, Function > > base_fun_2

	// attribute inherited from FiniteElement::StandAlone :
	// std::map < Cell::Core *, size_t > local_numbering_1
	// std::map < Function::Core *, size_t > basis_numbering_1

	// constructor

	inline TypeOne ( ) : FiniteElement::StandAlone ()  { }

	// get basis functions associated to vertices, to segments, etc
	inline Function basis_function ( Cell::Core * cll );
	inline Function basis_function ( Cell::Core * c1, Cell::Core * c2 );

	// dock_on  virtual from FiniteElement::Core
	void dock_on ( const Cell & cll ) = 0;
	void dock_on ( const Cell & cll, const tag::Winding & ) = 0;

	// two versions of 'pre_compute' stay pure virtual from FiniteElement::Core

	class Segment;  class Triangle;  class Quadrangle;
	
};  // end of  class FiniteElement::StandAlone::TypeOne

//-----------------------------------------------------------------------------------------//


class FiniteElement::StandAlone::TypeOne::Triangle : public FiniteElement::StandAlone::TypeOne

// triangular finite elements, no master element

{	public :

	// attributes inherited from FiniteElement::Core :
	// Integrator integr
	// std::map < Cell::Core *, Function > base_fun_1
	// std::map < Cell::Core *, std::map < Cell::Core *, Function > > base_fun_2

	// constructor

	inline Triangle ( ) : FiniteElement::StandAlone::TypeOne()  { }
	// define three basis function, mere symbols
	
	// dock_on  virtual from FiniteElement::Core
	void dock_on ( const Cell & cll );
	void dock_on ( const Cell & cll, const tag::Winding & );

	//  pre_compute  virtual from FiniteElement::Core
	void pre_compute ( const Function bf, std::vector < Function > result );
	void pre_compute ( const Function bf1, const Function bf2, std::vector < Function > result );

};  // end of  class FiniteElement::StandAlone::TypeOne::Triangle

//-----------------------------------------------------------------------------------------//


inline FiniteElement::FiniteElement  // no master, stand-alone
( const tag::Triangle &, const tag::lagrange &, const tag::OfDegree &, size_t deg )
	
:	core { nullptr }

{	assert ( deg == 1 );

	this->core = new FiniteElement::StandAlone::TypeOne::Triangle();  }


inline Integrator FiniteElement::set_integrator ( const tag::UFL_FFC & )
	
{ FiniteElement::StandAlone::TypeOne * this_core = tag::Util::assert_cast
		< FiniteElement::Core *, FiniteElement::StandAlone::TypeOne * > ( this->core );
	this_core->integr .core =
		new Integrator::UFL_FFC ( tag::from_finite_element, *this );
	return this_core->integr;                                                         }


}  // end of  namespace maniFEM

#endif
// ifndef MANIFEM_FINITE_ELEM_H
