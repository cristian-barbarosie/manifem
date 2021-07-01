
// finite-elem.h 2021.04.10

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
                          quad_4, quad_9 };
	struct FromFiniteElementWithMaster { };
	static const FromFiniteElementWithMaster from_finite_element_with_master;
	struct ThroughDockedFiniteElement { };
	static const ThroughDockedFiniteElement through_docked_finite_element;
}

class FiniteElement;

//-----------------------------------------------------------------------------------------//

class Integrator

// a wrapper for Gauss quadrature
// and perhaps other integration methods (e.g. symbolic integration)

{	public :

	class Core;  class Gauss;
	
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
	
};  // end of  class Integrator::Core

//-----------------------------------------------------------------------------------------//

inline double Integrator::operator()
( const Function & f, const tag::ThroughDockedFiniteElement &, const FiniteElement & fe )

{	return this->core->action ( f, fe );  }

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
	
};  // end of  class Integrator::Gauss

//-----------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------//

// Again ...
// What is a finite element ?

class FiniteElement

// wrapper for different types of finite elements	
	
{	public :

	class Core;  class WithMaster;  class Lagrange;
	
	FiniteElement::Core * core;

	// constructor

	inline FiniteElement ( const tag::WithMaster &, const tag::Segment &,
                         const tag::lagrange &, const tag::OfDegree &, size_t deg );
	inline FiniteElement ( const tag::WithMaster &, const tag::Triangle &,
                         const tag::lagrange &, const tag::OfDegree &, size_t deg );
	inline FiniteElement ( const tag::WithMaster &, const tag::Quadrangle &,
                         const tag::lagrange &, const tag::OfDegree &, size_t deg );

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

	inline void dock_on ( const Cell & cll );

	inline double integrate ( const Function & );
	
};  // end of  class FiniteElement
	
//-----------------------------------------------------------------------------------------//

inline Integrator::Integrator ( const tag::gauss &, const tag::gauss_quadrature & q,
                                const tag::FromFiniteElementWithMaster &, FiniteElement & fe )
:	core { new Integrator::Gauss ( q, tag::from_finite_element_with_master, fe ) }  { }

//-----------------------------------------------------------------------------------------//

class FiniteElement::Core

// a base class for FiniteElements

{	public :

	Integrator integr;

	// constructor

	inline Core () : integr ( tag::non_existent ) { };

	// destructor

	virtual ~Core()  { };

	virtual void dock_on ( const Cell & cll ) = 0;
	
};  // end of  class FiniteElement::Core

//-----------------------------------------------------------------------------------------//

inline FiniteElement::~FiniteElement ()  {  delete this->core;  }

inline void FiniteElement::dock_on ( const Cell & cll )
{	this->core->dock_on ( cll );  }

//-----------------------------------------------------------------------------------------//

class FiniteElement::WithMaster : public FiniteElement::Core

// finite elements which use a master element

{	public :

	Manifold master_manif;
	Function transf;
	std::map < Cell::Core *, Function > base_fun_1;
	std::map < Cell::Core *, std::map < Cell::Core *, Function > > base_fun_2;

	// at the beginning, 'docked_on' is a non-existent cell
	// later, we shall dock the finite element on a physical cell
	Cell docked_on;

	// constructor

	inline WithMaster ( Manifold m )
	: FiniteElement::Core (), master_manif (m), transf ( 0. ), docked_on ( tag::non_existent )
	{ }

	// get basis functions associated to vertices, to segments, etc
	inline Function basis_function ( Cell::Core * cll );
	inline Function basis_function ( Cell::Core * c1, Cell::Core * c2 );

	void dock_on ( const Cell & cll ) = 0;  // virtual from FiniteElement::Core

	class Segment;  class Triangle;  class Quadrangle;
	
};  // end of  class FiniteElement::withMaster

//-----------------------------------------------------------------------------------------//

inline Function FiniteElement::WithMaster::basis_function ( Cell::Core * cll )
{	std::map < Cell::Core *, Function > :: iterator it =
		this->base_fun_1.find ( cll );
	assert ( it != this->base_fun_1.end() );
  return it->second;                                           }
	
inline Function FiniteElement::WithMaster::basis_function ( Cell::Core * c1, Cell::Core * c2 )
{	std::map < Cell::Core *, std::map < Cell::Core *, Function > >
		:: iterator it = this->base_fun_2.find ( c1 );
	assert ( it != this->base_fun_2.end() );
	std::map < Cell::Core *, Function > :: iterator itt =
		it->second.find ( c2 );
	assert ( itt != it->second.end() );
  return itt->second;                                                          }

inline Function FiniteElement::basis_function ( const Cell cll )
{	FiniteElement::WithMaster * fe_core =
		dynamic_cast < FiniteElement::WithMaster * > ( this->core );
	assert ( fe_core );
	return fe_core->basis_function ( cll.core );                        }

inline Function FiniteElement::basis_function ( const Cell c1, const Cell c2 )
{	FiniteElement::WithMaster * fe_core =
		dynamic_cast < FiniteElement::WithMaster * > ( this->core );
	assert ( fe_core );
	return fe_core->basis_function ( c1.core, c2.core );                }

inline double FiniteElement::integrate ( const Function & f )
{	// assert that 'this' is already docked :
	FiniteElement::WithMaster * fe_core =
		dynamic_cast < FiniteElement::WithMaster * > ( this->core );
	assert ( fe_core );
	assert ( fe_core->docked_on.exists() );
	return this->core->integr ( f, tag::through_docked_finite_element, *this );  }

//-----------------------------------------------------------------------------------------//

class FiniteElement::WithMaster::Segment : public FiniteElement::WithMaster

// triangular finite elements which use a master element

{	public :

	// constructor

	inline Segment ( Manifold m ) : FiniteElement::WithMaster ( m )  { }
	
	void dock_on ( const Cell & cll );
	// virtual from FiniteElement::Core

};  // end of  class FiniteElement::WithMaster::Segment

//-----------------------------------------------------------------------------------------//

class FiniteElement::WithMaster::Triangle : public FiniteElement::WithMaster

// triangular finite elements which use a master element

{	public :

	// constructor

	inline Triangle ( Manifold m ) : FiniteElement::WithMaster ( m )  { }
	
	void dock_on ( const Cell & cll );
	// virtual from FiniteElement::Core

};  // end of  class FiniteElement::withMaster::Triangle

//-----------------------------------------------------------------------------------------//

class FiniteElement::WithMaster::Quadrangle : public FiniteElement::WithMaster

// quadrangular finite elements which use a master element

{	public :

	// constructor

	inline Quadrangle ( Manifold m ) : FiniteElement::WithMaster ( m )  { }
	
	void dock_on ( const Cell & cll );
	// virtual from FiniteElement::Core

};  // end of  class FiniteElement::withMaster::Quadrangle

//-----------------------------------------------------------------------------------------//


inline FiniteElement::FiniteElement
(	const tag::WithMaster &, const tag::Segment &,
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
(	const tag::WithMaster &, const tag::Triangle &,
	const tag::lagrange &, const tag::OfDegree &, size_t deg )
:	core { nullptr }

{	assert ( deg == 1 );

	// we keep the working manifold and restaure it at the end
	Manifold work_manif = Manifold::working;
	
	Manifold RR2_master ( tag::Euclid, tag::of_dim, 2 );
	Function xi_eta = RR2_master.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	// we should take advantage of the memory space already reserved for x and y
	// Function xi = xi_eta[0], eta = xi_eta[1];

	this->core = new FiniteElement::WithMaster::Triangle ( RR2_master );
	
	work_manif.set_as_working_manifold();                                                       }


inline FiniteElement::FiniteElement
(	const tag::WithMaster &, const tag::Quadrangle &,
	const tag::lagrange &, const tag::OfDegree &, size_t deg )
:	core { nullptr }

{	assert ( deg == 1 );

	// we keep the working manifold and restaure it at the end
	Manifold work_manif = Manifold::working;
	
	Manifold RR2_master ( tag::Euclid, tag::of_dim, 2 );
	Function xi_eta = RR2_master.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	// we should take advantage of the memory space already reserved for x and y
	// Function xi = xi_eta[0], eta = xi_eta[1];

	this->core = new FiniteElement::WithMaster::Quadrangle ( RR2_master );
	
	work_manif.set_as_working_manifold();                                                       }


inline Integrator FiniteElement::set_integrator
( const tag::gauss &, const tag::gauss_quadrature & q )
	
{ FiniteElement::WithMaster * this_core =
		dynamic_cast < FiniteElement::WithMaster * > ( this->core );
	assert ( this_core );
	this_core->integr.core =
		new Integrator::Gauss ( q, tag::from_finite_element_with_master, *this );
	return this_core->integr;                                                    }


}  // end of  namespace maniFEM

#endif
// ifndef MANIFEM_FINITE_ELEM_H
