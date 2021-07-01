
// manifold.h 2021.04.17

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

#ifndef MANIFEM_MANIFOLD_H
#define MANIFEM_MANIFOLD_H

#include <iostream>
#include <vector>
#include <forward_list>
#include "assert.h"

#include "mesh.h"
#include "field.h"
#include "function.h"


namespace maniFEM {

namespace tag {
	struct euclid { };  static const euclid Euclid;
	struct lagrange { };  static const lagrange Lagrange;
	struct Implicit { };  static const Implicit implicit;
	struct Parametric { };  static const Parametric parametric;
	struct DoNotSetAsWorking { };  static const DoNotSetAsWorking do_not_set_as_working;
}


class Manifold

{	public :

	class Core;
	
	Manifold::Core * core;

	inline Manifold ( const tag::WhoseCoreIs &, Manifold::Core * c )
	:	core ( c )
	{	assert ( c );
		Manifold::working = *this;  }

	inline Manifold ( const Manifold & m ) : core ( m.core ) { }

	inline Manifold ( const Manifold && m ) : core ( m.core ) { }
	// is this the right way to do it ?

	inline Manifold ( const tag::NonExistent & ) : core ( nullptr ) { }

	inline Manifold ( const tag::euclid &, const tag::OfDimension &, const size_t dim );

	inline Manifold ( const tag::euclid &, const tag::OfDimension &, const size_t dim,
                    const tag::DoNotSetAsWorking &                                   );

	inline Manifold ( const tag::Implicit &, const Manifold &, const Function & );

	inline Manifold ( const tag::Implicit &,
	                  const Manifold &, const Function &, const Function & );
	
	inline Manifold ( const tag::Parametric &, const Manifold &, const Function::Equality & );

	inline Manifold ( const tag::Parametric &,
	       const Manifold &, const Function::Equality &, const Function::Equality & );
	
	inline Manifold ( const tag::Parametric &, const Manifold &,
	       const Function::Equality &, const Function::Equality &, const Function::Equality & );

	inline ~Manifold ( )
	{} //	std::cout << "destructor Manifold" << std::endl;  }
	
	inline Manifold & operator= ( const Manifold & m )
	{	core = m.core;
		return *this;   }
	
	inline Manifold & operator= ( const Manifold && m )
	{	core = m.core;  // is this the right way to do it ?
		return *this;   }
	
	inline Function coordinates ( ) const;
	inline void set_coordinates ( const Function co );
	// inline void set_coordinates ( Function & ) const;

	inline Function build_coordinate_system
	( const tag::lagrange &, const tag::OfDegree &, size_t d );

	// a non-existent manifold has null core
	inline bool exists ( ) const { return core != nullptr; }

	// metric in the manifold (an inner product on the tangent space)
	inline double inner_prod ( const Cell & P, const std::vector<double> & v,
                                             const std::vector<double> & w ) const;

	// square of the distance, computed from the inner product
	// only meaningful for short distances
	inline double dist_sq ( const Cell & A, const Cell & B ) const
	{	const Function & coord = coordinates();
		std::vector < double > vA = coord ( A ), vB = coord ( B ),
			delta ( coord.nb_of_components() );
		for ( size_t i = 0; i < coord.nb_of_components(); i++ )
			delta[i] = vB[i] - vA[i];
		return inner_prod ( A, delta, delta );                        }
	
	// distance, computed from the inner product
	inline double distance ( const Cell & A, const Cell & B ) const
	{	return std::sqrt ( dist_sq ( A, B ) );  }
	
	// metric in the manifold (an inner product on the tangent space) (metric not used)
	static double default_inner_prod ( const Cell & P, const std::vector<double> & v,
                                     const std::vector<double> & w, const Function & metric );
	
	// metric in the manifold (an inner product on the tangent space) not constant
	static double zoom_inner_prod ( const Cell & P, const std::vector<double> & v,
                                     const std::vector<double> & w, const Function & metric );
	
	// metric in the manifold (an inner product on the tangent space) not constant, anisotropic
	static double matrix_inner_prod ( const Cell & P, const std::vector<double> & v,
                                     const std::vector<double> & w, const Function & metric );

	inline void set_metric ( const Function & metric );
	
	// P = sA + sB,  s+t == 1
	inline void interpolate
	( const Cell & P, double s, const Cell & A, double t, const Cell & B ) const;

	// P = sA + sB + uC + vD,  s+t+u+v == 1
	inline void interpolate ( const Cell & P, double s, const Cell & A,
		double t, const Cell & B, double u, const Cell & C, double v, const Cell & D ) const;

	// P = sA + sB + uC + vD + wE + zF,  s+t+u+v+w+z == 1
	inline void interpolate ( const Cell & P, double s, const Cell & A,
		double t, const Cell & B, double u, const Cell & C, double v, const Cell & D,
		double w, const Cell & E, double z, const Cell & F ) const;

	// P = sum c_k P_k,  sum c_k == 1
	inline void interpolate
		( const Cell & P, std::vector < double > & coefs, std::vector < Cell > & points ) const;

	inline void project ( const Cell & ) const;
	
	inline Manifold implicit ( const Function::Equality eq ) const;
	inline Manifold implicit ( const Function::Equality eq1,
	                           const Function::Equality eq2 ) const;
	inline Manifold parametric ( const Function::Equality eq ) const;
	inline Manifold parametric ( const Function::Equality eq1,
	                             const Function::Equality eq2 ) const;
	inline Manifold parametric ( const Function::Equality eq1,
        	const Function::Equality eq2, const Function::Equality eq3 ) const;

	static Manifold working;

	inline void set_as_working_manifold ( )
	{	Manifold::working = * this;  }

	class Euclid;  class Implicit;  class Parametric;

};  // end of  class Manifold


class Manifold::Core

{	public :

	size_t dim;

	Function metric;  // may be a matrix function

	double (*inner_prod) ( const Cell & P, const std::vector<double> & v,
												 const std::vector<double> & w, const Function & metric );

	inline Core ()
	:	metric (1.), inner_prod { & Manifold::default_inner_prod }
	{ }

	virtual Function build_coord_func
		( const tag::lagrange &, const tag::OfDegree &, size_t d ) = 0;
	
	virtual Function get_coord_func ( ) const = 0;
	
	virtual void set_coords ( const Function co ) = 0;

	virtual void project ( Cell::Positive::Vertex * ) const = 0;

	// P = sA + sB,  s+t == 1
	virtual void interpolate ( Cell::Positive::Vertex * P,
		double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B ) const = 0;

	// P = sA + sB + uC + vD,  s+t+u+v == 1
	virtual void interpolate ( Cell::Positive::Vertex * P,
		double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
		double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D ) const = 0;

	// P = sA + sB + uC + vD + wE + zF,  s+t+u+v+w+z == 1
	virtual void interpolate ( Cell::Positive::Vertex * P,
		double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
		double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D,
		double w, Cell::Positive::Vertex * E, double z, Cell::Positive::Vertex * F ) const = 0;
	
	// P = sum c_k P_k,  sum c_k == 1
	virtual void interpolate ( Cell::Positive::Vertex * P, std::vector < double > & coefs,
														 std::vector < Cell::Positive::Vertex * > & points ) const = 0;

}; // end of class Manifold::Core


inline Function Manifold::coordinates ( ) const
{	assert ( this->core );
	return this->core->get_coord_func();  }


inline Function Manifold::build_coordinate_system
( const tag::lagrange &, const tag::OfDegree &, size_t d )
{	return this->core->build_coord_func ( tag::Lagrange, tag::of_degree, d );  }


inline void Manifold::set_coordinates ( const Function co )
{	assert ( this->core );
	this->core->set_coords ( co );  }


// metric in the manifold (an inner product on the tangent space)
inline double Manifold::inner_prod
( const Cell & P, const std::vector<double> & v, const std::vector<double> & w ) const
{	assert ( v.size() == w.size() );
	assert ( P.dim() == 0 );
	return this->core->inner_prod ( P, v, w, this->core->metric );  }


inline void Manifold::set_metric ( const Function & m )
{	size_t dim = this->coordinates().nb_of_components();
	this->core->metric = m;
	if ( m.nb_of_components() == 1 )
		this->core->inner_prod = Manifold::zoom_inner_prod;
	else
	{	assert ( m.nb_of_components() == dim*dim );
		this->core->inner_prod = Manifold::matrix_inner_prod;  }  }

// P = sA + sB,  s+t == 1
inline void Manifold::interpolate
( const Cell & P, double s, const Cell & A, double t, const Cell & B ) const
{	assert ( P.dim() == 0 );  assert ( A.dim() == 0 );  assert ( B.dim() == 0 );
	assert ( P.is_positive() );  assert ( A.is_positive() );  assert ( B.is_positive() );
	assert ( s >= 0. );  assert ( s <= 1. );
	assert ( t >= 0. );  assert ( t <= 1. );
	// cannot assert the sum is 1. due to round-off errors
	this->core->interpolate ( ( Cell::Positive::Vertex * ) P.core,
		s, ( Cell::Positive::Vertex * ) A.core, t, ( Cell::Positive::Vertex * ) B.core );  }


// P = sA + sB + uC + vD,  s+t+u+v == 1
inline void Manifold::interpolate
( const Cell & P, double s, const Cell & A, double t, const Cell & B,
                  double u, const Cell & C, double v, const Cell & D ) const
{	assert ( P.dim() == 0 );
	assert ( A.dim() == 0 );  assert ( B.dim() == 0 );
	assert ( C.dim() == 0 );  assert ( D.dim() == 0 );
	assert ( P.is_positive() );
	assert ( A.is_positive() );  assert ( B.is_positive() );
	assert ( C.is_positive() );  assert ( D.is_positive() );
	assert ( s >= 0. );  assert ( s <= 1. );
	assert ( t >= 0. );  assert ( t <= 1. );
	assert ( u >= 0. );  assert ( u <= 1. );
	assert ( v >= 0. );  assert ( v <= 1. );
	// cannot assert the sum is 1. due to round-off errors
	this->core->interpolate ( ( Cell::Positive::Vertex * ) P.core,
		s, ( Cell::Positive::Vertex * ) A.core, t, ( Cell::Positive::Vertex * ) B.core,
		u, ( Cell::Positive::Vertex * ) C.core, v, ( Cell::Positive::Vertex * ) D.core );  }


// P = sA + sB + uC + vD + wE + zF,  s+t+u+v+w+z == 1
inline void Manifold::interpolate
( const Cell & P, double s, const Cell & A, double t, const Cell & B,
                  double u, const Cell & C, double v, const Cell & D,
                  double w, const Cell & E, double z, const Cell & F ) const
{	assert ( P.dim() == 0 );
	assert ( A.dim() == 0 );  assert ( B.dim() == 0 );
	assert ( C.dim() == 0 );  assert ( D.dim() == 0 );
	assert ( E.dim() == 0 );  assert ( F.dim() == 0 );
	assert ( P.is_positive() );
	assert ( A.is_positive() );  assert ( B.is_positive() );
	assert ( C.is_positive() );  assert ( D.is_positive() );
	assert ( E.is_positive() );  assert ( F.is_positive() );
	assert ( s >= 0. );  assert ( s <= 1. );
	assert ( t >= 0. );  assert ( t <= 1. );
	assert ( u >= 0. );  assert ( u <= 1. );
	assert ( v >= 0. );  assert ( v <= 1. );
	assert ( w >= 0. );  assert ( w <= 1. );
	assert ( z >= 0. );  assert ( z <= 1. );
	// cannot assert the sum is 1. due to round-off errors
	this->core->interpolate ( ( Cell::Positive::Vertex * ) P.core,
		s, ( Cell::Positive::Vertex * ) A.core, t, ( Cell::Positive::Vertex * ) B.core,
		u, ( Cell::Positive::Vertex * ) C.core, v, ( Cell::Positive::Vertex * ) D.core,
		w, ( Cell::Positive::Vertex * ) E.core, z, ( Cell::Positive::Vertex * ) F.core );  }


// P = sum c_k P_k,  sum c_k == 1
inline void Manifold::interpolate
( const Cell & P, std::vector < double > & coefs, std::vector < Cell > & points ) const
	
{	assert ( points.size() == coefs.size() );
	for ( size_t i = 0; i < points.size(); i++ )  // if debug !
	{	assert ( points[i].is_positive() );
		assert ( points[i].dim() == 0 );
		assert ( coefs[i] >= 0. );   assert ( coefs[i] <= 1. );  }
		// cannot assert the sum is 1. due to round-off errors
	std::vector < Cell::Positive::Vertex * > cores ( coefs.size() );
	for ( size_t i = 0; i < coefs.size(); i++ )
		cores[i] = ( Cell::Positive::Vertex * ) points[i].core;
	this->core->interpolate ( ( Cell::Positive::Vertex * ) P.core, coefs, cores );   }

					 
inline void Manifold::project ( const Cell & cll ) const
{	assert ( cll.is_positive () );
	assert ( cll.dim() == 0 );
	this->core->project ( ( Cell::Positive::Vertex * ) cll.core );  }


inline Manifold Manifold::implicit ( const Function::Equality eq ) const
{	return Manifold ( tag::implicit, *this, eq.lhs - eq.rhs );  }

inline Manifold Manifold::implicit
( const Function::Equality eq1, const Function::Equality eq2 ) const
{	return Manifold ( tag::implicit, *this, eq1.lhs - eq1.rhs, eq2.lhs - eq2.rhs );  }


inline Manifold Manifold::parametric ( const Function::Equality eq ) const
{	return Manifold ( tag::parametric, *this, eq );  }

inline Manifold Manifold::parametric
( const Function::Equality eq1, const Function::Equality eq2 ) const
{	return Manifold ( tag::parametric, *this, eq1, eq2 );  }

inline Manifold Manifold::parametric ( const Function::Equality eq1,
        const Function::Equality eq2, const Function::Equality eq3 ) const
{	return Manifold ( tag::parametric, *this, eq1, eq2, eq3 );  }

//------------------------------------------------------------------------------------------

class Manifold::Euclid : public Manifold::Core

{	public :

	size_t dim;
	
	Field::Core * coord_field { nullptr };
	Function coord_func { tag::non_existent };

	inline Euclid ( size_t d )
	:	Manifold::Core(), dim { d }
	{	assert ( d > 0 );  }

	// P = sA + sB,  s+t == 1     virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
		double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B ) const;
	void pretty_interpolate
		( const Cell & P, double s, const Cell & A, double t, const Cell & B ) const;

	// P = sA + sB + uC + vD,  s+t+u+v == 1     virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
		double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
		double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D  ) const;
	void pretty_interpolate
	(	const Cell & P, double s, const Cell & A, double t, const Cell & B,
		                double u, const Cell & C, double v, const Cell & D  ) const;

	// P = sA + sB + uC + vD + wE + zF,  s+t+u+v+w+z == 1     virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
		double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
		double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D,
		double w, Cell::Positive::Vertex * E, double z, Cell::Positive::Vertex * F ) const;
	void pretty_interpolate ( const Cell & P, double s, const Cell & A,
		double t, const Cell & B, double u, const Cell & C, double v, const Cell & D,
		double w, const Cell & E, double z, const Cell & F ) const;

	// P = sum c_k P_k,  sum c_k == 1     virtual from Manifold::Core
  void interpolate ( Cell::Positive::Vertex * P,
		std::vector < double > & coefs, std::vector < Cell::Positive::Vertex * > & points ) const;
	void pretty_interpolate ( const Cell & P,
		std::vector < double > & coefs, std::vector < Cell > & points ) const;
	
  Function build_coord_func ( const tag::lagrange &, const tag::OfDegree &, size_t d ) ;
  // virtual from Manifold::Core
	
  Function get_coord_func ( ) const;  // virtual from Manifold::Core
	
	void set_coords ( const Function co );  // virtual from Manifold::Core

	void project ( Cell::Positive::Vertex * ) const;
	// virtual from Manifold::Core, here execution forbidden

	class SqDist;
	// a callable object returning the square of the distance between two points
	// used for MetricTree, see paragraph 10.15 in the manual

};  // end of class Manifold::Euclid

//-----------------------------------------------------------------------------------------

class Manifold::Euclid::SqDist

// a callable object returning the square of the distance between two points
// used for MetricTree, see paragraph 10.15 in the manual

{	public :

	inline double operator() ( const Cell & A, const Cell & B )
	{	double res = 0.;
		const size_t nc = Manifold::working.coordinates().nb_of_components();
		for ( size_t i = 0; i < nc; i++ )
		{	Function x = Manifold::working.coordinates()[i];
			double tmp = x(B) - x(A);
			res += tmp*tmp;                                  }
		return res;                                                            }

};  // end of  class Manifold::Euclid::SqDist

// when we first declare a Manifold::Euclid, coord_field and coord_func will be nullptr
// immediately after declaring the manifold (in main), we must declare the coordinates
// to be of type Lagrange degree 1 for instance
//    	Manifold RR2 ( tag::Euclid, 2);
//      Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );

inline Manifold::Manifold
( const tag::euclid &, const tag::OfDimension &, const size_t d )
:	Manifold ( tag::whose_core_is, new Manifold::Euclid ( d ) )
{	assert ( d > 0 );  }

inline Manifold::Manifold
( const tag::euclid &, const tag::OfDimension &, const size_t d,
  const tag::DoNotSetAsWorking & )
:	core { new Manifold::Euclid ( d ) }
{	assert ( this->core );
	assert ( d > 0 );       }


class Manifold::Implicit : public Manifold::Core

// a submanifold of a Manifold::Euclid defined by one or more equations

{	public :

	Manifold surrounding_space;

	inline Implicit ( )	:	Manifold::Core(), surrounding_space ( tag::non_existent ) { }
		
	// the projection will be done by means of the Newton method
	static const short int steps_for_Newton = 10;
	
	// P = sA + sB,  s+t == 1     virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B ) const;

	// P = sA + sB + uC + vD,  s+t+u+v == 1     virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	  double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D  ) const;

	// P = sA + sB + uC + vD + wE + zF,  s+t+u+v+w+z == 1     virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	  double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D,
	  double w, Cell::Positive::Vertex * E, double z, Cell::Positive::Vertex * F ) const;

	// P = sum c_k P_k,  sum c_k == 1     virtual from Manifold::Core
  void interpolate ( Cell::Positive::Vertex * P,
		std::vector < double > & coefs, std::vector < Cell::Positive::Vertex * > & points ) const;
	
	Function build_coord_func ( const tag::lagrange &, const tag::OfDegree &, size_t d );
	//   virtual from Manifold::Core, here execution forbidden
	
	Function get_coord_func ( ) const;  // virtual from Manifold::Core

	void set_coords ( const Function co );  // virtual from Manifold::Core

	// void project ( Cell::Positive::Vertex * ) const  stays pure virtual from Manifold::Core

	class OneEquation; class TwoEquations;
	
};  // end of class Manifold::Implicit


class Manifold::Implicit::OneEquation : public Manifold::Implicit

// a submanifold of a Manifold::Euclid defined by one equation

{	public :

	// Manifold surrounding_space  in Manifold::Implicit
	Function level_function, grad_lev_func;
	
	inline OneEquation ( const Manifold & s, const Function & f );
	
	// void interpolate (different overloaded versions) defined by Manifold::Implicit
	
	void project ( Cell::Positive::Vertex * ) const;
	// virtual from Manifold::Core, through Manifold::Implicit

	// Function build_coord_func ( const tag::lagrange &, const tag::OfDegree &, size_t d )
	//   defined by Manifold::Implicit (execution forbidden)

	// Function get_coord_func ( ) const  defined by Manifold::Implicit

	// void set_coords ( const Function co )  defined by Manifold::Implicit

};  // end of class Manifold::Implicit::OneEquation


class Manifold::Implicit::TwoEquations : public Manifold::Implicit

// a submanifold of a Manifold::Euclid defined by two equations

{	public :

	// Manifold surrounding_space  in Manifold::Implicit
	Function level_function_1, grad_lev_func_1;
	Function level_function_2, grad_lev_func_2;
	
	inline TwoEquations ( const Manifold & s, const Function & f );
	inline TwoEquations ( const Manifold & s, const Function & f1, const Function & f2 );

	// void interpolate (different overloaded versions) defined by Manifold::Implicit
	
	// Function build_coord_func ( const tag::lagrange &, const tag::OfDegree &, size_t d )
	//   defined by Manifold::Implicit (execution forbidden)
	
	// Function get_coord_func ( ) const  defined by Manifold::Implicit

	// void set_coords ( const Function co )  defined by Manifold::Implicit

	void project ( Cell::Positive::Vertex * ) const;
	// virtual from Manifold::Core, through Manifold::Implicit

};  // end of class Manifold::Implicit::TwoEquations


inline Manifold::Manifold ( const tag::Implicit &, const Manifold & m, const Function & f )
:	Manifold ( tag::non_existent )  // temporarily empty manifold
{	// if m is Manifold::Euclid, we want a Manifold::Implicit::OneEquation
	// if m is Manifold::Implicit::OneEquation, we want a Manifold::Implicit::TwoEquations
	Manifold::Euclid * m_euclid = dynamic_cast<Manifold::Euclid*> ( m.core );
	if ( m_euclid )
		this->core = new Manifold::Implicit::OneEquation ( m, f );
	else
	{	Manifold::Implicit::OneEquation * m_one_eq =
			dynamic_cast<Manifold::Implicit::OneEquation*> ( m.core );
		assert ( m_one_eq );
		this->core = new Manifold::Implicit::TwoEquations ( m, f );  }
	Manifold::working = *this;                                                }
	

inline Manifold::Manifold
( const tag::Implicit &, const Manifold & m, const Function & f1, const Function & f2 )
:	Manifold ( tag::non_existent )  // temporarily empty manifold
{	Manifold::Euclid * m_euclid = dynamic_cast<Manifold::Euclid*> ( m.core );
	assert ( m_euclid );
	this->core = new Manifold::Implicit::TwoEquations ( m, f1, f2 );
	Manifold::working = *this;                                                }
	

inline Manifold::Implicit::OneEquation::OneEquation
( const Manifold & m, const Function & f )
: level_function ( f ),
	grad_lev_func ( 0. )  // temporarily zero gradient

{	this->surrounding_space = m;
	Manifold::Euclid * m_euclid = dynamic_cast<Manifold::Euclid*> ( m.core );
	assert ( m_euclid );
	Function coord = m_euclid->coord_func;
	size_t n = coord.nb_of_components();
	Function::Aggregate * grad = new Function::Aggregate ( tag::reserve_size, n );
	for ( size_t i = 0; i < n; i++ ) // grad->components[i] = f.deriv(coord[i]);
		grad->components.emplace_back ( f.deriv(coord[i]) );
	this->grad_lev_func = Function ( tag::whose_core_is, grad );                   }


inline Manifold::Implicit::TwoEquations::TwoEquations
( const Manifold & m, const Function & f )
: level_function_1 ( 0. ),  // temporarily zero function
	grad_lev_func_1 ( 0. ),  // temporarily zero gradient
	level_function_2 ( f ),
	grad_lev_func_2 ( 0. )  // temporarily zero gradient

{	Manifold::Implicit::OneEquation * m_one_eq =
		dynamic_cast < Manifold::Implicit::OneEquation * > ( m.core );
	assert ( m_one_eq );
	this->surrounding_space = m_one_eq->surrounding_space;
	this->level_function_1 = m_one_eq->level_function;
	this->grad_lev_func_1 = m_one_eq->grad_lev_func;
	Function coord = this->surrounding_space.coordinates();
	size_t n = coord.nb_of_components();
	Function::Aggregate * grad = new Function::Aggregate ( tag::reserve_size, n );
	for ( size_t i = 0; i < n; i++ ) // grad->components[i] = f.deriv(coord[i]);
		grad->components.emplace_back ( f.deriv(coord[i]) );
	this->grad_lev_func_2 = Function ( tag::whose_core_is, grad );                   }


inline Manifold::Implicit::TwoEquations::TwoEquations
( const Manifold & m, const Function & f1, const Function & f2 )
: level_function_1 ( f1 ), 
	grad_lev_func_1 ( 0. ),  // temporarily zero gradient
	level_function_2 ( f2 ),
	grad_lev_func_2 ( 0. )  // temporarily zero gradient

{	this->surrounding_space = m;
	Manifold::Euclid * m_euclid = dynamic_cast<Manifold::Euclid*> ( m.core );
	assert ( m_euclid );
	Function coord = m_euclid->coord_func;
	size_t n = coord.nb_of_components();
	Function::Aggregate * grad = new Function::Aggregate ( tag::reserve_size, n );
	for ( size_t i = 0; i < n; i++ ) // grad->components[i] = f1.deriv(coord[i]);
		grad->components.emplace_back ( f1.deriv(coord[i]) );
	this->grad_lev_func_1 = Function ( tag::whose_core_is, grad );
	grad = new Function::Aggregate ( tag::reserve_size, n );
	for ( size_t i = 0; i < n; i++ ) // grad->components[i] = f2.deriv(coord[i]);
		grad->components.emplace_back ( f2.deriv(coord[i]) );
	this->grad_lev_func_2 = Function ( tag::whose_core_is, grad );                 }


class Manifold::Parametric : public Manifold::Core

// a submanifold of a Manifold::Euclid defined by one or more explicit equations

{	public :

	Manifold surrounding_space;

	std::map < Function::Core *, Function::Core * > equations;

	inline Parametric ( )	:	Manifold::Core(), surrounding_space ( tag::non_existent ) { }

	inline Parametric ( const Manifold &, const Function::Equality & );

	inline Parametric ( const Manifold &, const Function::Equality &, const Function::Equality & );

	inline Parametric ( const Manifold &, const Function::Equality &,
                     const Function::Equality &, const Function::Equality & );

	// P = sA + sB,  s+t == 1     virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
		double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B ) const;

	// P = sA + sB + uC + vD,  s+t+u+v == 1     virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
		double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
		double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D  ) const;

	// P = sA + sB + uC + vD + wE + zF,  s+t+u+v+w+z == 1     virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
		double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
		double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D,
		double w, Cell::Positive::Vertex * E, double z, Cell::Positive::Vertex * F ) const;

	// P = sum c_k P_k,  sum c_k == 1     virtual from Manifold::Core
	virtual void interpolate ( Cell::Positive::Vertex * P,
		std::vector < double > & coefs, std::vector < Cell::Positive::Vertex * > & points ) const;
	
	Function build_coord_func ( const tag::lagrange &, const tag::OfDegree &, size_t d );
	//   virtual from Manifold::Core, here execution forbidden
	
	Function get_coord_func ( ) const;  // virtual from Manifold::Core

	void set_coords ( const Function co );  // virtual from Manifold::Core

	void project ( Cell::Positive::Vertex * ) const;
	// virtual from Manifold::Core

};  // end of class Manifold::Parametric



inline Manifold::Manifold ( const tag::Parametric &, const Manifold & m, const Function::Equality & f_eq )
:	Manifold ( tag::non_existent )  // temporary empty manifold
{	Manifold::Euclid * m_euclid = dynamic_cast<Manifold::Euclid*> ( m.core );
	assert ( m_euclid );
	this->core = new Manifold::Parametric ( m, f_eq );
	Manifold::working = *this;                                                 }

inline Manifold::Manifold ( const tag::Parametric &, const Manifold & m,
														const Function::Equality & f_eq_1, const Function::Equality & f_eq_2 )
:	Manifold ( tag::non_existent )  // temporary empty manifold
{	Manifold::Euclid * m_euclid = dynamic_cast<Manifold::Euclid*> ( m.core );
	assert ( m_euclid );
	this->core = new Manifold::Parametric ( m, f_eq_1, f_eq_2 );
	Manifold::working = *this;                                                 }

inline Manifold::Manifold ( const tag::Parametric &, const Manifold & m,
const Function::Equality & f_eq_1, const Function::Equality & f_eq_2, const Function::Equality & f_eq_3 )
:	Manifold ( tag::non_existent )  // temporary empty manifold
{	Manifold::Euclid * m_euclid = dynamic_cast<Manifold::Euclid*> ( m.core );
	assert ( m_euclid );
	this->core = new Manifold::Parametric ( m, f_eq_1, f_eq_2, f_eq_3 );
	Manifold::working = *this;                                                           }


inline Manifold::Parametric::Parametric ( const Manifold & m, const Function::Equality & f_eq )
:	Parametric()
{	this->surrounding_space = m;
	Manifold::Euclid * m_euclid = dynamic_cast<Manifold::Euclid*> ( m.core );
	assert ( m_euclid );
	this->equations [ f_eq.lhs.core ] = f_eq.rhs.core;                         }

inline Manifold::Parametric::Parametric
( const Manifold & m, const Function::Equality & f_eq_1, const Function::Equality & f_eq_2 )
:	Parametric()
{	this->surrounding_space = m;
	Manifold::Euclid * m_euclid = dynamic_cast<Manifold::Euclid*> ( m.core );
	assert ( m_euclid );
	this->equations [ f_eq_1.lhs.core ] = f_eq_1.rhs.core;
	this->equations [ f_eq_2.lhs.core ] = f_eq_2.rhs.core;                         }

inline Manifold::Parametric::Parametric ( const Manifold & m, const Function::Equality & f_eq_1,
  const Function::Equality & f_eq_2, const Function::Equality & f_eq_3 )
:	Parametric()
{	this->surrounding_space = m;
	Manifold::Euclid * m_euclid = dynamic_cast<Manifold::Euclid*> ( m.core );
	assert ( m_euclid );
	this->equations [ f_eq_1.lhs.core ] = f_eq_1.rhs.core;
	this->equations [ f_eq_2.lhs.core ] = f_eq_2.rhs.core;
	this->equations [ f_eq_3.lhs.core ] = f_eq_3.rhs.core;                         }


inline void Cell::project () const
{	this->project ( tag::onto, Manifold::working );  }
	
inline void Cell::project ( const tag::Onto &, const Manifold m ) const
{	assert ( m.core );
	assert ( this->core->get_dim() == 0 );
	assert ( this->core->is_positive() );
	Cell::Positive::Vertex * cll = ( Cell::Positive::Vertex * ) this->core;
	m.core->project ( cll );                                                }
		

void inline Mesh::baricenter ( const Cell & ver, const Cell & seg_ini )

// 'ver' is a vertex in 'this' mesh, 'seg' is a segment pointing towards 'ver'

{	assert ( ver.dim() == 0 );
	assert ( seg_ini.dim() == 1 );
	std::vector < Cell > neighbours;  // vertices
	Cell seg = seg_ini;
	size_t n = 0;
	while ( true )
	{	assert ( seg.tip() == ver );
		Cell tri = this->cell_in_front_of ( seg, tag::may_not_exist );
		if ( not tri.exists() ) return;
		// 'ver' is on the boundary of mesh so we do nothing
		n++;  neighbours.push_back ( seg.base().reverse() );
		seg = tri.boundary().cell_behind(ver);
		if ( seg == seg_ini ) break;                                    }
	std::vector < double > coefs ( n, 1./n );
	Manifold::working.interpolate ( ver, coefs, neighbours );             }


}  // end of  namespace maniFEM

#endif
// ifndef MANIFEM_MANIFOLD_H

