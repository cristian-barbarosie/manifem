
// manifold.h 2022.06.25

//   This file is part of maniFEM, a C++ library for meshes and finite elements on manifolds.

//   Copyright 2019 -- 2022 Cristian Barbarosie cristian.barbarosie@gmail.com

//   http://manifem.rd.ciencias.ulisboa.pt/
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

namespace tag
{	struct euclid { };  static const euclid Euclid;
	struct lagrange { };  static const lagrange Lagrange;
	struct Implicit { };  static const Implicit implicit;
	struct Intersect { };  static const Intersect intersect;
	struct Parametric { };  static const Parametric parametric;
	struct DoNotSetAsWorking { };  static const DoNotSetAsWorking do_not_set_as_working;
	struct StartWithInconsistentMesh { };
	static const StartWithInconsistentMesh start_with_inconsistent_mesh;                }


//------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------//

// class below will allow us to easily manipulate the metric of a manifold

class tag::Util::Metric

// abstract class, specialized in tag::Util::Metric::Trivial, Isotropic::***, Anisotropic::***

{	public :

	class Trivial;  class Isotropic;  class Anisotropic;

	inline   Metric() { }
	virtual ~Metric() { }

	// metamorphosis : produce a modified copy of self (a new, scaled, metric)
	virtual tag::Util::Metric * scale ( const Function & f ) = 0;
	virtual tag::Util::Metric * scale ( const double f ) = 0;

	// inner product
	virtual double inner_prod
	( const Cell & P, const std::vector < double > & v, const std::vector < double > & w ) = 0;	
	virtual double inner_prod
	( const Cell & A, const Cell & B,
	  const std::vector < double > & v, const std::vector < double > & w ) = 0;	

	// compute (square of) distance between two points
	inline double sq_dist ( const Cell & A, const Cell & B );
	inline double sq_dist ( const Cell & A, const Cell & B, const std::vector < double > & w );
	
};  // end of  class tag::Util::Metric

	
inline double tag::Util::Metric::sq_dist
( const Cell & A, const Cell & B, const std::vector < double > & w )
{	return this->inner_prod ( A, B, w, w );  }
		
	
class tag::Util::Metric::Trivial : public tag::Util::Metric

{	public :

	inline   Trivial() : tag::Util::Metric()  { }
	virtual ~Trivial() { }

	// two 'scale' methods are virtual from tag::Util::Metric
	tag::Util::Metric * scale ( const Function & f );
	tag::Util::Metric * scale ( const double f );

	// two 'inner_prod' methods are virtual from tag::Util::Metric
	double inner_prod
	( const Cell & P, const std::vector < double > & v, const std::vector < double > & w );	
	double inner_prod
	( const Cell & A, const Cell & B,
	  const std::vector < double > & v, const std::vector < double > & w );	

	// two inlined 'sq_dist' methods defined by tag::Util::Metric

};  // end of  class tag::Util::Metric::Trivial


class tag::Util::Metric::Isotropic : public tag::Util::Metric

// abstract class, specialized in tag::Util::Metric::Isotropic::Constant, Variable

{	public :

	inline   Isotropic() : tag::Util::Metric()  { }
	virtual ~Isotropic() { }

	// two 'scale' methods stay pure virtual from tag::Util::Metric
	
	// two inlined 'sq_dist' methods defined by tag::Util::Metric

	class Constant; class Variable;

};  // end of  class tag::Util::Metric::Isotropic

	
class tag::Util::Metric::Isotropic::Constant : public tag::Util::Metric::Isotropic

{	public :

	double zoom;
	
	inline   Constant ( const double z ) : tag::Util::Metric::Isotropic(), zoom {z}  { }
	virtual ~Constant() { }

	// two 'scale' methods are virtual from tag::Util::Metric, defined in frontal.cpp
	tag::Util::Metric * scale ( const Function & f );
	tag::Util::Metric * scale ( const double f );

	// inner product, virtual from tag::Util::Metric
	virtual double inner_prod
	( const Cell & P, const std::vector < double > & v, const std::vector < double > & w );	
	virtual double inner_prod
	( const Cell & A, const Cell & B,
	  const std::vector < double > & v, const std::vector < double > & w );	

	// two inlined 'sq_dist' methods defined by tag::Util::Metric

};  // end of  class tag::Util::Metric:Isotropic::Constant;


class tag::Util::Metric::Isotropic::Variable : public tag::Util::Metric::Isotropic

{	public :

	Function zoom;
	
	inline   Variable ( const Function & z ) : tag::Util::Metric::Isotropic(), zoom {z}  { }
	virtual ~Variable() { }

	// two 'scale' methods are virtual from tag::Util::Metric, defined in frontal.cpp
	tag::Util::Metric * scale ( const Function & f );
	tag::Util::Metric * scale ( const double f );

	// inner product, virtual from tag::Util::Metric
	virtual double inner_prod
	( const Cell & P, const std::vector < double > & v, const std::vector < double > & w );	
	virtual double inner_prod
	( const Cell & A, const Cell & B,
	  const std::vector < double > & v, const std::vector < double > & w );	

	// two inlined 'sq_dist' methods defined by tag::Util::Metric

};  // end of  class tag::Util::Metric::Isotropic::Variable


class tag::Util::Metric::Anisotropic : public tag::Util::Metric

// abstract class, specialized in tag::Util::Metric::Anisotropic::***

{	public :

	inline   Anisotropic() : tag::Util::Metric() { }
	virtual ~Anisotropic() { }

	// two 'scale' methods stay pure virtual from tag::Util::Metric

	// two 'inner_product' methods stay pure virtual from tag::Util::Metric

	// two inlined 'sq_dist' methods defined by tag::Util::Metric

	class Matrix;  class Rayleigh;

};  // end of  class tag::Util::Metric::Anisotropic

	
class tag::Util::Metric::Anisotropic::Matrix : public tag::Util::Metric::Anisotropic

// abstract class, specialized in tag::Util::Metric::Anisotropic::Matrix::Constant, Variable

{	public :

	inline   Matrix() : tag::Util::Metric::Anisotropic()  { }
	virtual ~Matrix() { }

	// two 'scale' methods stay pure virtual from tag::Util::Metric

	// two 'inner_product' methods stay pure virtual from tag::Util::Metric

	// two inlined 'sq_dist' methods defined by tag::Util::Metric

	class Constant; class Variable;

};  // end of  class tag::Util::Metric::Anisotropic::Matrix


class tag::Util::Metric::Anisotropic::Matrix::Constant
:	public tag::Util::Metric::Anisotropic::Matrix

{	public :

	std::vector < double > matrix;  // n^2 components, definite positive
	
	inline Constant ( const double a11, const double a12, const double a22 )
	:	tag::Util::Metric::Anisotropic::Matrix(), matrix { a11, a12, a22 }
	{ }
	
	inline Constant ( const double a11, const double a12, const double a13,
	                  const double a22, const double a23, const double a33 )
	:	tag::Util::Metric::Anisotropic::Matrix(), matrix { a11, a12, a13, a22, a23, a33 }
	{ }

	virtual ~Constant() { };

	// two 'scale' methods are virtual from tag::Util::Metric, defined in frontal.cpp
	tag::Util::Metric * scale ( const Function & f );
	tag::Util::Metric * scale ( const double f );

	// inner product, virtual from tag::Util::Metric
	virtual double inner_prod
	( const Cell & P, const std::vector < double > & v, const std::vector < double > & w );	
	virtual double inner_prod
	( const Cell & A, const Cell & B,
	  const std::vector < double > & v, const std::vector < double > & w );	

	// two inlined 'sq_dist' methods defined by tag::Util::Metric

};  // end of  class tag::Util::Metric::Anisotropic::Matrix::Constant


class tag::Util::Metric::Anisotropic::Matrix::Variable
:	public tag::Util::Metric::Anisotropic::Matrix

{	public :

	Function matrix;  // n^2 components, definite positive
	
	inline   Variable() : tag::Util::Metric::Anisotropic::Matrix() { };
	virtual ~Variable() { };

	// two 'scale' methods are virtual from tag::Util::Metric, defined in frontal.cpp
	tag::Util::Metric * scale ( const Function & f );
	tag::Util::Metric * scale ( const double f );

	// inner product, virtual from tag::Util::Metric
	virtual double inner_prod
	( const Cell & P, const std::vector < double > & v, const std::vector < double > & w );	
	virtual double inner_prod
	( const Cell & A, const Cell & B,
	  const std::vector < double > & v, const std::vector < double > & w );	

	// two inlined 'sq_dist' methods defined by tag::Util::Metric

};  // end of  class tag::Util::Metric::Anisotropic::Matrix::Variable


class tag::Util::Metric::Anisotropic::Rayleigh : public tag::Util::Metric::Anisotropic

// abstract class, specialized in tag::Util::Metric::Anisotropic::Rayleigh::Constant, Variable

{	public :

	Function rayl;  // scalar, positive
	Function dev;   // n^2 components, semi-definite positive
	// matrix = rayl * Id + dev
	
	inline   Rayleigh() { }
	virtual ~Rayleigh() { }

	// two 'scale' methods stay pure virtual from tag::Util::Metric

	// two inlined 'sq_dist' methods defined by tag::Util::Metric

	class Constant; class Variable;

};  // end of  class tag::Util::Metric::Anisotropic::Rayleigh


//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------

// there will be just a few Manifold objects in a program,
// and they will usually be destroyed only at the end of the program,
// so we do not use here the mechanism of inheriting from tag::Util::Wrapper and tag::Util::Core
	
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------


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
	                  const tag::DoNotSetAsWorking &                                  );

	inline Manifold ( const tag::Implicit &, const Manifold &, const Function & );
	inline Manifold ( const tag::Implicit &,
	                  const Manifold &, const Function &, const Function & );
	inline Manifold ( const tag::Implicit &,
	                  const Manifold &, const Function &, const Function &, const Function & );
	
	inline Manifold ( const tag::Intersect &, const Manifold &, const Manifold & );
	inline Manifold
	( const tag::Intersect &, const Manifold &, const Manifold &, const Manifold & );
	
	inline Manifold ( const tag::Parametric &, const Manifold &, const Function::Equality & );
	inline Manifold ( const tag::Parametric &,
	       const Manifold &, const Function::Equality &, const Function::Equality & );
	inline Manifold ( const tag::Parametric &, const Manifold &,
	       const Function::Equality &, const Function::Equality &, const Function::Equality & );

	inline ~Manifold ( )
	{} //	std::cout << "destructor Manifold" << std::endl;
	
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

	inline void build_intersection_of_two_manif ( const Manifold & m1, const Manifold & m2 );
		
	// measure of the entire manifold (length for 1D, area for 2D, volume for 3D)
	// produces run-time error for Euclidian manifolds
	// and for other manifolds whose measure is infinite (e.g. cylinder)
	// and for manifolds whose measure is too difficult to compute (e.g. implicit manifolds)
	// significant for torus
 	inline double measure ( ) const;
	
	// metric in the manifold (an inner product on the tangent space)
	inline double inner_prod ( const Cell & P, const std::vector<double> & v,
                                             const std::vector<double> & w ) const;

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
	
	typedef tag::Util::Action Action;  // aka  class Function::Action
	// we define it in function.h because we need it for Function::MultiValued
	// but we prefer the user to see it as an attribute of class Manifold

	// P = sA + sB,  s+t == 1
	inline void interpolate
	( const Cell & P, double s, const Cell & A, double t, const Cell & B ) const;

	// P = sA + sB,  s+t == 1
	inline void interpolate
	( const Cell & P, double s, const Cell & A, double t, const Cell & B,
	  const tag::Winding &, const Manifold::Action & exp            ) const;

	// P = sA + sB + uC + vD,  s+t+u+v == 1
	inline void interpolate ( const Cell & P, double s, const Cell & A,
	  double t, const Cell & B, double u, const Cell & C, double v, const Cell & D ) const;

	// P = sA + sB + uC,  s+t+u == 1
	inline void interpolate ( const Cell & P, double s, const Cell & A,
	  double t, const Cell & B, double u, const Cell & C               ) const;

	// P = sA + sB + uC + vD,  s+t+u+v == 1
	inline void interpolate ( const Cell & P, double s, const Cell & A,
	  double t, const Cell & B, const tag::Winding &, const Manifold::Action &,
	  double u, const Cell & C, const tag::Winding &, const Manifold::Action &,
	  double v, const Cell & D, const tag::Winding &, const Manifold::Action & ) const;

	// P = sA + sB + uC,  s+t+u == 1
	inline void interpolate ( const Cell & P, double s, const Cell & A,
	  double t, const Cell & B, const tag::Winding &, const Manifold::Action &,
	  double u, const Cell & C, const tag::Winding &, const Manifold::Action & ) const;

	// P = sA + sB + uC + vD + wE + zF,  s+t+u+v+w+z == 1
	inline void interpolate ( const Cell & P, double s, const Cell & A,
	  double t, const Cell & B, double u, const Cell & C, double v, const Cell & D,
	  double w, const Cell & E, double z, const Cell & F                           ) const;

	// P = sA + sB + uC + vD + wE + zF,  s+t+u+v+w+z == 1
	inline void interpolate ( const Cell & P, double s, const Cell & A,
	  double t, const Cell & B, const tag::Winding &, const Manifold::Action &,
	  double u, const Cell & C, const tag::Winding &, const Manifold::Action &,
	  double v, const Cell & D, const tag::Winding &, const Manifold::Action &,
	  double w, const Cell & E, const tag::Winding &, const Manifold::Action &,
	  double z, const Cell & F, const tag::Winding &, const Manifold::Action & ) const;

	// P = sum c_k P_k,  sum c_k == 1
	inline void interpolate
	( const Cell & P, const std::vector < double > & coefs,
	  const std::vector < Cell > & points                  ) const;

	// P = sum c_k P_k,  sum c_k == 1
	inline void interpolate
	( const Cell & P, const std::vector < double > & coefs,
	  const std::vector < Cell > & points, const tag::Winding &,
	  const std::vector < Manifold::Action >  ) const;

	inline void project ( const Cell & ) const;
	
	inline Manifold implicit ( const Function::Equality eq ) const;
	inline Manifold implicit ( const Function::Equality eq1,
	                           const Function::Equality eq2 ) const;
	inline Manifold parametric ( const Function::Equality eq ) const;
	inline Manifold parametric ( const Function::Equality eq1,
	                             const Function::Equality eq2 ) const;
	inline Manifold parametric ( const Function::Equality eq1,
        	const Function::Equality eq2, const Function::Equality eq3 ) const;

	inline Manifold quotient ( const Manifold::Action & a1 );
	inline Manifold quotient ( const Manifold::Action & a1, const Manifold::Action & a2 );
	inline Manifold quotient
	( const Manifold::Action & a1, const Manifold::Action & a2, const Manifold::Action & a3 );

	static Manifold working;

	inline void set_as_working_manifold ( )
	{	Manifold::working = * this;  }

	class Euclid;  class Implicit;  class Parametric;  class Quotient;

	struct Type  // used in frontal.cpp
	{	class Euclidian;  class Quotient;  };

};  // end of  class Manifold

//------------------------------------------------------------------------------------------------------//


inline double tag::Util::Metric::sq_dist ( const Cell & A, const Cell & B )

{	Manifold & space = Manifold::working;
	assert ( space .exists() );  // we use the current manifold
	const size_t n = space .coordinates() .nb_of_components();
	std::vector < double > w ( n );
	for ( size_t i = 0; i < n; i++ )
	{	const Function & x = space .coordinates() [i];
		w[i] = x(B) - x(A);                            }
	return this->sq_dist ( A, B, w );                         }

//------------------------------------------------------------------------------------------------------//


inline Cell::Cell ( const tag::Vertex &, const tag::OfCoordinates &, const std::vector < double > & v,
                    const tag::IsPositive & ispos                                                     )
// by default, ispos = tag::is_positive, so may be called with only three arguments
:	Cell ( tag::whose_core_is, new Cell::Positive::Vertex ( tag::one_dummy_wrapper ),
         tag::freshly_created                                                       )
{	Manifold & space = Manifold::working;
	assert ( space .exists() );  // we use the current manifold
	Function coords = space .coordinates();
	assert ( v .size() == coords .nb_of_components() );
	coords ( *this ) = v;                              }


inline Cell::Cell ( const tag::Vertex &, const tag::OfCoordinates &, const std::vector < double > & v,
                    const tag::Project &, const tag::IsPositive & ispos                               )
// by default, ispos = tag::is_positive, so may be called with only three arguments
:	Cell ( tag::whose_core_is, new Cell::Positive::Vertex ( tag::one_dummy_wrapper ),
         tag::freshly_created                                                       )
{	Manifold & space = Manifold::working;
	assert ( space .exists() );  // we use the current manifold
	Function coords = space .coordinates();
	assert ( v .size() == coords .nb_of_components() );
	coords ( *this ) = v;
	space .project ( *this );                           }

//------------------------------------------------------------------------------------------------------//


class Manifold::Core

{	public :

	size_t dim;

	tag::Util::Metric * metric;

	inline Core ( )
	:	metric ( new tag::Util::Metric::Trivial )
	{ }

	virtual ~Core ( ) { };

	virtual Function build_coord_func
	( const tag::lagrange &, const tag::OfDegree &, size_t d ) = 0;
	
	virtual Function get_coord_func ( ) const = 0;
	
	virtual void set_coords ( const Function co ) = 0;

	virtual void project ( Cell::Positive::Vertex * ) const = 0;

	// measure of the entire manifold (length for 1D, area for 2D, volume for 3D)
	// produces run-time error for Euclidian manifolds
	// and for other manifolds whose measure is infinite (e.g. cylinder)
	// and for manifolds whose measure is too difficult to compute (e.g. implicit manifolds)
	// significant for quotient manifolds (torus)
	
	virtual double measure ( ) const = 0;
	
	// P = sA + sB,  s+t == 1
	virtual void interpolate ( Cell::Positive::Vertex * P,
	 double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B ) const = 0;

	// P = sA + sB,  s+t == 1
	virtual void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	  const tag::Winding &, const Manifold::Action & exp                         ) const = 0;

	// P = sA + sB + uC,  s+t+u == 1
	virtual void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	  double u, Cell::Positive::Vertex * C                                       ) const = 0;

	// P = sA + sB + uC,  s+t+u == 1
	virtual void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A,
	  double t, Cell::Positive::Vertex * B,
	  const tag::Winding &, const Manifold::Action &,
	  double u, Cell::Positive::Vertex * C,
	  const tag::Winding &, const Manifold::Action & ) const = 0;

	// P = sA + sB + uC + vD,  s+t+u+v == 1
	virtual void interpolate ( Cell::Positive::Vertex * P,
		double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
		double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D ) const = 0;

	// P = sA + sB + uC + vD,  s+t+u+v == 1
	virtual void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A,
	  double t, Cell::Positive::Vertex * B,
	  const tag::Winding &, const Manifold::Action &,
	  double u, Cell::Positive::Vertex * C,
	  const tag::Winding &, const Manifold::Action &,
	  double v, Cell::Positive::Vertex * D,
	  const tag::Winding &, const Manifold::Action & ) const = 0;

	// P = sA + sB + uC + vD + wE + zF,  s+t+u+v+w+z == 1
	virtual void interpolate ( Cell::Positive::Vertex * P,
		double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
		double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D,
		double w, Cell::Positive::Vertex * E, double z, Cell::Positive::Vertex * F ) const = 0;
	
	// P = sA + sB + uC + vD + wE + zF,  s+t+u+v+w+z == 1
	virtual void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A,
	  double t, Cell::Positive::Vertex * B,
	  const tag::Winding &, const Manifold::Action &,
	  double u, Cell::Positive::Vertex * C,
	  const tag::Winding &, const Manifold::Action &,
	  double v, Cell::Positive::Vertex * D,
	  const tag::Winding &, const Manifold::Action &,
	  double w, Cell::Positive::Vertex * E,
	  const tag::Winding &, const Manifold::Action &,
	  double z, Cell::Positive::Vertex * F,
	  const tag::Winding &, const Manifold::Action & ) const = 0;

	// P = sum c_k P_k,  sum c_k == 1
	virtual void interpolate ( Cell::Positive::Vertex * P, const std::vector < double > & coefs,
														 const std::vector < Cell::Positive::Vertex * > & points ) const = 0;

	virtual void frontal_method  // defined in frontal.cpp, overridden by Manifold::Quotient
	( Mesh & msh, const tag::StartWithInconsistentMesh &,
	  const tag::StartAt &, const Cell & start,
	  const tag::Towards &, std::vector<double> tangent,
	  const tag::StopAt &,  const Cell & stop            );
	
}; // end of class Manifold::Core

//------------------------------------------------------------------------------------------------------//


inline Function Manifold::coordinates ( ) const
{	assert ( this->core );
	return this->core->get_coord_func();  }


inline Function Manifold::build_coordinate_system
( const tag::lagrange &, const tag::OfDegree &, size_t d )
{	return this->core->build_coord_func ( tag::Lagrange, tag::of_degree, d );  }


inline void Manifold::set_coordinates ( const Function co )
{	assert ( this->core );
	this->core->set_coords ( co );  }


inline double Manifold::measure ( ) const
{	return this->core->measure();  }


// P = sA + sB,  s+t == 1
inline void Manifold::interpolate
( const Cell & P, double s, const Cell & A, double t, const Cell & B ) const

{	assert ( P .dim() == 0 );  assert ( A .dim() == 0 );  assert ( B .dim() == 0 );
	assert ( P .is_positive() );  assert ( A .is_positive() );  assert ( B .is_positive() );
	assert ( s >= 0. );  assert ( s <= 1. );
	assert ( t >= 0. );  assert ( t <= 1. );

	double sum = s + t;
	// cannot assert the sum is 1. due to round-off errors
	assert ( sum > 0.99 );  assert ( sum < 1.01 );

	Cell::Positive::Vertex * P_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( P.core );
	Cell::Positive::Vertex * A_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	Cell::Positive::Vertex * B_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	this->core->interpolate ( P_c, s, A_c, t, B_c );                                       }


// P = sA + sB,  s+t == 1
inline void Manifold::interpolate
( const Cell & P, double s, const Cell & A, double t, const Cell & B,
  const tag::Winding &, const Manifold::Action & exp_AB              ) const

{	assert ( P .dim() == 0 );  assert ( A .dim() == 0 );  assert ( B .dim() == 0 );
	assert ( P .is_positive() );  assert ( A .is_positive() );  assert ( B .is_positive() );
	assert ( s >= 0. );  assert ( s <= 1. );
	assert ( t >= 0. );  assert ( t <= 1. );

	double sum = s + t;
	// cannot assert the sum is 1. due to round-off errors
	assert ( sum > 0.99 );  assert ( sum < 1.01 );

	Cell::Positive::Vertex * P_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( P.core );
	Cell::Positive::Vertex * A_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	Cell::Positive::Vertex * B_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	this->core->interpolate ( P_c, s, A_c, t, B_c, tag::winding, exp_AB );                    }


// P = sA + sB + uC,  s+t+u == 1
inline void Manifold::interpolate
( const Cell & P, double s, const Cell & A, double t, const Cell & B,
                  double u, const Cell & C                           ) const

{	assert ( P .dim() == 0 );
	assert ( A .dim() == 0 );  assert ( B .dim() == 0 );
	assert ( C .dim() == 0 );
	assert ( P .is_positive() );
	assert ( A .is_positive() );  assert ( B .is_positive() );
	assert ( C .is_positive() );
	assert ( s >= 0. );  assert ( s <= 1. );
	assert ( t >= 0. );  assert ( t <= 1. );
	assert ( u >= 0. );  assert ( u <= 1. );

	double sum = s + t + u;
	// cannot assert the sum is 1. due to round-off errors
	assert ( sum > 0.99 );  assert ( sum < 1.01 );

	Cell::Positive::Vertex * P_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( P.core );
	Cell::Positive::Vertex * A_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	Cell::Positive::Vertex * B_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	Cell::Positive::Vertex * C_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( C.core );
	this->core->interpolate ( P_c, s, A_c, t, B_c, u, C_c );  }


// P = sA + sB + uC,  s+t+u == 1
inline void Manifold::interpolate
( const Cell & P, double s, const Cell & A,
  double t, const Cell & B, const tag::Winding &, const Manifold::Action & exp_AB,
  double u, const Cell & C, const tag::Winding &, const Manifold::Action & exp_AC ) const

{	assert ( P .dim() == 0 );
	assert ( A .dim() == 0 );  assert ( B .dim() == 0 );
	assert ( C .dim() == 0 );
	assert ( P .is_positive() );
	assert ( A .is_positive() );  assert ( B .is_positive() );
	assert ( C .is_positive() );
	assert ( s >= 0. );  assert ( s <= 1. );
	assert ( t >= 0. );  assert ( t <= 1. );
	assert ( u >= 0. );  assert ( u <= 1. );

	double sum = s + t + u;
	// cannot assert the sum is 1. due to round-off errors
	assert ( sum > 0.99 );  assert ( sum < 1.01 );

	Cell::Positive::Vertex * P_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( P.core );
	Cell::Positive::Vertex * A_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	Cell::Positive::Vertex * B_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	Cell::Positive::Vertex * C_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( C.core );
	this->core->interpolate ( P_c, s, A_c, t, B_c, tag::winding, exp_AB,
              u, C_c, tag::winding, exp_AC                            );  }


// P = sA + sB + uC + vD,  s+t+u+v == 1
inline void Manifold::interpolate
( const Cell & P, double s, const Cell & A, double t, const Cell & B,
                  double u, const Cell & C, double v, const Cell & D ) const

{	assert ( P .dim() == 0 );
	assert ( A .dim() == 0 );  assert ( B .dim() == 0 );
	assert ( C .dim() == 0 );  assert ( D .dim() == 0 );
	assert ( P .is_positive() );
	assert ( A .is_positive() );  assert ( B .is_positive() );
	assert ( C .is_positive() );  assert ( D .is_positive() );
	assert ( s >= 0. );  assert ( s <= 1. );
	assert ( t >= 0. );  assert ( t <= 1. );
	assert ( u >= 0. );  assert ( u <= 1. );
	assert ( v >= 0. );  assert ( v <= 1. );

	double sum = s + t + u + v;
	// cannot assert the sum is 1. due to round-off errors
	assert ( sum > 0.99 );  assert ( sum < 1.01 );

	Cell::Positive::Vertex * P_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( P.core );
	Cell::Positive::Vertex * A_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	Cell::Positive::Vertex * B_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	Cell::Positive::Vertex * C_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( C.core );
	Cell::Positive::Vertex * D_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( D.core );
	this->core->interpolate ( P_c, s, A_c, t, B_c, u, C_c, v, D_c );  }


// P = sA + sB + uC + vD,  s+t+u+v == 1
inline void Manifold::interpolate
( const Cell & P, double s, const Cell & A,
  double t, const Cell & B, const tag::Winding &, const Manifold::Action & exp_AB,
  double u, const Cell & C, const tag::Winding &, const Manifold::Action & exp_AC,
  double v, const Cell & D, const tag::Winding &, const Manifold::Action & exp_AD ) const

{	assert ( P .dim() == 0 );
	assert ( A .dim() == 0 );  assert ( B .dim() == 0 );
	assert ( C .dim() == 0 );  assert ( D .dim() == 0 );
	assert ( P .is_positive() );
	assert ( A .is_positive() );  assert ( B .is_positive() );
	assert ( C .is_positive() );  assert ( D .is_positive() );
	assert ( s >= 0. );  assert ( s <= 1. );
	assert ( t >= 0. );  assert ( t <= 1. );
	assert ( u >= 0. );  assert ( u <= 1. );
	assert ( v >= 0. );  assert ( v <= 1. );

	double sum = s + t + u + v;
	// cannot assert the sum is 1. due to round-off errors
	assert ( sum > 0.99 );  assert ( sum < 1.01 );

	Cell::Positive::Vertex * P_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( P.core );
	Cell::Positive::Vertex * A_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	Cell::Positive::Vertex * B_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	Cell::Positive::Vertex * C_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( C.core );
	Cell::Positive::Vertex * D_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( D.core );
	this->core->interpolate ( P_c, s, A_c, t, B_c, tag::winding, exp_AB,
              u, C_c, tag::winding, exp_AC, v, D_c, tag::winding, exp_AD );  }


// P = sA + sB + uC + vD + wE + zF,  s+t+u+v+w+z == 1
inline void Manifold::interpolate
( const Cell & P, double s, const Cell & A, double t, const Cell & B,
                  double u, const Cell & C, double v, const Cell & D,
                  double w, const Cell & E, double z, const Cell & F ) const

{	assert ( P .dim() == 0 );
	assert ( A .dim() == 0 );  assert ( B .dim() == 0 );
	assert ( C .dim() == 0 );  assert ( D .dim() == 0 );
	assert ( E .dim() == 0 );  assert ( F .dim() == 0 );
	assert ( P .is_positive() );
	assert ( A .is_positive() );  assert ( B .is_positive() );
	assert ( C .is_positive() );  assert ( D .is_positive() );
	assert ( E .is_positive() );  assert ( F .is_positive() );
	assert ( s >= 0. );  assert ( s <= 1. );
	assert ( t >= 0. );  assert ( t <= 1. );
	assert ( u >= 0. );  assert ( u <= 1. );
	assert ( v >= 0. );  assert ( v <= 1. );
	assert ( w >= 0. );  assert ( w <= 1. );
	assert ( z >= 0. );  assert ( z <= 1. );

	double sum = s + t + u + v + w + z;
	// cannot assert the sum is 1. due to round-off errors
	assert ( sum > 0.99 );  assert ( sum < 1.01 );

	Cell::Positive::Vertex * P_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( P.core );
	Cell::Positive::Vertex * A_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	Cell::Positive::Vertex * B_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	Cell::Positive::Vertex * C_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( C.core );
	Cell::Positive::Vertex * D_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( D.core );
	Cell::Positive::Vertex * E_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( E.core );
	Cell::Positive::Vertex * F_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( F.core );
	this->core->interpolate ( P_c, s, A_c, t, B_c, u, C_c, v, D_c, w, E_c, z, F_c );  }


// P = sA + sB + uC + vD + wE + zF,  s+t+u+v+w+z == 1
inline void Manifold::interpolate
( const Cell & P, double s, const Cell & A,
  double t, const Cell & B, const tag::Winding &, const Manifold::Action & exp_AB,
  double u, const Cell & C, const tag::Winding &, const Manifold::Action & exp_AC,
  double v, const Cell & D, const tag::Winding &, const Manifold::Action & exp_AD,
  double w, const Cell & E, const tag::Winding &, const Manifold::Action & exp_AE,
  double z, const Cell & F, const tag::Winding &, const Manifold::Action & exp_AF ) const
{	assert ( P .dim() == 0 );
	assert ( A .dim() == 0 );  assert ( B .dim() == 0 );
	assert ( C .dim() == 0 );  assert ( D .dim() == 0 );
	assert ( E .dim() == 0 );  assert ( F .dim() == 0 );
	assert ( P .is_positive() );
	assert ( A .is_positive() );  assert ( B .is_positive() );
	assert ( C .is_positive() );  assert ( D .is_positive() );
	assert ( E .is_positive() );  assert ( F .is_positive() );
	assert ( s >= 0. );  assert ( s <= 1. );
	assert ( t >= 0. );  assert ( t <= 1. );
	assert ( u >= 0. );  assert ( u <= 1. );
	assert ( v >= 0. );  assert ( v <= 1. );
	assert ( w >= 0. );  assert ( w <= 1. );
	assert ( z >= 0. );  assert ( z <= 1. );

	double sum = s + t + u + v + w + z;
	// cannot assert the sum is 1. due to round-off errors
	assert ( sum > 0.99 );  assert ( sum < 1.01 );

	Cell::Positive::Vertex * P_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( P.core );
	Cell::Positive::Vertex * A_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	Cell::Positive::Vertex * B_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	Cell::Positive::Vertex * C_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( C.core );
	Cell::Positive::Vertex * D_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( D.core );
	Cell::Positive::Vertex * E_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( E.core );
	Cell::Positive::Vertex * F_c = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( F.core );
	this->core->interpolate ( P_c, s, A_c, t, B_c, tag::winding, exp_AB,
              u, C_c, tag::winding, exp_AC, v, D_c, tag::winding, exp_AD,
              w, E_c, tag::winding, exp_AE, z, F_c, tag::winding, exp_AF );  }


// P = sum c_k P_k,  sum c_k == 1
inline void Manifold::interpolate
( const Cell & P, const std::vector < double > & coefs, const std::vector < Cell > & points ) const
	
{	assert ( points.size() == coefs .size() );
	for ( size_t i = 0; i < points .size(); i++ )  // if debug !!
	{	assert ( points [i] .is_positive() );
		assert ( points [i] .dim() == 0 );
		assert ( coefs [i] >= 0. );   assert ( coefs [i] <= 1. );  }
		// cannot assert the sum is 1. due to round-off errors
	std::vector < Cell::Positive::Vertex * > cores ( coefs .size() );
	for ( size_t i = 0; i < coefs.size(); i++ )
		cores [i] = ( Cell::Positive::Vertex * ) points [i] .core;
	this->core->interpolate ( ( Cell::Positive::Vertex * ) P .core, coefs, cores );   }

					 
inline void Manifold::project ( const Cell & cll ) const
{	assert ( cll .is_positive () );
	assert ( cll .dim() == 0 );
	this->core->project ( tag::Util::assert_cast
                        < Cell::Core*, Cell::Positive::Vertex* > ( cll .core ) );  }


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

inline Manifold Manifold::parametric
( const Function::Equality eq1, const Function::Equality eq2, const Function::Equality eq3 ) const
{	return Manifold ( tag::parametric, *this, eq1, eq2, eq3 );  }

//------------------------------------------------------------------------------------------------------//


class Manifold::Euclid : public Manifold::Core

{	public :

	size_t dim;
	
	Field::Core * coord_field { nullptr };
	Function coord_func { tag::non_existent };

	inline Euclid ( size_t d )
	:	Manifold::Core(), dim { d }
	{	assert ( d > 0 );  }

	// P = sA + sB,  s + t == 1     virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B ) const;
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	  const tag::Winding &, const Manifold::Action & exp                         ) const ;
	void pretty_interpolate
		( const Cell & P, double s, const Cell & A, double t, const Cell & B ) const;

	// P = sA + sB + uC,  s+t+u == 1
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	  double u, Cell::Positive::Vertex * C                                       ) const;

	// P = sA + sB + uC,  s+t+u == 1
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A,
	  double t, Cell::Positive::Vertex * B,
	  const tag::Winding &, const Manifold::Action &,
	  double u, Cell::Positive::Vertex * C,
	  const tag::Winding &, const Manifold::Action & ) const;

	// P = sA + sB + uC + vD,  s + t + u + v == 1     virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
		double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
		double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D  ) const;

	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A,
	  double t, Cell::Positive::Vertex * B,
	  const tag::Winding &, const Manifold::Action &,
	  double u, Cell::Positive::Vertex * C,
	  const tag::Winding &, const Manifold::Action &,
	  double v, Cell::Positive::Vertex * D,
	  const tag::Winding &, const Manifold::Action & ) const;

	void pretty_interpolate
	( const Cell & P, double s, const Cell & A, double t, const Cell & B,
	                  double u, const Cell & C, double v, const Cell & D  ) const;

	// P = sA + sB + uC + vD + wE + zF,  s + t + u + v + w + z == 1    virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	  double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D,
	  double w, Cell::Positive::Vertex * E, double z, Cell::Positive::Vertex * F ) const;
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A,
	  double t, Cell::Positive::Vertex * B,
	  const tag::Winding &, const Manifold::Action &,
	  double u, Cell::Positive::Vertex * C,
	  const tag::Winding &, const Manifold::Action &,
	  double v, Cell::Positive::Vertex * D,
	  const tag::Winding &, const Manifold::Action &,
	  double w, Cell::Positive::Vertex * E,
	  const tag::Winding &, const Manifold::Action &,
	  double z, Cell::Positive::Vertex * F,
	  const tag::Winding &, const Manifold::Action & ) const;
	void pretty_interpolate ( const Cell & P, double s, const Cell & A,
	  double t, const Cell & B, double u, const Cell & C, double v, const Cell & D,
	  double w, const Cell & E, double z, const Cell & F ) const;

	// P = sum c_k P_k,  sum c_k == 1     virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
	  const std::vector < double > & coefs, const std::vector < Cell::Positive::Vertex * > & points ) const;
	void pretty_interpolate ( const Cell & P,
	  const std::vector < double > & coefs, const std::vector < Cell > & points ) const;
	
  Function build_coord_func ( const tag::lagrange &, const tag::OfDegree &, size_t d ) ;
  // virtual from Manifold::Core
	
  Function get_coord_func ( ) const;  // virtual from Manifold::Core
	
	void set_coords ( const Function co );  // virtual from Manifold::Core

	// measure of the entire manifold, here produces run-time error
	// (a Euclidian maifold has infinite measure)
	double measure ( ) const;  // virtual from Manifold::Core

	void project ( Cell::Positive::Vertex * ) const;
	// virtual from Manifold::Core, here does nothing

	// several versions of  frontal_method  defined by Manifold::Core

};  // end of class Manifold::Euclid

//-----------------------------------------------------------------------------------------


// when we first declare a Manifold::Euclid, coord_field and coord_func will be nullptr
// immediately after declaring the manifold (in main),
// we must declare the coordinates to be of type Lagrange degree 1 for instance
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

//-----------------------------------------------------------------------------------------


class Manifold::Implicit : public Manifold::Core

// a submanifold of a Manifold::Euclid defined by one or more equations

{	public :

	Manifold surrounding_space;

	inline Implicit ( )	:	Manifold::Core(), surrounding_space ( tag::non_existent ) { }
		
	// the projection will be done by means of the Newton method
	static const short int steps_for_Newton = 10;
	
	// P = sA + sB,  s + t == 1     virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B ) const;
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	  const tag::Winding &, const Manifold::Action & exp                    ) const ;

	// P = sA + sB + uC,  s+t+u == 1
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	  double u, Cell::Positive::Vertex * C                                       ) const;

	// P = sA + sB + uC,  s+t+u == 1
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A,
	  double t, Cell::Positive::Vertex * B,
	  const tag::Winding &, const Manifold::Action &,
	  double u, Cell::Positive::Vertex * C,
	  const tag::Winding &, const Manifold::Action & ) const;

	// P = sA + sB + uC + vD,  s + t + u + v == 1     virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	  double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D  ) const;
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A,
	  double t, Cell::Positive::Vertex * B,
	  const tag::Winding &, const Manifold::Action &,
	  double u, Cell::Positive::Vertex * C,
	  const tag::Winding &, const Manifold::Action &,
	  double v, Cell::Positive::Vertex * D,
	  const tag::Winding &, const Manifold::Action & ) const;

	// P = sA + sB + uC + vD + wE + zF,  s + t + u + v + w + z == 1     virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	  double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D,
	  double w, Cell::Positive::Vertex * E, double z, Cell::Positive::Vertex * F ) const;
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A,
	  double t, Cell::Positive::Vertex * B,
	  const tag::Winding &, const Manifold::Action &,
	  double u, Cell::Positive::Vertex * C,
	  const tag::Winding &, const Manifold::Action &,
	  double v, Cell::Positive::Vertex * D,
	  const tag::Winding &, const Manifold::Action &,
	  double w, Cell::Positive::Vertex * E,
	  const tag::Winding &, const Manifold::Action &,
	  double z, Cell::Positive::Vertex * F,
	  const tag::Winding &, const Manifold::Action & ) const;

	// P = sum c_k P_k,  sum c_k == 1     virtual from Manifold::Core
  void interpolate ( Cell::Positive::Vertex * P,
	  const std::vector < double > & coefs,
	  const std::vector < Cell::Positive::Vertex * > & points ) const;
	
	Function build_coord_func ( const tag::lagrange &, const tag::OfDegree &, size_t d );
	//   virtual from Manifold::Core, here execution forbidden
	
	Function get_coord_func ( ) const;  // virtual from Manifold::Core

	void set_coords ( const Function co );  // virtual from Manifold::Core

	// void project ( Cell::Positive::Vertex * ) const  stays pure virtual from Manifold::Core

	// measure of the entire manifold, here produces run-time error
	double measure ( ) const;  // virtual from Manifold::Core

	class OneEquation;  class TwoEquations;  class ThreeEquationsOrMore;
	
	// several versions of  frontal_method  defined by Manifold::Core
	
};  // end of class Manifold::Implicit

//-----------------------------------------------------------------------------------------


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

	// double measure ( ) const  virtual from Manifold::Core
	// defined by Manifold::Implicit, execution forbidden

	// several versions of  frontal_method  defined by Manifold::Core
	
};  // end of class Manifold::Implicit::OneEquation

//-----------------------------------------------------------------------------------------


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

	// double measure ( ) const  virtual from Manifold::Core
	// defined by Manifold::Implicit, execution forbidden

	// several versions of  frontal_method  defined by Manifold::Core
	
};  // end of class Manifold::Implicit::TwoEquations

//-----------------------------------------------------------------------------------------


class Manifold::Implicit::ThreeEquationsOrMore : public Manifold::Implicit

// a submanifold of a Manifold::Euclid defined by three or more equations

{	public :

	// Manifold surrounding_space  in Manifold::Implicit
	std::vector < Function > level_function, grad_lev_func;
	
	inline ThreeEquationsOrMore ( const Manifold & s, const Function & f );
	inline ThreeEquationsOrMore ( const Manifold & s, const Function & f1, const Function & f2 );
	inline ThreeEquationsOrMore
	( const Manifold & s, const Function & f1, const Function & f2, const Function & f3 );

	// void interpolate (different overloaded versions) defined by Manifold::Implicit
	
	// Function build_coord_func ( const tag::lagrange &, const tag::OfDegree &, size_t d )
	//   defined by Manifold::Implicit (execution forbidden)
	
	// Function get_coord_func ( ) const  defined by Manifold::Implicit

	// void set_coords ( const Function co )  defined by Manifold::Implicit

	void project ( Cell::Positive::Vertex * ) const;
	// virtual from Manifold::Core, through Manifold::Implicit

	// double measure ( ) const  virtual from Manifold::Core
	// defined by Manifold::Implicit, execution forbidden

	// several versions of  frontal_method  defined by Manifold::Core
	
};  // end of class Manifold::Implicit::ThreeEquationsOrMore

//-----------------------------------------------------------------------------------------


inline Manifold::Manifold ( const tag::Implicit &, const Manifold & m, const Function & f )

:	Manifold ( tag::non_existent )  // temporarily empty manifold

{	// if m is Manifold::Euclid, we want a Manifold::Implicit::OneEquation
	// if m is Manifold::Implicit::OneEquation, we want a Manifold::Implicit::TwoEquations
	// if m is Manifold::Implicit::TwoEquations, we want a Manifold::Implicit::ThreeEquationsOrMore
	// if m is Manifold::Implicit::ThreeEquationsOrMore, we want a Manifold::Implicit::ThreeEquationsOrMore

	Manifold::Euclid * m_euclid = dynamic_cast < Manifold::Euclid* > ( m .core );
	if ( m_euclid )
	{	this->core = new Manifold::Implicit::OneEquation ( m, f );
		Manifold::working = *this;
		return;                                                    }

	Manifold::Implicit::OneEquation * m_one_eq =
		dynamic_cast < Manifold::Implicit::OneEquation* > ( m .core );
	if ( m_one_eq )
	{	this->core = new Manifold::Implicit::TwoEquations ( m, f );
		Manifold::working = *this;
		return;                                                     }
	
	Manifold::Implicit::TwoEquations * m_two_eq =
		dynamic_cast < Manifold::Implicit::TwoEquations* > ( m .core );
	if ( m_two_eq )
	{	this->core = new Manifold::Implicit::ThreeEquationsOrMore ( m, f );
		Manifold::working = *this;
		return;                                                             }
	
	Manifold::Implicit::ThreeEquationsOrMore * m_three_eq =
		dynamic_cast < Manifold::Implicit::ThreeEquationsOrMore* > ( m .core );
	if ( m_three_eq )
	{	this->core = new Manifold::Implicit::ThreeEquationsOrMore ( m, f );
		Manifold::working = *this;
		return;                                                             }
	
	assert ( false );
	
}  // end of  Manifold constructor
	

inline Manifold::Manifold
( const tag::Implicit &, const Manifold & m, const Function & f1, const Function & f2 )

:	Manifold ( tag::non_existent )  // temporarily empty manifold

{	// if m is Manifold::Euclid, we want a Manifold::Implicit::TwoEquations
	// if m is Manifold::Implicit::{OneEquation,TwoEquations,ThreeEquationsOrMore},
	//      we want a Manifold::Implicit::TwoEquations

	Manifold::Euclid * m_euclid = dynamic_cast < Manifold::Euclid* > ( m .core );
	if ( m_euclid )
	{	this->core = new Manifold::Implicit::TwoEquations ( m, f1, f2 );
		Manifold::working = *this;
		return;                                                          }

	Manifold::Implicit::OneEquation * m_one_eq =
		dynamic_cast < Manifold::Implicit::OneEquation* > ( m .core );
	if ( m_one_eq )
	{	this->core = new Manifold::Implicit::ThreeEquationsOrMore ( m, f1, f2 );
		Manifold::working = *this;
		return;                                                                  }
	
	Manifold::Implicit::TwoEquations * m_two_eq =
		dynamic_cast < Manifold::Implicit::TwoEquations* > ( m .core );
	if ( m_two_eq )
	{	this->core = new Manifold::Implicit::ThreeEquationsOrMore ( m, f1, f2 );
		Manifold::working = *this;
		return;                                                                  }
	
	Manifold::Implicit::ThreeEquationsOrMore * m_three_eq =
		dynamic_cast < Manifold::Implicit::ThreeEquationsOrMore* > ( m .core );
	if ( m_three_eq )
	{	this->core = new Manifold::Implicit::ThreeEquationsOrMore ( m, f1, f2 );
		Manifold::working = *this;
		return;                                                                  }
	
	assert ( false );
	
}  // end of  Manifold constructor
		

inline Manifold::Manifold ( const tag::Implicit &, const Manifold & m,
                            const Function & f1, const Function & f2, const Function & f3 )

:	Manifold ( tag::non_existent )  // temporarily empty manifold

{	this->core = new Manifold::Implicit::ThreeEquationsOrMore ( m, f1, f2, f3 );
	Manifold::working = *this;                                                   }
		

inline Manifold::Implicit::OneEquation::OneEquation
( const Manifold & m, const Function & f )

: level_function ( f ),
	grad_lev_func ( 0. )  // temporarily zero gradient

{	this->surrounding_space = m;
	Manifold::Euclid * m_euclid = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Euclid* > ( m .core );
	Function coord = m_euclid->coord_func;
	size_t n = coord .nb_of_components();
	Function::Aggregate * grad = new Function::Aggregate ( tag::reserve_size, n );
	for ( size_t i = 0; i < n; i++ )  // grad->components [i] = f .deriv ( coord [i] );
		grad->components .emplace_back ( f .deriv ( coord [i] ) );
	this->grad_lev_func = Function ( tag::whose_core_is, grad );                        }

//-----------------------------------------------------------------------------------------


inline void Manifold::build_intersection_of_two_manif ( const Manifold & m1, const Manifold & m2 )

{	Manifold::Euclid * m1_euclid = dynamic_cast < Manifold::Euclid* > ( m1 .core );
	if ( m1_euclid )
	{	Manifold::Implicit * m2_implicit = dynamic_cast < Manifold::Implicit* > ( m2 .core );
		assert ( m2_implicit );
		assert ( m2_implicit->surrounding_space .core == m1 .core );
		this->core = m2 .core;
		Manifold::working = *this;
		return;                                                    }

	Manifold::Euclid * m2_euclid = dynamic_cast < Manifold::Euclid* > ( m2 .core );
	if ( m2_euclid )
	{	Manifold::Implicit * m1_implicit = dynamic_cast < Manifold::Implicit* > ( m1 .core );
		assert ( m1_implicit );
		assert ( m1_implicit->surrounding_space .core == m2 .core );
		this->core = m1 .core;
		Manifold::working = *this;
		return;                                                    }

	Manifold::Implicit * m1_implicit = dynamic_cast < Manifold::Implicit* > ( m1 .core );
	assert ( m1_implicit );
	Manifold::Implicit * m2_implicit = dynamic_cast < Manifold::Implicit* > ( m2 .core );
	assert ( m2_implicit );
	assert ( m1_implicit->surrounding_space .core == m2_implicit->surrounding_space .core );

	Manifold::Implicit::OneEquation * m1_implicit_one =
		dynamic_cast < Manifold::Implicit::OneEquation* > ( m1 .core );
	if ( m1_implicit_one )
	{	Manifold::Implicit::OneEquation * m2_implicit_one =
			dynamic_cast < Manifold::Implicit::OneEquation* > ( m2 .core );
		if ( m2_implicit_one )
		{	this->core = new Manifold::Implicit::TwoEquations ( m1, m2_implicit_one->level_function );
			Manifold::working = *this;
			return;                                                                                    }
		Manifold::Implicit::OneEquation * m2_implicit_two =
			dynamic_cast < Manifold::Implicit::OneEquation* > ( m2 .core );
		assert ( m2_implicit_two );
		this->core = new Manifold::Implicit::TwoEquations ( m2, m1_implicit_one->level_function );
		Manifold::working = *this;
		return;                                                                                         }

	Manifold::Implicit::TwoEquations * m1_implicit_two =
		dynamic_cast < Manifold::Implicit::TwoEquations* > ( m1 .core );
	assert ( m1_implicit_two );
	Manifold::Implicit::OneEquation * m2_implicit_one =
		dynamic_cast < Manifold::Implicit::OneEquation* > ( m2 .core );
	assert ( m2_implicit_one );
	this->core = new Manifold::Implicit::TwoEquations ( m1, m2_implicit_one->level_function );
	Manifold::working = *this;
	return;

}  // end of  build_intersection_two_manif
	

inline Manifold::Manifold ( const tag::Intersect &, const Manifold & m1, const Manifold & m2 )

:	Manifold ( tag::non_existent )  // temporarily empty manifold

{	this->build_intersection_of_two_manif ( m1, m2 );  }


inline Manifold::Manifold
( const tag::Intersect &, const Manifold & m1, const Manifold & m2, const Manifold & m3 )

:	Manifold ( tag::non_existent )  // temporarily empty manifold

{	Manifold::Euclid * m1_euclid = dynamic_cast < Manifold::Euclid* > ( m1 .core );
	if ( m1_euclid )
	{	this->build_intersection_of_two_manif ( m2, m3 );
		return;                                           }

	Manifold::Euclid * m2_euclid = dynamic_cast < Manifold::Euclid* > ( m2 .core );
	if ( m2_euclid )
	{	this->build_intersection_of_two_manif ( m1, m3 );
		return;                                           }
	
	Manifold::Euclid * m3_euclid = dynamic_cast < Manifold::Euclid* > ( m3 .core );
	if ( m3_euclid )
	{	this->build_intersection_of_two_manif ( m1, m2 );
		return;                                           }

	Manifold::Implicit::OneEquation * m1_implicit_one =
		dynamic_cast < Manifold::Implicit::OneEquation* > ( m1 .core );
	assert ( m1_implicit_one );
	Manifold::Implicit::OneEquation * m2_implicit_one =
		dynamic_cast < Manifold::Implicit::OneEquation* > ( m2 .core );
	assert ( m2_implicit_one );
	Manifold::Implicit::OneEquation * m3_implicit_one =
		dynamic_cast < Manifold::Implicit::OneEquation* > ( m3 .core );
	assert ( m3_implicit_one );

	this->core = new Manifold::Implicit::ThreeEquationsOrMore
		( m1, m2_implicit_one->level_function, m3_implicit_one->level_function );
	Manifold::working = *this;
	return;                                                                         }
	
	
//-----------------------------------------------------------------------------------------


inline Manifold::Implicit::TwoEquations::TwoEquations
( const Manifold & m, const Function & f )

: level_function_1 ( 0. ),  // temporarily zero function
	grad_lev_func_1 ( 0. ),  // temporarily zero gradient
	level_function_2 ( f ),
	grad_lev_func_2 ( 0. )  // temporarily zero gradient

{	Manifold::Implicit::OneEquation * m_one_eq = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Implicit::OneEquation * > ( m .core );
	this->surrounding_space = m_one_eq->surrounding_space;
	this->level_function_1 = m_one_eq->level_function;
	this->grad_lev_func_1 = m_one_eq->grad_lev_func;
	Function coord = this->surrounding_space .coordinates();
	size_t n = coord .nb_of_components();
	Function::Aggregate * grad = new Function::Aggregate ( tag::reserve_size, n );
	for ( size_t i = 0; i < n; i++ ) // grad->components[i] = f .deriv ( coord [i] );
		grad->components .emplace_back ( f .deriv ( coord [i] ) );
	this->grad_lev_func_2 = Function ( tag::whose_core_is, grad );                    }


inline Manifold::Implicit::TwoEquations::TwoEquations
( const Manifold & m, const Function & f1, const Function & f2 )

: level_function_1 ( f1 ), 
	grad_lev_func_1 ( 0. ),  // temporarily zero gradient
	level_function_2 ( f2 ),
	grad_lev_func_2 ( 0. )  // temporarily zero gradient

{	this->surrounding_space = m;
	Manifold::Euclid * m_euclid = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Euclid* > ( m .core );
	Function coord = m_euclid->coord_func;
	size_t n = coord .nb_of_components();
	Function::Aggregate * grad = new Function::Aggregate ( tag::reserve_size, n );
	for ( size_t i = 0; i < n; i++ ) // grad->components[i] = f1.deriv(coord[i]);
		grad->components .emplace_back ( f1 .deriv ( coord [i] ) );
	this->grad_lev_func_1 = Function ( tag::whose_core_is, grad );
	grad = new Function::Aggregate ( tag::reserve_size, n );
	for ( size_t i = 0; i < n; i++ ) // grad->components [i] = f2 .deriv ( coord [i] );
		grad->components .emplace_back ( f2 .deriv ( coord [i] ) );
	this->grad_lev_func_2 = Function ( tag::whose_core_is, grad );                      }


inline Manifold::Implicit::ThreeEquationsOrMore::ThreeEquationsOrMore
( const Manifold & m, const Function & f )

// level_function and grad_lev_function initialized by default as empty vectors

{	Manifold::Implicit::TwoEquations * m_two_eq = dynamic_cast
		< Manifold::Implicit::TwoEquations* > ( m .core );
	if ( m_two_eq )
	{	this->surrounding_space = m_two_eq->surrounding_space;
		this->level_function .push_back ( m_two_eq->level_function_1 );
		this->level_function .push_back ( m_two_eq->level_function_2 );
		this->grad_lev_func .push_back ( m_two_eq->grad_lev_func_1 );
		this->grad_lev_func .push_back ( m_two_eq->grad_lev_func_2 );    }
	else
	{	Manifold::Implicit::ThreeEquationsOrMore * m_three_eq = tag::Util::assert_cast
			< Manifold::Core*, Manifold::Implicit::ThreeEquationsOrMore* > ( m .core );
		this->surrounding_space = m_three_eq->surrounding_space;
		this->level_function = m_three_eq->level_function;
		this->grad_lev_func = m_three_eq->grad_lev_func;                                }
	this->level_function .push_back ( f );
	Function coord = this->surrounding_space .coordinates();
	size_t n = coord .nb_of_components();
	Function::Aggregate * grad = new Function::Aggregate ( tag::reserve_size, n );
	for ( size_t i = 0; i < n; i++ ) // grad->components[i] = f .deriv ( coord [i] );
		grad->components .emplace_back ( f .deriv ( coord [i] ) );
	// this->grad_lev_func .emplace_back ( std::piecewise_construct,
	//                                     std::forward_as_tuple ( tag::whose_core_is ),
	//                                     std::forward_as_tuple ( grad )               );
	this->grad_lev_func .push_back ( Function ( tag::whose_core_is, grad ) );           }


inline Manifold::Implicit::ThreeEquationsOrMore::ThreeEquationsOrMore
( const Manifold & m, const Function & f1, const Function & f2 )

// level_function and grad_lev_function initialized by default as empty vectors

{	Manifold::Implicit::OneEquation * m_one_eq = dynamic_cast
		< Manifold::Implicit::OneEquation* > ( m .core );
	Manifold::Implicit::TwoEquations * m_two_eq = dynamic_cast
		< Manifold::Implicit::TwoEquations* > ( m .core );
	if ( m_one_eq )
	{	this->surrounding_space = m_one_eq->surrounding_space;
		this->level_function .push_back ( m_one_eq->level_function );
		this->grad_lev_func .push_back ( m_one_eq->grad_lev_func );   }
	else if ( m_two_eq )
	{	this->surrounding_space = m_two_eq->surrounding_space;
		this->level_function .push_back ( m_two_eq->level_function_1 );
		this->level_function .push_back ( m_two_eq->level_function_2 );
		this->grad_lev_func .push_back ( m_two_eq->grad_lev_func_1 );
		this->grad_lev_func .push_back ( m_two_eq->grad_lev_func_2 );    }
	else
	{	Manifold::Implicit::ThreeEquationsOrMore * m_three_eq = tag::Util::assert_cast
			< Manifold::Core*, Manifold::Implicit::ThreeEquationsOrMore* > ( m .core );
		this->surrounding_space = m_three_eq->surrounding_space;
		this->level_function = m_three_eq->level_function;
		this->grad_lev_func = m_three_eq->grad_lev_func;                                }
	this->level_function .push_back ( f1 );
	this->level_function .push_back ( f2 );
	Function coord = this->surrounding_space .coordinates();
	size_t n = coord .nb_of_components();
	Function::Aggregate * grad = new Function::Aggregate ( tag::reserve_size, n );
	for ( size_t i = 0; i < n; i++ ) // grad->components[i] = f1 .deriv ( coord [i] );
		grad->components .emplace_back ( f1 .deriv ( coord [i] ) );
	this->grad_lev_func .push_back ( Function ( tag::whose_core_is, grad ) );
	grad = new Function::Aggregate ( tag::reserve_size, n );
	for ( size_t i = 0; i < n; i++ ) // grad->components[i] = f2 .deriv ( coord [i] );
		grad->components .emplace_back ( f2 .deriv ( coord [i] ) );
	this->grad_lev_func .push_back ( Function ( tag::whose_core_is, grad ) );           }


inline Manifold::Implicit::ThreeEquationsOrMore::ThreeEquationsOrMore
( const Manifold & m, const Function & f1, const Function & f2, const Function & f3 )

// level_function and grad_lev_function initialized by default as empty vectors

{	Manifold::Euclid * m_euclid = dynamic_cast < Manifold::Euclid* > ( m .core );
	Manifold::Implicit::OneEquation * m_one_eq = dynamic_cast
		< Manifold::Implicit::OneEquation* > ( m .core );
	Manifold::Implicit::TwoEquations * m_two_eq = dynamic_cast
		< Manifold::Implicit::TwoEquations* > ( m .core );
	if ( m_euclid ) this->surrounding_space = m;
	else if ( m_one_eq )
	{	this->surrounding_space = m_one_eq->surrounding_space;
		this->level_function .push_back ( m_one_eq->level_function );
		this->grad_lev_func .push_back ( m_one_eq->grad_lev_func );   }
	else if ( m_two_eq )
	{	this->surrounding_space = m_two_eq->surrounding_space;
		this->level_function .push_back ( m_two_eq->level_function_1 );
		this->level_function .push_back ( m_two_eq->level_function_2 );
		this->grad_lev_func .push_back ( m_two_eq->grad_lev_func_1 );
		this->grad_lev_func .push_back ( m_two_eq->grad_lev_func_2 );    }
	else
	{	Manifold::Implicit::ThreeEquationsOrMore * m_three_eq = tag::Util::assert_cast
			< Manifold::Core*, Manifold::Implicit::ThreeEquationsOrMore* > ( m .core );
		this->surrounding_space = m_three_eq->surrounding_space;
		this->level_function = m_three_eq->level_function;
		this->grad_lev_func = m_three_eq->grad_lev_func;                                }
	this->level_function .push_back ( f1 );
	this->level_function .push_back ( f2 );
	this->level_function .push_back ( f3 );
	Function coord = this->surrounding_space .coordinates();
	size_t n = coord .nb_of_components();
	Function::Aggregate * grad = new Function::Aggregate ( tag::reserve_size, n );
	for ( size_t i = 0; i < n; i++ ) // grad->components[i] = f1 .deriv ( coord [i] );
		grad->components .emplace_back ( f1 .deriv ( coord [i] ) );
	this->grad_lev_func .push_back ( Function ( tag::whose_core_is, grad ) );
	grad = new Function::Aggregate ( tag::reserve_size, n );
	for ( size_t i = 0; i < n; i++ ) // grad->components[i] = f2 .deriv ( coord [i] );
		grad->components .emplace_back ( f2 .deriv ( coord [i] ) );
	this->grad_lev_func .push_back ( Function ( tag::whose_core_is, grad ) );
	grad = new Function::Aggregate ( tag::reserve_size, n );
	for ( size_t i = 0; i < n; i++ ) // grad->components[i] = f3 .deriv ( coord [i] );
		grad->components .emplace_back ( f3 .deriv ( coord [i] ) );
	this->grad_lev_func .push_back ( Function ( tag::whose_core_is, grad ) );           }

//-----------------------------------------------------------------------------------------


class Manifold::Parametric : public Manifold::Core

// a submanifold of a Manifold::Euclid defined by one or more explicit equations

{	public :

	Manifold surrounding_space;

	std::map < Function, Function, bool (*) ( const Function &, const Function & ) > equations;
	//  decltype(Function::less_for_map)*  equals  bool (*) ( const Function &, const Function & )

	inline Parametric ( )
	: Manifold::Core(),
		surrounding_space ( tag::non_existent ),
		equations ( & Function::less_for_map )
	{ }

	inline Parametric ( const Manifold &, const Function::Equality & );

	inline Parametric ( const Manifold &, const Function::Equality &, const Function::Equality & );

	inline Parametric ( const Manifold &, const Function::Equality &,
                      const Function::Equality &, const Function::Equality & );

	// P = sA + sB,  s + t == 1     virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B ) const;
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	  const tag::Winding &, const Manifold::Action & exp                    ) const ;

	// P = sA + sB + uC,  s+t+u == 1
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	  double u, Cell::Positive::Vertex * C                                       ) const;

	// P = sA + sB + uC,  s+t+u == 1
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A,
	  double t, Cell::Positive::Vertex * B,
	  const tag::Winding &, const Manifold::Action &,
	  double u, Cell::Positive::Vertex * C,
	  const tag::Winding &, const Manifold::Action & ) const;

	// P = sA + sB + uC + vD,  s + t + u + v == 1     virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	  double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D  ) const;
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A,
	  double t, Cell::Positive::Vertex * B,
	  const tag::Winding &, const Manifold::Action &,
	  double u, Cell::Positive::Vertex * C,
	  const tag::Winding &, const Manifold::Action &,
	  double v, Cell::Positive::Vertex * D,
	  const tag::Winding &, const Manifold::Action & ) const;

	// P = sA + sB + uC + vD + wE + zF,  s + t + u + v + w + z == 1    virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	  double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D,
	  double w, Cell::Positive::Vertex * E, double z, Cell::Positive::Vertex * F ) const;
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A,
	  double t, Cell::Positive::Vertex * B,
	  const tag::Winding &, const Manifold::Action &,
	  double u, Cell::Positive::Vertex * C,
	  const tag::Winding &, const Manifold::Action &,
	  double v, Cell::Positive::Vertex * D,
	  const tag::Winding &, const Manifold::Action &,
	  double w, Cell::Positive::Vertex * E,
	  const tag::Winding &, const Manifold::Action &,
	  double z, Cell::Positive::Vertex * F,
	  const tag::Winding &, const Manifold::Action & ) const;

	// P = sum c_k P_k,  sum c_k == 1     virtual from Manifold::Core
	virtual void interpolate ( Cell::Positive::Vertex * P,
	  const std::vector < double > & coefs, const std::vector < Cell::Positive::Vertex * > & points ) const;
	
	Function build_coord_func ( const tag::lagrange &, const tag::OfDegree &, size_t d );
	//   virtual from Manifold::Core, here execution forbidden
	
	Function get_coord_func ( ) const;  // virtual from Manifold::Core

	void set_coords ( const Function co );  // virtual from Manifold::Core

	void project ( Cell::Positive::Vertex * ) const;
	// virtual from Manifold::Core

	// measure of the entire manifold, here produces run-time error
	double measure ( ) const;  // virtual from Manifold::Core

};  // end of class Manifold::Parametric

//-----------------------------------------------------------------------------------------


inline Manifold::Manifold
( const tag::Parametric &, const Manifold & m, const Function::Equality & f_eq )

:	Manifold ( tag::non_existent )  // temporary empty manifold

{	tag::Util::assert_cast < Manifold::Core*, Manifold::Euclid* > ( m.core );
	this->core = new Manifold::Parametric ( m, f_eq );
	Manifold::working = *this;                                                 }


inline Manifold::Manifold ( const tag::Parametric &, const Manifold & m,
                            const Function::Equality & f_eq_1, const Function::Equality & f_eq_2 )

:	Manifold ( tag::non_existent )  // temporary empty manifold

{	tag::Util::assert_cast < Manifold::Core*, Manifold::Euclid* > ( m.core );
	this->core = new Manifold::Parametric ( m, f_eq_1, f_eq_2 );
	Manifold::working = *this;                                                 }


inline Manifold::Manifold
( const tag::Parametric &, const Manifold & m, const Function::Equality & f_eq_1,
  const Function::Equality & f_eq_2, const Function::Equality & f_eq_3           )

:	Manifold ( tag::non_existent )  // temporary empty manifold

{	tag::Util::assert_cast < Manifold::Core*, Manifold::Euclid* > ( m.core );
	this->core = new Manifold::Parametric ( m, f_eq_1, f_eq_2, f_eq_3 );
	Manifold::working = *this;                                                           }


inline Manifold::Parametric::Parametric ( const Manifold & m, const Function::Equality & f_eq )

:	Parametric()

{	this->surrounding_space = m;
	tag::Util::assert_cast < Manifold::Core*, Manifold::Euclid* > ( m .core );
	assert ( this->equations .find ( f_eq .lhs ) == this->equations .end() );
	this->equations .insert ( std::pair ( f_eq .lhs, f_eq .rhs ) );               }
//	this->equations [ f_eq.lhs ] = f_eq .rhs;


inline Manifold::Parametric::Parametric
( const Manifold & m, const Function::Equality & f_eq_1, const Function::Equality & f_eq_2 )

:	Parametric()

{	this->surrounding_space = m;
	tag::Util::assert_cast < Manifold::Core*, Manifold::Euclid* > ( m .core );
	assert ( this->equations .find ( f_eq_1 .lhs ) == this->equations .end() );
	this->equations .insert ( std::pair ( f_eq_1 .lhs, f_eq_1 .rhs ) );
	assert ( this->equations.find ( f_eq_2 .lhs ) == this->equations .end() );
	this->equations .insert ( std::pair ( f_eq_2 .lhs, f_eq_2 .rhs ) );           }
//	this->equations [ f_eq_1 .lhs ] = f_eq_1 .rhs;
//	this->equations [ f_eq_2 .lhs ] = f_eq_2 .rhs;

inline Manifold::Parametric::Parametric
( const Manifold & m, const Function::Equality & f_eq_1,
  const Function::Equality & f_eq_2, const Function::Equality & f_eq_3 )

:	Parametric()

{	this->surrounding_space = m;
	tag::Util::assert_cast < Manifold::Core*, Manifold::Euclid* > ( m .core );
	assert ( this->equations .find ( f_eq_1.lhs ) == this->equations .end() );
	this->equations .insert ( std::pair ( f_eq_1 .lhs, f_eq_1 .rhs ) );
	assert ( this->equations .find ( f_eq_2.lhs ) == this->equations .end() );
	this->equations .insert ( std::pair ( f_eq_2 .lhs, f_eq_2 .rhs ) );
	assert ( this->equations .find ( f_eq_3.lhs ) == this->equations .end() );
	this->equations .insert ( std::pair ( f_eq_3 .lhs, f_eq_3 .rhs ) );            }
//	this->equations [ f_eq_1 .lhs ] = f_eq_1 .rhs;
//	this->equations [ f_eq_2 .lhs ] = f_eq_2 .rhs;
//	this->equations [ f_eq_3 .lhs ] = f_eq_3 .rhs;

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class Manifold::Quotient : public Manifold::Core

// a Euclidian manifold divided by a group of actions

{	public :

	Manifold base_space;  // wrap around a Manifold::Euclid
	
	Function coord_func { tag::non_existent };

	std::vector < Function::ActionGenerator > actions;  // set of generators for a discrete group
	
	std::vector < Field::ShortInt > winding_nbs;  // a jump (exponent) per action

	inline Quotient ( Manifold b, const Manifold::Action & g1 );

	inline Quotient ( Manifold b, const Manifold::Action & g1, const Manifold::Action & g2 );

	inline Quotient ( Manifold b,
	  const Manifold::Action & g1, const Manifold::Action & g2, const Manifold::Action & g3 );

	Function build_coord_func ( const tag::lagrange &, const tag::OfDegree &, size_t d );
	//   virtual from Manifold::Core, here execution forbidden
	
	Function get_coord_func ( ) const;  // virtual from Manifold::Core

	void set_coords ( const Function co );  // virtual from Manifold::Core

	void project ( Cell::Positive::Vertex * ) const;  // virtual from Manifold::Core

	// P = sA + sB,  s + t == 1     virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B ) const;
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	  const tag::Winding &, const Manifold::Action & exp                         ) const ;

	// P = sA + sB + uC,  s+t+u == 1
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	  double u, Cell::Positive::Vertex * C                                       ) const;

	// P = sA + sB + uC,  s+t+u == 1
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A,
	  double t, Cell::Positive::Vertex * B,
	  const tag::Winding &, const Manifold::Action &,
	  double u, Cell::Positive::Vertex * C,
	  const tag::Winding &, const Manifold::Action & ) const;

	// P = sA + sB + uC + vD,  s + t + u + v == 1     virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	  double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D  ) const;
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A,
	  double t, Cell::Positive::Vertex * B,
	  const tag::Winding &, const Manifold::Action &,
	  double u, Cell::Positive::Vertex * C,
	  const tag::Winding &, const Manifold::Action &,
	  double v, Cell::Positive::Vertex * D,
	  const tag::Winding &, const Manifold::Action & ) const;

	// P = sA + sB + uC + vD + wE + zF,  s + t + u + v + w + z == 1    virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	  double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D,
	  double w, Cell::Positive::Vertex * E, double z, Cell::Positive::Vertex * F ) const;
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A,
	  double t, Cell::Positive::Vertex * B,
	  const tag::Winding &, const Manifold::Action &,
	  double u, Cell::Positive::Vertex * C,
	  const tag::Winding &, const Manifold::Action &,
	  double v, Cell::Positive::Vertex * D,
	  const tag::Winding &, const Manifold::Action &,
	  double w, Cell::Positive::Vertex * E,
	  const tag::Winding &, const Manifold::Action &,
	  double z, Cell::Positive::Vertex * F,
	  const tag::Winding &, const Manifold::Action & ) const;

	// P = sum c_k P_k,  sum c_k == 1     virtual from Manifold::Core
	virtual void interpolate ( Cell::Positive::Vertex * P,
	  const std::vector < double > & coefs,
	  const std::vector < Cell::Positive::Vertex * > & points ) const;
	
	// measure of the entire manifold
	double measure ( ) const;  // virtual from Manifold::Core
	
	void frontal_method  // virtual from Manifold::Core
	// defined in frontal.cpp, overrides definition by Manifold::Core
	( Mesh & msh, const tag::StartWithInconsistentMesh &,
	  const tag::StartAt &, const Cell & start,
	  const tag::Towards &, std::vector<double> tangent,
	  const tag::StopAt &,  const Cell & stop            ) override;
	
};  // end of class Manifold::Quotient

//-----------------------------------------------------------------------------------------


inline Manifold::Quotient::Quotient ( Manifold b, const Manifold::Action & a1 )

// a Field::ShortInt is created and kept in this->winding_nbs
// this increments Cell::Positive::short_int_heap_size [1]
// thus, future segments will have space reserved for a winding number
	
: Manifold::Core(), base_space ( b ), actions { }

{	assert ( a1.index_map.size() == 1 );
	std::map < Function::ActionGenerator, short int > ::const_iterator
		it = a1 .index_map .begin();
	assert ( it->second == 1 );
	const Function::ActionGenerator & g1 = it->first;
	this->actions = { g1 };

	assert ( g1 .coords .core == b .coordinates() .core );
	assert ( this->coord_func .core == nullptr );
	this->winding_nbs .emplace_back ( tag::lives_on_positive_cells, tag::of_dim, 1 );
	assert ( this->coord_func .core == nullptr );

	if ( b .coordinates() .nb_of_components() == 1 )
	{	std::vector < double > v1 =
			Function::Scalar::MultiValued::JumpIsSum::analyse_linear_expression
			( g1 .transf, g1 .coords );
		assert ( v1 .size() == 1 );
		this->coord_func .core = new Function::Scalar::MultiValued::JumpIsSum
			( tag::associated_with, b .coordinates(), { g1 }, v1 );             }
	else	
	{	assert ( b.coordinates().nb_of_components() > 1 );
		std::vector < double > v1 =
			Function::Vector::MultiValued::JumpIsSum::analyse_linear_expression
			( g1.transf, g1.coords );
		if ( v1.size() > 0 )
			this->coord_func.core = new Function::Vector::MultiValued::JumpIsSum
				( tag::associated_with, b.coordinates(), { g1 }, { v1 } );
		else
		{	std::pair < std::vector < std::vector < double > >, std::vector < double > >
			p = Function::Vector::MultiValued::JumpIsLinear::analyse_linear_expression
			    ( g1.transf, g1.coords );
			this->coord_func.core = new Function::Vector::MultiValued::JumpIsLinear
				( tag::associated_with, b.coordinates(), { g1 }, { p.first }, { p.second } );  }  }
	assert ( this->coord_func.core );
	this->coord_func.core->nb_of_wrappers = 1;                                                }


//-----------------------------------------------------------------------------------------


inline Manifold::Quotient::Quotient
( Manifold b, const Manifold::Action & a1, const Manifold::Action & a2 )

// two Field::ShortInt are created and kept in this->winding_nbs
// this increments by two Cell::Positive::short_int_heap_size [1]
// thus, future segments will have space reserved for two winding numbers
	
: Manifold::Core(), base_space ( b ), actions { }

{	assert ( a1 .index_map .size() == 1 );
	std::map < Function::ActionGenerator, short int > ::const_iterator
		it = a1 .index_map .begin();
	assert ( it->second == 1 );
	const Function::ActionGenerator & g1 = it->first;
	assert ( a2 .index_map .size() == 1 );
	it = a2 .index_map .begin();
	assert ( it->second == 1 );
	const Function::ActionGenerator & g2 = it->first;
	this->actions = { g1, g2 };

	assert ( g1 .coords .core == b .coordinates() .core );
	assert ( g2 .coords .core == b .coordinates() .core );
	this->winding_nbs .reserve ( 2 );
	this->winding_nbs .emplace_back ( tag::lives_on_positive_cells, tag::of_dim, 1 );
	this->winding_nbs .emplace_back ( tag::lives_on_positive_cells, tag::of_dim, 1 );
	assert ( b .coordinates() .nb_of_components() >= 2 );

	assert ( this->coord_func.core == nullptr );
	std::vector < double > v1 =
		Function::Vector::MultiValued::JumpIsSum::analyse_linear_expression
		( g1 .transf, g1 .coords );
	std::vector < double > v2 =
		Function::Vector::MultiValued::JumpIsSum::analyse_linear_expression
		( g2 .transf, g2 .coords );
	if ( v1.size() > 0 )
	{	assert ( v2.size() > 0 );
		this->coord_func .core = new Function::Vector::MultiValued::JumpIsSum
			( tag::associated_with, b .coordinates(), { g1, g2 }, { v1, v2 } );  }
	else
	{	std::pair < std::vector < std::vector < double > >, std::vector < double > >
		p1 = Function::Vector::MultiValued::JumpIsLinear::analyse_linear_expression
		     ( g1 .transf, g1 .coords ),
		p2 = Function::Vector::MultiValued::JumpIsLinear::analyse_linear_expression
		     ( g2 .transf, g2 .coords );
		this->coord_func .core = new Function::Vector::MultiValued::JumpIsLinear
			( tag::associated_with, b.coordinates(), { g1, g2 },
			  { p1 .first, p2 .first }, { p1 .second, p2 .second } );                  }
	assert ( this->coord_func .core );
	this->coord_func .core->nb_of_wrappers = 1;                                       }


//-----------------------------------------------------------------------------------------


inline Manifold::Quotient::Quotient ( Manifold b,
  const Manifold::Action & a1, const Manifold::Action & a2, const Manifold::Action & a3 )

// three Field::ShortInt are created and kept in this->winding_nbs
// this increments by three Cell::Positive::short_int_heap_size [1]
// thus, future segments will have space reserved for three winding numbers
	
: Manifold::Core(), base_space ( b ), actions { }

{	assert ( a1 .index_map .size() == 1 );
	std::map < Function::ActionGenerator, short int > ::const_iterator
		it = a1 .index_map .begin();
	assert ( it->second == 1 );
	const Function::ActionGenerator & g1 = it->first;
	assert ( a2 .index_map .size() == 1 );
	it = a2 .index_map .begin();
	assert ( it->second == 1 );
	const Function::ActionGenerator & g2 = it->first;
	assert ( a3 .index_map .size() == 1 );
	it = a3 .index_map .begin();
	assert ( it->second == 1 );
	const Function::ActionGenerator & g3 = it->first;
	this->actions = { g1, g2, g3 };

	assert ( g1 .coords .core == b .coordinates() .core );
	assert ( g2 .coords .core == b .coordinates() .core );
	assert ( g3 .coords .core == b .coordinates() .core );
	this->winding_nbs .reserve ( 3 );
	this->winding_nbs .emplace_back ( tag::lives_on_positive_cells, tag::of_dim, 1 );
	this->winding_nbs .emplace_back ( tag::lives_on_positive_cells, tag::of_dim, 1 );
	this->winding_nbs .emplace_back ( tag::lives_on_positive_cells, tag::of_dim, 1 );
	std::vector < double > v1 =
		Function::Vector::MultiValued::JumpIsSum::analyse_linear_expression
		( g1 .transf, g1 .coords );
	assert ( v1 .size() > 0 );
	std::vector < double > v2 =
		Function::Vector::MultiValued::JumpIsSum::analyse_linear_expression
		( g2 .transf, g2 .coords );
	assert ( v2 .size() > 0 );
	std::vector < double > v3 =
		Function::Vector::MultiValued::JumpIsSum::analyse_linear_expression
		( g3 .transf, g3 .coords );
	assert ( v3 .size() > 0 );
	assert ( b .coordinates() .nb_of_components() >= 3 );
	assert ( this->coord_func .core == nullptr );
	this->coord_func .core = new Function::Vector::MultiValued::JumpIsSum
		( tag::associated_with, b .coordinates(), { g1, g2, g3 }, { v1, v2, v3 } );
	this->coord_func .core->nb_of_wrappers = 1;                                    }

//-----------------------------------------------------------------------------------------


inline Manifold Manifold::quotient ( const Manifold::Action & a1 )
{	return Manifold ( tag::whose_core_is, new Manifold::Quotient ( *this, a1 ) );  }

inline Manifold Manifold::quotient ( const Manifold::Action & a1, const Manifold::Action & a2 )
{	return Manifold ( tag::whose_core_is, new Manifold::Quotient ( *this, a1, a2 ) );  }

inline Manifold Manifold::quotient
( const Manifold::Action & a1, const Manifold::Action & a2, const Manifold::Action & a3 )
{	return Manifold ( tag::whose_core_is, new Manifold::Quotient ( *this, a1, a2, a3 ) );  }

//-----------------------------------------------------------------------------------------


inline Cell::Winding::operator Manifold::Action ( )

{	Manifold::Action res ( 0 );
	Manifold::Quotient * manif_q = dynamic_cast
		< Manifold::Quotient* > ( Manifold::working.core );
	assert ( manif_q );

	assert ( this->cll->get_dim() == 1 );  // usually is a segment

	size_t n = manif_q->actions.size();
	assert ( n == manif_q->winding_nbs.size() );
	Cell::Core * pos_cll = this->cll->get_positive().core;
	for ( size_t i = 0; i < n; i++ )
	{	short int exp = manif_q->winding_nbs[i].on_cell(pos_cll);
		if ( exp == 0 ) continue;
		if ( not this->cll->is_positive() ) exp = -exp;
		Function::ActionGenerator & g = manif_q->actions[i];
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map<Function::ActionGenerator,short int>::iterator itt =
			res.index_map.lower_bound ( g );
		assert ( ( itt == res.index_map.end() ) or
	           ( res.index_map.key_comp()(g,itt->first) ) );
		res.index_map.emplace_hint ( itt, std::piecewise_construct,
			std::forward_as_tuple(g), std::forward_as_tuple(exp) );       }
	return res;                                                          }


inline Manifold::Action Cell::Winding::operator= ( const Manifold::Action & a )
	
{	Manifold::Quotient * manif_q = dynamic_cast
		< Manifold::Quotient* > ( Manifold::working.core );
	assert ( manif_q );

	assert ( this->cll->get_dim() == 1 );  // usually is a segment

	size_t n = manif_q->actions.size();
	assert ( n == manif_q->winding_nbs.size() );
	Cell::Core * pos_cll = this->cll->get_positive().core;
	for ( size_t i = 0; i < n; i++ )
	{	Function::ActionGenerator & g = manif_q->actions[i];
		std::map<Function::ActionGenerator,short int>::const_iterator itt =
			a.index_map.find ( g );
		if ( itt == a.index_map.end() )
		{	manif_q->winding_nbs[i] .on_cell ( pos_cll ) = 0;
			continue;                                   }
		short int exp = itt->second;
		assert ( exp != 0 );
		if ( this->cll->is_positive() ) manif_q->winding_nbs[i].on_cell(pos_cll) = exp;
		else manif_q->winding_nbs[i].on_cell(pos_cll) = -exp;                           }
	return *this;                                                                  }


inline Manifold::Action Cell::Winding::operator+= ( const Manifold::Action & a )
	
{	Manifold::Quotient * manif_q = dynamic_cast
		< Manifold::Quotient* > ( Manifold::working.core );
	assert ( manif_q );

	assert ( this->cll->get_dim() == 1 );  // usually is a segment

	size_t n = manif_q->actions.size();
	assert ( n == manif_q->winding_nbs.size() );
	Cell::Core * pos_cll = this->cll->get_positive().core;
	for ( size_t i = 0; i < n; i++ )
	{	Function::ActionGenerator & g = manif_q->actions[i];
		std::map<Function::ActionGenerator,short int>::const_iterator itt =
			a.index_map.find ( g );
		if ( itt == a.index_map.end() ) continue;
		short int exp = itt->second;
		assert ( exp != 0 );
		if ( not this->cll->is_positive() ) exp = -exp;
		manif_q->winding_nbs[i].on_cell(pos_cll) += exp;                           }
	return *this;                                                             }


inline Manifold::Action Cell::Winding::operator-= ( const Manifold::Action & a )
	
{	Manifold::Quotient * manif_q = dynamic_cast
		< Manifold::Quotient* > ( Manifold::working.core );
	assert ( manif_q );

	assert ( this->cll->get_dim() == 1 );  // usually is a segment

	size_t n = manif_q->actions.size();
	assert ( n == manif_q->winding_nbs.size() );
	Cell::Core * pos_cll = this->cll->get_positive().core;
	for ( size_t i = 0; i < n; i++ )
	{	Function::ActionGenerator & g = manif_q->actions[i];
		std::map<Function::ActionGenerator,short int>::const_iterator itt =
			a.index_map.find ( g );
		if ( itt == a.index_map.end() ) continue;
		short int exp = itt->second;
		assert ( exp != 0 );
		if ( not this->cll->is_positive() ) exp = -exp;
		manif_q->winding_nbs[i].on_cell(pos_cll) -= exp;                           }
	return *this;                                                             }


inline Cell::Winding Cell::winding ( ) const
{ return Cell::Winding ( *this );  }

//-----------------------------------------------------------------------------------------


inline void Cell::project () const
{	this->project ( tag::onto, Manifold::working );  }
	
inline void Cell::project ( const tag::Onto &, const Manifold m ) const
{	assert ( m.core );
	assert ( this->core->get_dim() == 0 );
	assert ( this->core->is_positive() );
	Cell::Positive::Vertex * cll = ( Cell::Positive::Vertex * ) this->core;
	m.core->project ( cll );                                                }
		
//-------------------------------------------------------------------------------------------------


inline Mesh::Mesh ( const tag::Frontal &, const tag::EntireManifold &, Manifold manif,
                    const tag::Orientation &, const tag::OrientationChoice & oc,
                    const tag::DesiredLength &, const Function & length               )
:	Mesh ( tag::frontal, tag::entire_manifold, manif,
	       tag::desired_length, length, tag::orientation, oc )
{	}

//-------------------------------------------------------------------------------------------------

}  // end of  namespace maniFEM

#endif
// ifndef MANIFEM_MANIFOLD_H

