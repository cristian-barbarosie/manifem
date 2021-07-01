
// function.h 2021.04.18

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

#ifndef MANIFEM_FUNCTION_H
#define MANIFEM_FUNCTION_H

#include <iostream>
#include <vector>
#include <forward_list>
#include <memory>
#include "math.h"
#include "assert.h"

#include "mesh.h"
#include "field.h"


namespace maniFEM {


namespace tag
{	struct Diffeomorphism { };  static const Diffeomorphism diffeomorphism;
	struct Immersion { };  static const Immersion immersion;
	struct BuildJacobian { };  static const BuildJacobian build_jacobian;
	struct ComposedWith { };  static const ComposedWith composed_with;
	struct OneDim { };  static const OneDim one_dim;
	struct HighDim { };  static const HighDim high_dim;
	struct Iff { };  static const Iff iff;
	struct PreviouslyNonExistent { };
	  static const PreviouslyNonExistent previously_non_existent;
	struct LessThan { };  static const LessThan less_than;
	struct IfLessThan { };  static const IfLessThan if_less_than;
	struct Otherwise { };  static const Otherwise otherwise;                  }
	

class Manifold;
	

class Function

// just a a thin wrapper around Function::Core
// I guess this is what they call "delegation" ...
// useful for arithmetic operators and for vector functions

{	public :

	class Core;

	static unsigned int total_cores;

	#ifndef NDEBUG
	static std::map < const Function::Core*, std::string > name;
	#endif

	Function::Core * core;

	inline Function ( const tag::NonExistent & ) : core ( nullptr ) { };

	inline Function ( const tag::WhoseCoreIs &, Function::Core * c );

	inline Function ( const Function & f );

	inline Function ( const tag::HasSize &, size_t s );

	inline Function ( double c );

	inline Function ( const tag::Diffeomorphism &, const tag::OneDim &,
	                  const Function & geom_coords, const Function & master_coords,
	                  const Function & geom_back_coords                             );

	inline Function ( const tag::Diffeomorphism &, const tag::HighDim &,
	                  const Function & geom_coords, const Function & master_coords,
	                  const Function & geom_back_coords                             );

	inline Function ( const tag::Diffeomorphism &, const Function & geom_coords,
                    const Function & master_coords, const Function & geom_back_coords );

	inline Function ( const tag::Immersion &, const Function & geom_coords,
                    const Function & master_coords, const Function & geom_back_coords );

	inline Function ( const Function & expr,
                    const tag::ComposedWith &, const Function & map );
	// an 'expr'ession involving master coordinates ( e.g.  1. - xi - eta )
	// composed with a map (diffeomorphism or immersion) sending it in the physical space

	inline ~Function();

	inline Function & operator= ( const Function & m );
	
	inline size_t nb_of_components ( ) const;

	inline Function operator[] ( size_t ) const;

	inline void conditionally_dispose_core ( );
	inline void set_core ( Function::Core *, const tag::PreviouslyNonExistent & );
	inline void set_core_to_null ( );
	inline void change_core_to ( Function::Core * );
	
	class TakenOnCell;
	inline Function::TakenOnCell operator() ( const Cell & cll ) const;

	inline Function deriv ( const Function & x ) const;  // derivative with respect to x

	inline Function replace ( const Function & x, const Function & y ) const;
	// in an expression, replace x by y

	#ifndef NDEBUG
	inline std::string repr ( ) const;
	enum From { from_void, from_sum, from_product, from_power, from_function };
	#endif

	class Scalar;  class ArithmeticExpression;  class Constant;
	class Vector;  class Aggregate;  class CoupledWithField;
	class Sum;  class Product;  class Power;  class Sqrt;  class Sin;  class Cos;  class Step;
	class Map;  class Diffeomorphism;  class Immersion;  class Composition;
	class Equality;

};  // end of  class Function


//-----------------------------------------------------------------------------------------//

class Function::Core

// just a base for classes like Constant, Sum, Product and many others
// we use dynamic_cast

{	public :

	// Manifold * manifold;

	unsigned int nb_of_wrappers { 0 };

	inline Core ( ) { Function::total_cores++; };
	
	virtual ~Core ( );

	inline Core ( const Function::Core & ) = delete;
	inline Core ( Function::Core && ) = delete;
	
	inline bool dispose ( )
	{	assert ( nb_of_wrappers > 0 );
		nb_of_wrappers--;
		return ( nb_of_wrappers == 0 );  }

	virtual size_t nb_of_components ( ) const = 0;

	virtual Function component ( size_t i ) = 0;

	virtual Function deriv ( Function ) const = 0;

	virtual Function replace ( const Function & x, const Function & y ) = 0;
	// in an expression, replace x by y

	#ifndef NDEBUG	
	virtual std::string repr ( const Function::From & from = Function::from_void ) const = 0;
	#endif
};

//-----------------------------------------------------------------------------------------//

inline Function::Function ( const tag::WhoseCoreIs &, Function::Core * c )
:	core ( c )
{	assert ( c );  c->nb_of_wrappers++;  }

inline Function::Function ( const Function & f )
:	Function ( tag::whose_core_is, f.core )
{	}
	
inline void Function::set_core ( Function::Core * c, const tag::PreviouslyNonExistent & )
{	assert ( this->core == nullptr );
	assert ( c );
	this->core = c;
	c->nb_of_wrappers++;               }
	
inline void Function::conditionally_dispose_core ( )
{	if ( this->core )
		if ( this->core->dispose() )
			delete this->core;          }
	
inline void Function::set_core_to_null ( )
{	this->conditionally_dispose_core ();
	this->core = nullptr;                 }
	
inline void Function::change_core_to ( Function::Core * c )
{	this->conditionally_dispose_core ();
	assert ( c );
	this->core = c;
	c->nb_of_wrappers++;                  }
	
inline Function & Function::operator= ( const Function & m )
{	this->change_core_to ( m.core );
	return *this;                     }
	
inline Function::~Function()
{	this->conditionally_dispose_core ();  }

inline Function Function::replace ( const Function & x, const Function & y ) const
{	return this->core->replace ( x, y );  }

//-----------------------------------------------------------------------------------------//

class Function::Scalar : public Function::Core
	
{	public :

	inline Scalar ( ) { };
	
	inline Scalar ( const Function::Scalar & ) = delete;
	inline Scalar ( Function::Scalar && ) = delete;

	size_t nb_of_components ( ) const;  // virtual from Function::Core, here returns 1

	Function component ( size_t i ); // virtual from Function::Core

	virtual double set_value_on_cell ( Cell::Core *, const double & ) = 0;
	// assign a numeric value to the function on the cell and return that value
	
	virtual double get_value_on_cell ( Cell::Core * ) const = 0;

	// Function deriv ( Function )
	// Function replace ( const Function & x, const Function & y )
	//    stay pure virtual from Function::Core
	
	#ifndef NDEBUG	
	// string repr ( const Function::From & from = Function::from_void )
	//   stays pure virtual from Function::Core
	#endif
};


//-----------------------------------------------------------------------------------------//

class Function::ArithmeticExpression : public Function::Scalar
	
// base for classes like Constant, Sum, Product
// any function for which 'set_value_on_cell' does not make sense

{	public :

	inline ArithmeticExpression ( ) { };

	inline ArithmeticExpression ( const Function::ArithmeticExpression & ) = delete;
	inline ArithmeticExpression ( Function::ArithmeticExpression && ) = delete;
	
	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1
	// Function component ( size_t i ) defined by Function::Scalar, returns self

	double set_value_on_cell ( Cell::Core *, const double & );  // virtual from Function::Scalar
	// here execution forbidden
	
	// double get_value_on_cell ( Cell::Core * )  stays pure virtual from Function::Scalar

	// Function deriv ( Function )
	// Function replace ( const Function & x, const Function & y )
	//   stay pure virtual from Function::Core, through Function::Scalar
	
	#ifndef NDEBUG	
	// string repr ( const Function::From & from = Function::from_void )
	//   stays pure virtual from Function::Core, through Function::Scalar
	#endif
};

//-----------------------------------------------------------------------------------------//

class Function::Constant : public Function::ArithmeticExpression
	
{	public :

	double val;

	inline Constant ( double c )
	:	val { c } { }

	inline Constant ( const Function::Constant & ) = delete;
	inline Constant ( Function::Constant && ) = delete;
	
	inline Function::Constant operator= ( const Function::Constant & ) = delete;
	inline Function::Constant operator= ( Function::Constant && ) = delete;

	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1
	// Function component ( size_t i ) defined by Function::Scalar, returns self

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   defined by Function::ArithmeticExpression (execution forbidden)

	double get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Scalar, through Function::ArithmeticExpression

	Function deriv ( Function ) const;
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	#endif
};

//-----------------------------------------------------------------------------------------//

inline Function::Function ( double c )
:	Function ( tag::whose_core_is, new Function::Constant ( c ) )
{	}

//-----------------------------------------------------------------------------------------//

class Function::Sum : public Function::ArithmeticExpression
	
{	public :

	std::forward_list < Function > terms;  // we use wrappers as a pointers

	inline Sum ( ) { };
	
	inline Sum ( const Function::Sum & ) = delete;
	inline Sum ( Function::Sum && ) = delete;
	
	inline Function::Sum operator= ( const Function::Sum & ) = delete;
	inline Function::Sum operator= ( Function::Sum && ) = delete;

	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1
	// Function component ( size_t i ) defined by Function::Scalar, returns self

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   defined by Function::ArithmeticExpression (execution forbidden)

	double get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Scalar, through Function::ArithmeticExpression

	Function deriv ( Function ) const;
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	#endif
};

//-----------------------------------------------------------------------------------------//

class Function::Product : public Function::ArithmeticExpression
	
{	public :

	std::forward_list < Function > factors;  // we use wrappers as a pointers

	inline Product ( ) { };

	inline Product ( const Function::Product & ) = delete;
	inline Product ( Function::Product && ) = delete;
	
	inline Function::Product operator= ( const Function::Product & ) = delete;
	inline Function::Product operator= ( Function::Product && ) = delete;

	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1
	// Function component ( size_t i ) defined by Function::Scalar, returns self

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   defined by Function::ArithmeticExpression (execution forbidden)

	double get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Scalar, through Function::ArithmeticExpression

	Function deriv ( Function ) const;
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	#endif
};

//-----------------------------------------------------------------------------------------//

class Function::Power : public Function::ArithmeticExpression
	
{	public :

	Function base;  // we use wrapper as a pointer
	double exponent;

	inline Power ( Function b, double e )
	:	base { b }, exponent { e }
	{	assert ( dynamic_cast < Function::Scalar* > ( b.core ) );  }

	inline Power ( const Function::Power & ) = delete;
	inline Power ( Function::Power && ) = delete;
	
	inline Function::Power operator= ( const Function::Power & ) = delete;
	inline Function::Power operator= ( Function::Power && ) = delete;

	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1
	// Function component ( size_t i ) defined by Function::Scalar, returns self

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   defined by Function::ArithmeticExpression (execution forbidden)

	double get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Scalar, through Function::ArithmeticExpression

	Function deriv ( Function ) const;
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	#endif
};

//-----------------------------------------------------------------------------------------//

class Function::Sqrt : public Function::ArithmeticExpression
	
{	public :

	Function base;  // we use wrapper as a pointer

	inline Sqrt ( const Function & b ) : base { b }
	{	assert ( dynamic_cast < Function::Scalar* > ( b.core ) );  }

	inline Sqrt ( const Function::Sqrt & ) = delete;
	inline Sqrt ( Function::Sqrt && ) = delete;
	
	inline Function::Sqrt operator= ( const Function::Sqrt & ) = delete;
	inline Function::Sqrt operator= ( Function::Sqrt && ) = delete;
	
	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1
	// Function component ( size_t i ) defined by Function::Scalar, returns self

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   defined by Function::ArithmeticExpression (execution forbidden)

	double get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Scalar, through Function::ArithmeticExpression

	Function deriv ( Function ) const;
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	#endif
};

//-----------------------------------------------------------------------------------------//

inline Function sqrt ( const Function & f )
{	return Function ( tag::whose_core_is, new Function::Sqrt ( f ) );  }

//-----------------------------------------------------------------------------------------//

class Function::Sin : public Function::ArithmeticExpression
	
{	public :

	Function base;  // we use wrapper as a pointer

	inline Sin ( const Function & b ) : base { b }
	{	assert ( dynamic_cast < Function::Scalar* > ( b.core ) );  }

	inline Sin ( const Function::Sin & ) = delete;
	inline Sin ( Function::Sin && ) = delete;
	
	inline Function::Sin operator= ( const Function::Sin & ) = delete;
	inline Function::Sin operator= ( Function::Sin && ) = delete;
	
	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1
	// Function component ( size_t i ) defined by Function::Scalar, returns self

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   defined by Function::ArithmeticExpression (execution forbidden)

	double get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Scalar, through Function::ArithmeticExpression

	Function deriv ( Function ) const;
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	#endif
};

//-----------------------------------------------------------------------------------------//

inline Function sin ( const Function & f )
{ return Function ( tag::whose_core_is, new Function::Sin ( f ) );  }

//-----------------------------------------------------------------------------------------//

class Function::Cos : public Function::ArithmeticExpression
	
{	public :

	Function base;  // we use wrapper as a pointer

	inline Cos ( const Function & b ) : base { b }
	{	assert ( dynamic_cast < Function::Scalar* > ( b.core ) );  }

	inline Cos ( const Function::Cos & ) = delete;
	inline Cos ( Function::Cos && ) = delete;
	
	inline Function::Cos operator= ( const Function::Cos & ) = delete;
	inline Function::Cos operator= ( Function::Cos && ) = delete;

	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1
	// Function component ( size_t i ) defined by Function::Scalar, returns self

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   defined by Function::ArithmeticExpression (execution forbidden)

	double get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Scalar, through Function::ArithmeticExpression

	Function deriv ( Function ) const;
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	#endif
};

//-----------------------------------------------------------------------------------------//

inline Function cos ( const Function & f )
{	return Function ( tag::whose_core_is, new Function::Cos ( f ) );  }

//-----------------------------------------------------------------------------------------//

class Function::Step : public Function::ArithmeticExpression

// piecewise constant functions 

{	public :

	Function arg;  // we use wrappers as pointers
	// 'arg' is used for testing inequalities

	std::vector < Function > values;
	std::vector < double > cuts;  //  values.size() == cuts.size() + 1
	// 'cuts' contains values ordered increasingly
	// for instance, suppose arg==x, values=={x,1/x,3} and cuts=={-1,1}
	// this means this(x)=x for x<-1, this(x)=1/x for -1<x<1, this(x)=3 for x>1

	inline Step ( const Function & v1, const tag::Iff, const Function & x,
                const tag::LessThan &, double c,
	              const Function & v2, const tag::Otherwise &             )
	:	arg { x }, values { v1, v2 }, cuts { c }
	{	assert ( dynamic_cast < Function::Scalar* > ( x.core ) );
		assert ( dynamic_cast < Function::Scalar* > ( v1.core ) );
		assert ( dynamic_cast < Function::Scalar* > ( v2.core ) );  }
		
	inline Step ( const Function & v1, const tag::Iff, const Function & x,
                const tag::LessThan &, double c1,
                const Function & v2, const tag::IfLessThan, double c2,
	              const Function & v3, const tag::Otherwise &             )
	:	arg { x }, values { v1, v2, v3 }, cuts { c1, c2 }
	{	assert ( c1 < c2 );
		assert ( dynamic_cast < Function::Scalar* > ( x.core ) );
		assert ( dynamic_cast < Function::Scalar* > ( v1.core ) );
		assert ( dynamic_cast < Function::Scalar* > ( v2.core ) );
		assert ( dynamic_cast < Function::Scalar* > ( v3.core ) );  }
	
	inline Step ( const Function & v1, const tag::Iff, const Function & x,
                const tag::LessThan &, double c1,
                const Function & v2, const tag::IfLessThan, double c2,
                const Function & v3, const tag::IfLessThan, double c3,
	              const Function & v4, const tag::Otherwise &             )
	:	arg { x }, values { v1, v2, v3, v4 }, cuts { c1, c2, c3 }
	{	assert ( c1 < c2 );  assert ( c2 < c3 );
		assert ( dynamic_cast < Function::Scalar* > ( x.core ) );
		assert ( dynamic_cast < Function::Scalar* > ( v1.core ) );
		assert ( dynamic_cast < Function::Scalar* > ( v2.core ) );
		assert ( dynamic_cast < Function::Scalar* > ( v3.core ) );
		assert ( dynamic_cast < Function::Scalar* > ( v4.core ) );  }
	
	inline Step ( const Function & v1, const tag::Iff, const Function & x,
                const tag::LessThan &, double c1,
                const Function & v2, const tag::IfLessThan, double c2,
                const Function & v3, const tag::IfLessThan, double c3,
                const Function & v4, const tag::IfLessThan, double c4,
	              const Function & v5, const tag::Otherwise &             )
	:	arg { x }, values { v1, v2, v3, v4, v5 }, cuts { c1, c2, c3, c4 }
	{	assert ( c1 < c2 );  assert ( c1 < c2 );  assert ( c3 < c4 );
		assert ( dynamic_cast < Function::Scalar* > ( x.core ) );
		assert ( dynamic_cast < Function::Scalar* > ( v1.core ) );
		assert ( dynamic_cast < Function::Scalar* > ( v2.core ) );
		assert ( dynamic_cast < Function::Scalar* > ( v3.core ) );
		assert ( dynamic_cast < Function::Scalar* > ( v4.core ) );
		assert ( dynamic_cast < Function::Scalar* > ( v5.core ) );     }

	inline Step ( const Function & aarg, std::vector < Function > & vals,
                std::vector < double > cts )
	: arg { aarg }, values { vals }, cuts { cts }
	{	assert ( vals.size() == cts.size() + 1 );  }  // add asserts !!

	inline Step ( const Function::Step & ) = delete;
	inline Step ( Function::Step && ) = delete;
	
	inline Function::Step operator= ( const Function::Step & ) = delete;
	inline Function::Step operator= ( Function::Step && ) = delete;

	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1
	// Function component ( size_t i ) defined by Function::Scalar, returns self

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   defined by Function::ArithmeticExpression (execution forbidden)

	double get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Scalar, through Function::ArithmeticExpression

	Function deriv ( Function ) const;
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	#endif

}; // end of  class Function::Step

//-----------------------------------------------------------------------------------------//

class Function::Vector : public Function::Core
	
{	public :

	inline Vector ( ) { };

	inline Vector ( const Function::Vector & ) = delete;
	inline Vector ( Function::Vector && ) = delete;
	
	// size_t nb_of_components ( )  stays pure virtual from Function::Core
	
	// Function component ( )  stays pure virtual from Function::Core

	virtual std::vector<double> set_value_on_cell
	( Cell::Core *, const std::vector<double> & ) = 0;
	// assign a numeric vector to the function on the cell and return that vector
	
	virtual std::vector<double> get_value_on_cell ( Cell::Core * ) const = 0;

	Function deriv ( Function ) const;
	//  virtual from Function::Core, here execution forbidden, to change

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, here execution forbidden
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	//  virtual from Function::Core, here forbids execution
	#endif
};

//-----------------------------------------------------------------------------------------//

class Function::Aggregate : public Function::Vector
	
{	public :

	std::vector < Function > components;
	
	inline Aggregate ( const tag::ReserveSize &, size_t s )
	:	Function::Vector()
	{	components.reserve(s);  }

	inline Aggregate ( const Function::Aggregate & ) = delete;
	inline Aggregate ( Function::Aggregate && ) = delete;
	
	inline Function::Aggregate operator= ( const Function::Aggregate & ) = delete;
	inline Function::Aggregate operator= ( Function::Aggregate && ) = delete;

	// 'nb_of_components' and 'component' are virtual from Function::Core, through Function::Vector
	// overridden by Function::CoupledWithField::Vector
	size_t nb_of_components ( ) const;
	Function component ( size_t i );
	
	std::vector<double> set_value_on_cell ( Cell::Core *, const std::vector<double> & );
	// virtual from Function::Vector

	std::vector<double> get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Vector

	// Function deriv ( Function )
	//    defined by Function::Vector (execution forbidden), to change

	// Function replace ( const Function & x, const Function & y );
	//    defined by Function::Vector (execution forbidden), to change
	
	#ifndef NDEBUG	
	// std::string repr ( const Function::From & from = Function::from_void )
	//    defined by Function::Vector (execution forbidden)
	#endif
};


//-----------------------------------------------------------------------------------------//

inline Function operator&& ( Function f, Function g )

{	size_t nf = f.nb_of_components(), ng = g.nb_of_components();
	Function::Aggregate * res = new Function::Aggregate ( tag::reserve_size, nf + ng ) ;
	for ( size_t i = 0; i < nf; i++ ) res->components.emplace_back ( f[i] );
	for ( size_t i = 0; i < ng; i++ ) res->components.emplace_back ( g[i] );
	return Function ( tag::whose_core_is,	res );                                         }

//-----------------------------------------------------------------------------------------//

inline bool operator< ( const Function & f, const Function & g )
{	return f.core < g.core;  }
// needed for map 'jacobian' in class Function::Diffeomorphism::HighDim

//-----------------------------------------------------------------------------------------//


class Function::Map

// used for Function::Diffeomorphism::OneDim and Function::Immersion
// Functions inheriting from this class may be used for Function::Composition

{	public :

	Function geom_coords, master_coords, back_geom_coords;
	// back_geom_coords are expressions involving the master coordinates
	// mathematically, they correspond to the same function geom_coords

	// for immersions, 'jacobian' contains the derivatives of
	// back_geom_coords with respect to master_coords

	// for Function::Diffeomorphism::***Dim, 'jacobian' contains these derivatives inverted
	// that is, we keep the inverse matrix, which can be seen as the
	// derivatives of the master_coords with respect to geom_coords
	// they should be composed with 'this' map in the calling code

	std::map < Function, Function > jacobian;
	
	Function det;  // positive scalar, dilation coefficient

	inline Map ( const Function & gc, const Function & mc, const Function & bgc )
	:	geom_coords ( gc ), master_coords ( mc ), back_geom_coords ( bgc ),
		jacobian ( ), det( tag::non_existent )
	{	}
	
};  // end of Function::Map


class Function::Diffeomorphism

// useful only for asserting dynamic_cast

{	public :
	class OneDim; class HighDim; };

//-----------------------------------------------------------------------------------------//


class Function::Diffeomorphism::OneDim
: public Function::Scalar, public Function::Map, public Function::Diffeomorphism

// there are one-dimensional diffeomorphisms and we must treat them separately
// because they inherit from Function::Scalar, not from Function::Vector

{	public :

	// Function geom_coords, master_coords, back_geom_coords
	// inherited from Function::Map
	// back_geom_coords are expressions involving the master coordinates
	// mathematically, they correspond to the same function geom_coords

	// std::map < Function, Function > jacobian  inherited from Function::Map
	// contains one over the derivative of back_geom_coords
	//          with respect to master_coords  ( 1./det )
	
	// 'jacobian' is an expression involving master_cords
	// if we want to  look at them as functions of 'geom_coords',
	// in particular, if we want derivatives with respect to 'geom_coords',
	// they should be composed with 'this' map in the calling code

	// Function det  inherited from Function::Map
  // positive scalar, dilation coefficient
	// can be used for integration
	// here, 'det' is equal to the inverse of the only component of 'jacobian'
	// thus, 'det' is	the derivative of back_geom_coords with respect to master_coords
	
	OneDim ( const Function & gc, const Function & mc, const Function & bgc,
           const tag::BuildJacobian &                                      );
	
	inline OneDim ( const Function::Diffeomorphism::OneDim & ) = delete;
	inline OneDim ( Function::Diffeomorphism::OneDim && ) = delete;
	
	inline Function::Diffeomorphism::OneDim operator=
		( const Function::Diffeomorphism::OneDim & ) = delete;
	inline Function::Diffeomorphism::OneDim operator=
		( Function::Diffeomorphism::OneDim && ) = delete;
	
	// Function component ( size_t i ) defined by Function::Scalar, returns self
	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1

	double set_value_on_cell ( Cell::Core *, const double & );
	// virtual from Function::Vector

	double get_value_on_cell ( Cell::Core * ) const;  // virtual from Function::Scalar

	Function deriv ( Function ) const;
	//  virtual from Function::Core, through Function::Scalar

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, through Function::Scalar
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar
	#endif

};  // end of class Function::Diffeomorphism::OneDim
	
//-----------------------------------------------------------------------------------------//


inline Function::Function ( const tag::Diffeomorphism &, const tag::OneDim &,
                            const Function & geom_coords, const Function & master_coords,
                            const Function & back_geom_coords                             )
: Function ( tag::non_existent )
{	assert ( geom_coords.nb_of_components() == back_geom_coords.nb_of_components() );
	assert ( master_coords.nb_of_components() == geom_coords.nb_of_components() );
	assert ( geom_coords.nb_of_components() == 1 );
	this->set_core ( new Function::Diffeomorphism::OneDim
	   ( geom_coords, master_coords, back_geom_coords, tag::build_jacobian ),
	   tag::previously_non_existent                                           );       }
	

//-----------------------------------------------------------------------------------------//

class Function::Immersion
: public Function::Vector, public Function::Map

// a map from a master manifold to a geometric manifold
// geometric dimension striclty higher than master dimension

{	public :

	// Function geom_coords, master_coords, back_geom_coords
	// inherited from Function::Map
	// back_geom_coords are expressions involving the master coordinates
	// mathematically, they correspond to the same function geom_coords

	// std::map < Function, Function > jacobian  inherited from Function::Map
	// contains the derivatives of back_geom_coords with respect to master_coords
	// quite different from Function::Diffeomorphism::***Dim !

	// Function det  inherited from Function::Map
  // positive scalar, dilation coefficient, can be used for integration

	// if master_coords is 1D then
	// det = |D Phi| (Phi being back_geom_coords)
	// D Phi is a vector in [back_]geom_coords manifold
	// (more precisely, it is tangent to the manifold)

	// if master_coords is 2D
	// let (e,f) be an orthonormal basis in master_coords
	// let v = D Phi e, w = D Phi f  (Phi being back_geom_coords)
	// v and w belong to [back_]geom_coords manifold
	// (more precisely, they are tangent to the manifold)
	// then this dilation coef 'det' is square root of  v.v w.w - v.w v.w
	// that is, square root of  |v|^2 |w|^2 - (v.w)^2
	// that is, square root of  vi vi wj wj - vi vj wi wj
	// is there any way to avoid computing a square root ? probably no

	// components of 'jacobian', as well as 'det', are expressions involving master_cords

	inline Immersion ( const Function & gc, const Function & mc, const Function & bgc )
	:	Function::Vector ( ), Function::Map ( gc, mc, bgc )
	{	}  // 'jacobian' initialized as empty map, 'det' as non_existent

	Immersion ( const Function & gc, const Function & mc, const Function & bgc,
              const tag::BuildJacobian &                                      );

	inline Immersion ( const Function::Immersion & ) = delete;
	inline Immersion ( Function::Immersion && ) = delete;
	
	inline Function::Immersion operator= ( const Function::Immersion & ) = delete;
	inline Function::Immersion operator= ( Function::Immersion && ) = delete;

	size_t nb_of_components ( ) const;  // virtual from Function::Core, through Function::Vector

	Function component ( size_t i );  // virtual from Function::Core, through Function::Vector
	
	std::vector<double> set_value_on_cell ( Cell::Core *, const std::vector<double> & );
	// virtual from Function::Vector

	std::vector<double> get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Vector

	// Function deriv ( Function )
	//    defined by Function::Vector (execution forbidden), to change

	// Function replace ( const Function & x, const Function & y );
	//    defined by Function::Vector (execution forbidden), to change
	
	#ifndef NDEBUG	
	// std::string repr ( const Function::From & from = Function::from_void )
	//    defined by Function::Vector (execution forbidden)
	#endif

};  // end of class Function::Immersion

//-----------------------------------------------------------------------------------------//

inline Function::Function ( const tag::Immersion &, const Function & geom_coords,
                            const Function & master_coords, const Function & back_geom_coords )
:	Function ( tag::whose_core_is, new Function::Immersion
						 ( geom_coords, master_coords, back_geom_coords, tag::build_jacobian ) )
{	}

//-----------------------------------------------------------------------------------------//


class Function::Diffeomorphism::HighDim
: public Function::Immersion, public Function::Diffeomorphism

// the inheritance from Function::Immersion is quite misleading !
// here, geometric dimension is equal master dimension
// the interpretation of attribute 'jacobian' is different from Function::Immersion
	
{	public :

	// Function geom_coords, master_coords, back_geom_coords
	// inherited from Function::Map
	// back_geom_coords are expressions involving the master coordinates
	// mathematically, they correspond to the same function geom_coords

	// std::map < Function, Function > jacobian  inherited from Function::Map

	// here, the meaning of the 'jacobian' attribute is different from Function::Immersion
	// we keep here the matrix already inverted
	// conceptually, this 'jacobian' contains the derivatives of master_coords
	// with respect to geom_coords
	// this is implemented backwards, by computing the derivatives of
	// back_geom_coords with respect to master_coords
	// then taking the inverse matrix

	// Function det  inherited from Function::Immersion
  // determinant of the jacobian matrix - not inverted !
	// should be positive, so may be used directly for integration

	// components of 'jacobian', as well as 'det', are expressions involving master_cords
	// if we want to  look at them as functions of 'geom_coords',
	// in particular if we want derivatives with respect to 'geom_coords',
	// they should be composed with 'this' map in the calling code

	HighDim ( const Function & gc, const Function & mc, const Function & bgc,
                   const tag::BuildJacobian &                               );

	inline HighDim ( const Function::Diffeomorphism::HighDim & ) = delete;
	inline HighDim ( Function::Diffeomorphism::HighDim && ) = delete;
	
	inline Function::Diffeomorphism::HighDim operator=
		( const Function::Diffeomorphism::HighDim & ) = delete;
	inline Function::Diffeomorphism::HighDim operator=
		( Function::Diffeomorphism::HighDim && ) = delete;

	// size_t nb_of_components ( ) const   and
	// Function component ( size_t i )    defined in Function::Immersion, execution forbidden
	
	// the two methods below defined by Function::ArithmeticExpression (execution forbidden)
	// std::vector<double> set_value_on_cell ( Cell::Core *, const std::vector<double> & )
	// std::vector<double> get_value_on_cell ( Cell::Core * ) const

	// Function deriv ( Function )
	//    defined by Function::Vector (execution forbidden)

	// Function replace ( const Function & x, const Function & y );
	//    defined by Function::Vector (execution forbidden)
	
	#ifndef NDEBUG	
	// std::string repr ( const Function::From & from = Function::from_void )
	//    defined by Function::Vector (execution forbidden)
	#endif

};  // end of class Function::Diffeomorphism::HighDim

//-----------------------------------------------------------------------------------------//


inline Function::Function ( const tag::Diffeomorphism &, const tag::HighDim &,
                            const Function & geom_coords, const Function & master_coords,
                            const Function & back_geom_coords                             )
: Function ( tag::non_existent )
{	assert ( geom_coords.nb_of_components() == back_geom_coords.nb_of_components() );
	assert ( master_coords.nb_of_components() == geom_coords.nb_of_components() );
	assert ( geom_coords.nb_of_components() >= 2 );
	this->set_core ( new Function::Diffeomorphism::HighDim
	   ( geom_coords, master_coords, back_geom_coords, tag::build_jacobian ),
	   tag::previously_non_existent                                           );       }

	
inline Function::Function ( const tag::Diffeomorphism &, const Function & geom_coords,
                            const Function & master_coords, const Function & back_geom_coords )
: Function ( tag::non_existent )
{	assert ( geom_coords.nb_of_components() == back_geom_coords.nb_of_components() );
	assert ( master_coords.nb_of_components() == geom_coords.nb_of_components() );
	if ( geom_coords.nb_of_components() == 1 )
	{	this->set_core ( new Function::Diffeomorphism::OneDim
		   ( geom_coords, master_coords, back_geom_coords, tag::build_jacobian ),
		   tag::previously_non_existent                                           );
		return;                                                                       }
	this->set_core ( new Function::Diffeomorphism::HighDim
	   ( geom_coords, master_coords, back_geom_coords, tag::build_jacobian ),
	   tag::previously_non_existent                                           );       }
	

//-----------------------------------------------------------------------------------------//

class Function::Composition : public Function::Scalar
	
// an expression 'base' involving master coordinates ( e.g.  1. - xi - eta )
// composed with an 'transf'ormation sending it in the physical space
// note that 'transf' could be an Immersion or a Diffeomorphism::OneDim or ::HighDim

{	public :

	Function base, transf;
	
	inline Composition ( const Function & b, const Function & tr )
	:	Function::Scalar(), base (b), transf (tr)
	{ }

	inline Composition ( const Function::Composition & ) = delete;
	inline Composition ( Function::Composition && ) = delete;
	
	inline Function::Composition operator= ( const Function::Composition & ) = delete;
	inline Function::Composition operator= ( Function::Composition && ) = delete;

	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1
	// Function component ( size_t i ) defined by Function::Scalar, returns self

	double set_value_on_cell ( Cell::Core *, const double & );  // virtual from Function::Scalar
	double get_value_on_cell ( Cell::Core * ) const;  // virtual from Function::Scalar

	Function deriv ( Function ) const;
	//  virtual from Function::Core, through Function::Scalar

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, through Function::Scalar
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar
	#endif
};

//-----------------------------------------------------------------------------------------//


inline Function::Function
(	const Function & expr, const tag::ComposedWith &, const Function & f_map )

// an 'expr'ession involving master coordinates ( e.g.  1.- xi - eta )
// composed with a map 'f_map' sending it in the physical space

: Function ( tag::non_existent )

{	assert ( dynamic_cast < Function::Map* > ( f_map.core ) );
  Function::Constant * expr_c = dynamic_cast < Function::Constant* > ( expr.core );
	if ( expr_c ) this->set_core ( expr.core, tag::previously_non_existent );
	else this->set_core ( new Function::Composition ( expr, f_map ),
                        tag::previously_non_existent               );                }

//-----------------------------------------------------------------------------------------//


class Function::CoupledWithField
	
{	public :

  Field::Core * field;

	inline CoupledWithField ( Field::Core * f ) : field { f }
	{	assert ( f );  }
		
	class Scalar;  class Vector;
};

//-----------------------------------------------------------------------------------------//

class Function::CoupledWithField::Scalar
: public Function::Scalar,
	public Function::CoupledWithField
	
{	public :

	inline Scalar ( Field::Scalar * f )
	:	Function::Scalar(),
		Function::CoupledWithField ( f )
	{	assert ( f->nb_of_components() == 1 );  }
		
	inline Scalar ( const Function::CoupledWithField::Scalar & ) = delete;
	inline Scalar ( Function::CoupledWithField::Scalar && ) = delete;
	
	inline Function::CoupledWithField::Scalar operator=
		( const Function::CoupledWithField::Scalar & ) = delete;
	inline Function::CoupledWithField::Scalar operator=
		( Function::CoupledWithField::Scalar && ) = delete;

	// Function component ( size_t i ) defined by Function::Scalar, returns self
	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1

	double set_value_on_cell ( Cell::Core *, const double & );  // virtual from Function::Scalar
	double get_value_on_cell ( Cell::Core * ) const;  // virtual from Function::Scalar

	Function deriv ( Function ) const;
	//  virtual from Function::Core, through Function::Scalar

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, through Function::Scalar
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar
	#endif
};

//-----------------------------------------------------------------------------------------//

class Function::CoupledWithField::Vector
: public Function::Aggregate,
	public Function::CoupledWithField

// inheriting from Function::Aggregate means simply that there is a 'components' member

// quando calculamos as componentes duma funcao associada a um campo vectorial,
// nao precisamos de criar as componentes do campo vectorial
// podemos ter uma funcao escalar associada a um indice dentro dum campo vectorial
	
{	public :

	inline Vector ( Field::Block * f )
	:	Function::Aggregate ( tag::reserve_size, f->nb_of_components() ),
		Function::CoupledWithField ( f )
	{	size_t n = f->nb_of_components();
		for ( size_t j = 0; j < n; j++ )
//		components[j] = Function ( tag::whose_core_is,
			components.emplace_back ( tag::whose_core_is,
				new Function::CoupledWithField::Scalar ( field->component(j) ) );  }
	
	inline Vector ( const Function::CoupledWithField::Vector & ) = delete;
	inline Vector ( Function::CoupledWithField::Vector && ) = delete;
	
	inline Function::CoupledWithField::Vector operator=
		( const Function::CoupledWithField::Vector & ) = delete;
	inline Function::CoupledWithField::Vector operator=
		( Function::CoupledWithField::Vector && ) = delete;

	virtual size_t nb_of_components ( ) const override;
	// virtual from Function::Core through Function::Vector, Function::Aggregate
	
	virtual Function component ( size_t i ) override;
	// virtual from Function::Core, through Function::Vector, Function::Aggregate
	
	std::vector<double> set_value_on_cell ( Cell::Core *, const std::vector<double> & );
	// virtual from Function::Vector, through Function::Aggregate
	
	std::vector<double> get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Vector, through Function::Aggregate

	// Function deriv ( Function )
	//    defined by Function::Vector (execution forbidden), may change

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, through Function::Vector
	
	#ifndef NDEBUG	
	// std::string repr ( const Function::From & from = Function::from_void );
	//    defined by Function::Vector (execution forbidden)
	#endif

};  // end of class Function::CoupledWithField::Vector

//-----------------------------------------------------------------------------------------//


#ifndef NDEBUG
inline std::string Function::repr ( ) const  { return core->repr();  }
#endif

//-----------------------------------------------------------------------------------------//

Function operator+ ( const Function & f, const Function & g );
Function operator* ( const Function & f, const Function & g );
Function power ( const Function & f, double e );

inline Function operator- ( const Function & f )
{	Function minus_one ( -1. );
	return minus_one * f;         }

inline Function operator- ( const Function & f, const Function & g )
{	return f + (-g);  }

inline Function operator+= ( Function & f, const Function & g )
{	return f = f + g;  }

inline Function operator*= ( Function & f, const Function & g )
{	return f = f * g;  }

inline Function operator/ ( const Function & f, const Function & g )
{	return f * power ( g, -1. );  }

inline Function operator/= ( Function & f, const Function & g )
{	return f = f / g;  }

inline Function Function::deriv ( const Function & x ) const  // derivative with respect to x
{	return this->core->deriv ( x );  }

//-----------------------------------------------------------------------------------------//


namespace tag  {  struct Threshold { };  static const Threshold threshold;  }

inline Function sign_of ( const Function & f )
{	return Function ( tag::whose_core_is, new Function::Step
      ( Function (-1.), tag::iff, f, tag::less_than, 0., Function(1.), tag::otherwise ) );  }

inline Function smooth_abs ( const Function & f, const tag::Threshold &, double d )
// a smooth (C1) approximation of the absolute value
// 'd' is a threshold; if |f| >= d then smooth_abs(f) == |f|
// below 'd', the function is smooth but far from the absolute value (never less than d/2)
{	assert ( d > 0. );
	return Function ( tag::whose_core_is, new Function::Step
	      ( -f, tag::iff, f, tag::less_than, -d,
          ( f*f + d*d ) / (2.*d), tag::if_less_than, d, f, tag::otherwise ) );  }
//	Function abs_f = sign_of ( f ) * f;
//	return abs_f * abs_f / ( d + abs_f );  }

inline Function smooth_min
( const Function & f, const Function & g, const tag::Threshold &, double d )
// a smooth (C1) approximation of the minimum between two functions
{	return 0.5 * ( f + g - smooth_abs ( f-g, tag::threshold, d ) );  }

inline Function smooth_min
( const Function & f, const Function & g, const Function & h, const tag::Threshold &, double d )
// a smooth (C1) approximation of the minimum between three functions
{	return smooth_min ( f, smooth_min ( g, h, tag::threshold, d ), tag::threshold, d );  }

inline Function smooth_min
( const Function & f, const Function & g, const Function & h, const Function & i,
  const tag::Threshold &, double d                                                )
// a smooth (C1) approximation of the minimum between four functions
{	return smooth_min ( smooth_min ( f, g, tag::threshold, d ),
                      smooth_min ( h, i, tag::threshold, d ), tag::threshold, d );  }

inline Function smooth_min
( const Function & f, const Function & g, const Function & h, const Function & i,
  const Function & j ,const tag::Threshold &, double d                            )
// a smooth (C1) approximation of the minimum between five functions
{	return smooth_min ( smooth_min ( f, g, h, tag::threshold, d ),
	                    smooth_min ( i, j, tag::threshold, d ), tag::threshold, d );  }

inline Function smooth_min
( const Function & f, const Function & g, const Function & h, const Function & i,
  const Function & j, const Function & k, const tag::Threshold &, double d        )
// a smooth (C1) approximation of the minimum between six functions
{	return smooth_min ( smooth_min ( f, g, tag::threshold, d ),
											smooth_min ( h, i, tag::threshold, d ),
	                    smooth_min ( j, k, tag::threshold, d ), tag::threshold, d );  }

inline Function smooth_max
( const Function & f, const Function & g, const tag::Threshold &, double d )
// a smooth (C1) approximation of the maximum between two functions
{	return 0.5 * ( f + g + smooth_abs ( f-g, tag::threshold, d ) );  }

inline Function smooth_max
( const Function & f, const Function & g, const Function & h, const tag::Threshold &, double d )
// a smooth (C1) approximation of the maximum between three functions
{	return smooth_max ( f, smooth_max ( g, h, tag::threshold, d ), tag::threshold, d );  }

inline Function smooth_max
( const Function & f, const Function & g, const Function & h, const Function & i,
  const tag::Threshold &, double d                                                )
// a smooth (C1) approximation of the maximum between four functions
{	return smooth_max ( smooth_max ( f, g, tag::threshold, d ),
	                    smooth_max ( h, i, tag::threshold, d ), tag::threshold, d );  }

inline Function smooth_max
( const Function & f, const Function & g, const Function & h, const Function & i,
  const Function & j ,const tag::Threshold &, double d                            )
// a smooth (C1) approximation of the maximum between five functions
{	return smooth_max ( smooth_max ( f, g, h, tag::threshold, d ),
	                    smooth_max ( i, j, tag::threshold, d ), tag::threshold, d );  }

inline Function smooth_max
( const Function & f, const Function & g, const Function & h, const Function & i,
  const Function & j, const Function & k, const tag::Threshold &, double d        )
// a smooth (C1) approximation of the maximum between six functions
{	return smooth_max ( smooth_max ( f, g, tag::threshold, d ),
	                    smooth_max ( h, i, tag::threshold, d ),
	                    smooth_max ( j, k, tag::threshold, d ), tag::threshold, d );  }
										
//-----------------------------------------------------------------------------------------//


class Function::TakenOnCell	

{	public :

  Function::Core * f;
	// Function::TakenOnCell should only be used as temporary objects
	// they should be immediately converted to a (reference to a) double or vector<double>
	// so an ordinary pointer is OK here
	
	Cell::Core * cll;

	TakenOnCell ( const Function::TakenOnCell & other ) = delete;
	TakenOnCell ( const Function::TakenOnCell && other ) = delete;
	
	Function::TakenOnCell & operator= ( const Function::TakenOnCell & other )
	{	if ( f->nb_of_components() == 1 )
		{	(*this) = static_cast <double> (other);
			return *this;                            }
		else
		{	(*this) = static_cast <std::vector<double>> (other);
			return *this;                                          }  }

	Function::TakenOnCell & operator= ( const Function::TakenOnCell && other )
	{	if ( f->nb_of_components() == 1 )
		{	(*this) = static_cast <double> (other);
			return *this;                            }
		else
		{	(*this) = static_cast <std::vector<double>> (other);
			return *this;                                          }  }
	
	inline operator double() const
	// can be used like in  double x = f(cll)  or  cout << f(cll)
	{	Function::Scalar * f_scalar =
			Mesh::assert_cast < Function::Core*, Function::Scalar* > ( f );
		return f_scalar->get_value_on_cell ( cll );                       }
		
	inline double operator= ( const double & x )
	// can be used like in  f(cll) = 2.0
	{	Function::Scalar * f_scalar =
			Mesh::assert_cast < Function::Core*, Function::Scalar* > ( f );
		return f_scalar->set_value_on_cell ( cll, x );                  }

	inline double operator+= ( const double & x )
	// can be used like in  f(cll) += 2.0
	{	Function::Scalar * f_scalar =
			Mesh::assert_cast < Function::Core*, Function::Scalar* > ( f );
		return f_scalar->set_value_on_cell ( cll, f_scalar->get_value_on_cell(cll) + x ); }

	inline double operator-= ( const double & x )
	// can be used like in  f(cll) -= 2.0
	{	Function::Scalar * f_scalar =
			Mesh::assert_cast < Function::Core*, Function::Scalar* > ( f );
		return f_scalar->set_value_on_cell ( cll, f_scalar->get_value_on_cell(cll) - x ); }

	inline double operator*= ( const double & x )
	// can be used like in  f(cll) *= 2.0
	{	Function::Scalar * f_scalar =
			Mesh::assert_cast < Function::Core*, Function::Scalar* > ( f );
		return f_scalar->set_value_on_cell ( cll, f_scalar->get_value_on_cell(cll) * x ); }

	inline double operator/= ( const double & x )
	// can be used like in  f(cll) /= 2.0
	{	Function::Scalar * f_scalar =
			Mesh::assert_cast < Function::Core*, Function::Scalar* > ( f );
		return f_scalar->set_value_on_cell ( cll, f_scalar->get_value_on_cell(cll) / x ); }

	inline operator std::vector<double>() const
	// can be used like in  vector<double> vec = f(cll)
	{	Function::Vector * f_vect =
			Mesh::assert_cast < Function::Core*, Function::Vector* > ( f );
	  return f_vect->get_value_on_cell ( cll );                         }
		
	inline std::vector<double> operator= ( const std::vector<double> & x )
	// can be used like in  f(cll) = vec
	{	Function::Vector * f_vect =
			Mesh::assert_cast < Function::Core*, Function::Vector* > ( f );
		return f_vect->set_value_on_cell ( cll, x );                      }

};  // end of class Function::TakenOnCell	

//-----------------------------------------------------------------------------------------//

class Function::Equality

// the result of an equality comparison between Functions

{	public :

	Function lhs, rhs;
	
};

//-----------------------------------------------------------------------------------------//

inline Function::Equality operator== ( const Function & f, const Function & g )
{	return Function::Equality { f, g };  }


inline size_t Function::nb_of_components ( ) const
{	return this->core->nb_of_components();  }


inline Function Function::operator[] ( size_t i ) const
{	assert ( i < this->nb_of_components() );
	Function::Scalar * f_scalar = dynamic_cast < Function::Scalar * > ( this->core );
	if ( f_scalar )
	{	assert ( this->nb_of_components() == 1 );
		return * this;                             }
  Function::Vector * f_vector = dynamic_cast < Function::Vector * > ( this->core );
	assert ( f_vector );
	return Function ( tag::whose_core_is, f_vector->component(i).core );               }


inline Function::TakenOnCell Function::operator() ( const Cell & cll ) const
{	return Function::TakenOnCell { this->core, cll.core };   }


}  // namespace maniFEM


#endif
// ifndef MANIFEM_FUNCTION_H
