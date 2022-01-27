
// function.h 2022.01.26

//   This file is part of maniFEM, a C++ library for meshes and finite elements on manifolds.

//   Copyright 2019, 2020, 2021, 2022 Cristian Barbarosie cristian.barbarosie@gmail.com

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

#ifndef MANIFEM_FUNCTION_H
#define MANIFEM_FUNCTION_H

#include <iostream>
#include <vector>
#include <forward_list>
#include <memory>
#include "math.h"
#include "assert.h"

#include "mesh.h"
#include "iterator.h"
#include "field.h"


namespace maniFEM {


namespace tag
{	struct Diffeomorphism { };  static const Diffeomorphism diffeomorphism;
	struct Immersion { };  static const Immersion immersion;
	struct BuildJacobian { };  static const BuildJacobian build_jacobian;
	struct ComposedWith { };  static const ComposedWith composed_with;
	struct Iff { };  static const Iff iff;
	struct PreviouslyNonExistent { };
	  static const PreviouslyNonExistent previously_non_existent;
	struct BasisFunction { };  static const BasisFunction basis_function;
	struct LessThan { };  static const LessThan less_than;
	struct IfLessThan { };  static const IfLessThan if_less_than;
	struct Otherwise { };  static const Otherwise otherwise;
	struct MereSymbol { };  static const MereSymbol mere_symbol;
	struct AssociatedWith { };  static const AssociatedWith associated_with;
	struct Through { };  static const Through through;
	struct Becomes { };  static const Becomes becomes;
	struct Transforms { };  static const Transforms transforms;
	struct Into { };  static const Into into;                                 }
	

class Manifold;  class FiniteElement;
	

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

	inline Function ( ) { assert ( false );  };

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

	inline Function ( const tag::BasisFunction &, const tag::Within &, FiniteElement & );

	inline Function ( const tag::MereSymbol & );

	class ActionGenerator;  // a generator of a discrete group
	// an action will act on functions, particularly on coordinates of a quotient manifold
	typedef tag::Util::Action Action;  // aka  class Manifold::Action
	// we define it here because we need it for Function::MultiValued
	// but we prefer the user to see it as an attribute of class Manifold
	
	inline ~Function();

	inline Function & operator= ( const Function & );
	
	inline size_t nb_of_components ( ) const;

	inline bool exists ( )
	{	return this->core;  }

	inline Function operator[] ( size_t ) const;

	inline void conditionally_dispose_core ( );
	inline void set_core ( Function::Core *, const tag::PreviouslyNonExistent & );
	inline void set_core_to_null ( );
	inline void change_core_to ( Function::Core * );

	inline Function multivalued ( const tag::Through &, const Function::ActionGenerator &,
	                              const tag::Becomes &, const Function &         );
 	inline Function multivalued ( const tag::Through &, const Function::ActionGenerator &,
	                              const tag::Becomes &, const Function &,
	                              const tag::Through &, const Function::ActionGenerator &,
	                              const tag::Becomes &, const Function &         );
	inline Function multivalued ( const tag::Through &, const Function::ActionGenerator &,
	                              const tag::Becomes &, const Function &,
	                              const tag::Through &, const Function::ActionGenerator &,
	                              const tag::Becomes &, const Function &,
	                              const tag::Through &, const Function::ActionGenerator &,
	                              const tag::Becomes &, const Function &         );

	class Jump;
	inline Jump jump();
	
	class TakenOnCell;  class TakenOnWindingCell;
	inline Function::TakenOnCell operator() ( const Cell & cll ) const;
	inline Function::TakenOnWindingCell operator()
	( const Cell & cll, const tag::Winding &, const Function::Action & exp ) const;

	inline Function deriv ( const Function & x ) const;  // derivative with respect to x

	inline Function replace ( const Function & x, const Function & y ) const;
	// in an expression, replace x by y

	#ifndef NDEBUG
	inline std::string repr ( ) const;
	enum From { from_void, from_sum, from_product, from_power, from_function };
	#endif

	static bool less_for_map ( const Function & f, const Function & g );
	// needed for map 'jacobian' in class Function::Map
	// and for map 'equations' in class Manifold::Parametric
	// other solution :
	//namespace std
	//{	template <>
	//	struct less<Function>
	//	{  inline bool operator() ( const Function &, const Function & );  };  }

	class Scalar;  class ArithmeticExpression;  class Constant;
	class Vector;  class Aggregate;  class CoupledWithField;
	class Sum;  class Product;  class Power;  class Sqrt;  class Sin;  class Cos;  class Step;
	class Map;  class Diffeomorphism;  class Immersion;  class Composition;
	class MultiValued;  class MereSymbol;  class DelayedDerivative;
	class Equality;
	struct Inequality
	{ class LessThanZero;  class LessThan;  class GreaterThan;
		typedef tag::Util::InequalitySet Set;                    };

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

	virtual Function deriv ( Function ) = 0;

	virtual Function replace ( const Function & x, const Function & y ) = 0;
	// in an expression, replace x by y

	virtual Function::Jump jump ( );
	// here execution forbidden
	// later overridden by Function::***::MultiValued::JumpIsSum

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
	// never actually used for objects in this class because Function::operator[] returns self

	virtual void set_value ( double ) = 0;  // only for constants

	virtual double get_value_on_cell ( Cell::Core * ) const = 0;
	virtual double get_value_on_cell
	( Cell::Core *, const tag::Winding &, const Function::Action & exp ) const = 0;

	virtual double set_value_on_cell ( Cell::Core *, const double & ) = 0;
	// assign a numeric value to the function on the cell and return that value
	
	// Function deriv ( Function )
	// Function replace ( const Function & x, const Function & y )
	//    stay pure virtual from Function::Core

	// Function::Jump jump ( )  virtual, defined by Function::Core, execution forbidden
	// later overridden by Function::Scalar::MultiValued::JumpIsSum
	
	#ifndef NDEBUG	
	// string repr ( const Function::From & from = Function::from_void )
	//   stays pure virtual from Function::Core
	#endif

	class MultiValued;  

};  // end of class Function::Scalar


//-----------------------------------------------------------------------------------------//

class Function::ArithmeticExpression : public Function::Scalar
	
// base for classes like Constant, Sum, Product
// any function for which 'set_value_on_cell' does not make sense

{	public :

	inline ArithmeticExpression ( ) { };

	inline ArithmeticExpression ( const Function::ArithmeticExpression & ) = delete;
	inline ArithmeticExpression ( Function::ArithmeticExpression && ) = delete;
	
	// size_t nb_of_components ( )
	//   virtual from Function::Core, defined by Function::Scalar, returns 1

	// Function component ( size_t i )  virtual from Function::Core,
	//   defined by Function::Scalar, returns self, never actually used

	// void set_value  stays pure virtual from Function::Scalar
	
	// two versions of
	// double get_value_on_cell ( Cell::Core * )  stay pure virtual from Function::Scalar

	double set_value_on_cell ( Cell::Core *, const double & );  // virtual from Function::Scalar
	// here execution forbidden

	// Function deriv ( Function )
	// Function replace ( const Function & x, const Function & y )
	//   stay pure virtual from Function::Core
	
	// Function::Jump jump ( )  virtual, defined by Function::Core, execution forbidden
	
	#ifndef NDEBUG	
	// string repr ( const Function::From & from = Function::from_void )
	//   stays pure virtual from Function::Core
	#endif

};  // end of  class Function::ArithmeticExpression

//-----------------------------------------------------------------------------------------//


class Function::Constant : public Function::ArithmeticExpression
	
{	public :

	double value;

	inline Constant ( double c )
	:	value { c } { }

	inline Constant ( const Function::Constant & ) = delete;
	inline Constant ( Function::Constant && ) = delete;
	
	inline Function::Constant operator= ( const Function::Constant & ) = delete;
	inline Function::Constant operator= ( Function::Constant && ) = delete;

	// size_t nb_of_components ( )
	//   virtual from Function::Core, defined by Function::Scalar, returns 1

	// Function component ( size_t i )  virtual from Function::Core,
	//   defined by Function::Scalar, returns self, never actually used

	void set_value ( double );  // virtual from Function::Scalar

	double get_value_on_cell ( Cell::Core * ) const;
	double get_value_on_cell
	( Cell::Core *, const tag::Winding &, const Function::Action & exp ) const;
	// virtual from Function::Scalar

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   virtual from Function::Scalar
	//   defined by Function::ArithmeticExpression (execution forbidden)

	Function deriv ( Function );  // virtual from Function::Core

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core
	
	// Function::Jump jump ( )  virtual, defined by Function::Core, execution forbidden
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core
	#endif

};  // end of  class Function::Constant

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

	// size_t nb_of_components ( )
	//   virtual from Function::Core, defined by Function::Scalar, returns 1

	// Function component ( size_t i )  virtual from Function::Core,
	//   defined by Function::Scalar, returns self, never actually used

	void set_value ( double );  // virtual from Function::Scalar, execution forbidden

	double get_value_on_cell ( Cell::Core * ) const;
	double get_value_on_cell
	( Cell::Core *, const tag::Winding &, const Function::Action & exp ) const;
	// virtual from Function::Scalar

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   virtual from Function::Scalar
	//   defined by Function::ArithmeticExpression (execution forbidden)

	Function deriv ( Function );  // virtual from Function::Core

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core
	
	// Function::Jump jump ( )  virtual, defined by Function::Core, execution forbidden
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	//  virtual from Function::Core
	#endif

};  // end of  class Function::Sum

//-----------------------------------------------------------------------------------------//

class Function::Product : public Function::ArithmeticExpression
	
{	public :

	std::forward_list < Function > factors;  // we use wrappers as a pointers

	inline Product ( ) { };

	inline Product ( const Function::Product & ) = delete;
	inline Product ( Function::Product && ) = delete;
	
	inline Function::Product operator= ( const Function::Product & ) = delete;
	inline Function::Product operator= ( Function::Product && ) = delete;

	// size_t nb_of_components ( )
	//   virtual from Function::Core, defined by Function::Scalar, returns 1

	// Function component ( size_t i )  virtual from Function::Core,
	//   defined by Function::Scalar, returns self, never actually used

	void set_value ( double );  // virtual from Function::Scalar, execution forbidden

	double get_value_on_cell ( Cell::Core * ) const;
	double get_value_on_cell
	( Cell::Core *, const tag::Winding &, const Function::Action & exp ) const;
	// virtual from Function::Scalar

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   virtual from Function::Scalar
	//   defined by Function::ArithmeticExpression (execution forbidden)

	Function deriv ( Function );  // virtual from Function::Core

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core
	
	// Function::Jump jump ( )  virtual, defined by Function::Core, execution forbidden
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	//  virtual from Function::Core
	#endif

};  // end of  class Function::Product

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

	// size_t nb_of_components ( )
	//   virtual from Function::Core, defined by Function::Scalar, returns 1

	// Function component ( size_t i )  virtual from Function::Core,
	//   defined by Function::Scalar, returns self, never actually used

	void set_value ( double );  // virtual from Function::Scalar, execution forbidden

	double get_value_on_cell ( Cell::Core * ) const;
	double get_value_on_cell
	( Cell::Core *, const tag::Winding &, const Function::Action & exp ) const;
	// virtual from Function::Scalar

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   virtual from Function::Scalar
	//   defined by Function::ArithmeticExpression (execution forbidden)

	Function deriv ( Function );  // virtual from Function::Core

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core
	
	// Function::Jump jump ( )  virtual, defined by Function::Core, execution forbidden
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	//  virtual from Function::Core
	#endif

};  // end of  class Function::Power

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
	
	// size_t nb_of_components ( )
	//   virtual from Function::Core, defined by Function::Scalar, returns 1

	// Function component ( size_t i )  virtual from Function::Core,
	//   defined by Function::Scalar, returns self, never actually used

	void set_value ( double );  // virtual from Function::Scalar, execution forbidden

	double get_value_on_cell ( Cell::Core * ) const;
	double get_value_on_cell
	( Cell::Core *, const tag::Winding &, const Function::Action & exp ) const;
	// virtual from Function::Scalar

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   virtual from Function::Scalar
	//   defined by Function::ArithmeticExpression (execution forbidden)

	Function deriv ( Function );  // virtual from Function::Core

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core
	
	// Function::Jump jump ( )  virtual, defined by Function::Core, execution forbidden
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	//  virtual from Function::Core
	#endif

};  // end of  class Function::Sqrt

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
	
	// size_t nb_of_components ( )
	//   virtual from Function::Core, defined by Function::Scalar, returns 1

	// Function component ( size_t i )  virtual from Function::Core,
	//   defined by Function::Scalar, returns self, never actually used

	void set_value ( double );  // virtual from Function::Scalar, execution forbidden

	double get_value_on_cell ( Cell::Core * ) const;
	double get_value_on_cell
	( Cell::Core *, const tag::Winding &, const Function::Action & exp ) const;
	// virtual from Function::Scalar

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   virtual from Function::Scalar
	//   defined by Function::ArithmeticExpression (execution forbidden)

	Function deriv ( Function );  // virtual from Function::Core

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core
	
	// Function::Jump jump ( )  virtual, defined by Function::Core, execution forbidden
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	//  virtual from Function::Core
	#endif

};  // end of  class Function::Sin

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

	// size_t nb_of_components ( )
	//   virtual from Function::Core, defined by Function::Scalar, returns 1

	// Function component ( size_t i )  virtual from Function::Core,
	//   defined by Function::Scalar, returns self, never actually used

	void set_value ( double );  // virtual from Function::Scalar, execution forbidden

	double get_value_on_cell ( Cell::Core * ) const;
	double get_value_on_cell
	( Cell::Core *, const tag::Winding &, const Function::Action & exp ) const;
	// virtual from Function::Scalar

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   virtual from Function::Scalar
	//   defined by Function::ArithmeticExpression (execution forbidden)

	Function deriv ( Function );  // virtual from Function::Core

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core
	
	// Function::Jump jump ( )  virtual, defined by Function::Core, execution forbidden
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	//  virtual from Function::Core
	#endif

};  // end of  class Function::Cos

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
	{	assert ( vals .size() == cts .size() + 1 );  }  // add asserts !!

	inline Step ( const Function::Step & ) = delete;
	inline Step ( Function::Step && ) = delete;
	
	inline Function::Step operator= ( const Function::Step & ) = delete;
	inline Function::Step operator= ( Function::Step && ) = delete;

	// size_t nb_of_components ( )
	//   virtual from Function::Core, defined by Function::Scalar, returns 1

	// Function component ( size_t i )  virtual from Function::Core,
	//   defined by Function::Scalar, returns self, never actually used

	void set_value ( double );  // virtual from Function::Scalar, execution forbidden

	double get_value_on_cell ( Cell::Core * ) const;
	double get_value_on_cell
	( Cell::Core *, const tag::Winding &, const Function::Action & exp ) const;
	// virtual from Function::Scalar

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   virtual from Function::Scalar
	//   defined by Function::ArithmeticExpression (execution forbidden)

	Function deriv ( Function );  // virtual from Function::Core

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core
	
	// Function::Jump jump ( )  virtual, defined by Function::Core, execution forbidden
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	//  virtual from Function::Core
	#endif

}; // end of  class Function::Step

//-----------------------------------------------------------------------------------------//

class Function::MereSymbol : public Function::Scalar

// used by some finite elements as (slack) basis function
	
{	public :

	// no data

	inline MereSymbol ( ) { };

	inline MereSymbol ( const Function::MereSymbol & ) = delete;
	inline MereSymbol ( Function::MereSymbol && ) = delete;
	
	inline Function::MereSymbol operator= ( const Function::MereSymbol & ) = delete;
	inline Function::MereSymbol operator= ( Function::MereSymbol && ) = delete;

	// size_t nb_of_components ( )
	//   virtual from Function::Core, defined by Function::Scalar, returns 1

	// Function component ( size_t i )  virtual from Function::Core,
	//   defined by Function::Scalar, returns self, never actually used

	void set_value ( double );  // virtual from Function::Scalar, here execution forbidden

	double get_value_on_cell ( Cell::Core * ) const;
	double get_value_on_cell
	( Cell::Core *, const tag::Winding &, const Function::Action & exp ) const;
	// virtual from Function::Scalar, here execution forbidden

	double set_value_on_cell ( Cell::Core *, const double & );
	//   virtual from Function::Scalar, here execution forbidden

	Function deriv ( Function );  // virtual from Function::Core

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core
	
	// Function::Jump jump ( )  virtual, defined by Function::Core, execution forbidden
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core
	#endif

};  // end of  class Function::MereSymbol

//-----------------------------------------------------------------------------------------//

inline Function::Function
( const tag::BasisFunction &, const tag::Within &, FiniteElement & )
:	Function ( tag::whose_core_is, new Function::MereSymbol )
{	}

inline Function::Function ( const tag::MereSymbol & )
:	Function ( tag::whose_core_is, new Function::MereSymbol )
{	}

//-----------------------------------------------------------------------------------------//


class Function::DelayedDerivative : public Function::Scalar

// a function 'base' differentiated with respect to 'variable'

{	public :

	Function base, variable;  // we use wrappers as a pointers

	inline DelayedDerivative ( const Function & b, const Function & v )
	:	base { b }, variable { v }
	{	}

	inline DelayedDerivative ( const Function::DelayedDerivative & ) = delete;
	inline DelayedDerivative ( Function::DelayedDerivative && ) = delete;
	
	inline Function::DelayedDerivative operator= ( const Function::DelayedDerivative & ) = delete;
	inline Function::DelayedDerivative operator= ( Function::DelayedDerivative && ) = delete;

	// size_t nb_of_components ( )
	//   virtual from Function::Core, defined by Function::Scalar, returns 1

	// Function component ( size_t i )  virtual from Function::Core,
	//   defined by Function::Scalar, returns self, never actually used

	void set_value ( double );  // virtual from Function::Scalar, here execution forbidden

	double get_value_on_cell ( Cell::Core * ) const;
	double get_value_on_cell
	( Cell::Core *, const tag::Winding &, const Function::Action & exp ) const;
	// virtual from Function::Scalar, here execution forbidden

	double set_value_on_cell ( Cell::Core *, const double & );
	//   virtual from Function::Scalar, here execution forbidden

	Function deriv ( Function );  // virtual from Function::Core

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core
	
	// Function::Jump jump ( )  virtual, defined by Function::Core, execution forbidden
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core
	#endif

};  // end of  class Function::DelayedDerivative

//-----------------------------------------------------------------------------------------//


class Function::Vector : public Function::Core
	
{	public :

	inline Vector ( ) { };

	inline Vector ( const Function::Vector & ) = delete;
	inline Vector ( Function::Vector && ) = delete;
	
	// size_t nb_of_components ( )  stays pure virtual from Function::Core
	
	// Function component ( )  stays pure virtual from Function::Core

	virtual std::vector<double> get_value_on_cell ( Cell::Core * ) const = 0;
	virtual std::vector<double> get_value_on_cell
	( Cell::Core *, const tag::Winding &, const Function::Action & exp ) const = 0;

	virtual std::vector<double> set_value_on_cell
	( Cell::Core *, const std::vector<double> & ) = 0;
	// assign a numeric vector to the function on the cell and return that vector
	
	Function deriv ( Function );
	//  virtual from Function::Core, here execution forbidden, to change

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, here execution forbidden
	
	// Function::Jump jump ( )  virtual, defined by Function::Core, execution forbidden
	// later overridden by Function::Vector::MultiValued::JumpIsSum
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	//  virtual from Function::Core, here forbids execution
	#endif

	class MultiValued;
	
};  // end of class Function::Vector

//-----------------------------------------------------------------------------------------//

class Function::Aggregate : public Function::Vector
	
// inheriting from this class means simply that there is a 'components' member

{	public :

	std::vector < Function > components;
	
	inline Aggregate ( const tag::ReserveSize &, size_t s )
	:	Function::Vector()
	{	components.reserve(s);  }

	inline Aggregate ( const Function::Aggregate & ) = delete;
	inline Aggregate ( Function::Aggregate && ) = delete;
	
	inline Function::Aggregate operator= ( const Function::Aggregate & ) = delete;
	inline Function::Aggregate operator= ( Function::Aggregate && ) = delete;

	// 'nb_of_components' and 'component' are virtual from Function::Core
	// later overridden by Function::CoupledWithField::Vector
	size_t nb_of_components ( ) const;
	Function component ( size_t i );
	
	std::vector<double> get_value_on_cell ( Cell::Core * ) const;
	std::vector<double> get_value_on_cell
	( Cell::Core *, const tag::Winding &, const Function::Action & exp ) const;
	// virtual from Function::Vector

	std::vector<double> set_value_on_cell ( Cell::Core *, const std::vector<double> & );
	// virtual from Function::Vector

	// Function deriv ( Function )
	//    defined by Function::Vector (execution forbidden), to change

	// Function replace ( const Function & x, const Function & y );
	//    defined by Function::Vector (execution forbidden), to change
	
	// Function::Jump jump ( )  virtual, defined by Function::Core, execution forbidden
	
	#ifndef NDEBUG	
	// std::string repr ( const Function::From & from = Function::from_void )
	//    defined by Function::Vector (execution forbidden)
	#endif

};  // end of  class Function::Aggregate


//-----------------------------------------------------------------------------------------//

inline Function operator&& ( Function f, Function g )

{	size_t nf = f .nb_of_components(), ng = g .nb_of_components();
	Function::Aggregate * res = new Function::Aggregate ( tag::reserve_size, nf + ng ) ;
	for ( size_t i = 0; i < nf; i++ ) res->components .emplace_back ( f [i] );
	for ( size_t i = 0; i < ng; i++ ) res->components .emplace_back ( g [i] );
	return Function ( tag::whose_core_is,	res );                                         }

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

	std::map < Function, Function, bool (*) ( const Function &, const Function & ) > jacobian;
	//  decltype(Function::less_for_map)*  equals  bool (*) ( const Function &, const Function & )
	
	Function det;  // positive scalar, dilation coefficient

	inline Map ( const Function & gc, const Function & mc, const Function & bgc )
	:	geom_coords ( gc ), master_coords ( mc ), back_geom_coords ( bgc ),
		jacobian ( & Function::less_for_map ), det( tag::non_existent )
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
	
	// size_t nb_of_components ( )
	//   virtual from Function::Core, defined by Function::Scalar, returns 1

	// Function component ( size_t i )  virtual from Function::Core,
	//   defined by Function::Scalar, returns self, never actually used

	void set_value ( double );  // virtual from Function::Scalar, execution forbidden

	double get_value_on_cell ( Cell::Core * ) const;
	double get_value_on_cell
	( Cell::Core *, const tag::Winding &, const Function::Action & exp ) const;
	// virtual from Function::Scalar

	double set_value_on_cell ( Cell::Core *, const double & );
	// virtual from Function::Vector

	Function deriv ( Function );  // virtual from Function::Core

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core
	
	// Function::Jump jump ( )  virtual, defined by Function::Core, execution forbidden
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	//  virtual from Function::Core
	#endif

};  // end of class Function::Diffeomorphism::OneDim
	
//-----------------------------------------------------------------------------------------//


inline Function::Function ( const tag::Diffeomorphism &, const tag::OneDim &,
                            const Function & geom_coords, const Function & master_coords,
                            const Function & back_geom_coords                             )
: Function ( tag::non_existent )
{	assert ( geom_coords .nb_of_components() == back_geom_coords .nb_of_components() );
	assert ( master_coords .nb_of_components() == geom_coords .nb_of_components() );
	assert ( geom_coords .nb_of_components() == 1 );
	this->set_core ( new Function::Diffeomorphism::OneDim
	   ( geom_coords, master_coords, back_geom_coords, tag::build_jacobian ),
	   tag::previously_non_existent                                           );       }
	

//-----------------------------------------------------------------------------------------//

class Function::Immersion : public Function::Vector, public Function::Map

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

	size_t nb_of_components ( ) const;  // virtual from Function::Core

	Function component ( size_t i );  // virtual from Function::Core
	
	std::vector<double> get_value_on_cell ( Cell::Core * ) const;
	std::vector<double> get_value_on_cell
	( Cell::Core *, const tag::Winding &, const Function::Action & exp ) const;
	// virtual from Function::Vector

	std::vector<double> set_value_on_cell ( Cell::Core *, const std::vector<double> & );
	// virtual from Function::Vector

	// Function deriv ( Function )  defined by Function::Vector (execution forbidden), to change

	// Function replace ( const Function & x, const Function & y );
	//    defined by Function::Vector (execution forbidden), to change
	
	// Function::Jump jump ( )  virtual, defined by Function::Core, execution forbidden
	
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
	// Function component ( size_t i )   defined by Function::Immersion, execution forbidden
	
	// the two methods below defined by Function::Immersion (execution forbidden)
	// std::vector<double> set_value_on_cell ( Cell::Core *, const std::vector < double > & )
	// std::vector<double> get_value_on_cell ( Cell::Core * ) const

	// Function deriv ( Function )
	//    defined by Function::Vector (execution forbidden)

	// Function replace ( const Function & x, const Function & y );
	//    defined by Function::Vector (execution forbidden)
	
	// Function::Jump jump ( )  virtual, defined by Function::Core, execution forbidden
	
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
{	assert ( geom_coords .nb_of_components() == back_geom_coords .nb_of_components() );
	assert ( master_coords .nb_of_components() == geom_coords .nb_of_components() );
	assert ( geom_coords .nb_of_components() >= 2 );
	this->set_core ( new Function::Diffeomorphism::HighDim
	   ( geom_coords, master_coords, back_geom_coords, tag::build_jacobian ),
	   tag::previously_non_existent                                           );       }

	
inline Function::Function ( const tag::Diffeomorphism &, const Function & geom_coords,
                            const Function & master_coords, const Function & back_geom_coords )
: Function ( tag::non_existent )
{	assert ( geom_coords .nb_of_components() == back_geom_coords .nb_of_components() );
	assert ( master_coords .nb_of_components() == geom_coords .nb_of_components() );
	if ( geom_coords .nb_of_components() == 1 )
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

	// size_t nb_of_components ( )
	//   virtual from Function::Core, defined by Function::Scalar, returns 1

	// Function component ( size_t i )  virtual from Function::Core,
	//   defined by Function::Scalar, returns self, never actually used

	void set_value ( double );  // virtual from Function::Scalar, execution forbidden

	double get_value_on_cell ( Cell::Core * ) const;
	double get_value_on_cell
	( Cell::Core *, const tag::Winding &, const Function::Action & exp ) const;
	double set_value_on_cell ( Cell::Core *, const double & );
	// virtual from Function::Scalar

	Function deriv ( Function );  // virtual from Function::Core

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core
	
	// Function::Jump jump ( )  virtual, defined by Function::Core, execution forbidden
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	//  virtual from Function::Core
	#endif

};  // end of  class Function::Composition

//-----------------------------------------------------------------------------------------//


inline Function::Function
(	const Function & expr, const tag::ComposedWith &, const Function & f_map )

// an 'expr'ession involving master coordinates ( e.g.  1.- xi - eta )
// composed with a map 'f_map' sending it in the physical space

: Function ( tag::non_existent )

{	assert ( dynamic_cast < Function::Map* > ( f_map .core ) );
  Function::Constant * expr_c = dynamic_cast < Function::Constant* > ( expr .core );
	if ( expr_c ) this->set_core ( expr .core, tag::previously_non_existent );
	else this->set_core ( new Function::Composition ( expr, f_map ),
                        tag::previously_non_existent               );                }

//-----------------------------------------------------------------------------------------//


class Function::CoupledWithField
	
{	public :

	Field::Double::Core * field;

	inline CoupledWithField ( Field::Double::Core * f ) : field { f }
	{	assert ( f );  }
		
	class Scalar;  class Vector;
};

//-----------------------------------------------------------------------------------------//

class Function::CoupledWithField::Scalar
: public Function::Scalar,
	public Function::CoupledWithField
	
{	public :

	inline Scalar ( Field::Double::Scalar * f )
	:	Function::Scalar(),
		Function::CoupledWithField ( f )
	{	assert ( f->nb_of_components() == 1 );  }
		
	inline Scalar ( const Function::CoupledWithField::Scalar & ) = delete;
	inline Scalar ( Function::CoupledWithField::Scalar && ) = delete;
	
	inline Function::CoupledWithField::Scalar operator=
		( const Function::CoupledWithField::Scalar & ) = delete;
	inline Function::CoupledWithField::Scalar operator=
		( Function::CoupledWithField::Scalar && ) = delete;

	// size_t nb_of_components ( )
	//   virtual from Function::Core, defined by Function::Scalar, returns 1

	// Function component ( size_t i )  virtual from Function::Core,
	//   defined by Function::Scalar, returns self, never actually used

	void set_value ( double );  // virtual from Function::Scalar, execution forbidden

	double get_value_on_cell ( Cell::Core * ) const;
	double get_value_on_cell
	( Cell::Core *, const tag::Winding &, const Function::Action & exp ) const;
	double set_value_on_cell ( Cell::Core *, const double & );
	// virtual from Function::Scalar

	Function deriv ( Function );  // virtual from Function::Core

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core
	
	// Function::Jump jump ( )  virtual, defined by Function::Core, execution forbidden
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	//  virtual from Function::Core
	#endif

};  // end of  class Function::CoupledWithField::Scalar

//-----------------------------------------------------------------------------------------//

class Function::CoupledWithField::Vector
: public Function::Aggregate,
	public Function::CoupledWithField

// inheriting from Function::Aggregate means simply that there is a 'components' member

// quando calculamos as componentes duma funcao associada a um campo vectorial,
// nao precisamos de criar as componentes do campo vectorial
// podemos ter uma funcao escalar associada a um indice dentro dum campo vectorial Block
	
{	public :

	inline Vector ( Field::Double::Block * f )
	:	Function::Aggregate ( tag::reserve_size, f->nb_of_components() ),
		Function::CoupledWithField ( f )
	{	size_t n = f->nb_of_components();
		for ( size_t j = 0; j < n; j++ )
//		this->components [j] = Function ( tag::whose_core_is,
			this->components .emplace_back ( tag::whose_core_is,
				new Function::CoupledWithField::Scalar ( field->component (j) ) );  }
	
	inline Vector ( const Function::CoupledWithField::Vector & ) = delete;
	inline Vector ( Function::CoupledWithField::Vector && ) = delete;
	
	inline Function::CoupledWithField::Vector operator=
		( const Function::CoupledWithField::Vector & ) = delete;
	inline Function::CoupledWithField::Vector operator=
		( Function::CoupledWithField::Vector && ) = delete;

	virtual size_t nb_of_components ( ) const override;
	// virtual from Function::Core
	
	virtual Function component ( size_t i ) override;
	// virtual from Function::Core
	
	std::vector<double> get_value_on_cell ( Cell::Core * ) const;
	std::vector<double> get_value_on_cell
	( Cell::Core *, const tag::Winding &, const Function::Action & exp ) const;
	std::vector<double> set_value_on_cell ( Cell::Core *, const std::vector<double> & );
	// virtual from Function::Vector

	// Function deriv ( Function )
	//    defined by Function::Vector (execution forbidden), may change

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, through Function::Vector
	
	// Function::Jump jump ( )  virtual, defined by Function::Core, execution forbidden
	
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
	return minus_one * f;       }

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

inline Function abs ( const Function & f )
{	return Function ( tag::whose_core_is, new Function::Step
      ( -f, tag::iff, f, tag::less_than, 0., f, tag::otherwise ) );  }

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
//-----------------------------------------------------------------------------------------//


class Function::ActionGenerator

// a generator of a discrete group

{	public :

	static size_t counter;
	size_t id;
	
	Function coords, transf;
	// here we momentarily keep the coordinates function together with the composed ones
	// just prior to the declaration of the respective quotient manifold
	// the quotient manifold will then keep these coordinates
	// after that, we can forget about them

	inline ActionGenerator ( )
	:	id { Function::ActionGenerator::counter },
		coords ( tag::non_existent ), transf ( tag::non_existent )
	{	Function::ActionGenerator::counter++;  }

	inline ActionGenerator
	( const tag::Transforms &, const Function f, const tag::Into &, const Function g )
	:	id { Function::ActionGenerator::counter }, coords ( f ), transf ( g )
	{	Function::ActionGenerator::counter++;  }

	inline ActionGenerator ( const Function::ActionGenerator & a )
	:	id { a .id }, coords ( a .coords ), transf ( a .transf )
	{	}

	inline ActionGenerator operator= ( const Function::ActionGenerator & a )
	{	this->id = a .id;  return *this;  }
	
	inline operator Function::Action() const;

	struct Applied { class ToFunction;  };	

};  // end of class Function::ActionGenerator


inline bool operator== ( const Function::ActionGenerator & a, const Function::ActionGenerator & b )
{	return a .id == b .id;  }


inline bool operator< ( const Function::ActionGenerator & f, const Function::ActionGenerator & g )
{	return f .id < g .id;  }
// needed for Function::Action::index_map

//-----------------------------------------------------------------------------------------//


class tag::Util::Action  //  aka  class Function::Action, aka class Manfold::Action
		// we define it here because we need it for Function::MultiValued
		// but we prefer the user to see it as an attribute of class Manifold

// a composition of actions, thus an element of the group
// essentially, a multi-index, more precisely a map < Function::ActionGenerator, short int >

// the final user will compose actions in an additive notation,
// e.g.  2*g1 - g2  rather than  g1^2 g2^-1

{	public :

	std::map < Function::ActionGenerator, short int > index_map;

	inline Action ( ) { };
		
	inline Action ( short int i )
	{	assert ( i == 0 );  }
	
	inline Action ( const Function::ActionGenerator & g, short int i )
	:	index_map { std::pair < Function::ActionGenerator, short int > ( g, i ) }
	{	}

	inline Action ( const Function::ActionGenerator & g, short int i,
                          const Function::ActionGenerator & h, short int j )
	:	index_map { std::pair < Function::ActionGenerator, short int > ( g, i ),
		            std::pair < Function::ActionGenerator, short int > ( h, j )  }
	{	}

	inline Action ( const Function::ActionGenerator & g1, short int i,
                          const Function::ActionGenerator & g2, short int j,
                          const Function::ActionGenerator & g3, short int k )
	:	index_map { std::pair < Function::ActionGenerator, short int > ( g1, i ),
		            std::pair < Function::ActionGenerator, short int > ( g2, j ),
		            std::pair < Function::ActionGenerator, short int > ( g3, k )  }
	{	}

	inline Action ( const tag::Transforms &, const Function f, const tag::Into &, const Function g )
	{	this->index_map.emplace ( std::piecewise_construct,
		  std::forward_as_tuple( tag::transforms, f, tag::into, g ), std::forward_as_tuple(1) );  }
		//	{	Function::ActionGenerator a ( tag::transforms, f, tag::into, g );
		//		this->index_map [a] = 1;                                           }

	inline Function::Action operator= ( const short int zero );

};  // end of class Function::Action


inline bool operator== ( const Function::Action & a, const Function::Action & b )

{	{ // just a block for hiding names
	std::map <Function::ActionGenerator, short int > ::const_iterator it_a =
		a .index_map .begin();
	for ( ; it_a != a .index_map .end(); it_a++ )
	{	Function::ActionGenerator aa = it_a->first;
		std::map < Function::ActionGenerator, short int > ::const_iterator it_b =
			b .index_map .find ( aa );
		if ( it_b == b .index_map .end() ) return false;
		if ( it_a->second != it_b->second ) return false;                              }
	} { // just a block for hiding names
	std::map < Function::ActionGenerator, short int > ::const_iterator it_b =
		b .index_map .begin();
	for ( ; it_b != b .index_map .end(); it_b++ )
	{	Function::ActionGenerator bb = it_b->first;
		std::map < Function::ActionGenerator, short int > ::const_iterator it_a =
			a .index_map .find ( bb );
		if ( it_a == a .index_map .end() ) return false;
		assert ( it_a->second == it_b->second );                                       }
	} // just a block for hiding names
	return true;                                                                       }


inline bool operator!=
( const Function::Action & a, const Function::Action & b )
{	return not ( a == b );  }


inline bool operator== ( const Function::Action & a, const short int zero )
// shorthand for quering the identity action :  id == 0
{	assert ( zero == 0 );
	return a .index_map .size() == 0;  }


inline bool operator!=
( const Function::Action & a, const short int zero )
{	return not ( a == zero );  }
	

inline Function::Action Function::Action::operator=
( const short int zero )
// shorthand for quering the identity action :  id = 0
{	assert ( zero == 0 );
	this->index_map .clear();
	return *this;             }

	
inline Function::ActionGenerator::operator Function::Action() const
{	return Function::Action ( *this, 1 );  }


inline Function::Action operator+
( const Function::Action & a, const Function::Action & b )
{	Function::Action res = a;
	std::map < Function::ActionGenerator, short int > ::const_iterator it =
		b .index_map .begin();
	for ( ; it != b .index_map .end(); it++ )
	{	const Function::ActionGenerator & g = it->first;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Function::ActionGenerator, short int > ::iterator itt =
			res .index_map .lower_bound (g);
		if ( ( itt == res .index_map .end() ) or
		     ( res .index_map .key_comp() ( g, itt->first ) ) )
			// new action
			res .index_map .emplace_hint ( itt, std::piecewise_construct,
	      std::forward_as_tuple (g), std::forward_as_tuple ( it->second ) ); 
		else  // action already there
		{	itt->second += it->second;  // could be zero
			if ( itt->second == 0 ) res .index_map .erase ( itt );  }            }
	return res;                                                                 }


inline Function::Action operator+=
( Function::Action & a, const Function::Action & b )
{	std::map < Function::ActionGenerator, short int > ::const_iterator it =
		b .index_map .begin();
	for ( ; it != b .index_map .end(); it++ )
	{	const Function::ActionGenerator & g = it->first;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Function::ActionGenerator, short int > ::iterator itt =
			a .index_map .lower_bound (g);
		if ( ( itt == a .index_map .end() ) or
	       ( a .index_map .key_comp() ( g, itt->first ) ) )
			// new action
			a .index_map .emplace_hint ( itt, std::piecewise_construct,
	      std::forward_as_tuple (g), std::forward_as_tuple ( it->second ) ); 
		else  // action already there
		{	itt->second += it->second;  // could be zero
			if ( itt->second == 0 ) a .index_map .erase ( itt );  }              }
	return a;                                                                   }


inline Function::Action operator-
( const Function::Action & a, const Function::Action & b )
{	Function::Action res = a;
	std::map < Function::ActionGenerator, short int > ::const_iterator it =
		b .index_map .begin();
	for ( ; it != b .index_map .end(); it++ )
	{	const Function::ActionGenerator & g = it->first;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Function::ActionGenerator, short int > ::iterator itt =
			res .index_map .lower_bound (g);
		if ( ( itt == res .index_map .end() ) or
	       ( res .index_map .key_comp() ( g, itt->first ) ) )
			// new action
			res .index_map .emplace_hint ( itt, std::piecewise_construct,
	      std::forward_as_tuple (g), std::forward_as_tuple ( -it->second ) ); 
		else  // action already there
		{	itt->second -= it->second;  // could be zero
			if ( itt->second == 0 ) res .index_map .erase ( itt );  }              }
	return res;                                                                  }


inline Function::Action operator-=
( Function::Action & a, const Function::Action & b )
{	std::map < Function::ActionGenerator, short int > ::const_iterator it =
		b .index_map .begin();
	for ( ; it != b .index_map .end(); it++ )
	{	const Function::ActionGenerator & g = it->first;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Function::ActionGenerator, short int > ::iterator itt =
			a .index_map .lower_bound (g);
		if ( ( itt == a .index_map .end() ) or
	       ( a .index_map .key_comp() ( g, itt->first ) ) )
			// new action
			a .index_map .emplace_hint ( itt, std::piecewise_construct,
	      std::forward_as_tuple (g), std::forward_as_tuple ( -it->second ) ); 
		else  // action already there
		{	itt->second -= it->second;  // could be zero
			if ( itt->second == 0 ) a .index_map .erase ( itt );  }               }
	return a;                                                                   }


inline Function::Action operator*
( const short int k, const Function::Action & a )
{	if ( k == 0 ) return Function::Action ( 0 );
	Function::Action res = a;
	std::map < Function::ActionGenerator, short int > ::iterator it =
		res .index_map .begin();
	for ( ; it != res .index_map .end(); it++ ) it->second *= k;
	return res;                                                        }


inline Function::Action operator-
( const Function::Action & a )
{	return (-1) * a;  }

//---------------------------------------------------------------------------------------

	
class Cell::Winding  // inutil ?

// temporary object returned by Cell::winding() and used for assingment of winding numbers

// only positive segments have space reserved for holding a winding number
// however, the arithmetic operators of this class are crafted in a manner
// which allows getting and setting the winding number of a negative segment, too
// reversed arithmetic operations are performed on the reverse, positive, segment

{	public :

	Cell::Core * cll;  // usually a segment
	// Cell::Winding should only be used as temporary objects
	// they should be immediately converted to a (reference to a) double or vector<double>
	// so an ordinary pointer is OK here

	inline Winding ( const Cell & c )
	:	cll { c.core }
	{ assert ( this->cll );  }
	
	inline Function::Action operator=
	( const Function::Action & a );  // defined in manifold.h
	inline Function::Action operator+=
	( const Function::Action & a );  // defined in manifold.h
	inline Function::Action operator-=
	( const Function::Action & a );  // defined in manifold.h
	
	inline operator Function::Action ( );

};

//-----------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------//


class Function::MultiValued

// classes Function::Scalar::MultiValued and Function::Vector::MultiValued inherit from here

{	public :

	Function base;

	std::vector < Function::ActionGenerator > actions;

	inline MultiValued ( const Function & b, std::vector < Function::ActionGenerator > a )
	:	base { b }, actions { a }
	{	// assert ( Manifold::working.actions == a );
	}

};  // end of class Function::Scalar::Multivalued

//-----------------------------------------------------------------------------------------


class Function::Scalar::MultiValued : public Function::MultiValued, public Function::Scalar

// here (finally) the method get_value_on_cell with tag::winding is meaningful
// suppose 'exp' is a pair of short integers (i,j)
// then the above refered method checks that the 'actions' match
// those of the current working manifold
// then takes the value of 'base' on the cell, applies the first action 'i' times
// then applies the second action 'j' times (recall the group should be commutative)
// the vector of 'actions' is only used for the above refered checking operation

// abstract class, specialized in Function::Scalar::MultiValued::JumpIsSum and JumpIsLinear

{	public :

	// members inherited from Function::MultiValued :
	// Function base  -- here must be Function::Scalar
	// std::vector < Function::ActionGenerator > actions

	inline MultiValued ( const tag::AssociatedWith &, const Function & b,
	                     std::vector < Function::ActionGenerator > a               )
	:	Function::MultiValued ( b, a )
	{	}

	inline MultiValued ( const Function::Scalar::MultiValued & ) = delete;
	inline MultiValued ( Function::Scalar::MultiValued && ) = delete;
	
	// size_t nb_of_components ( )
	//   virtual from Function::Core, defined by Function::Scalar, returns 1

	// Function component ( size_t i )  virtual from Function::Core,
	//   defined by Function::Scalar, returns self, never actually used

	void set_value ( double );  // virtual from Function::Scalar, delegates to base

	double get_value_on_cell ( Cell::Core * ) const;
	//  virtual from Function::Scalar, delegates to base

	// double get_value_on_cell  with tag::winding  stays pure virtual from Function::Scalar

	double set_value_on_cell ( Cell::Core *, const double & );
	//   virtual from Function::Scalar, here execution forbidden

	Function deriv ( Function );  // virtual from Function::Core

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core
	
	// Function::Jump jump ( )  virtual, defined by Function::Core, execution forbidden
	// later overridden by Function::Scalar::MultiValued::JumpIsSum
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	//  virtual from Function::Core
	#endif

	class JumpIsSum;  class JumpIsLinear;

};  // end of class Function::Scalar::MultiValued

//-----------------------------------------------------------------------------------------//


class Function::Scalar::MultiValued::JumpIsSum : public Function::Scalar::MultiValued

// here the actions are mere translations on RR (sums)

{	public :

	// members inherited from Function::MultiValued :
	// Function base  -- here must be Function::Scalar
	// std::vector < Function::ActionGenerator > actions
	
	std::vector < double > beta;
	// upon each action [i], the value v of the function becomes  v + beta [i]
	// use a Function::Jump::Sum::Scalar instead !
	
	inline JumpIsSum ( const tag::AssociatedWith &, const Function & f,
										 std::vector < Function::ActionGenerator > ac, std::vector < double > be )
	:	Function::Scalar::MultiValued ( tag::associated_with, f, ac ), beta { be }
	{	assert ( ac .size() == be .size() );  }

	inline JumpIsSum ( const Function::Scalar::MultiValued::JumpIsSum & ) = delete;
	inline JumpIsSum ( Function::Scalar::MultiValued::JumpIsSum && ) = delete;
	
	inline Function::Scalar::MultiValued::JumpIsSum operator=
		( const Function::Scalar::MultiValued::JumpIsSum & ) = delete;
	inline Function::Scalar::MultiValued::JumpIsSum operator=
		( Function::Scalar::MultiValued::JumpIsSum && ) = delete;

	// size_t nb_of_components ( )
	//   virtual from Function::Core, defined by Function::Scalar, returns 1

	// Function component ( size_t i )  virtual from Function::Core,
	//   defined by Function::Scalar, returns self, never actually used

	// void set_value ( double )
	//   defined by Function::Scalar::Multivalued, delegates to base

	// double get_value_on_cell ( Cell::Core * ) const;
	//   defined by Function::Scalar::MultiValued, delegates to base

	double get_value_on_cell
	( Cell::Core *, const tag::Winding &, const Function::Action & exp ) const;
	//  virtual from Function::Scalar
	
	// double set_value_on_cell ( Cell::Core *, const double & )
	//   defined by Function::Scalar::MultiValued, execution forbidden

	Function::Jump jump ( );
  // virtual, defined by Function::Core, execution forbidden, here overridden
	
	// the return value of 'analyse_linear_expression' should be double
	// however, we use std::vector < double > instead
	// a zero-length vector means not succeeded
	// success is represented by a one-length vector
	inline static std::vector < double > analyse_linear_expression
		( Function expression, Function base );
	
};  // end of class Function::Scalar::MultiValued::JumpIsSum

//-----------------------------------------------------------------------------------------//


class Function::Scalar::MultiValued::JumpIsLinear : public Function::Scalar::MultiValued

// here the actions are linear (affine) maps on RR

{	public :

	// members inherited from Function::MultiValued :
	// Function base  -- here must be Function::Scalar
	// std::vector < Function::ActionGenerator > actions
	
	std::vector < double > alpha;
	double gamma;
	// upon each action [i], the value v of the function becomes  alpha[i] * v + beta[i]
	// group is commutative only if  beta[i] / (alpha[i]-1.)  does not depend on i
	// we call  gamma  the common value, so  beta[i] = gamma * (alpha[i]-1.)
	// upon action[i], the value v of the function becomes  alpha[i] * ( v - gamma ) + gamma
	// it is the user's responsibility to ensure that each action is invertible
	//   that is,  alpha [i] != 0.
	// and that they commute, that is,  beta[i] / (alpha[i]-1.)  does not depend on i
	//   often, beta [i] are all zero so this is not a problem
	
	inline JumpIsLinear ( const Function::Scalar::MultiValued::JumpIsLinear & ) = delete;
	inline JumpIsLinear ( Function::Scalar::MultiValued::JumpIsLinear && ) = delete;
	
	inline Function::Scalar::MultiValued::JumpIsLinear operator=
		( const Function::Scalar::MultiValued::JumpIsLinear & ) = delete;
	inline Function::Scalar::MultiValued::JumpIsLinear operator=
		( Function::Scalar::MultiValued::JumpIsLinear && ) = delete;

	// size_t nb_of_components ( )
	//   virtual from Function::Core, defined by Function::Scalar, returns 1

	// Function component ( size_t i )  virtual from Function::Core,
	//   defined by Function::Scalar, returns self, never actually used

	// void set_value ( double )
	//   defined by Function::Scalar::MultiValued, delegates to base

	// double get_value_on_cell ( Cell::Core * ) const;
	//   defined by Function::Scalar::MultiValued, delegates to base

	double get_value_on_cell
	( Cell::Core *, const tag::Winding &, const Function::Action & exp ) const;
	//  virtual from Function::Scalar
	
	// double set_value_on_cell ( Cell::Core *, const double & )
	//   defined by Function::Scalar::MultiValued, execution forbidden

	// Function::Jump jump ( )  virtual, defined by Function::Core, execution forbidden
	
	static inline std::pair < std::vector < double >, double > analyse_linear_expression
	( Function expression, Function base );

};  // end of class Function::Scalar::MultiValued::JumpIsLinear

//-----------------------------------------------------------------------------------------//


class Function::Vector::MultiValued : public Function::MultiValued, public Function::Aggregate

// same as Function::Scalar::MultiValued, here with vector values

// abstract class, specialized in Function::Vector::MultiValued::JumpIsSum and JumpIsLinear

// inheriting from Function::Aggregate means simply that there is a 'components' member
							 
{	public :

	// members inherited from Function::MultiValued :
	// Function base  -- here must be Function::Vector
	// std::vector < Function::ActionGenerator > actions

	inline MultiValued ( const tag::AssociatedWith &, const Function & b,
	                     std::vector < Function::ActionGenerator > a     )
	:	Function::MultiValued ( b, a ),
		Function::Aggregate ( tag::reserve_size, b.nb_of_components() )
	{	}
	
	inline MultiValued ( const Function::Vector::MultiValued & ) = delete;
	inline MultiValued ( Function::Vector::MultiValued && ) = delete;
	
	size_t nb_of_components ( ) const;  // virtual from Function::Core, delegates to base

	// Function component ( size_t i )
	//   virtual from Function::Core, defined by Function::Aggregate

	void set_value ( std::vector < double > );
	// virtual from Function::Vector, here delegates to base

	std::vector < double > get_value_on_cell ( Cell::Core * ) const;

	// std::vector < double > get_value_on_cell  with tag::winding
	//   stays pure virtual from Function::Vector

	std::vector < double > set_value_on_cell ( Cell::Core *, const std::vector < double > & );
	// virtual from Function::Vector

	Function deriv ( Function );  // virtual from Function::Core

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core
	
	// Function::Jump jump ( )  virtual, defined by Function::Core, execution forbidden
	// later overridden by Function::Vector::MultiValued::JumpIsSum
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	//  virtual from Function::Core
	#endif

	class JumpIsSum;  class JumpIsLinear;

};  // end of class Function::Vector::MultiValued

//-----------------------------------------------------------------------------------------//


class Function::Vector::MultiValued::JumpIsSum : public Function::Vector::MultiValued

// here the actions are mere translations on RR^n (sums)

{	public :

	// members inherited from Function::MultiValued :
	// Function base  -- here must be Function::Vector
	// std::vector < Function::ActionGenerator > actions
	
	std::vector < std::vector < double > > beta;
	// upon each action[i], the value v of the function becomes  v + beta[i]
	// use a Function::Jump::Sum::Vector instead !

	inline JumpIsSum ( const tag::AssociatedWith &, const Function & f,
	                   std::vector < Function::ActionGenerator > ac,
	                   std::vector < std::vector < double > > be        )
	:	Function::Vector::MultiValued ( tag::associated_with, f, ac ), beta { be }
	{	size_t n_act = ac .size();
		assert ( n_act == be .size() );
		size_t n_comp = f .nb_of_components();
		std::vector < double > beta_scalar ( n_act, 0. );
		for ( size_t j = 0; j < n_comp; j++ )
		{	for ( size_t i = 0; i < n_act; i++ ) beta_scalar [i] = be [i][j];
//		this->components [j] = Function ( tag::whose_core_is,
			this->components.emplace_back ( tag::whose_core_is,
			  new Function::Scalar::MultiValued::JumpIsSum
	                    ( tag::associated_with, f[j], ac, beta_scalar ) );  }  }

	inline JumpIsSum ( const Function::Vector::MultiValued::JumpIsSum & ) = delete;
	inline JumpIsSum ( Function::Vector::MultiValued::JumpIsSum && ) = delete;
	
	inline Function::Vector::MultiValued::JumpIsSum operator=
		( const Function::Vector::MultiValued::JumpIsSum & ) = delete;
	inline Function::Vector::MultiValued::JumpIsSum operator=
		( Function::Vector::MultiValued::JumpIsSum && ) = delete;

	// size_t nb_of_components ( )  virtual from Function::Core
	//   defined by Function::Vector::MultiValued, delegates to base

	// Function component ( size_t i )
	//   virtual from Function::Core, defined by Function::Aggregate

	// void set_value ( std::vector < double > )
	//   defined by Function::Vector::MultiValued, delegates to base

	// std::vector < double > get_value_on_cell ( Cell::Core * ) const;
	//   defined by Function::Vector::MultiValued, delegates to base

	std::vector < double > get_value_on_cell
	( Cell::Core *, const tag::Winding &, const Function::Action & exp ) const;
	//  virtual from Function::Vector
	
	// std::vector < double > set_value_on_cell ( Cell::Core *, const std::vector < double > & )
	//   defined by Function::Vector::MultiValued

	Function::Jump jump ( );
  // virtual, defined by Function::Core, here overridden
	
	inline static std::vector < double > analyse_linear_expression
	( Function expression, Function base );
	
};  // end of class Function::Vector::MultiValued::JumpIsSum

//-----------------------------------------------------------------------------------------//


class Function::Vector::MultiValued::JumpIsLinear : public Function::Vector::MultiValued

// here the actions are linear (affine) maps on RR^n

{	public :

	// members inherited from Function::MultiValued :
	// Function base  -- here must be Function::Vector
	// std::vector < Function::ActionGenerator > actions
	
	std::vector < std::vector < std::vector < double > > > A;
	std::vector < std::vector < double > > b;
	// upon each action [i], the value v of the function becomes  A[i]*v + b[i]
	// it is the user's responsibility to ensure that each action is invertible
	//   (i.e. the matrix A[i] is invertible) and that they commute
	// often, b[i] are all zero so it suffices that the matrices A[i] commute
	//   A[i1] A[i2] == A[i2] A[i1]

	// we keep also the inverse matrices
	std::vector < std::vector < std::vector < double > > > Ainv;

	inline JumpIsLinear ( const tag::AssociatedWith &, const Function & f,
	                      std::vector < Function::ActionGenerator > ac,
	                      std::vector < std::vector < std::vector < double > > > AA,
	                      std::vector < std::vector < double > > bb                 )
	:	Function::Vector::MultiValued ( tag::associated_with, f, ac ), A { AA },  b { bb }
	{	assert ( ac .size() == AA .size() );
		assert ( ac .size() == bb .size() );
		// we need to compute the inverse matrices, we do this only for 2x2 matrices
		std::vector < std::vector < std::vector < double > > > ::const_iterator it;
		for ( it = this->A .begin(); it != this->A .end(); it++ )
		{	const std::vector < std::vector < double > > & mat = *it;
			assert ( mat.size() == 2 );
			double det = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
			this->Ainv .push_back ( { {   mat[1][1] / det, - mat[0][1] / det },
			                          { - mat[1][0] / det,   mat[0][0] / det } } );  }  }
	
	inline JumpIsLinear ( const Function::Vector::MultiValued::JumpIsLinear & ) = delete;
	inline JumpIsLinear ( Function::Vector::MultiValued::JumpIsLinear && ) = delete;
	
	inline Function::Vector::MultiValued::JumpIsLinear operator=
		( const Function::Vector::MultiValued::JumpIsLinear & ) = delete;
	inline Function::Vector::MultiValued::JumpIsLinear operator=
		( Function::Vector::MultiValued::JumpIsLinear && ) = delete;

	// size_t nb_of_components ( )  virtual from Function::Core
	//   defined by Function::Vector::MultiValued, delegates to base

	Function component ( size_t i ) override;   // virtual from Function::Core,
	//  virtual from Function::Core, defined by Function::Aggregate, here overridden
	//  execution forbidden for now, this is a difficult case

	// void set_value ( std::vector < double > )
	//   defined by Function::Vector::MultiValued, delegates to base

	// std::vector < double > get_value_on_cell ( Cell::Core * ) const;
	//   defined by Function::Vector::MultiValued, delegates to base

	std::vector < double > get_value_on_cell
	( Cell::Core *, const tag::Winding &, const Function::Action & exp ) const;
	//  virtual from Function::Vector
	
	// std::vector < double > set_value_on_cell ( Cell::Core *, const std::vector < double > & )
	//   defined by Function::Vector::MultiValued

	// Function::Jump jump ( )  virtual, defined by Function::Core, execution forbidden
	
	inline static std::pair < std::vector < std::vector < double > >, std::vector < double > >
	analyse_linear_expression ( Function expression, Function base );

};  // end of class Function::Vector::MultiValued::JumpIsLinear

//-----------------------------------------------------------------------------------------//


class Function::Jump

// a Function::Jump object describes the behaviour of a multi-function
// when the point goes "around" the manifold a given number of times
// a Cell::Winding describes this number of times

// so a Function::Jump can be applied to a Cell::Winding, producing as result
// a double for a Function::Scalar::MultiValued
// or a vector of doubles for a Function::Vector::MultiValued


{	public :

  std::vector < Function::ActionGenerator > actions;
	
	std::vector < double > ju;
	
	double operator() ( const Function::Action & a ) const;

};  // end of class Function::Jump::Sum::Scalar


//-----------------------------------------------------------------------------------------//


Function::Jump operator+ ( const Function::Jump & j1, const Function::Jump & j2 );

Function::Jump operator+= ( Function::Jump & j1, const Function::Jump & j2 );

Function::Jump operator- ( const Function::Jump & j1, const Function::Jump & j2 );

Function::Jump operator-= ( Function::Jump & j1, const Function::Jump & j2 );

Function::Jump operator* ( double a, const Function::Jump & j );

inline Function::Jump operator* ( const Function::Jump & j, double a )
{	return a*j;  }

Function::Jump operator*= ( Function::Jump & j, double a );

inline Function::Jump operator/ ( const Function::Jump & j, double a )
{	return (1./a)*j;  }

inline Function::Jump operator/= ( Function::Jump & j, double a )
{	return operator*= ( j, 1./a );  }

//-----------------------------------------------------------------------------------------//


inline std::vector < double >  // static
Function::Scalar::MultiValued::JumpIsSum::analyse_linear_expression
( Function expression, Function base )

// 'base' is a variable, say x
// 'expression' is an expression in x (a sum x+c)
// we want to identify c

// the return value should be double
// however, we use std::vector < double > instead
// a zero-length vector means not succeeded
// success is represented by a vector of length one

{	size_t n_coord = base .nb_of_components();
	size_t dim_expr = expression .nb_of_components();
	assert ( n_coord == 1 );
	assert ( dim_expr == 1 );
	Function::CoupledWithField * cf =
		dynamic_cast < Function::CoupledWithField* > ( base .core );
	assert ( cf );
	Field::Double::Scalar * fi =
		dynamic_cast < Field::Double::Scalar * > ( cf->field );
	assert ( fi );
	Function::CoupledWithField * cfe =
		dynamic_cast < Function::CoupledWithField* > ( expression .core );
	if ( cfe )
	{	Field::Double::Scalar * fie =
			dynamic_cast < Field::Double::Scalar * > ( cfe->field );
		if ( fie == nullptr ) return {};
		if ( fi->index_in_heap != fie->index_in_heap ) return {};
		return { 0. };                                             }
	Function::Sum * sum = dynamic_cast < Function::Sum * > ( expression .core );
	if ( sum == nullptr ) return {};
	std::forward_list < Function > terms = sum->terms;
	std::forward_list<Function>::iterator it = terms .begin();
	assert ( it != terms .end() );
	Function x = *it;
	if ( x.core != base .core )
	{	Function c = x;
		Function::Constant * cc = dynamic_cast < Function::Constant * > ( c.core );
		if ( cc == nullptr ) return {};
		it++;
		assert ( it != terms .end() );
		Function xx = *it;
		if ( xx.core != base .core ) return {};
		return { cc->value };                                                           }
	it++;
	assert ( it != terms .end() );
	Function c = *it;
	it++;
	if ( it != terms .end() ) return {};
	Function::Constant * cc = dynamic_cast < Function::Constant * > ( c .core );
	if ( cc == nullptr ) return {};
	return { cc->value };                                                                }

//-----------------------------------------------------------------------------------------//


inline std::vector < double >  // static
Function::Vector::MultiValued::JumpIsSum::analyse_linear_expression
( Function expression, Function base )

// 'base' is a set of variables, say xyz
// 'expression' is a vector expression in xyz (a sum : x+a, y+b, z+c)
// we want to identify a, b, c
// a zero-length return vector means not succeeded

{	size_t n_coord = base .nb_of_components();
	size_t dim_expr = expression .nb_of_components();
	assert ( n_coord == dim_expr );
	std::vector < double > res;
	res.reserve ( n_coord );
	for ( size_t i = 0; i < n_coord; i++ )
	{	std::vector < double > res_i =
			Function::Scalar::MultiValued::JumpIsSum::analyse_linear_expression
			( expression [i], base [i] );
		if ( res_i .size() == 0 ) return {};
		assert ( res_i .size() == 1 );
		res .push_back ( res_i [0] );                                          }
	return res;                                                                }

//-----------------------------------------------------------------------------------------//


inline std::pair < std::vector < double >, double >  // static
Function::Scalar::MultiValued::JumpIsLinear::analyse_linear_expression
( Function expression, Function base )

// 'base' is a set of variables, say X
// 'expression' is a scalar expression in X  ( a1 x1 + a2 x2 + ... + b )
// we want to identify ai and b

{	size_t n_coord = base .nb_of_components();
	size_t dim_expr = expression .nb_of_components();
	assert ( dim_expr == 1 );
	std::pair < std::vector < double >, double > res
		{ std::vector < double > ( n_coord, 0. ), 0. };
	// often, 'expression' will be a sum
	Function::Sum * expr_sum = dynamic_cast < Function::Sum* > ( expression .core );
	if ( expr_sum )
	{	std::forward_list < Function > tl = expr_sum->terms;
		for ( std::forward_list < Function > ::const_iterator
					it = tl .begin(); it != tl .end(); it++           )
		{	Function term = *it;
			// often, 'term' will be a product
			Function::Product * term_prod = dynamic_cast < Function::Product* > ( term.core );
			if ( term_prod )  // we assume this product has the form  constant * variable
			{	std::forward_list < Function > fl = term_prod->factors;
				std::forward_list < Function > ::const_iterator itt = fl .begin();
				assert ( itt != fl .end() );
				Function c = *itt;
				Function::Constant * cc =
					tag::Util::assert_cast < Function::Core*, Function::Constant* > ( c .core );
				itt++;
				assert ( itt != fl .end() );
				Function var = *itt;
				itt++;
				assert ( itt == fl .end() );
				Function::CoupledWithField::Scalar * var_cf =
					tag::Util::assert_cast < Function::Core*,
					                         Function::CoupledWithField::Scalar* > ( var .core );
				Field::Double::Scalar * var_field =
					tag::Util::assert_cast < Field::Core*, Field::Double::Scalar * > ( var_cf->field );
				bool found_var = false;
				size_t i = 0;
				for ( ; i < n_coord; i++ )
				{	Function::CoupledWithField::Scalar * base_i_cf =
						tag::Util::assert_cast < Function::Core*, Function::CoupledWithField::Scalar* >
						( base [i] .core );
					Field::Double::Scalar * base_i_field =
						tag::Util::assert_cast < Field::Core*, Field::Double::Scalar * >
						( base_i_cf->field );
					if ( base_i_field->index_in_heap == var_field->index_in_heap )
					{	found_var = true;  break;  }                                                   }
				assert ( found_var );
				assert ( res .first [i] == 0. );
				res .first [i] = cc->value;                                                          }
			else  // term is not a product
			{	// could a variable xi or the constant b
				Function::CoupledWithField * term_var =
					dynamic_cast < Function::CoupledWithField* > ( term .core );
				if ( term_var )
				{	Field::Double::Scalar * var_field =
						tag::Util::assert_cast < Field::Core*, Field::Double::Scalar * >
						( term_var->field );
					bool found_var = false;
					size_t i = 0;
					for ( ; i < n_coord; i++ )
					{	Function::CoupledWithField::Scalar * base_i_cf =
							tag::Util::assert_cast < Function::Core*, Function::CoupledWithField::Scalar* >
							( base [i] .core );
						Field::Double::Scalar * base_i_field =
							tag::Util::assert_cast < Field::Core*, Field::Double::Scalar * >
							( base_i_cf->field );
						if ( base_i_field->index_in_heap == var_field->index_in_heap )
							{	found_var = true;  break;  }                                            }
					assert ( found_var );
					assert ( res .first [i] == 0. );
					res .first [i] = 1.;                                                               }
				else  // term must be the constant b
				{	Function::Constant * cc =
						tag::Util::assert_cast < Function::Core*, Function::Constant* > ( term .core );
					assert ( res .second == 0. );
					res .second = cc->value;
				}  }  }  }
	else  // expression is not a sum
	{	// could be a product
		Function::Product * expr_prod = dynamic_cast < Function::Product* > ( expression .core );
		if ( expr_prod )  // we assume this product has the form  constant * variable
		{	std::forward_list < Function > fl = expr_prod->factors;
			std::forward_list < Function > ::const_iterator itt = fl .begin();
			assert ( itt != fl .end() );
			Function c = *itt;
			Function::Constant * cc =
				tag::Util::assert_cast < Function::Core*, Function::Constant* > ( c .core );
			itt++;
			assert ( itt != fl .end() );
			Function var = *itt;
			itt++;
			assert ( itt == fl .end() );
			Function::CoupledWithField::Scalar * var_cf =
				tag::Util::assert_cast < Function::Core*,
				                         Function::CoupledWithField::Scalar* > ( var .core );
			Field::Double::Scalar * var_field =
				tag::Util::assert_cast < Field::Core*, Field::Double::Scalar * > ( var_cf->field );
			bool found_var = false;
			size_t i = 0;
			for ( ; i < n_coord; i++ )
			{	Function::CoupledWithField::Scalar * base_i_cf =
					tag::Util::assert_cast < Function::Core*, Function::CoupledWithField::Scalar* >
					( base [i] .core );
				Field::Double::Scalar * base_i_field =
					tag::Util::assert_cast < Field::Core*, Field::Double::Scalar * >
					( base_i_cf->field );
				if ( base_i_field->index_in_heap == var_field->index_in_heap )
				{	found_var = true;  break;  }                                            }
			assert ( found_var );
			assert ( res.first[i] == 0. );
			res .first [i] = cc->value;                                                       }
		else  // expression is not a sum, not a product
		{	// must be a variable xi
			Function::CoupledWithField::Scalar * expr_var =
				tag::Util::assert_cast < Function::Core*, Function::CoupledWithField::Scalar* >
				( expression .core );
			Field::Double::Scalar * var_field =
				tag::Util::assert_cast < Field::Core*, Field::Double::Scalar * >
				( expr_var->field );
			bool found_var = false;
			size_t i = 0;
			for ( ; i < n_coord; i++ )
			{	Function::CoupledWithField::Scalar * base_i_cf =
					tag::Util::assert_cast < Function::Core*, Function::CoupledWithField::Scalar* >
					( base [i] .core );
				Field::Double::Scalar * base_i_field =
					tag::Util::assert_cast < Field::Core*, Field::Double::Scalar * >
					( base_i_cf->field );
				if ( base_i_field->index_in_heap == var_field->index_in_heap )
				{	found_var = true;  break;  }                                            }
			assert ( found_var );
			assert ( res .first [i] == 0. );
			res .first [i] = 1.;                                                               }  }
	return res;
	
}  // end of Function::Scalar::MultiValued::JumpIsLinear::analyse_linear_expression

//-----------------------------------------------------------------------------------------//


inline std::pair < std::vector < std::vector < double > >, std::vector < double > >
Function::Vector::MultiValued::JumpIsLinear::analyse_linear_expression  // static
( Function expression, Function base )

// 'base' is a set of variables, say X
// 'expression' is a vector expression in X :  A X + b
// we want to identify A and b

{	size_t n_coord = base .nb_of_components();
	size_t dim_expr = expression .nb_of_components();
	assert ( n_coord == dim_expr );
	std::pair < std::vector < std::vector < double > >, std::vector < double > > res;
	for ( size_t i = 0; i < n_coord; i++ )
		{	std::pair < std::vector < double >, double > res_i =
			Function::Scalar::MultiValued::JumpIsLinear::analyse_linear_expression
			( expression [i], base );
		assert ( res_i .first .size() == n_coord );
		res .first .push_back ( std::move ( res_i .first ) );
		res .second .push_back ( res_i.second );                                   }
	return res;                                                                    }

//-----------------------------------------------------------------------------------------//


inline Function::Jump Function::jump ( )
{	return this->core->jump();  }

//-----------------------------------------------------------------------------------------//
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
	
	inline Function::TakenOnCell & operator= ( const Function::TakenOnCell & other )
	{	if ( this->f->nb_of_components() == 1 )
		{	(*this) = static_cast < double > (other);
			return *this;                            }
		else
		{	(*this) = static_cast < std::vector < double > > (other);
			return *this;                                             }  }

	inline Function::TakenOnCell & operator= ( const Function::TakenOnCell && other )
	{	if ( this->f->nb_of_components() == 1 )
		{	(*this) = static_cast < double > (other);
			return *this;                            }
		else
		{	(*this) = static_cast < std::vector < double > > (other);
			return *this;                                             }  }
	
	inline operator double() const
	// can be used like in  double x = f ( cll )  or  cout << f ( cll )
	{	Function::Scalar * f_scalar =
			tag::Util::assert_cast < Function::Core*, Function::Scalar* > ( this->f );
		return f_scalar->get_value_on_cell ( this->cll );                            }
		
	inline double operator= ( const double & x )
	// can be used like in  f ( cll ) = 2.0
	{	Function::Scalar * f_scalar =
			tag::Util::assert_cast < Function::Core*, Function::Scalar* > ( this->f );
		return f_scalar->set_value_on_cell ( this->cll, x );                         }

	inline double operator+= ( const double & x )
	// can be used like in  f ( cll ) += 2.0
	{	Function::Scalar * f_scalar =
			tag::Util::assert_cast < Function::Core*, Function::Scalar* > ( this->f );
		return f_scalar->set_value_on_cell
			( this->cll, f_scalar->get_value_on_cell ( this->cll ) + x );              }

	inline double operator-= ( const double & x )
	// can be used like in  f ( cll ) -= 2.0
	{	Function::Scalar * f_scalar =
			tag::Util::assert_cast < Function::Core*, Function::Scalar* > ( this->f );
		return f_scalar->set_value_on_cell
			( this->cll, f_scalar->get_value_on_cell ( this->cll ) - x );              }

	inline double operator*= ( const double & x )
	// can be used like in  f ( cll ) *= 2.0
	{	Function::Scalar * f_scalar =
			tag::Util::assert_cast < Function::Core*, Function::Scalar* > ( this->f );
		return f_scalar->set_value_on_cell
			( this->cll, f_scalar->get_value_on_cell ( this->cll ) * x );              }

	inline double operator/= ( const double & x )
	// can be used like in  f ( cll ) /= 2.0
	{	Function::Scalar * f_scalar =
			tag::Util::assert_cast < Function::Core*, Function::Scalar* > ( this->f );
		return f_scalar->set_value_on_cell
			( this->cll, f_scalar->get_value_on_cell ( this->cll ) / x );              }

	inline operator std::vector<double>() const
	// can be used like in  vector<double> vec = f(cll)
	{	Function::Vector * f_vect =
			tag::Util::assert_cast < Function::Core*, Function::Vector* > ( this->f );
	  return f_vect->get_value_on_cell ( this->cll );                              }
		
	inline std::vector<double> operator= ( const std::vector < double > & x )
	// can be used like in  f ( cll ) = vec
	{	Function::Vector * f_vect =
			tag::Util::assert_cast < Function::Core*, Function::Vector* > ( this->f );
		return f_vect->set_value_on_cell ( this->cll, x );                           }

};  // end of class Function::TakenOnCell	

//-----------------------------------------------------------------------------------------//


class Function::TakenOnWindingCell

{	public :

  Function::Core * f;
	const Function::Action & winding;
	
	// Function::TakenOnWindingCell should only be used as temporary objects
	// they should be immediately converted to a (reference to a) double or vector<double>
	// so an ordinary pointer is OK here
	
	Cell::Core * const cll;

	inline TakenOnWindingCell ( const Function & ff, const Cell & c,
                               const Function::Action & exp )
	:	f { ff .core }, winding { exp }, cll { c .core }
	{	}

	TakenOnWindingCell ( const Function::TakenOnWindingCell & other ) = delete;
	TakenOnWindingCell ( const Function::TakenOnWindingCell && other ) = delete;
	
	inline Function::TakenOnWindingCell & operator= ( const Function::TakenOnWindingCell & other )
	{	if ( this->f->nb_of_components() == 1 )
		{	(*this) = static_cast < double > ( other );
			return *this;                            }
		else
		{	(*this) = static_cast < std::vector < double > > ( other );
			return *this;                                          }  }

	inline Function::TakenOnWindingCell & operator= ( const Function::TakenOnWindingCell && other )
	{	if ( this->f->nb_of_components() == 1 )
		{	(*this) = static_cast < double > ( other );
			return *this;                             }
		else
		{	(*this) = static_cast < std::vector < double > > ( other );
			return *this;                                          }  }
	
	inline operator double() const
	// can be used like in  double x = f ( cll )  or  cout << f ( cll )
	{	Function::Scalar * f_scalar =
			tag::Util::assert_cast < Function::Core*, Function::Scalar* > ( this->f );
		return f_scalar->get_value_on_cell ( this->cll, tag::winding, this->winding );     }
		
	inline double operator= ( const double & x )
	// can be used like in  f ( cll ) = 2.0
	{	Function::Scalar * f_scalar =
			tag::Util::assert_cast < Function::Core*, Function::Scalar* > ( this->f );
		return f_scalar->set_value_on_cell ( this->cll, x );                         }

	inline double operator+= ( const double & x )
	// can be used like in  f ( cll ) += 2.0
	{	Function::Scalar * f_scalar =
			tag::Util::assert_cast < Function::Core*, Function::Scalar* > ( this->f );
		return f_scalar->set_value_on_cell
			( this->cll, f_scalar->get_value_on_cell ( this->cll ) + x );              }

	inline double operator-= ( const double & x )
	// can be used like in  f ( cll ) -= 2.0
	{	Function::Scalar * f_scalar =
			tag::Util::assert_cast < Function::Core*, Function::Scalar* > ( this->f );
		return f_scalar->set_value_on_cell
			( this->cll, f_scalar->get_value_on_cell ( this->cll ) - x );              }

	inline double operator*= ( const double & x )
	// can be used like in  f ( cll ) *= 2.0
	{	Function::Scalar * f_scalar =
			tag::Util::assert_cast < Function::Core*, Function::Scalar* > ( this->f );
		return f_scalar->set_value_on_cell
			( this->cll, f_scalar->get_value_on_cell ( this->cll ) * x );              }

	inline double operator/= ( const double & x )
	// can be used like in  f ( cll ) /= 2.0
	{	Function::Scalar * f_scalar =
			tag::Util::assert_cast < Function::Core*, Function::Scalar* > ( this->f );
		return f_scalar->set_value_on_cell
			( this->cll, f_scalar->get_value_on_cell ( this->cll ) / x );              }

	inline operator std::vector<double>() const
	// can be used like in  vector<double> vec = f ( cll )
	{	Function::Vector * f_vect =
			tag::Util::assert_cast < Function::Core*, Function::Vector* > ( this->f );
	  return f_vect->get_value_on_cell ( this->cll, tag::winding, this->winding );  }
		
	inline std::vector<double> operator= ( const std::vector<double> & x )
	// can be used like in  f ( cll ) = vec
	{	Function::Vector * f_vect =
			tag::Util::assert_cast < Function::Core*, Function::Vector* > ( this->f );
		return f_vect->set_value_on_cell ( this->cll, x );                           }

};  // end of class Function::TakenOnWindingCell

//-----------------------------------------------------------------------------------------//

class Function::Equality

// the result of an equality comparison between Functions

{	public :

	Function lhs, rhs;
	
};

//-----------------------------------------------------------------------------------------//
//  classes below allow for boolean epressions like
//  f < g < h  --  equivalent to  ( f < g ) < h
//  and   f > g > h  --  equivalent to  ( f > g ) > h
//-----------------------------------------------------------------------------------------//

class Function::Inequality::LessThanZero

// the result of an inequality comparison between a Function and zero

{	public :

	Function expr;

	inline operator Function::Inequality::Set() const;
	
};

//-----------------------------------------------------------------------------------------//

class Function::Inequality::LessThan

// the result of an inequality comparison between Functions

{	public :

	Function low, high;
	
	inline operator Function::Inequality::LessThanZero() const;
	inline operator Function::Inequality::Set() const;

};

//-----------------------------------------------------------------------------------------//

class Function::Inequality::GreaterThan

// the result of an inequality comparison between Functions

{	public :

	Function low, high;
	
	inline operator Function::Inequality::LessThanZero() const;
	inline operator Function::Inequality::Set() const;

};

//-----------------------------------------------------------------------------------------//

class tag::Util::InequalitySet

// a vector of inequalities

{	public :

	std::vector < Function::Inequality::LessThanZero > vec;

	inline bool on_cell ( const Cell & ) const;
	inline bool on_cell
	( const Cell &, const tag::Winding &, const Function::Action & ) const;
	
};

//-----------------------------------------------------------------------------------------//

inline Function::Inequality::LessThan::operator Function::Inequality::LessThanZero() const
{	return Function::Inequality::LessThanZero { this->low - this->high };  }

inline Function::Inequality::LessThan::operator Function::Inequality::Set() const
{	return Function::Inequality::Set
		{ { Function::Inequality::LessThanZero { this->low - this->high } } };  }

inline Function::Inequality::GreaterThan::operator Function::Inequality::LessThanZero() const
{	return Function::Inequality::LessThanZero { this->low - this->high };  }

inline Function::Inequality::GreaterThan::operator Function::Inequality::Set() const
{	return Function::Inequality::Set
		{ { Function::Inequality::LessThanZero { this->low - this->high } } };  }

inline Function::Inequality::LessThanZero::operator Function::Inequality::Set() const
{	return Function::Inequality::Set { { *this } };  }


inline Function::Inequality::Set operator&&
( const Function::Inequality::Set & f, const Function::Inequality::Set & g )

{	Function::Inequality::Set res;
	for ( std::vector < Function::Inequality::LessThanZero > ::const_iterator
        it = f .vec .begin(); it != f .vec .end(); it++                     )
		res .vec .push_back ( *it );
	for ( std::vector < Function::Inequality::LessThanZero > ::const_iterator
        it = g .vec .begin(); it != g .vec .end(); it++                     )
		res .vec .push_back ( *it );
	return res;                                                                   }
// optimize !!


inline bool Function::Inequality::Set::on_cell ( const Cell & cll ) const

{	for ( std::vector < Function::Inequality::LessThanZero > ::const_iterator
        it = this->vec .begin(); it != this->vec .end(); it++               )
		if ( ( it->expr ) ( cll ) > 0. ) return false;
	return true;                                                               }


inline bool Function::Inequality::Set::on_cell
( const Cell & cll, const tag::Winding &, const Function::Action & a ) const

{	for ( std::vector < Function::Inequality::LessThanZero > ::const_iterator
        it = this->vec .begin(); it != this->vec .end(); it++               )
		if ( ( it->expr ) ( cll, tag::winding, a ) > 0. ) return false;
	return true;                                                             }


inline Function::Equality operator== ( const Function & f, const Function & g )
{	return Function::Equality { f, g };  }


inline Function::Inequality::LessThan operator<= ( const Function & f, const Function & g )
{	return Function::Inequality::LessThan { f, g };  }

inline Function::Inequality::LessThan operator< ( const Function & f, const Function & g )
{	return Function::Inequality::LessThan { f, g };  }

inline Function::Inequality::GreaterThan operator>= ( const Function & f, const Function & g )
{	return Function::Inequality::GreaterThan { g, f };  }

inline Function::Inequality::GreaterThan operator> ( const Function & f, const Function & g )
{	return Function::Inequality::GreaterThan { g, f };  }

inline Function::Inequality::Set operator<=   //  f <= g <= h
( const Function::Inequality::LessThan & fg, const Function & h )
{	Function::Inequality::Set res;
	res.vec.push_back ( fg );
	res.vec.push_back ( fg .high < h );
	return res;                         }

inline Function::Inequality::Set operator<   //  f < g < h
( const Function::Inequality::LessThan & fg, const Function & h )
{	Function::Inequality::Set res;
	res.vec.push_back ( fg );
	res.vec.push_back ( fg .high < h );
	return res;                         }

inline Function::Inequality::Set operator>=   //  f >= g >= h
( const Function::Inequality::GreaterThan & fg, const Function & h )
{	Function::Inequality::Set res;
	res.vec.push_back ( fg );
	res.vec.push_back ( fg .low > h );
	return res;                        }

inline Function::Inequality::Set operator>   //  f > g > h
( const Function::Inequality::GreaterThan & fg, const Function & h )
{	Function::Inequality::Set res;
	res.vec.push_back ( fg );
	res.vec.push_back ( fg .low > h );
	return res;                        }


inline size_t Function::nb_of_components ( ) const
{	return this->core->nb_of_components();  }


inline Function Function::operator[] ( size_t i ) const
{	assert ( i < this->nb_of_components() );
	Function::Scalar * f_scalar = dynamic_cast < Function::Scalar * > ( this->core );
	if ( f_scalar )
	{	assert ( this->nb_of_components() == 1 );
		assert ( i == 0 );
		return *this;                             }
  Function::Vector * f_vector = dynamic_cast < Function::Vector * > ( this->core );
	assert ( f_vector );
	return Function ( tag::whose_core_is, f_vector->component(i).core );               }


inline Function::TakenOnCell Function::operator() ( const Cell & cll ) const
{	return Function::TakenOnCell { this->core, cll.core };   }


inline Function::TakenOnWindingCell Function::operator()
( const Cell & cll, const tag::Winding &, const Function::Action & exp ) const
{	return Function::TakenOnWindingCell ( *this, cll, exp );   }

//----------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------//


inline void Mesh::draw_ps ( std::string file_name,
         const tag::Unfold &, const std::vector < Function::Action > & v,
         const tag::OverRegion &, const tag::Util::InequalitySet & c1,
                                  const tag::Util::InequalitySet & c2    ) const
{	this->draw_ps ( file_name, tag::unfold, v, tag::over_region, c1 && c2 );  }

//----------------------------------------------------------------------------------//
														
inline void Mesh::draw_ps ( std::string file_name,
         const tag::Unfold &, const std::vector < Function::Action > & v,
         const tag::OverRegion &, const tag::Util::InequalitySet & c1,
                                  const tag::Util::InequalitySet & c2,
                                  const tag::Util::InequalitySet & c3     ) const
{	this->draw_ps ( file_name, tag::unfold, v, tag::over_region, c1 && c2 && c3 );  }

//----------------------------------------------------------------------------------//
														
inline void Mesh::draw_ps ( std::string file_name,
         const tag::Unfold &, const std::vector < Function::Action > & v,
         const tag::OverRegion &, const tag::Util::InequalitySet & c1,
                                  const tag::Util::InequalitySet & c2,
                                  const tag::Util::InequalitySet & c3,
                                  const tag::Util::InequalitySet & c4     ) const
{	this->draw_ps ( file_name, tag::unfold, v, tag::over_region, c1 && c2 && c3 && c4 );  }

														
inline void Mesh::draw_ps ( std::string file_name,
                            const tag::Unfold &, const tag::OverRegion &,
                            const tag::Util::InequalitySet & c1,
                            const tag::Util::InequalitySet & c2           ) const
{	this->draw_ps ( file_name, tag::unfold, tag::over_region, c1 && c2 );  }


inline void Mesh::draw_ps ( std::string file_name,
                            const tag::Unfold &, const tag::OverRegion &,
                            const tag::Util::InequalitySet & c1,
                            const tag::Util::InequalitySet & c2,
                            const tag::Util::InequalitySet & c3           ) const
{	this->draw_ps ( file_name, tag::unfold, tag::over_region, c1 && c2 && c3 );  }


inline void Mesh::draw_ps ( std::string file_name,
                            const tag::Unfold &, const tag::OverRegion &,
                            const tag::Util::InequalitySet & c1,
                            const tag::Util::InequalitySet & c2,
                            const tag::Util::InequalitySet & c3,
                            const tag::Util::InequalitySet & c4           ) const
{	this->draw_ps ( file_name, tag::unfold, tag::over_region, c1 && c2 && c3 && c4 );  }

//----------------------------------------------------------------------------------//


inline Mesh Mesh::unfold ( const tag::OverRegion &,
                           const tag::Util::InequalitySet & c1,
                           const tag::Util::InequalitySet & c2 ) const
{	return this->unfold ( tag::over_region, c1 && c2 );  }

inline Mesh Mesh::unfold ( const tag::OverRegion &,
                           const tag::Util::InequalitySet & c1,
                           const tag::Util::InequalitySet & c2,
                           const tag::Util::InequalitySet & c3 ) const
{	return this->unfold ( tag::over_region, c1 && c2 && c3 );  }

inline Mesh Mesh::unfold ( const tag::OverRegion &,
                           const tag::Util::InequalitySet & c1,
                           const tag::Util::InequalitySet & c2,
                           const tag::Util::InequalitySet & c3,
                           const tag::Util::InequalitySet & c4 ) const
{	return this->unfold ( tag::over_region, c1 && c2 && c3 && c4 );  }


inline Mesh Mesh::unfold ( const tag::OverRegion &,
                           const tag::Util::InequalitySet & c1,
                           const tag::Util::InequalitySet & c2,
         const tag::ReturnMapBetween &, const tag::CellsOfDim &,
         size_t dim, std::map < Cell, std::pair < Cell, tag::Util::Action > > & mapping ) const
{	return this->unfold ( tag::over_region, c1 && c2,
	                      tag::return_map_between, tag::cells_of_dim, dim, mapping );  }

inline Mesh Mesh::unfold ( const tag::OverRegion &,
                           const tag::Util::InequalitySet & c1,
                           const tag::Util::InequalitySet & c2,
                           const tag::Util::InequalitySet & c3,
         const tag::ReturnMapBetween &, const tag::CellsOfDim &,
         size_t dim, std::map < Cell, std::pair < Cell, tag::Util::Action > > & mapping ) const
{	return this->unfold ( tag::over_region, c1 && c2 && c3,
	                      tag::return_map_between, tag::cells_of_dim, dim, mapping );  }

inline Mesh Mesh::unfold ( const tag::OverRegion &,
                           const tag::Util::InequalitySet & c1,
                           const tag::Util::InequalitySet & c2,
                           const tag::Util::InequalitySet & c3,
                           const tag::Util::InequalitySet & c4,
         const tag::ReturnMapBetween &, const tag::CellsOfDim &,
         size_t dim, std::map < Cell, std::pair < Cell, tag::Util::Action > > & mapping ) const
{	return this->unfold ( tag::over_region, c1 && c2 && c3 && c4,
	                      tag::return_map_between, tag::cells_of_dim, dim, mapping );  }

inline Mesh Mesh::unfold ( const std::vector < tag::Util::Action > & aa, const tag::OverRegion &,
                           const tag::Util::InequalitySet & c1,
                           const tag::Util::InequalitySet & c2           ) const
{	return this->unfold ( aa, tag::over_region, c1 && c2 );  }

inline Mesh Mesh::unfold ( const std::vector < tag::Util::Action > & aa, const tag::OverRegion &,
                           const tag::Util::InequalitySet & c1,
                           const tag::Util::InequalitySet & c2,
                           const tag::Util::InequalitySet & c3           ) const
{	return this->unfold ( aa, tag::over_region, c1 && c2 && c3 );  }

inline Mesh Mesh::unfold ( const std::vector < tag::Util::Action > & aa, const tag::OverRegion &,
                           const tag::Util::InequalitySet & c1,
                           const tag::Util::InequalitySet & c2,
                           const tag::Util::InequalitySet & c3,
                           const tag::Util::InequalitySet & c4           ) const
{	return this->unfold ( aa, tag::over_region, c1 && c2 && c3 && c4 );  }


inline Mesh Mesh::unfold ( const std::vector < tag::Util::Action > & aa, const tag::OverRegion &,
                           const tag::Util::InequalitySet & c1,
                           const tag::Util::InequalitySet & c2,
         const tag::ReturnMapBetween &, const tag::CellsOfDim &,
         size_t dim, std::map < Cell, std::pair < Cell, tag::Util::Action > > & mapping ) const
{	return this->unfold ( aa, tag::over_region, c1 && c2,
	                      tag::return_map_between, tag::cells_of_dim, dim, mapping );  }

inline Mesh Mesh::unfold ( const std::vector < tag::Util::Action > & aa, const tag::OverRegion &,
                           const tag::Util::InequalitySet & c1,
                           const tag::Util::InequalitySet & c2,
                           const tag::Util::InequalitySet & c3,
         const tag::ReturnMapBetween &, const tag::CellsOfDim &,
         size_t dim, std::map < Cell, std::pair < Cell, tag::Util::Action > > & mapping ) const
{	return this->unfold ( aa, tag::over_region, c1 && c2 && c3,
	                      tag::return_map_between, tag::cells_of_dim, dim, mapping );  }

inline Mesh Mesh::unfold ( const std::vector < tag::Util::Action > & aa, const tag::OverRegion &,
                           const tag::Util::InequalitySet & c1,
                           const tag::Util::InequalitySet & c2,
                           const tag::Util::InequalitySet & c3,
                           const tag::Util::InequalitySet & c4,
         const tag::ReturnMapBetween &, const tag::CellsOfDim &,
         size_t dim, std::map < Cell, std::pair < Cell, tag::Util::Action > > & mapping ) const
{	return this->unfold ( aa, tag::over_region, c1 && c2 && c3 && c4,
	                      tag::return_map_between, tag::cells_of_dim, dim, mapping );  }


}  // namespace maniFEM


#endif
// ifndef MANIFEM_FUNCTION_H
