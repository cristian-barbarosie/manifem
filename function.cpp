
// maniFEM function.cpp 2019.12.28

//    This file is part of maniFEM, a C++ library for meshes on manifolds and finite elements.

//    ManiFEM is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Lesser General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.

//    ManiFEM is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.

//    You should have received a copy of the GNU Lesser General Public License
//    along with maniFEM.  If not, see <https://www.gnu.org/licenses/>.

//    Copyright 2019, 2020 Cristian Barbarosie cristian.barbarosie@gmail.com
//    https://github.com/cristian-barbarosie/manifem

#include "math.h"
#include <sstream>

#include "function.h"

using namespace maniFEM;

size_t Function::Scalar::nb_of_components ( ) const  // virtual from Function::Core
{	return 1;  }

size_t Function::Aggregate::nb_of_components ( ) const
// virtual from Function::Core, through Function::Vector
{	return this->components.size();  }

size_t Function::CoupledWithField::Vector::nb_of_components ( ) const
// virtual from Function::Core through Function::Vector
{	return this->field->nb_of_components();  }


Function Function::Aggregate::component ( size_t i )
// virtual from Function::Vector
{	assert ( i < this->components.size() );
	return components[i];                    }

Function Function::CoupledWithField::Vector::component ( size_t i )
// virtual from Function::Vector
{	size_t n = this->field->nb_of_components();
	assert ( i < n );
	assert ( this->components.size() == n );
	return this->components[i];                   }

	
double Function::ArithmeticExpression::set_value_on_cell
( Cell::Core * cll, const double & )  // virtual from Function::Scalar
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
            << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot assign to an arithmetic expression." << std::endl;
	exit ( 1 );                                                              }
	

double Function::Constant::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Scalar, through Function::ArithmeticExpression
{ return this->val;  }

double Function::Step::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Scalar, through Function::ArithmeticExpression
{	Function::Scalar * arg_scalar = Function::core_to_scalar ( this->arg.core.get() );
	double arg_v = arg_scalar->get_value_on_cell(cll);
	for ( size_t i = 0; i < this->cuts.size(); i++ )
		if ( arg_v < cuts[i] )
		{	Function::Scalar * val_i_scalar = Function::core_to_scalar ( this->values[i].core.get() );
			return val_i_scalar->get_value_on_cell(cll);                                                }
	Function::Scalar * val_scalar = Function::core_to_scalar ( this->values.back().core.get() );
	return val_scalar->get_value_on_cell(cll);                                                        }
	
double Function::Sum::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Scalar, through Function::ArithmeticExpression
{ std::forward_list<Function>::const_iterator it = this->terms.begin();
	double sum = 0.;
	for ( ; it != this->terms.end(); it++ )
	{	Function::Scalar * term_scalar = Function::core_to_scalar ( it->core.get() );
		sum += term_scalar->get_value_on_cell ( cll );                               }
	return sum;                                                                      }

double Function::Product::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Scalar, through Function::ArithmeticExpression
{ std::forward_list<Function>::const_iterator it = this->factors.begin();
	double prod = 1.;
	for ( ; it != this->factors.end(); it++ )
	{	Function::Scalar * fact_scalar = Function::core_to_scalar ( it->core.get() );
		prod *= fact_scalar->get_value_on_cell ( cll );                              }
	return prod;                                                                     }
	
double Function::Power::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Scalar, through Function::ArithmeticExpression
{	Function::Scalar * base_scalar = Function::core_to_scalar ( this->base.core.get() );
	return std::pow ( base_scalar->get_value_on_cell(cll), this->exponent );              }
	
double Function::Sin::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Scalar, through Function::ArithmeticExpression
{	Function::Scalar * base_scalar = Function::core_to_scalar ( this->base.core.get() );
	return std::sin ( base_scalar->get_value_on_cell(cll) );                              }
	
double Function::Cos::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Scalar, through Function::ArithmeticExpression
{	Function::Scalar * base_scalar = Function::core_to_scalar ( this->base.core.get() );
	return std::cos ( base_scalar->get_value_on_cell(cll) );                              }
	
double Function::JustTesting::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Scalar, through Function::ArithmeticExpression
{ assert ( false );   }
	
std::vector<double> Function::Aggregate::set_value_on_cell
( Cell::Core * cll, const std::vector<double> & x )  // virtual from Function::Vector
{ size_t n = this->nb_of_components();
	assert ( n == x.size() );
	for ( size_t i = 0; i < n; i++ )
	{	Function::Scalar * comp_scalar =
			Function::core_to_scalar ( this->components[i].core.get() );
		comp_scalar->set_value_on_cell ( cll, x[i] );                  }
	return x;                                                          }

std::vector<double> Function::Aggregate::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Vector
{ size_t n = this->nb_of_components();
	std::vector<double> result ( n );
	for ( size_t i = 0; i < n; i++ )
	{	Function::Scalar * comp_scalar =
			Function::core_to_scalar ( this->components[i].core.get() );
		result[i] = comp_scalar->get_value_on_cell ( cll );            }
	return result;                                                      }

double Function::CoupledWithField::Scalar::set_value_on_cell
( Cell::Core * cll, const double & x )  // virtual from Function::Scalar
{ return this->field->on_cell(cll) = x;  }

double Function::CoupledWithField::Scalar::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Scalar
{ return this->field->on_cell(cll);  }
// { return this->field->on_cell(cll).reference();  }
	
std::vector<double> Function::CoupledWithField::Vector::set_value_on_cell
( Cell::Core * cll, const std::vector<double> & x )  // virtual from Function::Vector
{ return this->field->on_cell(cll) = x;  }

std::vector<double> Function::CoupledWithField::Vector::get_value_on_cell
( Cell::Core * cll ) const  // virtual from Function::Vector
{ return this->field->on_cell(cll);  }
	

#ifndef NDEBUG	

std::string Function::Constant::repr ( const Function::From & from ) const
{	std::stringstream ss;
	ss << this->val; 
	std::string s = ss.str();
	if ( this->val < 0. ) s = "(" + s + ")";
	return s;                         }

std::string Function::Sum::repr ( const Function::From & from ) const
{ std::forward_list<Function>::const_iterator it = this->terms.begin();
	assert ( it != terms.end() );
	std::string s = it->core->repr ( Function::from_sum );
	for ( it++; it != terms.end(); it++ )
		s = s + '+' + it->core->repr ( Function::from_sum );
	if ( ( from == Function::from_power ) or ( from == Function::from_product ) )
		s = "(" + s + ")";
	return s;                                                                      }

std::string Function::Product::repr ( const Function::From & from ) const
{ std::forward_list<Function>::const_iterator it = this->factors.begin();
	assert ( it != factors.end() );
	std::string s = it->core->repr ( Function::from_product );
	for ( it++; it != factors.end(); it++ )
		s = s + '*' + it->core->repr ( Function::from_product );
	if ( from == Function::from_power ) s = "(" + s + ")";
	return s;                                                    }

std::string Function::Power::repr ( const Function::From & from ) const
{	std::string s = this->base.core->repr ( Function::from_power ) + "^";
	std::stringstream ss;
	if ( exponent >= 0. )  ss << exponent;
	else  ss << "(" << exponent << ")";
	return s+ss.str();                      }

std::string Function::Sin::repr ( const Function::From & from ) const
{	return "sin" + this->base.core->repr ( Function::from_function );  }

std::string Function::Cos::repr ( const Function::From & from ) const
{	return "cos" + this->base.core->repr ( Function::from_function );  }

std::string Function::Step::repr ( const Function::From & from ) const
{	return "step";  }

std::string Function::JustTesting::repr ( const Function::From & from ) const
{	return this->name;  }

std::string Function::Vector::repr ( const Function::From & from ) const
{	assert ( false );  }

std::string Function::CoupledWithField::Scalar::repr ( const Function::From & from ) const
{	return "scalar";  }

#endif // DEBUG


Function maniFEM::operator+ ( const Function & f, const Function & g )

{	// both should be scalar :
	std::shared_ptr < Function::Scalar > f_scalar =
		std::dynamic_pointer_cast < Function::Scalar > ( f.core );
	std::shared_ptr < Function::Scalar > g_scalar =
		std::dynamic_pointer_cast < Function::Scalar > ( g.core );
	assert ( f_scalar );  assert ( g_scalar );

	// if one of them is zero :
	std::shared_ptr < Function::Constant > f_const =
		std::dynamic_pointer_cast < Function::Constant > ( f.core );
	if ( f_const )
		if ( f_const->val == 0. ) return g;
	std::shared_ptr < Function::Constant > g_const =
		std::dynamic_pointer_cast < Function::Constant > ( g.core );
	if ( g_const )
		if ( g_const->val == 0. ) return f;
	
	// if one of them is a sum, or both :
	std::shared_ptr < Function::Sum > f_sum =
		std::dynamic_pointer_cast < Function::Sum > ( f.core );
	std::shared_ptr < Function::Sum > g_sum =
		std::dynamic_pointer_cast < Function::Sum > ( g.core );

  Function::Sum * result = new Function::Sum;  // empty sum
	if ( g_sum )  // g is a sum
	{	std::forward_list<Function>::iterator it_g;
		for ( it_g = g_sum->terms.begin(); it_g != g_sum->terms.end(); it_g++ )
			result->terms.push_front ( *it_g );                                   }
	else  result->terms.push_front ( g );
	if ( f_sum )  // f is a sum
	{	std::forward_list<Function>::iterator it_f;
		for ( it_f = f_sum->terms.begin(); it_f != f_sum->terms.end(); it_f++ )
			result->terms.push_front ( *it_f );                                   }
	else  result->terms.push_front ( f );

	return Function ( tag::whose_core_is, result );                              }


Function maniFEM::operator* ( const Function & f, const Function & g )

{	// both should be scalar
	std::shared_ptr < Function::Scalar > f_scalar =
		std::dynamic_pointer_cast < Function::Scalar > ( f.core );
	std::shared_ptr < Function::Scalar > g_scalar =
		std::dynamic_pointer_cast < Function::Scalar > ( g.core );
	assert ( f_scalar );  assert ( g_scalar );
	
	// if any one of them is zero or one :
	std::shared_ptr < Function::Constant > f_const =
		std::dynamic_pointer_cast < Function::Constant > ( f.core );
	if ( f_const )
	{	if ( f_const->val == 0. ) return f;
		if ( f_const->val == 1. ) return g;  }
	std::shared_ptr < Function::Constant > g_const =
		std::dynamic_pointer_cast < Function::Constant > ( g.core );
	if ( g_const )
	{	if ( g_const->val == 0. ) return g;
		if ( g_const->val == 1. ) return f;  }
	
	// if one of them is a product, or both :
	std::shared_ptr < Function::Product > f_prod =
		std::dynamic_pointer_cast < Function::Product > ( f.core );
	std::shared_ptr < Function::Product > g_prod =
		std::dynamic_pointer_cast < Function::Product > ( g.core );

	Function::Product * result = new Function::Product;  // empty product
	if ( g_prod )  // g is a product
	{	std::forward_list<Function>::iterator it_g;
		for ( it_g = g_prod->factors.begin(); it_g != g_prod->factors.end(); it_g++ )
			result->factors.push_front ( *it_g );                             }
	else  result->factors.push_front ( g );
	if ( f_prod )  // f is a product
	{	std::forward_list<Function>::iterator it_f;
		for ( it_f = f_prod->factors.begin(); it_f != f_prod->factors.end(); it_f++ )
			result->factors.push_front ( *it_f );                             }
	else  result->factors.push_front ( f );

	return Function ( tag::whose_core_is, result );                              }


Function maniFEM::power ( const Function & f, double e )

{	// if f is a product, return a product, otherwise return a power
	// more simplifications can be done, but for the moment we don't go much deeper

	std::shared_ptr < Function::Scalar > f_scalar =
		std::dynamic_pointer_cast < Function::Scalar > ( f.core );
	assert ( f_scalar );
	
	std::shared_ptr < Function::Constant > f_const =
		std::dynamic_pointer_cast < Function::Constant > ( f.core );
	if ( f_const )
	{	if ( f_const->val == 0. ) return Function ( 0. );
		if ( f_const->val == 1. ) return Function ( 1. );
		return Function ( pow ( f_const->val, e ) );   }

	std::shared_ptr < Function::Product > f_prod =
		std::dynamic_pointer_cast < Function::Product > ( f.core );

	if ( f_prod )  // f is a product
	{ Function::Product * result = new Function::Product;  // empty product
		std::forward_list<Function>::iterator it;
		for ( it = f_prod->factors.begin(); it != f_prod->factors.end(); it++ )
		{	Function g = *it;
			result->factors.push_front ( power ( g, e ) );  }
		return Function ( tag::whose_core_is, result );                              }

	return Function ( tag::whose_core_is, new Function::Power ( f, e ) );              }


Function Function::Constant::deriv ( Function ) const
//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
{	return Function ( 0. );  }


Function Function::Step::deriv ( Function x ) const
//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
{	std::vector < Function > derivs;
	for ( size_t i = 0; i < this->values.size(); i++ )
		derivs.push_back ( this->values[i].deriv ( x ) );
	return Function ( tag::whose_core_is, new Function::Step ( this->arg, derivs, this->cuts ) );  }


Function Function::Sum::deriv ( Function x ) const
//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
{	std::forward_list<Function>::const_iterator it = this->terms.begin();
	Function result = 0.;
	for ( ; it != this->terms.end(); it++ )
		result += it->core->deriv(x);
	return result;                                                    }

Function Function::Product::deriv ( Function x ) const
//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
{	Function result = 0.;
	std::forward_list<Function>::const_iterator it1, it2;
	size_t c1, c2;
	for ( it1 = this->factors.begin(), c1 = 0;
        it1 != this->factors.end(); it1++, c1++ )
	{ Function partial_res = 1.;
		for ( it2 = this->factors.begin(), c2 = 0;
	        it2 != this->factors.end(); it2++, c2++ )
			if ( c1 == c2 )  // later we can eliminate c1 and c2
			{	assert ( it1 == it2 );
				partial_res *= it2->core->deriv(x);  }
			else				
			{	assert ( it1 != it2 );
				partial_res *= *it2;  }
		result += partial_res;                              }
	return result;                                            }

Function Function::Power::deriv ( Function x ) const
//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
{	Function result = this->exponent;
	result *= power ( this->base, this->exponent - 1. );
	result *= this->base.core->deriv ( x );
	return result;                                                           }

Function Function::Sin::deriv ( Function x ) const
//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
{	return cos ( this->base ) * this->base.deriv(x);  }

Function Function::Cos::deriv ( Function x ) const
//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
{	return - sin ( this->base ) * this->base.deriv(x);  }

Function Function::JustTesting::deriv ( Function x ) const
{ if ( this == x.core.get() ) return Function ( 1. );
	return Function ( 0. );                             }

Function Function::CoupledWithField::Scalar::deriv ( Function x ) const
{ if ( this == x.core.get() ) return Function ( 1. );
	return Function ( 0. );                        }

Function Function::Vector::deriv ( Function x ) const
{ assert ( false );  }

