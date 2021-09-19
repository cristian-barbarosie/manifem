
// function.cpp 2021.09.17

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

#include "math.h"
#include <sstream>

#include "function.h"
#include "manifold.h"

using namespace maniFEM;

//-----------------------------------------------------------------------------------------//

unsigned int Function::total_cores { 0 };
size_t Function::Action::counter { 0 };

#ifndef NDEBUG
std::map < const Function::Core*, std::string > Function::name;
#endif

//-----------------------------------------------------------------------------------------//

bool Function::less_for_map ( const Function & f, const Function & g )  //static
{	return f.core < g.core;  }
// needed for map 'jacobian' in class Function::Map
// and for map 'equations' in class Manifold::Parametric

//-----------------------------------------------------------------------------------------//

size_t Function::Scalar::nb_of_components ( ) const  // virtual from Function::Core
{	return 1;  }

size_t Function::Aggregate::nb_of_components ( ) const
// virtual from Function::Core, through Function::Vector
{	return this->components.size();  }

size_t Function::Immersion::nb_of_components ( ) const
// virtual from Function::Core, through Function::Vector
{ assert ( false );  }

size_t Function::CoupledWithField::Vector::nb_of_components ( ) const
// virtual from Function::Core, defined in Function::Aggregate, here overridden
{	assert ( this->components.size() == this->field->nb_of_components() );
	return this->field->nb_of_components();                          }

size_t Function::Vector::MultiValued::nb_of_components ( ) const
// virtual from Function::Core
{	return this->base.nb_of_components();  }

//-----------------------------------------------------------------------------------------//


Function Function::Scalar::component ( size_t i )
// virtual from Function::Core
// never actually used because Function::operator[] returns self
{	assert ( i == 0 );
	return Function ( tag::whose_core_is, this );  }

Function Function::Aggregate::component ( size_t i )
// virtual from Function::Core
{	assert ( i < this->components.size() );
	return this->components[i];             }

Function Function::Immersion::component ( size_t i )
// virtual from Function::Core
{	assert ( false );  }

Function Function::CoupledWithField::Vector::component ( size_t i )
// virtual from Function::Core, defined in Function::Aggregate, here overridden
{	size_t n = this->field->nb_of_components();
	assert ( i < n );
	assert ( this->components.size() == n );
	return this->components[i];                   }

Function Function::Vector::MultiValued::JumpIsSum::component ( size_t i )
// virtual from Function::Core
// build a new Function::Scalar::MultiValued::JumpIsSum
{	assert ( false );  }

Function Function::Vector::MultiValued::JumpIsLinear::component ( size_t i )
// difficult case ! it's impossible to compute separately the value of
// a component of such a multifunction, the jump depends on all other components
{	assert ( false );  }


//-----------------------------------------------------------------------------------------//
	

void Function::Constant::set_value ( double v )
// virtual from Function::Scalar
{ this->value = v;  }

double Function::Constant::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Scalar
{ return this->value;  }

double Function::Constant::get_value_on_cell
( Cell::Core *, const tag::Spin &, const Function::CompositionOfActions & exp ) const
// virtual from Function::Scalar
{ return this->value;  }

void Function::Sum::set_value ( double v )
// virtual from Function::Scalar
{ assert ( false );  }

double Function::Sum::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Scalar
{ std::forward_list<Function>::const_iterator it = this->terms.begin();
	double sum = 0.;
	for ( ; it != this->terms.end(); it++ )
	{	Function::Scalar * term_scalar = tag::Util::assert_cast
		< Function::Core*, Function::Scalar* > ( it->core );
		sum += term_scalar->get_value_on_cell ( cll );       }
	return sum;                                                           }

double Function::Sum::get_value_on_cell
	( Cell::Core * cll, const tag::Spin &, const Function::CompositionOfActions & exp ) const
// virtual from Function::Scalar
{ std::forward_list<Function>::const_iterator it = this->terms.begin();
	double sum = 0.;
	for ( ; it != this->terms.end(); it++ )
	{	Function::Scalar * term_scalar = tag::Util::assert_cast
		< Function::Core*, Function::Scalar* > ( it->core );
		sum += term_scalar->get_value_on_cell ( cll, tag::spin, exp );  }
	return sum;                                                           }

void Function::Product::set_value ( double v )
// virtual from Function::Scalar
{ assert ( false );  }

double Function::Product::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Scalar
{ std::forward_list<Function>::const_iterator it = this->factors.begin();
	double prod = 1.;
	for ( ; it != this->factors.end(); it++ )
	{	Function::Scalar * fact_scalar = tag::Util::assert_cast
			< Function::Core*, Function::Scalar* > ( it->core );
		prod *= fact_scalar->get_value_on_cell ( cll );        }
	return prod;                                                             }
	
double Function::Product::get_value_on_cell
( Cell::Core * cll, const tag::Spin &, const Function::CompositionOfActions & exp ) const
// virtual from Function::Scalar
{ std::forward_list<Function>::const_iterator it = this->factors.begin();
	double prod = 1.;
	for ( ; it != this->factors.end(); it++ )
	{	Function::Scalar * fact_scalar = tag::Util::assert_cast
			< Function::Core*, Function::Scalar* > ( it->core );
		prod *= fact_scalar->get_value_on_cell ( cll, tag::spin, exp );  }
	return prod;                                                            }
	
void Function::Power::set_value ( double v )
// virtual from Function::Scalar
{ assert ( false );  }

double Function::Power::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Scalar
{	Function::Scalar * base_scalar = tag::Util::assert_cast
		< Function::Core*, Function::Scalar* > ( this->base.core );
	return std::pow ( base_scalar->get_value_on_cell(cll), this->exponent );  }
	
double Function::Power::get_value_on_cell
( Cell::Core * cll, const tag::Spin &, const Function::CompositionOfActions & exp ) const
// virtual from Function::Scalar
{	Function::Scalar * base_scalar = tag::Util::assert_cast
		< Function::Core*, Function::Scalar* > ( this->base.core );
	return std::pow ( base_scalar->get_value_on_cell ( cll, tag::spin, exp ), this->exponent );  }
	
void Function::Sqrt::set_value ( double v )
// virtual from Function::Scalar
{ assert ( false );  }

double Function::Sqrt::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Scalar
{	Function::Scalar * base_scalar = tag::Util::assert_cast
		< Function::Core*, Function::Scalar* > ( this->base.core );
	return std::sqrt ( base_scalar->get_value_on_cell(cll) );      }
	
double Function::Sqrt::get_value_on_cell
( Cell::Core * cll, const tag::Spin &, const Function::CompositionOfActions & exp ) const
// virtual from Function::Scalar
{	Function::Scalar * base_scalar = tag::Util::assert_cast
		< Function::Core*, Function::Scalar* > ( this->base.core );
	return std::sqrt ( base_scalar->get_value_on_cell ( cll, tag::spin, exp ) );  }
	
void Function::Sin::set_value ( double v )
// virtual from Function::Scalar
{ assert ( false );  }

double Function::Sin::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Scalar
{	Function::Scalar * base_scalar = tag::Util::assert_cast
		< Function::Core*, Function::Scalar* > ( this->base.core );
	return std::sin ( base_scalar->get_value_on_cell(cll) );      }
	
double Function::Sin::get_value_on_cell
( Cell::Core * cll, const tag::Spin &, const Function::CompositionOfActions & exp ) const
// virtual from Function::Scalar
{	Function::Scalar * base_scalar = tag::Util::assert_cast
		< Function::Core*, Function::Scalar* > ( this->base.core );
	return std::sin ( base_scalar->get_value_on_cell ( cll, tag::spin, exp ) );  }
	
void Function::Cos::set_value ( double v )
// virtual from Function::Scalar
{ assert ( false );  }

double Function::Cos::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Scalar
{	Function::Scalar * base_scalar = tag::Util::assert_cast
		< Function::Core*, Function::Scalar* > ( this->base.core );
	return std::cos ( base_scalar->get_value_on_cell(cll) );      }
	
double Function::Cos::get_value_on_cell
( Cell::Core * cll, const tag::Spin &, const Function::CompositionOfActions & exp ) const
// virtual from Function::Scalar
{	Function::Scalar * base_scalar = tag::Util::assert_cast
		< Function::Core*, Function::Scalar* > ( this->base.core );
	return std::cos ( base_scalar->get_value_on_cell ( cll, tag::spin, exp ) );  }
	
void Function::Step::set_value ( double v )
// virtual from Function::Scalar
{ assert ( false );  }

double Function::Step::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Scalar
{	Function::Scalar * arg_scalar = tag::Util::assert_cast
			< Function::Core*, Function::Scalar* > ( this->arg.core );
	double arg_v = arg_scalar->get_value_on_cell(cll);
	for ( size_t i = 0; i < this->cuts.size(); i++ )
		if ( arg_v < cuts[i] )
		{	Function::Scalar * val_i_scalar = tag::Util::assert_cast
				< Function::Core*, Function::Scalar* > ( this->values[i].core );
			return val_i_scalar->get_value_on_cell(cll);                        }
	Function::Scalar * val_scalar = tag::Util::assert_cast
		< Function::Core*, Function::Scalar* > ( this->values.back().core );
	return val_scalar->get_value_on_cell(cll);                                }
	
double Function::Step::get_value_on_cell
( Cell::Core * cll, const tag::Spin &, const Function::CompositionOfActions & exp ) const
// virtual from Function::Scalar
{	Function::Scalar * arg_scalar = tag::Util::assert_cast
			< Function::Core*, Function::Scalar* > ( this->arg.core );
	double arg_v = arg_scalar->get_value_on_cell ( cll, tag::spin, exp );
	for ( size_t i = 0; i < this->cuts.size(); i++ )
		if ( arg_v < cuts[i] )
		{	Function::Scalar * val_i_scalar = tag::Util::assert_cast
				< Function::Core*, Function::Scalar* > ( this->values[i].core );
			return val_i_scalar->get_value_on_cell ( cll, tag::spin, exp );    }
	Function::Scalar * val_scalar = tag::Util::assert_cast
		< Function::Core*, Function::Scalar* > ( this->values.back().core );
	return val_scalar->get_value_on_cell ( cll, tag::spin, exp );             }
	
std::vector<double> Function::Aggregate::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Vector
{ size_t n = this->nb_of_components();
	std::vector<double> result ( n );
	for ( size_t i = 0; i < n; i++ )
	{	Function::Scalar * comp_scalar = tag::Util::assert_cast
			< Function::Core*, Function::Scalar* > ( this->components[i].core );
		result[i] = comp_scalar->get_value_on_cell ( cll );                      }
	return result;                                                                }

std::vector<double> Function::Aggregate::get_value_on_cell
( Cell::Core * cll, const tag::Spin &, const Function::CompositionOfActions & exp ) const
// virtual from Function::Vector
{ assert ( false );  }

void Function::CoupledWithField::Scalar::set_value ( double v )
// virtual from Function::Scalar
{ assert ( false );  }

double Function::CoupledWithField::Scalar::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Scalar
{ return this->field->on_cell(cll);  }
// { return this->field->on_cell(cll).reference();  }
	
double Function::CoupledWithField::Scalar::get_value_on_cell
( Cell::Core * cll, const tag::Spin &, const Function::CompositionOfActions & exp ) const
// virtual from Function::Scalar
{ return this->field->on_cell(cll);  }
// { return this->field->on_cell(cll).reference();  }

void Function::Diffeomorphism::OneDim::set_value ( double v )
// virtual from Function::Scalar
{ assert ( false );  }
	
double Function::Diffeomorphism::OneDim::get_value_on_cell ( Cell::Core * cll ) const
{ assert ( false );  }  // virtual from Function::Scalar
	
double Function::Diffeomorphism::OneDim::get_value_on_cell
( Cell::Core * cll, const tag::Spin &, const Function::CompositionOfActions & exp ) const
{ assert ( false );  }  // virtual from Function::Scalar
	
std::vector<double> Function::CoupledWithField::Vector::get_value_on_cell
( Cell::Core * cll ) const  // virtual from Function::Vector
{ return this->field->on_cell(cll);  }
// { return this->field->on_cell(cll).reference();  }
	
std::vector<double> Function::CoupledWithField::Vector::get_value_on_cell
( Cell::Core * cll, const tag::Spin &, const Function::CompositionOfActions & exp ) const
{ return this->field->on_cell(cll);  }  // virtual from Function::Vector
// { return this->field->on_cell(cll).reference();  }
	
std::vector<double> Function::Immersion::get_value_on_cell ( Cell::Core * cll ) const
{ assert ( false );  }  // virtual from Function::Vector

std::vector<double> Function::Immersion::get_value_on_cell
( Cell::Core * cll, const tag::Spin &, const Function::CompositionOfActions & exp ) const
{ assert ( false );  }  // virtual from Function::Vector

void Function::Composition::set_value ( double v )  // virtual from Function::Scalar
{ assert ( false );  }

double Function::Composition::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Scalar
{	Function::Scalar * base_scalar = tag::Util::assert_cast
		< Function::Core*, Function::Scalar* > ( this->base.core );
	// we assume here that cll lives in master coordinates
	return base_scalar->get_value_on_cell(cll);                   }

double Function::Composition::get_value_on_cell  // virtual from Function::Scalar
( Cell::Core * cll, const tag::Spin &, const Function::CompositionOfActions & exp ) const
{ assert ( false );  }
	
//-----------------------------------------------------------------------------------------//


void Function::Scalar::MultiValued::set_value ( double v )  // virtual from Function::Scalar
{	Function::Scalar * base_scalar = tag::Util::assert_cast
		< Function::Core*, Function::Scalar* > ( this->base.core );
	base_scalar->set_value ( v );                                 }

//-----------------------------------------------------------------------------------------//


double Function::Scalar::MultiValued::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Scalar
{	Function::Scalar * base_scalar = tag::Util::assert_cast
		< Function::Core*, Function::Scalar* > ( this->base.core );
	return base_scalar->get_value_on_cell(cll);                   }

//-----------------------------------------------------------------------------------------//


double Function::Scalar::MultiValued::JumpIsSum::get_value_on_cell  // virtual from Function::Scalar
( Cell::Core * cll, const tag::Spin &, const Function::CompositionOfActions & exp ) const

{	Function::Scalar * base_scalar = tag::Util::assert_cast
		< Function::Core*, Function::Scalar* > ( this->base.core );
	Manifold::Quotient * manif = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( Manifold::working.core );
	size_t n = this->actions.size();
	assert ( manif->actions == this->actions );
	assert ( this->beta.size() == n );
  double val = base_scalar->get_value_on_cell ( cll );
	for ( size_t i = 0; i < n; i++ )
	{	std::map<Function::Action,short int>::const_iterator it =
			exp.index_map.find ( this->actions[i] );
		assert ( it != exp.index_map.end() );
		val += it->second * this->beta[i];                         }
	return val;                                                          }

//-----------------------------------------------------------------------------------------//


double Function::Scalar::MultiValued::JumpIsLinear::get_value_on_cell
// virtual from Function::Scalar
( Cell::Core * cll, const tag::Spin &, const Function::CompositionOfActions & exp ) const

{	Function::Scalar * base_scalar = tag::Util::assert_cast
		< Function::Core*, Function::Scalar* > ( this->base.core );
	Manifold::Quotient * manif = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( Manifold::working.core );
	assert ( manif->actions == this->actions );
	size_t n = this->alpha.size();
	assert ( this->actions.size() == n );
  double val = base_scalar->get_value_on_cell ( cll ) - this->gamma;
	for ( size_t i = 0; i < n; i++ )
	{	std::map<Function::Action,short int>::const_iterator it =
			exp.index_map.find ( this->actions[i] );
		if ( it == exp.index_map.end() ) continue;
		short int exp_i = it->second;
		assert ( exp_i != 0 );
		double alpha_i = this->alpha[i];
		if ( exp_i > 0 )
		{	size_t abs_exp_i = exp_i;
			for ( size_t j = 0; j < abs_exp_i; j++ ) val *= alpha_i;  }
		if ( exp_i < 0 )
		{	size_t abs_exp_i = -exp_i;
			for ( size_t j = 0; j < abs_exp_i; j++ ) val /= alpha_i;  }  }
	return val + this->gamma;                                             }

//-----------------------------------------------------------------------------------------//


std::vector<double> Function::Vector::MultiValued::get_value_on_cell
( Cell::Core * cll ) const  // virtual from Function::Vector
{	Function::Vector * base_vector = tag::Util::assert_cast
		< Function::Core*, Function::Vector* > ( this->base.core );
	return base_vector->get_value_on_cell(cll);                   }

//-----------------------------------------------------------------------------------------//


std::vector<double> Function::Vector::MultiValued::JumpIsSum::get_value_on_cell
// virtual from Function::Vector
( Cell::Core * cll, const tag::Spin &, const Function::CompositionOfActions & exp ) const
	
{	Function::Vector * base_vector = tag::Util::assert_cast
		< Function::Core*, Function::Vector* > ( this->base.core );
	Manifold::Quotient * manif = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( Manifold::working.core );
	size_t n = this->actions.size();
	size_t dim = base_vector->nb_of_components();
	assert ( manif->actions == this->actions );
	assert ( this->beta.size() == n );
	std::vector < double > val = base_vector->get_value_on_cell ( cll );
	for ( size_t i = 0; i < n; i++ )
	{	assert ( beta[i].size() == dim );
		std::map<Function::Action,short int>::const_iterator it =
			exp.index_map.find ( this->actions[i] );
		if ( it == exp.index_map.end() ) continue;
		short int exp_i = it->second;
		assert ( exp_i != 0 );
		for ( size_t j = 0; j < dim; j++ )
			val[j] += exp_i * this->beta[i][j];                       }
	return val;                                                            }

//-----------------------------------------------------------------------------------------//


std::vector<double> Function::Vector::MultiValued::JumpIsLinear::get_value_on_cell
// virtual from Function::Vector
( Cell::Core * cll, const tag::Spin &, const Function::CompositionOfActions & exp ) const
	
{	Function::Vector * base_vector = tag::Util::assert_cast
		< Function::Core*, Function::Vector* > ( this->base.core );
	Manifold::Quotient * manif = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( Manifold::working.core );
	size_t n = this->actions.size();
	size_t dim = base_vector->nb_of_components();
	assert ( manif->actions == this->actions );
	assert ( this->A.size() == n );
	assert ( this->b.size() == n );
	std::vector < double > val = base_vector->get_value_on_cell ( cll );
	for ( size_t i = 0; i < n; i++ )
	{	assert ( this->A[i].size() == dim );
		std::map<Function::Action,short int>::const_iterator it =
			exp.index_map.find ( this->actions[i] );
		if ( it == exp.index_map.end() ) continue;
		short int exp_i = it->second;
		if ( exp_i > 0 )
		{	size_t abs_exp_i = exp_i;
			for ( size_t j = 0; j < abs_exp_i; j++ )
			{	std::vector < double > temp ( dim, 0. );
				for ( size_t k = 0; k < dim; k++ )
				for ( size_t l = 0; l < dim; l++ )
					temp[k] += this->A[i][k][l]*val[l];
				for ( size_t k = 0; k < dim; k++ ) val[k] = temp[k] + this->b[i][k];  }  }
		if ( exp_i < 0 )
		{	size_t abs_exp_i = -exp_i;
			for ( size_t j = 0; j < abs_exp_i; j++ )
				{	std::vector < double > temp ( dim, 0. );
				for ( size_t k = 0; k < dim; k++ )
				for ( size_t l = 0; l < dim; l++ )
					temp[k] += this->Ainv[i][k][l]*(val[l]-this->b[i][l]);
				val = std::move ( temp );                                }  }              }
	return val;                                                                         }

//-----------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------//
	

double Function::ArithmeticExpression::set_value_on_cell
( Cell::Core * cll, const double & )  // virtual from Function::Scalar
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
            << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot assign to an arithmetic expression." << std::endl;
	exit ( 1 );                                                              }
	
std::vector<double> Function::Aggregate::set_value_on_cell
( Cell::Core * cll, const std::vector<double> & x )  // virtual from Function::Vector
// there should be a faster way !!
{ size_t n = this->nb_of_components();
	assert ( n == x.size() );
	for ( size_t i = 0; i < n; i++ )
	{	Function::Scalar * comp_scalar = tag::Util::assert_cast
			< Function::Core*, Function::Scalar* > ( this->components[i].core );
		comp_scalar->set_value_on_cell ( cll, x[i] );               	         }
	return x;     		                                                         }

double Function::CoupledWithField::Scalar::set_value_on_cell
( Cell::Core * cll, const double & x )  // virtual from Function::Scalar
{ return this->field->on_cell(cll) = x;  }

double Function::Diffeomorphism::OneDim::set_value_on_cell
( Cell::Core * cll, const double & x )  // virtual from Function::Scalar
{ assert ( false );  }

std::vector<double> Function::CoupledWithField::Vector::set_value_on_cell
( Cell::Core * cll, const std::vector<double> & x )  // virtual from Function::Vector
{ return this->field->on_cell(cll) = x;  }

std::vector<double> Function::Immersion::set_value_on_cell
( Cell::Core * cll, const std::vector<double> & x )  // virtual from Function::Vector
{ assert ( false );  }

double Function::Composition::set_value_on_cell
( Cell::Core * cll, const double & x )  // virtual from Function::Vector
{ assert ( false );  }

double Function::Scalar::MultiValued::set_value_on_cell
( Cell::Core * cll, const double & x )  // virtual from Function::Scalar
{ assert ( false );  }

std::vector < double > Function::Vector::MultiValued::set_value_on_cell
( Cell::Core * cll, const std::vector<double> & x )  // virtual from Function::Scalar
{ assert ( false );  }

//-----------------------------------------------------------------------------------------//

#ifndef NDEBUG	

std::string Function::Constant::repr ( const Function::From & from ) const
{	std::stringstream ss;
	ss << this->value; 
	std::string s = ss.str();
	if ( ( this->value < 0. ) and ( from != Function::from_void ) ) s = "(" + s + ")";
	return s;                         }

std::string Function::Sum::repr ( const Function::From & from ) const
{ std::forward_list<Function>::const_iterator it = this->terms.begin();
	assert ( it != terms.end() );
	std::string s = it->core->repr ( Function::from_sum );
	for ( it++; it != terms.end(); it++ )
		s = s + '+' + it->core->repr ( Function::from_sum );
	if ( ( from == Function::from_power ) or ( from == Function::from_product )
       or ( from == Function::from_function ) )  s = "(" + s + ")";
	return s;                                                                      }

std::string Function::Product::repr ( const Function::From & from ) const
{ std::forward_list<Function>::const_iterator it = this->factors.begin();
	assert ( it != factors.end() );
	std::string s = it->core->repr ( Function::from_product );
	for ( it++; it != factors.end(); it++ )
		s = s + '*' + it->core->repr ( Function::from_product );
	if ( ( from == Function::from_power ) or ( from == Function::from_function ) )
		s = "(" + s + ")";
	return s;                                                    }

std::string Function::Power::repr ( const Function::From & from ) const
{	std::string s = this->base.core->repr ( Function::from_power ) + "^";
	std::stringstream ss;
	if ( exponent >= 0. )  ss << exponent;
	else  ss << "(" << exponent << ")";
	return s+ss.str();                      }

std::string Function::Sqrt::repr ( const Function::From & from ) const
{	return "sqrt" + this->base.core->repr ( Function::from_function );  }

std::string Function::Sin::repr ( const Function::From & from ) const
{	return "sin" + this->base.core->repr ( Function::from_function );  }

std::string Function::Cos::repr ( const Function::From & from ) const
{	return "cos" + this->base.core->repr ( Function::from_function );  }

std::string Function::Step::repr ( const Function::From & from ) const
{	return "step";  }

std::string Function::Vector::repr ( const Function::From & from ) const
{	return "vector";  }

std::string Function::CoupledWithField::Scalar::repr ( const Function::From & from ) const
{	if ( Function::name.find(this) != Function::name.end() )
		return Function::name[this];
	std::stringstream ss;
	ss << this->field;
	return "scalar"+ss.str();                                      }

std::string Function::Composition::repr ( const Function::From & from ) const
{	return this->base.core->repr ( Function::from_product ) + "º";  }

std::string Function::Diffeomorphism::OneDim::repr ( const Function::From & from ) const
{ assert ( false );  }

std::string Function::Scalar::MultiValued::repr ( const Function::From & from ) const
{	return "multi" + this->base.core->repr ( Function::from_function );  }

std::string Function::Vector::MultiValued::repr ( const Function::From & from ) const
{ assert ( false );  }

#endif // DEBUG

//-----------------------------------------------------------------------------------------//


Function maniFEM::operator+ ( const Function & f, const Function & g )

{	// both should be scalar
	assert ( dynamic_cast < Function::Scalar* > ( f.core ) );
	assert ( dynamic_cast < Function::Scalar* > ( g.core ) );

	// if one of them is zero :
  Function::Constant * f_const = dynamic_cast < Function::Constant * > ( f.core );
	if ( f_const )  if ( f_const->value == 0. )  return g;
  Function::Constant * g_const = dynamic_cast < Function::Constant * > ( g.core );
	if ( g_const )  if ( g_const->value == 0. )  return f;
	
	// if one of them is a sum, or both :
  Function::Sum * f_sum = dynamic_cast < Function::Sum * > ( f.core );
  Function::Sum * g_sum = dynamic_cast < Function::Sum * > ( g.core );

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

//-----------------------------------------------------------------------------------------//


Function maniFEM::operator* ( const Function & f, const Function & g )

{	// both should be scalar
	assert ( dynamic_cast < Function::Scalar* > ( f.core ) );
	assert ( dynamic_cast < Function::Scalar* > ( g.core ) );
	
	// if any one of them is zero or one :
  Function::Constant * f_const = dynamic_cast < Function::Constant * > ( f.core );
	if ( f_const )
	{	if ( f_const->value == 0. ) return f;
		if ( f_const->value == 1. ) return g;
		// is g a product ?
		Function::Product * g_prod = dynamic_cast < Function::Product * > ( g.core );
		if ( g_prod )
		{	Function::Product * result = new Function::Product;  // empty product
			bool constant_not_yet_used = true;
			// we make a copy of the list of factors so to keep the original order
			std::forward_list < Function > list_fact;
			std::forward_list<Function>::iterator it_g;
			for ( it_g = g_prod->factors.begin(); it_g != g_prod->factors.end(); it_g++ )
				list_fact.push_front ( *it_g );
			for ( it_g = list_fact.begin(); it_g != list_fact.end(); it_g++ )
			{	Function fact = *it_g;
				if ( constant_not_yet_used )
				{	Function::Constant * fact_const =
						dynamic_cast < Function::Constant * > ( fact.core );
					if ( fact_const )
					{	constant_not_yet_used = false;
						result->factors.push_front ( f * fact );  }
					else result->factors.push_front ( fact );              }
				else result->factors.push_front ( fact );                   }
			return Function ( tag::whose_core_is, result );                               }  }
  Function::Constant * g_const = dynamic_cast < Function::Constant * > ( g.core );
	if ( g_const )
	{	if ( g_const->value == 0. ) return g;
		if ( g_const->value == 1. ) return f;
		// is f a product ?
		Function::Product * f_prod = dynamic_cast < Function::Product * > ( f.core );
		if ( f_prod )
		{	Function::Product * result = new Function::Product;  // empty product
			bool constant_not_yet_used = true;
			// we make a copy of the list of factors so to keep the original order
			std::forward_list < Function > list_fact;
			std::forward_list<Function>::iterator it_f;
			for ( it_f = f_prod->factors.begin(); it_f != f_prod->factors.end(); it_f++ )
				list_fact.push_front ( *it_f );
			for ( it_f = list_fact.begin(); it_f != list_fact.end(); it_f++ )
			{	Function fact = *it_f;
				if ( constant_not_yet_used )
				{	Function::Constant * fact_const =
						dynamic_cast < Function::Constant * > ( fact.core );
					if ( fact_const )
					{	constant_not_yet_used = false;
						result->factors.push_front ( g * fact );  }
					else result->factors.push_front ( fact );              }
				else result->factors.push_front ( fact );                   }
			return Function ( tag::whose_core_is, result );                               }  }

	// if both are constant :
	if ( f_const and g_const ) return Function ( f_const->value * g_const->value );
	
	// if one of them is a product, or both :
  Function::Product * f_prod = dynamic_cast < Function::Product * > ( f.core );
  Function::Product * g_prod = dynamic_cast < Function::Product * > ( g.core );

	Function::Product * result = new Function::Product;  // empty product
	if ( g_prod )  // g is a product
	{	std::forward_list<Function>::iterator it_g;
		for ( it_g = g_prod->factors.begin(); it_g != g_prod->factors.end(); it_g++ )
			result->factors.push_front ( *it_g );                                       }
	else  result->factors.push_front ( g );
	if ( f_prod )  // f is a product
	{	std::forward_list<Function>::iterator it_f;
		for ( it_f = f_prod->factors.begin(); it_f != f_prod->factors.end(); it_f++ )
			result->factors.push_front ( *it_f );                                       }
	else  result->factors.push_front ( f );

	return Function ( tag::whose_core_is, result );                                   }

//-----------------------------------------------------------------------------------------//


Function maniFEM::power ( const Function & f, double e )

{	assert ( dynamic_cast < Function::Scalar* > ( f.core ) );

	if ( e == 0. ) return Function ( 1. );
	if ( e == 1. ) return f;
	if ( e == 2. ) return f*f;
	
  Function::Constant * f_const = dynamic_cast < Function::Constant * > ( f.core );
	if ( f_const )
	{	if ( f_const->value == 0. ) return Function ( 0. );
		if ( f_const->value == 1. ) return Function ( 1. );
		return Function ( pow ( f_const->value, e ) );       }

  Function::Power * f_pow = dynamic_cast < Function::Power * > ( f.core );
	if ( f_pow ) return power ( f_pow->base, f_pow->exponent * e );

  Function::Product * f_prod = dynamic_cast < Function::Product * > ( f.core );
	if ( f_prod )  // f is a product
	{ Function::Product * result = new Function::Product;  // empty product
		std::forward_list<Function>::iterator it;
		for ( it = f_prod->factors.begin(); it != f_prod->factors.end(); it++ )
		{	Function g = *it;
			result->factors.push_front ( power ( g, e ) );  }
		return Function ( tag::whose_core_is, result );                              }

	return Function ( tag::whose_core_is, new Function::Power ( f, e ) );              }

//-----------------------------------------------------------------------------------------//


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
		{	if ( c1 == c2 )  // later we can eliminate c1 and c2
			{	assert ( it1 == it2 );
				partial_res *= it2->core->deriv(x); }
			else				
			{	assert ( it1 != it2 );
				partial_res *= *it2;   } }
		result += partial_res;                              }
	return result;                                            }


Function Function::Power::deriv ( Function x ) const
//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
{	Function result = this->exponent;
	result *= power ( this->base, this->exponent - 1. );
	result *= this->base.core->deriv ( x );
	return result;                                                           }

Function Function::Sqrt::deriv ( Function x ) const
//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
{	return 0.5 * this->base.deriv(x) / sqrt ( this->base );    }

Function Function::Sin::deriv ( Function x ) const
//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
{	return cos ( this->base ) * this->base.deriv(x);  }

Function Function::Cos::deriv ( Function x ) const
//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
{	return - sin ( this->base ) * this->base.deriv(x);  }

Function Function::CoupledWithField::Scalar::deriv ( Function x ) const
{	if ( this == x.core ) return Function ( 1. );
	return Function ( 0. );                        }

Function Function::Vector::deriv ( Function x ) const
{ assert ( false );  }

Function Function::Diffeomorphism::OneDim::deriv ( Function x ) const
{ assert ( false );  }

Function Function::Scalar::MultiValued::deriv ( Function x ) const
{ assert ( false );  }

Function Function::Vector::MultiValued::deriv ( Function x ) const
{ assert ( false );  }

//-------------------------------------------------------------------------------------


Function::Diffeomorphism::OneDim::OneDim
( const Function & gc, const Function & mc, const Function & bgc, const tag::BuildJacobian & )
	
:	Function::Scalar ( ), Function::Map ( gc, mc, bgc ), Function::Diffeomorphism ( )
// 'jacobian' initialized as empty map, 'det' as non_existent

{	assert ( dynamic_cast < Function::Scalar* > ( gc.core ) );
	assert ( dynamic_cast < Function::Scalar* > ( mc.core ) );
	assert ( dynamic_cast < Function::Scalar* > ( bgc.core ) );
	assert ( gc.nb_of_components() == 1 );
	assert ( mc.nb_of_components() == 1 );
	assert ( bgc.nb_of_components() == 1 );

	// 'det' should be good for use in integration, see e.g. Integrator::Gauss::action
	this->det = bgc.deriv(mc);

	// code below should match the conventions in Function::Composition::deriv
	this->jacobian.insert ( std::pair < Function, Function > ( gc, 1. / this->det ) );  }


Function::Immersion::Immersion
( const Function & gc, const Function & mc, const Function & bgc, const tag::BuildJacobian & )

:	Function::Immersion ( gc, mc, bgc )
// 'jacobian' initialized as empty map, 'det' as non_existent

{	assert ( dynamic_cast < Function::Vector* > ( gc.core ) );
	assert ( dynamic_cast < Function::Vector* > ( bgc.core ) );
	size_t geom_dim = gc.nb_of_components();
	size_t master_dim = mc.nb_of_components();
	assert ( geom_dim > master_dim );
	assert ( geom_dim == bgc.nb_of_components() );
	for ( size_t k = 0; k < master_dim; k++ )
	{	Function t = mc[k];
		Function jacob = bgc[0].deriv(t);
		for ( size_t i = 1; i < geom_dim; i++ ) jacob = jacob && bgc[i].deriv(t);
		this->jacobian.insert ( std::pair < Function, Function > ( t, jacob ) );   }
	// this->jacobian is a non-square matrix, so it makes no sense to think about its inverse
	// this is very different from Function::Diffeomorphism::HighDim !

	// 'det' should be good for use in integration, see e.g. Integrator::Gauss::action
	if ( master_dim == 1 )
	{	assert ( dynamic_cast < Function::Scalar* > ( mc.core ) );
		Function v = this->jacobian.find(mc)->second;
		Function norm_v_2 ( 0. );
		for ( size_t i = 0; i < geom_dim; i++ ) norm_v_2 += v[i]*v[i];
		this->det = sqrt ( norm_v_2 );                                  }
	else
	{	assert ( master_dim == 2 );  // sorry, no 3d for now ...
		assert ( dynamic_cast < Function::Vector* > ( mc.core ) );
		Function v = this->jacobian.find(mc[0])->second;
		Function w = this->jacobian.find(mc[1])->second;
		// det = square root of  |v|^2 |w|^2 - (v.w)^2
		Function norm_v_2 ( 0. ), norm_w_2 ( 0. ), vw ( 0. );
		for ( size_t i = 0; i < geom_dim; i++ )
		{	norm_v_2 += v[i]*v[i];
			norm_w_2 += w[i]*w[i];
			vw       += v[i]*w[i];  }
		this->det = sqrt ( norm_v_2 * norm_w_2 - vw * vw );    }                        }


Function::Diffeomorphism::HighDim::HighDim
( const Function & gc, const Function & mc, const Function & bgc, const tag::BuildJacobian & )
	
:	Function::Immersion ( gc, mc, bgc ), Function::Diffeomorphism ( )
// 'jacobian' initialized as empty map, 'det' as non_existent

{	assert ( dynamic_cast < Function::Vector* > ( gc.core ) );
	assert ( dynamic_cast < Function::Vector* > ( mc.core ) );
	assert ( dynamic_cast < Function::Vector* > ( bgc.core ) );
	assert ( gc.nb_of_components() == mc.nb_of_components() );
	assert ( gc.nb_of_components() == bgc.nb_of_components() );
  assert ( gc.nb_of_components() >= 2 );
	assert ( mc.nb_of_components() == 2 );
	assert ( gc.nb_of_components() == 2 );  // sorry, no 3d for now ...
	Function dx_dxi = bgc[0].deriv(mc[0]);
	Function dx_deta = bgc[0].deriv(mc[1]);
	Function dy_dxi = bgc[1].deriv(mc[0]);
	Function dy_deta = bgc[1].deriv(mc[1]);

	// 'det' should be good for use in integration, see e.g. Integrator::Gauss::action
	this->det = dx_dxi * dy_deta - dx_deta * dy_dxi;
	
	// now we compute the inverse of the matrix
	// this is quite different from Function::Immersion !
	// code below should match the conventions in Function::Composition::deriv
	Function dxi_dx =   dy_deta / this->det,  deta_dx = - dy_dxi / this->det,
	         dxi_dy = - dx_deta / this->det,  deta_dy =   dx_dxi / this->det;
	this->jacobian.insert ( std::pair < Function, Function > ( gc[0], dxi_dx && deta_dx ) );
	this->jacobian.insert ( std::pair < Function, Function > ( gc[1], dxi_dy && deta_dy ) );
	assert ( this->jacobian.find(geom_coords[0]) != this->jacobian.end() );
	assert ( this->jacobian.find(geom_coords[1]) != this->jacobian.end() );                   }


Function Function::Composition::deriv ( Function x ) const
	
// there are two kinds of derivatives for a Function::Composition
// if x is a master coordinate, then we return the usual arithmetic derivative
// if x is a geometric coordinate, we apply the chain rule

// in both cases, we compose the result just before return
// since the calling code may differentiate the result again
	
{	Function::Map * transf_c = dynamic_cast < Function::Map* > ( this->transf.core );
	assert ( transf_c );
	// in the above we must use dynamic_cast
	// below, with -DNDEBUG, assert_cast calls a static_cast which does not compile
	// classes Function::Core and Function::Map are not directly related
	// perhaps we should have used virtual inheritance ?
	// Function::Map * transf_c = tag::Util::assert_cast
	// 	< Function::Core*, Function::Map* > ( this->transf.core );

	// we have a map, which may be a Function::Diffeomorphism
	// or a Function::Immersion
	// in the latter case we cannot differentiate wrt a geometric coordinate
	// only with respect to a master coordinate
	
	if ( dynamic_cast < Function::Diffeomorphism* > ( this->transf.core ) )

	{	// 'this' may be Function::Diffeomorphism::OneDim or Function::Diffeomorphism::HighDim
		// in either case, transf_c->jacobian derivatives of master coordinates
		// with respect to geometric coordinates
		std::map<Function,Function>::const_iterator it = transf_c->jacobian.find(x);
		if ( it != transf_c->jacobian.end() )
			// x is a geometric coordinate, so we apply the chain rule
		{	const Function & jacob = it->second;
			// 'jacob' contains derivatives of the master coordinates with respect to x
			Function result = 0.;
			for ( size_t i = 0; i < transf_c->master_coords.nb_of_components(); i++ )
				result += this->base.deriv(transf_c->master_coords[i]) * jacob[i];
			return Function ( result, tag::composed_with, this->transf );             }
		else  // x is a master coordinate
			return Function ( this->base.deriv(x), tag::composed_with, this->transf );      }

	else  //  not a Diffeomorphism, but an Immersion
		
	{	assert ( dynamic_cast < Function::Immersion* > ( this->transf.core ) );
		std::map<Function,Function>::const_iterator
			it = transf_c->jacobian.find(x); // DEBUG mode only !! [?]
		assert ( it != transf_c->jacobian.end() );  // x must be a master coordinate
		return Function ( this->base.deriv(x), tag::composed_with, this->transf );   }
	
}  // end of Function::Composition::deriv	

//-----------------------------------------------------------------------------------------//


Function Function::Constant::replace ( const Function & x, const Function & y )
//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
{ return Function ( tag::whose_core_is, this );  }

Function Function::Sum::replace ( const Function & x, const Function & y )
//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
{ if ( this == x.core ) return y;
	std::forward_list<Function>::const_iterator it = this->terms.begin();
	Function result = 0.;
	for ( ; it != this->terms.end(); it++ )
	{ result += it->replace ( x, y );          }
  return result;                                                                 }

Function Function::Product::replace ( const Function & x, const Function & y )
//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
{ if ( this == x.core ) return y;
	std::forward_list<Function>::const_iterator it = this->factors.begin();
	Function result = 1.;
	for ( ; it != this->factors.end(); it++ ) result *= it->replace ( x, y );
	return result;                                                                   }

Function Function::Power::replace ( const Function & x, const Function & y )
//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
{	if ( this == x.core ) return y;
	return power ( this->base.replace ( x, y ), this->exponent );  }

Function Function::Sqrt::replace ( const Function & x, const Function & y )
//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
{	if ( this == x.core ) return y;
	return sqrt ( this->base.replace ( x, y ) );  }

Function Function::Sin::replace ( const Function & x, const Function & y )
//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
{	if ( this == x.core ) return y;
	return sin ( this->base.replace ( x, y ) );  }

Function Function::Cos::replace ( const Function & x, const Function & y )
//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
{	if ( this == x.core ) return y;
	return cos ( this->base.replace ( x, y ) );  }

Function Function::Step::replace ( const Function & x, const Function & y )
//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
{	assert ( false );  }

Function Function::Vector::replace ( const Function & x, const Function & y )
//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
{	assert ( false );  }

Function Function::CoupledWithField::Scalar::replace
( const Function & x, const Function & y )
//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
{	if ( this == x.core ) return y;
	return Function ( tag::whose_core_is, this );  }

Function Function::CoupledWithField::Vector::replace
( const Function & x, const Function & y )
//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
{	assert ( false );  }

Function Function::Diffeomorphism::OneDim::replace
( const Function & x, const Function & y )
//  virtual from Function::Vector
{	assert ( false );  }

Function Function::Composition::replace ( const Function & x, const Function & y )
//  virtual from Function::Core, through Function::Vector
{	return Function ( tag::whose_core_is, this );  }

Function Function::Scalar::MultiValued::replace ( const Function & x, const Function & y )
{	assert ( false );  }  //  virtual from Function::Scalar

Function Function::Vector::MultiValued::replace ( const Function & x, const Function & y )
{	assert ( false );  }  //  virtual from Function::Vector

//--------------------------------------------------------------------------

Function::Core::~Core ( )
{	assert ( Function::total_cores > 0 );
	Function::total_cores--; }


