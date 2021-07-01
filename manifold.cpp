
// manifold.cpp 2021.04.10

//   Copyright 2019, 2020, 2021 Cristian Barbarosie cristian.barbarosie@gmail.com
//   https://github.com/cristian-barbarosie/manifem

//   This file is part of maniFEM, a C++ library for meshes and finite elements on manifolds.

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

#include "manifold.h"

using namespace maniFEM;

	
Manifold Manifold::working ( tag::non_existent );
// anything would do, the user must set this variable before anything else
// by simply declaring any Manifold object (constructor assigns to Manifold::working)

	
Function Manifold::Euclid::build_coord_func
( const tag::lagrange &, const tag::OfDegree &, size_t deg )
// virtual from Manifold::Core
	
{	assert ( deg == 1 );
	Field::Core * field;
	Function::Core * func;
	if ( this->dim == 1 )
	{	Field::Scalar * field_scalar = new Field::Scalar
			( tag::lives_on_positive_cells, tag::of_dim, 0 );
		field = field_scalar;
		func = new Function::CoupledWithField::Scalar ( field_scalar );  }
	else
	{	assert ( this->dim > 1 );
		Field::Block * field_block = new Field::Block
			( tag::lives_on_positive_cells, tag::of_dim, 0, tag::has_size, this->dim );
		field = field_block;
	  func = new Function::CoupledWithField::Vector ( field_block );                }
	this->coord_field = field;
	this->coord_func = Function ( tag::whose_core_is, func );
	return this->coord_func;                                                            }


Function Manifold::Implicit::build_coord_func
( const tag::lagrange &, const tag::OfDegree &, size_t deg )
// virtual from Manifold::Core
{	assert ( false );  }

Function Manifold::Parametric::build_coord_func
( const tag::lagrange &, const tag::OfDegree &, size_t deg )
// virtual from Manifold::Core
{	assert ( false );  }

	
Function Manifold::Euclid::get_coord_func ( ) const  // virtual from Manifold::Core
{	return this->coord_func;  }
	
Function Manifold::Implicit::get_coord_func ( ) const  // virtual from Manifold::Core
{	return this->surrounding_space.coordinates();  }

Function Manifold::Parametric::get_coord_func ( ) const  // virtual from Manifold::Core
{	return this->surrounding_space.coordinates();  }


void Manifold::Euclid::set_coords ( const Function co )  // virtual from Manifold::Core
{	this->coord_func = co;  }

void Manifold::Implicit::set_coords ( const Function co )  // virtual from Manifold::Core
{	Manifold m = this->surrounding_space;
	Manifold::Euclid * m_euclid = dynamic_cast<Manifold::Euclid*> ( m.core );
	assert ( m_euclid );
	m_euclid->coord_func = co;                                                 }

void Manifold::Parametric::set_coords ( const Function co )  // virtual from Manifold::Core
{	Manifold m = this->surrounding_space;
	Manifold::Euclid * m_euclid = dynamic_cast<Manifold::Euclid*> ( m.core );
	assert ( m_euclid );
	m_euclid->coord_func = co;                                                 }


// metric in the manifold (an inner product on the tangent space)
double Manifold::default_inner_prod  // static
(	const Cell & P, const std::vector<double> & v,
	const std::vector<double> & w, const Function & )

{	size_t dim = v.size();
	assert ( dim == w.size() );
	double res = 0.;
	for ( size_t i = 0; i < v.size(); i++ )  res += v[i]*w[i];
	return res;                                                  }


// metric in the manifold (an inner product on the tangent space) not constant
double Manifold::zoom_inner_prod  // static
(	const Cell & P, const std::vector<double> & v,
	const std::vector<double> & w, const Function & metric )

{	size_t dim = v.size();
	assert ( dim == w.size() );
	assert ( metric.nb_of_components() == 1 );
	// std::cout << "zoom (" << metric(P) <<") ";
	double res = 0.;
	for ( size_t i = 0; i < v.size(); i++ )  res += v[i]*w[i];
	return res * metric(P);                                     }


// metric in the manifold (an inner product on the tangent space) not constant, anisotropic
double Manifold::matrix_inner_prod  // static
(	const Cell & P, const std::vector<double> & v,
	const std::vector<double> & w, const Function & metric )

{	size_t dim = v.size();
	assert ( dim == w.size() );
	assert ( metric.nb_of_components() == dim*dim );
	double res = 0.;
	for ( size_t i = 0; i < dim; i++ )
	for ( size_t j = 0; j < dim; j++ )
		res += v[i]*w[j] * metric[dim*i+j](P);
	return res;                                        }


// P = sA + sB,  s+t == 1

void Manifold::Euclid::pretty_interpolate
( const Cell & P, double s, const Cell & A, double t, const Cell & B ) const

{	Function coord = this->get_coord_func();
	size_t n = coord.nb_of_components();
	for ( size_t i = 0; i < n; i++ )
		coord[i](P) = s * coord[i](A) + t * coord[i](B);  }	

// code below does the same as code above
// above is pretty, slightly slower - below is ugly, slightly faster


void Manifold::Euclid::interpolate ( Cell::Positive::Vertex * P,
  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B ) const
//  virtual from Manifold::Core

// we could inline these, as interpolate_euclid, to gain speed	
	
{	Function coord = this->get_coord_func();
	Function::Scalar * coord_scalar = dynamic_cast < Function::Scalar * > ( coord.core );
	if ( coord_scalar )
	{	assert ( coord.nb_of_components() == 1 );
	  coord_scalar->set_value_on_cell
			( P, s * coord_scalar->get_value_on_cell(A) +
			     t * coord_scalar->get_value_on_cell(B)   );
		return;                                            } 
	Function::Vector * coord_vector = dynamic_cast < Function::Vector * > ( coord.core );
	assert ( coord_vector );
	size_t n = coord.nb_of_components();
	for ( size_t i = 0; i < n; i++ )
	{	Function::Scalar * coord_i = Mesh::assert_cast
			< Function::Core*, Function::Scalar* > ( coord_vector->component(i).core );
		coord_i->set_value_on_cell
			( P, s * coord_i->get_value_on_cell(A) +
			     t * coord_i->get_value_on_cell(B) );             }                               }


// P = sA + sB + uC + vD,  s+t+u+v == 1

void Manifold::Euclid::pretty_interpolate
( const Cell & P, double s, const Cell & A, double t, const Cell & B,
                  double u, const Cell & C, double v, const Cell & D ) const

{	Function coord = this->get_coord_func();
	size_t n = coord.nb_of_components();
	for ( size_t i = 0; i < n; i++ )
		coord[i](P) = s * coord[i](A) + t * coord[i](B) +
		              u * coord[i](C) + v * coord[i](D);   }

// code below does the same as code above
// above is pretty, slightly slower - below is ugly, slightly faster


void Manifold::Euclid::interpolate ( Cell::Positive::Vertex * P,
	double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D ) const
//  virtual from Manifold::Core

// we could inline these, as interpolate_euclid, to gain speed	
	
{	Function coord = this->get_coord_func();
	Function::Scalar * coord_scalar = dynamic_cast < Function::Scalar * > ( coord.core );
	if ( coord_scalar )
	{	assert ( coord.nb_of_components() == 1 );
	  coord_scalar->set_value_on_cell
			( P, s * coord_scalar->get_value_on_cell(A) +
			     t * coord_scalar->get_value_on_cell(B) +
			     u * coord_scalar->get_value_on_cell(C) +
			     v * coord_scalar->get_value_on_cell(D)   );
		return;                                                        } 
	Function::Vector * coord_vector = dynamic_cast < Function::Vector * > ( coord.core );
	assert ( coord_vector );
	size_t n = coord.nb_of_components();
	for ( size_t i = 0; i < n; i++ )
	{	Function::Scalar * coord_i = Mesh::assert_cast
			< Function::Core*, Function::Scalar* > ( coord_vector->component(i).core );
		coord_i->set_value_on_cell
			( P, s * coord_i->get_value_on_cell(A) +
			     t * coord_i->get_value_on_cell(B) +
			     u * coord_i->get_value_on_cell(C) +
			     v * coord_i->get_value_on_cell(D)   );           }                             }


// P = sA + sB + uC + vD + wE + zF,  s+t+u+v+w+z == 1

void Manifold::Euclid::pretty_interpolate
( const Cell & P, double s, const Cell & A, double t, const Cell & B,
                  double u, const Cell & C, double v, const Cell & D,
                  double w, const Cell & E, double z, const Cell & F ) const

{	Function coord = this->get_coord_func();
	size_t n = coord.nb_of_components();
	for ( size_t i = 0; i < n; i++ )
		coord[i](P) = s * coord[i](A) + t * coord[i](B) +
		              u * coord[i](C) + v * coord[i](D) +
		              w * coord[i](E) + z * coord[i](F);   }

// code below does the same as code above
// above is pretty, slightly slower - below is ugly, slightly faster


void Manifold::Euclid::interpolate ( Cell::Positive::Vertex * P,
	double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D,
	double w, Cell::Positive::Vertex * E, double z, Cell::Positive::Vertex * F ) const
//  virtual from Manifold::Core

// we could inline these, as interpolate_euclid, to gain speed	
	
{	Function coord = this->get_coord_func();
	Function::Scalar * coord_scalar = dynamic_cast < Function::Scalar * > ( coord.core );
	if ( coord_scalar )
	{	assert ( coord.nb_of_components() == 1 );
	  coord_scalar->set_value_on_cell
			( P, s * coord_scalar->get_value_on_cell(A) +
			     t * coord_scalar->get_value_on_cell(B) +
			     u * coord_scalar->get_value_on_cell(C) +
			     v * coord_scalar->get_value_on_cell(D) +
			     w * coord_scalar->get_value_on_cell(E) +
			     z * coord_scalar->get_value_on_cell(F)   );
		return;                                             } 
	Function::Vector * coord_vector = dynamic_cast < Function::Vector * > ( coord.core );
	assert ( coord_vector );
	size_t n = coord.nb_of_components();
	for ( size_t i = 0; i < n; i++ )
	{	Function::Scalar * coord_i = Mesh::assert_cast
			< Function::Core*, Function::Scalar* > ( coord_vector->component(i).core );
		coord_i->set_value_on_cell
			( P, s * coord_i->get_value_on_cell(A) +
			     t * coord_i->get_value_on_cell(B) +
				   u * coord_i->get_value_on_cell(C) +
			     v * coord_i->get_value_on_cell(D) +
			     w * coord_i->get_value_on_cell(E) +
			     z * coord_i->get_value_on_cell(F)   );            }                              }


// P = sA + sB,  s+t == 1

void Manifold::Euclid::pretty_interpolate ( const Cell & P,
	std::vector < double > & coefs, std::vector < Cell > & points ) const

{	Function coord = this->get_coord_func();
	size_t n = coord.nb_of_components();
	size_t m = points.size();  // == coefs.size()
	for ( size_t i = 0; i < n; i++ )
	{	double v = 0.;
		for ( size_t j = 0; j < m; j++ )  v += coefs[j] * coord[i]( points[j] );
		coord[i](P) = v;                                                          }	 }

// code below does the same as code above
// above is pretty, slightly slower - below is ugly, slightly faster


void Manifold::Euclid::interpolate ( Cell::Positive::Vertex * P,
	std::vector < double > & coefs, std::vector < Cell::Positive::Vertex * > & points ) const
//  virtual from Manifold::Core

// we could inline these, as interpolate_euclid, to gain speed	
	
{	Function coord = this->get_coord_func();
	Function::Scalar * coord_scalar = dynamic_cast < Function::Scalar * > ( coord.core );
	if ( coord_scalar )
	{	assert ( coord.nb_of_components() == 1 );
		double v = 0.;
		size_t m = points.size();  // m== coefs.size()
		for ( size_t j = 0; j < m; j++ )
			v += coefs[j] * coord_scalar->get_value_on_cell ( points[j] );
		coord_scalar->set_value_on_cell ( P, v );
		return;                                            } 
	Function::Vector * coord_vector = dynamic_cast < Function::Vector * > ( coord.core );
	assert ( coord_vector );
	size_t n = coord.nb_of_components(), m = points.size();  // m== coefs.size()
	for ( size_t i = 0; i < n; i++ )
	{	Function::Scalar * coord_i = Mesh::assert_cast
			< Function::Core*, Function::Scalar* > ( coord_vector->component(i).core );
		double v = 0.;
		for ( size_t j = 0; j < m; j++ )
			v += coefs[j] * coord_i->get_value_on_cell ( points[j] );
		coord_i->set_value_on_cell ( P, v );                         }                       }


// P = sA + sB,  s+t == 1     virtual from Manifold::Core
void Manifold::Implicit::interpolate ( Cell::Positive::Vertex * P,
	double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B ) const

{	this->surrounding_space.core->interpolate ( P, s, A, t, B );
	this->project ( P );                                          }


// P = sA + sB + uC + vD,  s+t == 1     virtual from Manifold::Core
void Manifold::Implicit::interpolate ( Cell::Positive::Vertex * P,
  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
  double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D ) const

{	this->surrounding_space.core->interpolate ( P, s, A, t, B, u, C, v, D );
	this->project ( P );                                                      }


// P = sA + sB + uC + vD + wE + zF,  s+t+u+v+w+z == 1     virtual from Manifold::Core
void Manifold::Implicit::interpolate ( Cell::Positive::Vertex * P,
  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
  double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D,
  double w, Cell::Positive::Vertex * E, double z, Cell::Positive::Vertex * F ) const

{	this->surrounding_space.core->interpolate ( P, s, A, t, B, u, C, v, D, w, E, z, F );
	this->project ( P );                                                                  }


// P = sum c_k P_k,  sum c_k == 1     virtual from Manifold::Core
void Manifold::Implicit::interpolate ( Cell::Positive::Vertex * P,
	std::vector < double > & coefs, std::vector < Cell::Positive::Vertex * > & points ) const

{	this->surrounding_space.core->interpolate ( P, coefs, points );
	this->project ( P );                                            }


	
// P = sA + sB,  s+t == 1     virtual from Manifold::Core
void Manifold::Parametric::interpolate ( Cell::Positive::Vertex * P,
  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B ) const

{	this->surrounding_space.core->interpolate ( P, s, A, t, B );
	this->project ( P );                                          }


// P = sA + sB + uC + vD,  s+t == 1     virtual from Manifold::Core
void Manifold::Parametric::interpolate ( Cell::Positive::Vertex * P,
  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
  double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D ) const

{	this->surrounding_space.core->interpolate ( P, s, A, t, B, u, C, v, D );
	this->project ( P );                                                      }


// P = sA + sB + uC + vD + wE + zF,  s+t+u+v+w+z == 1     virtual from Manifold::Core
void Manifold::Parametric::interpolate ( Cell::Positive::Vertex * P,
	double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D,
	double w, Cell::Positive::Vertex * E, double z, Cell::Positive::Vertex * F ) const
//  virtual from Manifold::Core

{	this->surrounding_space.core->interpolate ( P, s, A, t, B, u, C, v, D, w, E, z, F );
	this->project ( P );                                                                  }


// P = sum c_k P_k,  sum c_k == 1     virtual from Manifold::Core
void Manifold::Parametric::interpolate ( Cell::Positive::Vertex * P,
	std::vector < double > & coefs, std::vector < Cell::Positive::Vertex * > & points ) const

{	this->surrounding_space.core->interpolate ( P, coefs, points );
	this->project ( P );                                            }


void Manifold::Euclid::project ( Cell::Positive::Vertex * P_c ) const  { }
	

void Manifold::Implicit::OneEquation::project ( Cell::Positive::Vertex * P_c ) const

// just a few steps of Newton's method for under-determined equations

// we could inline this, as interpolate_implicit_one, to gain speed	
	
{	const Function & coord = this->get_coord_func();
	size_t n = coord.nb_of_components();
	const Function & lev_func = this->level_function;
	assert ( lev_func.nb_of_components() == 1 );
	const Function & grad_lev = this->grad_lev_func;
	assert ( grad_lev.nb_of_components() == n );
	const Cell P ( tag::whose_core_is, P_c );  // temporary wrapper for P_c
	for ( short int k = 0; k < Manifold::Implicit::steps_for_Newton; k++ )
	// we move doubles around a lot
	// how to do it faster ?
	// somehow bind references to coord_at_P to 'coord'
	{	std::vector<double> coord_at_P = coord ( P );
		double lev_at_P = lev_func ( P );
		std::vector<double> grad_lev_at_P = grad_lev ( P );
		double norm2 = 0.;
		for ( size_t i = 0; i < n; i++ )
		{	double & tmp = grad_lev_at_P[i];
			norm2 += tmp*tmp;                 }
		double coef = lev_at_P / norm2;
		for ( size_t i = 0; i < n; i++ )
			coord_at_P[i] -= coef * grad_lev_at_P[i];
		coord ( P ) = coord_at_P;                               }                    }


void Manifold::Implicit::TwoEquations::project ( Cell::Positive::Vertex * P_c ) const

// just a few steps of Newton's method for under-determined systems of equations

// we could inline this, as interpolate_implicit_two, to gain speed	
	
{	const Function & coord = this->get_coord_func();
	size_t n = coord.nb_of_components();
	const Function & lev_func_1 = this->level_function_1;
	assert ( lev_func_1.nb_of_components() == 1 );
	const Function & grad_lev_1 = this->grad_lev_func_1;
	assert ( grad_lev_1.nb_of_components() == n );
	const Function & lev_func_2 = this->level_function_2;
	assert ( lev_func_2.nb_of_components() == 1 );
	const Function & grad_lev_2 = this->grad_lev_func_2;
	assert ( grad_lev_2.nb_of_components() == n );
	const Cell P ( tag::whose_core_is, P_c );  // temporary wrapper for P_c
	for ( short int k = 0; k < Manifold::Implicit::steps_for_Newton; k++ )
	// we move doubles around a lot
	// how to do it faster ?
	// somehow bind references to coord_at_P to 'coord'
	{	std::vector<double> coord_at_P = coord ( P );
		double lev_1_at_P = lev_func_1 ( P );
		std::vector<double> grad_lev_1_at_P = grad_lev_1 ( P );
		double lev_2_at_P = lev_func_2 ( P );
		std::vector<double> grad_lev_2_at_P = grad_lev_2 ( P );
		double a11 = 0., a12 = 0., a22 = 0.;
		for ( size_t i = 0; i < n; i++ )
		{	double & tmp1 = grad_lev_1_at_P[i];
			double & tmp2 = grad_lev_2_at_P[i];
			a11 += tmp1*tmp1;
			a12 += tmp1*tmp2;
			a22 += tmp2*tmp2;                    }
		double det = a11*a22 - a12*a12;
		double l1 = ( - a22 * lev_1_at_P + a12 * lev_2_at_P ) / det;
		double l2 = (   a12 * lev_1_at_P - a11 * lev_2_at_P ) / det;
		for ( size_t i = 0; i < n; i++ )
			coord_at_P[i] += l1 * grad_lev_1_at_P[i] + l2 * grad_lev_2_at_P[i];
		coord ( P ) = coord_at_P;                                               }   }


void Manifold::Parametric::project ( Cell::Positive::Vertex * P_c ) const

{	std::map< Function::Core *, Function::Core * >::const_iterator it;
	for ( it = this->equations.begin(); it != this->equations.end(); it++ )
	{	Function::Scalar * coord_scalar = Mesh::assert_cast
			< Function::Core*, Function::Scalar* > ( it->first );
		Function::Scalar * expr_scalar = Mesh::assert_cast
			< Function::Core*, Function::Scalar* > ( it->second );
		coord_scalar->set_value_on_cell ( P_c, expr_scalar->get_value_on_cell(P_c) );  }  }


