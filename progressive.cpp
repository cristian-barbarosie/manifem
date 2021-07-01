
// manifem/progressive.cpp 2020.02.02

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

#include <stack>
#include "math.h"
#include <random>

#include "maniFEM.h"
#include "metric-tree.h"

namespace maniFEM { namespace tag

{	struct OrthogonalTo { };  static const OrthogonalTo orthogonal_to;         }  }

using namespace maniFEM;


// global variables and functions for this file, not visible for other object files
namespace {  // equivalent to 'static' (internal linkage)

std::map < Cell, MetricTree<Cell,Manifold::Euclid::SqDist>::Node * > node_in_cloud;
std::map < Cell, std::vector < double > > normals;
// use Cell::hook !

Cell temporary_vertex ( tag::non_existent );

size_t progress_nb_of_coords;  // dimension of the surrounding Euclidian space
// also known as "geometric dimension"

const double one_plus_tolerance = 1.15;
const double sqrt_of_075 = std::sqrt ( 0.75 );
const double half_of_sqrt_of_075 = sqrt_of_075 / 2.;
double desired_length, desired_length_sq,
       progress_long_dist, progress_long_dist_sq;

Mesh mesh_under_constr ( tag::of_dimension_one );
Mesh progress_interface ( tag::of_dimension_one );
// empty temporary meshes, dimension doesn't matter
// these variables will be set when we start the meshing process,
// in the body of the progressive Mesh constructor

bool correctly_oriented ( const Mesh msh );
bool correctly_oriented_complicated ( const Mesh msh );
void switch_orientation ( Mesh msh );

//-------------------------------------------------------------------------------------------------

inline double approx_sqrt ( double arg, const tag::Around &, double centre, double sqrt_of_centre )

// we assume all segments have length around the same value, here called 'sqrt_of_centre'
// so the result should not be too far from 'sqrt_of_centre'
	
{	double res = 0.5 / sqrt_of_centre * ( arg + centre );
	#ifndef NDEBUG
	if ( ( res < 0.25*sqrt_of_centre ) or ( res > 4.*sqrt_of_centre ) )
	{	std::cout << "bad approximation of square root" << std::endl;
		std::cout << arg << " " << centre << std::endl;
		std::cout << res << " " << sqrt_of_centre << std::endl;
		exit(1);                                                         }
	#endif  // debug
	return res;
}

//-------------------------------------------------------------------------------------------------

std::vector < double > compute_tangent_vec
( Cell start, bool check_orth, std::vector < double > given_vec )

// computes a vector tangent to Manifold::working at point 'start'

// if third argument is true, candidates will be projected onto the space orthogonal to given_vec
// given_vec must be tangent to Manifold::working at point 'start'
// and must have length approximately equal to desired_length
	
{	// Manifold::Implicit * m_impl =  dynamic_cast<Manifold::Implicit*> ( Manifold::working.core );
	// assert ( m_impl );
	// Manifold::Euclid * m_euclid = dynamic_cast<Manifold::Euclid*> ( m_impl->surrounding_space.core );
	// assert ( m_euclid );
	std::vector < double > best_tangent;  double longest_projection = 0.;
	const double len_sq = desired_length * desired_length;
	const size_t nc = Manifold::working.coordinates().nb_of_components();
	for ( size_t n = 1; n <= nc; n++ )
	{	const double coef = desired_length / std::sqrt(n);
		// we make sums of n vectors in the canonical basis with both signs
		std::vector < size_t > indices ( n+1 );
		for ( size_t i = 0; i < n; i++ )  indices[i] = i;
		indices[n] = nc;
		while ( true )
		{	std::vector < short int > signs ( n, 1. );
			while ( true )
			{	std::vector < double > tangent ( nc, 0. );
				for ( size_t i = 0; i < n; i++ )  tangent[indices[i]] = signs[i];
				// we normalize 'tangent'
				for ( size_t i = 0; i < nc; i++ ) tangent[i] *= coef;
				// we project
				for ( size_t i = 0; i < nc; i++ )
				{	Function x = Manifold::working.coordinates()[i];
					x ( temporary_vertex ) = x ( start ) + tangent[i];  }
				Manifold::working.project ( temporary_vertex );
				for ( size_t i = 0; i < nc; i++ )
				{	Function x = Manifold::working.coordinates()[i];
					tangent[i] = x ( temporary_vertex ) - x ( start );  }
				if ( check_orth )
				{	double prod = Manifold::working.inner_prod ( start, tangent, given_vec );
					double lambd = prod / len_sq;
					for ( size_t i = 0; i < nc; i++ )  tangent[i] -= lambd * given_vec[i];    }
				// we choose the longest projection
				double n2 = Manifold::working.inner_prod ( start, tangent, tangent );
				if ( n2 > longest_projection )
				{	best_tangent = tangent;  longest_projection = n2;  }
				// now change signs
				bool found = false;
				for ( int i = n-1; i >= 0; i-- )
					if ( signs[i] == 1 )
					{	found = true;  signs[i] = -1;
						for ( size_t j = i+1; j < n; j++ ) signs[j] = 1;
						break;                                            }
				if ( not found ) break;                                                  }
 			// now change indices
			bool found = false;
			for ( int i = n-1; i >= 0; i-- )
				if ( indices[i] < indices[i+1] - 1 )
				{	found = true;
					indices[i]++;
					for ( size_t j = i+1; j<n; j++ )
					{	indices[j] = indices[j-1] + 1;
						assert ( indices[j] < nc );     }
					break;                                                            }
			if ( not found ) break;                                                  }
	}  // end of  for n
	// normalize best_tangent
	double n2 = Manifold::working.inner_prod ( start, best_tangent, best_tangent );
	double norm = approx_sqrt ( n2, tag::around, len_sq, desired_length );
	double coef = desired_length / approx_sqrt ( n2, tag::around, norm*norm, norm );
	for ( size_t i = 0; i < nc; i++ )  best_tangent[i] *= coef;
	return best_tangent;

}  // end of  compute_tangent_vec

//-------------------------------------------------------------------------------------------------

inline void progress_add_point
( const Cell & P, MetricTree<Cell,Manifold::Euclid::SqDist> & cloud )

{	assert ( P.dim() == 0 );
	node_in_cloud[P] = cloud.add ( P );  }
// use hook !

//-------------------------------------------------------------------------------------------------

inline bool positive_orientation
( const Cell & A, const Cell & B, const Cell & AB, const Cell & BC )
	
{	assert ( A == AB.base().reverse() );
	assert ( B == AB.tip() );
	assert ( B == BC.base().reverse() );
	// code below is identical to part of progress_cos_sq_120
	// if you change anything, please change both; keep them identical
	assert ( normals.find(BC) != normals.end() );
	std::vector < double > e ( progress_nb_of_coords ),
		& f = normals[BC];
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		e[i] = x ( B ) - x ( A );              }
	double prod = 0.;  // scalar product e.f
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )  prod += e[i] * f[i];
	// we could have used the Riemannian product, but the sign should be the same
	// code above is identical to part of progress_cos_sq_120
	return prod < 0.;                                                                      }

//-------------------------------------------------------------------------------------------------

inline double progress_cos_sq_60
( const Cell & A, const Cell & B, const Cell & C, const Cell & AB, const Cell & BC )

// return cosine square of 180 - ABC (or 0. if wrong orientation)
// check that  cos_sq > 0.03  to ensure angle ABC is below 80 deg

{	assert ( A == AB.base().reverse() );
	assert ( B == AB.tip() );
	assert ( B == BC.base().reverse() );
	assert ( C == BC.tip() );
	// first, check orientation
	std::vector < double > e1 ( progress_nb_of_coords ),
	                       e2 ( progress_nb_of_coords );
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		e1[i] = x ( B ) - x ( A );
		e2[i] = x ( C ) - x ( B );                        }
	double prod = 0.;  // scalar product e1.e2
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )  prod += e1[i] * e2[i];
	// use Riemannian product !
	if ( prod > 0. ) return 0.;
	// use Riemannian metric below
	double norm1 = 0., norm2 = 0.;
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	norm1 += e1[i] * e1[i];
		norm2 += e2[i] * e2[i];  }
	return prod * prod / norm1 / norm2;  // cosine square
}

//-------------------------------------------------------------------------------------------------
			
inline double progress_cos_sq_120
( const Cell & A, const Cell & B, const Cell & C, const Cell & AB, const Cell & BC )

// return cosine square of 180 - ABC (or 2. if wrong orientation)
// for instance, calling function should check cos_sq < 0.671 to ensure angle ABC is below 145 deg
// or check cos_sq < -0.17365 to ensure angle ABC is below 80 deg

{	assert ( A == AB.base().reverse() );
	assert ( B == AB.tip() );
	assert ( B == BC.base().reverse() );
	assert ( C == BC.tip() );
	// first, check orientation
	// code below is identical to part of 'positive_orientation'
	// if you change anything, please change both; keep them identical
	assert ( normals.find(BC) != normals.end() );
	std::vector < double > e ( progress_nb_of_coords ),
		& f = normals[BC];
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		e[i] = x ( B ) - x ( A );                         }
	double prod = 0.;  // scalar product e.f
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )  prod += e[i] * f[i];
	// we could have used the Riemannian product, but the sign should be the same
	// code above is identical to part of 'positive_orientation'
	if ( prod > 0. ) return 2.;
	std::vector < double > e2 ( progress_nb_of_coords );
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		e2[i] = x ( C ) - x ( B );                       }
	// use Riemannian metric below
	double norm1 = 0., norm2 = 0.;  prod = 0.;
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	norm1 += e [i] * e [i];
		norm2 += e2[i] * e2[i];
		prod  += e [i] * e2[i];  }
	return prod * prod / norm1 / norm2;  // cosine square
}

//-------------------------------------------------------------------------------------------------
			
inline void improve_normal ( const Cell & A, std::vector < double > & nor )

// project 'nor' onto working manifold

{	// we normalize 'nor'
	double n2 = Manifold::working.inner_prod ( A, nor, nor );
	double norm = approx_sqrt ( n2, tag::around, desired_length_sq, desired_length );
	norm = approx_sqrt ( n2, tag::around, norm*norm, norm );
	norm = desired_length / norm;
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )  nor[i] *= norm;
	// we project on Manifold::working
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		x ( temporary_vertex ) = x(A) + nor[i];          }
	Manifold::working.project ( temporary_vertex );
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		nor[i] = x ( temporary_vertex ) - x(A);          }                               }

//-------------------------------------------------------------------------------------------------

inline void improve_normal
( const Cell & A, const Cell & B, std::vector < double > & AB_coord, std::vector < double > & nor )

// project 'nor' onto working manifold and normalize it relatively to segment AB

{	// we normalize 'nor'
	double n2 = Manifold::working.inner_prod ( A, nor, nor );
	double norm = approx_sqrt 
		( n2, tag::around, desired_length_sq, desired_length );
	norm = approx_sqrt ( n2, tag::around, norm*norm, norm );
	norm = desired_length / norm;
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )  nor[i] *= norm;
	// we project on Manifold::working
	std::vector < double > mid_seg ( progress_nb_of_coords );
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		mid_seg[i] = ( x(A) + x(B) ) / 2.;
		x ( temporary_vertex ) = mid_seg[i] + nor[i];  }
	Manifold::working.project ( temporary_vertex );
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
		nor[i] = Manifold::working.coordinates()[i] ( temporary_vertex ) - mid_seg[i];
	// we ensure again that 'nor' is orthogonal to AB_coord
	double prod = Manifold::working.inner_prod ( A, AB_coord, nor );
	n2 = Manifold::working.inner_prod ( A, AB_coord, AB_coord );
	prod /= n2;
	for ( size_t i = 0; i < progress_nb_of_coords; i++ ) nor[i] -= prod * AB_coord[i];
	// we normalize 'nor'
	n2 = Manifold::working.inner_prod ( A, nor, nor );
	norm = approx_sqrt 
		( n2, tag::around, desired_length_sq, desired_length );
	norm = approx_sqrt ( n2, tag::around, norm*norm, norm );
	norm = desired_length / norm;
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )  nor[i] *= norm;
}

//-------------------------------------------------------------------------------------------------

inline bool opposite_signs ( const double a, const double b )
{	if ( a < 0. )  return b >= 0.;
	if ( a == 0. )  return true;
	return b <= 0.;                  }

//-------------------------------------------------------------------------------------------------

inline Cell search_start_ver_c1 ( )

// search for a starting point in a manifold of co-dimension one
// e.g. a curve in the plane or a surface in 3D

{	Manifold::Implicit::OneEquation * m_impl =
		dynamic_cast<Manifold::Implicit::OneEquation*> ( Manifold::working.core );
	assert ( m_impl );
	// Manifold::Euclid * m_euclid = dynamic_cast<Manifold::Euclid*> ( m_impl->surrounding_space.core );
	// assert ( m_euclid );
	// Function m_impl->level_function;

	const size_t nc = Manifold::working.coordinates().nb_of_components();
	const double len_sq = desired_length * desired_length;
	
	Cell tmp_ver ( tag::vertex );
	Cell tmp_ver_1 ( tag::vertex );
	Cell tmp_ver_2 ( tag::vertex );

	std::default_random_engine random_generator;
	size_t size_of_cube = 5;
	while ( true )
	{	double s = size_of_cube * desired_length;
		std::uniform_real_distribution<double> distr ( -s, s );
		size_t nb = 1;
		for ( size_t j = 0; j < nc; j++ ) nb *= size_of_cube;
		for ( size_t j = 0; j < nb; j++ )
		{	for ( size_t i = 0; i < nc; i++ )
			{	Function x = Manifold::working.coordinates()[i];
				x ( tmp_ver_1 ) = distr(random_generator);
				x ( tmp_ver_2 ) = distr(random_generator);        }
			double v1 = m_impl->level_function ( tmp_ver_1 ),
			       v2 = m_impl->level_function ( tmp_ver_2 );
			if ( opposite_signs ( v1, v2 ) )
				// refine by applying bissection algorithm
				while ( true )
				{	if ( Manifold::working.dist_sq ( tmp_ver_1, tmp_ver_2 ) < len_sq )
					{	tmp_ver.dispose();  tmp_ver_2.dispose();
						Manifold::working.project ( tmp_ver_1 );
						return tmp_ver_1;                          }
					m_impl->surrounding_space.interpolate ( tmp_ver, 0.5, tmp_ver_1, 0.5, tmp_ver_2 );
					double v = m_impl->level_function ( tmp_ver );
					if ( opposite_signs ( v, v2 ) )
					{	for ( size_t i = 0; i < nc; i++ )
						{	Function x = Manifold::working.coordinates()[i];
							x ( tmp_ver_1 ) = x ( tmp_ver );                 }
						v1 = v;                                               }
					else if ( opposite_signs ( v, v1 ) )
					{	for ( size_t i = 0; i < nc; i++ )
						{	Function x = Manifold::working.coordinates()[i];
							x ( tmp_ver_2 ) = x ( tmp_ver );                 }
						v2 = v;                                               }
					else assert( false );                                                }  }
		size_of_cube *= 2;                                                                                 }
}

//-------------------------------------------------------------------------------------------------

inline double ext_prod_R2 ( const double vx, const double vy, const double wx, const double wy )
{	return vx*wy - wx*vy;  }


//-------------------------------------------------------------------------------------------------

inline bool origin_outside ( const double & Ax, const double & Ay,
                             const double & Bx, const double & By,
                             const double & Cx, const double & Cy )

{	return
	   opposite_signs ( - ext_prod_R2 ( Bx-Ax, By-Ay, Ax, Ay ), ext_prod_R2 ( Bx-Ax, By-Ay, Cx-Ax, Cy-Ay ) )
	or opposite_signs ( - ext_prod_R2 ( Cx-Bx, Cy-By, Bx, By ), ext_prod_R2 ( Cx-Bx, Cy-By, Ax-Bx, Ay-By ) )
	or opposite_signs ( - ext_prod_R2 ( Ax-Cx, Ay-Cy, Cx, Cy ), ext_prod_R2 ( Ax-Cx, Ay-Cy, Bx-Cx, By-Cy ) ); }

//-------------------------------------------------------------------------------------------------

inline Cell search_start_ver_c2 ( )

// search for a starting point in a manifold of co-dimension two (a curve in 3D)

{	Manifold::Implicit::TwoEquations * m_impl =
		dynamic_cast<Manifold::Implicit::TwoEquations*> ( Manifold::working.core );
	assert ( m_impl );
	// Manifold::Euclid * m_euclid = dynamic_cast<Manifold::Euclid*> ( m_impl->surrounding_space.core );
	// assert ( m_euclid );

	const size_t nc = Manifold::working.coordinates().nb_of_components();
	const double len_sq = desired_length * desired_length;
	
	Cell tmp_A ( tag::vertex );
	Cell tmp_B ( tag::vertex );
	Cell tmp_C ( tag::vertex );
	Cell tmp_AB ( tag::vertex );
	Cell tmp_BC ( tag::vertex );
	Cell tmp_CA ( tag::vertex );

	std::default_random_engine random_generator;
	size_t counter = 0;
	size_t size_of_cube = 5;
	while ( true )
	{	double s = size_of_cube * desired_length;
		std::uniform_real_distribution<double> distr ( -s, s );
		size_t nb = 1;
		restart :
		for ( size_t j = 0; j < nc; j++ ) nb *= size_of_cube;
		for ( size_t j = 0; j < nb; j++ )
		{	for ( size_t i = 0; i < nc; i++ )
			{	Function x = Manifold::working.coordinates()[i];
				x ( tmp_A ) = distr(random_generator);
				x ( tmp_B ) = distr(random_generator);
				x ( tmp_C ) = distr(random_generator);             }
			counter++;
			double vA1 = m_impl->level_function_1 ( tmp_A ),
			       vA2 = m_impl->level_function_2 ( tmp_A );
			double vB1 = m_impl->level_function_1 ( tmp_B ),
			       vB2 = m_impl->level_function_2 ( tmp_B );
			double vC1 = m_impl->level_function_1 ( tmp_C ),
			       vC2 = m_impl->level_function_2 ( tmp_C );
			if ( not origin_outside ( vA1, vA2, vB1, vB2, vC1, vC2 ) )
				// refine by repeatedly cutting the triangle
				while ( true )
				{	if ( ( Manifold::working.dist_sq ( tmp_A, tmp_B ) < len_sq ) and
					     ( Manifold::working.dist_sq ( tmp_B, tmp_C ) < len_sq ) and
					     ( Manifold::working.dist_sq ( tmp_C, tmp_A ) < len_sq )      )
					{	tmp_B.dispose();  tmp_C.dispose();
						tmp_AB.dispose();  tmp_BC.dispose();  tmp_CA.dispose();
						Manifold::working.project ( tmp_A );
						return tmp_A;                                            }
					m_impl->surrounding_space.interpolate ( tmp_AB, 0.5, tmp_A, 0.5, tmp_B );
					double vAB1 = m_impl->level_function_1 ( tmp_AB ),
					       vAB2 = m_impl->level_function_2 ( tmp_AB );
				  m_impl->surrounding_space.interpolate ( tmp_BC, 0.5, tmp_B, 0.5, tmp_C );
					double vBC1 = m_impl->level_function_1 ( tmp_BC ),
					       vBC2 = m_impl->level_function_2 ( tmp_BC );
				  m_impl->surrounding_space.interpolate ( tmp_CA, 0.5, tmp_C, 0.5, tmp_A );
					double vCA1 = m_impl->level_function_1 ( tmp_CA ),
					       vCA2 = m_impl->level_function_2 ( tmp_CA );
					if ( not origin_outside ( vA1, vA2, vAB1, vAB2, vCA1, vCA2 ) )
					{	for ( size_t i = 0; i < nc; i++ )
						{	Function x = Manifold::working.coordinates()[i];
							x ( tmp_B ) = x ( tmp_AB );  x ( tmp_C ) = x ( tmp_CA );  }
						vB1 = vAB1;  vB2 = vAB2;  vC1 = vCA1;  vC2 = vCA2;            }
					else if ( not origin_outside ( vB1, vB2, vAB1, vAB2, vBC1, vBC2 ) )
					{	for ( size_t i = 0; i < nc; i++ )
						{	Function x = Manifold::working.coordinates()[i];
							x ( tmp_A ) = x ( tmp_AB );  x ( tmp_C ) = x ( tmp_BC );  }
						vA1 = vAB1;  vA2 = vAB2;  vC1 = vBC1;  vC2 = vBC2;            }
					else if ( not origin_outside ( vC1, vC2, vBC1, vBC2, vCA1, vCA2 ) )
					{	for ( size_t i = 0; i < nc; i++ )
						{	Function x = Manifold::working.coordinates()[i];
							x ( tmp_A ) = x ( tmp_CA );  x ( tmp_B ) = x ( tmp_BC );  }
						vA1 = vCA1;  vA2 = vCA2;  vB1 = vBC1;  vB2 = vBC2;            }
					else if ( not origin_outside ( vAB1, vAB2, vBC1, vBC2, vCA1, vCA2 ) )
					{	for ( size_t i = 0; i < nc; i++ )
						{	Function x = Manifold::working.coordinates()[i];
							x ( tmp_A ) = x ( tmp_BC );  x ( tmp_B ) = x ( tmp_CA );
							x ( tmp_C ) = x ( tmp_AB );                               }
						vA1 = vBC1;  vA2 = vBC2;  vB1 = vCA1;  vB2 = vCA2;
						vC1 = vAB1;  vC2 = vAB2;                                                 }
					else  // nasty nonlinear level functions ...
						goto restart;                                                               }  }
		size_of_cube *= 2;                                                                                 }
}

//-------------------------------------------------------------------------------------------------

inline Cell search_start_ver ( const double & length )

// search for a starting point
// current working manifold may have co-dimension one
// e.g. a curve in the plane or a surface in 3D
// of co-dimension two (a curve in 3D)

{	desired_length = length;
	Manifold::Implicit::OneEquation * m_impl =
		dynamic_cast<Manifold::Implicit::OneEquation*> ( Manifold::working.core );
	if ( m_impl ) return search_start_ver_c1 ( );  // co-dimension one
	return search_start_ver_c2 ( );  // co-dimension two
}

//-------------------------------------------------------------------------------------------------

inline void redistribute_vertices ( const Mesh & msh,
  const Cell & start, const Cell & stop, double last_length, size_t n )

// chain of n segments, desired length d, last length d'
// move p1 with   (d'-d) / n
//    alpha x0 + beta x2 == x1 + (d'-d) / n   alpha+beta == 1
//    x2 == x0 + 2d   x1 == x0 + d
//    2 d beta == d + (d'-d)/ n 
// move p2 with 2*(d'-d) / n
//    alpha x1 + beta x3 == x2 + 2(d'-d)/n   alpha+beta == 1
//    x1 == x0 + d + (d'-d)/n
//    x3 == x0 + 3d
//    alpha (d'-d)/n + 2 d beta == d + 2(d'-d)/n
//    ( 2d - (d'-d)/n ) beta == d + (d'-d)/n
// move p2 with 3*(d'-d) / n
//    ( 2d - 2(d'-d)/n ) beta == d + (d'-d)/n

// how to do this if the metric is not uniform ?
// use sums of lengths os segments, mimiking geodesics
	
{	Cell A = msh.cell_behind(stop).base().reverse();
	for ( size_t i = 1; i < n; i++ )
	{	if ( A == start )  {  n = i;  break;  }
		A = msh.cell_behind(A).base().reverse();  }
	assert ( n > 1 );
	Cell B = msh.cell_in_front_of(A).tip();
	Cell C = msh.cell_in_front_of(B).tip();
	double epsilon = ( last_length - desired_length ) / n;
	double v1 = desired_length + epsilon,  v2 = 2*desired_length;
	while ( C != stop )
	{	double beta = v1/v2;
		Manifold::working.interpolate ( B, 1.-beta, A, beta, C );
		v2 -= epsilon;  assert ( v2 > 0. );
		A = B;  B = C;
		C = msh.cell_in_front_of(B).tip();                            }
	Manifold::working.interpolate ( B, 0.5, A, 0.5, C );               }

//-------------------------------------------------------------------------------------------------

inline double get_z_baric ( const Cell & tri )

{	assert ( Manifold::working.coordinates().nb_of_components() == 3 );
	CellIterator it = tri.boundary().iter_over ( tag::vertices );
	Function z = Manifold::working.coordinates()[2];
	double zz = 0.;
	size_t counter = 0;
	for ( it.reset(); it.in_range(); it++, counter++ )  zz += z(*it);
	assert ( counter == 3 );
	return  zz/3.;                                                               }

//-------------------------------------------------------------------------------------------------

inline bool tri_correctly_oriented ( const Cell & tri )

{	assert ( tri.dim() == 2 );

	CellIterator it = tri.boundary().iter_over ( tag::segments );
	it.reset();  assert ( it.in_range() );
	Cell AB = *it;
	it++;  assert ( it.in_range() );
	Cell BC = *it;
	it++;  assert ( it.in_range() );
	// Cell CA = *it;
	it++;  assert ( not it.in_range() );

	Cell A = AB.base().reverse();
	Cell B = AB.tip();
	assert ( B == BC.base().reverse() );
	Cell C = BC.tip();

	Function x = Manifold::working.coordinates()[0];
	Function y = Manifold::working.coordinates()[1];
	double  xAB = x(B) - x(A),  yAB = y(B) - y(A),
	        xBC = x(C) - x(B),  yBC = y(C) - y(B);
	return  xAB * yBC > yAB * xBC;;                                  }

//-------------------------------------------------------------------------------------------------

bool correctly_oriented ( const Mesh msh )

// tells whether 'msh's orientation is consistent with the orientation of the
// surrounding Euclidian space

{	Manifold::Implicit::OneEquation * m_impl =
		dynamic_cast<Manifold::Implicit::OneEquation*> ( Manifold::working.core );
	assert ( m_impl );  // more explicit error message : did you forget the tag::random_orientation ?
	Manifold::Euclid * m_euclid = dynamic_cast<Manifold::Euclid*> ( m_impl->surrounding_space.core );
	assert ( m_euclid );
	const size_t nc = Manifold::working.coordinates().nb_of_components();

	// for surfaces, we should search the vertex with zmax and check the orientation
	// of all surrounding triangles (we need an iterator over cells above)
	if ( msh.dim() != 1 )
	{	assert ( msh.dim() == 2 );
		assert ( nc == 3 );
		CellIterator it = msh.iter_over ( tag::cells_of_dim, 2 );
		it.reset();  assert ( it.in_range() );
		Cell trimax = *it;
		double zmax = get_z_baric ( trimax );
		for ( it++; it.in_range(); it++ )
		{	double zz = get_z_baric ( *it );
			if ( zz > zmax )
			{	zmax = zz;  trimax = *it;  }   }
		return tri_correctly_oriented ( trimax );                   }

	assert ( msh.dim() == 1 );
	Function x = Manifold::working.coordinates()[0];
	Function y = Manifold::working.coordinates()[1];

	CellIterator it = msh.iter_over ( tag::vertices );
	it.reset();  assert ( it.in_range() );
	Cell ver = *it;
	double ymax = y(ver);
	for ( it++; it.in_range(); it++ )
	{	Cell other_ver = *it;
		double other_y = y(other_ver);
		if ( other_y > ymax )
		{	ymax = other_y;  ver = other_ver;  }  }
	Cell prev_seg = msh.cell_behind ( ver );
	Cell next_seg = msh.cell_in_front_of ( ver );
	assert ( prev_seg.tip() == ver );
	assert ( next_seg.base().reverse() == ver );
	Cell A = prev_seg.base().reverse();
	Cell C = next_seg.tip();
	bool prev_orient = ( x(ver) < x(A) );
	bool next_orient = ( x(C) < x(ver) );
	if ( prev_orient == next_orient )  return prev_orient;
	double  xAB = x(ver) - x(A),  yAB = y(ver) - y(A),
	        xBC = x(C) - x(ver),  yBC = y(C) - y(ver);
	return  xAB * yBC > yAB * xBC;;
	
}

//-------------------------------------------------------------------------------------------------

bool correctly_oriented_complicated ( const Mesh msh )

// tells whether 'msh's orientation is consistent with the orientation of the
// surrounding Euclidian space

{	Manifold::Implicit::OneEquation * m_impl =
		dynamic_cast<Manifold::Implicit::OneEquation*> ( Manifold::working.core );
	assert ( m_impl );
	Manifold::Euclid * m_euclid = dynamic_cast<Manifold::Euclid*> ( m_impl->surrounding_space.core );
	assert ( m_euclid );
	const size_t nc = Manifold::working.coordinates().nb_of_components();

	if ( nc != 2 )
	{	std::cout << "for the moment I can only check the orientation of "
	            << "(closed) curves in the plane - sorry" << std::endl;
		exit ( 1 );                                                         }

	Function x = Manifold::working.coordinates()[0];
	Function y = Manifold::working.coordinates()[1];
	CellIterator it = msh.iter_over ( tag::segments );
	it.reset();  assert ( it.in_range() );
	Cell seg = *it;
	assert ( seg.dim() == 1 );
	bool dx_pos = ( x (seg.tip()) - x (seg.base().reverse()) ) > 0.;
	bool dy_pos = ( y (seg.tip()) - y (seg.base().reverse()) ) > 0.;
	int counter = 0;
	for ( it++; it.in_range(); it++ )
	{	seg = *it;
		bool dx_pos_now = ( x(seg.tip()) - x (seg.base().reverse()) ) > 0.;
		bool dy_pos_now = ( y(seg.tip()) - y (seg.base().reverse()) ) > 0.;
		if ( dx_pos_now != dx_pos )  //  dx has changed sign
		{	if ( dy_pos_now != dy_pos )  // dy has changed sign, too
			{	std::cout << "I cannot check the orientation if the curve "
			            << "has too sharp angles - sorry" << std::endl;
				// we can do better here
				exit ( 1 );                                                  }
			if ( dx_pos == dy_pos ) counter++;
			else counter--;                                                    }
		dx_pos = dx_pos_now;  dy_pos = dy_pos_now;                              }
	assert ( ( counter == 2 ) or ( counter == -2 ) );
	return counter == 2;                                                         }

//-------------------------------------------------------------------------------------------------

inline void switch_orientation ( Cell cll )

// this is always called from switch_orientation ( Mesh )
// we should deal with segments separately, as well as with negative cells !
	
{	Mesh msh = cll.boundary();
	std::vector < Cell > vec_of_cells;
	vec_of_cells.reserve ( msh.number_of ( tag::cells_of_dim, msh.dim() ) );
	CellIterator itt = msh.iter_over ( tag::cells_of_dim, msh.dim(), tag::force_positive );
	for ( itt.reset(); itt.in_range(); itt++ )
		vec_of_cells.push_back ( *itt );
	for ( std::vector<Cell>::iterator it = vec_of_cells.begin(); it != vec_of_cells.end(); it++ )
		it->remove_from ( msh );
	for ( std::vector<Cell>::iterator it = vec_of_cells.begin(); it != vec_of_cells.end(); it++ )
		it->reverse().add_to ( msh );                                                                 }

//-------------------------------------------------------------------------------------------------

void switch_orientation ( Mesh msh )

// do not use reverse !
// call switch_orientation on each cell

{	std::vector < Cell > vec_of_cells;
	vec_of_cells.reserve ( msh.number_of ( tag::cells_of_dim, msh.dim() ) );
	CellIterator itt = msh.iter_over ( tag::cells_of_dim, msh.dim() );
	for ( itt.reset(); itt.in_range(); itt++ )
		vec_of_cells.push_back ( *itt );
	for ( std::vector<Cell>::iterator it = vec_of_cells.begin(); it != vec_of_cells.end(); it++ )
		it->remove_from ( msh );
	for ( std::vector<Cell>::iterator it = vec_of_cells.begin(); it != vec_of_cells.end(); it++ )
		it->reverse().add_to ( msh );                                                                 }

//-------------------------------------------------------------------------------------------------

inline void build_one_normal ( Cell & B, Cell & C, Cell & new_seg )

// builds the normal vector for 'new_seg', based on information from previous segment
	
{	assert ( B == new_seg.base().reverse() );
	assert ( C == new_seg.tip() );
	Cell AB = progress_interface.cell_behind ( B );
	
	assert ( AB.tip() == B );
	Cell A = AB.base().reverse();

	std::vector < double > vA = Manifold::working.coordinates() ( A );
	std::vector < double > vB = Manifold::working.coordinates() ( B );
	std::vector < double > old_e ( progress_nb_of_coords );
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )  old_e[i] = vB[i] - vA[i];
	assert ( normals.find(AB) != normals.end() );
	std::vector < double > & old_f = normals[AB];
	// 'e' is the vector of the segment, 'f' is orthogonal
	// they are all of approximately the same length, equal to desired_length

	// code below is identical to part of 'build_each_normal'
	// if you change anything, please change both; keep them identical
	std::vector < double > vC = Manifold::working.coordinates() ( C );
	std::vector < double > new_e ( progress_nb_of_coords );
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )  new_e[i] = vC[i] - vB[i];
	// scalar products :
	double with_e = 0.,  with_f = 0., norm_e_sq = 0., norm_f_sq = 0.;
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )  // use Riemannian metric !
	{	with_e += new_e[i]*old_e[i];     with_f += new_e[i]*old_f[i];
		norm_e_sq += old_e[i]*old_e[i];  norm_f_sq += old_f[i]*old_f[i];  }
	with_e /= norm_e_sq;
	with_f /= norm_f_sq;
	std::vector < double > new_f ( progress_nb_of_coords );
	// we rotate 'new_e' with 90 degrees, in the same sense as 'f' is rotated from 'e'
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
		new_f[i] = - with_f * old_e[i] + with_e * old_f[i];
	// project and normalize :
	improve_normal ( B, C, new_e, new_f );
	normals[new_seg] = new_f;  // optimize
	// code above is identical to part of 'build_each_normal'
}  // end of build_one_normal

//-------------------------------------------------------------------------------------------------

inline void build_each_normal
( Cell & B, Cell & C, Cell & new_seg, std::vector < double > & old_e, std::vector < double > & old_f )

// 'e' is the vector of the segment, 'f' is orthogonal
// they are all of approximately the same length, equal to desired_length
// the computed normal is 'new_f'
// at the end the old vectors 'e' and 'f' are replaced by the new ones

{	assert ( B == new_seg.base().reverse() );
	assert ( C == new_seg.tip() );
	std::vector < double > vB = Manifold::working.coordinates() ( B );
	// code below is identical to part of 'build_one_normal'
	// if you change anything, please change both; keep them identical
	std::vector < double > vC = Manifold::working.coordinates() ( C );
	std::vector < double > new_e ( progress_nb_of_coords );
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )  new_e[i] = vC[i] - vB[i];
	// scalar products :
	double with_e = 0.,  with_f = 0.;
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )  // use Riemannian metric !
	{	with_e += new_e[i]*old_e[i];  with_f += new_e[i]*old_f[i];  }
	with_e /= desired_length_sq;
	with_f /= desired_length_sq;
	std::vector < double > new_f ( progress_nb_of_coords );
	// we rotate 'new_e' with 90 degrees, in the same sense as 'f' is rotated from 'e'
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
		new_f[i] = - with_f * old_e[i] + with_e * old_f[i];
	// project and normalize :
	improve_normal ( B, C, new_e, new_f );
	normals[new_seg] = new_f;  // optimize
	// code above is identical to part of 'build_one_normal'
	old_e = new_e;  old_f = new_f;
}  // end of build_each_normal

//-------------------------------------------------------------------------------------------------

inline Cell build_normals ( const Cell & start )

// from a cell 'start', propagate normals along progress_interface
// (will only cover the connected component containing 'start')
// return the first segment which already has a normal
// see paragraph 9.6 in the manual
// 'normal' should have norm approximately equal to desired_length

{	assert ( start.belongs_to ( progress_interface, tag::oriented ) );
	std::cout << "building normals" << std::endl;
	Cell seg = start;
	Cell A = seg.base().reverse(),  B = seg.tip();
	std::vector < double > va = Manifold::working.coordinates() (A),
		vb = Manifold::working.coordinates() (B);
	assert ( normals.find ( seg ) != normals.end() );
	std::vector < double > e ( progress_nb_of_coords ),
		f = normals[seg];  // assert key found
	for ( size_t i = 0; i < progress_nb_of_coords; i++ ) e[i] = vb[i] - va[i];
	// 'e' and 'f' form an oriented basis in the two-dimensional space
	// tangent to the manifold at the current point
	// they will be used to build further vectors pointing outwards
	// (on the correct side of progress_interface)
	while ( true )
	// progress_interface may be disconnected, so we cannot use CellIterators
	// this loop will only cover its current connected component
	{ Cell new_seg = progress_interface.cell_in_front_of ( B );
		if ( normals.find(new_seg) != normals.end() )  return new_seg;
		assert ( new_seg != start );
		Cell C = new_seg.tip();
		build_each_normal ( B, C, new_seg, e, f );
		// 'e' and 'f' get updated within 'build_each_normal'
		seg = new_seg;  B = C;                                                       }
	
}  // end of build_normals

//-------------------------------------------------------------------------------------------------

inline void progress_fill_60
(	Cell & AB, Cell & BC, const Cell & CA, const Cell & B,
	MetricTree<Cell,Manifold::Euclid::SqDist> & cloud     )

{	AB.remove_from ( progress_interface );
	BC.remove_from ( progress_interface );
	assert ( normals.find(AB) != normals.end() );
	assert ( normals.find(BC) != normals.end() );
	normals.erase ( AB );
	normals.erase ( BC );
	cloud.remove ( node_in_cloud[B] );
	Cell new_tri ( tag::triangle, AB, BC, CA );
	
	std::vector < double > vA = Manifold::working.coordinates() ( CA.tip() );
	std::vector < double > vB = Manifold::working.coordinates() ( B );
	std::vector < double > vC = Manifold::working.coordinates() ( BC.tip() );

	new_tri.add_to ( mesh_under_constr );
	mesh_under_constr.baricenter ( B, AB );
}

//-------------------------------------------------------------------------------------------------

inline void glue_two_segs_common
( Cell & A, Cell & B, Cell & AB, Cell & C, Cell & D, Cell & CD, Cell & AD, Cell & BC )

{	AB.remove_from ( progress_interface );
	CD.remove_from ( progress_interface );
	Cell CB = BC.reverse();
	AD.add_to ( progress_interface );
	CB.add_to ( progress_interface );
	assert ( normals.find(AB) != normals.end() );
	normals.erase ( AB );  // assert that key was found
	assert ( normals.find ( CD ) != normals.end() );
	normals.erase ( CD );
	// build normals for two newly added segments AD and CB
	build_one_normal ( A, D, AD );  // from previous segment
	build_one_normal ( C, B, CB );  // from previous segment
}  // end of void glue_two_segs_common

//-------------------------------------------------------------------------------------------------

inline Cell glue_two_segs_S
(	Cell & A, Cell & B, Cell & AB, Cell & C, Cell & D, Cell & CD,
	MetricTree<Cell,Manifold::Euclid::SqDist> & cloud             )

// see paragraph 9.8 in the manual

// progress_interface.cell_in_front_of(B) may have tip C
// that is, BC may belong already to 'progress_progress_interface'

{ Cell AD ( tag::segment, A.reverse(), D );
	Cell AC ( tag::segment, A.reverse(), C );
	Cell BC = progress_interface.cell_in_front_of(B);
	if ( BC.tip() == C )
	{	progress_fill_60 ( AB, BC, AC.reverse(), B, cloud );
		AC.add_to ( progress_interface );
		build_one_normal ( A, C, AC );  // based on previous segment
		progress_fill_60 ( AC, CD, AD.reverse(), C, cloud );
		AD.add_to ( progress_interface );
		build_one_normal ( A, D, AD );  }  // based on previous segment
	else
	{	BC = Cell ( tag::segment, B.reverse(), C );
		Cell ABC ( tag::triangle, AB, BC, AC.reverse() );
		Cell ACD ( tag::triangle, AC, CD, AD.reverse() );
		ABC.add_to ( mesh_under_constr );
		ACD.add_to ( mesh_under_constr );
		glue_two_segs_common ( A, B, AB, C, D, CD, AD, BC );  }
	return AD;                                                                    }

//-------------------------------------------------------------------------------------------------

inline Cell glue_two_segs_Z
(	Cell & A, Cell & B, Cell & AB, Cell & C, Cell & D, Cell & CD,
	MetricTree<Cell,Manifold::Euclid::SqDist> & cloud      )

// see paragraph 9.8 in the manual

// progress_interface.cell_behind(A) may have base D
// that is, DA may belong already to 'progress_interface'

{	Cell BC ( tag::segment, B.reverse(), C );
	Cell DB ( tag::segment, D.reverse(), B );
	Cell DA = progress_interface.cell_behind(A);
	if ( DA.base().reverse() == D )
	{	progress_fill_60 ( DA, AB, DB.reverse(), A, cloud );
		DB.add_to ( progress_interface );
		build_one_normal ( D, B, DB );  // based on previous segment
		progress_fill_60 ( CD, DB, BC, D, cloud );
		BC.reverse().add_to ( progress_interface );
		build_one_normal ( B, C, BC );  }  // based on previous segment
	else
	{	DA = Cell ( tag::segment, D.reverse(), A );
		Cell ABD ( tag::triangle, AB, DB.reverse(), DA );
		Cell BCD ( tag::triangle, BC, CD, DB );
		ABD.add_to ( mesh_under_constr );
		BCD.add_to ( mesh_under_constr );
		Cell AD = DA.reverse();
		glue_two_segs_common ( A, B, AB, C, D, CD, AD, BC );  }
	return BC.reverse();                                                 }

//-------------------------------------------------------------------------------------------------

inline void progress_fill_last_triangle
(	const Cell & A, const Cell & B, const Cell & C, Cell & AB, Cell & BC, Cell & CA,
	MetricTree<Cell,Manifold::Euclid::SqDist> & cloud                         )

{	progress_fill_60 ( AB, BC, CA, B, cloud );
	CA.remove_from ( progress_interface );
	assert ( normals.find(CA) != normals.end() );
	normals.erase ( CA );
	cloud.remove ( node_in_cloud [ A ] );
	cloud.remove ( node_in_cloud [ C ] );
	mesh_under_constr.baricenter ( A, CA );
	mesh_under_constr.baricenter ( C, BC );                   }

//-------------------------------------------------------------------------------------------------

inline std::vector < double > compute_tangent_vec
(	const tag::StartAt &, Cell start, const tag::DesiredLength &, double length,
	const tag::OrthogonalTo &, std::vector < double > given_vec                  )

// computes a vector tangent to Manifold::working at point 'start', normal to given_vec

// given_vec must be tangent to Manifold::working at point 'start'
// and must have length approximately equal to desired_length

{	desired_length = length;
	return compute_tangent_vec ( start, true, given_vec );  }
	// 'true' as second argument means "do check orthogonality"
	
//-------------------------------------------------------------------------------------------------

inline std::vector < double > compute_tangent_vec
( const tag::StartAt &, Cell start, const tag::DesiredLength &, double length )

// computes a vector tangent to Manifold::working at point 'start'

{	desired_length = length;
	return compute_tangent_vec ( start, false, std::vector<double>() );  }
	// 'false' as second argument means "do not check orthogonality"

//-------------------------------------------------------------------------------------------------

void progress_relocate
( const Cell & P, size_t n, std::vector<double> & normal_dir,
  const size_t & nc, std::set<Cell> & set_of_ver,
  MetricTree<Cell,Manifold::Euclid::SqDist> & cloud             )

// re-compute the placement of a newly created vertex

// vertex has been located according to two segments, from angles_120 :   n == 2
// or according to only one segment, if built from a brand new triangle : n == 1
	
// we compute here 'set_of_ver' (which is the set of all vertices in the cloud
// close enough to 'ver') and keep it for future use in 'check_touching'

{	// make a list (using a set) of nearby points
	// build a vector of segments from it
	// relocate point P by averaging all normals

	std::list < Cell > list_of_ver =
		cloud.find_close_neighbours_of ( P, progress_long_dist );
	// P has not been added to the cloud yet, so it will not show up in 'list_of_ver'
	set_of_ver.clear();
	for ( std::list<Cell>::const_iterator it = list_of_ver.begin();
        it != list_of_ver.end(); it++                                    )
		set_of_ver.insert ( *it );
	// std::cout << "relocate, n = " << n << ", found " << list_of_ver.size() << " neighbours, ";

	std::vector < Cell > vector_of_seg;
	// vector_of_seg will contain all segments whose both extremities belong to 'set_of-ver'
	// but not the segments adjacent to P (since P does not belong to 'set_of_ver')
	for ( std::set<Cell>::iterator it = set_of_ver.begin(); it != set_of_ver.end(); it++ )
	{	Cell A = *it;
		Cell AB = progress_interface.cell_in_front_of ( A );
		if ( set_of_ver.find ( AB.tip() ) != set_of_ver.end() )
		{	vector_of_seg.push_back ( AB );                        }  }

	// std::cout << vector_of_seg.size() << " segments" << std::endl;
	if ( vector_of_seg.empty() ) return;

	// there are two cases : 
	// the piece of the interface that we just encountered may have normals or not
	size_t counter = 0;
	Cell kept_seg ( tag::non_existent );
	for ( size_t i = 0; i < vector_of_seg.size(); i++ )
	{	Cell seg_p = vector_of_seg[i];
		if ( normals.find ( seg_p ) == normals.end() )
		{	counter++; kept_seg = seg_p;  }                  }
	if ( counter > 0 )  // there are 'counter' segments with no normal
	{	assert ( counter == 1 );
		// build normal of 'kept_seg' from 'normal_dir'
		Cell A = kept_seg.base().reverse();
		Cell B = kept_seg.tip();
		std::vector < double > tangent_dir ( nc );
		for ( size_t i = 0; i < nc; i++ )
		{	Function x = Manifold::working.coordinates()[i];
			tangent_dir[i] = x(B) - x(A);
			normal_dir[i] *= -1.;                             }
	  improve_normal ( A, B, tangent_dir, normal_dir );  // modifies normal_dir
		normals[kept_seg] = normal_dir;
		Cell ret = build_normals ( kept_seg );
		assert ( ret == kept_seg );                                           }

	std::vector < double > pos = Manifold::working.coordinates() ( P );
	if ( n != 1 )
	{	assert ( n == 2 );
		for ( size_t i = 0; i < nc; i++ )  pos[i] *= n;  }
	for ( size_t j = 0; j < vector_of_seg.size(); j++ )
	{	Cell AB = vector_of_seg[j];
	  std::map < Cell, std::vector<double> > :: iterator it =
			normals.find ( AB );
		assert ( it != normals.end() );
		std::vector < double > & nor = it->second;
		// std::vector < double > & nor = normals [ AB ];
		Cell A = AB.base().reverse();
		Cell B = AB.tip();
		for ( size_t i = 0; i < nc; i++ )
		{	Function x = Manifold::working.coordinates()[i];
			pos[i] += ( x(A) + x(B) ) / 2. + nor[i] * sqrt_of_075;  }  }
	n += vector_of_seg.size();
	for ( size_t i = 0; i < nc; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		x(P) = pos[i] / n;                                 }
	
	Manifold::working.project(P);

	// check P is on the right side of kept_seg

}  // end of  progress_relocate

//-------------------------------------------------------------------------------------------------

inline bool check_touching
(	Cell & ver, std::set<Cell> & set_of_ver,
	MetricTree<Cell,Manifold::Euclid::SqDist> & cloud )

// analyse position of recently built vertex 'ver' relatively to other vertices on the interface
// occasionally, different connected components of the interface touch and merge
// or, the current connected component of the interface may touch itself and split in two

// return true if a touch was detected and the corresponding pieces have been merged

// we take advantage of 'set_of_ver' which is the set of all vertices in the cloud
// close enough to 'ver', previously computed in 'relocate'
// we can destroy it here, it won't be used anymore
	
// see paragraph 9.8 in the manual

{	if ( not ver.exists() )  return false;  // no touch
	if ( not ver.belongs_to ( progress_interface, tag::not_oriented ) )  return false;
	//  because 'ver' might have been left behind in the meanwhile

	assert ( set_of_ver.find(ver) == set_of_ver.end() );
	Cell prev_seg = progress_interface.cell_behind(ver);
	Cell prev_ver = prev_seg.base().reverse();
	Cell next_seg = progress_interface.cell_in_front_of(ver);
	Cell next_ver = next_seg.tip();
	set_of_ver.erase ( prev_ver );
	set_of_ver.erase ( next_ver );
	if ( set_of_ver.empty() )  return false;  // no touch
	if ( set_of_ver.size() == 1 )  return false;  // almost touch, no merge
	// do we need to check orientations ? probably not

	// we deal with two vertices or three vertices
	std::set<Cell>::iterator iit = set_of_ver.begin();
	assert ( iit != set_of_ver.end() );
	Cell one = *iit;
	assert ( prev_ver != one );  assert ( next_ver != one );
	iit++;  assert ( iit != set_of_ver.end() );
	Cell two = *iit;
	assert ( prev_ver != two );  assert ( next_ver != two );
	iit++;
	if ( iit != set_of_ver.end() )  // three vertices
	{	Cell three = *iit;
		iit++;  assert ( iit == set_of_ver.end() );
		// we want one-two-three to be in order along the interface
		Cell seg = progress_interface.cell_behind (one);
		Cell tmp = seg.base().reverse();
		if ( set_of_ver.find ( tmp ) != set_of_ver.end() )  // 'one' is not the first one
		{	if ( tmp == two )  { two = one;  one = tmp;  }
			else  {	assert ( tmp == three );  three = one;  one = tmp;  }  }
		seg = progress_interface.cell_behind ( one );
		tmp = seg.base().reverse();
		if ( set_of_ver.find ( tmp ) != set_of_ver.end() )
		// 'one' is still not the first one
		{	if ( tmp == two )
			{ two = one;  one = tmp;  }
			else  {	assert ( tmp == three );  three = one;  one = tmp;  }  }
		// now 'one' must be the first one
		seg = progress_interface.cell_behind ( one );
		assert ( set_of_ver.find ( seg.base().reverse() ) == set_of_ver.end() );
		Cell one_two = progress_interface.cell_in_front_of(one);
		two = one_two.tip();
		Cell two_three = progress_interface.cell_in_front_of(two);
		three = two_three.tip();
		if ( Manifold::working.dist_sq ( one, next_ver ) < progress_long_dist_sq )
		{	// std::cout << "touching interface S fill" << std::endl << std::flush;
			//  progress_interface.cell_in_front_of(next_ver) may have tip 'one'
			Cell ver_two = glue_two_segs_S
				( ver, next_ver, next_seg, one, two, one_two, cloud );
			Cell ver_three ( tag::segment, ver.reverse(), three );
			progress_fill_60 ( ver_two, two_three, ver_three.reverse(), two, cloud );
			ver_three.add_to ( progress_interface );
			build_one_normal ( ver, three, ver_three );  // based on previous segment
			return true;                                                                }
		if ( Manifold::working.dist_sq ( two, prev_ver ) < progress_long_dist_sq )
		{	// std::cout << "touching interface Z fill" << std::endl << std::flush;
			//  progress_interface.cell_behind(prev_ver) may have base 'three'
			Cell two_ver = glue_two_segs_Z
				( prev_ver, ver, prev_seg, two, three, two_three, cloud );
			Cell one_ver ( tag::segment, one.reverse(), ver );
			progress_fill_60 ( one_two, two_ver, one_ver.reverse(), two, cloud );
			one_ver.add_to ( progress_interface );
			build_one_normal ( one, ver, one_ver );  // based on previous segment
			return true;                                                            }    }
	else  // two vertices
	{	// if 'two' is not next to 'one' within the interface, switch them
		Cell seg = progress_interface.cell_in_front_of (one);
		if ( seg.tip() != two )
		{	// in some cases, 'one' and 'two' are not even adjacent to each other
			Cell other_seg = progress_interface.cell_behind ( one );
			if ( other_seg.base().reverse() != two ) return false;  // no merge
			Cell tmp = one;  one = two;
			two = tmp;                            }
		Cell one_two = progress_interface.cell_in_front_of(one);
		assert ( one_two.tip() == two );
		// sometimes vertices are all adjacent
		if ( progress_interface.cell_in_front_of(two).tip() == prev_ver ) return false;  // no merge
		if ( progress_interface.cell_behind(one).base().reverse() == next_ver )
			return false;  // no merge
		if ( Manifold::working.dist_sq ( one, next_ver ) < progress_long_dist_sq )
		{	// std::cout << "touching interface S" << std::endl << std::flush;
			//  progress_interface.cell_in_front_of(next_ver) may have tip 'one'
			glue_two_segs_S ( ver, next_ver, next_seg, one, two, one_two, cloud );
			return true;                                                             }
		if ( Manifold::working.dist_sq ( two, prev_ver ) < progress_long_dist_sq )
		{	// std::cout << "touching interface Z" << std::endl << std::flush;
			//  progress_interface.cell_behind(prev_ver) may have base 'two'
			glue_two_segs_Z ( prev_ver, ver, prev_seg, one, two, one_two, cloud );
			return true;                                                            }  }

	return false;  // almost touch, no merge
	
}  // end of check_touching

//-------------------------------------------------------------------------------------------------

void progressive_construct ( Mesh & msh,
	const tag::StartAt &, const Cell & start, const tag::Towards &, std::vector<double> & tangent,
	const tag::StopAt &, const Cell & stop, const tag::DesiredLength &, double length )

// builds a one-dimensional mesh (a curve)
	
// 'start' and 'stop' are vertices (may be one and the same)
	
{	const size_t nc = Manifold::working.coordinates().nb_of_components();
	assert ( start.dim() == 0 );
	assert ( stop.dim() == 0 );

	desired_length = length;
	double len_sq = desired_length * desired_length;
	double augm_length = desired_length * 1.618034,  // golden number
	       augm_len_sq = augm_length * augm_length;
	size_t counter = 1;
	size_t max_counter = 0;  // std::cin >> max_counter;
	Cell A = start;
	while ( true )
	{	double d = Manifold::working.dist_sq ( A, stop );  // long range distance !
		if ( d < augm_len_sq )
		{	// below we could use the Riemannian metric
			// however, the sign should be the same
			std::vector < double > e ( nc );
			double prod = 0.;
			for ( size_t i = 0; i < nc; i++ )
			{	Function x = Manifold::working.coordinates()[i];
				e[i] = x(stop) - x(A);
				prod += tangent[i] * e[i];                        }
			if ( prod > 0. )
			{	Cell last ( tag::segment, A.reverse(), stop );
				last.add_to ( msh );
				// redistribute vertices
				double n2 = Manifold::working.inner_prod ( A, e, e );
				double norm = approx_sqrt ( n2, tag::around, len_sq, desired_length );
				norm = approx_sqrt ( n2, tag::around, norm*norm, norm );
				redistribute_vertices ( msh, start, stop, norm, 15 );
				break;                                                                 }  }
		Cell B ( tag::vertex );
		for ( size_t i = 0; i < nc; i++ )
		{	Function x = Manifold::working.coordinates()[i];
			x(B) = x(A) + tangent[i];                         }
		if ( counter == max_counter ) return;
		Manifold::working.project ( B );
		for ( size_t i = 0; i < nc; i++ )
		{	Function x = Manifold::working.coordinates()[i];
			tangent[i] = x(B) - x(A);                         }
		double n2 = 0.;  // Manifold::working.inner_prod ( A, tangent, tangent ); ?!!
		for ( size_t i = 0; i < nc; i++ )
		{	double temp = tangent[i];  temp *= temp;  n2 += temp;  }
		n2 = approx_sqrt ( n2, tag::around, len_sq, desired_length );
		n2 = desired_length / n2;
		for ( size_t i = 0; i < nc; i++ )  tangent[i] *= n2;
		Cell AB ( tag::segment, A.reverse(), B );
		AB.add_to ( msh );  counter++;  A = B;                                            }

} // end of  progressive_construct

//-------------------------------------------------------------------------------------------------

void progressive_construct ( Mesh & msh, const tag::StartAt &, const Cell & start,
	const tag::StopAt &, const Cell & stop, const tag::DesiredLength &, double length )
	
// builds a one-dimensional mesh (a curve)
	
// no tangent vector provided, thus no way of knowing which way to go
// the only solution is to start walking in both directions simultaneously
// the first one to reach 'stop' wins

// 'start' and 'stop' are vertices (may be one and the same)

{	const size_t nc = Manifold::working.coordinates().nb_of_components();
	assert ( start.dim() == 0 );
	assert ( stop.dim() == 0 );

	std::vector < double > best_tangent = compute_tangent_vec
		( tag::start_at, start, tag::desired_length, length );

	// start walking along the manifold from 'start' in the direction of best_tangent
	// and, simultaneously, in the opposite direction, given by -best_tangent
	double augm_length = length * 1.5,
	       augm_len_sq = augm_length * augm_length;
	std::vector < double > tan1 = best_tangent, tan2 = best_tangent;
	for ( size_t i = 0; i < nc; i++ ) tan2[i] *= -1.;
	Cell ver1 ( tag::vertex ), ver2 ( tag::vertex );
	for ( size_t i = 0; i < nc; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		x ( ver1 ) = x ( start );  x ( ver2 ) = x ( start );  }
	int winner;  //  will be 1 or -1
	while ( true )
	{	for ( size_t i = 0; i < nc; i++ )
		{	Function x = Manifold::working.coordinates()[i];
			x ( temporary_vertex ) = x ( ver1 ) + tan1[i];  }
		Manifold::working.project ( temporary_vertex );
		for ( size_t i = 0; i < nc; i++ )
		{	Function x = Manifold::working.coordinates()[i];
			tan1[i] = x ( temporary_vertex ) - x ( ver1 );
			x ( ver1 ) = x ( temporary_vertex );            }
		double d = Manifold::working.dist_sq ( ver1, stop );
		if ( d < augm_len_sq )
		{	double prod = 0.;
			for ( size_t i = 0; i < nc; i++ )
			{	Function x = Manifold::working.coordinates()[i];
				prod += tan1[i] * ( x(stop) - x(ver1) );          }
			if ( prod > 0. )  { winner = 1;  break;  }           }
		for ( size_t i = 0; i < nc; i++ )
		{	Function x = Manifold::working.coordinates()[i];
			x ( temporary_vertex ) = x ( ver2 ) + tan2[i];  }
		Manifold::working.project ( temporary_vertex );
		for ( size_t i = 0; i < nc; i++ )
		{	Function x = Manifold::working.coordinates()[i];
			tan2[i] = x ( temporary_vertex ) - x ( ver2 );
			x ( ver2 ) = x ( temporary_vertex );            }
		d = Manifold::working.dist_sq ( ver2, stop );
		if ( d < augm_len_sq )
		{	double prod = 0.;
			for ( size_t i = 0; i < nc; i++ )
			{	Function x = Manifold::working.coordinates()[i];
				prod += tan2[i] * ( x(stop) - x(ver2) );          }
			if ( prod > 0. )  { winner = -1;  break;  }           }
	}  // end of  while true

	ver1.dispose();  ver2.dispose();

	for ( size_t i = 0; i < nc; i++ ) best_tangent[i] *= winner;
	progressive_construct ( msh, tag::start_at, start, tag::towards, best_tangent,
	                        tag::stop_at, stop, tag::desired_length, length   );
}

//-------------------------------------------------------------------------------------------------

void progressive_construct
( Mesh & msh, const tag::StartAt &, const Cell & start, const tag::StopAt &, const Cell & stop,
  const tag::DesiredLength &, double length, const tag::InherentOrientation & )

// 'start' and 'stop' are vertices (may be one and the same)

{	const size_t nc = Manifold::working.coordinates().nb_of_components();
	assert ( start.dim() == 0 );
	assert ( stop.dim() == 0 );

	std::vector < double > best_tangent = compute_tangent_vec
		( tag::start_at, start, tag::desired_length, length );

	Mesh msh1 ( tag::of_dimension_one );
	progressive_construct ( msh1, tag::start_at, start, tag::towards, best_tangent,
                          tag::stop_at, stop, tag::desired_length, length );
	for ( size_t i = 0; i < nc; i++ )  best_tangent[i] *= -1.;
	Mesh msh2 ( tag::of_dimension_one );
	progressive_construct ( msh2, tag::start_at, start, tag::towards, best_tangent,
                          tag::stop_at, stop, tag::desired_length, length );
	switch_orientation ( msh2 );
	Mesh whole ( tag::join, msh1, msh2 );

	if ( correctly_oriented ( whole ) )
	{	msh2.dispose();  msh = msh1;  }
	else
	{	msh1.dispose();  switch_orientation ( msh2 );  msh = msh2;  }

	whole.dispose();
}
	
//-------------------------------------------------------------------------------------------------

void progressive_construct
( Mesh & msh, const tag::StartAt &, Cell start,
  const tag::Towards &, std::vector<double> & normal,
  const tag::Boundary &, Mesh bdry,
  const tag::DesiredLength, double length        )

// for two-dimensional meshes (arbitrary geometric dimension)
	
// 'start' is a vertex or segment belonging to 'bdry'
// 'normal' is a vector tangent to the working manifold, orthogonal to 'start'

{	desired_length = length;
	desired_length_sq = desired_length * desired_length;
	progress_long_dist = desired_length * one_plus_tolerance;
	progress_long_dist_sq = progress_long_dist * progress_long_dist;

	// we don't want to change 'bdry' so we make a copy of it
	// Mesh interface ( tag::deep_copy_of, bdry );  // wait, we want to switch to a list of interfaces ...
	progress_interface = bdry;
	mesh_under_constr = msh;
	const size_t nc = Manifold::working.coordinates().nb_of_components();
	progress_nb_of_coords = nc;
	Cell vertex_recently_built ( tag::non_existent );
	std::set < Cell > set_of_nearby_vertices;
	// vertices close to vertex_recently_built

	Manifold::Euclid::SqDist square_dist;
	MetricTree<Cell,Manifold::Euclid::SqDist> cloud ( square_dist, desired_length, 6. );
	// first argument : a callable object returning the square of the distance
	// between two points, measured in the surrounding, Euclidian, space
	// second argument : distance for rank zero nodes, which is a mere hint
	// about how to initialize the tree; the tree can change a lot later
	// third argument : ratio between distances of successive ranks

	{ // just a block of code for hiding variables
	// 'progress_interface' is a one-dimensional mesh, not necessarily connected
	// so we cannot use a CellIterator - perhaps an unstructured one ?
	std::list < Cell::Core* > & l = progress_interface.core->cells[0];
	std::list<Cell::Core*>::iterator it = l.begin();
	for ( ; it != l.end(); it++ )
	{	Cell::Core * ver_p = *it;
		Cell ver ( tag::whose_core_is, ver_p );
		progress_add_point ( ver, cloud );       }
	} // just a block of code for hiding variables

	if ( start.dim() != 1 )
	{	assert ( start.dim() == 0 );
		start = progress_interface.cell_in_front_of ( start );                    }
	assert ( start.dim() == 1 );
	assert ( bdry.dim() == 1 );
	assert ( start.belongs_to ( progress_interface, tag::oriented ) );
	
	{ // just a block of code for hiding variables
	double n2 = 0.;  // ensure that normal has the right norm
	for ( size_t i = 0; i < nc; i++ )
	{	double temp = normal[i];  temp *= temp;  n2 += temp;  }
	double coef = desired_length / std::sqrt(n2);
	for ( size_t i = 0; i < nc; i++ )  normal[i] *= coef;
	Cell start_base = start.base().reverse();
	Cell start_tip = start.tip();
	std::vector < double > vec ( nc );
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		vec[i] = x ( start_tip ) - x ( start_base );      }
	improve_normal ( start_base, start_tip, vec, normal );
	// ensures again the norm is right, projects on the tangent space,
	// ensures orthogonality with 'start' segment
	} // just a block of code for hiding variables

	normals[start] = normal;
	Cell ret = build_normals ( start );
	assert ( ret == start );

	Cell point_60 = start.tip();
	
	int stopping_criterion = 0;
	// std::cout << "stopping criterion : ";  std::cin >> stopping_criterion;
	int current_name = 1;

restart:
	
	Cell stop_point_60 = point_60, point_120 = point_60, stop_point_120 = point_60;

angles_60 :

	// starting at 'point_60', we go along this connected component of 'progress_interface'
	{ // just a block of code for hiding variables
	Cell prev_seg = progress_interface.cell_behind ( point_60, tag::surely_exists );
	Cell A = prev_seg.base().reverse();
	while ( true )
	{	Cell next_seg = progress_interface.cell_in_front_of ( point_60, tag::surely_exists );
		Cell B = next_seg.tip();  assert ( A != B );
		double d = Manifold::working.dist_sq ( A, B );
		if ( d < 1.1 * progress_long_dist_sq )  // we may want to form a triangle
		// but first we must make sure it is on the correct side
		if ( positive_orientation ( A, point_60, prev_seg, next_seg ) )
		if ( ( progress_cos_sq_60 ( A, point_60, B, prev_seg, next_seg) > 0.02 )
		  or ( d < desired_length )                                              )
		// triangle waiting to be filled; see paragraph 9.7 in the manual
		{	Cell seg_next_to_B = progress_interface.cell_in_front_of(B);
			Cell ver_next_to_B = seg_next_to_B.tip();
			set_of_nearby_vertices.erase ( point_60 );
			if ( ver_next_to_B == A )  // this is the last triangle in this piece of progress_interface
			{	progress_fill_last_triangle
					( A, point_60, B, prev_seg, next_seg, seg_next_to_B, cloud );
				std::cout << "shrinking triangle " << ++current_name << std::endl;
				if ( current_name == stopping_criterion ) return;			
				goto search_for_start;  	                                                        }
			Cell AB ( tag::segment, A.reverse(), B );
			progress_fill_60 ( prev_seg, next_seg, AB.reverse(), point_60, cloud );
			AB.add_to ( progress_interface );
			build_one_normal ( A, B, AB );  // based on previous segment
			std::cout << "found angle around 60 deg " << ++current_name << std::endl;
			if ( current_name == stopping_criterion ) return;
			if ( stop_point_120 == B )
			{	if ( stop_point_120 == point_120 )
				{	point_120 = A;  stop_point_120 = A;  }
				else  stop_point_120 = ver_next_to_B;    }
			if ( stop_point_120 == point_60 )  // we have all the loop to cover for 120 deg
			{	assert ( point_120 == stop_point_120 );
				stop_point_120 = A;                     }
			if ( point_120 == point_60 ) point_120 = A;
			if ( stop_point_60 == B )  stop_point_60 = ver_next_to_B;
			if ( stop_point_60 == point_60 )  // we have all the loop to cover for 60 deg
				stop_point_60 = A;
			point_60 = A;  goto angles_60;                                                    }
		if ( B == stop_point_60 ) break;
		A = point_60;  prev_seg = next_seg;  point_60 = B;                                            }
	} // just a block of code for hiding variables

check_touching :

	{ // just a block of code for hiding 'touch'
	bool touch = check_touching ( vertex_recently_built, set_of_nearby_vertices, cloud );
	vertex_recently_built = Cell ( tag::non_existent );
	if ( touch )
	{	assert ( point_120.belongs_to ( progress_interface, tag::not_oriented ) );
		std::cout << "touch " << ++current_name << std::endl;
		if ( current_name == stopping_criterion ) return;			
		point_60 = point_120;  stop_point_60 = point_60;  stop_point_120 = point_60;
		goto angles_60;                                                                     }
	} // just a block of code for hiding 'touch'

// look for angles around 120 deg :

	// starting at 'point_120', we go along this connected component of 'progress_interface'
	{  // just a block of code for hiding prev_seg and A
	Cell prev_seg = progress_interface.cell_behind ( point_120, tag::surely_exists );
	Cell A = prev_seg.base().reverse();
	while ( true )
	{	Cell next_seg = progress_interface.cell_in_front_of ( point_120, tag::surely_exists );
		Cell  B = next_seg.tip();
		if ( progress_cos_sq_120 ( A, point_120, B, prev_seg, next_seg) < 0.55 )  // 0.67
		// angle around 120 deg, we want to form two triangles; see paragraph 9.7 in the manual
		{	// we don't build a new vertex yet, we want to check for a quadrangle first
			Cell seg_prev_to_A = progress_interface.cell_behind ( A );
			Cell seg_next_to_B = progress_interface.cell_in_front_of ( B );
			Cell ver_prev_to_A = seg_prev_to_A.base().reverse();
			Cell ver_next_to_B = seg_next_to_B.tip();
			if ( ver_prev_to_A == ver_next_to_B )  // quadrangle
			{	// check orientations : correct side ?
				std::cout << "shrinking quadrangle" << std::endl;
				// choose the shortest diagonal, add two triangles to the mesh under construction
				if ( Manifold::working.dist_sq ( A, B ) <
						 Manifold::working.dist_sq ( point_120, ver_next_to_B ) )
				{	Cell AB ( tag::segment, A.reverse(), B );
					progress_fill_60 ( prev_seg, next_seg, AB.reverse(), point_120, cloud );
					progress_fill_60 ( seg_next_to_B, seg_prev_to_A, AB, ver_prev_to_A, cloud ); }
				else
				{	Cell seg ( tag::segment, point_120.reverse(), ver_next_to_B );
					progress_fill_60 ( next_seg, seg_next_to_B, seg.reverse(), B, cloud );
					progress_fill_60 ( seg_prev_to_A, prev_seg, seg, A, cloud );            }
				goto search_for_start;                                                             }
			Cell P ( tag::vertex );  vertex_recently_built = P;
			// now we want to place this new vertex accordingly
			assert ( normals.find(prev_seg) != normals.end() );
			assert ( normals.find(next_seg) != normals.end() );
			std::vector < double > & nor_a = normals [ prev_seg ];
			std::vector < double > & nor_b = normals [ next_seg ];
			std::vector < double > sum_of_nor ( nc );
			for ( size_t i = 0; i < nc; i++ )  sum_of_nor[i] = nor_a[i] + nor_b[i];
			for ( size_t i = 0; i < nc; i++ )
			{	Function x = Manifold::working.coordinates()[i];
				x(P) = x(point_120) / 2. + ( x(A) + x(B) ) / 4.
					+ sum_of_nor[i] * half_of_sqrt_of_075;         }
			Manifold::working.project(P);
			Cell AP ( tag::segment, A.reverse(), P );
			Cell BP ( tag::segment, B.reverse(), P );
			Cell PB = BP.reverse();
			Cell sP ( tag::segment, point_120.reverse(), P );
			Cell tri1 ( tag::triangle, prev_seg, sP, AP.reverse() );
			Cell tri2 ( tag::triangle, next_seg, BP, sP.reverse() );
			tri1.add_to ( msh );
			tri2.add_to ( msh );
			prev_seg.remove_from ( progress_interface );
			next_seg.remove_from ( progress_interface );
			normals.erase ( prev_seg );
			normals.erase ( next_seg );
			AP.add_to ( progress_interface );
			PB.add_to ( progress_interface );
			cloud.remove ( node_in_cloud[point_120] );
			build_one_normal ( A, P, AP );  // based on previous segment
			build_one_normal ( P, B, PB );  // based on previous segment
			progress_relocate ( P, 2, sum_of_nor, nc, set_of_nearby_vertices, cloud );
			// find more vertices close to P and take them all into account; modifies sum_of_nor
			assert ( prev_seg.tip() == point_120 );
			msh.baricenter ( point_120, prev_seg );
			progress_add_point ( P, cloud );
			std::cout << "found angle around 120 deg " << ++current_name << std::endl;
			if ( current_name == stopping_criterion ) return;			
			if ( stop_point_120 == point_120 )  // we have all the loop to cover
				stop_point_120 = A;
			if ( stop_point_120 == B )  stop_point_120 = ver_next_to_B;
			point_120 = A;
			point_60 = A;
			stop_point_60 = ver_next_to_B;
			goto angles_60;                                                        }
		if ( B == stop_point_120 ) break;
		A = point_120;  prev_seg = next_seg;  point_120 = B;
	}  // end of while
	}  // just a block of code for hiding prev_seg and A

// build a brand new triangle :
	
	{ // just a block of code for hiding variables
	// we use point_120, any other point would do
	Cell next_seg = progress_interface.cell_in_front_of ( point_120, tag::surely_exists );
	Cell B = next_seg.tip();
	Cell P ( tag::vertex );  vertex_recently_built = P;
	assert ( normals.find(next_seg) != normals.end() );
	std::vector < double > f = normals [ next_seg ];
	for ( size_t i = 0; i < nc; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		x(P) = ( x(point_120) + x(B) ) / 2. + f[i] * sqrt_of_075;  }
	Manifold::working.project(P);
	Cell AP ( tag::segment, point_120.reverse(), P );
	Cell BP ( tag::segment, B.reverse(), P );
	Cell PB = BP.reverse();
	Cell tri ( tag::triangle, next_seg, BP, AP.reverse() );
	tri.add_to ( msh );
	next_seg.remove_from ( progress_interface );
	normals.erase ( next_seg );
	AP.add_to ( progress_interface );
	PB.add_to ( progress_interface );
	std::cout << "building brand new triangle " << ++current_name << std::endl;
	if ( current_name == stopping_criterion ) return;			
	build_one_normal ( point_120, P, AP );  // based on previous segment
	build_one_normal ( P, B, PB );  // based on previous segment
	progress_relocate ( P, 1, f, nc, set_of_nearby_vertices, cloud );
	// find more vertices close to P and take them all into account; modifies f
	progress_add_point ( P, cloud );
	stop_point_120 = progress_interface.cell_in_front_of(B).tip();
	goto check_touching;
	} // just a block of code for hiding variables

search_for_start :
// execution only reaches this point through 'goto'

	{ // just a block of code for hiding variables
	// we look for a segment in 'progress_interface' which has a normal
	// 'progress_interface' is a one-dimensional mesh, not necessarily connected
	// so we cannot use a CellIterator - perhaps an unstructured one ?
	std::list < Cell::Core* > & l = progress_interface.core->cells[1];
	if ( l.size() == 0 ) return;  // empty interface, meshing process ended
	std::cout << "search for start " << current_name << std::endl;
	Cell start_seg ( tag::non_existent );
	assert ( l.size() >= 2 );
	std::list<Cell::Core*>::iterator it = l.begin();
	for ( ; it != l.end(); it++ )
	{	Cell tmp ( tag::whose_core_is, *it );
		if ( normals.find ( tmp ) != normals.end() )
		{	start_seg = tmp;  break;  }                }
	assert ( start_seg.exists() );
	point_60 = start_seg.tip();
	// any point on this connected component would do
	goto restart;
	} // just a block of code for hiding variables

}  // end of progressive_construct

//-------------------------------------------------------------------------------------------------

void progressive_construct
( Mesh & msh, const tag::DesiredLength &, double length, bool check_and_switch )
	
// last argument tells whether to check the orientation of the resulting mesh
// and switch it if necessary

{	const size_t nc = Manifold::working.coordinates().nb_of_components();
	
	// call to 'search_start_ver' does not depend on the dimension of the mesh
	Cell start = search_start_ver ( length );

	// call to 'compute_tangent_vec' does not depend on the dimension of the mesh
	std::vector < double > tangent = compute_tangent_vec
		( tag::start_at, start, tag::desired_length, length );

	// now we branch, depending on the dimension
	if ( msh.dim() == 1 )
		progressive_construct ( msh, tag::start_at, start, tag::towards, tangent,
		                        tag::stop_at, start, tag::desired_length, length );
	else
	{	assert ( msh.dim() == 2 );  // no 3D for now
		assert ( nc == 3 );
		Cell B ( tag::vertex );
		for ( size_t i = 0; i < nc; i++ )
		{	Function x = Manifold::working.coordinates()[i];
			x ( B ) = x ( start ) + tangent[i];              }
		Cell AB ( tag::segment, start.reverse(), B );
		std::vector < double > normal = compute_tangent_vec ( tag::start_at, start,
		  tag::desired_length, length, tag::orthogonal_to, tangent );
		Cell C ( tag::vertex );
		for ( size_t i = 0; i < nc; i++ )
		{	Function x = Manifold::working.coordinates()[i];
			x ( C ) = x ( start ) + 0.5 * tangent[i] - 0.866 * normal[i];  }
		Cell BC ( tag::segment, B.reverse(), C );
		Cell CA ( tag::segment, C.reverse(), start );
		Cell tri ( tag::triangle, AB, BC, CA );
		tri.add_to ( msh );
		Mesh interf ( tag::of_dimension_one );
		AB.reverse().add_to(interf);
		BC.reverse().add_to(interf);
		CA.reverse().add_to(interf);
		progressive_construct ( msh, tag::start_at, AB.reverse(), tag::towards, normal,
		                        tag::boundary, interf, tag::desired_length, length );
	}  // end of  else  with  msh.dim() == 2

	if ( not check_and_switch ) return;
	if ( not correctly_oriented ( msh ) )  switch_orientation ( msh );
	assert ( correctly_oriented ( msh ) );                                                }

//-------------------------------------------------------------------------------------------------
	
void progressive_construct
( Mesh & msh, const tag::StartAt &, const Cell & start,
  const tag::Boundary &, Mesh interface,
  const tag::DesiredLength, double length,
	const tag::InherentOrientation, bool check_and_switch )

{	assert ( Manifold::working.coordinates().nb_of_components() == 3 );
	assert ( msh.dim() == 2 );
	assert ( start.core->belongs_to ( interface.core, tag::oriented ) );
	Mesh interf_rev ( tag::of_dimension_one );
	CellIterator it = interface.iter_over ( tag::cells_of_dim, 1 );
	for ( it.reset(); it.in_range(); it++ )
		(*it).reverse().add_to ( interf_rev );
	assert ( start.reverse().core->belongs_to ( interf_rev.core, tag::oriented ) );

	// compute a normal vector, on an arbitrary side of 'start'
	Cell A = start.base().reverse();
	Cell B = start.tip();
	Function x = Manifold::working.coordinates()[0];
	Function y = Manifold::working.coordinates()[1];
	Function z = Manifold::working.coordinates()[2];
	std::vector < double > tan { x(B) - x(A), y(B) - y(A), z(B) - z(A) };
	std::vector < double > nor = compute_tangent_vec
		( tag::start_at, A, tag::desired_length, length, tag::orthogonal_to, tan );

	progressive_construct ( msh, tag::start_at, start, tag::towards, nor,
	                        tag::boundary, interface, tag::desired_length, length );	

	if ( not check_and_switch ) return;

	// build the mesh on the other side of 'interface'
	for ( size_t i = 0; i < 3; i++ )  nor[i] *= -1.;
	Mesh msh2 ( tag::of_dimension, 2, tag::greater_than_one );
	progressive_construct ( msh2, tag::start_at, start.reverse(), tag::towards, nor,
	                        tag::boundary, interf_rev, tag::desired_length, length );	
	// join everything to get a mesh on the entire manifold
	std::cout << "wait a minute ..." << std::endl;
	// we should switch_orientation ( msh2 )
	Mesh glob ( tag::join, msh, msh2 );

	if ( not correctly_oriented ( glob ) )
	{	switch_orientation ( msh2 );  msh = msh2;  msh.dispose();  }
	else  msh2.dispose();
	glob.dispose();
}

//-------------------------------------------------------------------------------------------------
	
void progressive_construct
( Mesh & msh, const tag::StartAt &, const Cell & start,
  const tag::Boundary &, Mesh interface,
  const tag::DesiredLength, double length        )

// for two-dimensional meshes in RR^2 (intrinsic orientation)
//   or in a 2D submanifold of RR^3 (random orientation)
	
// 'start' is a segment belonging to 'interface'
// no normal vector provided, we need to build our own

{	if ( Manifold::working.coordinates().nb_of_components() == 3 )
	// surface in RR^3
	{	progressive_construct ( msh, tag::start_at, start, tag::boundary, interface,
		        tag::desired_length, length, tag::inherent_orientation, false );
		// last argument 'false' means do not check orientation, leave it random
		return;                                                                            }
	
	// domain in the plane RR^2
	Cell A = start.base().reverse();
	Cell B = start.tip();
	Function x = Manifold::working.coordinates()[0];
	Function y = Manifold::working.coordinates()[1];
	std::vector < double > nor { y(A) - y(B), x(B) - x(A) };   // rotate with 90 deg

	progressive_construct ( msh, tag::start_at, start, tag::towards, nor,
	                        tag::boundary, interface, tag::desired_length, length );	
}
	
//-------------------------------------------------------------------------------------------------

void progressive_construct
( Mesh & msh, const tag::Boundary &, Mesh interface,
  const tag::DesiredLength, double length    )

// for two-dimensional meshes in RR^2 (intrinsic orientation)
//   or in a 2D submanifold of RR^3 (random orientation)
	
{	// we search for a starting point
	// 'interface' is a one-dimensional mesh, not necessarily connected
	// so we cannot use a CellIterator - perhaps an unstructured one ?
	std::list < Cell::Core* > & l = interface.core->cells[1];
	std::list<Cell::Core*>::iterator it = l.begin();
	Cell start ( tag::whose_core_is, *it );

	progressive_construct ( msh, tag::start_at, start, tag::boundary, interface,
	                        tag::desired_length, length                  );       }

//-------------------------------------------------------------------------------------------------

void progressive_construct
( Mesh & msh, const tag::Boundary &, Mesh interface,
  const tag::DesiredLength, double length,
	const tag::InherentOrientation, bool check_and_switch )

// if last argument is true, compute inherent orientation
// otherwise, random orientation

{	// we search for a starting point
	// 'interface' is a one-dimensional mesh, not necessarily connected
	// so we cannot use a CellIterator - perhaps an unstructured one ?
	std::list < Cell::Core* > & l = interface.core->cells[1];
	std::list<Cell::Core*>::iterator it = l.begin();
	Cell start ( tag::whose_core_is, *it );

	progressive_construct ( msh, tag::start_at, start, tag::boundary, interface,
	                        tag::desired_length, length,
                          tag::inherent_orientation, check_and_switch );               }

//-------------------------------------------------------------------------------------------------

inline size_t get_topological_dim ( )

{	Manifold::Implicit::OneEquation * m_impl_1 =
		dynamic_cast<Manifold::Implicit::OneEquation*> ( Manifold::working.core );
	if ( m_impl_1 )
	{	Manifold::Euclid * m_euclid = dynamic_cast<Manifold::Euclid*>
			( m_impl_1->surrounding_space.core );
		assert ( m_euclid );
		return m_euclid->coord_func.nb_of_components() - 1;              }
	Manifold::Implicit::TwoEquations * m_impl_2 =
		dynamic_cast<Manifold::Implicit::TwoEquations*> ( Manifold::working.core );
	if ( m_impl_2 )
	{	Manifold::Euclid * m_euclid = dynamic_cast<Manifold::Euclid*>
			( m_impl_2->surrounding_space.core );
		assert ( m_euclid );
		return m_euclid->coord_func.nb_of_components() - 2;              }
	assert ( false );                                                               } 

//-------------------------------------------------------------------------------------------------

}  // end of anonymous namespace

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::DesiredLength &, double length )

:	Mesh ( tag::of_dimension, get_topological_dim(), tag::might_be_one )

{	temporary_vertex = Cell ( tag::vertex );

	progressive_construct ( *this,  tag::desired_length, length, true );
	// last argument true means : check orientation, switch it if necessary

	temporary_vertex.dispose();                                            }

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::EntireManifold, Manifold manif,
             const tag::DesiredLength &, double length                            )

:	Mesh ( tag::of_dimension, get_topological_dim(), tag::might_be_one )

{	Manifold tmp_manif = Manifold::working;
	Manifold::working = manif;
	
	temporary_vertex = Cell ( tag::vertex );

	progressive_construct ( *this,  tag::desired_length, length, true );
	// last argument true means : check orientation, switch it if necessary

	temporary_vertex.dispose();
	Manifold::working = tmp_manif;                                           }

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::DesiredLength &, double length,
             const tag::RandomOrientation &                                        )

:	Mesh ( tag::of_dimension, get_topological_dim(), tag::might_be_one )

{	temporary_vertex = Cell ( tag::vertex );

	progressive_construct ( *this,  tag::desired_length, length, false );
	// last argument false means : do not check orientation, leave it random

	temporary_vertex.dispose();                                             }

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::EntireManifold, Manifold manif,
             const tag::DesiredLength &, double length, const tag::RandomOrientation & )

:	Mesh ( tag::of_dimension, get_topological_dim(), tag::might_be_one )

{	Manifold tmp_manif = Manifold::working;
	Manifold::working = manif;
	
	temporary_vertex = Cell ( tag::vertex );

	progressive_construct ( *this,  tag::desired_length, length, false );
	// last argument false means : do not check orientation, leave it random

	temporary_vertex.dispose();
	Manifold::working = tmp_manif;                                          }

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::DesiredLength &, double length,
             const tag::InherentOrientation &                                     )

:	Mesh ( tag::of_dimension, get_topological_dim(), tag::might_be_one )

{	assert ( Manifold::working.coordinates().nb_of_components() == this->dim() + 1 );
	temporary_vertex = Cell ( tag::vertex );

	progressive_construct ( *this,  tag::desired_length, length, true );
	// last argument true means : check orientation, switch it if necessary

	temporary_vertex.dispose();                                                        }

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::EntireManifold, Manifold manif,
             const tag::DesiredLength &, double length, const tag::InherentOrientation & )

:	Mesh ( tag::of_dimension, get_topological_dim(), tag::might_be_one )

{	assert ( Manifold::working.coordinates().nb_of_components() == this->dim() + 1 );
	Manifold tmp_manif = Manifold::working;
	Manifold::working = manif;
	temporary_vertex = Cell ( tag::vertex );

	progressive_construct ( *this,  tag::desired_length, length, true );
	// last argument true means : check orientation, switch it if necessary

	temporary_vertex.dispose();
	Manifold::working = tmp_manif;                                                        }

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
             const tag::Towards &, std::vector<double> & tangent,
             const tag::StopAt &, const Cell & stop,
             const tag::DesiredLength &, double length                            )

// 'start' and 'stop' may be the same cell

:	Mesh ( tag::of_dimension_one )  // positive, by default

{	temporary_vertex = Cell ( tag::vertex );

	const size_t nc = Manifold::working.coordinates().nb_of_components();	
	double n2 = 0.;  // ensure that normal has the right norm
	for ( size_t i = 0; i < nc; i++ )
	{	double temp = tangent[i];  temp *= temp;  n2 += temp;  }
	double coef = length / std::sqrt(n2);
	for ( size_t i = 0; i < nc; i++ )  tangent[i] *= coef;
	desired_length = length;

	improve_normal ( start, tangent );
	// ensures again the norm is right, then projects on the tangent space

	progressive_construct ( *this, tag::start_at, start, tag::towards, tangent,
                          tag::stop_at, stop, tag::desired_length, length );

	temporary_vertex.dispose();                                                  }

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
             const tag::Towards &, std::vector<double> & tangent,
             const tag::DesiredLength &, double length                           )

// 'start' and 'stop' may be the same cell

:	Mesh ( tag::of_dimension_one )  // positive, by default

{	temporary_vertex = Cell ( tag::vertex );

	const size_t nc = Manifold::working.coordinates().nb_of_components();	
	double n2 = 0.;  // ensure that normal has the right norm
	for ( size_t i = 0; i < nc; i++ )
	{	double temp = tangent[i];  temp *= temp;  n2 += temp;  }
	double coef = length / std::sqrt(n2);
	for ( size_t i = 0; i < nc; i++ )  tangent[i] *= coef;
	desired_length = length;

	improve_normal ( start, tangent );
	// ensures again the norm is right, then projects on the tangent space

	progressive_construct ( *this, tag::start_at, start, tag::towards, tangent,
                          tag::stop_at, start, tag::desired_length, length );

	temporary_vertex.dispose();                                                  }

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
             const tag::StopAt &, const Cell & stop,
             const tag::DesiredLength &, double length                            )

:	Mesh ( tag::of_dimension_one )  // positive, by default

{	temporary_vertex = Cell ( tag::vertex );
	
	progressive_construct ( *this, tag::start_at, start, tag::stop_at, stop,
	                        tag::desired_length, length              );

	temporary_vertex.dispose();                                              }

//-------------------------------------------------------------------------------------------------

// below we should use the inherent orientation if the co-dimension is 1
// also, we should define a similar constructor ending in tag::random_orientation
// a similar constructor ending in tag::inherent_orientation is defined inline in mesh.h

Mesh::Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
             const tag::DesiredLength &, double length                           )

// dimension may be 2 !
:	Mesh ( tag::of_dimension_one )  // positive, by default

{	temporary_vertex = Cell ( tag::vertex );
	
	progressive_construct ( *this, tag::start_at, start, tag::stop_at, start,
	                        tag::desired_length, length                       );

	temporary_vertex.dispose();                                                  }

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
             const tag::StopAt &, const Cell & stop,
             const tag::DesiredLength &, double length, const tag::InherentOrientation & )

:	Mesh ( tag::of_dimension_one )  // positive, by default

{	temporary_vertex = Cell ( tag::vertex );
	
	progressive_construct ( *this, tag::start_at, start, tag::stop_at, stop,
	                        tag::desired_length, length, tag::inherent_orientation );

	temporary_vertex.dispose();                                                       }

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::Boundary &, Mesh interface,
             const tag::DesiredLength &, double length                        )

// for now, only works for two-dimensional meshes (either in RR2 or in RR3)
// should be adapted for three-dimensional meshes

:	Mesh ( tag::of_dimension, 2, tag::greater_than_one )  // positive, by default

{	temporary_vertex = Cell ( tag::vertex );

	progressive_construct ( *this, tag::boundary, interface, tag::desired_length, length );
	
	temporary_vertex.dispose();                                                              }

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::Boundary &, Mesh interface,
             const tag::DesiredLength &, double length, const tag::IntrinsicOrientation & )

// for now, only works for two-dimensional meshes in RR2
// should be adapted for three-dimensional meshes

:	Mesh ( tag::of_dimension, 2, tag::greater_than_one )  // positive, by default

{	temporary_vertex = Cell ( tag::vertex );

	assert ( Manifold::working.coordinates().nb_of_components() == 2 );
	progressive_construct ( *this, tag::boundary, interface, tag::desired_length, length );
	
	temporary_vertex.dispose();                                                              }

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::Boundary &, Mesh interface,
             const tag::DesiredLength &, double length, const tag::InherentOrientation & )

// for two-dimensional meshes in RR^3

:	Mesh ( tag::of_dimension, 2, tag::greater_than_one )  // positive, by default

{	temporary_vertex = Cell ( tag::vertex );

	assert ( Manifold::working.coordinates().nb_of_components() == 3 );
	progressive_construct ( *this, tag::boundary, interface,
                          tag::desired_length, length, tag::inherent_orientation, true );
	
	temporary_vertex.dispose();                                                        }

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::Boundary &, Mesh interface,
             const tag::StartAt &, const Cell & start,
             const tag::Towards &, std::vector<double> & normal,             
             const tag::DesiredLength &, double length                        )
	
// 'start' is a vertex or segment belonging to 'interface'

:	Mesh ( tag::of_dimension, 2, tag::greater_than_one )  // positive, by default

{	temporary_vertex = Cell ( tag::vertex );

	progressive_construct ( *this, tag::start_at, start, tag::towards, normal,
	                        tag::boundary, interface, tag::desired_length, length );
	
	temporary_vertex.dispose();                                                       }

//-------------------------------------------------------------------------------------------------

int main5 ( )

// skeleton of a tetrahedron
	
{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz[0],  y = xyz[1],  z = xyz[2];

	// A ( 1, 0, 1 )
	// B ( 0, 1,-1 )
	// C (-1, 0, 1 )
	// D ( 0,-1,-1 )
	Function lambda_AB = ( x - y + 2.*z - 3. ) / 6.;
	Function dist_AB_sq = ( x - 1. - lambda_AB ) * ( x - 1. - lambda_AB )
		+ ( y + lambda_AB ) * ( y + lambda_AB )
		+ ( z - 1. - 2.*lambda_AB ) * ( z - 1. - 2.*lambda_AB );
	Function dist_AC_sq = y*y + (z-1.)*(z-1.);
	Function lambda_AD = ( x + y + 2.*z - 3. ) / 6.;
	Function dist_AD_sq = ( x - 1. - lambda_AD ) * ( x - 1. - lambda_AD )
		+ ( y - lambda_AD ) * ( y - lambda_AD )
		+ ( z - 1. - 2.*lambda_AD ) * ( z - 1. - 2.*lambda_AD );
	Function lambda_BC = ( - x - y + 2.*z + 3. ) / 6.;
	Function dist_BC_sq = ( x + lambda_BC ) * ( x + lambda_BC )
		+ ( y - 1. + lambda_BC ) * ( y - 1. + lambda_BC )
		+ ( z + 1. - 2.*lambda_BC ) * ( z + 1. - 2.*lambda_BC );
	Function dist_BD_sq = x*x + (z+1.)*(z+1.);
	Function lambda_CD = ( - x + y + 2.*z - 3. ) / 6.;
	Function dist_CD_sq = ( x + 1. + lambda_CD ) * ( x + 1. + lambda_CD )
		+ ( y - lambda_CD ) * ( y - lambda_CD )
		+ ( z - 1. - 2.*lambda_CD ) * ( z - 1. - 2.*lambda_CD );
	Function smd = smooth_min ( dist_AB_sq, dist_AC_sq, dist_AD_sq,
	                            dist_BC_sq, dist_BD_sq, dist_CD_sq, tag::threshold, 0.02 );
	Manifold tetra = RR3.implicit ( smd == 0.02 );
	// the product of distances could be a good example for non-uniform mesh generation
	// ( dist_AB_sq * dist_AC_sq * dist_AD_sq * dist_BC_sq * dist_BD_sq * dist_CD_sq == 0.2 );

	Cell A ( tag::vertex );  x(A) =  0.11;  y(A) = 0.   ;  z(A) = 1.2;
	Cell B ( tag::vertex );  x(B) =  0.09;  y(B) = 0.   ;  z(B) = 1.2;
	Cell C ( tag::vertex );  x(C) =  0.1 ;  y(C) = 0.017;  z(C) = 1.197;
	tetra.project(A);  tetra.project(B);  tetra.project(C);
	Cell AB ( tag::segment, A.reverse(), B );
	Cell BC ( tag::segment, B.reverse(), C );
	Cell CA ( tag::segment, C.reverse(), A );
	Mesh chain ( tag::of_dim_one );
	AB.add_to ( chain );
	BC.add_to ( chain );
	CA.add_to ( chain );

	std::vector < double > N { 0, -0.02, 0. };
	Mesh msh ( tag::progressive, tag::boundary, chain,
	           tag::start_at, AB, tag::towards, N,
	           tag::desired_length, 0.02              );

	msh.export_msh ("msh-new.msh");
	std::cout << "reached end, " << chain.core->cells[1].size() << " segments on the interface" << std::endl;

	return 0;
}

void Cell::print_coords ( )
{	assert ( this->dim() == 0 );
	for ( size_t i = 0; i < 2; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		std::cout << x(*this) << " ";                       }  }

