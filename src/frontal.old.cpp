
// frontal.cpp 2022.05.02

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

// this is not the most elegant programming style
// I use some global variables, static_cast and also some goto's ...

#include <stack>
#include "math.h"
#include <memory>
#include <random>

#include "maniFEM.h"
#include "metric-tree.h"

namespace maniFEM {

namespace tag

{	struct OrthogonalTo { };  static const OrthogonalTo orthogonal_to;
	struct StartWithNonExistentMesh { };
	static const StartWithNonExistentMesh start_with_non_existent_mesh;
	struct AtPoint { };  static const AtPoint at_point ;                 }

}  // end of  namespace maniFEM

using namespace maniFEM;

//-----------------------------------------------------------------------------------------


// we want to deal with Euclidian manifolds and implicit submanifolds, on one side
// quotient manifolds on the other side
// not easy to reach both goals while keeping some of the common code ... common

// for quotient manifolds, we use a dirty trick
// we provide as Point not just a Cell but a Cell together with a winding number
// each call to SqDist will set the "winning" winding of the second Point B
// relatively to the first one A
// (the winding number of a future segment AB where the minimum distance is achieved)

// perhaps change the template of MetricTree like
// template < typename Point, typename SqDist, typname RichPoint = Point >
// where RichPoint contains a winding number while Point does not
// we will have to provide two versions of SqDist :
//    SqDist ( const Point & A, const Point & B )
//    SqDist ( const Point & A, RichPoint & B )
// the latter sets the winding of B to that winding which achieves the minimum distance

// for computing this winding number,
// implement a discrete version of the steepest descent method
// rather than a blind search as in draw_ps with tag::windng

// or is it desirable to keep a winding number for each node in the cloud ?
// e.g. relative to its parent ? or to the root ?
// if yes, no need to introduce RichPoint
	
//-----------------------------------------------------------------------------------------


class Manifold::Type::Euclidian

{	public :

	class sq_dist;
	typedef Cell winding_cell;
	typedef MetricTree < Cell, sq_dist > metric_tree;

};
	
//-------------------------------------------------------------------------------------------------

	
class Manifold::Type::Euclidian::sq_dist

// a callable object returning the square of the distance between two points
// used for MetricTree, see paragraphs 12.10 and 12.11 in the manual

{	public :

	inline double operator() ( const Cell & A, const Cell & B )
	{	double res = 0.;
		const size_t nc = Manifold::working .coordinates() .nb_of_components();
		for ( size_t i = 0; i < nc; i++ )
		{	Function x = Manifold::working .coordinates() [i];
			double tmp = x(B) - x(A);
			res += tmp*tmp;                                  }
		return res;                                                            }

};  // end of  class Manifold::Type::Euclidian::sq_dist

//-----------------------------------------------------------------------------------------


class Manifold::Type::Quotient

{	public :

	class sq_dist;	
	typedef std::pair < Cell, Manifold::Action > winding_cell;
	typedef MetricTree < winding_cell, sq_dist > metric_tree;
};

// in the above, we don't really need a MetricTree of windng_cells
// MetricTree < Cell, sq_dist >  would do as well
// we only need the winding for searching close neighbours of a given vertex
// in that situation, we need 'sq_dist' to keep the "winning" winding, see below
	
//-------------------------------------------------------------------------------------------------

	
class Manifold::Type::Quotient::sq_dist

// a callable object returning the square of the distance between two points
// used for MetricTree, see paragraphs 12.10 and 12.11 in the manual

{	public :

	inline double dist2
	( const Cell & A, const Cell & B, const Manifold::Action & winding,
	  const Function & coords_Eu, const Function & coords_q            )
	{	double res = 0.;
		std::vector < double > coord_A = coords_Eu ( A ),
		                       coord_B = coords_q ( B, tag::winding, winding );
		const size_t nc = coord_A .size();
		assert ( nc == coord_B .size() );
		for ( size_t i = 0; i < nc; i++ )
		{	double tmp = coord_B [i] - coord_A [i];
			res += tmp*tmp;                        }
		return res;                                                            }

	// we provide as Point not just a Cell but a Cell together with a winding number
	// each call to SqDist will set the "winning" winding of the second Point B
	// relatively to the first one A
	// (the winding number of a future segment AB where the minimum distance is achieved)

	inline double operator() ( std::pair < Cell, Manifold::Action > A,
	                           std::pair < Cell, Manifold::Action > & B )
	{	Manifold space = Manifold::working;
		assert ( space .exists() );  // we use the current (quotient) manifold
		Manifold::Quotient * manif_q = tag::Util::assert_cast
			< Manifold::Core*, Manifold::Quotient* > ( space .core );
		assert ( manif_q );
		Function coords_q = space.coordinates();
		Manifold mani_Eu = manif_q->base_space;  // underlying Euclidian manifold
		Function coords_Eu = mani_Eu .coordinates();
		assert ( coords_Eu .nb_of_components() == 2 );
		Function x = coords_Eu [0],  y = coords_Eu [1];

		// the action group may have one or two generators
		assert ( manif_q->actions .size() == 2 );
		assert ( manif_q->winding_nbs .size() == 2 );
		Function::ActionGenerator g1 = manif_q->actions [0], g2 = manif_q->actions [1];
	
		std::vector < std::vector < short int > > directions
			{ { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 } };
		// declare global, here and in global.cpp draw_ps
		
		size_t unsuccessful_tries = 0;
		// we describe a sort of spiral
		// at the end, we stop after 10 unsuccessful rounds
		size_t size_of_round = 0;
		short int ii = 0, jj = 0;
		// we choose a high value of 'dist', it will be overridden for ii == jj == 0
		double dist = this->dist2 ( A .first, B .first, 0, coords_Eu, coords_q ) + 1.;
		while ( unsuccessful_tries < 10 )
		{	size_of_round++;
			for ( size_t d = 0; d < 4; d++ )
			{	if ( d == 2 ) size_of_round++;
				for ( size_t i = 0; i < size_of_round; i++ )
				{	Manifold::Action s = ii*g1 + jj*g2;
					double di = this->dist2 ( A .first, B .first, s, coords_Eu, coords_q );
					if ( di < dist )
					{	dist = di;
						B .second = s;
						unsuccessful_tries = 0;  }
					ii += directions [d] [0];
					jj += directions [d] [1];                                               }  }
			unsuccessful_tries++;
		}  // end of  while unsuccessful_tries < 10

	return dist;                                                                          }

};  // end of  class Manifold::Type::Quotient::sq_dist

//-----------------------------------------------------------------------------------------


// global variables and functions for this file, not visible for other object files
namespace {  // anonymous namespace, mimics static linkage

// std::map < Cell, MetricTree<Cell,Manifold::Euclid::SqDist>::Node * > node_in_cloud;

Cell temporary_vertex ( tag::non_existent );

size_t progress_nb_of_coords;  // dimension of the surrounding Euclidian space
// also known as "geometric dimension"

const double one_plus_tolerance = 1.15;
const double sqrt_of_075 = std::sqrt ( 0.75 );
const double half_of_sqrt_of_075 = sqrt_of_075 / 2.;
Function desired_length = 0.;
double progress_long_dist, progress_sq_long_dist;
double desired_len_at_point, sq_desired_len_at_point;

Mesh mesh_under_constr ( tag::non_existent );
Mesh progress_interface ( tag::non_existent );
// empty temporary meshes
// these variables will be set when we start the meshing process,
// in the body of the progressive Mesh constructor

//-----------------------------------------------------------------------------------------------


inline double approx_sqrt ( double arg, const tag::Around &, double centre, double ini )
// hidden in anonymous namespace
	
// we assume all segments have length around the same value, here called 'ini'
// so the result should not be too far from 'ini'

// centre == ini*ini
	
{	double res = 0.5 * ( ini + arg / ini );
	res = 0.5 * ( res + arg / res );
	#ifndef NDEBUG
	if ( ( res < 0.25*ini ) or ( res > 4.*ini ) )
	{	std::cout << "bad approximation of square root" << std::endl;
		std::cout << arg << " " << ini*ini << std::endl;
		std::cout << res << " " << ini << std::endl;
		exit(1);                                                         }
	#endif  // debug
	return res;
}

//-----------------------------------------------------------------------------------------------


inline size_t get_topological_dim ( )       // hidden in anonymous namespace
	
{	Manifold::Implicit::OneEquation * m_impl_1 =
		dynamic_cast < Manifold::Implicit::OneEquation* > ( Manifold::working .core );
	if ( m_impl_1 )
	{	Manifold::Euclid * m_euclid = dynamic_cast < Manifold::Euclid* >
			( m_impl_1->surrounding_space .core );
		assert ( m_euclid );
		return m_euclid->coord_func .nb_of_components() - 1;              }
	Manifold::Implicit::TwoEquations * m_impl_2 =
		dynamic_cast < Manifold::Implicit::TwoEquations* > ( Manifold::working .core );
	if ( m_impl_2 )
		{	Manifold::Euclid * m_euclid = tag::Util::assert_cast
				< Manifold::Core*, Manifold::Euclid* > ( m_impl_2->surrounding_space .core );
		return m_euclid->coord_func .nb_of_components() - 2;                              }
	assert ( false );  // we return zero just to avoid compilation errors
	return 0;                                                                             }

//-----------------------------------------------------------------------------------------------


std::vector < double > compute_tangent_vec      // hidden in anonymous namespace
(	Cell start, bool check_orth, std::vector < double > given_vec )
	
// computes a vector tangent to Manifold::working at point 'start'

// if third argument is true, candidates will be projected onto the space orthogonal to given_vec
// given_vec must be tangent to Manifold::working at point 'start'
// and must have length approximately equal to desired_length
	
{	// Manifold::Implicit * m_impl =  dynamic_cast<Manifold::Implicit*> ( Manifold::working.core );
	// assert ( m_impl );
	// Manifold::Euclid * m_euclid =
	// 	dynamic_cast<Manifold::Euclid*> ( m_impl->surrounding_space.core );
	// assert ( m_euclid );
	std::vector < double > best_tangent;  double longest_projection = 0.;
	for ( size_t n = 1; n <= progress_nb_of_coords; n++ )
	{	const double coef = desired_len_at_point / std::sqrt(n);
		// we make sums of n vectors in the canonical basis with both signs
		std::vector < size_t > indices ( n+1 );
		for ( size_t i = 0; i < n; i++ )  indices [i] = i;
		indices [n] = progress_nb_of_coords;
		while ( true )
		{	std::vector < short int > signs ( n, 1. );
			while ( true )
			{	std::vector < double > tangent ( progress_nb_of_coords, 0. );
				for ( size_t i = 0; i < n; i++ )  tangent [ indices[i] ] = signs [i];
				// we normalize 'tangent'
				for ( size_t i = 0; i < progress_nb_of_coords; i++ ) tangent [i] *= coef;
				// we project
				for ( size_t i = 0; i < progress_nb_of_coords; i++ )
				{	Function x = Manifold::working .coordinates() [i];
					x ( temporary_vertex ) = x ( start ) + tangent [i];  }
				Manifold::working.project ( temporary_vertex );
				for ( size_t i = 0; i < progress_nb_of_coords; i++ )
				{	Function x = Manifold::working .coordinates() [i];
					tangent [i] = x ( temporary_vertex ) - x ( start );  }
				if ( check_orth )
				{	double prod = Manifold::working.inner_prod ( start, tangent, given_vec );
					double lambd = prod / sq_desired_len_at_point;
					for ( size_t i = 0; i < progress_nb_of_coords; i++ )
						tangent [i] -= lambd * given_vec [i];                                      }
				// we choose the longest projection
				double n2 = Manifold::working.inner_prod ( start, tangent, tangent );
				if ( n2 > longest_projection )
				{	best_tangent = tangent;  longest_projection = n2;  }
				// now change signs
				bool found = false;
				for ( int i = n-1; i >= 0; i-- )
					if ( signs [i] == 1 )
					{	found = true;  signs [i] = -1;
						for ( size_t j = i+1; j < n; j++ ) signs [j] = 1;
						break;                                            }
				if ( not found ) break;                                                  }
 			// now change indices
			bool found = false;
			for ( int i = n-1; i >= 0; i-- )
				if ( indices [i] < indices [i+1] - 1 )
				{	found = true;
					indices [i] ++;
					for ( size_t j = i+1; j<n; j++ )
					{	indices [j] = indices [j-1] + 1;
						assert ( indices [j] < progress_nb_of_coords );     }
					break;                                                            }
			if ( not found ) break;                                                  }
	}  // end of  for n
	// normalize best_tangent
	double n2 = Manifold::working.inner_prod ( start, best_tangent, best_tangent );
	double norm = approx_sqrt ( n2, tag::around, sq_desired_len_at_point, desired_len_at_point );
	double coef = desired_len_at_point / approx_sqrt ( n2, tag::around, norm*norm, norm );
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )  best_tangent [i] *= coef;
	return best_tangent;

}  // end of  compute_tangent_vec

//-----------------------------------------------------------------------------------------------


inline std::vector < double > compute_tangent_vec
(	const tag::AtPoint &, Cell start, const tag::OrthogonalTo &, std::vector < double > given_vec )
// hidden in anonymous namespace

// computes a vector tangent to Manifold::working at point 'start', normal to given_vec

// given_vec must be tangent to Manifold::working at point 'start'
// and must have length approximately equal to desired_length

{	return compute_tangent_vec ( start, true, given_vec );  }
	// 'true' as second argument means "do check orthogonality"
	
//-----------------------------------------------------------------------------------------------


inline std::vector < double > compute_tangent_vec ( const tag::AtPoint &, Cell start )
// hidden in anonymous namespace

// computes a vector tangent to Manifold::working at point 'start'

{	return compute_tangent_vec ( start, false, std::vector<double>() );  }
	// 'false' as second argument means "do not check orthogonality"

//----------------------------------------------------------------------------------------------


template < class manif_type >
inline void progress_add_point ( const typename manif_type::winding_cell & P,
                                 typename manif_type::metric_tree & cloud      )
// hidden in anonymous namespace

{	assert ( P.dim() == 0 );
	P.core->hook[tag::node_in_cloud] = static_cast < void * > ( cloud.add ( P ) );  }
// optimize !

//-----------------------------------------------------------------------------------------------


inline bool positive_orientation
(	const Cell & A, const Cell & B, const Cell & AB, const Cell & BC )
// hidden in anonymous namespace
	
{	assert ( A == AB.base().reverse() );
	assert ( B == AB.tip() );
	assert ( B == BC.base().reverse() );
	// code below is identical to part of progress_cos_sq_120
	// if you change anything, please change both; keep them identical
	assert ( BC.core->hook.find(tag::normal_vector) != BC.core->hook.end() );
	std::vector < double > & f =
		* static_cast < std::vector < double > * > ( BC.core->hook[tag::normal_vector] );
	std::vector < double >  e ( progress_nb_of_coords );
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		e[i] = x ( B ) - x ( A );          }  // recover tangent from hook !
	double prod = 0.;  // scalar product e.f
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )  prod += e[i] * f[i];
	// we could have used the Riemannian product, but the sign should be the same
	// code above is identical to part of progress_cos_sq_120
	return prod < 0.;                                                                   }

//-----------------------------------------------------------------------------------------------


inline double progress_cos_sq_60
(	const Cell & A, const Cell & B, const Cell & C, const Cell & AB, const Cell & BC )
// hidden in anonymous namespace

// return cosine square of 180 - ABC (or 0. if wrong orientation)
// calling function should check  cos_sq > 0.03  to ensure angle ABC is below 80 deg

{	assert ( A == AB.base().reverse() );
	assert ( B == AB.tip() );
	assert ( B == BC.base().reverse() );
	assert ( C == BC.tip() );
	// first, check orientation
	std::vector < double > e1 ( progress_nb_of_coords ),
	                       e2 ( progress_nb_of_coords );
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		e1[i] = x ( B ) - x ( A );  // recover tangents from hook !
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

//---------------------------------------------------------------------------------------------


inline double progress_cos_sq_120
(	const Cell & A, const Cell & B, const Cell & C, const Cell & AB, const Cell & BC )
// hidden in anonymous namespace

// return cosine square of 180 - ABC (or 2. if wrong orientation)
// calling function should check  cos_sq < 0.671  to ensure angle ABC is below 145 deg
// or check  cos_sq < -0.17365  to ensure angle ABC is below 80 deg [?!]

{	assert ( A == AB.base().reverse() );
	assert ( B == AB.tip() );
	assert ( B == BC.base().reverse() );
	assert ( C == BC.tip() );
	// first, check orientation
	// code below is identical to part of 'positive_orientation'
	// if you change anything, please change both; keep them identical
	assert ( BC.core->hook.find(tag::normal_vector) != BC.core->hook.end() );
	std::vector < double > & f =
		* static_cast < std::vector < double > * > ( BC.core->hook[tag::normal_vector] );
	std::vector < double >  e ( progress_nb_of_coords );
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		e[i] = x ( B ) - x ( A );          }  // recover tangent from hook !
	double prod = 0.;  // scalar product e.f
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )  prod += e[i] * f[i];
	// we could have used the Riemannian product, but the sign should be the same
	// code above is identical to part of 'positive_orientation'
	if ( prod > 0. ) return 2.;
	std::vector < double > e2 ( progress_nb_of_coords );
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	Function x = Manifold::working.coordinates()[i];  // recover tangent from hook !
		e2[i] = x ( C ) - x ( B );                       }
	// use Riemannian metric below
	double norm1 = 0., norm2 = 0.;  prod = 0.;
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	norm1 += e [i] * e [i];
		norm2 += e2[i] * e2[i];
		prod  += e [i] * e2[i];  }
	return prod * prod / norm1 / norm2;  // cosine square
}

//-----------------------------------------------------------------------------------------------
			

inline bool opposite_signs ( const double a, const double b )
// hidden in anonymous namespace
{	if ( a < 0. )  return b >= 0.;
	if ( a == 0. )  return true;
	return b <= 0.;                  }

//----------------------------------------------------------------------------------------------


inline Cell search_start_ver_c1 ( )    // hidden in anonymous namespace

// search for a starting point in a manifold of co-dimension one
// e.g. a curve in the plane or a surface in 3D

{	Manifold::Implicit::OneEquation * m_impl =
		dynamic_cast<Manifold::Implicit::OneEquation*> ( Manifold::working.core );
	assert ( m_impl );
	// Manifold::Euclid * m_euclid =
	// 	dynamic_cast<Manifold::Euclid*> ( m_impl->surrounding_space.core );
	// assert ( m_euclid );
	// Function m_impl->level_function;

	// this function is not called at some point in the manifold
	// (thats's what it is here for, to find a starting point)
	// so we have no desired_len_at_point, we must invent some length
	// we could take desired_length at the origin of the space
	// but in order to reduce the risk of falling into a singularity
	// we prefer to introduce some randomness
	std::default_random_engine random_generator;
	std::uniform_real_distribution<double> distr11 ( -1., 1. );
	Cell tmp_ver ( tag::vertex );
	Cell tmp_ver_1 ( tag::vertex );
	Cell tmp_ver_2 ( tag::vertex );
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		x ( tmp_ver ) = distr11(random_generator);        }
	desired_len_at_point = desired_length ( tmp_ver );
	sq_desired_len_at_point = desired_len_at_point * desired_len_at_point;
	
	size_t size_of_cube = 5;
	while ( true )
	{	double s = size_of_cube * desired_len_at_point;
		std::uniform_real_distribution<double> distr ( -s, s );
		size_t nb = 1;
		for ( size_t j = 0; j < progress_nb_of_coords; j++ ) nb *= size_of_cube;
		for ( size_t j = 0; j < nb; j++ )
		{	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
			{	Function x = Manifold::working.coordinates()[i];
				x ( tmp_ver_1 ) = distr(random_generator);
				x ( tmp_ver_2 ) = distr(random_generator);        }
			double v1 = m_impl->level_function ( tmp_ver_1 ),
			       v2 = m_impl->level_function ( tmp_ver_2 );
			if ( opposite_signs ( v1, v2 ) )
				// refine by applying bissection algorithm
				while ( true )
				{	if ( Manifold::working.dist_sq ( tmp_ver_1, tmp_ver_2 ) < sq_desired_len_at_point )
					{	Manifold::working.project ( tmp_ver_1 );
						return tmp_ver_1;                        }
					m_impl->surrounding_space.interpolate ( tmp_ver, 0.5, tmp_ver_1, 0.5, tmp_ver_2 );
					double v = m_impl->level_function ( tmp_ver );
					if ( opposite_signs ( v, v2 ) )
					{	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
						{	Function x = Manifold::working.coordinates()[i];
							x ( tmp_ver_1 ) = x ( tmp_ver );                 }
						v1 = v;                                               }
					else if ( opposite_signs ( v, v1 ) )
					{	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
						{	Function x = Manifold::working.coordinates()[i];
							x ( tmp_ver_2 ) = x ( tmp_ver );                 }
						v2 = v;                                               }
					else assert( false );                                               }  }
		size_of_cube *= 2;                                                                          }
}

//----------------------------------------------------------------------------------------------


inline double ext_prod_R2 ( const double vx, const double vy, const double wx, const double wy )
// hidden in anonymous namespace
{	return vx*wy - wx*vy;  }


//----------------------------------------------------------------------------------------------


inline bool origin_outside ( const double & Ax, const double & Ay,
                             const double & Bx, const double & By,
                             const double & Cx, const double & Cy )
// hidden in anonymous namespace

{	return
	   opposite_signs ( - ext_prod_R2 ( Bx-Ax, By-Ay, Ax, Ay ),
	                      ext_prod_R2 ( Bx-Ax, By-Ay, Cx-Ax, Cy-Ay ) )
	or opposite_signs ( - ext_prod_R2 ( Cx-Bx, Cy-By, Bx, By ),
	                      ext_prod_R2 ( Cx-Bx, Cy-By, Ax-Bx, Ay-By ) )
	or opposite_signs ( - ext_prod_R2 ( Ax-Cx, Ay-Cy, Cx, Cy ),
	                      ext_prod_R2 ( Ax-Cx, Ay-Cy, Bx-Cx, By-Cy ) ); }

//-----------------------------------------------------------------------------------------------


inline Cell search_start_ver_c2 ( )
// hidden in anonymous namespace

// search for a starting point in a manifold of co-dimension two (a curve in 3D)

{	Manifold::Implicit::TwoEquations * m_impl =
		dynamic_cast<Manifold::Implicit::TwoEquations*> ( Manifold::working.core );
	assert ( m_impl );
	// Manifold::Euclid * m_euclid =
	// 	dynamic_cast<Manifold::Euclid*> ( m_impl->surrounding_space.core );
	// assert ( m_euclid );

	// this function is not called at some point in the manifold
	// (thats's what it is here for, to find a starting point)
	// so we have no desired_len_at_point, we must invent some length
	// we could take desired_length at the origin of the space
	// but in order to reduce the risk of falling into a singularity
	// we prefer to introduce some randomness
	std::default_random_engine random_generator;
	std::uniform_real_distribution<double> distr11 ( -1., 1. );
	Cell tmp_A ( tag::vertex );
	Cell tmp_B ( tag::vertex );
	Cell tmp_C ( tag::vertex );
	Cell tmp_AB ( tag::vertex );
	Cell tmp_BC ( tag::vertex );
	Cell tmp_CA ( tag::vertex );
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		x ( tmp_A ) = distr11(random_generator);        }
	desired_len_at_point = desired_length ( tmp_A );
	sq_desired_len_at_point = desired_len_at_point * desired_len_at_point;
	
	size_t size_of_cube = 5;
	while ( true )
	{	double s = size_of_cube * desired_len_at_point;
		std::uniform_real_distribution<double> distr ( -s, s );
		size_t nb = 1;
		restart :
		for ( size_t j = 0; j < progress_nb_of_coords; j++ ) nb *= size_of_cube;
		for ( size_t j = 0; j < nb; j++ )
		{	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
			{	Function x = Manifold::working.coordinates()[i];
				x ( tmp_A ) = distr(random_generator);
				x ( tmp_B ) = distr(random_generator);
				x ( tmp_C ) = distr(random_generator);             }
			double vA1 = m_impl->level_function_1 ( tmp_A ),
			       vA2 = m_impl->level_function_2 ( tmp_A );
			double vB1 = m_impl->level_function_1 ( tmp_B ),
			       vB2 = m_impl->level_function_2 ( tmp_B );
			double vC1 = m_impl->level_function_1 ( tmp_C ),
			       vC2 = m_impl->level_function_2 ( tmp_C );
			if ( not origin_outside ( vA1, vA2, vB1, vB2, vC1, vC2 ) )
				// refine by repeatedly cutting the triangle
				while ( true )
				{	if ( ( Manifold::working.dist_sq ( tmp_A, tmp_B ) < sq_desired_len_at_point ) and
					     ( Manifold::working.dist_sq ( tmp_B, tmp_C ) < sq_desired_len_at_point ) and
					     ( Manifold::working.dist_sq ( tmp_C, tmp_A ) < sq_desired_len_at_point )      )
					{	Manifold::working.project ( tmp_A );
						return tmp_A;                        }
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
					{	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
						{	Function x = Manifold::working.coordinates()[i];
							x ( tmp_B ) = x ( tmp_AB );  x ( tmp_C ) = x ( tmp_CA );  }
						vB1 = vAB1;  vB2 = vAB2;  vC1 = vCA1;  vC2 = vCA2;            }
					else if ( not origin_outside ( vB1, vB2, vAB1, vAB2, vBC1, vBC2 ) )
					{	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
						{	Function x = Manifold::working.coordinates()[i];
							x ( tmp_A ) = x ( tmp_AB );  x ( tmp_C ) = x ( tmp_BC );  }
						vA1 = vAB1;  vA2 = vAB2;  vC1 = vBC1;  vC2 = vBC2;            }
					else if ( not origin_outside ( vC1, vC2, vBC1, vBC2, vCA1, vCA2 ) )
					{	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
						{	Function x = Manifold::working.coordinates()[i];
							x ( tmp_A ) = x ( tmp_CA );  x ( tmp_B ) = x ( tmp_BC );  }
						vA1 = vCA1;  vA2 = vCA2;  vB1 = vBC1;  vB2 = vBC2;            }
					else if ( not origin_outside ( vAB1, vAB2, vBC1, vBC2, vCA1, vCA2 ) )
					{	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
						{	Function x = Manifold::working.coordinates()[i];
							x ( tmp_A ) = x ( tmp_BC );  x ( tmp_B ) = x ( tmp_CA );
							x ( tmp_C ) = x ( tmp_AB );                               }
						vA1 = vBC1;  vA2 = vBC2;  vB1 = vCA1;  vB2 = vCA2;
						vC1 = vAB1;  vC2 = vAB2;                                                 }
					else  // nasty nonlinear level functions ...
						goto restart;                                                               }  }
		size_of_cube *= 2;                                                                          }
}

//-----------------------------------------------------------------------------------------------


inline Cell search_start_ver ( )
// hidden in anonymous namespace

// search for a starting point
// current working manifold may have co-dimension one
// e.g. a curve in the plane or a surface in 3D
// of co-dimension two (a curve in 3D)

{	Manifold::Implicit::OneEquation * m_impl =
		dynamic_cast<Manifold::Implicit::OneEquation*> ( Manifold::working.core );
	if ( m_impl ) return search_start_ver_c1 ( );  // co-dimension one
	return search_start_ver_c2 ( );  // co-dimension two
}

//-----------------------------------------------------------------------------------------------


inline void redistribute_vertices                 // line 710
( const Mesh & msh, const Cell & start, const Cell & stop, size_t n )
// hidden in anonymous namespace    // called only once

// just make some baricenters
	
{	Cell A = stop;
	// just in case n is too large, or the curve is too short, we look for 'start'
	// the statement below is important for closed loops
	A = msh .cell_behind (A) .base() .reverse();
	for ( size_t i = 2; i < n; i++ )
	{	if ( A == start )  {  n = i;  break;  }
		A = msh .cell_behind ( A, tag::surely_exists ) .base() .reverse();  }
	assert ( n > 2 );
	A = stop;
	Cell B = msh .cell_behind ( A, tag::surely_exists ) .base() .reverse();
	Cell C = msh .cell_behind ( B, tag::surely_exists ) .base() .reverse();
	for ( size_t i = 2; i < n; i++ )
	{	Manifold::working .interpolate ( B, 0.3, A, 0.4, B, 0.3, C );
		if ( C == start ) break;
		A = B;  B = C;
		C = msh .cell_behind ( B, tag::surely_exists ) .base() .reverse(); }    }

//-----------------------------------------------------------------------------------------------


inline double get_z_baric ( const Cell & tri )
// hidden in anonymous namespace

{	assert ( progress_nb_of_coords == 3 );
	Mesh::Iterator it = tri.boundary().iterator ( tag::over_vertices );
	Function z = Manifold::working.coordinates()[2];
	double zz = 0.;
	size_t counter = 0;
	for ( it.reset(); it.in_range(); it++, counter++ )  zz += z(*it);
	assert ( counter == 3 );
	return  zz/3.;                                                               }

//-----------------------------------------------------------------------------------------------


inline bool tri_correctly_oriented ( const Cell & tri )
// hidden in anonymous namespace

{	assert ( tri.dim() == 2 );

	Mesh::Iterator it = tri.boundary().iterator ( tag::over_segments, tag::require_order );
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
	return  xAB * yBC > yAB * xBC;;                                                         }

//-------------------------------------------------------------------------------------------------


bool correctly_oriented    // hidden in anonymous namespace
( const Mesh msh, const tag::Orientation &, const tag::OrientationChoice & oc )

// tells whether 'msh's orientation is consistent with the orientation of the
// surrounding Euclidian space

{	std::string error_message = "????";
	if ( oc == tag::inherent )
		error_message = "inherent orientation only makes sense if co-dimension is one";
	else
	{	assert ( oc == tag::not_provided );
		error_message = "please provide tag::orientation, tag::random";  }
	Manifold::Implicit::OneEquation * m_impl =
		dynamic_cast < Manifold::Implicit::OneEquation* > ( Manifold::working.core );
	if ( m_impl == nullptr ) 
	{	std::cout << error_message << std::endl;
		exit (1);                               }
	Manifold::Euclid * m_euclid =
		dynamic_cast < Manifold::Euclid* > ( m_impl->surrounding_space.core );
	if ( m_euclid == nullptr ) 
	{	std::cout << error_message << std::endl;
		exit (1);                               }

	// for surfaces, we search the vertex with zmax and check the orientation
	// of all surrounding triangles (we need an iterator over cells around that vertex)
	if ( msh.dim() != 1 )
	{	assert ( msh.dim() == 2 );
		assert ( progress_nb_of_coords == 3 );
		Mesh::Iterator it = msh.iterator ( tag::over_cells_of_max_dim );
		it.reset();  assert ( it.in_range() );
		Cell trimax = *it;
		double zmax = get_z_baric ( trimax );
		for ( it++; it.in_range(); it++ )
		{	double zz = get_z_baric ( *it );
			if ( zz > zmax )
			{	zmax = zz;  trimax = *it;  }   }
		return tri_correctly_oriented ( trimax );                     }

	Function x = Manifold::working.coordinates()[0];
	Function y = Manifold::working.coordinates()[1];

	Mesh::Iterator it = msh.iterator ( tag::over_vertices );
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
	
}  // end of correctly_oriented

//-------------------------------------------------------------------------------------------------

// a more complicated version of correctly_oriented can be found at
// https://github.com/cristian-barbarosie/attic - manifem.cpp

//-------------------------------------------------------------------------------------------------


inline void switch_orientation_direct ( Mesh & msh )
// hidden in anonymous namespace
	
{	std::vector < Cell > vec_of_cells;
	vec_of_cells .reserve ( msh .number_of ( tag::cells_of_max_dim ) );
	Mesh::Iterator itt = msh .iterator ( tag::over_cells_of_max_dim );
	for ( itt .reset(); itt .in_range(); itt++ )
		vec_of_cells .push_back ( *itt );
	for ( std::vector<Cell>::iterator it = vec_of_cells.begin(); it != vec_of_cells.end(); it++ )
		it->remove_from_mesh ( msh, tag::do_not_bother );
	for ( std::vector<Cell>::iterator it = vec_of_cells.begin(); it != vec_of_cells.end(); it++ )
		it->reverse() .add_to_mesh ( msh, tag::do_not_bother );                                     }
// the meaning of tag::do_not_bother is explained at the end of paragraph 11.6 in the manual
// no need to call update_info_connected_one_dim because we are at hands with a closed loop
// nb_of_segs remains the same, as well as first_ver and last_ver

//-------------------------------------------------------------------------------------------------


void switch_orientation_of_each_cell ( Mesh & msh )
// hidden in anonymous namespace

// 'msh' is a closed mesh (has no boundary)
	
// since 'msh' has just been created, we choose not to use 'reverse'
// instead, call switch_orientation_direct on boundary of each cell

// does not work for one-dimensional meshes
// because switch_orientation_direct does not work for segments
// should be adapted

{	assert ( msh .dim() == 2 );
	std::vector < Cell > vec_of_cells;
	vec_of_cells .reserve ( msh .number_of ( tag::cells_of_max_dim ) );
	Mesh::Iterator itt = msh .iterator ( tag::over_cells_of_max_dim );
	for ( itt .reset(); itt .in_range(); itt++ )
		vec_of_cells .push_back ( *itt );
	for ( std::vector<Cell>::iterator it = vec_of_cells.begin(); it != vec_of_cells.end(); it++ )
		it->remove_from_mesh ( msh, tag::do_not_bother );
	for ( std::vector<Cell>::iterator it = vec_of_cells.begin(); it != vec_of_cells.end(); it++ )
	{	Mesh bdry = it->boundary();
		switch_orientation_direct ( bdry );  }
	for ( std::vector<Cell>::iterator it = vec_of_cells.begin(); it != vec_of_cells.end(); it++ )
		it->add_to_mesh ( msh, tag::do_not_bother );                                                }
// the meaning of tag::do_not_bother is explained at the end of paragraph 11.6 in the manual
// no need to call update_info_connected_one_dim because we are at hands with a closed loop
// nb_of_segs remains the same, as well as first_ver and last_ver

//-------------------------------------------------------------------------------------------------


inline void improve_tangent ( const Cell & A, std::vector < double > & nor )
// hidden in anonymous namespace

// project 'nor' onto working manifold

{	// we normalize 'nor'
	double n2 = Manifold::working.inner_prod ( A, nor, nor );
	double norm = approx_sqrt ( n2, tag::around, sq_desired_len_at_point, desired_len_at_point );
	norm = approx_sqrt ( n2, tag::around, norm*norm, norm );
	norm = desired_len_at_point / norm;
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
// hidden in anonymous namespace

// project 'nor' onto working manifold and normalize it relatively to segment AB

{	// we normalize 'nor'
	double n2 = Manifold::working.inner_prod ( A, nor, nor );
	double norm = approx_sqrt ( n2, tag::around, sq_desired_len_at_point, desired_len_at_point );
	norm = approx_sqrt ( n2, tag::around, norm*norm, norm );
	norm = desired_len_at_point / norm;
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
	norm = approx_sqrt ( n2, tag::around, sq_desired_len_at_point, desired_len_at_point );
	norm = approx_sqrt ( n2, tag::around, norm*norm, norm );
	norm = desired_len_at_point / norm;
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )  nor[i] *= norm;                  }

//-------------------------------------------------------------------------------------------------


inline void build_one_normal ( Cell & B, Cell & C, Cell & new_seg )
// hidden in anonymous namespace

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
	// recover tangent from hook !
	assert ( AB.core->hook.find(tag::normal_vector) != AB.core->hook.end() );
	std::vector < double > & old_f =
		* static_cast < std::vector < double > * > ( AB.core->hook[tag::normal_vector] );
	// 'e' is the vector of the segment, 'f' is orthogonal
	// they are all of approximately the same length, equal to desired_length

	// code below is identical to part of 'build_each_normal'
	// if you change anything, please change both; keep them identical
	std::vector < double > vC = Manifold::working.coordinates() ( C );
	std::vector < double > new_e ( progress_nb_of_coords );
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )  new_e[i] = vC[i] - vB[i];
	// recover tangent from hook !
	// scalar products :
	double with_e = 0.,  with_f = 0., norm_e_sq = 0., norm_f_sq = 0.;
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )  // use Riemannian metric !
	{	with_e += new_e[i]*old_e[i];     with_f += new_e[i]*old_f[i];
		norm_e_sq += old_e[i]*old_e[i];  norm_f_sq += old_f[i]*old_f[i];  }
	with_e /= norm_e_sq;
	with_f /= norm_f_sq;
	std::vector < double > & new_f = * new std::vector < double > ( progress_nb_of_coords );
	// we rotate 'new_e' with 90 degrees, in the same sense as 'f' is rotated from 'e'
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
		new_f[i] = - with_f * old_e[i] + with_e * old_f[i];
	// project and normalize :
	improve_normal ( B, C, new_e, new_f );
	new_seg.core->hook[tag::normal_vector] = static_cast < void * > ( & new_f );  // optimize
	// code above is identical to part of 'build_each_normal'
}  // end of build_one_normal

//-------------------------------------------------------------------------------------------------


inline void build_each_normal     // hidden in anonymous namespace
( Cell & B, Cell & C, Cell & new_seg,
  std::vector < double > & old_e, std::vector < double > & old_f )

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
	// recover tangent from hetero_info !
	// scalar products :
	double with_e = 0.,  with_f = 0.;
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )  // use Riemannian metric !
	{	with_e += new_e[i]*old_e[i];  with_f += new_e[i]*old_f[i];  }
	with_e /= sq_desired_len_at_point;
	with_f /= sq_desired_len_at_point;
	std::vector < double > & new_f = * new std::vector < double > ( progress_nb_of_coords );
	// we rotate 'new_e' with 90 degrees, in the same sense as 'f' is rotated from 'e'
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
		new_f[i] = - with_f * old_e[i] + with_e * old_f[i];
	// project and normalize :
	improve_normal ( B, C, new_e, new_f );
	new_seg.core->hook[tag::normal_vector] = static_cast < void * > ( & new_f );  // optimize
	// code above is identical to part of 'build_one_normal'
	old_e = new_e;  old_f = new_f;
}  // end of build_each_normal

//-------------------------------------------------------------------------------------------------


inline Cell build_normals ( const Cell & start )
// hidden in anonymous namespace

// from a cell 'start', propagate normals along progress_interface
// (will only cover the connected component containing 'start')
// return the first segment which already has a normal
// see paragraph 12.6 in the manual
// 'normal' should have norm approximately equal to desired_length

{	assert ( start.belongs_to ( progress_interface, tag::same_dim, tag::oriented ) );
	#ifndef NDEBUG
	// std::cout << "building normals" << std::endl;
	#endif

	// we keep local copies of desired_len_at_point and sq_desired_len_at_point
	double dlp = desired_len_at_point, sdlp = sq_desired_len_at_point;

	Cell seg = start;
	Cell A = seg.base().reverse(),  B = seg.tip();
	std::vector < double > va = Manifold::working.coordinates() (A),
		vb = Manifold::working.coordinates() (B);
	assert ( seg.core->hook.find(tag::normal_vector) != seg.core->hook.end() );
	std::vector < double > f =
		* static_cast < std::vector < double > * > ( seg.core->hook[tag::normal_vector] );
	std::vector < double > e ( progress_nb_of_coords );
	for ( size_t i = 0; i < progress_nb_of_coords; i++ ) e[i] = vb[i] - va[i];
	// recover tangent from hook !
	// 'e' and 'f' form an oriented basis in the two-dimensional space
	// tangent to the manifold at the current point
	// they will be used to build further vectors pointing outwards
	// (on the correct side of progress_interface)
	while ( true )
	// progress_interface may be disconnected, so we cannot use Mesh::Iterators
	// this loop will only cover its current connected component
	{	desired_len_at_point = desired_length ( B );
		sq_desired_len_at_point = desired_len_at_point * desired_len_at_point;
 		Cell new_seg = progress_interface.cell_in_front_of ( B );
		if ( new_seg.core->hook.find(tag::normal_vector) != new_seg.core->hook.end() )
		{	desired_len_at_point = dlp;
			sq_desired_len_at_point = sdlp;
			return new_seg;                 }
		assert ( new_seg != start );
		Cell C = new_seg.tip();
		build_each_normal ( B, C, new_seg, e, f );
		// 'e' and 'f' get updated within 'build_each_normal'
		seg = new_seg;  B = C;                                                       }

}  // end of build_normals

//-------------------------------------------------------------------------------------------------


template < class manif_type >
inline void progress_fill_60
(	Cell & AB, Cell & BC, const Cell & CA, const Cell & B,
	typename manif_type::metric_tree & cloud               )
// hidden in anonymous namespace

{	AB.remove_from_mesh ( progress_interface );
	BC.remove_from_mesh ( progress_interface );
	// optimize map access in statements below !!
	assert ( AB.core->hook.find(tag::normal_vector) != AB.core->hook.end() );
	assert ( BC.core->hook.find(tag::normal_vector) != BC.core->hook.end() );
	delete static_cast < std::vector < double > * > ( AB.core->hook[tag::normal_vector] );
	delete static_cast < std::vector < double > * > ( BC.core->hook[tag::normal_vector] );
	AB.core->hook.erase ( tag::normal_vector );
	BC.core->hook.erase ( tag::normal_vector );
	cloud.remove ( static_cast < typename manif_type::metric_tree::Node * >
	               ( B.core->hook[tag::node_in_cloud] )                    );
	B.core->hook.erase ( tag::node_in_cloud );  // optimize !
	Cell new_tri ( tag::triangle, AB, BC, CA );
	
	new_tri.add_to_mesh ( mesh_under_constr );
	if ( B.is_inner_to ( mesh_under_constr ) ) mesh_under_constr.baricenter ( B );         }

//-------------------------------------------------------------------------------------------------


inline void glue_two_segs_common
( Cell & A, Cell & B, Cell & AB, Cell & C, Cell & D, Cell & CD, Cell & AD, Cell & BC )
// hidden in anonymous namespace

{	AB.remove_from_mesh ( progress_interface );
	CD.remove_from_mesh ( progress_interface );
	Cell CB = BC.reverse();
	AD.add_to_mesh ( progress_interface );
	CB.add_to_mesh ( progress_interface );
	// optimize map access in statements below !!
	assert ( AB.core->hook.find(tag::normal_vector) != AB.core->hook.end() );
	assert ( CD.core->hook.find(tag::normal_vector) != CD.core->hook.end() );
	delete static_cast < std::vector < double > * > ( AB.core->hook[tag::normal_vector] );
	delete static_cast < std::vector < double > * > ( CD.core->hook[tag::normal_vector] );
	AB.core->hook.erase ( tag::normal_vector );
	CD.core->hook.erase ( tag::normal_vector );
	// build normals for two newly added segments AD and CB
	build_one_normal ( A, D, AD );  // from previous segment
	build_one_normal ( C, B, CB );  // from previous segment
}  // end of void glue_two_segs_common

//-------------------------------------------------------------------------------------------------


template < class manif_type >
inline Cell glue_two_segs_S
(	Cell & A, Cell & B, Cell & AB, Cell & C, Cell & D, Cell & CD,
	typename manif_type::metric_tree & cloud                       )
// hidden in anonymous namespace

// see paragraph 12.8 in the manual

// progress_interface.cell_in_front_of(B) may have tip C
// that is, BC may belong already to 'progress_interface'

{	Cell AD ( tag::segment, A.reverse(), D );  // winding !!
	Cell AC ( tag::segment, A.reverse(), C );  // winding !!
	Cell BC = progress_interface.cell_in_front_of(B);
	if ( BC.tip() == C )
	{	progress_fill_60 < manif_type > ( AB, BC, AC.reverse(), B, cloud );
		AC.add_to_mesh ( progress_interface );
		build_one_normal ( A, C, AC );  // based on previous segment
		progress_fill_60 < manif_type > ( AC, CD, AD.reverse(), C, cloud );
		AD.add_to_mesh ( progress_interface );
		build_one_normal ( A, D, AD );  }  // based on previous segment
	else
	{	BC = Cell ( tag::segment, B.reverse(), C );  // winding !!
		Cell ABC ( tag::triangle, AB, BC, AC.reverse() );
		Cell ACD ( tag::triangle, AC, CD, AD.reverse() );
		ABC.add_to_mesh ( mesh_under_constr );
		ACD.add_to_mesh ( mesh_under_constr );
		glue_two_segs_common ( A, B, AB, C, D, CD, AD, BC );  }
	return AD;                                                                    }

//-------------------------------------------------------------------------------------------------


template < class manif_type >
inline Cell glue_two_segs_Z
(	Cell & A, Cell & B, Cell & AB, Cell & C, Cell & D, Cell & CD,
	typename manif_type::metric_tree & cloud                      )
// hidden in anonymous namespace

// see paragraph 12.8 in the manual

// progress_interface.cell_behind(A) may have base D
// that is, DA may belong already to 'progress_interface'

{	Cell BC ( tag::segment, B.reverse(), C );  // winding !!
	Cell DB ( tag::segment, D.reverse(), B );  // winding !!
	Cell DA = progress_interface.cell_behind(A);
	if ( DA.base().reverse() == D )
	{	progress_fill_60 < manif_type > ( DA, AB, DB.reverse(), A, cloud );
		DB.add_to_mesh ( progress_interface );
		build_one_normal ( D, B, DB );  // based on previous segment
		progress_fill_60 < manif_type > ( CD, DB, BC, D, cloud );
		BC.reverse().add_to_mesh ( progress_interface );
		build_one_normal ( B, C, BC );  }  // based on previous segment
	else
	{	DA = Cell ( tag::segment, D.reverse(), A );  // winding !!
		Cell ABD ( tag::triangle, AB, DB.reverse(), DA );
		Cell BCD ( tag::triangle, BC, CD, DB );
		ABD.add_to_mesh ( mesh_under_constr );
		BCD.add_to_mesh ( mesh_under_constr );
		Cell AD = DA.reverse();
		glue_two_segs_common ( A, B, AB, C, D, CD, AD, BC );  }
	return BC.reverse();                                                 }

//-------------------------------------------------------------------------------------------------


template < class manif_type >
inline void progress_fill_last_triangle
(	const Cell & A, const Cell & B, const Cell & C, Cell & AB, Cell & BC, Cell & CA,
	typename manif_type::metric_tree & cloud                                         )
// hidden in anonymous namespace

{	progress_fill_60 < manif_type > ( AB, BC, CA, B, cloud );
	CA.remove_from_mesh ( progress_interface );
	// optimize map access in statements below !!
	assert ( CA.core->hook.find(tag::normal_vector) != CA.core->hook.end() );
	delete static_cast < std::vector < double > * > ( CA.core->hook[tag::normal_vector] );
	CA.core->hook.erase ( tag::normal_vector );
	cloud.remove ( static_cast < typename manif_type::metric_tree::Node * >
	               ( A.core->hook[tag::node_in_cloud] )                     );
	A.core->hook.erase ( tag::node_in_cloud );  // optimize !
	cloud.remove ( static_cast < typename manif_type::metric_tree::Node * >
	               ( C.core->hook[tag::node_in_cloud] )                     );
	C.core->hook.erase ( tag::node_in_cloud );  // optimize !
	if ( A.is_inner_to ( mesh_under_constr ) ) mesh_under_constr.baricenter ( A );
	if ( C.is_inner_to ( mesh_under_constr ) ) mesh_under_constr.baricenter ( C );         }

//-------------------------------------------------------------------------------------------------


template < class manif_type >
void progress_relocate
(	const Cell & P, size_t n, std::vector<double> & normal_dir,
	std::set < typename manif_type::winding_cell > & set_of_ver,
	typename manif_type::metric_tree & cloud                      )
// hidden in anonymous namespace

// modifies normal_dir
	
// re-compute the placement of a newly created vertex

// vertex has been located according to two segments, from angles_120 :   n == 2
// or according to only one segment, if built from a brand new triangle : n == 1
	
// we compute here 'set_of_ver'
// (which is the set of all vertices in the cloud close enough to 'ver')
// and keep it for future use in 'check_touching'

{	// make a list (using a set) of nearby points
	// build a vector of segments from it
	// relocate point P by averaging all normals

	std::list < typename manif_type::winding_cell > list_of_ver =
		cloud .find_close_neighbours_of ( P, progress_long_dist );
	// P has not been added to the cloud yet, so it will not show up in 'list_of_ver'
	set_of_ver .clear();
	for ( typename std::list < typename  manif_type::winding_cell >
					::const_iterator it = list_of_ver .begin();
        it != list_of_ver .end(); it++                            )
		set_of_ver .insert ( *it );

	std::vector < Cell > vector_of_seg;
	// vector_of_seg will contain all segments whose both extremities belong to 'set_of-ver'
	// but not the segments adjacent to P (since P does not belong to 'set_of_ver')
	for ( std::set < Cell > ::iterator it = set_of_ver .begin(); it != set_of_ver .end(); it++ )
	{	Cell A = *it;
		Cell AB = progress_interface .cell_in_front_of ( A );
		if ( set_of_ver .find ( AB .tip() ) != set_of_ver .end() )
			vector_of_seg .push_back ( AB );                        }

	if ( vector_of_seg .empty() ) return;

	// there are two cases : 
	// the piece of the interface that we just encountered may have normals or not
	size_t counter = 0;
	Cell kept_seg ( tag::non_existent );
	for ( size_t i = 0; i < vector_of_seg .size(); i++ )
	{	Cell seg_p = vector_of_seg [i];
		if ( seg_p .core->hook .find ( tag::normal_vector ) == seg_p .core->hook .end() )
		{	counter++;  kept_seg = seg_p;  }                                               }
	if ( counter > 0 )  // there are 'counter' segments with no normal
	{	// assert ( counter == 1 );
		// build normal of 'kept_seg' from 'normal_dir'
		Cell A = kept_seg .base() .reverse();
		Cell B = kept_seg .tip();
		std::vector < double > tangent_dir ( progress_nb_of_coords );
		for ( size_t i = 0; i < progress_nb_of_coords; i++ )
		{	Function x = Manifold::working.coordinates()[i];
			tangent_dir [i] = x(B) - x(A);  // recover tangent from hook !
			normal_dir [i] *= -1.;                             }
	  improve_normal ( A, B, tangent_dir, normal_dir );  // modifies normal_dir
		kept_seg .core->hook [ tag::normal_vector ] = static_cast < void * >
			( new std::vector < double > { normal_dir } );  // optimize !!
		Cell ret = build_normals ( kept_seg );
		assert ( ret == kept_seg );                                           }

	std::vector < double > pos = Manifold::working.coordinates() ( P );
	if ( n != 1 )
	{	assert ( n == 2 );
		for ( size_t i = 0; i < progress_nb_of_coords; i++ )  pos[i] *= n;  }
	for ( size_t j = 0; j < vector_of_seg.size(); j++ )
	{	Cell AB = vector_of_seg[j];
		assert ( AB.core->hook.find(tag::normal_vector) != AB.core->hook.end() );
	  std::vector < double > & nor =
			* static_cast < std::vector < double > * > ( AB.core->hook[tag::normal_vector] );
		Cell A = AB.base().reverse();
		Cell B = AB.tip();
		for ( size_t i = 0; i < progress_nb_of_coords; i++ )
		{	Function x = Manifold::working.coordinates()[i];
			pos[i] += ( x(A) + x(B) ) / 2. + nor[i] * sqrt_of_075;  }  }
	n += vector_of_seg.size();
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		x(P) = pos[i] / n;                                 }
	
	Manifold::working.project(P);

	// check P is on the right side of kept_seg

}  // end of  progress_relocate

//-------------------------------------------------------------------------------------------------


template < class manif_type >
inline bool check_touching       // hidden in anonymous namespace
(	Cell & ver, std::set < typename manif_type::winding_cell > & set_of_ver,
	Cell & point_120, Cell & stop_point_120,
	typename manif_type::metric_tree & cloud                                  )

// analyse position of recently built vertex 'ver' relatively to other vertices on the interface
// occasionally, different connected components of the interface touch and merge
// or, the current connected component of the interface may touch itself and split in two

// return true if a touch was detected and the corresponding pieces have been merged

// we take advantage of 'set_of_ver'
// which is the set of all vertices in the cloud close enough to 'ver',
// previously computed in 'progress_relocate'
// we can destroy it here, it won't be used anymore
	
// see paragraph 12.8 in the manual

{	if ( not ver.exists() )  return false;  // no touch
	if ( not ver.belongs_to ( progress_interface, tag::cell_has_low_dim, tag::not_oriented ) )
		return false;  //  because 'ver' might have been left behind in the meanwhile

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

	{ // just a block of code for hiding variables
	start_again :
	Cell next_next_seg = progress_interface.cell_in_front_of(next_ver);
	Cell next_next_ver = next_next_seg.tip();
	Cell prev_prev_seg = progress_interface.cell_behind(prev_ver);
	Cell prev_prev_ver = prev_prev_seg.base().reverse();
	for ( std::set<Cell>::iterator it = set_of_ver.begin(); it != set_of_ver.end(); it++ )
	{	Cell new_ver = *it;
		if ( new_ver == next_next_ver )
		{	Cell new_seg ( tag::segment, ver.reverse(), new_ver );
		  progress_fill_60 < manif_type >
				( next_seg, next_next_seg, new_seg.reverse(), next_ver, cloud );
			new_seg.add_to_mesh ( progress_interface );
			build_one_normal ( ver, new_ver, new_seg );  // based on previous segment
			assert ( point_120 != next_ver );
			if ( stop_point_120 == next_ver ) stop_point_120 = next_next_ver;
			next_ver = next_next_ver;  next_seg = new_seg;
			set_of_ver.erase ( it );  goto start_again;                                       }
		if ( new_ver == prev_prev_ver )
		{	Cell new_seg ( tag::segment, new_ver.reverse(), ver );
			progress_fill_60 < manif_type >
				( prev_prev_seg, prev_seg, new_seg.reverse(), prev_ver, cloud );
			new_seg.add_to_mesh ( progress_interface );
			build_one_normal ( new_ver, ver, new_seg );  // based on previous segment
			if ( point_120 == prev_ver )
			{	if ( stop_point_120 == prev_ver ) stop_point_120 = prev_prev_ver;
				point_120 = prev_prev_ver;                                        }
			prev_ver = prev_prev_ver;  prev_seg = new_seg;
			set_of_ver.erase ( it );  goto start_again;                                       }  }
	} // just a block of code for hiding variables

	// we deal with two vertices or three vertices
	std::set<Cell>::iterator it = set_of_ver.begin();
	if ( it == set_of_ver.end() )  return false;
	Cell one = *it;
	assert ( prev_ver != one );  assert ( next_ver != one );
	it++;  if ( it == set_of_ver.end() )  return false;
	Cell two = *it;
	assert ( prev_ver != two );  assert ( next_ver != two );
	it++;
	if ( it != set_of_ver.end() )  // three vertices
	{	Cell three = *it;
		it++;  assert ( it == set_of_ver.end() );
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
		if ( Manifold::working.dist_sq ( one, next_ver ) < progress_sq_long_dist )
		{	// std::cout << "touching interface S fill" << std::endl << std::flush;
			//  progress_interface.cell_in_front_of(next_ver) may have tip 'one'
			Cell ver_two = glue_two_segs_S < manif_type >
				( ver, next_ver, next_seg, one, two, one_two, cloud );
			Cell ver_three ( tag::segment, ver.reverse(), three );
			progress_fill_60 < manif_type >
				( ver_two, two_three, ver_three.reverse(), two, cloud );
			ver_three.add_to_mesh ( progress_interface );
			build_one_normal ( ver, three, ver_three );  // based on previous segment
			return true;                                                                }
		if ( Manifold::working.dist_sq ( two, prev_ver ) < progress_sq_long_dist )
		{	// std::cout << "touching interface Z fill" << std::endl << std::flush;
			//  progress_interface.cell_behind(prev_ver) may have base 'three'
			Cell two_ver = glue_two_segs_Z < manif_type >
				( prev_ver, ver, prev_seg, two, three, two_three, cloud );
			Cell one_ver ( tag::segment, one.reverse(), ver );
			progress_fill_60 < manif_type > ( one_two, two_ver, one_ver.reverse(), two, cloud );
			one_ver.add_to_mesh ( progress_interface );
			build_one_normal ( one, ver, one_ver );  // based on previous segment
			return true;                                                                         }  }
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
		// we should do something different in this situation ! e' como um beco quadrado
		if ( progress_interface.cell_in_front_of(two).tip() == prev_ver ) return false;  // no merge
		if ( progress_interface.cell_behind(one).base().reverse() == next_ver )
			return false;  // no merge
		if ( Manifold::working.dist_sq ( one, next_ver ) < progress_sq_long_dist )
		{	// std::cout << "touching interface S" << std::endl << std::flush;
			//  progress_interface.cell_in_front_of(next_ver) may have tip 'one'
			glue_two_segs_S < manif_type > ( ver, next_ver, next_seg, one, two, one_two, cloud );
			return true;                                                                         }
		if ( Manifold::working.dist_sq ( two, prev_ver ) < progress_sq_long_dist )
		{	// std::cout << "touching interface Z" << std::endl << std::flush;
			//  progress_interface.cell_behind(prev_ver) may have base 'two'
			glue_two_segs_Z < manif_type > ( prev_ver, ver, prev_seg, one, two, one_two, cloud );
			return true;                                                                         }  }

	return false;  // almost touch, no merge
	
}  // end of check_touching

//-------------------------------------------------------------------------------------------------


inline void update_info_connected_one_dim ( const Mesh msh, const Cell start, const Cell stop )
// hidden in anonymous namespace

// 'start' and 'stop' are positive vertices (may be one and the same)

{	assert ( start.dim() == 0 );
	assert ( stop.dim() == 0 );
	assert ( start.is_positive() );
	assert ( stop.is_positive() );

	Mesh::Connected::OneDim * msh_core = tag::Util::assert_cast
		< Mesh::Core*, Mesh::Connected::OneDim* > ( msh.core );
	msh_core->first_ver = start.reverse();
	msh_core->last_ver = stop;
	// now we can use an iterator

	Mesh::Iterator it = msh.iterator ( tag::over_segments, tag::require_order );
	size_t n = 0;
	for ( it.reset(); it.in_range(); it++ ) n++;
	msh_core->nb_of_segs = n;                                                   }
	
//-------------------------------------------------------------------------------------------------


// manif_type is useful for distinguishing between
// a quotient manifold and other types of manifold

template < class manif_type >
void progressive_construct           // hidden in anonymous namespace  // line 1506
(	Mesh & msh, const tag::StartAt &, Cell start,
	const tag::Towards &, std::vector<double> & normal,
	const tag::Boundary &, Mesh bdry                   )

// for two-dimensional meshes (arbitrary geometric dimension)
	
// 'start' is a vertex or segment belonging to 'bdry'
// 'normal' is a vector tangent to the working manifold, orthogonal to 'start'

{	// we don't want to change 'bdry' so we make a copy of it
	// one more reason : bdry may be Mesh::Connected::OneDim,
	// we want a Mesh::Fuzzy interface to play with
	// in the future, we will want a Mesh::STSI

	{ // just a block of code for hiding 'it', 'interface'
	Mesh interface ( tag::fuzzy, tag::of_dim, 1 );
	Mesh::Iterator it = bdry .iterator ( tag::over_cells_of_max_dim );
	for ( it .reset(); it .in_range(); it++ ) (*it) .add_to_mesh ( interface );
	progress_interface = interface;
	} // just a block of code 
	mesh_under_constr = msh;
	Cell vertex_recently_built ( tag::non_existent );
	std::set < typename manif_type::winding_cell > set_of_nearby_vertices;
	// vertices close to vertex_recently_built

	if ( start.dim() != 1 )
	{	assert ( start.dim() == 0 );
		start = progress_interface.cell_in_front_of ( start );                    }
	assert ( start.dim() == 1 );
	assert ( bdry.dim() == 1 );
	assert ( msh.dim() == 2 );
	assert ( start.belongs_to ( progress_interface, tag::same_dim, tag::oriented ) );
	
	typename manif_type::sq_dist square_dist;
	desired_len_at_point = desired_length ( start.tip() );
	sq_desired_len_at_point = desired_len_at_point * desired_len_at_point;
	typename manif_type::metric_tree cloud ( square_dist, desired_len_at_point, 6. );
	// first argument : a callable object returning the square of the distance
	// between two points, measured in the surrounding, Euclidian, space
	// second argument : distance for rank zero nodes, which is a mere hint
	// about how to initialize the tree (the tree can change a lot later)
	// third argument : ratio between distances of successive ranks
	// see paragraphs 12.10 and 12.11 in the manual

	{ // just a block of code for hiding variables
	size_t n = 0;
  Mesh::Iterator it = progress_interface.iterator ( tag::over_vertices );
	for ( it.reset(); it.in_range(); it++, n++ )
		progress_add_point < manif_type > ( *it, cloud );
	} // just a block of code for hiding variables

	{ // just a block of code for hiding variables
	double n2 = 0.;  // ensure that normal has the right norm
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	double temp = normal[i];  temp *= temp;  n2 += temp;  }
	double coef = desired_len_at_point / std::sqrt(n2);
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )  normal[i] *= coef;
	Cell start_base = start.base().reverse();
	Cell start_tip = start.tip();
	std::vector < double > vec ( progress_nb_of_coords );
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		vec[i] = x ( start_tip ) - x ( start_base );      }
	// recover tangent from hook !
	improve_normal ( start_base, start_tip, vec, normal );
	// ensures again the norm is right, projects on the tangent space,
	// ensures orthogonality with 'start' segment
	} // just a block of code for hiding variables

	assert ( start.core->hook.find(tag::normal_vector) == start.core->hook.end() );
  std::vector < double > * ff = new std::vector < double > { normal };
	start.core->hook[tag::normal_vector] = static_cast < void * > ( ff );
	// optimize !
	Cell ret = build_normals ( start );
	assert ( ret == start );

	Cell point_60 = start.tip();

	#ifndef NDEBUG
	int stopping_criterion = 0;
	// std::cout << "stopping criterion : ";  std::cin >> stopping_criterion;
	int current_name = 1;
	#endif

restart:
	
	Cell stop_point_60 = point_60, point_120 = point_60, stop_point_120 = point_60;

angles_60 :

	// starting at 'point_60', we go along this connected component of 'progress_interface'
	// we stop at 'stop_point_60', not including 'stop_point_60'
	// if point_60 == stop_point_60, we want to treat every vertex in the
	// current connected component of 'progress_interface'
	{ // just a block of code for hiding variables
	Cell prev_seg = progress_interface.cell_behind ( point_60, tag::surely_exists );
	Cell A = prev_seg.base().reverse();
	while ( true )
	{	desired_len_at_point = desired_length ( point_60 );
		sq_desired_len_at_point = desired_len_at_point * desired_len_at_point;
		progress_long_dist = desired_len_at_point * one_plus_tolerance;
		progress_sq_long_dist = progress_long_dist * progress_long_dist;
		Cell next_seg = progress_interface.cell_in_front_of ( point_60, tag::surely_exists );
		Cell B = next_seg.tip();  assert ( A != B );
		double d2 = Manifold::working.dist_sq ( A, B );
		if ( d2 < progress_sq_long_dist )  // we may want to form a triangle
		// but first we must make sure it is on the correct side
		if ( positive_orientation ( A, point_60, prev_seg, next_seg ) )
		if ( ( progress_cos_sq_60 ( A, point_60, B, prev_seg, next_seg) > 0.02 )
		  or ( d2 < sq_desired_len_at_point )                                      )
		// triangle waiting to be filled; see paragraph 12.7 in the manual
		{	Cell seg_next_to_B = progress_interface.cell_in_front_of(B);
			Cell ver_next_to_B = seg_next_to_B.tip();
			set_of_nearby_vertices.erase ( point_60 );
			if ( ver_next_to_B == A )  // this is the last triangle in this piece of progress_interface
			{	progress_fill_last_triangle < manif_type >
					( A, point_60, B, prev_seg, next_seg, seg_next_to_B, cloud );
				#ifndef NDEBUG
				// std::cout << "shrinking triangle " << ++current_name << std::endl;
				if ( current_name == stopping_criterion ) return;
				#endif
				goto search_for_start;  	                                                        }
			Cell AB ( tag::segment, A.reverse(), B );
			progress_fill_60 < manif_type > ( prev_seg, next_seg, AB.reverse(), point_60, cloud );
			AB.add_to_mesh ( progress_interface );
			build_one_normal ( A, B, AB );  // based on previous segment
			#ifndef NDEBUG
			// std::cout << "found angle around 60 deg " << ++current_name << std::endl;
			if ( current_name == stopping_criterion ) return;
			#endif
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

// check_touching :

	{ // just a block of code for hiding 'touch'
	bool touch = check_touching < manif_type >
		( vertex_recently_built, set_of_nearby_vertices,
      point_120, stop_point_120, cloud              );
	vertex_recently_built = Cell ( tag::non_existent );
	if ( touch )
	{	// if ( current_name == 296 ) return;
		assert ( point_120.belongs_to
			( progress_interface, tag::cell_has_low_dim, tag::not_oriented ) );
		#ifndef NDEBUG
		// std::cout << "touch " << ++current_name << std::endl;
		if ( current_name == stopping_criterion ) return;
		#endif
		point_60 = point_120;  stop_point_60 = point_60;  stop_point_120 = point_60;
		goto angles_60;                                                                     }
	} // just a block of code for hiding 'touch'

// look for angles around 120 deg :

	// starting at 'point_120', we go along this connected component of 'progress_interface'
	// we stop at 'stop_point_120', not including 'stop_point_120'
	// if point_120 == stop_point_120, we want to treat every vertex in the
	// current connected component of 'progress_interface'
	{  // just a block of code for hiding prev_seg and A
	Cell prev_seg = progress_interface.cell_behind ( point_120, tag::surely_exists );
	Cell A = prev_seg.base().reverse();
	while ( true )
	{	desired_len_at_point = desired_length ( point_120 );
		sq_desired_len_at_point = desired_len_at_point * desired_len_at_point;
		progress_long_dist = desired_len_at_point * one_plus_tolerance;
		progress_sq_long_dist = progress_long_dist * progress_long_dist;
		Cell next_seg = progress_interface.cell_in_front_of ( point_120, tag::surely_exists );
		Cell  B = next_seg.tip();

		if ( progress_cos_sq_120 ( A, point_120, B, prev_seg, next_seg) < 0.55 )  // 0.67
		// angle around 120 deg, we want to form two triangles; see paragraph 12.7 in the manual
		{	// we don't build a new vertex yet, we want to check for a quadrangle first
			Cell seg_prev_to_A = progress_interface.cell_behind ( A );
			Cell seg_next_to_B = progress_interface.cell_in_front_of ( B );
			Cell ver_prev_to_A = seg_prev_to_A.base().reverse();
			Cell ver_next_to_B = seg_next_to_B.tip();
			if ( ver_prev_to_A == ver_next_to_B )  // quadrangle
			{	// check orientations : correct side ?
				#ifndef NDEBUG
				// std::cout << "shrinking quadrangle" << std::endl;
				#endif
				// choose the shortest diagonal, add two triangles to the mesh under construction
				if ( Manifold::working.dist_sq ( A, B ) <
						 Manifold::working.dist_sq ( point_120, ver_next_to_B ) )
				{	Cell AB ( tag::segment, A.reverse(), B );
					progress_fill_60 < manif_type >
						( prev_seg, next_seg, AB.reverse(), point_120, cloud );
					progress_fill_60 < manif_type >
						( seg_next_to_B, seg_prev_to_A, AB, ver_prev_to_A, cloud );  }
				else
				{	Cell seg ( tag::segment, point_120.reverse(), ver_next_to_B );
					progress_fill_60 < manif_type >
						( next_seg, seg_next_to_B, seg.reverse(), B, cloud );
					progress_fill_60 < manif_type > ( seg_prev_to_A, prev_seg, seg, A, cloud );  }
				goto search_for_start;                                                             }
			Cell P ( tag::vertex );  vertex_recently_built = P;
			// now we want to place this new vertex accordingly
			assert ( prev_seg.core->hook.find(tag::normal_vector) !=
		           prev_seg.core->hook.end()                        );
		  std::vector < double > & nor_a = * static_cast < std::vector < double > * >
				( prev_seg.core->hook[tag::normal_vector] );
			assert ( next_seg.core->hook.find(tag::normal_vector) !=
		           next_seg.core->hook.end()                        );
			std::vector < double > & nor_b = * static_cast < std::vector < double > * >
				( next_seg.core->hook[tag::normal_vector] );
			std::vector < double > sum_of_nor ( progress_nb_of_coords );
			for ( size_t i = 0; i < progress_nb_of_coords; i++ )  sum_of_nor[i] = nor_a[i] + nor_b[i];
			for ( size_t i = 0; i < progress_nb_of_coords; i++ )
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
			tri1.add_to_mesh ( msh );
			tri2.add_to_mesh ( msh );
			prev_seg.remove_from_mesh ( progress_interface );
			next_seg.remove_from_mesh ( progress_interface );
			delete static_cast < std::vector < double > * > ( prev_seg.core->hook[tag::normal_vector] );
			delete static_cast < std::vector < double > * > ( next_seg.core->hook[tag::normal_vector] );
			// optimize !
			prev_seg.core->hook.erase ( tag::normal_vector );
			next_seg.core->hook.erase ( tag::normal_vector );
			AP.add_to_mesh ( progress_interface );
			PB.add_to_mesh ( progress_interface );
			cloud.remove ( static_cast < typename manif_type::metric_tree::Node * >
			               ( point_120.core->hook[tag::node_in_cloud] )            );
			point_120.core->hook.erase ( tag::node_in_cloud );  // optimize !
			build_one_normal ( A, P, AP );  // based on previous segment
			build_one_normal ( P, B, PB );  // based on previous segment
			progress_relocate < manif_type > ( P, 2, sum_of_nor, set_of_nearby_vertices, cloud );
			// find more vertices close to P and take them all into account; modifies sum_of_nor
			assert ( prev_seg.tip() == point_120 );
			if ( point_120 .is_inner_to ( mesh_under_constr ) ) msh.baricenter ( point_120 );
			progress_add_point < manif_type > ( P, cloud );
			#ifndef NDEBUG
			// std::cout << "found angle around 120 deg " << ++current_name << std::endl << std::flush;
			if ( current_name == stopping_criterion ) return;
			#endif
			if ( stop_point_120 == point_120 )  // we have all the loop to cover
				stop_point_120 = A;
			if ( stop_point_120 == B ) stop_point_120 = ver_next_to_B;
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
	// we use point_120, the desired distances have been computed there
	Cell next_seg = progress_interface.cell_in_front_of ( point_120, tag::surely_exists );

	Cell B = next_seg.tip();
	Cell P ( tag::vertex );  vertex_recently_built = P;
	assert ( next_seg.core->hook.find(tag::normal_vector) !=
	         next_seg.core->hook.end()                       );
	std::vector < double > f = * static_cast < std::vector < double > * >
		( next_seg.core->hook[tag::normal_vector] );
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		x(P) = ( x(point_120) + x(B) ) / 2. + f[i] * sqrt_of_075;  }
	Manifold::working.project(P);

	Cell AP ( tag::segment, point_120.reverse(), P );
	Cell BP ( tag::segment, B.reverse(), P );
	Cell PB = BP.reverse();
	Cell tri ( tag::triangle, next_seg, BP, AP.reverse() );
	tri.add_to_mesh ( msh );
	next_seg.remove_from_mesh ( progress_interface );
	delete static_cast < std::vector < double > * >
		( next_seg.core->hook[tag::normal_vector] );  // optimize !
	next_seg.core->hook.erase ( tag::normal_vector );
	AP.add_to_mesh ( progress_interface );
	PB.add_to_mesh ( progress_interface );
	build_one_normal ( point_120, P, AP );  // based on previous segment
	build_one_normal ( P, B, PB );  // based on previous segment
	progress_relocate < manif_type > ( P, 1, f, set_of_nearby_vertices, cloud );
	// find more vertices close to P and take them all into account; modifies f
	progress_add_point < manif_type > ( P, cloud );
	stop_point_120 = progress_interface.cell_in_front_of(B).tip();

	#ifndef NDEBUG
	// std::cout << "building brand new triangle " << ++current_name << std::endl;
	if ( current_name == stopping_criterion ) return;
	#endif

	point_60 = point_120;
	Cell BC = progress_interface.cell_in_front_of ( B, tag::surely_exists );
	stop_point_60 = BC.tip();
	stop_point_120 = stop_point_60;
	goto angles_60;
	// goto check_touching;
	} // just a block of code for hiding variables

search_for_start :  // execution only reaches this point through 'goto'

	{ // just a block of code for hiding variables
	// we look for a segment in 'progress_interface' which has a normal
	// 'progress_interface' is a one-dimensional mesh, not necessarily connected
	// so we cannot use a Mesh::Iterator - perhaps an unstructured one ?
	if ( progress_interface.number_of ( tag::segments ) == 0 ) return;
	// empty interface, meshing process ended
	#ifndef NDEBUG
	// std::cout << "search for start " << current_name << std::endl;
	#endif
	Mesh::Iterator it = progress_interface.iterator ( tag::over_segments );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell tmp = *it;
		if ( ( tmp.core->hook.find(tag::normal_vector) ) != tmp.core->hook.end() )
		{	point_60 = (*it).tip();  // any point on this connected component would do
			goto restart;            }                                                }
	assert ( false );
	} // just a block of code for hiding variables

}  // end of progressive_construct

//-----------------------------------------------------------------------------------------------


void progressive_construct ( Mesh & msh,
	const tag::StartAt &, const Cell & start,
	const tag::Towards &, std::vector < double > & tangent,
	const tag::StopAt &, const Cell & stop                 )
// hidden in anonymous namespace

// builds a one-dimensional mesh (a curve)
// orientation given by 'tangent'
	
// 'start' and 'stop' are vertices (may be one and the same)
	
{	assert ( start.dim() == 0 );
	assert ( stop.dim() == 0 );

	size_t counter = 1;
	size_t max_counter = 0;  //  std::cin >> max_counter;
	Cell A = start;
	desired_len_at_point = desired_length (A);
	sq_desired_len_at_point = desired_len_at_point * desired_len_at_point;
	while ( true )
	{	double d = Manifold::working.dist_sq ( A, stop );  // long range distance !
		double augm_length = desired_len_at_point * 1.618034,  // golden number
		       augm_len_sq = augm_length * augm_length;
		if ( d < augm_len_sq )
		{	// below we could use the Riemannian metric
			// however, the sign should be the same
			std::vector < double > e ( progress_nb_of_coords );
			double prod = 0.;
			for ( size_t i = 0; i < progress_nb_of_coords; i++ )
			{	Function x = Manifold::working .coordinates() [i];
				e[i] = x(stop) - x(A);  // recover tangent from hook !
				prod += tangent[i] * e[i];                        }
			if ( prod > 0. )
			{	Cell last ( tag::segment, A.reverse(), stop );
				last .add_to_mesh ( msh, tag::do_not_bother );
				// the meaning of tag::do_not_bother is explained
				// at the end of paragraph 11.6 in the manual
				redistribute_vertices ( msh, start, stop, 6 );    // line 710
				break;                                                                       }  }
		Cell B ( tag::vertex );
		for ( size_t i = 0; i < progress_nb_of_coords; i++ )
		{	Function x = Manifold::working.coordinates()[i];
			x(B) = x(A) + tangent[i];                         }
		if ( counter == max_counter ) return;
		Manifold::working.project ( B );
		desired_len_at_point = desired_length (B);
		sq_desired_len_at_point = desired_len_at_point * desired_len_at_point;
		for ( size_t i = 0; i < progress_nb_of_coords; i++ )
		{	Function x = Manifold::working.coordinates()[i];  // recover tangent from hook !
			tangent[i] = x(B) - x(A);                         }
		double n2 = 0.;  // Manifold::working.inner_prod ( A, tangent, tangent ); ?!!
		for ( size_t i = 0; i < progress_nb_of_coords; i++ )
		{	double temp = tangent[i];  temp *= temp;  n2 += temp;  }
		n2 = approx_sqrt ( n2, tag::around, sq_desired_len_at_point, desired_len_at_point );
		n2 = desired_len_at_point / n2;
		for ( size_t i = 0; i < progress_nb_of_coords; i++ )  tangent[i] *= n2;
		Cell AB ( tag::segment, A.reverse(), B );
		AB.add_to_mesh ( msh, tag::do_not_bother );
		// the meaning of tag::do_not_bother is explained at the end of paragraph 11.6 in the manual
		counter++;  A = B;                                                                      }
	
		update_info_connected_one_dim ( msh, start, stop );
		
} // end of  progressive_construct

//-------------------------------------------------------------------------------------------------


void progressive_construct     // hidden in anonymous namespace
( Mesh & msh, const tag::StartAt &, const Cell & start,
  const tag::Orientation &, const tag::OrientationChoice & oc )

// meshing process starts and stops at the same vertex
// if orientation not specified, we choose inherent orientation

{	assert ( start .dim() == 0 );
	assert ( start .is_positive() );
	assert ( msh .dim() == 1 );
	assert ( ( oc == tag::inherent ) or ( oc == tag::random ) or ( oc == tag::not_provided ) );

	desired_len_at_point = desired_length ( start );
	sq_desired_len_at_point = desired_len_at_point * desired_len_at_point;
	std::vector < double > best_tangent = compute_tangent_vec ( tag::at_point, start );

	progressive_construct ( msh, tag::start_at, start, tag::towards, best_tangent,
                          tag::stop_at, start                                   );

	if ( oc != tag::random )
		if ( not correctly_oriented ( msh, tag::orientation, oc ) )
			switch_orientation_direct ( msh );                                                       }
	
//-------------------------------------------------------------------------------------------------


void progressive_construct ( Mesh & msh, const tag::StartAt &, const Cell & start,   // line 1947
                             const tag::StopAt &, const Cell & stop,
                             const tag::Orientation &, const tag::OrientationChoice & oc )
// hidden in anonymous namespace

// 'start' and 'stop' are positive vertices (may be one and the same)

{	assert ( start .dim() == 0 );
	assert ( stop .dim() == 0 );
	assert ( start .is_positive() );
	assert ( stop .is_positive() );
	assert ( msh .dim() == 1 );

	desired_len_at_point = desired_length ( start );
	sq_desired_len_at_point = desired_len_at_point * desired_len_at_point;
	std::vector < double > best_tangent = compute_tangent_vec ( tag::at_point, start );

	if ( oc == tag::not_provided )
	{	std::cout << "when starting and stopping points are provided," << std::endl;
		std::cout << "maniFEM needs to know how to choose the orientation of the curve;" << std::endl;
		std::cout << "please specify either tag::orientation or tag::shortest_path" << std::endl;
		exit (1);                                                                                    }
	
	if ( oc == tag::geodesic )   // shortest path

	{	assert ( start != stop );
		
		// start walking along the manifold from 'start' in the direction of best_tangent
		// and, simultaneously, in the opposite direction, given by -best_tangent
		std::vector < double > tan1 = best_tangent, tan2 = best_tangent;
		for ( size_t i = 0; i < progress_nb_of_coords; i++ ) tan2 [i] *= -1.;
		Cell ver1 ( tag::vertex ), ver2 ( tag::vertex );
		for ( size_t i = 0; i < progress_nb_of_coords; i++ )
		{	Function x = Manifold::working .coordinates() [i];
			x ( ver1 ) = x ( start );  x ( ver2 ) = x ( start );  }
		int winner = 0;  //  will be 1 or -1
		while ( true )
		{	double augm_length = desired_length(ver1) * 1.5,
			       augm_len_sq = augm_length * augm_length;
			for ( size_t i = 0; i < progress_nb_of_coords; i++ )
			{	Function x = Manifold::working .coordinates() [i];
				x ( temporary_vertex ) = x ( ver1 ) + tan1 [i];  }
			Manifold::working .project ( temporary_vertex );
			for ( size_t i = 0; i < progress_nb_of_coords; i++ )
			{	Function x = Manifold::working .coordinates() [i];
				tan1[i] = x ( temporary_vertex ) - x ( ver1 );
				x ( ver1 ) = x ( temporary_vertex );            }
			double d = Manifold::working .dist_sq ( ver1, stop );
			if ( d < augm_len_sq )
			{	double prod = 0.;
				for ( size_t i = 0; i < progress_nb_of_coords; i++ )
				{	Function x = Manifold::working .coordinates() [i];
					prod += tan1[i] * ( x (stop) - x (ver1) );         }
				if ( prod > 0. )  { winner = 1;  break;  }             }
			for ( size_t i = 0; i < progress_nb_of_coords; i++ )
			{	Function x = Manifold::working .coordinates() [i];
				x ( temporary_vertex ) = x ( ver2 ) + tan2 [i];  }
			Manifold::working .project ( temporary_vertex );
			for ( size_t i = 0; i < progress_nb_of_coords; i++ )
			{	Function x = Manifold::working .coordinates() [i];
				tan2[i] = x ( temporary_vertex ) - x ( ver2 );
				x ( ver2 ) = x ( temporary_vertex );            }
			d = Manifold::working .dist_sq ( ver2, stop );
			if ( d < augm_len_sq )
			{	double prod = 0.;
				for ( size_t i = 0; i < progress_nb_of_coords; i++ )
				{	Function x = Manifold::working .coordinates() [i];
					prod += tan2 [i] * ( x (stop) - x (ver2) );          }
				if ( prod > 0. )  { winner = -1;  break;  }           }
		}  // end of  while true

		assert ( ( winner == 1 ) or ( winner == -1 ) );
		for ( size_t i = 0; i < progress_nb_of_coords; i++ ) best_tangent [i] *= winner;
		progressive_construct ( msh, tag::start_at, start, tag::towards, best_tangent,
		                        tag::stop_at, stop                                    );
		return;                                                                            }

	if ( ( oc == tag::inherent ) or ( oc == tag::random ) )

	{	if ( start == stop )
		{	progressive_construct ( msh, tag::start_at, start, tag::orientation, oc );
			return;                                                                   }

		progressive_construct ( msh, tag::start_at, start, tag::towards, best_tangent,
 		                        tag::stop_at, stop );

		if ( oc == tag::random ) return;

		for ( size_t i = 0; i < progress_nb_of_coords; i++ )  best_tangent[i] *= -1.;
		// the number of segments does not count, and we don't know it yet
		Mesh msh2 ( tag::whose_core_is,
		    new Mesh::Connected::OneDim ( tag::with, 1, tag::segments, tag::one_dummy_wrapper ),
	  	  tag::freshly_created, tag::is_positive                                              );
		// the number of segments does not count, and we don't know it yet
		// we compute it after the mesh is built, by counting segments
		// but we count segments using an iterator, and the iterator won't work
		// if this->msh->nb_of_segs == 0, so we set nb_of_segs to 1 (dirty trick)
		// see Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder::reset
		// in iterator.cpp
		progressive_construct ( msh2, tag::start_at, start, tag::towards, best_tangent,
		                        tag::stop_at, stop                                     );

		switch_orientation_direct ( msh2 );
		update_info_connected_one_dim ( msh2, stop, start );
		Mesh whole ( tag::join, msh, msh2 );

		if ( not correctly_oriented ( whole, tag::orientation, oc ) )
		{	switch_orientation_direct ( msh2 );  msh = msh2;  }

		return;                                                                                  }

	assert ( oc == tag::intrinsic );
	std::cout << "intrinsic orientation does not make sense for one-dimensional manifolds"
	          << std::endl;
	exit (1);

}  // end of  progressive_construct

//-------------------------------------------------------------------------------------------------


inline void progressive_construct     // hidden in anonymous namespace  // line 2075
( Mesh & msh, const tag::StartWithNonExistentMesh &,
  const tag::StartAt, const Cell & start,
  const tag::Orientation &, const tag::OrientationChoice & oc )

// last argument may be  tag::random  or  tag::inherent	 or  tag::not_provided

{	assert ( not msh .exists() );

	// call to 'compute_tangent_vec' does not depend on the dimension of the mesh
	desired_len_at_point = desired_length ( start );
	sq_desired_len_at_point = desired_len_at_point * desired_len_at_point;
	std::vector < double > tangent = compute_tangent_vec ( tag::at_point, start );

	// now we branch, depending on the dimension
	if ( get_topological_dim() == 1 )
	{	msh .core = new Mesh::Connected::OneDim
			( tag::with, 1, tag::segments, tag::one_dummy_wrapper );
		// the number of segments does not count, and we don't know it yet
		// we compute it after the mesh is built, by counting segments
		// but we count segments using an iterator, and the iterator won't work
		// if this->msh->nb_of_segs == 0, so we set nb_of_segs to 1 (dirty trick)
		// see Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder::reset
		// in iterator.cpp

		if ( oc == tag::intrinsic )
		{	std::cout << "intrinsic orientation does not apply to one-dimensional meshes"
			          << std::endl;
			exit (1);                                                                    }

		progressive_construct ( msh, tag::start_at, start, tag::towards, tangent,
		                        tag::stop_at, start                              );

		if ( oc == tag::random ) return;

		if ( progress_nb_of_coords > 2 )
		{	std::cout << "co-dimension 2 or larger, "
			          << "please specify tag::orientation, tag::random" << std::endl;
			exit (1);                                                                 }

		assert ( progress_nb_of_coords == 2 );
		assert ( ( oc == tag::inherent ) or ( oc == tag::not_provided ) );
		// here we interpret "not provided" as "inherent"
		if ( not correctly_oriented ( msh, tag::orientation, oc ) )
			switch_orientation_direct ( msh );
		assert ( correctly_oriented ( msh, tag::orientation, oc ) );
		return;                                                                            }

	// else :  top dim > 1
	{	assert ( get_topological_dim() == 2 );  // no 3D for now
		assert ( progress_nb_of_coords == 3 );

		if ( oc == tag::intrinsic )
		{	std::cout << "intrinsic orientation makes no sense here, "
			          << "did you mean inherent ?" << std::endl;
			exit (1);                                                  }

		msh.core = new Mesh::Fuzzy ( tag::of_dim, 3, tag::minus_one, tag::one_dummy_wrapper );
		Cell B ( tag::vertex );
		for ( size_t i = 0; i < progress_nb_of_coords; i++ )
		{	Function x = Manifold::working .coordinates() [i];
			x ( B ) = x ( start ) + tangent [i];               }
		Cell AB ( tag::segment, start .reverse(), B );
		std::vector < double > normal =
			compute_tangent_vec ( tag::at_point, start, tag::orthogonal_to, tangent );
		Cell C ( tag::vertex );
		for ( size_t i = 0; i < progress_nb_of_coords; i++ )
		{	Function x = Manifold::working.coordinates()[i];
			x ( C ) = x ( start ) + 0.5 * tangent [i] - sqrt_of_075 * normal [i];  }
		Cell BC ( tag::segment, B .reverse(), C );
		Cell CA ( tag::segment, C .reverse(), start );
		Cell tri ( tag::triangle, AB, BC, CA );
		tri.add_to_mesh ( msh );
		
		// Mesh interf ( tag::of_dimension_one );
		Mesh interf ( tag::whose_core_is,
		    new Mesh::Connected::OneDim ( tag::with, 3, tag::segments, tag::one_dummy_wrapper ),
		    tag::freshly_created, tag::is_positive                                              );
		AB .reverse() .add_to_mesh ( interf, tag::do_not_bother );
		BC .reverse() .add_to_mesh ( interf, tag::do_not_bother );
		CA .reverse() .add_to_mesh ( interf, tag::do_not_bother );
		// the meaning of tag::do_not_bother is explained at the end of paragraph 11.6 in the manual
		update_info_connected_one_dim ( interf, B, B );
		progressive_construct < Manifold::Type::Euclidian >
			( msh, tag::start_at, AB.reverse(), tag::towards, normal, tag::boundary, interf );

		if ( oc == tag::random ) return;
		assert ( ( oc == tag::inherent ) or ( oc == tag::not_provided ) );
		// here we interpret "not provided" as "inherent"
		if ( not correctly_oriented ( msh, tag::orientation, oc ) )
			switch_orientation_of_each_cell ( msh );
		assert ( correctly_oriented ( msh, tag::orientation, oc ) );
		return;                                                                                 }

	assert ( false );

}  // end of  progressive_construct

//-------------------------------------------------------------------------------------------------
	

inline void progressive_construct         // hidden in anonymous namespace  // line 2175
( Mesh & msh, const tag::StartAt &, const Cell & start,
  const tag::Boundary &, Mesh interface,
  const tag::Orientation &, const tag::OrientationChoice & oc )

// for two-dimensional meshes in RR^2 (intrinsic orientation)
//   or in a 2D submanifold of RR^3 (random or inherent orientation)
	
// 'start' is a segment belonging to 'interface'
// no normal vector provided, we need to build our own

{	if ( progress_nb_of_coords == 2 )
		
	{	// domain in the plane RR^2, intrinsic orientation
		
		if ( oc == tag::random )
		{	std::cout << "it is unsafe to try to mesh a plane region with random orientation;"
			          << std::endl;
			std::cout << "the region may be unbounded" << std::endl;
			exit (1);                                                                          }
		if ( oc == tag::inherent )
		{	std::cout << "inherent orientation does not apply to a plane region;" << std::endl;
			std::cout << "did you mean 'intrinsic' ?" << std::endl;
			exit (1);                                                                           }
		assert ( ( oc == tag::intrinsic ) or ( oc == tag::not_provided ) );
		// here we interpret "not provided" as "intrinsic"
		
		Cell A = start .base() .reverse();
		Cell B = start .tip();
		Function x = Manifold::working.coordinates() [0];
		Function y = Manifold::working.coordinates() [1];
		std::vector < double > nor { y(A) - y(B), x(B) - x(A) };   // rotate with 90 deg

		progressive_construct < Manifold::Type::Euclidian >
			( msh, tag::start_at, start, tag::towards, nor, tag::boundary, interface );
		return;                                                                        }
	
	assert ( progress_nb_of_coords == 3 );
	assert ( msh .dim() == 2 );
	assert ( start .core->belongs_to ( interface.core, tag::same_dim, tag::oriented ) );
	
	// compute a normal vector, on an arbitrary side of 'start'
	Cell A = start .base() .reverse();
	Cell B = start .tip();
	Function x = Manifold::working .coordinates() [0];
	Function y = Manifold::working .coordinates() [1];
	Function z = Manifold::working .coordinates() [2];
	std::vector < double > tan { x(B) - x(A), y(B) - y(A), z(B) - z(A) };
	desired_len_at_point = desired_length ( A );
	sq_desired_len_at_point = desired_len_at_point * desired_len_at_point;
	std::vector < double > nor =
		compute_tangent_vec ( tag::at_point, A, tag::orthogonal_to, tan );

	progressive_construct < Manifold::Type::Euclidian >
		( msh, tag::start_at, start, tag::towards, nor, tag::boundary, interface );

	if ( oc == tag::random ) return;
	assert ( ( oc == tag::inherent ) or ( oc == tag::not_provided ) );
	// here we interpret "not provided" as "inherent"

	Mesh interf_rev ( tag::whose_core_is,   // of dimesion one
		new Mesh::Fuzzy ( tag::of_dim, 2, tag::minus_one, tag::one_dummy_wrapper ),
		tag::freshly_created, tag::is_positive                                      );
	Mesh::Iterator it = interface .iterator ( tag::over_segments );
	for ( it .reset(); it .in_range(); it++ )
		(*it) .reverse() .add_to_mesh ( interf_rev );
	assert ( start .reverse() .core->belongs_to
			( interf_rev .core, tag::same_dim, tag::oriented ) );

	// std::cout << "wait a minute ..." << std::endl;
	// build the mesh on the other side of 'interface'
	for ( size_t i = 0; i < 3; i++ )  nor [i] *= -1.;
	Mesh msh2 ( tag::whose_core_is,   // of dimesion two
		new Mesh::Fuzzy ( tag::of_dim, 3, tag::minus_one, tag::one_dummy_wrapper ),
		tag::freshly_created, tag::is_positive                                      );
	progressive_construct < Manifold::Type::Euclidian >
		( msh2, tag::start_at, start .reverse(), tag::towards, nor, tag::boundary, interf_rev );	
	// join everything to get a mesh on the entire manifold
	Mesh whole ( tag::join, msh, msh2 );
	
	if ( not correctly_oriented ( whole, tag::orientation, oc ) )
	{	switch_orientation_of_each_cell ( msh2 );  msh = msh2;  }

}  // end of  progressive_construct

//-------------------------------------------------------------------------------------------------
	
}  // end of anonymous namespace

//-------------------------------------------------------------------------------------------------


Mesh::Mesh ( const tag::Frontal &, const tag::DesiredLength &, const Function & length )

// since no boundary is provided, we assume the user wants the entire working manifold
// we assume the working manifold is compact
// we also assume the user wants the inherent orientation
// a special case : mesh a quotient manifold entirely
	
:	Mesh ( tag::non_existent )
// we don't know yet the dimension, so we postpone the constructor

{	temporary_vertex = Cell ( tag::vertex );
	progress_nb_of_coords = Manifold::working .coordinates() .nb_of_components();
	desired_length = length;

	// call to 'search_start_ver' does not depend on the dimension of the mesh
	Cell start = search_start_ver ( );

	progressive_construct ( *this, tag::start_with_non_existent_mesh,   // line 2075
	                        tag::start_at, start, tag::orientation, tag::not_provided );  }
	// last argument is equivalent, in this case, to tag::inherent

//-------------------------------------------------------------------------------------------------


Mesh::Mesh ( const tag::Frontal &, const tag::EntireManifold &, Manifold manif,
             const tag::DesiredLength &, const Function &  length                  )

// we assume 'manif' is a compact manifold
// we also assume the user wants the inherent orientation
// a special case : mesh a quotient manifold entirely
	
:	Mesh ( tag::non_existent )
// we don't know yet the dimension, so we postpone the constructor

{	Manifold tmp_manif = Manifold::working;
	Manifold::working = manif;
	
	temporary_vertex = Cell ( tag::vertex );
	progress_nb_of_coords = Manifold::working .coordinates() .nb_of_components();
	desired_length = length;

	// call to 'search_start_ver' does not depend on the dimension of the mesh
	Cell start = search_start_ver ( );

	progressive_construct ( *this, tag::start_with_non_existent_mesh,   // line 2075
	                        tag::start_at, start, tag::orientation, tag::not_provided );
	// last argument is equivalent, in this case, to tag::inherent

	Manifold::working = tmp_manif;                                                   }

//-------------------------------------------------------------------------------------------------


Mesh::Mesh ( const tag::Frontal &, const tag::DesiredLength &, const Function & length,
             const tag::Orientation &, const tag::OrientationChoice & oc                   )

// since no boundary is provided, we assume the user wants the entire working manifold
// we assume the working manifold is compact
	
:	Mesh ( tag::non_existent )
// we don't know yet the dimension, so we postpone the constructor

{	temporary_vertex = Cell ( tag::vertex );
	progress_nb_of_coords = Manifold::working .coordinates().nb_of_components();
	desired_length = length;

	// call to 'search_start_ver' does not depend on the dimension of the mesh
	Cell start = search_start_ver ( );

	progressive_construct ( *this, tag::start_with_non_existent_mesh,   // line 2075
                          tag::start_at, start, tag::orientation, oc );        }

//-------------------------------------------------------------------------------------------------


Mesh::Mesh ( const tag::Frontal &, const tag::EntireManifold &, Manifold manif,
             const tag::DesiredLength &, const Function & length,
             const tag::Orientation &, const tag::OrientationChoice & oc           )

// we assume 'manif' is a compact manifold
	
:	Mesh ( tag::non_existent )
// we don't know yet the dimension, so we postpone the constructor

{	Manifold tmp_manif = Manifold::working;
	Manifold::working = manif;
	
	temporary_vertex = Cell ( tag::vertex );
	progress_nb_of_coords = Manifold::working .coordinates() .nb_of_components();
	desired_length = length;

	// call to 'search_start_ver' does not depend on the dimension of the mesh
	Cell start = search_start_ver ( );

	progressive_construct ( *this, tag::start_with_non_existent_mesh,   // line 2075
                          tag::start_at, start, tag::orientation, oc );

	Manifold::working = tmp_manif;                                                 }

//-------------------------------------------------------------------------------------------------


Mesh::Mesh ( const tag::Frontal &, const tag::Boundary &, Mesh interface,
             const tag::DesiredLength &, const Function & length              )

// for now, only works for two-dimensional meshes (either in RR2 or in RR3)
// should be adapted for three-dimensional meshes

:	Mesh ( tag::whose_core_is,
	       new Mesh::Fuzzy ( tag::of_dim, 3, tag::minus_one, tag::one_dummy_wrapper ),
	       tag::freshly_created, tag::is_positive                                     )

{	temporary_vertex = Cell ( tag::vertex );
	progress_nb_of_coords = Manifold::working .coordinates() .nb_of_components();
	desired_length = length;

	// we search for a starting point
	Mesh::Iterator it = interface.iterator ( tag::over_segments );
	it.reset();  assert ( it.in_range() );

	progressive_construct ( *this, tag::start_at, *it,         // line 2175
                          tag::boundary, interface, tag::orientation, tag::not_provided );  }
	
//-------------------------------------------------------------------------------------------------


Mesh::Mesh ( const tag::Frontal &, const tag::Boundary &, Mesh interface,
             const tag::DesiredLength &, const Function & length,
             const tag::Orientation &, const tag::OrientationChoice & oc     )

// for now, only works for two-dimensional meshes in RR2
// should be adapted for three-dimensional meshes

// 'oc' may be  inherent  intrinsic  random
	
:	Mesh ( tag::whose_core_is,
	       new Mesh::Fuzzy ( tag::of_dim, 3, tag::minus_one, tag::one_dummy_wrapper ),
	       tag::freshly_created, tag::is_positive                                     )

{	temporary_vertex = Cell ( tag::vertex );
	progress_nb_of_coords = Manifold::working .coordinates() .nb_of_components();
	desired_length = length;

	// we search for a starting point
	Mesh::Iterator it = interface.iterator ( tag::over_segments );
	it.reset();  assert ( it.in_range() );

	progressive_construct ( *this, tag::start_at, *it,         // line 2175
	                        tag::boundary, interface, tag::orientation, oc );          }
	
//-------------------------------------------------------------------------------------------------


Mesh::Mesh ( const tag::Frontal &, const tag::Boundary &, Mesh interface,
             const tag::StartAt &, const Cell & start,
             const tag::Towards &, std::vector<double> normal,             
             const tag::DesiredLength &, const Function & length             )
	
// 'start' is a vertex or segment belonging to 'interface'

:	Mesh ( tag::whose_core_is,
	       new Mesh::Fuzzy ( tag::of_dim, 3, tag::minus_one, tag::one_dummy_wrapper ),
	       tag::freshly_created, tag::is_positive                                     )

{	temporary_vertex = Cell ( tag::vertex );
	progress_nb_of_coords = Manifold::working .coordinates() .nb_of_components();
	desired_length = length;

	progressive_construct < Manifold::Type::Euclidian >   // line 1506
		( *this, tag::start_at, start, tag::towards, normal, tag::boundary, interface );  }
	
//-------------------------------------------------------------------------------------------------


Mesh::Mesh ( const tag::Frontal &, const tag::StartAt &, const Cell & start,
             const tag::Towards &, std::vector<double> tangent,
             const tag::StopAt &, const Cell & stop,
             const tag::DesiredLength &, const Function & length                )

// 'start' and 'stop' are positive vertices, may be one and the same

:	Mesh ( tag::whose_core_is,
	       new Mesh::Connected::OneDim ( tag::with, 1, tag::segments, tag::one_dummy_wrapper ),
	       tag::freshly_created, tag::is_positive                                              )
// the number of segments does not count, and we don't know it yet
// we compute it after the mesh is built, by counting segments
// but we count segments using an iterator, and the iterator won't work
// if this->msh->nb_of_segs == 0, so we set nb_of_segs to 1 (dirty trick)
// see Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder::reset
// in iterator.cpp

{	temporary_vertex = Cell ( tag::vertex );
	progress_nb_of_coords = Manifold::working .coordinates() .nb_of_components();
	desired_length = length;

	desired_len_at_point = desired_length ( start );
	sq_desired_len_at_point = desired_len_at_point * desired_len_at_point;

	double n2 = 0.;  // ensure that normal has the right norm
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	double temp = tangent [i];  temp *= temp;  n2 += temp;  }
	double coef = desired_len_at_point / std::sqrt(n2);
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )  tangent [i] *= coef;

	improve_tangent ( start, tangent );
	// ensures again the norm is right, then projects on the tangent space

	progressive_construct ( *this, tag::start_at, start, tag::towards, tangent,
	                        tag::stop_at, stop                                 );  }

//-------------------------------------------------------------------------------------------------


Mesh::Mesh ( const tag::Frontal &, const tag::StartAt &, const Cell & start,
             const tag::Towards &, std::vector<double> tangent,
             const tag::DesiredLength &, const Function & length                )

: Mesh ( tag::whose_core_is,
         new Mesh::Connected::OneDim ( tag::with, 1, tag::segments, tag::one_dummy_wrapper ),
         tag::freshly_created, tag::is_positive                                               )
// the number of segments does not count, and we don't know it yet
// we compute it after the mesh is built, by counting segments
// but we count segments using an iterator, and the iterator won't work
// if this->msh->nb_of_segs == 0, so we set nb_of_segs to 1 (dirty trick)
// see Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder::reset
// in iterator.cpp

{	temporary_vertex = Cell ( tag::vertex );
	progress_nb_of_coords = Manifold::working .coordinates() .nb_of_components();
	desired_length = length;

	desired_len_at_point = desired_length ( start );
	sq_desired_len_at_point = desired_len_at_point * desired_len_at_point;

	double n2 = 0.;  // ensure that normal has the right norm
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	double temp = tangent[i];  temp *= temp;  n2 += temp;  }
	double coef = desired_len_at_point / std::sqrt(n2);
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )  tangent[i] *= coef;
	desired_length = length;

	improve_tangent ( start, tangent );
	// ensures again the norm is right, then projects on the tangent space

	progressive_construct ( *this, tag::start_at, start, tag::towards, tangent,
	                        tag::stop_at, start                                );  }

//-------------------------------------------------------------------------------------------------


Mesh::Mesh ( const tag::Frontal &, const tag::StartAt &, const Cell & start,
             const tag::DesiredLength &, const Function & length                )

:	Mesh ( tag::non_existent )
// we don't know yet the dimension, so we postpone the constructor

{	temporary_vertex = Cell ( tag::vertex );
	progress_nb_of_coords = Manifold::working .coordinates() .nb_of_components();
	desired_length = length;

	if ( get_topological_dim() == 1 )
	{	this->core = new Mesh::Connected::OneDim
			( tag::with, 1, tag::segments, tag::one_dummy_wrapper );
		// the number of segments does not count, and we don't know it yet
		// we compute it after the mesh is built, by counting segments
		// but we count segments using an iterator, and the iterator won't work
		// if this->msh->nb_of_segs == 0, so we set nb_of_segs to 1 (dirty trick)
		// see Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder::reset
		// in iterator.cpp
		// no stopping point, we assume the user wants to mesh the entire working manifold
		// we assume the working manifold is compact
		if ( progress_nb_of_coords == 2 )
			progressive_construct ( *this, tag::start_at, start, tag::orientation, tag::inherent );
		else  // random orientation
			progressive_construct ( *this, tag::start_at, start,
			                        tag::stop_at, start, tag::orientation, tag::random );          }
	else
	{	assert ( get_topological_dim() == 2 );  // no 3D meshing for now
		progressive_construct ( *this, tag::start_with_non_existent_mesh,
		                        tag::start_at, start, tag::orientation, tag::not_provided );  }    }

//-------------------------------------------------------------------------------------------------


Mesh::Mesh ( const tag::Frontal &, const tag::StartAt &, const Cell & start,
             const tag::DesiredLength &, const Function & length,
             const tag::Orientation &, const tag::OrientationChoice & oc        )

:	Mesh ( tag::non_existent )
// we don't know yet the dimension, so we postpone the constructor
	
{	temporary_vertex = Cell ( tag::vertex );
	progress_nb_of_coords = Manifold::working .coordinates() .nb_of_components();
	desired_length = length;

	if ( get_topological_dim() == 1 )
	{	this->core = new Mesh::Connected::OneDim
			( tag::with, 1, tag::segments, tag::one_dummy_wrapper );
		// the number of segments does not count, and we don't know it yet
		// we compute it after the mesh is built, by counting segments
		// but we count segments using an iterator, and the iterator won't work
		// if this->msh->nb_of_segs == 0, so we set nb_of_segs to 1 (dirty trick)
		// see Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder::reset
		// in iterator.cpp

		// this mesh is inconsistent, its destructor may produce an error

		progressive_construct ( *this, tag::start_at, start,  // line 1947
		                        tag::stop_at, start, tag::orientation, oc );   }
	else
	{	assert ( get_topological_dim() == 2 );  // no 3D meshing for now
		progressive_construct ( *this, tag::start_with_non_existent_mesh,  // line 2075
	                          tag::start_at, start, tag::orientation, oc );  }        }

//-------------------------------------------------------------------------------------------------


Mesh::Mesh ( const tag::Frontal &, const tag::StartAt &, const Cell & start,
             const tag::StopAt &, const Cell & stop,
             const tag::DesiredLength &, const Function & length                )

// build a one-dimensional mesh, choose shortest path
	
:	Mesh ( tag::whose_core_is,
	       new Mesh::Connected::OneDim ( tag::with, 1, tag::segments, tag::one_dummy_wrapper ),
	       tag::freshly_created, tag::is_positive                                              )
// the number of segments does not count, and we don't know it yet
// we compute it after the mesh is built, by counting segments
// but we count segments using an iterator, and the iterator won't work
// if this->msh->nb_of_segs == 0, so we set nb_of_segs to 1 (dirty trick)
// see Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder::reset
// in iterator.cpp

// this mesh is inconsistent, its destructor may produce an error

{	temporary_vertex = Cell ( tag::vertex );
	progress_nb_of_coords = Manifold::working .coordinates() .nb_of_components();
	desired_length = length;

	progressive_construct ( *this, tag::start_at, start, tag::stop_at, stop,   // line 1947
	                        tag::orientation, tag::not_provided             );    }

//-------------------------------------------------------------------------------------------------


Mesh::Mesh ( const tag::Frontal &, const tag::StartAt &, const Cell & start,
             const tag::StopAt &, const Cell & stop,
             const tag::DesiredLength &, const Function & length,
             const tag::Orientation &, const tag::OrientationChoice & oc        )

: Mesh ( tag::whose_core_is,
         new Mesh::Connected::OneDim ( tag::with, 1, tag::segments, tag::one_dummy_wrapper ),
         tag::freshly_created, tag::is_positive                                              )
// the number of segments does not count, and we don't know it yet
// we compute it after the mesh is built, by counting segments
// but we count segments using an iterator, and the iterator won't work
// if this->msh->nb_of_segs == 0, so we set nb_of_segs to 1 (dirty trick)
// see Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder::reset
// in iterator.cpp

// this mesh is inconsistent, its destructor may produce an error

{	temporary_vertex = Cell ( tag::vertex );
	progress_nb_of_coords = Manifold::working .coordinates() .nb_of_components();
	desired_length = length;
	
	progressive_construct ( *this, tag::start_at, start, tag::stop_at, stop,   // line 1947
	                        tag::orientation, oc                            );   }

//-------------------------------------------------------------------------------------------------


Mesh::Mesh ( const tag::Frontal &, const tag::StartAt &, const Cell & start,
             const tag::StopAt &, const Cell & stop,
             const tag::DesiredLength &, const Function & length,
             const tag::ShortestPath &                                       )

: Mesh ( tag::whose_core_is,
         new Mesh::Connected::OneDim ( tag::with, 1, tag::segments, tag::one_dummy_wrapper ),
         tag::freshly_created, tag::is_positive                                              )
// the number of segments does not count, and we don't know it yet
// we compute it after the mesh is built, by counting segments
// but we count segments using an iterator, and the iterator won't work
// if this->msh->nb_of_segs == 0, so we set nb_of_segs to 1 (dirty trick)
// see Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder::reset
// in iterator.cpp

// this mesh is inconsistent, its destructor may produce an error

{	temporary_vertex = Cell ( tag::vertex );
	progress_nb_of_coords = Manifold::working .coordinates() .nb_of_components();
	desired_length = length;
	
	progressive_construct ( *this, tag::start_at, start, tag::stop_at, stop,   // line 1947
	                        tag::orientation, tag::geodesic                 );   }


//-------------------------------------------------------------------------------------------------


void print_winding ( Manifold::Action a )

{	Manifold::Quotient * manif_q = dynamic_cast
		< Manifold::Quotient* > ( Manifold::working .core );
	assert ( manif_q );
	Function xy = manif_q->base_space .coordinates();
	Function x = xy[0], y = xy[1];
	size_t n = manif_q->actions.size();
	assert ( n == manif_q->winding_nbs.size() );
	std::cout << "(";
	for ( size_t i = 0; i < n; i++ )
	{	Function::ActionGenerator & g = manif_q->actions[i];
		std::map<Function::ActionGenerator,short int>::const_iterator itt = a.index_map.find ( g );
		if ( itt == a.index_map.end() )
		{	std::cout << "0,"; continue;  }
		short int exp = itt->second;
		assert ( exp != 0 );
		std::cout << exp << ",";                                                            }
	std::cout << ")" << std::endl;                                                           }


void test_sq_dist ()

{	Manifold space = Manifold::working;
	assert ( space.exists() );  // we use the current (quotient) manifold
	Manifold::Quotient * manif_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space.core );
	assert ( manif_q );
	Function coords_q = space.coordinates();
	Manifold mani_Eu = manif_q->base_space;  // underlying Euclidian manifold
	Function coords_Eu = mani_Eu.coordinates();
	assert ( coords_Eu.nb_of_components() == 2 );
	Function x = coords_Eu[0],  y = coords_Eu[1];

	Cell A ( tag::vertex );  x(A) = 0.; y(A) = 0.;
	Cell B ( tag::vertex );  x(B) = 0.1;  y(B) = 0.7;

	Manifold::Type::Quotient::sq_dist sd;

	std::pair < Cell, Manifold::Action > AA { A, 0 };
	std::pair < Cell, Manifold::Action > BB { B, 0 };
	double dist2 = sd ( AA, BB );
	std::cout << dist2 << std::endl;
	print_winding ( BB.second );
}



