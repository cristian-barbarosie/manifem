
// frontal.cpp 2022.07.02

// almost total remake

//   This file is part of maniFEM, a C++ library for meshes and finite elements on manifolds.

//   Copyright 2019 -- 2022 Cristian Barbarosie cristian.barbarosie@gmail.com

//   http://manifem.rd.ciencias.ulisboa.pt/
//   https://github.com/cristian-barbarosie/manifem

//   ManiFEM is free software: you can redistribute it and/or modify it
//   under the terms of the GNU Lesser General Public License as published
//   by the Free Software Foundation, either version 3 of the License
//   or (at your option) any later version.

//   ManiFEM is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty
//   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//   See the GNU Lesser General Public License for more details.

//   You should have received a copy of the GNU Lesser General Public License
//   along with maniFEM.  If not, see <https://www.gnu.org/licenses/>.


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

//------------------------------------------------------------------------------------------------------//


// we want to deal with Euclidian manifolds and implicit submanifolds, on one side
// quotient manifolds on the other side
// not easy to reach both goals while keeping some of the common code ... common

// for quotient manifolds, we use a dirty trick
// we provide as Point not just a Cell but a Cell together with a winding number
// each call to SqDist will set the "winning" winding of the second Point B
// relatively to the first one A
// (the winding number of a future segment AB where the minimum distance is achieved)

// perhaps change the template of MetricTree like
// template < typename Point, typename SqDist, typename RichPoint = Point >
// where RichPoint contains a winding number while Point does not
// we will have to provide two versions of SqDist :
//    SqDist ( const Point & A, const Point & B )
//    SqDist ( const Point & A, RichPoint & B )
// the latter sets the winding of B to that winding which achieves the minimum distance


// there are two difficulties here

// first, we must distinguish between a "usual" manifold and a quotient one
// a "usual" manifold may be e.g. Euclidian or implicit
// the 'manif_type' will help us make this distinction
// it allows to have two types of cells (simple and rich, the latter having winding)

// second, we must distinguish between different types of metric :
// trivial, isotropic (constant zoom, variable zoom), anisotropic, Rayleigh
// virtual methods of the metric itself will allow us to treat differently these situations


// for computing this winding number,
// implement a discrete version of the steepest descent method
// rather than a blind spiral search as in draw_ps with tag::windng

//------------------------------------------------------------------------------------------------------//


// global variables and functions for this file, not visible for other object files
namespace {  // anonymous namespace, mimics static linkage

Cell temporary_vertex ( tag::non_existent );

size_t frontal_nb_of_coords;  // dimension of the surrounding Euclidian space
// also known as "geometric dimension"


//------------------------------------------------------------------------------------------------------//

inline void update_info_connected_one_dim ( const Mesh & msh, const Cell & start, const Cell & stop )
// hidden in anonymous namespace

// 'start' and 'stop' are positive vertices (may be one and the same)

// define this as a method of class Mesh !!

{	assert ( start .dim() == 0 );
	assert ( stop  .dim() == 0 );
	assert ( start .is_positive() );
	assert ( stop  .is_positive() );

	Mesh::Connected::OneDim * msh_core = tag::Util::assert_cast
		< Mesh::Core*, Mesh::Connected::OneDim* > ( msh .core );
	msh_core->first_ver = start .reverse();
	msh_core->last_ver = stop;
	// now we can use an iterator

	Mesh::Iterator it = msh .iterator ( tag::over_segments, tag::require_order );
	size_t n = 0;
	for ( it .reset(); it .in_range(); it++ ) n++;
	msh_core->nb_of_segs = n;                                                    }
	
//-------------------------------------------------------------------------------------------------


class ManifoldNoWinding  // hidden in anonymous namespace

// Euclidian manifold, implicit submanifold, parametric

{	public :

	class sq_dist;
	typedef Cell winding_cell;
	typedef std::vector < double > winding_vector;
	typedef MetricTree < Cell, sq_dist > metric_tree;

	inline static const Cell & remove_winding ( const Cell & P )  {  return P;  }
	inline static const std::vector < double > & remove_winding ( const std::vector < double > & v )
	{	return v;  }

	inline static std::vector < double > get_vector ( const Cell & A, const Cell & B );

	inline static Cell build_seg ( const Cell & A, const Cell & B, const winding_vector & e )
	{	return Cell ( tag::segment, A, B );  }  // A is negative
			
	inline static void redistribute_vertices
	( const Mesh & msh, const Cell & start, const Cell & stop, size_t n, const std::vector < double > & );
		
};  // end of  class ManifoldNoWinding
	
//------------------------------------------------------------------------------------------------------//


inline std::vector < double > ManifoldNoWinding::get_vector ( const Cell & A, const Cell & B )
// static
	
{	const Function & coord = Manifold::working .coordinates();
	const size_t n = coord .nb_of_components();
	std::vector < double > res (n);
	for ( size_t i = 0; i < n; i++ )
	{	const Function & x = coord [i];
		res[i] = x(B) - x(A);           }
	return res;                                                  }
			
//------------------------------------------------------------------------------------------------------//

	
inline void ManifoldNoWinding::redistribute_vertices  // static
( const Mesh & msh, const Cell & start, const Cell & stop, size_t n, const std::vector < double > & )

// just make some baricenters; last argument not used
	
{	Cell A = stop;
	// just in case n is too large, or the curve is too short, we look for 'start'
	// the statement below is important for closed loops, where start == stop
	A = msh .cell_behind (A) .base() .reverse();
	for ( size_t i = 2; i < n; i++ )
	{	if ( A == start )  {  n = i;  break;  }
		A = msh .cell_behind ( A, tag::surely_exists ) .base() .reverse();  }
	assert ( n > 2 );
	Cell C = stop;
	Cell B = msh .cell_behind ( A, tag::surely_exists ) .base() .reverse();
	A = msh .cell_behind ( B, tag::surely_exists ) .base() .reverse();
	for ( size_t i = 2; i < n; i++ )
	{	Manifold::working .interpolate ( B, 0.3, A, 0.4, B, 0.3, C );
		// the above method is not compatible with a non-uniform metric
		// so it will produce imperfect (but hopefully acceptable) results
		assert ( A != start );
		if ( A == start ) break;  // is this necessary ?
		C = B;  B = A;
		A = msh .cell_behind ( B, tag::surely_exists ) .base() .reverse(); }    }

//-----------------------------------------------------------------------------------------------


class ManifoldNoWinding::sq_dist

// a callable object returning the square of the distance between two points
// used for MetricTree, see paragraphs 12.10 and 12.11 in the manual
 
{	public :

	inline double operator() ( const Cell & A, const Cell & B )
	{	return Manifold::working .sq_dist  ( A, B );  }
		
};  // end of  class ManifoldNoWinding::sq_dist

//------------------------------------------------------------------------------------------------------//


class ManifoldQuotient  // hidden in anonymous namespace

{	public :

	class sq_dist;	
	typedef std::pair < Cell, Manifold::Action > winding_cell;
	typedef std::pair < std::vector < double >, Manifold::Action > winding_vector;
	typedef MetricTree < winding_cell, sq_dist > metric_tree;

	inline static const Cell & remove_winding ( const winding_cell & P )
	{	return P .first;  }
	inline static const std::vector < double > & remove_winding ( const winding_vector & v )
	{	return v .first;  }

	inline static winding_vector get_vector ( const Cell & A, const Cell & B );
			
	inline static Cell build_seg ( const Cell & A, const Cell & B, const winding_vector & e )
	{	Cell seg ( tag::segment, A, B );   // A is negative
		seg .winding() = e .second;
		return seg;                      }

	inline static void redistribute_vertices
	( const Mesh & msh, const Cell & start, const Cell & stop, size_t n,
	  const winding_vector & e                                          );
		
};  // end of  class ManifoldQuotient

// we only need the winding for searching close neighbours of a given vertex
// we need 'sq_dist' to keep the "winning" winding, see below
	

//------------------------------------------------------------------------------------------------------//


inline std::pair < std::vector < double >, Manifold::Action > ManifoldQuotient::get_vector
( const Cell & A, const Cell & B )  // static

{	assert ( false );  return std::pair < std::vector < double >, Manifold::Action > { };  }
			
//------------------------------------------------------------------------------------------------------//

	
inline void ManifoldQuotient::redistribute_vertices  // static
( const Mesh & msh, const Cell & start, const Cell & stop, size_t n,
  const ManifoldQuotient::winding_vector & e                        )

// just make some baricenters; last argument provides winding number
	
{	Cell A = stop;
	// just in case n is too large, or the curve is too short, we look for 'start'
	// the statement below is important for closed loops, where start == stop
	A = msh .cell_behind (A) .base() .reverse();
	for ( size_t i = 2; i < n; i++ )
	{	if ( A == start )  {  n = i;  break;  }
		A = msh .cell_behind ( A, tag::surely_exists ) .base() .reverse();  }
	assert ( n > 2 );
	Manifold::Action w = e .second;
	Cell C = stop;
	Cell B = msh .cell_behind ( A, tag::surely_exists ) .base() .reverse();
	A = msh .cell_behind ( B, tag::surely_exists ) .base() .reverse();
	for ( size_t i = 2; i < n; i++ )
	{	Manifold::working .interpolate ( B, 0.3, A, 0.4, B, tag::winding, 0, 0.3, C, tag::winding, w );
		// the above method is not compatible with a non-uniform metric
		// so it will produce imperfect (but hopefully acceptable) results
		w = 0;
		assert ( A != start );
		if ( A == start ) break;  // is this necessary ?
		C = B;  B = A;
		A = msh .cell_behind ( B, tag::surely_exists ) .base() .reverse(); }                             }

//------------------------------------------------------------------------------------------------------//


class ManifoldQuotient::sq_dist

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

//------------------------------------------------------------------------------------------------------//


inline double approx_sqrt ( const double x )  // hidden in anonymous namespace

// a good approximation of the square root of x, for x between 0.25 and 4.

{	assert ( ( x > 0.25 ) and ( x < 4. ) );
	constexpr double coef = 0.27;
	constexpr double coef1 = 0.5 - coef;
	constexpr double coef2 = 4.*coef;
	const double tmp = x + 1.;
	return coef1 * tmp + coef2 * x / tmp;   }


inline void project_tangent  // hidden in anonymous namespace
( const Cell & A, const Cell & B, std::vector < double > & tangent )

// modifies 'tangent', sets coordinates of B
// often, B == temporary_vertex
	
{	for ( size_t i = 0; i < frontal_nb_of_coords; i++ )
	{	const Function & x = Manifold::working .coordinates() [i];
		x (B) = x(A) + tangent [i];                               }
	Manifold::working .project (B);
	for ( size_t i = 0; i < frontal_nb_of_coords; i++ )
	{	const Function & x = Manifold::working .coordinates() [i];
		tangent[i] = x (B) - x(A);                                }  }

	
inline void improve_tangent  // hidden in anonymous namespace
( const Cell & A, std::vector < double > & tangent )
	
// ensure the norm is 1., project, ensure again the norm is 1.

// in this version we do not assume the norm of 'tangent' is close to 1.
// also, 'tangent' may be far from the tangent line or tangent plane
	
{	double n2 = Manifold::working .core->metric->inner_prod ( A, tangent, tangent );
	// tangent may have norm far away from 1., so we use true square root below
	double norm = std::sqrt ( n2 );
	for ( size_t i = 0; i < frontal_nb_of_coords; i++ )  tangent [i] /= norm;
	project_tangent ( A, temporary_vertex, tangent );  // modifies 'tangent', sets coordinates of temp_ver
  n2 = Manifold::working .core->metric->sq_dist ( A, temporary_vertex, tangent );
	// because the 'tangent' provided may not be tangent to the working manifold,
	// we use the true square root below
	norm = std::sqrt ( n2 );
	for ( size_t i = 0; i < frontal_nb_of_coords; i++ )  tangent[i] /= norm;         }
		
	
inline Cell project_vertex_forward  // hidden in anonymous namespace
( const Cell & A, std::vector < double > & tangent )

// similar to 'improve_tangent'
// we assume the norm of 'tangent' is close to 1.
// we also assume it is nearly tangent to the working manifold
	
{	double n2 = Manifold::working .core->metric->inner_prod ( A, tangent, tangent );
	double norm = approx_sqrt ( n2 );
	Cell B ( tag::vertex );
	for ( size_t i = 0; i < frontal_nb_of_coords; i++ )  tangent [i] /= norm;
	project_tangent ( A, B, tangent );  // modifies 'tangent', sets coordinates of B
  n2 = Manifold::working .core->metric->sq_dist ( A, B, tangent );
	norm = approx_sqrt ( n2 );
	for ( size_t i = 0; i < frontal_nb_of_coords; i++ )
	{	const Function & x = Manifold::working .coordinates() [i]; 
		x(B) = x(A) + tangent [i] / norm;                          }
	Manifold::working .project (B);
	return B;                                                                         }
	
//------------------------------------------------------------------------------------------------------//


template < class manif_type >                // line 357
void frontal_construct          // hidden in anonymous namespace
( Mesh & msh, const tag::StartWithInconsistentMesh &,
	const tag::StartAt &, const Cell & start,
	const tag::Towards &, std::vector < double > tangent,
	const tag::StopAt &, const Cell & stop               )

// manif_type could be ManifoldNoWinding or ManifoldQuotient (defined above)
	
// builds a one-dimensional mesh (a curve)
// orientation given by 'tangent'
	
// 'start' and 'stop' are positive vertices (may be one and the same)
	
{	assert ( start .dim() == 0 );
	assert ( stop  .dim() == 0 );
	assert ( start .is_positive() );
	assert ( stop  .is_positive() );
	assert ( msh .dim() == 1 );

	improve_tangent ( start, tangent );  // modifies 'tangent'

	size_t counter = 1;
	// size_t max_counter = 0;  //  std::cin >> max_counter;
	Cell A = start;
	// the manifold's metric has been scaled, so the desired distance is now one
	while ( true )
	{	typename manif_type::winding_vector e = manif_type::get_vector ( A, stop );
		const std::vector < double > & ee = manif_type::remove_winding (e);
		if ( Manifold::working .sq_dist ( A, stop, ee ) < 1.44 )
		// distance is less than 1.2, we may stop now
		{	// check that 'stop' is in front of us
			// because, at the beginning, it could be behind us (or, even worse, A == stop)
			// below we could use the Riemannian inner product
			// however, the sign should be the same so it does not really matter
			double prod = 0.;
			for ( size_t i = 0; i < frontal_nb_of_coords; i++ )
				prod += tangent[i] * ee[i];
			// assert ( ( std::abs ( prod - 1.) < 0.3 ) or ( std::abs ( prod + 1.) < 0.3 ) );
			if ( prod > 0. )  // yes, 'stop' is in front of us
			{	Cell last = manif_type::build_seg ( A.reverse(), stop, e );
				// for a quotient manifold, 'e' provides the winding number
				last .add_to_mesh ( msh, tag::do_not_bother );
				// the meaning of tag::do_not_bother is explained
				// at the end of paragraph 11.6 in the manual
				manif_type::redistribute_vertices ( msh, start, stop, 6, e );
				// for a quotient manifold, 'e' provides the winding number
				return;                                                       }        }
		Cell B = project_vertex_forward ( A, tangent );  // modifies 'tangent'
		Cell AB ( tag::segment, A.reverse(), B );        // zero winding
		AB .add_to_mesh ( msh, tag::do_not_bother );
		// the meaning of tag::do_not_bother is explained at the end of paragraph 11.6 in the manual
		counter++;  A = B;                                                                      }

	assert ( false );
	
} // end of  frontal_construct

//------------------------------------------------------------------------------------------------------//


std::vector < double > compute_tangent_vec      // hidden in anonymous namespace
(	const Cell & start, bool check_orth, std::vector < double > given_vec )
	
// computes a vector tangent to Manifold::working at point 'start'
// here the working manifold is not a quotient manifold

// if second argument is true, candidates will be projected onto the space orthogonal to given_vec
// given_vec must be tangent to Manifold::working at point 'start'
// and must have length (approximately) equal to 1.
	
{	// Manifold::Implicit * m_impl =  dynamic_cast<Manifold::Implicit*> ( Manifold::working.core );
	// assert ( m_impl );
	// Manifold::Euclid * m_euclid =
	// 	dynamic_cast<Manifold::Euclid*> ( m_impl->surrounding_space.core );
	// assert ( m_euclid );

	std::vector < double > best_tangent;
	double longest_projection = 0.;
	// 'direc' contains 8 directions for 2D, 26 directions for 3D
	const std::vector < std::vector < double > > & direc =
		tag::Util::directions [ frontal_nb_of_coords ];
	const size_t n_dir = direc .size();
	for ( size_t dir = 1; dir < n_dir; dir++ )
	{	std::vector < double > tangent = direc [dir];
		project_tangent ( start, temporary_vertex, tangent );  // modifies 'tangent'
		if ( check_orth )
		{	double prod = Manifold::working .inner_prod ( start, tangent, given_vec );
			assert ( false );
			for ( size_t i = 0; i < frontal_nb_of_coords; i++ )
				tangent [i] -= prod * given_vec [i];                                    }
		// we choose the longest projection
		double n2 = Manifold::working.inner_prod ( start, tangent, tangent );
		if ( n2 > longest_projection )
		{	best_tangent = tangent;  longest_projection = n2;  }                            }

	// normalize best_tangent
	double norm = Manifold::working.inner_prod ( start, best_tangent, best_tangent );
	norm = std::sqrt ( norm );
	for ( size_t i = 0; i < frontal_nb_of_coords; i++ )  best_tangent [i] /= norm;
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
	

inline std::vector < double > compute_tangent_vec ( const tag::AtPoint &, Cell start )
// hidden in anonymous namespace

// computes a vector tangent to Manifold::working at point 'start'

{	return compute_tangent_vec ( start, false, std::vector<double>() );  }
	// 'false' as second argument means "do not check orthogonality"

}  // anonymous namespace

//------------------------------------------------------------------------------------------------------//


void frontal_construct          // hidden in anonymous namespace
( Mesh & msh, const tag::StartWithInconsistentMesh &,
  const tag::StartAt &, const Cell & start,
  const tag::StopAt &, const Cell & stop,
  const tag::Orientation &, const tag::OrientationChoice & oc )

// 'start' and 'stop' are positive vertices (may be one and the same)

// here the working manifold is not a quotient manifold

// uses two new temporary vertices

{	assert ( start .dim() == 0 );
	assert ( stop  .dim() == 0 );
	assert ( start .is_positive() );
	assert ( stop  .is_positive() );
	assert ( msh .dim() == 1 );

	// here the working manifold is not a quotient manifold
	#ifndef NDEBUG
	{ // just a block of code for hiding 'm'
	Manifold::Quotient * m = dynamic_cast < Manifold::Quotient * > ( Manifold::working .core );
	assert ( m == nullptr );
	} // just a block of code
	#endif

	if ( oc == tag::not_provided )
	{	std::cout << "when starting and stopping points are provided," << std::endl;
		std::cout << "maniFEM needs to know how to choose the orientation of the curve;" << std::endl;
		std::cout << "please specify either tag::orientation or tag::shortest_path" << std::endl;
		exit (1);                                                                                      }
	
	std::vector < double > best_tangent = compute_tangent_vec ( tag::at_point, start );

	if ( oc == tag::geodesic )   // shortest path

	{	assert ( start != stop );
		
		// start walking along the manifold from 'start' in the direction of best_tangent
		// and, simultaneously, in the opposite direction, given by  - best_tangent
		std::vector < double > tan1 = best_tangent, tan2 = best_tangent;
		for ( size_t i = 0; i < frontal_nb_of_coords; i++ ) tan2 [i] *= -1.;
		Cell ver1 ( tag::vertex ), ver2 ( tag::vertex );
		for ( size_t i = 0; i < frontal_nb_of_coords; i++ )
		{	Function x = Manifold::working .coordinates() [i];
			x ( ver1 ) = x ( start );  x ( ver2 ) = x ( start );  }
		double winner = 0.;  //  will be 1. or -1.
		while ( true )
		{	for ( size_t i = 0; i < frontal_nb_of_coords; i++ )
			{	Function x = Manifold::working .coordinates() [i];
				x ( temporary_vertex ) = x ( ver1 ) + tan1 [i];    }
			Manifold::working .project ( temporary_vertex );
			for ( size_t i = 0; i < frontal_nb_of_coords; i++ )
			{	Function x = Manifold::working .coordinates() [i];
				tan1 [i] = x ( temporary_vertex ) - x ( ver1 );
				x ( ver1 ) = x ( temporary_vertex );               }
			double sd = Manifold::working .sq_dist ( ver1, stop );
			if ( sd < 1.44 )  // dist < 1.2
			{	// check that 'stop' is in front of us
				// below we could use the Riemannian inner product
				// however, the sign should be the same so it does not really matter
				double prod = 0.;
				for ( size_t i = 0; i < frontal_nb_of_coords; i++ )
				{	Function x = Manifold::working .coordinates() [i];
					prod += tan1[i] * ( x (stop) - x (ver1) );         }
				if ( prod > 0. )  { winner = 1.;  break;  }             }
			for ( size_t i = 0; i < frontal_nb_of_coords; i++ )
			{	Function x = Manifold::working .coordinates() [i];
				x ( temporary_vertex ) = x ( ver2 ) + tan2 [i];  }
			Manifold::working .project ( temporary_vertex );
			for ( size_t i = 0; i < frontal_nb_of_coords; i++ )
			{	Function x = Manifold::working .coordinates() [i];
				tan2 [i] = x ( temporary_vertex ) - x ( ver2 );
				x ( ver2 ) = x ( temporary_vertex );            }
			sd = Manifold::working .sq_dist ( ver2, stop );
			if ( sd < 1.44 )  // dist < 1.2
			{	// check that 'stop' is in front of us
				// below we could use the Riemannian inner product
				// however, the sign should be the same so it does not really matter
				double prod = 0.;
				for ( size_t i = 0; i < frontal_nb_of_coords; i++ )
				{	Function x = Manifold::working .coordinates() [i];
					prod += tan2 [i] * ( x (stop) - x (ver2) );        }
				if ( prod > 0. )  { winner = -1.;  break;  }            }
		}  // end of  while true

		assert ( ( winner == 1. ) or ( winner == -1. ) );
		for ( size_t i = 0; i < frontal_nb_of_coords; i++ ) best_tangent [i] *= winner;
		frontal_construct < ManifoldNoWinding >
		( msh, tag::start_with_inconsistent_mesh, tag::start_at, start,
		       tag::towards, best_tangent,        tag::stop_at, stop   );
		return;
	}  // end of  if

	if ( ( oc == tag::inherent ) or ( oc == tag::random ) )
	{	}  // see tmp.cpp

	assert ( oc == tag::intrinsic );
	assert ( false );
	std::cout << "intrinsic orientation does not make sense for one-dimensional manifolds"
	          << std::endl;   // ???
	exit (1);

} // end of  frontal_construct

//------------------------------------------------------------------------------------------------------//


void Manifold::Core::frontal_method  // virtual, overridden by Manifold::Quotient
( Mesh & msh, const tag::StartWithInconsistentMesh &,
  const tag::StartAt &, const Cell & start,
  const tag::Towards &, std::vector<double> tangent,
  const tag::StopAt &,  const Cell & stop            )

{	frontal_construct < ManifoldNoWinding >  // line 357
	( msh, tag::start_with_inconsistent_mesh,
	  tag::start_at, start, tag::towards, tangent,
	  tag::stop_at, stop                          );  }


void Manifold::Quotient::frontal_method  // virtual from Manifold::Core, here overridden
( Mesh & msh, const tag::StartWithInconsistentMesh &,
  const tag::StartAt &, const Cell & start,
  const tag::Towards &, std::vector<double> tangent,
  const tag::StopAt &,  const Cell & stop            )

{	frontal_construct < ManifoldQuotient >  // line 357
	( msh, tag::start_with_inconsistent_mesh,
	  tag::start_at, start, tag::towards, tangent,
	  tag::stop_at, stop                          );  }
	
//-------------------------------------------------------------------------------------------------


void Manifold::Core::frontal_method  // virtual, overridden by Manifold::Quotient
( Mesh & msh, const tag::StartWithInconsistentMesh &,
  const tag::StartAt &, const Cell & start,
  const tag::StopAt &,  const Cell & stop, const tag::ShortestPath & )

{	frontal_construct  // line ???
	( msh, tag::start_with_inconsistent_mesh,
	  tag::start_at, start, tag::stop_at, stop,
		tag::orientation, tag::geodesic          );  }


void Manifold::Quotient::frontal_method  // virtual from Manifold::Core, here overridden
( Mesh & msh, const tag::StartWithInconsistentMesh &,
  const tag::StartAt &, const Cell & start,
  const tag::StopAt &,  const Cell & stop, const tag::ShortestPath & )

{	assert ( false );  }
	
//-------------------------------------------------------------------------------------------------


Mesh::Mesh ( const tag::Frontal &, const tag::StartAt &, const Cell & start,
             const tag::Towards &, std::vector<double> tangent,
             const tag::StopAt &, const Cell & stop,
             const tag::DesiredLength &, const double length                )

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
	frontal_nb_of_coords = Manifold::working .coordinates() .nb_of_components();

	// rescale the working manifold's metric
	// thus we may work with a desired distance of 1.
	tag::Util::Metric * old_metric = Manifold::working .core->metric;
	Manifold::working .core->metric = old_metric->scale ( 1. / length );

	// method below is virtual, calls Euclidian or Quotient version
	// we could have used instead a dynamic_cast followed by an if
	Manifold::working .core->frontal_method
	( *this, tag::start_with_inconsistent_mesh, tag::start_at, start,
	  tag::towards, tangent, tag::stop_at, stop                      );

	delete ( Manifold::working .core->metric );
	Manifold::working .core->metric = old_metric;

	update_info_connected_one_dim ( *this, start, stop );                              }
		
//------------------------------------------------------------------------------------------------------//


Mesh::Mesh ( const tag::Frontal &, const tag::StartAt &, const Cell & start,
             const tag::Towards &, std::vector<double> tangent,
             const tag::StopAt &, const Cell & stop,
             const tag::DesiredLength &, const Function & length            )

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
	frontal_nb_of_coords = Manifold::working .coordinates() .nb_of_components();

	// rescale the working manifold's metric
	// thus we may work with a desired distance of 1.
	tag::Util::Metric * old_metric = Manifold::working .core->metric;
	Manifold::working .core->metric = old_metric->scale ( 1. / length );

	// method below is virtual, calls Euclidian or Quotient version
	// we could have used instead a dynamic_cast followed by an if
	Manifold::working .core->frontal_method
	( *this, tag::start_with_inconsistent_mesh, tag::start_at, start,
	  tag::towards, tangent, tag::stop_at, stop                      );

	delete ( Manifold::working .core->metric );
	Manifold::working .core->metric = old_metric;

	update_info_connected_one_dim ( *this, start, stop );                              }
		
//------------------------------------------------------------------------------------------------------//


Mesh::Mesh ( const tag::Frontal &, const tag::StartAt &, const Cell & start,
             const tag::StopAt &, const Cell & stop,
             const tag::DesiredLength &, const double length,
             const tag::ShortestPath &                                       )

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
	frontal_nb_of_coords = Manifold::working .coordinates() .nb_of_components();

	// rescale the working manifold's metric
	// thus we may work with a desired distance of 1.
	tag::Util::Metric * old_metric = Manifold::working .core->metric;
	Manifold::working .core->metric = old_metric->scale ( 1. / length );

	// method below is virtual, calls Euclidian or Quotient version
	// we could have used instead a dynamic_cast followed by an if
	Manifold::working .core->frontal_method
	( *this, tag::start_with_inconsistent_mesh, tag::start_at, start,
	  tag::stop_at, stop, tag::shortest_path                         );

	delete ( Manifold::working .core->metric );
	Manifold::working .core->metric = old_metric;

	update_info_connected_one_dim ( *this, start, stop );                              }
		
//------------------------------------------------------------------------------------------------------//


Mesh::Mesh ( const tag::Frontal &, const tag::StartAt &, const Cell & start,
             const tag::StopAt &, const Cell & stop,
             const tag::DesiredLength &, const Function & length,
             const tag::ShortestPath &                                       )

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
	frontal_nb_of_coords = Manifold::working .coordinates() .nb_of_components();

	// rescale the working manifold's metric
	// thus we may work with a desired distance of 1.
	tag::Util::Metric * old_metric = Manifold::working .core->metric;
	Manifold::working .core->metric = old_metric->scale ( 1. / length );

	// method below is virtual, calls Euclidian or Quotient version
	// we could have used instead a dynamic_cast followed by an if
	Manifold::working .core->frontal_method
	( *this, tag::start_with_inconsistent_mesh, tag::start_at, start,
	  tag::stop_at, stop, tag::shortest_path                         );

	delete ( Manifold::working .core->metric );
	Manifold::working .core->metric = old_metric;

	update_info_connected_one_dim ( *this, start, stop );                              }
		
//------------------------------------------------------------------------------------------------------//



