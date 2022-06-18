
// frontal.cpp 2022.06.18

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
// (no need to change anything in MetricTree)

// for computing this winding number,
// implement a discrete version of the steepest descent method
// rather than a blind search as in draw_ps with tag::windng

//------------------------------------------------------------------------------------------------------//


// global variables and functions for this file, not visible for other object files
namespace {  // anonymous namespace, mimics static linkage

Cell temporary_vertex ( tag::non_existent );

size_t progress_nb_of_coords;  // dimension of the surrounding Euclidian space
// also known as "geometric dimension"

}  // anonymous namespace

//------------------------------------------------------------------------------------------------------//


class Manifold::Type::Euclidian

{	public :

	class sq_dist;
	typedef Cell winding_cell;
	typedef MetricTree < Cell, sq_dist > metric_tree;

	class ConstInnProd;
};
	
//------------------------------------------------------------------------------------------------------//

	
class Manifold::Type::Euclidian::sq_dist

// a callable object returning the square of the distance between two points
// used for MetricTree, see paragraphs 12.10 and 12.11 in the manual

// the inner product is not constant, it may be different in A from B
// this is why we return the arithmetic mean between the two products
// however, we prefer not to divide by two in order to save computing time
// the calling code must adjust the desired distance accordingly

{	public :

	inline double operator() ( const Cell & A, const Cell & B )
	{	const Function & coord = Manifold::working .coordinates();
		std::vector < double > vA = coord ( A ), vB = coord ( B ),
			delta ( coord .nb_of_components() );
		for ( size_t i = 0; i < coord .nb_of_components(); i++ )
			delta[i] = vB[i] - vA[i];
		return Manifold::working .inner_prod ( A, delta, delta ) +
		       Manifold::working .inner_prod ( B, delta, delta )  ;  }
		// return ( Manifold::working .inner_prod ( A, delta, delta ) +
		//          Manifold::working .inner_prod ( B, delta, delta )  ) / 2.;  }

};  // end of  class Manifold::Type::Euclidian::sq_dist

//------------------------------------------------------------------------------------------------------//


class Manifold::Type::Euclidian::ConstInnProd

{	public :

	class sq_dist;
	typedef Cell winding_cell;
	typedef MetricTree < Cell, sq_dist > metric_tree;
};
	
//------------------------------------------------------------------------------------------------------//

	
class Manifold::Type::Euclidian::ConstInnProd::sq_dist

// a callable object returning the square of the distance between two points
// used for MetricTree, see paragraphs 12.10 and 12.11 in the manual

// here the inner product is constant, doesn't matter if we use A or B

{	public :

	inline double operator() ( const Cell & A, const Cell & B )
	{	const Function & coord = Manifold::working .coordinates();
		std::vector < double > vA = coord ( A ), vB = coord ( B ),
			delta ( coord .nb_of_components() );
		for ( size_t i = 0; i < coord .nb_of_components(); i++ )
			delta[i] = vB[i] - vA[i];
		return Manifold::working .inner_prod ( A, delta, delta );  }

};  // end of  class Manifold::Type::Euclidian::ConstInnProd::sq_dist

//------------------------------------------------------------------------------------------------------//


class Manifold::Type::Quotient

{	public :

	class sq_dist;	
	typedef std::pair < Cell, Manifold::Action > winding_cell;
	typedef MetricTree < winding_cell, sq_dist > metric_tree;
};

// we only need the winding for searching close neighbours of a given vertex
// we need 'sq_dist' to keep the "winning" winding, see below
	

//------------------------------------------------------------------------------------------------------//


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

//------------------------------------------------------------------------------------------------------//


// there are two aspects to be taken into account

// first, we must distinguish between a "usual" manifold and a quotient one
// a "usual" manifold may be e.g. Euclidian or implicit
// the 'manif_type' will help us make this distinction
// it allows to have two types of cells (simple and rich, the latter having winding)

// second, we must distinguish between different types of metric :
// trivial, isotropic (constant zoom, variable zoom), anisotropic, Rayleigh
// virtual methods of the metric itself will allow us to treat differently these situations


namespace { // anonymous namespace, mimics static linkage

inline double approx_sqrt ( const double x )  // hidden in anonymous namespace

// a good approximation of the square root of x
// for x between 	0.25 and 4.

{	assert ( ( x > 0.25 ) and ( x < 4. ) );
	constexpr double coef = 0.27;
	constexpr double coef1 = 0.5 - coef;
	constexpr double coef2 = 4.*coef;
	const double tmp = x + 1.;
	return coef1 * tmp + coef2 * x / tmp;   }


inline void improve_tangent  // hidden in anonymous namespace
( const Cell & A, std::vector < double > & tangent )

// ensure the norm is 1., project, ensure again the norm is 1.

// in this version we do not assume the norm of 'tangent' is close to 1.
// also, 'tangent' may be far from the tangent line or tangent plane
	
{	double n2 = Manifold::working .core->metric->inner_prod ( A, tangent, tangent );
	// tangent may have norm far away from 1., so we use true square root below
	double norm = std::sqrt ( n2 );
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	const Function & x = Manifold::working .coordinates() [i];
		x ( temporary_vertex ) = x(A) + tangent [i] / norm;        }
	Manifold::working .project ( temporary_vertex );
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	const Function & x = Manifold::working .coordinates() [i];
		tangent[i] = x ( temporary_vertex ) - x(A);                }
  n2 = Manifold::working .core->metric->sq_dist ( A, temporary_vertex, tangent );
	// because the 'tangent' provided may not be tangent to the working manifold,
	// we use the true square root below
	norm = std::sqrt ( n2 );
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )  tangent[i] /= norm;       }
		
	
inline Cell project_vertex_forward  // hidden in anonymous namespace
( const Cell & A, std::vector < double > & tangent )

// similar to 'improve_tangent'
// we assume the norm of 'tangent' is close to 1.
// we also assume it is nearly tangent to the working manifold
	
{	double n2 = Manifold::working .core->metric->inner_prod ( A, tangent, tangent );
	double norm = aprox_sqrt ( n2 );
	Cell B ( tag::vertex );
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	const Function & x = Manifold::working .coordinates() [i];
		x(B) = x(A) + tangent [i] / norm;                          }
	Manifold::working .project (B);
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	const Function & x = Manifold::working .coordinates() [i];
		tangent[i] = x(B) - x(A);                                  }
  n2 = Manifold::working .core->metric->sq_dist ( A, B, tangent );
	norm = aprox_sqrt ( n2 );
	for ( size_t i = 0; i < progress_nb_of_coords; i++ )
	{	const Function & x = Manifold::working .coordinates() [i]; 
		x(B) = x(A) + tangent [i] / norm;                          }
	Manifold::working .project (B);
	return B;                                                          }
	
//------------------------------------------------------------------------------------------------------//


inline void redistribute_vertices
( const Mesh & msh, const Cell & start, const Cell & stop, size_t n )
// hidden in anonymous namespace    // called only once

// just make some baricenters
	
{	Cell A = stop;
	// just in case n is too large, or the curve is too short, we look for 'start'
	// the statement below is important for closed loops, where start == stop
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
		// the above method is not compatible with a non-uniform metric
		if ( C == start ) break;
		A = B;  B = C;
		C = msh .cell_behind ( B, tag::surely_exists ) .base() .reverse(); }    }

//-----------------------------------------------------------------------------------------------


template < class manif_type >                // line 289
void progressive_construct          // hidden in anonymous namespace
( Mesh & msh, const tag::StartWithInconsistentMesh &,
	const tag::StartAt &, const Cell & start,
	const tag::Towards &, std::vector < double > tangent,
	const tag::StopAt &, const Cell & stop               )

// builds a one-dimensional mesh (a curve)
// orientation given by 'tangent'
	
// 'start' and 'stop' are vertices (may be one and the same)
	
{	assert ( start .dim() == 0 );
	assert ( stop  .dim() == 0 );

	improve_tangent ( start, tangent);  // modifies 'tangent'

	size_t counter = 1;
	size_t max_counter = 0;  //  std::cin >> max_counter;
	Cell A = start;
	// the manifold's metric has been scaled
	// so the desired distance is now one
	while ( true )
	{	std::vector < double > e = Manifold::working .core->get_vector ( A, stop );
		if ( Manifold::working .core->dist_less_than_one ( A, stop, e ) )
			// check that stop' is in front of us, not behind (or, even worse, A == stop)
			// below we could use the Riemannian metric
			// however, the sign should be the same so it does not really matter
		{	double prod = 0.;
			for ( size_t i = 0; i < progress_nb_of_coords; i++ )
				prod += tangent[i] * e[i];
			assert ( ( std::abs ( prod - 1.) < 0.3 ) or ( std::abs ( prod + 1.) < 0.3 ) );
			if ( prod > 0. )
			{	Cell last ( tag::segment, A.reverse(), stop );  // how about winding numbers ?
				last .add_to_mesh ( msh, tag::do_not_bother );
				// the meaning of tag::do_not_bother is explained
				// at the end of paragraph 11.6 in the manual
				redistribute_vertices ( msh, start, stop, 6 );  // how about winding numbers ?
				return;                                        }      }
		Cell B = project_vertex_forward ( A, tangent );  // modifies 'tangent'
		Cell AB ( tag::segment, A.reverse(), B );        // zero winding
		AB .add_to_mesh ( msh, tag::do_not_bother );
		// the meaning of tag::do_not_bother is explained at the end of paragraph 11.6 in the manual
		counter++;  A = B;                                                                      }

	assert ( false );
	
} // end of  progressive_construct

}  // anonymous namespace

//------------------------------------------------------------------------------------------------------//


void Manifold::Core::frontal_method  // virtual, overridden by Manifold::Quotient
( Mesh & msh, const tag::StartWithInconsistentMesh &,
  const tag::StartAt &, const Cell & start,
  const tag::Towards &, std::vector<double> tangent,
  const tag::StopAt &,  const Cell & stop            )

{	progressive_construct < Manifold::Type::Euclidian >  // line 289
	( msh, tag::start_with_inconsistent_mesh,
	  tag::start_at, start, tag::towards, tangent,
	  tag::stop_at, stop                          );  }
	
//------------------------------------------------------------------------------------------------------//


void Manifold::Quotient::frontal_method  // virtual from Manifold, here overridden
( Mesh & msh, const tag::StartWithInconsistentMesh &,
  const tag::StartAt &, const Cell & start,
  const tag::Towards &, std::vector<double> tangent,
  const tag::StopAt &,  const Cell & stop            )

{	progressive_construct < Manifold::Type::Quotient >  // line 289
	( msh, tag::start_with_inconsistent_mesh,
	  tag::start_at, start, tag::towards, tangent,
	  tag::stop_at, stop                          );  }
	
//-------------------------------------------------------------------------------------------------


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
	progress_nb_of_coords = Manifold::working .coordinates() .nb_of_components();

	// rescale the working manifold's metric
	// thus we may work with a desired distance of 1.
	tag::Util::Metric * old_metric = Manifold::working .metric;
	Manifold::working .metric = old_metric->scale ( 1. / length );

	// method below is virtual, calls Euclidian or Quotient version
	Manifold::working .core->frontal_method
	( *this, tag::start_with_inconsistent_mesh, tag::start_at, start,
	  tag::towards, tangent, tag::stop_at, stop                      );

	delete ( Manifold::working .metric );
	Manifold::working .metric = old_metric;

	update_info_connected_one_dim ( msh, start, stop );                              }
		
//------------------------------------------------------------------------------------------------------//


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
	progress_nb_of_coords = Manifold::working .coordinates() .nb_of_components();

	// rescale the working manifold's metric
	// thus we may work with a desired distance of 1.
	tag::Util::Metric * old_metric = Manifold::working .metric;
	Manifold::working .metric = old_metric->scale ( 1. / length );

	// method below is virtual, calls Euclidian or Quotient version
	Manifold::working .core->frontal_method
	( *this, tag::start_with_inconsistent_mesh, tag::start_at, start,
	  tag::towards, tangent, tag::stop_at, stop                      );

	delete ( Manifold::working .metric );
	Manifold::working .metric = old_metric;

	update_info_connected_one_dim ( msh, start, stop );                              }
		
//------------------------------------------------------------------------------------------------------//



