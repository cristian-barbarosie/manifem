
// frontal.cpp 2022.05.16

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
	{	return Manifold::working .sq_dist ( A, B );  }

};  // end of  class Manifold::Type::Euclidian::sq_dist

//-----------------------------------------------------------------------------------------


class Manifold::Type::Quotient

{	public :

	class sq_dist;	
	typedef std::pair < Cell, Manifold::Action > winding_cell;
	typedef MetricTree < winding_cell, sq_dist > metric_tree;
};

// we only need the winding for searching close neighbours of a given vertex
// we need 'sq_dist' to keep the "winning" winding, see below
	

//-----------------------------------------------------------------------------------------


template < class manif_type >
void progressive_construct ( Mesh & msh,
	const tag::StartAt &, const Cell & start,
	const tag::Towards &, std::vector < double > & tangent,
	const tag::StopAt &, const Cell & stop                 )
// hidden in anonymous namespace

// builds a one-dimensional mesh (a curve)
// orientation given by 'tangent'
	
// 'start' and 'stop' are vertices (may be one and the same)
	
{	assert ( start .dim() == 0 );
	assert ( stop  .dim() == 0 );

	size_t counter = 1;
	size_t max_counter = 0;  //  std::cin >> max_counter;
	Cell A = start;
	desired_len_at_point = desired_length (A);
	sq_desired_len_at_point = desired_len_at_point * desired_len_at_point;
	while ( true )
	{	double d = Manifold::working .dist_sq ( A, stop );  // long range distance !
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
				break;                                          }      }
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
		AB .add_to_mesh ( msh, tag::do_not_bother );
		// the meaning of tag::do_not_bother is explained at the end of paragraph 11.6 in the manual
		counter++;  A = B;                                                                      }
	
	update_info_connected_one_dim ( msh, start, stop );
		
} // end of  progressive_construct

//-------------------------------------------------------------------------------------------------


void Manifold::Core::frontal_method  // virtual, overridden by Manifold::Quotient
( Mesh & msh, const tag::StartWithInconsistentMesh &,
  const tag::StartAt &, const Cell & start,
  const tag::Towards &, std::vector<double> tangent,
  const tag::StopAt &, const Cell & stop             )

{	progressive_construct < Manifold::Type::Euclidian >
	( *this, tag::start_at, start, tag::towards, tangent,
	  tag::stop_at, stop                                 );  }
	
//-------------------------------------------------------------------------------------------------


void Manifold::Quotient::frontal_method  // virtual, overridden by Manifold::Quotient
( Mesh & msh, const tag::StartWithInconsistentMesh &,
  const tag::StartAt &, const Cell & start,
  const tag::Towards &, std::vector<double> tangent,
  const tag::StopAt &, const Cell & stop             )

{	progressive_construct < Manifold::Type::Quotient >
	( *this, tag::start_at, start, tag::towards, tangent,
	  tag::stop_at, stop                                 );  }
	
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

	// method below is virtual, calls Euclidian or Quotient version
	Manifold::working .core->frontal_method
	( *this, tag::start_with_inconsistent_mesh, tag::start_at, start,
	  tag::towards, tangent, tag::stop_at, stop                      );              }

//------------------------------------------------------------------------------------------------------//



