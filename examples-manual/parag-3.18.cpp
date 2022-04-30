
// example presented in paragraph 3.18 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// builds an infinity-shaped mesh

#include "maniFEM.h"
#include <set>
using namespace maniFEM;

void limit_number_of_neighbours ( Mesh msh );


int main ( )

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy [0], y = xy [1];

	double d = 0.15;
	
	Manifold small_circle_left = RR2 .implicit ( (x-2.)*(x-2.) + y*y == 1. );
	Mesh inner_left ( tag::frontal, tag::desired_length, d );
	Manifold small_circle_right = RR2 .implicit ( (x+1.)*(x+1.) + y*y == 1. );
	Mesh inner_right ( tag::frontal, tag::desired_length, d );

	Manifold big_circle_left  = RR2 .implicit ( (x-2.)*(x-2.) + y*y == 4. );
	Manifold big_circle_right = RR2 .implicit ( (x+1.)*(x+1.) + y*y == 4. );
	Manifold two_points ( tag::intersect, big_circle_left, big_circle_right );

	Cell A ( tag::vertex );  x(A) =  0.;  y(A) =  1.;  two_points .project (A);
	Cell B ( tag::vertex );  x(B) =  0.;  y(B) = -1.;  two_points .project (B);
  
	// alternative syntax :
	// Cell A ( tag::vertex, tag::of_coords, { 0.,  1.}, tag::project );
	// Cell B ( tag::vertex, tag::of_coords, { 0., -1.}, tag::project );

	// std::cout << x(A) << " " << y(A) << " " << x(B) << " " << y(B) << std::endl;

	big_circle_left .set_as_working_manifold();
	Mesh outer_left ( tag::frontal, tag::start_at, B, tag:: stop_at, A,
	                  tag::orientation, tag::inherent, tag::desired_length, d );
	big_circle_right .set_as_working_manifold();
	Mesh outer_right ( tag::frontal, tag::start_at, A, tag:: stop_at, B,
	                   tag::orientation, tag::inherent, tag::desired_length, d );

	Mesh bdry ( tag::join, inner_right .reverse(), inner_left .reverse(), outer_right, outer_left );

	RR2 .set_as_working_manifold();
	Mesh infinity ( tag::frontal, tag::boundary, bdry, tag::desired_length, d );
	limit_number_of_neighbours ( infinity );

	Mesh::Iterator it2 = infinity .iterator ( tag::over_vertices );
	for ( int i = 1; i < 20; i++ )
	for ( it2 .reset(); it2 .in_range(); it2++ )
	{	Cell ver = *it2;
		if ( ver .is_inner_to ( infinity ) ) infinity .baricenter ( ver );  }

	infinity .export_to_file ( tag::msh, "infinity.msh");
	std::cout << "produced file ininity.msh" << std::endl;
	
}  // end of main

//-----------------------------------------------------------------------------------------


inline bool flip_segment ( Mesh & msh, Cell & seg )

// flip 'seg' if it is inner to 'msh'
// equilibrates (baricenter) the four neighbour vertices

// return true if the segment has been flipped, false if not

// assumes there are only triangular cells
	
// assumes there is no higher-dimensional mesh "above" 'msh'

{	Cell tri2 = msh.cell_in_front_of ( seg, tag::may_not_exist );
	if ( not tri2 .exists() ) return false; 
	Cell tri1 = msh .cell_behind ( seg, tag::may_not_exist );
	if ( not tri1 .exists() ) return false;
	// or, equivalently :  if ( not seg .is_inner_to ( msh ) ) return false

	Cell A = seg .base() .reverse();
	Cell B = seg .tip();

	Cell BC = tri1 .boundary() .cell_in_front_of ( B, tag::surely_exists );
	Cell CA = tri1 .boundary() .cell_behind ( A, tag::surely_exists );
	Cell AD = tri2 .boundary() .cell_in_front_of ( A, tag::surely_exists );
	Cell DB = tri2 .boundary() .cell_behind ( B, tag::surely_exists );
	Cell C = BC .tip();
	assert ( CA .base() .reverse() == C );
	Cell D = AD .tip();
	assert ( DB .base() .reverse() == D );
	
	B .cut_from_bdry_of ( seg, tag::do_not_bother );
	A .reverse() .cut_from_bdry_of ( seg, tag::do_not_bother );
	CA .cut_from_bdry_of ( tri1, tag::do_not_bother );
	DB .cut_from_bdry_of ( tri2, tag::do_not_bother );
	C .reverse() .glue_on_bdry_of ( seg, tag::do_not_bother );
	D .glue_on_bdry_of ( seg, tag::do_not_bother );
	DB .glue_on_bdry_of ( tri1, tag::do_not_bother );
	CA .glue_on_bdry_of ( tri2, tag::do_not_bother );

	tri1 .boundary() .closed_loop ( B );
	tri2 .boundary() .closed_loop ( A );

	if ( A .is_inner_to ( msh ) ) msh .baricenter ( A );
	if ( B .is_inner_to ( msh ) ) msh .baricenter ( B );
	if ( C .is_inner_to ( msh ) ) msh .baricenter ( C );
	if ( D .is_inner_to ( msh ) ) msh .baricenter ( D );

	return true;

}  // end of  flip_segment
	
//-----------------------------------------------------------------------------------------


double length_square ( Cell AB )
	
{	Manifold space = Manifold::working;
	assert ( space .exists() );  // we use the current (quotient) manifold
	Function coords = space .coordinates();
	Cell A = AB .base() .reverse();
	Cell B = AB .tip();
	std::vector < double > A_co = coords (A);
	std::vector < double > B_co = coords (B);
	size_t n = A_co .size();
	assert ( n == B_co .size() );
	double len_AB_2 = 0.;
	for ( size_t i = 0; i < n; i++ )
		{	double v = B_co [i] - A_co [i];  len_AB_2 += v*v;  }
	return len_AB_2;                                        }
	

class compare_lenghts_of_segs

{ public :
	inline bool operator() ( Cell AB, Cell CD ) const
	{	return length_square ( AB ) > length_square ( CD );  }
};

//-----------------------------------------------------------------------------------//


void limit_number_of_neighbours ( Mesh msh )

// in a mesh of triangles, flip segments in order to prevent
// vertices with only four (or three) neighbours
// and vertices with eight neighbours or more
	
// assumes there is no higher-dimensional mesh "above" 'msh'

{	std::forward_list < Cell > has_few_neighbours, has_many_neighbours;

	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = msh .iterator ( tag::over_vertices );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell P = * it;
		if ( not P .is_inner_to ( msh ) ) continue;
		// how many neighbours does P have ?
		size_t counter = 0;
		Mesh::Iterator it_around_P = msh .iterator ( tag::over_segments, tag::around, P );
		for ( it_around_P .reset(); it_around_P .in_range(); it_around_P ++ )
			counter ++;
		if ( counter < 5 ) has_few_neighbours .push_front ( P );
		if ( counter > 7 ) has_many_neighbours .push_front ( P );                         }
	} // just a block of code for hiding 'it'

	// we use a map to order neighbour segments, we flip the longest one
	compare_lenghts_of_segs comp_len;
	for ( std::forward_list < Cell > ::iterator it = has_few_neighbours .begin();
				it != has_few_neighbours .end(); it ++                                 )
	{	Cell P = * it;
		// we count again the neighbours, configuration may have changed in the meanwhile
		size_t counter = 0;
		std::multiset < Cell, compare_lenghts_of_segs > ms ( comp_len );
		Mesh::Iterator it_around_P =
			msh .iterator ( tag::over_cells_of_dim, 2, tag::around, P );
		for ( it_around_P .reset(); it_around_P .in_range(); it_around_P ++ )
		{	Cell tri = * it_around_P;
			counter ++;
			assert ( tri.dim() == 2 );
			Cell SP = tri .boundary() .cell_behind ( P, tag::surely_exists );
			Cell PT = tri .boundary() .cell_in_front_of ( P, tag::surely_exists );
			Cell T = PT .tip();
			Cell TS = tri .boundary() .cell_in_front_of ( T, tag::surely_exists );
			assert ( TS .tip() == SP .base() .reverse() );
			ms .insert ( TS );                                                     }
		if ( counter > 4 ) continue;
		for ( std::multiset < Cell, compare_lenghts_of_segs > ::iterator itt = ms .begin();
					itt != ms .end(); itt ++                                                      )
		{	Cell seg = * itt;
			if ( flip_segment ( msh, seg ) )  break;  }                                         }
			// we choose to to flip only one segment

	for ( std::forward_list < Cell > ::iterator it = has_many_neighbours .begin();
				it != has_many_neighbours .end(); it ++                                 )
	{	Cell P = * it;
		// we count again the neighbours, configuration may have changed in the meanwhile
		size_t counter = 0;
		std::multiset < Cell, compare_lenghts_of_segs > ms ( comp_len );
		Mesh::Iterator it_around_P =
			msh .iterator ( tag::over_segments, tag::around, P );
		for ( it_around_P .reset(); it_around_P .in_range(); it_around_P ++ )
		{	Cell seg = * it_around_P;
			assert ( seg.tip() == P );
			counter ++;
			ms .insert ( seg );        }
		if ( counter < 8 ) continue;
		for ( std::multiset < Cell, compare_lenghts_of_segs > ::iterator itt = ms .begin();
					itt != ms .end(); itt ++                                                      )
		{	Cell seg = * itt;
			if ( flip_segment ( msh, seg ) )  break;  }                                         }
			// we choose to to flip only one segment

}  // end of  limit_number_of_neighbours
	
//-----------------------------------------------------------------------------------------

