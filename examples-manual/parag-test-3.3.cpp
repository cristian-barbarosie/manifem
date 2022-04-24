
// example presented in paragraph 3.3 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// builds an infinity-shaped overlaping mesh

#include "maniFEM.h"
#include <set>
using namespace maniFEM;

void limit_number_of_neighbours ( Mesh msh );


int main ( )

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy [0], y = xy [1];

	int n = 8;
	double d = 1./n;

	Cell O1 ( tag::vertex );  x(O1) = 0.;  y(O1) = 0.;
	Cell A1 ( tag::vertex );  x(A1) = 1.;  y(A1) = 0.;
	Mesh O1A1 ( tag::segment, O1 .reverse(), A1, tag::divided_in, n );
	Cell O2 ( tag::vertex );  x(O2) = 0.;  y(O2) = 0.;
	Cell A2 ( tag::vertex );  x(A2) = 1.;  y(A2) = 0.;
	Mesh O2A2 ( tag::segment, O2 .reverse(), A2, tag::divided_in, n );

	// although they occupy the same position in RR2,
	// the two segments O1A1 and O2A2 are entirely different
	// they have no common vertex and no common segment
	
	// you may want to give a slight curvature to one of the segments O1A1, O2A2
	// or you may want to curve one of them upwards and the other downwards
	// without changing the position of the extremities O1, O2, A1, A2
	// this may improve the look of the final mesh
	
	Cell B ( tag::vertex );  x(B) =  3.;  y(B) = 0.;
	Cell C ( tag::vertex );  x(C) =  4.;  y(C) = 0.;
	Mesh BC ( tag::segment, B .reverse(), C, tag::divided_in, n );
	Cell D ( tag::vertex );  x(D) = -2.;  y(D) = 0.;
	Cell E ( tag::vertex );  x(E) = -3.;  y(E) = 0.;
	Mesh DE ( tag::segment, D .reverse(), E, tag::divided_in, n );
	
	Manifold big_circle_left = RR2 .implicit ( (x-2.)*(x-2.) + y*y == 4. );
	Mesh outer_left_down ( tag::frontal, tag::start_at, O1, tag::towards, { 0., -1. },
	                      tag::stop_at, C, tag::desired_length, d                    );
	Mesh outer_left_up ( tag::frontal, tag::start_at, C, tag::towards, { 0., 1. },
	                     tag::stop_at, O2, tag::desired_length, d                 );
	Manifold small_circle_left = RR2 .implicit ( (x-2.)*(x-2.) + y*y == 1. );
	Mesh inner_left_down ( tag::frontal, tag::start_at, B, tag::towards, { 0., -1. },
	                      tag::stop_at, A1, tag::desired_length, d                  );
	Mesh inner_left_up ( tag::frontal, tag::start_at, A2, tag::towards, { 0., 1. },
	                     tag::stop_at, B, tag::desired_length, d                   );

	Manifold big_circle_right = RR2 .implicit ( (x+1.)*(x+1.) + y*y == 4. );
	Mesh outer_right_down ( tag::frontal, tag::start_at, E, tag::towards, { 0., -1. },
	                        tag::stop_at, A2, tag::desired_length, d                  );
	Mesh outer_right_up ( tag::frontal, tag::start_at, A1, tag::towards, { 0., 1. },
	                      tag::stop_at, E, tag::desired_length, d                   );
	Manifold small_circle_right = RR2 .implicit ( (x+1.)*(x+1.) + y*y == 1. );
	Mesh inner_right_down ( tag::frontal, tag::start_at, O2, tag::towards, { 0., -1. },
	                        tag::stop_at, D, tag::desired_length, d                    );
	Mesh inner_right_up ( tag::frontal, tag::start_at, D, tag::towards, { 0., 1. },
	                      tag::stop_at, O1, tag::desired_length, d                 );

	RR2 .set_as_working_manifold();
	Mesh bdry_left_up ( tag::join, outer_left_up, O2A2, inner_left_up, BC );
	Mesh left_up ( tag::frontal, tag::boundary, bdry_left_up, tag::desired_length, d );
	Mesh bdry_left_down ( tag::join, outer_left_down, BC .reverse(), inner_left_down, O1A1 .reverse() );
	Mesh left_down ( tag::frontal, tag::boundary, bdry_left_down, tag::desired_length, d );
	Mesh bdry_right_up ( tag::join, outer_right_up, DE .reverse(), inner_right_up, O1A1 );
	Mesh right_up ( tag::frontal, tag::boundary, bdry_right_up, tag::desired_length, d );
	Mesh bdry_right_down ( tag::join, outer_right_down, O2A2 .reverse(), inner_right_down, DE );
	Mesh right_down ( tag::frontal, tag::boundary, bdry_right_down, tag::desired_length, d );

	Mesh infinity ( tag::join, left_up, left_down, right_up, right_down );
	limit_number_of_neighbours ( infinity );

	Mesh::Iterator it2 = infinity .iterator ( tag::over_vertices );
	for ( int i = 1; i < 20; i++ )
	for ( it2 .reset(); it2 .in_range(); it2++ )
	{	Cell ver = *it2;
		if ( ver .is_inner_to ( infinity ) ) infinity .baricenter ( ver );  }

	infinity .export_to_file ( tag::msh, "infinity.msh");
	std::cout << "produced file ininity.msh" << std::end;

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

