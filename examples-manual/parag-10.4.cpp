
// example presented in paragraph 10.4 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// in a mesh of triangles, flip segments with the aim of controlling the number of neighbours

#include "maniFEM.h"
#include <set>
using namespace maniFEM;


void limit_number_of_neighbours ( Mesh msh );


int main ( )

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy [0],  y = xy [1];

	RR2 .implicit ( x*x + 0.7*y*y == 1. );

	Cell A ( tag::vertex );  x (A) = 1.;  y (A) = 0.;
	Mesh curve ( tag::progressive, tag::start_at, A, tag::desired_length, 0.33 );

	RR2 .set_as_working_manifold();
	Mesh ellipse ( tag::progressive, tag::boundary, curve, tag::desired_length, 0.33 );

	limit_number_of_neighbours ( ellipse );

	// in rare occasions, we must call this function twice :
	// limit_number_of_neighbours ( ellipse );

	ellipse .export_to_file ( tag::msh, "ellipse.msh");
	std::cout << "produced file ellipse.msh" << std::endl;

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

