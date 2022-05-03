
// example presented in paragraph 7.23 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// build an intricate mesh on a torus, having two V-shaped holes
// homogenized elastic tensor exhibits negative Poisson coefficient

// frontal mesh generation does not work yet on quotient manifolds
// we fake the result by folding a mesh whose exterior boundary is a polygonal line
// see https://github.com/cristian-barbarosie/attic, file five-V-holes.png

#include "maniFEM.h"
#include <set>

using namespace maniFEM;

void limit_number_of_neighbours ( Mesh msh );
	

int main ( )

{	// begin with the usual two-dimensional space
	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0], y = xy[1];

	Cell A ( tag::vertex );  x(A) =  0.  ;  y(A) = -0.35;
	Cell B ( tag::vertex );  x(B) =  0.  ;  y(B) = -1.05;
	Cell C ( tag::vertex );  x(C) =  1.15;  y(C) = -2.15;
	Cell D ( tag::vertex );  x(D) =  2.3 ;  y(D) = -1.05;
	Cell E ( tag::vertex );  x(E) =  2.3 ;  y(E) = -0.35;
	Cell F ( tag::vertex );  x(F) =  1.15;  y(F) = -1.45;
	Cell G ( tag::vertex );  x(G) =  1.15;  y(G) = -0.75;
	Cell H ( tag::vertex );  x(H) =  0.  ;  y(H) =  0.35;
	Cell I ( tag::vertex );  x(I) = -1.15;  y(I) = -0.75;
	Cell J ( tag::vertex );  x(J) = -1.15;  y(J) = -1.45;

	Mesh AB ( tag::segment, A .reverse(), B, tag::divided_in, 14 );
	Mesh BC ( tag::segment, B .reverse(), C, tag::divided_in, 32 );
	Mesh CD ( tag::segment, C .reverse(), D, tag::divided_in, 32 );
	Mesh DE ( tag::segment, D .reverse(), E, tag::divided_in, 14 );
	Mesh EF ( tag::segment, E .reverse(), F, tag::divided_in, 32 );
	Mesh FG ( tag::segment, F .reverse(), G, tag::divided_in, 14 );
	Mesh GH ( tag::segment, G .reverse(), H, tag::divided_in, 32 );
	Mesh HI ( tag::segment, H .reverse(), I, tag::divided_in, 32 );
	Mesh IJ ( tag::segment, I .reverse(), J, tag::divided_in, 14 );
	Mesh JA ( tag::segment, J .reverse(), A, tag::divided_in, 32 );

	Mesh zig_zag ( tag::join, { AB, BC, CD, DE, EF, FG, GH, HI, IJ, JA } );

	const double e = 1.5;
	Manifold curve1 = RR2 .implicit 
		( smooth_min ( 300.*power((x+y)*(x+y),e) + power((x-y-1.)*(x-y-1.),e),
		               300.*power((x-y)*(x-y),e) + power((x+y+1.)*(x+y+1.),e),
	                 tag::threshold, 1.5                     )  == 1. );
	Mesh loop1 ( tag::frontal, tag::desired_length, 0.05 );

	double a = 0.15, b = 0.82;
	Manifold curve2 = RR2 .implicit 
		( smooth_min ( 300.*power((x-a+y+b)*(x-a+y+b),e) + power((x-a-y-b-1.)*(x-a-y-b-1.),e),
		        300.*power((x-2.-a-y-b)*(x-2.-a-y-b),e) + power((x-2.-a+y+b+1.)*(x-2.-a+y+b+1.),e),
	                 tag::threshold, 1.5                     )  == 1. );
	Mesh loop2 ( tag::frontal, tag::desired_length, 0.05 );

	b = 1.4;
	Manifold curve3 = RR2 .implicit 
		( smooth_min ( 300.*power((x+y+b)*(x+y+b),e) + power((x-y-b-1.)*(x-y-b-1.),e),
		               300.*power((x-y-b)*(x-y-b),e) + power((x+y+b+1.)*(x+y+b+1.),e),
	                 tag::threshold, 1.5                     )  == 1. );
	Mesh loop3 ( tag::frontal, tag::desired_length, 0.05 );  // not used

	b = 0.82-1.4;
	Manifold curve4 = RR2 .implicit 
		( smooth_min ( 300.*power((x-a+y+b)*(x-a+y+b),e) + power((x-a-y-b-1.)*(x-a-y-b-1.),e),
		        300.*power((x-2.-a-y-b)*(x-2.-a-y-b),e) + power((x-2.-a+y+b+1.)*(x-2.-a+y+b+1.),e),
	                 tag::threshold, 1.5                     )  == 1. );
	Mesh loop4 ( tag::frontal, tag::desired_length, 0.05 );  // not used

	a -= 2.3;  b = 0.82;
	Manifold curve5 = RR2 .implicit 
		( smooth_min ( 300.*power((x-a+y+b)*(x-a+y+b),e) + power((x-a-y-b-1.)*(x-a-y-b-1.),e),
		        300.*power((x-2.-a-y-b)*(x-2.-a-y-b),e) + power((x-2.-a+y+b+1.)*(x-2.-a+y+b+1.),e),
	                 tag::threshold, 1.5                     )  == 1. );
	Mesh loop5 ( tag::frontal, tag::desired_length, 0.05 );  // not used

	Mesh loop ( tag::join, loop1 .reverse(), loop2 .reverse(), zig_zag );

	RR2 .set_as_working_manifold();
	Mesh domain ( tag::frontal, tag::boundary, loop, tag::desired_length, 0.05 );

	Mesh IJAB ( tag::join, IJ, JA, AB );
	Mesh GFED ( tag::join, DE .reverse(), EF .reverse(), FG .reverse() );

	Mesh VV = domain.fold ( tag::identify, IJAB, tag::with, GFED,
                          tag::identify, BC, tag::with, GH .reverse(),
													tag::identify, CD, tag::with, HI .reverse(),
													tag::use_existing_vertices                  );

	std::cout << "produced folded mesh, now drawing, please wait" << std::endl << std::flush;
	VV .draw_ps ( "VV.eps", tag::unfold, tag::over_region, -1.5 < x < 2.5, -2 < y < 0.5 );
	
	std::cout << "now smoothening ... " << std::flush;

	limit_number_of_neighbours ( VV );
	
	std::cout << "and drawing again, please wait" << std::endl << std::flush;

	VV .draw_ps ( "VV-smooth.eps", tag::unfold,
               tag::over_region, -1.5 < x < 2.5, -2 < y < 0.5 );

	std::cout << "produced files VV.eps and VV-smooth.eps - please edit before viewing" << std::endl;	

}  // end of  main

// looks like there is a problem in split_segment or in remove_short_segments
// shows up with MANIFEM_COLLECT_CM on


//-----------------------------------------------------------------------------------------

void remove_short_segments ( Mesh & msh, double threshold );

void flip_split_long_segments ( Mesh & msh, double threshold );

void baricenters ( Mesh & msh );

//-----------------------------------------------------------------------------------------


inline bool flip_segment ( Mesh & msh, Cell & seg )

// flip 'seg' if it is inner to 'msh'
// equilibrates (baricenter) the four neighbour vertices

// return true if the segment has been flipped, false if not

// assumes there are only triangular cells
	
// assumes there is no higher-dimensional mesh "above" 'msh'

// assumes that the current manifold is a quotient manifold (manipulates windings)

{	Cell tri2 = msh .cell_in_front_of ( seg, tag::may_not_exist );
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
	
	assert ( seg .winding() == - BC .winding() - CA .winding() );
	assert ( seg .winding() == AD .winding() + DB .winding() );
	seg .winding() = CA .winding() + AD .winding();
	assert ( seg .winding() == - BC .winding() - DB .winding() );

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

	if ( A .is_inner_to ( msh ) ) msh .baricenter ( A, tag::winding );
	if ( B .is_inner_to ( msh ) ) msh .baricenter ( B, tag::winding );
	if ( C .is_inner_to ( msh ) ) msh .baricenter ( C, tag::winding );
	if ( D .is_inner_to ( msh ) ) msh .baricenter ( D, tag::winding );

	return true;

}  // end of  flip_segment
	
//-----------------------------------------------------------------------------------------


double length_square ( Cell AB )
	
{	Manifold space = Manifold::working;
	assert ( space.exists() );  // we use the current (quotient) manifold
	Manifold::Quotient * mani_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space .core );
	Function coords_q = space .coordinates();
	Manifold mani_Eu = mani_q->base_space;  // underlying Euclidian manifold
	Function coords_Eu = mani_Eu .coordinates();
	Cell A = AB .base() .reverse();
	Cell B = AB .tip();
	Manifold::Action s = AB .winding();
	std::vector < double > A_co = coords_q ( A );  // same as coords_Eu ( A )
	std::vector < double > B_co = coords_q ( B, tag::winding, s );
	size_t n = A_co .size();
	assert ( n == B_co .size() );
	double len_AB_2 = 0.;
	for ( size_t i = 0; i < n; i++ )
		{	double v = B_co [i] - A_co [i];  len_AB_2 += v*v;  }
	return len_AB_2;                                              }
	

class compare_lenghts_of_segs

{ public :
	inline bool operator() ( Cell AB, Cell CD ) const
	{	return length_square ( AB ) > length_square ( CD );  }
};

//-----------------------------------------------------------------------------------//


inline bool split_segment ( Mesh & msh, Cell & seg )

// splits 'seg' in two if it is inner to 'msh'
// splits also the two neighbour triangles
// equilibrates (baricenter) the four neighbour vertices
// among the four segments "around" the newly created vertex,
// flips the longest

// returns true if the segment has been split, false if not
	
// assumes there is no higher-dimensional mesh "above" 'msh'

// assumes also that the current manifold is not implicit
// (for an implicit manifold, a projection operation should be added)

// assumes that the current manifold is a quotient manifold (manipulates windings)

{	Manifold space = Manifold::working;
	assert ( space .exists() );  // we use the current (quotient) manifold
	Manifold::Quotient * mani_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space .core );
	Function coords_q = space .coordinates();
	Manifold mani_Eu = mani_q->base_space;  // underlying Euclidian manifold
	Function coords_Eu = mani_Eu .coordinates();

	if ( not seg .belongs_to ( msh ) ) return false;

	Cell tri2 = msh .cell_in_front_of ( seg, tag::may_not_exist );
	if ( not tri2 .exists() ) return false;
	Cell tri1 = msh .cell_behind ( seg, tag::may_not_exist );
	if ( not tri1 .exists() ) return false;
	// or, equivalently :  if ( not seg .is_inner_to ( msh ) ) return false
	// refuse to split segments on the boundary

	Cell A = seg .base() .reverse();
	Cell B = seg .tip();
	Manifold::Action s = seg .winding();
	std::vector < double > A_co = coords_q ( A );  // same as coords_Eu ( A )
	std::vector < double > B_co = coords_q ( B, tag::winding, s );
	size_t n = A_co .size();

	Cell BC = tri1 .boundary() .cell_in_front_of ( B, tag::surely_exists );
	Cell CA = tri1 .boundary() .cell_behind ( A, tag::surely_exists );
	Cell AD = tri2 .boundary() .cell_in_front_of ( A, tag::surely_exists );
	Cell DB = tri2 .boundary() .cell_behind ( B, tag::surely_exists );
	Cell C = BC .tip();
	assert ( CA .base() .reverse() == C );
	Cell D = AD .tip();
	assert ( DB .base() .reverse() == D );

	assert ( seg .winding() == - BC .winding() - CA .winding() );
	assert ( seg .winding() == AD .winding() + DB .winding() );

	Cell E ( tag::vertex );   // put 'E' at the middle of 'seg' (aka AB)
	for ( size_t i = 0; i < n; i++ )
		coords_Eu [i] ( E ) = ( A_co [i] + B_co [i] ) / 2.;

	A .reverse() .cut_from_bdry_of ( seg, tag::do_not_bother );
	E .reverse() .glue_on_bdry_of ( seg, tag::do_not_bother );
	CA .cut_from_bdry_of ( tri1, tag::do_not_bother );
	AD .cut_from_bdry_of ( tri2, tag::do_not_bother );

	Cell AE ( tag::segment, A .reverse(), E );  // no winding
	Cell CE ( tag::segment, C .reverse(), E );
	CE .winding() = CA .winding();
	Cell DE ( tag::segment, D .reverse(), E );
	DE .winding() = - AD .winding();

	CE .glue_on_bdry_of ( tri1, tag::do_not_bother );
	DE .reverse() .glue_on_bdry_of ( tri2, tag::do_not_bother );

	tri1 .boundary() .closed_loop ( B );
	tri2 .boundary() .closed_loop ( B );

	Cell AEC ( tag::triangle, AE, CE .reverse(), CA );
	AEC .add_to_mesh ( msh );
	Cell ADE ( tag::triangle, AD, DE, AE .reverse() );
	ADE .add_to_mesh ( msh );

	// we want to sweep over the four segments "around" E : CA, AD, DB, BC
	// we want to flip the longest one (could be the two longest)
	// if the longest one is stuck (is on the boundary) we flip the next one

	// we use a map to order them
	compare_lenghts_of_segs comp_len;
	std::multiset < Cell, compare_lenghts_of_segs > ms ( comp_len );
	ms .insert ( CA );
	ms .insert ( AD );
	ms .insert ( DB );
	ms .insert ( BC );

	for ( std::multiset < Cell, compare_lenghts_of_segs > ::iterator it = ms .begin();
				it != ms .end(); it ++                                                      )
	{	Cell sseg = * it;
		if ( flip_segment ( msh, sseg ) )  break;  }
		// we choose to to flip only one segment

	if ( A .is_inner_to ( msh ) ) msh .baricenter ( A, tag::winding );
	if ( B .is_inner_to ( msh ) ) msh .baricenter ( B, tag::winding );
	if ( C .is_inner_to ( msh ) ) msh .baricenter ( C, tag::winding );
	if ( D .is_inner_to ( msh ) ) msh .baricenter ( D, tag::winding );
	if ( E .is_inner_to ( msh ) ) msh .baricenter ( E, tag::winding );

	return true;

}  // end of  split_segment


//-----------------------------------------------------------------------------------------

void limit_number_of_neighbours ( Mesh msh )

// only applies to meshes of triangular cells
	
// assumes there is no higher-dimensional mesh "above" 'msh'

// calls 'flip_segment' which assumes that the current manifold is a quotient manifold
// (manipulates windings)

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
			assert ( tri .dim() == 2 );
			Cell SP = tri .boundary() .cell_behind ( P, tag::surely_exists );
			Cell PT = tri .boundary() .cell_in_front_of ( P, tag::surely_exists );
			Cell T = PT .tip();
			Cell TS = tri .boundary() .cell_in_front_of ( T, tag::surely_exists );
			assert ( TS .tip() == SP .base() .reverse() );
			ms .insert ( TS );                                                      }
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
			assert ( seg .tip() == P );
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

	
void flip_split_long_segments ( Mesh & msh, double threshold )

// flip long segments, make baricenters on the four neighbour vertices
// if a segment is still long after flip, split it

// assumes there is no higher-dimensional mesh "above" 'msh'

// assumes also that the current manifold is not implicit
// (for an implicit manifold, a projection operation should be added)

// calls 'flip_segment' and 'split_segment' which assume that the current manifold
// is a quotient manifold (manipulate windings)

{
	double thr_sq = threshold * threshold;

	std::list < Cell > list_of_segments;

	Mesh::Iterator it = msh .iterator ( tag::over_segments );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell seg = *it;

	  Cell tri1 = msh .cell_in_front_of ( seg, tag::may_not_exist );
		if ( not tri1 .exists() ) continue;
	  Cell tri2 = msh .cell_behind ( seg, tag::may_not_exist );
		if ( not tri2 .exists() ) continue;
		// or, equivalently :  if ( not seg .is_inner_to ( msh ) ) continue
		// segments on the boundary cannot be flipped

		if ( length_square ( seg ) > thr_sq )  list_of_segments .push_back ( seg );   }

	for ( size_t ii = 0; ii < 2; ii++ )
	for ( std::list < Cell > ::iterator itt = list_of_segments .begin();
        itt != list_of_segments .end();                                )
	{	Cell seg = * itt;

		// 'seg' may have been eliminated from 'msh' in the meanwhile
		if ( not seg .belongs_to ( msh ) )
			{	itt = list_of_segments .erase ( itt );  continue;  }

		// the configuration around 'seg' may have changed in the meanwhile
		Cell tri2 = msh .cell_in_front_of ( seg, tag::may_not_exist );
		if ( not tri2 .exists() ) 
		{	itt = list_of_segments .erase ( itt );  continue;  }
		Cell tri1 = msh .cell_behind ( seg, tag::may_not_exist );
		if ( not tri1 .exists() )
		{	itt = list_of_segments .erase ( itt );  continue;  }
		// or, equivalently :  if ( not seg.is_inner_to ( msh ) ) continue
		// segments on the boundary cannot be flipped

		// length of 'seg' may have changed in the meanwhile
		if ( length_square ( seg ) < thr_sq ) 
		{	itt = list_of_segments .erase ( itt );  continue;  }

		flip_segment ( msh, seg );
		itt ++;                                                                   }

	// segments still on the list will be split

	for ( std::list < Cell > ::iterator itt = list_of_segments .begin();
        itt != list_of_segments .end(); itt ++                         )
	{	Cell seg = * itt;
		// the configuration around 'seg' may have changed in the meanwhile
	  Cell tri1 = msh .cell_in_front_of ( seg, tag::may_not_exist );
		if ( not tri1 .exists() ) continue;
	  Cell tri2 = msh .cell_behind ( seg, tag::may_not_exist );
		if ( not tri2 .exists() ) continue;
		// or, equivalently :  if ( not seg .is_inner_to ( msh ) ) continue
		if ( length_square ( seg ) > thr_sq )  split_segment ( msh, seg );  }

}

//-----------------------------------------------------------------------------------------


void remove_short_segments ( Mesh & msh, double threshold )

// remove segments shorter than the given threshold

// assumes there is no higher-dimensional mesh "above" 'msh'

// assumes also that the current manifold is not implicit
// (for an implicit manifold, a projection operation should be added)

// assumes that the current manifold is a quotient manifold
// (it manipulates windings)

{	Manifold space = Manifold::working;
	assert ( space.exists() );  // we use the current (quotient) manifold
	Manifold::Quotient * mani_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space .core );
	Function coords_q = space .coordinates();
	Manifold mani_Eu = mani_q->base_space;  // underlying Euclidian manifold
	Function coords_Eu = mani_Eu .coordinates();

	double thr_sq = threshold * threshold;

	std::forward_list < Cell > list_of_segments;

	Mesh::Iterator it = msh .iterator ( tag::over_segments );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell seg = *it;

	  Cell tri1 = msh .cell_in_front_of ( seg, tag::may_not_exist );
		if ( not tri1 .exists() ) continue;
	  Cell tri2 = msh .cell_behind ( seg, tag::may_not_exist );
		if ( not tri2 .exists() ) continue;
		// or, equivalently :  if ( not seg .is_inner_to ( msh ) ) continue
		// segments on the boundary are left unchanged

		Cell A = seg .base() .reverse();
		Cell B = seg .tip();
		if ( length_square ( seg ) > thr_sq ) continue;
		// long segments are left unchanged

		// method 'is_inner_to' says "vertex is not on the boundary of mesh"
		bool wing_A_free = A .is_inner_to (msh );
		bool wing_B_free = B .is_inner_to (msh );

		if ( wing_A_free or wing_B_free ) list_of_segments .push_front ( seg );  }

	for ( std::forward_list < Cell > ::iterator	itt = list_of_segments .begin();
	      itt != list_of_segments .end(); itt ++                                 )
	{	Cell seg = *itt;
		assert ( seg .exists() );

		// 'seg' may have been already eliminated
		if ( not seg .belongs_to ( msh ) ) continue;
		
	  Cell tri1 = msh .cell_in_front_of ( seg, tag::may_not_exist );
		if ( not tri1 .exists() ) continue;
	  Cell tri2 = msh .cell_behind ( seg, tag::may_not_exist );
		if ( not tri2 .exists() ) continue;
		// or, equivalently :  if ( not seg.is_inner_to ( msh ) ) continue
		// segments on the boundary are left unchanged

		Cell A = seg .base() .reverse();
		Cell B = seg .tip();
		Manifold::Action s = seg .winding();
		std::vector < double > A_co = coords_q ( A );
		std::vector < double > B_co = coords_q ( B, tag::winding, s );
		size_t n = A_co.size();

		Cell CB = tri1 .boundary() .cell_behind ( B, tag::surely_exists );
		Cell AC = tri1 .boundary() .cell_in_front_of ( A, tag::surely_exists );
		Cell DA = tri2 .boundary() .cell_behind ( A, tag::surely_exists );
		Cell BD = tri2 .boundary() .cell_in_front_of ( B, tag::surely_exists );
		assert ( CB .base() .reverse() == AC .tip() );
		assert ( DA .base() .reverse() == BD .tip() );
		Cell tri3 = msh .cell_in_front_of ( CB, tag::may_not_exist );
		Cell tri4 = msh .cell_in_front_of ( BD, tag::may_not_exist );
		Cell tri5 = msh .cell_in_front_of ( AC, tag::may_not_exist );
		Cell tri6 = msh .cell_in_front_of ( DA, tag::may_not_exist );

		// method 'is_inner_to' says "vertex is not on the boundary of mesh"
		bool wing_A_free = A .is_inner_to (msh );
		bool wing_B_free = B .is_inner_to (msh );

		// due to previous changes, wings may be not free any longer
		if ( ( not wing_A_free ) and ( not wing_B_free ) ) continue;

		// beware, if msh happens to be the boundary of some higher-dimensional cell,
		// and that higher-dim cell has neighbours in some higher-dim mesh,
		// statements below do not work as expected

		if ( wing_B_free and wing_A_free )
			// both wings free, we collapse one of them at random
			// we collapse wing_B, vertex B disappears,
			// A will occupy middle of segment seg

		{	// make a list of segments pointing towards B
			std::forward_list < Cell > list_of_segs;
			Mesh::Iterator it_around_B = msh .iterator ( tag::over_segments, tag::around, B );
			it_around_B .reset ( tag::start_at, BD .reverse() );
			assert ( it_around_B .in_range() );
			assert ( *it_around_B == BD .reverse() );
			it_around_B++;
			assert ( it_around_B .in_range() );
			assert ( *it_around_B == seg );
			it_around_B++;
			assert ( it_around_B .in_range() );
			assert ( *it_around_B == CB );
			// seg, BD and CB will disappear
			for ( it_around_B++; it_around_B .in_range(); it_around_B++ )
				list_of_segs .push_front ( *it_around_B );
			
			// change winding number of these segments
			if ( s != 0 )
				for ( std::forward_list < Cell > ::iterator it_list = list_of_segs .begin();
			        it_list != list_of_segs .end(); it_list ++                            )
				{	Cell segm = *it_list;      // 'segm' points towards B
				  segm.winding() -= s;     }

			tri1 .remove_from_mesh ( msh );
			tri2 .remove_from_mesh ( msh );
			assert ( tri3 .exists() );
			assert ( tri4 .exists() );
			assert ( tri5 .exists() );
			assert ( tri6 .exists() );
			CB .reverse() .cut_from_bdry_of ( tri3, tag::do_not_bother );
			BD .reverse() .cut_from_bdry_of ( tri4, tag::do_not_bother );
			AC .glue_on_bdry_of ( tri3, tag::do_not_bother );
			DA .glue_on_bdry_of ( tri4, tag::do_not_bother );

			// replace B by A in segments neighbour to B
			for ( std::forward_list < Cell > ::iterator it_list = list_of_segs .begin();
			      it_list != list_of_segs .end(); it_list ++                            )
			{	Cell segm = *it_list;
				assert ( segm.tip() == B );
				B .cut_from_bdry_of ( segm, tag::do_not_bother );
				A .glue_on_bdry_of ( segm, tag::do_not_bother );  }

			tri3 .boundary() .closed_loop ( A );
			tri4 .boundary() .closed_loop ( A );

			// change coordinates of A
			// if we are on an implicit manifold, we should project A
			for ( size_t i = 0; i < n; i++ ) A_co [i] = ( A_co [i] + B_co [i] ) / 2.;
			coords_Eu ( A ) = A_co;                                                        }

		else if ( wing_B_free )
			// we collapse wing_B, vertex B disappears, vertex A stays where it is

		{	// make a list of segments pointing towards B
			std::forward_list < Cell > list_of_segs;
			Mesh::Iterator it_around_B = msh .iterator ( tag::over_segments, tag::around, B );
			it_around_B .reset ( tag::start_at, BD .reverse() );
			assert ( it_around_B .in_range() );
			assert ( *it_around_B == BD .reverse() );
			it_around_B++;
			assert ( it_around_B .in_range() );
			assert ( *it_around_B == seg );
			it_around_B++;
			assert ( it_around_B .in_range() );
			assert ( *it_around_B == CB );
			// seg, BD and CB will disappear
			for ( it_around_B++; it_around_B .in_range(); it_around_B++ )
				list_of_segs .push_front ( *it_around_B );
			
			// change winding number of these segments
			if ( s != 0 )
				for ( std::forward_list < Cell > ::iterator it_list = list_of_segs .begin();
			        it_list != list_of_segs .end(); it_list ++                            )
				{	Cell segm = *it_list;      // 'segm' points towards B
				  segm .winding() -= s;     }

			tri1 .remove_from_mesh ( msh );
			tri2 .remove_from_mesh ( msh );
			assert ( tri3 .exists() );
			assert ( tri4 .exists() );
			CB .reverse() .cut_from_bdry_of ( tri3, tag::do_not_bother );
			BD .reverse() .cut_from_bdry_of ( tri4, tag::do_not_bother );
			AC .glue_on_bdry_of ( tri3, tag::do_not_bother );
			DA .glue_on_bdry_of ( tri4, tag::do_not_bother );

			// replace B by A in segments neighbour to B
			for ( std::forward_list < Cell > ::iterator it_list = list_of_segs .begin();
			      it_list != list_of_segs .end(); it_list ++                            )
			{	Cell segm = *it_list;
				assert ( segm.tip() == B );
				assert ( B.exists() );  assert ( segm.exists() );
				B .cut_from_bdry_of ( segm, tag::do_not_bother );
				A .glue_on_bdry_of ( segm, tag::do_not_bother );  }

			tri3 .boundary() .closed_loop ( A );
			tri4 .boundary() .closed_loop ( A );                                           }
		
		else if ( wing_A_free )
			// we collapse wing_A, vertex A disappears,, vertex B stays where it is

		{	// make a list of segments pointing towards A
			std::forward_list < Cell > list_of_segs;
			Mesh::Iterator it_around_A = msh .iterator ( tag::over_segments, tag::around, A );
			it_around_A .reset ( tag::start_at, AC .reverse() );
			assert ( it_around_A .in_range() );
			assert ( *it_around_A == AC .reverse() );
			it_around_A++;
			assert ( it_around_A .in_range() );
			assert ( *it_around_A == seg .reverse() );
			it_around_A++;
			assert ( it_around_A .in_range() );
			assert ( *it_around_A == DA );
			// seg, AC and DA will disappear
			for ( it_around_A++; it_around_A .in_range(); it_around_A++ )
				list_of_segs .push_front ( *it_around_A );
			
			// change winding number of these segments
			if ( s != 0 )
				for ( std::forward_list < Cell > ::iterator it_list = list_of_segs .begin();
			        it_list != list_of_segs .end(); it_list ++                            )
				{	Cell segm = *it_list;      // 'segm' points towards A
				  segm .winding() += s;     }

			tri1 .remove_from_mesh ( msh );
			tri2 .remove_from_mesh ( msh );
			assert ( tri5 .exists() );
			assert ( tri6 .exists() );
			AC .reverse() .cut_from_bdry_of ( tri5, tag::do_not_bother );
			DA .reverse() .cut_from_bdry_of ( tri6, tag::do_not_bother );
			CB .glue_on_bdry_of ( tri5, tag::do_not_bother );
			BD .glue_on_bdry_of ( tri6, tag::do_not_bother );

			// replace A by B in segments neighbour to A
			for ( std::forward_list < Cell > ::iterator it_list = list_of_segs .begin();
			      it_list != list_of_segs .end(); it_list ++                            )
			{	Cell segm = *it_list;
				assert ( segm .tip() == A );
				A .cut_from_bdry_of ( segm, tag::do_not_bother );
				B .glue_on_bdry_of ( segm, tag::do_not_bother );  }

			tri5 .boundary() .closed_loop ( B );
			tri6 .boundary() .closed_loop ( B );                                          }

		else assert ( false );    // at least one wing should be free
		
	}  // end of  for  over segments of msh
	
}  // end of  remove_short_segments
	

void baricenters ( Mesh & msh )

// segments on the boundary are not moved -- we use method  is_inner_to

{	Mesh::Iterator it = msh .iterator ( tag::over_vertices );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell P = *it;
		if ( P .is_inner_to ( msh ) ) msh .baricenter ( P, tag::winding );  }  }
	
