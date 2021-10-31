
// example presented in paragraph 7.17 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// build an intricate mesh on a torus, having two V-shaped holes
// homogenized elastic tensor exhibits negative Poisson coefficient

// progressive mesh generation does not work yet on quotient manifolds
// we fake the result by folding a mesh whose exterior boundary is a polygonal line
// see https://github.com/cristian-barbarosie/attic/blob/main/five-V-holes.png

#include "maniFEM.h"

using namespace maniFEM;

	
void remove_short_segments ( Mesh & msh, double threshold );

void flip_long_segments ( Mesh & msh, double threshold );

void baricenters ( Mesh & msh );


int main ( )

{	// begin with the usual two-dimensional space
	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
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

	Mesh AB ( tag::segment, A.reverse(), B, tag::divided_in, 14 );
	Mesh BC ( tag::segment, B.reverse(), C, tag::divided_in, 32 );
	Mesh CD ( tag::segment, C.reverse(), D, tag::divided_in, 32 );
	Mesh DE ( tag::segment, D.reverse(), E, tag::divided_in, 14 );
	Mesh EF ( tag::segment, E.reverse(), F, tag::divided_in, 32 );
	Mesh FG ( tag::segment, F.reverse(), G, tag::divided_in, 14 );
	Mesh GH ( tag::segment, G.reverse(), H, tag::divided_in, 32 );
	Mesh HI ( tag::segment, H.reverse(), I, tag::divided_in, 32 );
	Mesh IJ ( tag::segment, I.reverse(), J, tag::divided_in, 14 );
	Mesh JA ( tag::segment, J.reverse(), A, tag::divided_in, 32 );

	Mesh zig_zag ( tag::join, { AB, BC, CD, DE, EF, FG, GH, HI, IJ, JA } );

	const double e = 1.5;
	Manifold curve1 = RR2.implicit 
		( smooth_min ( 300.*power((x+y)*(x+y),e) + power((x-y-1.)*(x-y-1.),e),
		               300.*power((x-y)*(x-y),e) + power((x+y+1.)*(x+y+1.),e),
	                 tag::threshold, 1.5                     )  == 1. );
	Mesh loop1 ( tag::progressive, tag::desired_length, 0.05 );

	double a = 0.15, b = 0.82;
	Manifold curve2 = RR2.implicit 
		( smooth_min ( 300.*power((x-a+y+b)*(x-a+y+b),e) + power((x-a-y-b-1.)*(x-a-y-b-1.),e),
		               300.*power((x-2.-a-y-b)*(x-2.-a-y-b),e) + power((x-2.-a+y+b+1.)*(x-2.-a+y+b+1.),e),
	                 tag::threshold, 1.5                     )  == 1. );
	Mesh loop2 ( tag::progressive, tag::desired_length, 0.05 );

	b = 1.4;
	Manifold curve3 = RR2.implicit 
		( smooth_min ( 300.*power((x+y+b)*(x+y+b),e) + power((x-y-b-1.)*(x-y-b-1.),e),
		               300.*power((x-y-b)*(x-y-b),e) + power((x+y+b+1.)*(x+y+b+1.),e),
	                 tag::threshold, 1.5                     )  == 1. );
	Mesh loop3 ( tag::progressive, tag::desired_length, 0.05 );  // not used

	b = 0.82-1.4;
	Manifold curve4 = RR2.implicit 
		( smooth_min ( 300.*power((x-a+y+b)*(x-a+y+b),e) + power((x-a-y-b-1.)*(x-a-y-b-1.),e),
		               300.*power((x-2.-a-y-b)*(x-2.-a-y-b),e) + power((x-2.-a+y+b+1.)*(x-2.-a+y+b+1.),e),
	                 tag::threshold, 1.5                     )  == 1. );
	Mesh loop4 ( tag::progressive, tag::desired_length, 0.05 );  // not used

	a -= 2.3;  b = 0.82;
	Manifold curve5 = RR2.implicit 
		( smooth_min ( 300.*power((x-a+y+b)*(x-a+y+b),e) + power((x-a-y-b-1.)*(x-a-y-b-1.),e),
		               300.*power((x-2.-a-y-b)*(x-2.-a-y-b),e) + power((x-2.-a+y+b+1.)*(x-2.-a+y+b+1.),e),
	                 tag::threshold, 1.5                     )  == 1. );
	Mesh loop5 ( tag::progressive, tag::desired_length, 0.05 );  // not used

	Mesh loop ( tag::join, loop1.reverse(), loop2.reverse(), zig_zag );

	RR2.set_as_working_manifold();
	Mesh domain ( tag::progressive, tag::boundary, loop, tag::desired_length, 0.05 );

	Mesh IJAB ( tag::join, IJ, JA, AB );
	Mesh GFED ( tag::join, DE.reverse(), EF.reverse(), FG.reverse() );

	Mesh VV = domain.fold ( tag::identify, IJAB, tag::with, GFED,
                          tag::identify, BC, tag::with, GH.reverse(),
													tag::identify, CD, tag::with, HI.reverse(),
													tag::use_existing_vertices                 );

	// std::cout << "produced folded mesh, now drawing, please wait" << std::endl << std::flush;
	// VV.draw_ps ( "VV.eps", tag::unfold, tag::over_region, -1.5 < x < 2.5, -2 < y < 0.5 );
	
	std::cout << "now smoothening ... " << std::flush;

	double d = 0.05;
	
	baricenters ( VV );
	// re-define the position of each vertex as baricenter of its neighbours
	
	flip_long_segments ( VV, 1.5*d );
	// remove segments shorter than the given threshold
	
	remove_short_segments ( VV, 0.5*d );
	// remove segments shorter than the given threshold
	
	baricenters ( VV );
	// re-define the position of each vertex as baricenter of its neighbours
	
	flip_long_segments ( VV, 1.5*d );
	remove_short_segments ( VV, 0.5*d );
	baricenters ( VV );
	flip_long_segments ( VV, 1.5*d );
	remove_short_segments ( VV, 0.5*d );
	baricenters ( VV );
	flip_long_segments ( VV, 1.45*d );
	remove_short_segments ( VV, 0.55*d );
	baricenters ( VV );
	flip_long_segments ( VV, 1.45*d );
	remove_short_segments ( VV, 0.55*d );
	baricenters ( VV );
	flip_long_segments ( VV, 1.45*d );
	remove_short_segments ( VV, 0.55*d );
	baricenters ( VV );
	flip_long_segments ( VV, 1.4*d );
	remove_short_segments ( VV, 0.6*d );
	baricenters ( VV );
	flip_long_segments ( VV, 1.4*d );
	remove_short_segments ( VV, 0.6*d );
	baricenters ( VV );
	flip_long_segments ( VV, 1.4*d );
	remove_short_segments ( VV, 0.6*d );
	baricenters ( VV );
	flip_long_segments ( VV, 1.35*d );
	remove_short_segments ( VV, 0.65*d );
	baricenters ( VV );
	flip_long_segments ( VV, 1.35*d );
	remove_short_segments ( VV, 0.65*d );
	baricenters ( VV );
	flip_long_segments ( VV, 1.35*d );
	remove_short_segments ( VV, 0.65*d );
	baricenters ( VV );

	std::cout << "and drawing again, please wait" << std::endl << std::flush;

	VV.draw_ps ( "VV-smooth.eps", tag::unfold,
               tag::over_region, -1.5 < x < 2.5, -2 < y < 0.5 );

	std::cout << "produced files VV.eps and VV-smooth.eps - please edit before viewing" << std::endl;	
}



//-----------------------------------------------------------------------------------------


void flip_split_long_segments ( Mesh & msh, double threshold )

// flip long segments, make baricenters on the four neighbour vertices
// if a segment is still long after flip, split it

// this function assumes there is no higher-dimensional mesh "above" 'msh'

// assumes also that the current manifold is not implicit
// (for an implicit manifold, a projection operation should be added)

// assumes that the current manifold is a quotient manifold
// (it manipulates spins)

{	Manifold space = Manifold::working;
	assert ( space.exists() );  // we use the current (quotient) manifold
	Manifold::Quotient * mani_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space.core );
	Function coords_q = space.coordinates();
	Manifold mani_Eu = mani_q->base_space;  // underlying Euclidian manifold
	Function coords_Eu = mani_Eu.coordinates();

	double thr_sq = threshold * threshold;

	std::list < Cell > list_of_segments;

	CellIterator it = msh.iterator ( tag::over_segments );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell seg = *it;

	  Cell tri1 = msh.cell_in_front_of ( seg, tag::may_not_exist );
		if ( not tri1.exists() ) continue;
	  Cell tri2 = msh.cell_behind ( seg, tag::may_not_exist );
		if ( not tri2.exists() ) continue;
		// or, equivalently :  if ( not seg.is_inner_to ( msh ) ) continue
		// segments on the boundary cannot be flipped

		Cell A = seg.base().reverse();
		Cell B = seg.tip();
		Function::Action s = seg.spin();
		std::vector < double > A_co = coords_q ( A );
		std::vector < double > B_co = coords_q ( B, tag::spin, s );
		size_t n = A_co.size();
		double len_sq = 0.;
		for ( size_t i = 0; i < n; i++ )
		{	double d = B_co[i] - A_co[i];
			len_sq += d*d;                }
		if ( len_sq > thr_sq )  list_of_segments .push_back ( seg );   }

	for ( size_t i = 0; i < 2; i++ )
	for ( std::list < Cell > ::iterator itt = list_of_segments .begin();
        itt != list_of_segments .end();                                )
	{	Cell seg = * itt;

		// 'seg' may have been eliminated from 'msh' in the meanwhile
		if ( not seg .belongs_to ( msh ) )
			{	itt = list_of_segments .erase ( itt );  continue;  }

		// the configuration around 'seg' may have changed in the meanwhile
		Cell tri2 = msh.cell_in_front_of ( seg, tag::may_not_exist );
		if ( not tri2.exists() ) 
			{	itt = list_of_segments .erase ( itt );  continue;  }
		Cell tri1 = msh.cell_behind ( seg, tag::may_not_exist );
		if ( not tri1.exists() )
			{	itt = list_of_segments .erase ( itt );  continue;  }
		// or, equivalently :  if ( not seg.is_inner_to ( msh ) ) continue
		// segments on the boundary cannot be flipped

		// length of 'seg' may have changed in the meanwhile
		Cell A = seg.base().reverse();
		Cell B = seg.tip();
		Function::Action s = seg.spin();
		std::vector < double > A_co = coords_q ( A );
		std::vector < double > B_co = coords_q ( B, tag::spin, s );
		size_t n = A_co.size();
		double len_sq = 0.;
		for ( size_t i = 0; i < n; i++ )
			{	double d = B_co[i] - A_co[i];
				len_sq += d*d;                }
		if ( len_sq < thr_sq ) 
			{	itt = list_of_segments .erase ( itt );  continue;  }

		// flip
		Cell BC = tri1 .boundary() .cell_in_front_of ( B, tag::surely_exists );
		Cell CA = tri1 .boundary() .cell_behind ( A, tag::surely_exists );
		Cell AD = tri2 .boundary() .cell_in_front_of ( A, tag::surely_exists );
		Cell DB = tri2 .boundary() .cell_behind ( B, tag::surely_exists );
		Cell C = BC.tip();
		assert ( CA.base().reverse() == C );
		Cell D = AD.tip();
		assert ( DB.base().reverse() == D );

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

		if ( A .is_inner_to ( msh ) ) msh .baricenter ( A, tag::spin );
		if ( B .is_inner_to ( msh ) ) msh .baricenter ( B, tag::spin );
		if ( C .is_inner_to ( msh ) ) msh .baricenter ( C, tag::spin );
		if ( D .is_inner_to ( msh ) ) msh .baricenter ( D, tag::spin );
			
		itt ++;                                                                   }

	// segments still on the list will be split
			
	for ( std::list < Cell > ::iterator itt = list_of_segments .begin();
        itt != list_of_segments .end();                                )
	{	Cell seg = * itt;

	  Cell tri1 = msh.cell_in_front_of ( seg, tag::may_not_exist );
		assert ( tri1.exists() );
	  Cell tri2 = msh.cell_behind ( seg, tag::may_not_exist );
		assert ( tri2.exists() );
		// or, equivalently :  assert ( not seg.is_inner_to ( msh ) )

		Cell A = seg.base().reverse();
		Cell B = seg.tip();
		Function::Action s = seg.spin();
		std::vector < double > A_co = coords_q ( A );
		std::vector < double > B_co = coords_q ( B, tag::spin, s );
		size_t n = A_co.size();

		
}

//-----------------------------------------------------------------------------------------


void remove_short_segments ( Mesh & msh, double threshold )

// this function assumes there is no higher-dimensional mesh "above" 'msh'

// assumes also that the current manifold is not implicit
// (for an implicit manifold, a projection operation should be added)

// assumes that the current manifold is a quotient manifold
// (it manipulates spins)

{	Manifold space = Manifold::working;
	assert ( space.exists() );  // we use the current (quotient) manifold
	Manifold::Quotient * mani_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space.core );
	Function coords_q = space.coordinates();
	Manifold mani_Eu = mani_q->base_space;  // underlying Euclidian manifold
	Function coords_Eu = mani_Eu.coordinates();

	double thr_sq = threshold * threshold;

	std::forward_list < Cell > list_of_segments;

	CellIterator it = msh.iterator ( tag::over_segments );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell seg = *it;

	  Cell tri1 = msh.cell_in_front_of ( seg, tag::may_not_exist );
		if ( not tri1.exists() ) continue;
	  Cell tri2 = msh.cell_behind ( seg, tag::may_not_exist );
		if ( not tri2.exists() ) continue;
		// or, equivalently :  if ( not seg.is_inner_to ( msh ) ) continue
		// segments on the boundary are left unchanged

		Cell A = seg.base().reverse();
		Cell B = seg.tip();
		Function::Action s = seg.spin();
		std::vector < double > A_co = coords_q ( A );
		std::vector < double > B_co = coords_q ( B, tag::spin, s );
		size_t n = A_co.size();
		double len_sq = 0.;
		for ( size_t i = 0; i < n; i++ )
		{	double d = B_co[i] - A_co[i];
			len_sq += d*d;                }
		if ( len_sq > thr_sq ) continue;
		// long segments are left unchanged

		// method 'is_inner_to' says "vertex is not on the boundary of mesh"
		bool wing_A_free = A .is_inner_to (msh );
		bool wing_B_free = B .is_inner_to (msh );

		if ( wing_A_free or wing_B_free ) list_of_segments .push_front ( seg );  }

	for ( std::forward_list < Cell > ::iterator	itt = list_of_segments .begin();
	      itt != list_of_segments .end(); itt ++                                 )
	{	Cell seg = *itt;
		assert ( seg.exists() );

		// 'seg' may have been already eliminated
		if ( not seg .belongs_to ( msh ) ) continue;
		
	  Cell tri1 = msh.cell_in_front_of ( seg, tag::may_not_exist );
		if ( not tri1.exists() ) continue;
	  Cell tri2 = msh.cell_behind ( seg, tag::may_not_exist );
		if ( not tri2.exists() ) continue;
		// or, equivalently :  if ( not seg.is_inner_to ( msh ) ) continue
		// segments on the boundary are left unchanged

		Cell A = seg.base().reverse();
		Cell B = seg.tip();
		Function::Action s = seg.spin();
		std::vector < double > A_co = coords_q ( A );
		std::vector < double > B_co = coords_q ( B, tag::spin, s );
		size_t n = A_co.size();

		Cell CB = tri1 .boundary() .cell_behind ( B, tag::surely_exists );
		Cell AC = tri1 .boundary() .cell_in_front_of ( A, tag::surely_exists );
		Cell DA = tri2 .boundary() .cell_behind ( A, tag::surely_exists );
		Cell BD = tri2 .boundary() .cell_in_front_of ( B, tag::surely_exists );
		assert ( CB.base().reverse() == AC.tip() );
		assert ( DA.base().reverse() == BD.tip() );
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
			CellIterator it_around_B = msh.iterator ( tag::over_segments, tag::around, B );
			it_around_B .reset ( tag::start_at, BD.reverse() );
			assert ( it_around_B .in_range() );
			assert ( *it_around_B == BD.reverse() );
			it_around_B++;
			assert ( it_around_B .in_range() );
			assert ( *it_around_B == seg );
			it_around_B++;
			assert ( it_around_B .in_range() );
			assert ( *it_around_B == CB );
			// seg, BD and CB will disappear
			for ( it_around_B++; it_around_B .in_range(); it_around_B++ )
				list_of_segs .push_front ( *it_around_B );
			
			// change spin of these segments
			if ( s != 0 )
				for ( std::forward_list < Cell > ::iterator it_list = list_of_segs .begin();
			        it_list != list_of_segs .end(); it_list ++                            )
				{	Cell segm = *it_list;      // 'segm' points towards B
				  segm.spin() -= s;     }

			tri1 .remove_from_mesh ( msh );
			tri2 .remove_from_mesh ( msh );
			assert ( tri3.exists() );
			assert ( tri4.exists() );
			assert ( tri5.exists() );
			assert ( tri6.exists() );
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
			for ( size_t i = 0; i < n; i++ ) A_co[i] = ( A_co[i] + B_co[i] ) / 2.;
			coords_Eu ( A ) = A_co;                                                        }

		else if ( wing_B_free )
			// we collapse wing_B, vertex B disappears, vertex A stays where it is

		{	// make a list of segments pointing towards B
			std::forward_list < Cell > list_of_segs;
			CellIterator it_around_B = msh.iterator ( tag::over_segments, tag::around, B );
			it_around_B .reset ( tag::start_at, BD.reverse() );
			assert ( it_around_B .in_range() );
			assert ( *it_around_B == BD.reverse() );
			it_around_B++;
			assert ( it_around_B .in_range() );
			assert ( *it_around_B == seg );
			it_around_B++;
			assert ( it_around_B .in_range() );
			assert ( *it_around_B == CB );
			// seg, BD and CB will disappear
			for ( it_around_B++; it_around_B .in_range(); it_around_B++ )
				list_of_segs .push_front ( *it_around_B );
			
			// change spin of these segments
			if ( s != 0 )
				for ( std::forward_list < Cell > ::iterator it_list = list_of_segs .begin();
			        it_list != list_of_segs .end(); it_list ++                            )
				{	Cell segm = *it_list;      // 'segm' points towards B
				  segm.spin() -= s;     }

			tri1 .remove_from_mesh ( msh );
			tri2 .remove_from_mesh ( msh );
			assert ( tri3.exists() );
			assert ( tri4.exists() );
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
			CellIterator it_around_A = msh.iterator ( tag::over_segments, tag::around, A );
			it_around_A .reset ( tag::start_at, AC.reverse() );
			assert ( it_around_A .in_range() );
			assert ( *it_around_A == AC.reverse() );
			it_around_A++;
			assert ( it_around_A .in_range() );
			assert ( *it_around_A == seg.reverse() );
			it_around_A++;
			assert ( it_around_A .in_range() );
			assert ( *it_around_A == DA );
			// seg, AC and DA will disappear
			for ( it_around_A++; it_around_A .in_range(); it_around_A++ )
				list_of_segs .push_front ( *it_around_A );
			
			// change spin of these segments
			if ( s != 0 )
				for ( std::forward_list < Cell > ::iterator it_list = list_of_segs .begin();
			        it_list != list_of_segs .end(); it_list ++                            )
				{	Cell segm = *it_list;      // 'segm' points towards A
				  segm.spin() += s;     }

			tri1 .remove_from_mesh ( msh );
			tri2 .remove_from_mesh ( msh );
			assert ( tri5.exists() );
			assert ( tri6.exists() );
			AC .reverse() .cut_from_bdry_of ( tri5, tag::do_not_bother );
			DA .reverse() .cut_from_bdry_of ( tri6, tag::do_not_bother );
			CB .glue_on_bdry_of ( tri5, tag::do_not_bother );
			BD .glue_on_bdry_of ( tri6, tag::do_not_bother );

			// replace A by B in segments neighbour to A
			for ( std::forward_list < Cell > ::iterator it_list = list_of_segs .begin();
			      it_list != list_of_segs .end(); it_list ++                            )
			{	Cell segm = *it_list;
				assert ( segm.tip() == A );
				A .cut_from_bdry_of ( segm, tag::do_not_bother );
				B .glue_on_bdry_of ( segm, tag::do_not_bother );  }

			tri5 .boundary() .closed_loop ( B );
			tri6 .boundary() .closed_loop ( B );                                          }

		else assert ( false );    // at least one wing should be free
		
	}  // end of  for  over segments of msh
	
}  // end of  remove_short_segments
	

void baricenters ( Mesh & msh )

// segments on the boundary are not moved -- we use method  is_inner_to

{	CellIterator it = msh.iterator ( tag::over_vertices );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell P = *it;
		if ( P.is_inner_to ( msh ) ) msh.baricenter ( P, tag::spin );  }  }
	
