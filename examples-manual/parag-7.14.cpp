
// example presented in paragraph 7.14 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// folds a square around a cylinder

#include "maniFEM.h"
using namespace maniFEM;


void remove_short_segments ( Mesh & msh, double threshold );

void baricenters ( Mesh & msh );


int main ( )

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0], y = xy[1];

	size_t n = 20;
	double d = 2.6 / double(n);

	Cell A ( tag::vertex );  x(A) = -1.3;  y(A) = -1.3;
	Cell B ( tag::vertex );  x(B) =  1.3;  y(B) = -1.3;
	Cell C ( tag::vertex );  x(C) =  1.3;  y(C) =  1.3;
	Cell D ( tag::vertex );  x(D) = -1.3;  y(D) =  1.3;

	Mesh AB ( tag::segment, A.reverse(), B, tag::divided_in, n );
	Mesh BC ( tag::segment, B.reverse(), C, tag::divided_in, n );
	Mesh CD ( tag::segment, C.reverse(), D, tag::divided_in, n );
	Mesh DA ( tag::segment, D.reverse(), A, tag::divided_in, n );

	Manifold circle = RR2.implicit ( x*x + y*y == 1. );
	Mesh inner ( tag::progressive, tag::entire_manifold, circle, tag::desired_length, d );

	Mesh bdry ( tag::join, AB, BC, CD, DA, inner.reverse() );

	RR2.set_as_working_manifold();
	Mesh square ( tag::progressive, tag::boundary, bdry, tag::desired_length, d );

	Mesh cyl = square.fold
		( tag::identify, BC, tag::with, DA.reverse(), tag::use_existing_vertices );

	cyl.draw_ps ( "cylinder.eps", tag::unfold,
               tag::over_region, -2.1 < x < 4.3, -3.6 < y < 2.1 );

	remove_short_segments ( cyl, 0.5*d );
	// remove segments shorter than the given threshold
	
	baricenters ( cyl );
	// re-define the position of each vertex as baricenter of its neighbours
	
	cyl.draw_ps ( "cylinder-0.5.eps", tag::unfold,
               tag::over_region, -2.1 < x < 4.3, -3.6 < y < 2.1 );

	remove_short_segments ( cyl, 0.6*d );
	baricenters ( cyl );
	
	cyl.draw_ps ( "cylinder-0.6.eps", tag::unfold,
               tag::over_region, -2.1 < x < 4.3, -3.6 < y < 2.1 );

	remove_short_segments ( cyl, 0.65*d );
	baricenters ( cyl );
	
	cyl.draw_ps ( "cylinder-0.65.eps", tag::unfold,
               tag::over_region, -2.1 < x < 4.3, -3.6 < y < 2.1 );

	std::cout << "produced files cylinder-*.eps - please edit before viewing" << std::endl;	

} // end of main

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
	
