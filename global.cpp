
// global.cpp 2021.09.19

//   This file is part of maniFEM, a C++ library for meshes and finite elements on manifolds.

//   Copyright 2019, 2020, 2021 Cristian Barbarosie cristian.barbarosie@gmail.com
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

#include <fstream>
#include "maniFEM.h"

using namespace maniFEM;


void Mesh::build ( const tag::Segment &, const Cell & A, const Cell & B,
                   const tag::DividedIn &, size_t n                      )

// see paragraph 12.2 in the manual
	
{	Manifold space = Manifold::working;
	assert ( space.exists() );  // we use the current manifold

	assert ( not A.is_positive() );
	Cell posA = A.reverse();
	assert ( posA.is_positive() );
	assert ( B.is_positive() );
	assert ( A.dim() == 0 );
	assert ( B.dim() == 0 );

	Cell prev_point = A;
	for ( size_t i=1; i < n; ++i )
	{	Cell P ( tag::vertex );
		double frac = double(i)/double(n);
		space.interpolate ( P, 1.-frac, posA, frac, B );
		Cell seg ( tag::segment, prev_point, P );
		assert ( seg.exists() );
		seg.add_to_mesh ( *this, tag::do_not_bother );
		assert ( P.exists() );
		prev_point = P.reverse();                         }
	Cell seg ( tag::segment, prev_point, B );
	seg.add_to_mesh ( *this, tag::do_not_bother );

}  // end of Mesh::build with tag::segment

//----------------------------------------------------------------------------------//


void Mesh::build ( const tag::Segment &, const Cell & A, const Cell & B,
                   const tag::DividedIn &, size_t n,
                   const tag::Spin &, const tag::Util::CompositionOfActions & s )

// see paragraph 12.2 in the manual
// in this version, the segment may be a loop around a cylinder or torus
// beware, A.reverse() and B may be one and the same vertex !
	
{	Manifold space = Manifold::working;
	assert ( space.exists() );  // we use the current (quotient) manifold
	Manifold::Quotient * mani_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space.core );
	Function coords_q = space.coordinates();
	Manifold mani_Eu = mani_q->base_space;  // underlying Euclidian manifold
	Function coords_Eu = mani_Eu.coordinates();

	assert ( not A.is_positive() );
	Cell posA = A.reverse();
	assert ( posA.is_positive() );
	assert ( B.is_positive() );
	assert ( A.dim() == 0 );
	assert ( B.dim() == 0 );

	Cell shadow_of_B ( tag::vertex );
	if ( coords_Eu.nb_of_components() == 1 )
	{	double new_co = coords_q ( B, tag::spin, s );
		coords_Eu ( shadow_of_B ) = new_co;            }
	else
	{	assert ( coords_Eu.nb_of_components() > 1 );
		std::vector < double > new_co = coords_q ( B, tag::spin, s );
		coords_Eu ( shadow_of_B ) = new_co;                           }

	Cell prev_point = A;
	for ( size_t i=1; i < n; ++i )
	{	Cell P ( tag::vertex );
		double frac = double(i)/double(n);
		mani_Eu.interpolate ( P, 1.-frac, posA, frac, shadow_of_B );
		Cell seg ( tag::segment, prev_point, P );
		assert ( seg.exists() );
		seg.add_to_mesh ( *this, tag::do_not_bother );
		assert ( P.exists() );
		prev_point = P.reverse();                                     }

	Cell seg ( tag::segment, prev_point, B );
	// size_t nb_spins = mani_q->spins.size();
	// for ( size_t i = 0; i < nb_spins; i++ )
	// {	Field::ShortInt & sp = mani_q->spins[i];
	// 	std::map<Function::Action,short int>::const_iterator it =
	// 		s.index_map.lower_bound(mani_q->actions[i]);
	// 	sp.on_cell ( seg.core ) = it->second;                     }
	seg.spin() = s;
	seg.add_to_mesh ( *this, tag::do_not_bother );

	space.set_as_working_manifold();

}  // end of Mesh::build with tag::segment and tag::spin

//----------------------------------------------------------------------------------//


void Mesh::build ( const tag::Triangle &, const Mesh & AB, const Mesh & BC, const Mesh & CA )

// see paragraph 12.4 in the manual
	
{	Manifold space = Manifold::working;
	assert ( space.exists() );  // we use the current manifold

	// sides must be split in the same number of segments :
	size_t N = AB.number_of ( tag::segments );
	assert ( N == BC.number_of ( tag::segments ) );
	assert ( N == CA.number_of ( tag::segments ) );

	Cell A = AB.first_vertex().reverse();
	Cell B = BC.first_vertex().reverse();
	Cell C = CA.first_vertex().reverse();
	assert ( A == CA.last_vertex() );
	assert ( B == AB.last_vertex() );
	assert ( C == BC.last_vertex() );

	// we keep a list of horizontal segments (parallel to AB)
	// useful for the next layer of triangles
	std::list < Cell > ground, ceiling;
	{ // just a block of code for hiding 'it'
	CellIterator it = AB.iterator ( tag::over_segments, tag::require_order );
	for ( it.reset(); it.in_range(); it++ ) ground.push_back ( *it );
	} // just a block of code for hiding 'it'

	// we shall use six vertices, two on AB, two on BC, two on CA
	// like shadows of the point currently buing built
	Cell Q_AB_ini = A, P_BC = B, Q_CA = A;
	Cell seg_on_BC = BC.cell_in_front_of ( B );

	for ( size_t i = 1; i < N; i++ ) // "vertical" movement
	{	// advance one level upwards and slightly right (parallel to CA)
		std::list<Cell>::iterator it_ground = ground.begin();
		Cell ground_seg = *it_ground;
		assert ( ground_seg.base().reverse() == Q_CA );
		Cell ground_ver = ground_seg.tip();
		Cell seg_on_AB = AB.cell_in_front_of ( Q_AB_ini );
		Q_AB_ini = seg_on_AB.tip();
		Cell Q_AB = Q_AB_ini;
		assert ( seg_on_BC == BC.cell_in_front_of ( P_BC ) );
		P_BC = seg_on_BC.tip();
		Cell seg_on_CA = CA.cell_behind(Q_CA);
		Cell P_CA = Q_CA = seg_on_CA.base().reverse();
		Cell P_AB = A, Q_BC = C;
		// build the first triangle on this layer
		Cell previous_seg ( tag::segment, ground_ver.reverse(), P_CA );
		Cell tri ( tag::triangle, ground_seg, previous_seg, seg_on_CA );
		tri.add_to_mesh ( *this );  // 'this' is the mesh we are building
		Cell previous_ver = Q_CA;
		ceiling.clear();
		for ( size_t j = i+1; j <= N; j++ ) // "horizontal" movement
		{	// advance one step horizontally (parallel to AB)
			P_AB = AB.cell_in_front_of(P_AB).tip();
			Q_AB = AB.cell_in_front_of(Q_AB).tip();
			Q_BC = BC.cell_behind(Q_BC).base().reverse();
			P_CA = CA.cell_behind(P_CA).base().reverse();
			Cell S ( tag::non_existent );  // temporary non-existent cell
			if ( j == N ) S = P_BC;
			else
			{	// we prepare for building a new point S and we need fractions
				// distance to AB : i
				// distance to BC : N-j
				// distance to CA : j-i
				double frac_AB = 1. / double(i),
				       frac_BC = 1. / double(N-j),
				       frac_CA = 1. / double(j-i);
				double s = 2.* ( frac_AB + frac_BC + frac_CA );
				frac_AB /= s;  frac_BC /= s;  frac_CA /= s;
				S = Cell ( tag::vertex );
				space.interpolate ( S, frac_AB, P_AB,  frac_AB, Q_AB,
														   frac_BC, P_BC,  frac_BC, Q_BC,
														   frac_CA, P_CA,  frac_CA, Q_CA );  }
			Cell new_seg ( tag::segment, ground_ver.reverse(), S );
			Cell horizontal_seg ( tag::segment, S.reverse(), previous_ver );
		  Cell tri_1 ( tag::triangle, previous_seg.reverse(), new_seg, horizontal_seg );
		  tri_1.add_to_mesh ( *this );  // 'this' is the mesh we are building
			it_ground++;  assert ( it_ground != ground.end() );
			ground_seg = *it_ground;
			if ( j == N ) previous_seg = seg_on_BC;
			else
			{	ground_ver = ground_seg.tip();
				previous_seg = Cell ( tag::segment, ground_ver.reverse(), S );  }
		  Cell tri_2 ( tag::triangle, ground_seg, previous_seg, new_seg.reverse() );
			tri_2.add_to_mesh ( *this );  // 'this' is the mesh we are building
			previous_ver = S;
			// add horizontal_seg.reverse() to future ground
			ceiling.push_back ( horizontal_seg.reverse() );                              	  }
		assert ( seg_on_BC.tip() == P_BC );
		seg_on_BC = BC.cell_in_front_of ( P_BC );
		ground = ceiling;                                                                    }
		// improve by moving ceiling to ground, leaving ceiling empty !

	// last triangle
	Cell seg_on_CA = CA.cell_behind ( Q_CA );
	assert ( not ground.empty() );
	Cell ground_seg = *(ground.begin());
	Cell tri ( tag::triangle, seg_on_BC, seg_on_CA, ground_seg );
	tri.add_to_mesh ( *this );  // 'this' is the mesh we are building
		
} // end of Mesh::build with tag::triangle

//----------------------------------------------------------------------------------//


namespace {  // anonymous namespace, mimics static linkage

Cell find_common_vertex ( const Mesh & seg1, const Mesh & seg2 )

// we look for a common vertex, where seg1 ends and seg2 begins (in this order)
// this does not apply to closed loops of course, so we give them a different treatment

{	std::vector < Cell > vec;
	CellIterator it = seg1.iterator ( tag::over_vertices );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell V = *it;
		if ( V.belongs_to ( seg2 ) ) vec.push_back ( V );  }
	assert ( vec.size() > 0 );
	if ( vec.size() == 1 )  // one common vertex only, so we have no doubt
		return vec[0];
	assert ( vec.size() == 2 );  // it makes no sense to have more than two
	// we are looking for the one where seg1 ends and seg2 begins
	Cell V = vec[0];
	if ( seg1.cell_in_front_of ( V, tag::may_not_exist ) .exists() )
	{	// seg1 does not end in V, so this is not the vertex we are looking for
		// perhaps the other one ?
		V = vec[1];
		assert ( not seg1.cell_in_front_of ( V, tag::may_not_exist ) .exists() );
		// seg1 ends in V
		assert ( not seg2.cell_behind ( V, tag::may_not_exist ) .exists() );
		// seg2 begins in V
		return V;                                                                 }
	// else : seg1 ends in V
	assert ( not seg2.cell_behind ( V, tag::may_not_exist ) .exists() );
	// seg2 begins in V
	return V;                                                                      }

} // end of anonymous namespace


void Mesh::build ( const tag::Triangle &,
                   const Mesh & AB, const Mesh & BC, const Mesh & CA,
                   const tag::Spin &                                 )

// see paragraph 12.4 in the manual
// the tag:::spin tells maniFEM that we are on a quotient manifold
// and that the segments provided (AB, BC, CA) may have spin
	
// beware, sides may be closed loops
	
{	Manifold space = Manifold::working;
	assert ( space.exists() );  // we use the current manifold
	Manifold::Quotient * mani_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space.core );
	Function coords_q = space.coordinates();
	Manifold mani_Eu = mani_q->base_space;  // underlying Euclidian manifold
	Function coords_Eu = mani_Eu.coordinates();

	// sides must be split in the same number of segments :
	size_t N = AB.number_of ( tag::segments );
	assert ( N == BC.number_of ( tag::segments ) );
	assert ( N == CA.number_of ( tag::segments ) );

	// recover corners from the sides
	// the process is different from the one in 'build' without spin
	// here, sides may be closed loops and then methods 'first_vertex' and 'last_vertex'
	// become meaningless
	Cell A = find_common_vertex ( CA, AB );
	Cell B = find_common_vertex ( AB, BC );
	Cell C = find_common_vertex ( BC, CA );

	// we keep a list of horizontal segments (parallel to AB)
	// useful for the next layer of triangles
	std::list < Cell > ground, ceiling;
	{ // just a block of code for hiding 'it'
		CellIterator it = AB.iterator ( tag::over_segments, tag::require_order );
	for ( it.reset(); it.in_range(); it++ ) ground.push_back ( *it );
	} // just a block of code for hiding 'it'

	// we shall use six vertices, two on AB, two on BC, two on CA
	// like shadows of the point currently buing built
	Cell Q_AB_ini = A, P_BC = B, Q_CA = A;
	Cell seg_P_BC = BC.cell_in_front_of ( B );
	// we use six shadow vertices for interpolation
	Cell shadow_P_AB ( tag::vertex ), shadow_Q_AB ( tag::vertex ),
	     shadow_P_BC ( tag::vertex ), shadow_Q_BC ( tag::vertex ),
	     shadow_P_CA ( tag::vertex ), shadow_Q_CA ( tag::vertex );

  // we keep spins of vertices relative to A
	Function::CompositionOfActions spin_B = 0, spin_C = 0;
	{ // just a block of code for hiding 'it'
	CellIterator it = AB.iterator ( tag::over_segments );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell seg = *it;  spin_B += seg.spin();  }
	} { // just a block of code for hiding 'it'
	CellIterator it = CA.iterator ( tag::over_segments );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell seg = *it;  spin_C -= seg.spin();  }
	} { // just a block of code for hiding 'it'
	Function::CompositionOfActions spin_test = spin_B;
	CellIterator it = BC.iterator ( tag::over_segments );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell seg = *it;  spin_test += seg.spin();  }
	assert ( spin_test == spin_C );
	} // just a block of code for hiding 'it'
	
	Function::CompositionOfActions spin_Q_AB_ini = 0, spin_P_BC = spin_B, spin_Q_CA = 0;

	for ( size_t i = 1; i < N; i++ ) // "vertical" movement
	{	// advance one level upwards and slightly right (parallel to CA)
		std::list<Cell>::iterator it_ground = ground.begin();
		Cell ground_seg = *it_ground;
		assert ( ground_seg.base().reverse() == Q_CA );
		Cell ground_ver = ground_seg.tip();
		Function::CompositionOfActions spin_ground_ver = spin_Q_CA + ground_seg.spin();
		Cell seg_Q_AB = AB.cell_in_front_of ( Q_AB_ini );
		Q_AB_ini = seg_Q_AB.tip();
		spin_Q_AB_ini += seg_Q_AB.spin();
		Cell Q_AB = Q_AB_ini;
		Function::CompositionOfActions spin_Q_AB = spin_Q_AB_ini;
		assert ( seg_P_BC == BC.cell_in_front_of ( P_BC ) );
		P_BC = seg_P_BC.tip();
		spin_P_BC += seg_P_BC.spin();
		Cell seg_Q_CA = CA.cell_behind(Q_CA);
		Cell P_CA = Q_CA = seg_Q_CA.base().reverse();
		spin_Q_CA -= seg_Q_CA.spin();
		Function::CompositionOfActions spin_P_CA = spin_Q_CA;
		std::vector < double > v = coords_q ( Q_CA, tag::spin, spin_Q_CA );
		coords_Eu ( shadow_Q_CA ) = v;
		v = coords_q ( P_BC, tag::spin, spin_P_BC );
		coords_Eu ( shadow_P_BC ) = v;
		Cell P_AB = A, Q_BC = C;
		Function::CompositionOfActions spin_P_AB = 0, spin_Q_BC = spin_C;
		// build the first triangle on this layer
		Cell previous_seg ( tag::segment, ground_ver.reverse(), P_CA );
		previous_seg.spin() = - ground_seg.spin() - seg_Q_CA.spin();
		Cell tri ( tag::triangle, ground_seg, previous_seg, seg_Q_CA );
		tri.add_to_mesh ( *this );  // 'this' is the mesh we are building
		Cell previous_ver = Q_CA;
		Function::CompositionOfActions spin_prev_ver = spin_Q_CA;
		ceiling.clear();
		for ( size_t j = i+1; j <= N; j++ ) // "horizontal" movement
		{	// advance one step horizontally (parallel to AB)
			Cell seg_P_AB = AB.cell_in_front_of ( P_AB );
			P_AB = seg_P_AB.tip();
			spin_P_AB += seg_P_AB.spin();
			seg_Q_AB = AB.cell_in_front_of ( Q_AB );
			Q_AB = seg_Q_AB.tip();
			spin_Q_AB += seg_Q_AB.spin();
			Cell seg_Q_BC = BC.cell_behind ( Q_BC );
			Q_BC = seg_Q_BC.base().reverse();
			spin_Q_BC -= seg_Q_BC.spin();
			Cell seg_P_CA = CA.cell_behind ( P_CA );
			P_CA = seg_P_CA.base().reverse();
			spin_P_CA -= seg_P_CA.spin();
			Cell S ( tag::non_existent );  // temporary non-existent cell
			Function::CompositionOfActions spin_S = 0;
			if ( j == N )  { S = P_BC; spin_S = spin_P_BC;  }
			else
			{	// we prepare for building a new point S and we need fractions
				// distance to AB : i
				// distance to BC : N-j
				// distance to CA : j-i
				double frac_AB = 1. / double(i),
				       frac_BC = 1. / double(N-j),
				       frac_CA = 1. / double(j-i);
				double s = 2.* ( frac_AB + frac_BC + frac_CA );
				frac_AB /= s;  frac_BC /= s;  frac_CA /= s;
				S = Cell ( tag::vertex );
				v = coords_q ( P_AB, tag::spin, spin_P_AB );
				coords_Eu ( shadow_P_AB ) = v;
				v = coords_q ( Q_AB, tag::spin, spin_Q_AB );
				coords_Eu ( shadow_Q_AB ) = v;
				v = coords_q ( Q_BC, tag::spin, spin_Q_BC );
				coords_Eu ( shadow_Q_BC ) = v;
				v = coords_q ( P_CA, tag::spin, spin_P_CA );
				coords_Eu ( shadow_P_CA ) = v;
				mani_Eu.interpolate ( S, frac_AB, shadow_P_AB,  frac_AB, shadow_Q_AB,
													       frac_BC, shadow_P_BC,  frac_BC, shadow_Q_BC,
													       frac_CA, shadow_P_CA,  frac_CA, shadow_Q_CA );  }
			Cell new_seg ( tag::segment, ground_ver.reverse(), S );
			new_seg.spin() = spin_S - spin_ground_ver;
			Cell horizontal_seg ( tag::segment, S.reverse(), previous_ver );
			horizontal_seg.spin() = spin_prev_ver - spin_S;
			assert ( horizontal_seg.spin() + new_seg.spin() - previous_seg.spin() == 0 );
		  Cell tri_1 ( tag::triangle, previous_seg.reverse(), new_seg, horizontal_seg );
		  tri_1.add_to_mesh ( *this );  // 'this' is the mesh we are building
			it_ground++;  assert ( it_ground != ground.end() );
			ground_seg = *it_ground;
			if ( j == N ) previous_seg = seg_P_BC;
			else
			{	ground_ver = ground_seg.tip();
				previous_seg = Cell ( tag::segment, ground_ver.reverse(), S );
				previous_seg.spin() = spin_S - spin_ground_ver;                }
			assert ( ground_seg.spin() + previous_seg.spin() - new_seg.spin() == 0 );
		  Cell tri_2 ( tag::triangle, ground_seg, previous_seg, new_seg.reverse() );
			tri_2.add_to_mesh ( *this );  // 'this' is the mesh we are building
			previous_ver = S;  spin_prev_ver = spin_S;
			// add horizontal_seg.reverse() to future ground
			ceiling.push_back ( horizontal_seg.reverse() );                              	  }
		assert ( seg_P_BC.tip() == P_BC );
		seg_P_BC = BC.cell_in_front_of ( P_BC );
		ground = ceiling;                                                                    }
		// improve by moving ceiling to ground, leaving ceiling empty !

	// last triangle
	Cell seg_on_CA = CA.cell_behind ( Q_CA );
	assert ( not ground.empty() );
	Cell ground_seg = *(ground.begin());
	assert ( seg_P_BC.spin() + seg_on_CA.spin() + ground_seg.spin() == 0 );
	Cell tri ( tag::triangle, seg_P_BC, seg_on_CA, ground_seg );
	tri.add_to_mesh ( *this );  // 'this' is the mesh we are building
		
} // end of Mesh::build with tag::triangle and tag::spin

//----------------------------------------------------------------------------------//


void Mesh::build ( const tag::Quadrangle &, const Mesh & south, const Mesh & east,
                   const Mesh & north, const Mesh & west, bool cut_rectangles_in_half )

// see paragraph 12.3 in the manual
	
{	Manifold space = Manifold::working;
	assert ( space.exists() );  // we use the current manifold

	// recover corners from the sides
	Cell SW = south.first_vertex().reverse();
	assert ( SW == west.last_vertex() );
	Cell SE = east.first_vertex().reverse();
	assert ( SE == south.last_vertex() );
	Cell NE = north.first_vertex().reverse();
	assert ( NE == east.last_vertex() );
	Cell NW = west.first_vertex().reverse();
	assert ( NW == north.last_vertex() );
	size_t N_horiz = south.number_of ( tag::segments );
	assert ( N_horiz == north.number_of ( tag::segments ) );
	size_t N_vert = east.number_of ( tag::segments );
	assert ( N_vert == west.number_of ( tag::segments ) );

	// prepare horizon
	std::list <Cell> horizon;
	{ // just a block of code for hiding 'it'
	CellIterator it = south.iterator ( tag::over_segments, tag::require_order );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell seg = *it;  horizon.push_back ( seg );  }
	} // just a block of code for hiding 'it'

	// start mesh generation
	CellIterator it_east = east.iterator ( tag::over_vertices, tag::require_order );
	CellIterator it_west = west.iterator ( tag::over_vertices, tag::backwards );
	CellIterator it_south = south.iterator ( tag::over_vertices, tag::require_order );
	CellIterator it_north = north.iterator ( tag::over_vertices, tag::backwards );
	it_east.reset();  it_east++;
	it_west.reset();  it_west++;
	for ( size_t i = 1; i < N_vert; ++i )
	{	std::list<Cell>::iterator it = horizon.begin();
		Cell seg = *it;
		Cell A = seg.base().reverse();
		Cell DA = west.cell_behind ( A, tag::surely_exists );
		Cell D = DA.base().reverse();
		Cell ver_east = *it_east;
		Cell ver_west = *it_west;
		assert ( ver_west == D );
		double frac_N = double(i) / double(N_vert),  alpha = frac_N * (1-frac_N);
		alpha = alpha*alpha*alpha;
		it_south.reset();  it_south++;
		it_north.reset();  it_north++;
		for ( size_t j = 1; j < N_horiz; j++ )
		{	Cell AB = *it;  // 'it' points into the 'horizon' list of segments
			Cell B = AB.tip();
			Cell C ( tag::vertex );  // create a new vertex
			Cell ver_south = *it_south;
			Cell ver_north = *it_north;
			double frac_E = double(j) / double(N_horiz),  beta = frac_E * (1-frac_E);
			beta = beta*beta*beta;
			double sum = alpha + beta, aa = alpha/sum, bb = beta/sum;
			space.interpolate ( C, bb*(1-frac_N), ver_south, aa*frac_E,     ver_east,     
		                         bb*frac_N,     ver_north, aa*(1-frac_E), ver_west );
			Cell BC ( tag::segment, B.reverse(), C );  // create a new segment
			Cell CD ( tag::segment, C.reverse(), D );  // create a new segment
			if ( cut_rectangles_in_half )
			{	Cell BD ( tag::segment, B.reverse(), D );  // create a new segment
				Cell BCD ( tag::triangle, BD.reverse(), BC, CD );  // create a new triangle
				Cell ABD ( tag::triangle, BD, DA, AB );  // create a new triangle
				BCD.add_to_mesh (*this);  // 'this' is the mesh we are building
				ABD.add_to_mesh (*this);                          }
			else // with quadrilaterals
			{	Cell Q ( tag::rectangle, AB, BC, CD, DA );  // create a new rectangle
				Q.add_to_mesh (*this);                      }
			// CD is on the ceiling, we keep it in the 'horizon' list
			// it will be on the ground when we build the next layer of cells
			*it = CD.reverse(); // 'it' points into the 'horizon' list
			it++;
			D = C;
			DA = BC.reverse();
			it_south++;  assert ( it_south.in_range() );
			it_north++;  assert ( it_north.in_range() );
		} // end of for j
		it_south++;  assert ( not it_south.in_range() );
		it_north++;  assert ( not it_north.in_range() );
		// last rectangle of this row, east side already exists
		Cell AB = *it;
		Cell B = AB.tip();
		Cell BC = east.cell_in_front_of ( B, tag::surely_exists );
		Cell C = BC.tip();
		Cell CD ( tag::segment, C.reverse(), D );  // create a new segment
		if ( cut_rectangles_in_half )
		{	Cell BD ( tag::segment, B.reverse(), D );  // create a new segment
			Cell BCD ( tag::triangle, BD.reverse(), BC, CD );  // create a new triangle
			Cell ABD ( tag::triangle, BD, DA, AB );  // create a new triangle
			BCD.add_to_mesh (*this);
			ABD.add_to_mesh (*this);                          }
		else // with quadrilaterals
		{	Cell Q ( tag::rectangle, AB, BC, CD, DA );  // create a new rectangle
			Q.add_to_mesh (*this);                     }
		*it = CD.reverse();
		it_east++;  assert ( it_east.in_range() );
		it_west++;  assert ( it_west.in_range() );
		it++;  assert ( it == horizon.end() );
	} // end of for i
	it_east++;  assert ( not it_east.in_range() );
	it_west++;  assert ( not it_west.in_range() );
	// last row of rectangles is different, north sides already exist
	std::list<Cell>::iterator it = horizon.begin();
	Cell DA = west.cell_in_front_of ( NW, tag::surely_exists );
	Cell D = NW;
	for (size_t j=1; j < N_horiz; j++)
	{	Cell AB = *it;
		Cell B = AB.tip();
		Cell CD = north.cell_behind ( D );
		Cell C = CD.base().reverse();
		Cell BC ( tag::segment, B.reverse(), C );  // create a new segment
		if ( cut_rectangles_in_half )
		{	Cell BD ( tag::segment, B.reverse(), D );  // create a new segment
			Cell BCD ( tag::triangle, BD.reverse(), BC, CD );  // create a new triangle
			Cell ABD ( tag::triangle, BD, DA, AB );  // create a new triangle
			BCD.add_to_mesh (*this);  // 'this' is the mesh we are building
			ABD.add_to_mesh (*this);                           }
		else // with quadrilaterals
		{	Cell Q ( tag::rectangle, AB, BC, CD, DA );  // create a new rectangle
			Q.add_to_mesh (*this);                     }
		it++;
		D = C;
		DA = BC.reverse();                                          }
	// and the last rectangle of the last row
	Cell AB = *it;
	Cell B = AB.tip();
	Cell BC = east.cell_in_front_of (B);
	Cell C = BC.tip();
	assert ( C == NE );
	Cell CD = north.cell_behind ( D );
	if ( cut_rectangles_in_half )
	{	Cell BD ( tag::segment, B.reverse(), D );
		Cell BCD ( tag::triangle, BD.reverse(), BC, CD );
		Cell ABD ( tag::triangle, BD, DA, AB );
		BCD.add_to_mesh (*this);  // 'this' is the mesh we are building
		ABD.add_to_mesh (*this);                          }
	else // with quadrilaterals
	{	Cell Q ( tag::rectangle, AB, BC, CD, DA );
		Q.add_to_mesh (*this);                     }
	it++;  assert ( it == horizon.end() );

} // end of Mesh::build with tag::quadrangle

//----------------------------------------------------------------------------------//


void Mesh::build ( const tag::Quadrangle &, const Mesh & south, const Mesh & east,
                   const Mesh & north, const Mesh & west, bool cut_rectangles_in_half,
                   const tag::Spin &                                                   )

// see paragraph 12.3 in the manual
// the tag:::spin tells maniFEM that we are on a quotient manifold
// and that the segments provided (south, east, north, west) may have spin

// beware, south may be equal to north.reverse, east may be equal to west.reverse
// or they may be not equal but share the same vertices (and segments, reversed)
// beware, the correspondence may be not face-to-face,
// e.g. the first vertex of south may show up somewhere in the middle of north

// beware, sides may be closed loops
	
{	Manifold space = Manifold::working;
	assert ( space.exists() );  // we use the current manifold
	Manifold::Quotient * mani_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space.core );
	Function coords_q = space.coordinates();
	Manifold mani_Eu = mani_q->base_space;  // underlying Euclidian manifold
	Function coords_Eu = mani_Eu.coordinates();

	// recover corners from the sides
	// the process is different from the one in 'build' without spin
	// here, sides may be closed loops and then methods 'first_vertex' and 'last_vertex'
	// become meaningless
	Cell SW = find_common_vertex ( west, south );
	Cell SE = find_common_vertex ( south, east );
	Cell NE = find_common_vertex ( east, north );
	Cell NW = find_common_vertex ( north, west );
	
	size_t N_horiz = south.number_of ( tag::segments );
	assert ( N_horiz == north.number_of ( tag::segments ) );
	size_t N_vert = east.number_of ( tag::segments );
	assert ( N_vert == west.number_of ( tag::segments ) );

	// prepare horizon
	std::list <Cell> horizon;
	{ // just a block of code for hiding 'it'
	CellIterator it = south.iterator ( tag::over_segments, tag::require_order );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell seg = *it;  horizon.push_back ( seg );  }
	} // just a block of code for hiding 'it'

	// we have to deal with possible spins
	// we choose that, at each interpolation operation, i.e., for each new vertex,
	// the new vertex will have spin zero relatively to SW
	// we must keep track of the spins of ver_south, ver_east, ver_north, ver_west, B, D
	// (relatively to SW)
	// we use four shadow vertices for interpolation
	Cell shadow_south ( tag::vertex ), shadow_east ( tag::vertex );
	Cell shadow_north ( tag::vertex ), shadow_west ( tag::vertex );
	
	Function::CompositionOfActions spin_NW = 0, spin_SE = 0;
	// spin_SW is zero by our choice, spin_NE is not needed
	{ // just a block of code for hiding 'it'
	CellIterator it = south.iterator ( tag::over_segments );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell seg = *it;  spin_SE += seg.spin();  }
	} { // just a block of code for hiding 'it'
	CellIterator it = west.iterator ( tag::over_segments );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell seg = *it;  spin_NW -= seg.spin();  }
	} // just a block of code for hiding 'it'

	// start mesh generation
	// sides may be closed loops
	// have SW, SE, NW and NE been correctly defined ?
	Cell seg_east = east.cell_in_front_of ( SE, tag::surely_exists );
	Cell seg_west = west.cell_behind ( SW, tag::surely_exists );
	Function::CompositionOfActions spin_ver_west = 0;  // spin_SW is zero by our choice
	Function::CompositionOfActions spin_ver_east = spin_SE;
	Function::CompositionOfActions spin_B = 0;  // spin_SW is zero by our choice
	for ( size_t i = 1; i < N_vert; ++i )
	{	std::list<Cell>::iterator it = horizon.begin();
		Cell AB = *it;
		Cell A = AB.base().reverse();
	  Cell DA = seg_west;
		Cell D = DA.base().reverse();
		spin_ver_east += seg_east.spin();
		Cell ver_east = seg_east.tip();
		seg_east = east.cell_in_front_of ( ver_east, tag::surely_exists );
		spin_ver_west -= seg_west.spin();
		Cell ver_west = seg_west.base().reverse();
		seg_west = west.cell_behind ( ver_west, tag::surely_exists );
		assert ( ver_west == D );
		double frac_N = double(i) / double(N_vert),  alpha = frac_N * (1-frac_N);
		alpha = alpha*alpha*alpha;
		std::vector < double > v = coords_q ( ver_east, tag::spin, spin_ver_east );
		coords_Eu ( shadow_east ) = v;
		v = coords_q ( ver_west, tag::spin, spin_ver_west );
		coords_Eu ( shadow_west ) = v;
		Cell seg_south = south.cell_in_front_of ( SW, tag::surely_exists );
		Cell seg_north = north.cell_behind ( NW, tag::surely_exists );
		Function::CompositionOfActions spin_ver_south = 0;  // spin_SW is zero by our choice
		Function::CompositionOfActions spin_ver_north = spin_NW;
		Function::CompositionOfActions spin_D = spin_ver_west;
		for ( size_t j = 1; j < N_horiz; j++ )
		{	AB = *it;  // 'it' points into the 'horizon' list of segments
			Cell B = AB.tip();
			spin_B += AB.spin();
			spin_ver_south += seg_south.spin();
			Cell ver_south = seg_south.tip();
			seg_south = south.cell_in_front_of ( ver_south );
			spin_ver_north -= seg_north.spin();
			Cell ver_north = seg_north.base().reverse();
			seg_north = north.cell_behind ( ver_north );
			Cell C ( tag::vertex );  // create a new vertex
			double frac_E = double(j) / double(N_horiz),  beta = frac_E * (1-frac_E);
			beta = beta*beta*beta;
			double sum = alpha + beta, aa = alpha/sum, bb = beta/sum;
			v = coords_q ( ver_south, tag::spin, spin_ver_south );
			coords_Eu ( shadow_south ) = v;
			v = coords_q ( ver_north, tag::spin, spin_ver_north );
			coords_Eu ( shadow_north ) = v;
			mani_Eu.interpolate ( C, bb*(1-frac_N), shadow_south, aa*frac_E,     shadow_east,     
		                           bb*frac_N,     shadow_north, aa*(1-frac_E), shadow_west );
			Cell BC ( tag::segment, B.reverse(), C );  // create a new segment
			Cell CD ( tag::segment, C.reverse(), D );  // create a new segment
			BC.spin() = -spin_B;
			CD.spin() =  spin_D;
		  assert ( AB.spin() + BC.spin() + CD.spin() + DA.spin() == 0 );
		  spin_D = 0;
			if ( cut_rectangles_in_half )
			{	Cell BD ( tag::segment, B.reverse(), D );  // create a new segment
				BD.spin() = spin_D - spin_B;
				Cell BCD ( tag::triangle, BD.reverse(), BC, CD );  // create a new triangle
				Cell ABD ( tag::triangle, BD, DA, AB );  // create a new triangle
				BCD.add_to_mesh (*this);  // 'this' is the mesh we are building
				ABD.add_to_mesh (*this);                           }
			else // with quadrilaterals
			{	Cell Q ( tag::rectangle, AB, BC, CD, DA );  // create a new rectangle
				Q.add_to_mesh (*this);                     }
			// CD is on the ceiling, we keep it in the 'horizon' list
			// it will be on the ground when we build the next layer of cells
			*it = CD.reverse(); // 'it' points into the 'horizon' list
			it++;
			D = C;
			DA = BC.reverse();
		} // end of for j
		// last rectangle of this row, east side already exists
		AB = *it;
		Cell B = AB.tip();
		Cell BC = east.cell_in_front_of ( B, tag::surely_exists );
		Cell C = BC.tip();
		Cell CD ( tag::segment, C.reverse(), D );  // create a new segment
		CD.spin() = - DA.spin() - AB.spin() - BC.spin();
		assert ( spin_D == 0 );
		spin_B = spin_ver_west;
		if ( cut_rectangles_in_half )
		{	Cell BD ( tag::segment, B.reverse(), D );  // create a new segment
			BD.spin() = - AB.spin() - DA.spin();
			Cell BCD ( tag::triangle, BD.reverse(), BC, CD );  // create a new triangle
			Cell ABD ( tag::triangle, BD, DA, AB );  // create a new triangle
			BCD.add_to_mesh (*this);  // 'this' is the mesh we are building
			ABD.add_to_mesh (*this);                            }
		else // with quadrilaterals
		{	Cell Q ( tag::rectangle, AB, BC, CD, DA );  // create a new rectangle
			Q.add_to_mesh (*this);                     }
		*it = CD.reverse();
		it++;  assert ( it == horizon.end() );
	} // end of for i
	// last row of rectangles is different, north sides already exist
	std::list<Cell>::iterator it = horizon.begin();
	Cell DA = west.cell_in_front_of ( NW, tag::surely_exists );
	Cell D = NW;
	for (size_t j=1; j < N_horiz; j++)
	{	Cell AB = *it;
		Cell B = AB.tip();
		Cell CD = north.cell_behind ( D );
		Cell C = CD.base().reverse();
		Cell BC ( tag::segment, B.reverse(), C );  // create a new segment
		BC.spin() = - CD.spin() - DA.spin() - AB.spin();
		if ( cut_rectangles_in_half )
		{	Cell BD ( tag::segment, B.reverse(), D );  // create a new segment
			BD.spin() = - DA.spin() - AB.spin();
			Cell BCD ( tag::triangle, BD.reverse(), BC, CD );  // create a new triangle
			Cell ABD ( tag::triangle, BD, DA, AB );  // create a new triangle
			BCD.add_to_mesh (*this);  // 'this' is the mesh we are building
			ABD.add_to_mesh (*this);                          }
		else // with quadrilaterals
		{	Cell Q ( tag::rectangle, AB, BC, CD, DA );  // create a new rectangle
			Q.add_to_mesh (*this);                       }
		it++;
		D = C;
		DA = BC.reverse();                                          }
	// and the last rectangle of the last row
	Cell AB = *it;
	Cell B = AB.tip();
	Cell BC = east.cell_in_front_of (B);
	Cell C = BC.tip();
	assert ( C == NE );
	Cell CD = north.cell_behind ( D );
	assert ( AB.spin() + BC.spin() + CD.spin() + DA.spin() == 0 );
	if ( cut_rectangles_in_half )
	{	Cell BD ( tag::segment, B.reverse(), D );  // create a new segment
		BD.spin() = BC.spin() + CD.spin();
		Cell BCD ( tag::triangle, BD.reverse(), BC, CD );
		Cell ABD ( tag::triangle, BD, DA, AB );
		BCD.add_to_mesh (*this);  // 'this' is the mesh we are building
		ABD.add_to_mesh (*this);                          }
	else // with quadrilaterals
	{	Cell Q ( tag::rectangle, AB, BC, CD, DA );
		Q.add_to_mesh (*this);                     }
	it++;  assert ( it == horizon.end() );

} // end of Mesh::build with tag::quadrangle and tag::spin

//----------------------------------------------------------------------------------//


Mesh Mesh::fold ( const tag::Identify &, const Mesh & msh1,
                      const tag::With &, const Mesh & msh2,
                  const tag::BuildNewVertices &            )

// take a mesh having as external boundary a parallelogram
// and identify one pair of opposite sides

// a quotient manifold with one action generator will be built
	
{	// we use the current (Euclidian) manifold
	Manifold space = Manifold::working;
	assert ( space.exists() );
	Function coord = space.coordinates();
	assert ( coord.nb_of_components() == 2 );
	// we split 'coord' into its components
	Function x = coord[0], y = coord[1];

	// first we need to identify a translation which moves msh1 into msh2

	Cell A = msh1.first_vertex();
	Cell B = msh1.last_vertex();
	Cell C = msh2.first_vertex();
	Cell D = msh2.last_vertex();

	double dx = x(D) - x(A), dy = y(D) - y(A);
	double norm = std::sqrt ( dx*dx + dy*dy );
	assert ( std::abs ( dx - ( x(C) - x(B) ) ) < 1.e-4 * norm );
	assert ( std::abs ( dy - ( y(C) - y(B) ) ) < 1.e-4 * norm );

	// the desired translation is ( dx, dy )
	Function::Action g ( tag::transforms, coord, tag::into, (x+dx) && ( y+dy) );
	Manifold manif_q = space.quotient ( g );

	// we use a map -- for a faster code, we could use Cell::Core::hook
	std::map < Cell, std::pair < Cell, Function::CompositionOfActions > > corresp_ver;

	assert ( msh1.number_of ( tag::segments ) == msh2.number_of ( tag::segments ) );
	CellIterator it1 = msh1.iterator ( tag::over_vertices, tag::require_order );
	CellIterator it2 = msh2.iterator ( tag::over_vertices, tag::backwards );
	for ( it1.reset(), it2.reset(); it1.in_range(); it1++, it2++ )
	{	assert ( it2.in_range() );
		Cell V = *it1;  Cell W = *it2;
		assert ( std::abs ( dx - ( x(W) - x(V) ) ) < 1.e-4 * norm );
		assert ( std::abs ( dy - ( y(W) - y(V) ) ) < 1.e-4 * norm );
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Function::CompositionOfActions > >
			::iterator it_map = corresp_ver.lower_bound ( W );
		assert ( ( it_map == corresp_ver.end() ) or
		         ( corresp_ver.key_comp()(W,it_map->first) ) );
		corresp_ver.emplace_hint ( it_map, std::piecewise_construct,
		    std::forward_as_tuple ( W ), std::forward_as_tuple
		    ( std::pair < Cell, Function::CompositionOfActions > { V, g } ) );      }
		// corresp_ver [ W ] = { V, g };
	assert ( not it2.in_range() );
						 
	CellIterator it_ver = this->iterator ( tag::over_vertices );
	for ( it_ver.reset(); it_ver.in_range(); it_ver++ )
	{	Cell V = *it_ver;
		if ( V.belongs_to ( msh2 ) ) continue;
		Cell new_V ( tag::vertex );
		x ( new_V ) = x ( V );   y ( new_V ) = y ( V );
		std::pair < Cell, Function::CompositionOfActions > p { new_V, 0 };
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Function::CompositionOfActions > >
			::iterator it_map = corresp_ver.lower_bound ( V );
		assert ( ( it_map == corresp_ver.end() ) or
		         ( corresp_ver.key_comp()(V,it_map->first) ) );
		corresp_ver.emplace_hint ( it_map, std::piecewise_construct,
		    std::forward_as_tuple ( V ), std::forward_as_tuple
		    ( std::pair < Cell, Function::CompositionOfActions > { new_V, 0 } ) );  }
		// corresp_ver [ V ] = p;

	std::map < Cell, Cell > corresp_seg;
	CellIterator it_seg = this->iterator ( tag::over_segments );
	for ( it_seg.reset(); it_seg.in_range(); it_seg++ )
	{	Cell seg = *it_seg;
		std::map < Cell, std::pair < Cell, Function::CompositionOfActions > >
			::iterator it_base_rev = corresp_ver.find ( seg.base().reverse() );
		assert ( it_base_rev != corresp_ver.end() );
		std::map < Cell, std::pair < Cell, Function::CompositionOfActions > >
			::iterator it_tip = corresp_ver.find ( seg.tip() );
		assert ( it_tip != corresp_ver.end() );
		Cell new_seg ( tag::segment,
		               it_base_rev->second.first.reverse(), it_tip->second.first );
		new_seg.spin() = it_tip->second.second - it_base_rev->second.second;
		// new_seg.spin() = corresp_ver [ seg.tip()            ] .second -
		//	                corresp_ver [ seg.base().reverse() ] .second  ;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, Cell > ::iterator it_map = corresp_seg.lower_bound ( seg );
		assert ( ( it_map == corresp_seg.end() ) or
		         ( corresp_seg.key_comp()(seg,it_map->first) ) );
		corresp_seg.emplace_hint ( it_map, std::piecewise_construct,
		     std::forward_as_tuple ( seg ), std::forward_as_tuple ( new_seg ) );     }
		// corresp_seg [ seg ] = new_seg;               

	Mesh result ( tag::fuzzy, tag::of_dim, 2 );
	
	CellIterator it_cll = this->iterator ( tag::over_cells_of_dim, 2 );
	for ( it_cll.reset(); it_cll.in_range(); it_cll++ )
	{	Cell cll = *it_cll;
		Cell::Positive::HighDim * new_cll_ptr = new Cell::Positive::HighDim
			( tag::whose_boundary_is,
				Mesh ( tag::whose_core_is,
			         new Mesh::Connected::OneDim ( tag::with,
			           cll.boundary().number_of ( tag::segments ),
	               tag::segments, tag::one_dummy_wrapper       ),
	         tag::freshly_created                                 ),
				tag::one_dummy_wrapper                                     );
		Cell new_cll ( tag::whose_core_is, new_cll_ptr, tag::freshly_created );
		CellIterator it_bdry = cll.boundary() .iterator ( tag::over_segments );
		for ( it_bdry.reset(); it_bdry.in_range(); it_bdry++ )
		{	Cell seg = *it_bdry;
			// no need for glue_on_bdry_of :
			// new_cll has just been created, it has no meshes above
			// we call add_to_mesh instead (as if the mesh were not a boundary)
			corresp_seg [ seg ] .core->add_to_mesh
				( new_cll.boundary().core, tag::do_not_bother );  }
		assert ( not it_bdry.in_range() );
		it_bdry.reset();  assert ( it_bdry.in_range() );
		Cell seg = *it_bdry;
		Cell V = seg.tip();
		Mesh::Connected::OneDim * bdry = tag::Util::assert_cast
			< Mesh::Core*, Mesh::Connected::OneDim* > ( new_cll.boundary().core );
		bdry->first_ver = V.reverse();  // negative
		bdry->last_ver = V;
		new_cll.add_to_mesh ( result );
	} 

	return result;
	
}

//----------------------------------------------------------------------------------//


Mesh Mesh::fold ( const tag::Identify &, const Mesh & msh1,
                      const tag::With &, const Mesh & msh2,
                  const tag::Identify &, const Mesh & msh3,
                      const tag::With &, const Mesh & msh4,
                  const tag::BuildNewVertices &            )

// take a mesh having as external boundary a parallelogram
// and identify two pairs of opposite sides

// a quotient manifold with two action generators will be built
	
{	// we use the current manifold
	Manifold space = Manifold::working;
	assert ( space.exists() );
	Function coord = space.coordinates();
	assert ( coord.nb_of_components() == 2 );
	// we split 'coord' into its components
	Function x = coord[0],  y = coord[1];

	// we need to identify two translations, one which moves msh1 into msh2
	// and another one which moves msh3 into msh4


	// we use a map
	// for a faster code, we could use Cell::Core::hook
	std::map < Cell, Cell > corresp;

}

//----------------------------------------------------------------------------------//


// gs -q -dNOPAUSE -dBATCH -sDEVICE=png16 -sOUTPUTFILE=out.png in.ps


void Mesh::draw_ps ( std::string file_name )
	
{	// we use the current manifold
	Manifold space = Manifold::working;
	assert ( space.exists() );
	Function coord = space.coordinates();
	assert ( coord.nb_of_components() == 2 );
	// we split 'coord' into its components
	Function x = coord[0],  y = coord[1];

	double xmin, xmax, ymin, ymax, maxside;

	{ // just a block for hiding variables
	CellIterator it = this->iterator ( tag::over_vertices );
	it.reset();
	assert( it.in_range() );
	Cell Vfirst = *it;
	xmin = xmax = x(Vfirst);
	ymin = ymax = y(Vfirst);	
	for ( it++ ; it.in_range(); it++ )
	{ Cell V = *it; 
		double xV = x(V), yV = y(V);
		if ( xV < xmin ) xmin = xV;
	  if ( xV > xmax ) xmax = xV;
	  if ( yV < ymin ) ymin = yV;
	  if ( yV > ymax ) ymax = yV;      }
	} // just a block for hiding variables 
	if ( xmax-xmin < ymax-ymin ) maxside = ymax-ymin;
	else maxside = xmax-xmin;
	double border = 0.02*maxside;
	double scale_factor = 500/maxside;
	double translation_x = -scale_factor*xmin;
	double translation_y = -scale_factor*ymin;

	std::ofstream file_ps ( file_name );
	file_ps << "%!PS-Adobe-3.0 EPSF-3.0" << std::endl;
	file_ps << "%%Title:                     maniFEM" << std::endl;
	file_ps << "%%BoundingBox:  0 0 " << " " << scale_factor*(xmax-xmin+2*border)
	        << "   " << scale_factor*(ymax-ymin+2*border) << std::endl;
	file_ps << "%%EndComments" << std::endl;
	file_ps << "%%BeginSetup" << std::endl;
	file_ps << "<< /PageSize [" << scale_factor*(xmax-xmin+2*border) << " "
	        << scale_factor*(ymax-ymin+2*border) << "] >> setpagedevice" << std::endl;
	file_ps << "%%EndSetup" << std::endl << std::endl;
	
	file_ps << "gsave" << std::endl;
	// file_ps << "/m{moveto}def" << std::endl;
	// file_ps << "/l{lineto}def" << std::endl;
	// file_ps << "/s{stroke}def" << std::endl;
	file_ps << translation_x + scale_factor*border << " "
	        << translation_y + scale_factor*border << " translate" << std::endl;
	file_ps << scale_factor << " dup scale" << std::endl << std::endl;

	file_ps << "gsave " << 1.5 / scale_factor << " setlinewidth" << std::endl;
	
	{ // just a block for hiding variables
	CellIterator it = this->iterator ( tag::over_segments );
	for ( it.reset() ; it.in_range(); it++ )
	{	Cell seg = *it;
		Cell base = seg.base().reverse();
		Cell tip  = seg.tip();
		file_ps << x(base) << " " << y(base) << " moveto" << std::endl;
		file_ps << x(tip) << " " << y(tip) << " lineto stroke" << std::endl;  }
	} // just a block for hiding variables
																					
	file_ps << "grestore" << std::endl;
	file_ps << "grestore" << std::endl << std::endl;
	file_ps << "showpage" << std::endl;
	file_ps << "%%Trailer" << std::endl;
	file_ps << "%EOF" << std::endl;

	if ( ! file_ps.good() )
	{	std::cerr << "error writing postscript file" << std::endl;
		exit (1);                                                   }

} // end of  Mesh::draw_ps

//----------------------------------------------------------------------------------//


void Mesh::draw_ps ( std::string file_name, const tag::Unfold &,
                     const tag::OverRegion &, const Function::Inequality::Set & constraints )

// we draw several translations (or, more generally, transformations)
// of each segment, within the region described by the constraints
	
{	Manifold space = Manifold::working;
	assert ( space.exists() );  // we use the current (quotient) manifold
	Manifold::Quotient * manif_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space.core );
	assert ( manif_q );

	// the action group may have one or two generators
	size_t n = manif_q->actions.size();
	assert ( n == manif_q->spins.size() );
	if ( n == 1 ) this->draw_ps ( file_name, tag::unfold, tag::one_generator,
	                              tag::over_region, constraints               );
	else
	{	assert ( n == 2 );
		this->draw_ps ( file_name, tag::unfold, tag::two_generators,
		                tag::over_region, constraints               );  }

} // end of  Mesh::draw_ps with tag::unfold
	
//----------------------------------------------------------------------------------//


void Mesh::draw_ps ( std::string file_name, const tag::Unfold &, const tag::OneGenerator &,
                     const tag::OverRegion &, const Function::Inequality::Set & constraints )

{	assert ( false );  }

//----------------------------------------------------------------------------------//


void Mesh::draw_ps ( std::string file_name, const tag::Unfold &, const tag::TwoGenerators &,
                     const tag::OverRegion &, const Function::Inequality::Set & constraints )

// we draw several translations (or, more generally, transformations)
// of each segment, within the region described by the constraints
	
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

	// the action group may have one or two generators
	assert ( manif_q->actions.size() == 2 );
	assert ( manif_q->spins.size() == 2 );
	Function::Action g1 = manif_q->actions[0], g2 = manif_q->actions[1];
	
	std::ofstream file_ps ( file_name );
	file_ps << "please copy here the preamble from the end of file - after %EOF " << std::endl;
	file_ps	<< "you can move the shadow by simply changing the initial point" << std::endl;
	file_ps	<< "you may also want to define a different contour path" << std::endl;
	file_ps	<< "please erase the preamble from the end of file; erase also these four lines";
	file_ps << std::endl << std::endl;
						 
	double xmin = 1.e8, xmax = -1.e8, ymin = 1.e8, ymax = -1.e8, maxside;

	Cell shadow ( tag::vertex );
	std::vector < double > coords_base, coords_tip;
	std::vector < std::vector < short int > > directions
		{ { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 } };
	// declare global, here and in progressive.cpp Manifold::Type::Quotient::sq_dist
	// and perhaps in other places

	CellIterator it = this->iterator ( tag::over_segments );
  size_t first_unsuccessful_tries = 1, last_unsuccessful_tries = 0;
	for ( it.reset() ; it.in_range(); it++ )
	{	Cell seg = *it;
		Cell base = seg.base().reverse();
		Cell tip  = seg.tip();
		// we describe a sort of spiral
		// if the first tries are out of the region, we give up after 50 unsuccsessful rounds
		// at the end, we stop after 10 unsuccessful rounds
		size_t size_of_round = 0;
		short int ii = 0, jj = 0;
		while ( true )
		{	size_of_round++;
			bool successful_round = false;
			for ( size_t d = 0; d < 4; d++ )
			{	if ( d == 2 ) size_of_round++;
				for ( size_t i = 0; i < size_of_round; i++ )
				{	Function::CompositionOfActions a = ii*g1 + jj*g2;
					bool touches_region = false;
					coords_base = coords_q ( base, tag::spin, a );
					coords_Eu ( shadow ) = coords_base;
					touches_region = touches_region or constraints.on_cell ( shadow );
					a += seg.spin();
					coords_tip = coords_q ( tip, tag::spin, a );
					coords_Eu ( shadow ) = coords_tip;
					touches_region = touches_region or constraints.on_cell ( shadow );
					if ( touches_region )
					{	successful_round = true;
						double xx = coords_base[0], yy = coords_base[1];
						file_ps << xx << " " << yy << " moveto" << std::endl;
						if ( xx < xmin ) xmin = xx;
						if ( xx > xmax ) xmax = xx;
						if ( yy < ymin ) ymin = yy;
						if ( yy > ymax ) ymax = yy;
						xx = coords_tip[0];  yy = coords_tip[1];
						file_ps << xx << " "  << yy << " lineto stroke" << std::endl;
						if ( xx < xmin ) xmin = xx;
						if ( xx > xmax ) xmax = xx;
						if ( yy < ymin ) ymin = yy;
						if ( yy > ymax ) ymax = yy;                                           }  
					ii += directions[d][0];
					jj += directions[d][1];                                                   }  }
			if ( successful_round )
			{	first_unsuccessful_tries = 0;
				last_unsuccessful_tries = 0;   }
			else  // either we have not started yet, or we are approaching the end
			{	if ( first_unsuccessful_tries > 0 )  // not started yet
				{	first_unsuccessful_tries++;
					if ( first_unsuccessful_tries > 50 ) goto give_up;  }
				else  // approaching end
				{	last_unsuccessful_tries++;
					if ( last_unsuccessful_tries > 10 ) goto give_up;   }   }
		}  // while true
	give_up: ;  // next segment, please
	}  // end of loop over segments

	if ( xmax-xmin < ymax-ymin ) maxside = ymax-ymin;
	else maxside = xmax-xmin;
	double border = 0.02*maxside;
	double scale_factor = 500/maxside;
	double translation_x = -scale_factor*xmin;
	double translation_y = -scale_factor*ymin;

	file_ps << std::endl << "0.5 setgray 0.03 setlinewidth" << std::endl;
	file_ps << "newpath contour stroke" << std::endl << std::endl;
	file_ps << "grestore" << std::endl << std::endl;
	file_ps << "showpage" << std::endl;
	file_ps << "%%Trailer" << std::endl;
	file_ps << "%EOF" << std::endl << std::endl;

	file_ps << "%!PS-Adobe-3.0 EPSF-3.0" << std::endl;
	file_ps << "%%Title:                     maniFEM" << std::endl;
	file_ps << "%%BoundingBox:  0 0 " << " " << scale_factor*(xmax-xmin+2*border)
	        << "   " << scale_factor*(ymax-ymin+2*border) << std::endl;
	file_ps << "%%EndComments" << std::endl;
	file_ps << "%%BeginSetup" << std::endl;
	file_ps << "<< /PageSize [" << scale_factor*(xmax-xmin+2*border) << " "
	        << scale_factor*(ymax-ymin+2*border) << "] >> setpagedevice" << std::endl;
	file_ps << "%%EndSetup" << std::endl << std::endl;
													
	file_ps << "/contour" << std::endl << "{ ";
	file_ps << "/xmin {" << xmin << "} def  /xmax {" << xmax
          << "} def  /ymin {" << ymin << "} def  /ymax {"
          << ymax << "} def  /border {" << border << "}  def" << std::endl;
	file_ps << "xmin border add ymin border add border 180 270 arc" << std::endl;
	file_ps << "xmax border sub ymin lineto" << std::endl;
	file_ps << "xmax border sub ymin border add border -90 0 arc" << std::endl;
	file_ps << "xmax ymax border sub lineto" << std::endl;
	file_ps << "xmax border sub ymax border sub border 0 90 arc" << std::endl;
	file_ps << "xmin border add ymax lineto" << std::endl;
	file_ps << "xmin border add ymax border sub border 90 180 arc" << std::endl;
	file_ps << "xmin ymin border add lineto" << std::endl;
	file_ps << "} def" << std::endl << std::endl;

	file_ps << "/shadow" << std::endl << "{ 0 0 moveto" << std::endl << "  ";
	coords_Eu ( shadow ) = { 0., 0. };
	coords_base = coords_q ( shadow, tag::spin, g1 );
	file_ps << coords_base[0] << " " << coords_base[1] << " rlineto ";
	coords_base = coords_q ( shadow, tag::spin, g2 );
	file_ps << coords_base[0] << " " << coords_base[1] << " rlineto ";
	coords_base = coords_q ( shadow, tag::spin, -g1 );
	file_ps << coords_base[0] << " " << coords_base[1] << " rlineto ";
	file_ps << "} def" << std::endl << std::endl;

	file_ps << "gsave" << std::endl;
	// file_ps << "/m{moveto}def" << std::endl;
	// file_ps << "/l{lineto}def" << std::endl;
	// file_ps << "/s{stroke}def" << std::endl;
	file_ps << translation_x + scale_factor*border << " "
	        << translation_y + scale_factor*border << " translate" << std::endl;
	file_ps << scale_factor << " dup scale" << std::endl << std::endl;

  file_ps << "0 setgray contour clip" << std::endl;
  file_ps << "newpath 0.6 0.8 1. setrgbcolor shadow fill" << std::endl;

  file_ps << 1.5 / scale_factor << " setlinewidth 0 setgray" << std::endl;
	
	if ( ! file_ps.good() )
	{	std::cerr << "error writing postscript file" << std::endl;
		exit (1);                                                   }

} // end of  Mesh::draw_ps with tag::unfold and tag::two_generators

//----------------------------------------------------------------------------------//


void Mesh::draw_ps ( std::string file_name,
         const tag::Unfold &, const std::vector < Function::CompositionOfActions > & v,
         const tag::OverRegion &, const Function::Inequality::Set & constraints        )

// we draw several translations (or, more generally, transformations)
// of each segment, within the region described by the constraints
	
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

	std::ofstream file_ps ( file_name );
	file_ps << "please copy here the preamble from the end of file - after %EOF " << std::endl;
	file_ps	<< "you can move the shadow by simply changing the initial point" << std::endl;
	file_ps	<< "you may also want to define a different contour path" << std::endl;
	file_ps	<< "please erase the preamble from the end of file; erase also these four lines";
	file_ps << std::endl << std::endl;
						 
	double xmin = 1.e8, xmax = -1.e8, ymin = 1.e8, ymax = -1.e8, maxside;

	Cell shadow ( tag::vertex );
	std::vector < double > coords_base, coords_tip;

	CellIterator it = this->iterator ( tag::over_segments );
	for ( it.reset() ; it.in_range(); it++ )
	{	Cell seg = *it;
		Cell base = seg.base().reverse();
		Cell tip  = seg.tip();
		for ( std::vector < Function::CompositionOfActions > ::const_iterator
						it_v = v.begin(); it_v != v.end(); it_v++                     )
		{	Function::CompositionOfActions a = *it_v;
			bool touches_region = false;
			coords_base = coords_q ( base, tag::spin, a );
			coords_Eu ( shadow ) = coords_base;
			touches_region = touches_region or constraints.on_cell ( shadow );
			a += seg.spin();
			coords_tip = coords_q ( tip, tag::spin, a );
			coords_Eu ( shadow ) = coords_tip;
			touches_region = touches_region or constraints.on_cell ( shadow );
			if ( touches_region )
			{	double xx = coords_base[0], yy = coords_base[1];
				file_ps << xx << " " << yy << " moveto" << std::endl;
				if ( xx < xmin ) xmin = xx;
				if ( xx > xmax ) xmax = xx;
				if ( yy < ymin ) ymin = yy;
				if ( yy > ymax ) ymax = yy;
				xx = coords_tip[0];  yy = coords_tip[1];
				file_ps << xx << " "  << yy << " lineto stroke" << std::endl;
				if ( xx < xmin ) xmin = xx;
				if ( xx > xmax ) xmax = xx;
				if ( yy < ymin ) ymin = yy;
				if ( yy > ymax ) ymax = yy;                                           }  }  }

	if ( xmax-xmin < ymax-ymin ) maxside = ymax-ymin;
	else maxside = xmax-xmin;
	double border = 0.02*maxside;
	double scale_factor = 500/maxside;
	double translation_x = -scale_factor*xmin;
	double translation_y = -scale_factor*ymin;

	file_ps << std::endl << "0.5 setgray 0.03 setlinewidth" << std::endl;
	file_ps << "newpath contour stroke" << std::endl << std::endl;
	file_ps << "grestore" << std::endl << std::endl;
	file_ps << "showpage" << std::endl;
	file_ps << "%%Trailer" << std::endl;
	file_ps << "%EOF" << std::endl << std::endl;

	file_ps << "%!PS-Adobe-3.0 EPSF-3.0" << std::endl;
	file_ps << "%%Title:                     maniFEM" << std::endl;
	file_ps << "%%BoundingBox:  0 0 " << " " << scale_factor*(xmax-xmin+2*border)
	        << "   " << scale_factor*(ymax-ymin+2*border) << std::endl;
	file_ps << "%%EndComments" << std::endl;
	file_ps << "%%BeginSetup" << std::endl;
	file_ps << "<< /PageSize [" << scale_factor*(xmax-xmin+2*border) << " "
	        << scale_factor*(ymax-ymin+2*border) << "] >> setpagedevice" << std::endl;
	file_ps << "%%EndSetup" << std::endl << std::endl;
													
	file_ps << "/contour" << std::endl << "{ ";
	file_ps << "/xmin {" << xmin << "} def  /xmax {" << xmax
          << "} def  /ymin {" << ymin << "} def  /ymax {"
          << ymax << "} def  /border {" << border << "}  def" << std::endl;
	file_ps << "xmin border add ymin border add border 180 270 arc" << std::endl;
	file_ps << "xmax border sub ymin lineto" << std::endl;
	file_ps << "xmax border sub ymin border add border -90 0 arc" << std::endl;
	file_ps << "xmax ymax border sub lineto" << std::endl;
	file_ps << "xmax border sub ymax border sub border 0 90 arc" << std::endl;
	file_ps << "xmin border add ymax lineto" << std::endl;
	file_ps << "xmin border add ymax border sub border 90 180 arc" << std::endl;
	file_ps << "xmin ymin border add lineto" << std::endl;
	file_ps << "} def" << std::endl << std::endl;

	file_ps << "/shadow" << std::endl << "{ 0 0 moveto" << std::endl << "  ";
	coords_Eu ( shadow ) = { 0., 0. };
	coords_base = coords_q ( shadow, tag::spin, 0 );
	file_ps << coords_base[0] << " " << coords_base[1] << " rlineto ";
	coords_base = coords_q ( shadow, tag::spin, 0 );
	file_ps << coords_base[0] << " " << coords_base[1] << " rlineto ";
	coords_base = coords_q ( shadow, tag::spin, 0 );
	file_ps << coords_base[0] << " " << coords_base[1] << " rlineto ";
	file_ps << "} def" << std::endl << std::endl;

	file_ps << "gsave" << std::endl;
	// file_ps << "/m{moveto}def" << std::endl;
	// file_ps << "/l{lineto}def" << std::endl;
	// file_ps << "/s{stroke}def" << std::endl;
	file_ps << translation_x + scale_factor*border << " "
	        << translation_y + scale_factor*border << " translate" << std::endl;
	file_ps << scale_factor << " dup scale" << std::endl << std::endl;

  file_ps << "0 setgray contour clip" << std::endl;
  file_ps << "newpath 0.6 0.8 1. setrgbcolor shadow fill" << std::endl;

  file_ps << 1.5 / scale_factor << " setlinewidth 0 setgray" << std::endl;
	
	if ( ! file_ps.good() )
	{	std::cerr << "error writing postscript file" << std::endl;
		exit (1);                                                   }

} // end of  Mesh::draw_ps with tag::unfold and tag::two_generators

//----------------------------------------------------------------------------------//


//  method below relies on some postscript macros for (very) rudimentary 3d drawings
//  available at https://github.com/cristian-barbarosie/manifem/blob/main/3d.ps

void Mesh::draw_ps_3d ( std::string file_name )
	
{	// we use the current manifold
	Manifold space = Manifold::working;
	assert ( space.exists() );
	Function coord = space.coordinates();
	assert ( coord.nb_of_components() == 3 );
	// we split 'coord' into its components
	Function x = coord[0],  y = coord[1],  z = coord[2];

	double xmin, xmax, ymin, ymax, zmin, zmax, maxside;
	
	{ // just a block for hiding variables
	CellIterator it = this->iterator ( tag::over_vertices );
	it.reset();  assert( it.in_range() );
	Cell Vfirst = *it;
	xmin = xmax = x(Vfirst);
	ymin = ymax = y(Vfirst);	
	zmin = zmax = z(Vfirst);	
	for ( it++ ; it.in_range(); it++ )
	{	Cell V = *it; 
		double xV = x(V), yV = y(V), zV = z(V);
		if ( xV < xmin ) xmin = xV;
		if ( xV > xmax ) xmax = xV;
		if ( yV < ymin ) ymin = yV;
		if ( yV > ymax ) ymax = yV;
		if ( zV < zmin ) zmin = zV;
		if ( zV > zmax ) zmax = zV;      }
	} // just a block for hiding variables 
	// we look at the object along the y axis, so values of y do not count
	if ( xmax-xmin < zmax-zmin ) maxside = zmax-zmin;
	else maxside = xmax-xmin;
	double border = 0.02*maxside;
	double scale_factor = 500/maxside;
	double translation_x = -scale_factor*xmin;
	double translation_y = -scale_factor*zmin;

	std::ofstream file_ps ( file_name );

	file_ps << "%!PS-Adobe-3.0 EPSF-3.0" << std::endl;
	file_ps << "%%Title:                     maniFEM" << std::endl;
	file_ps << "%%BoundingBox:  0 0 " << " " << scale_factor*(xmax-xmin+2*border)
	        << "   " << scale_factor*(zmax-zmin+2*border) << std::endl;
	file_ps << "%%EndComments" << std::endl;
	file_ps << "%%BeginSetup" << std::endl;
	file_ps << "<< /PageSize [" << scale_factor*(xmax-xmin+2*border) << " "
	        << scale_factor*(zmax-zmin+2*border) << "] >> setpagedevice" << std::endl;
	file_ps << "%%EndSetup" << std::endl << std::endl;

	std::ifstream file_3d ("3d.ps");
	while ( true )
	{	std::string line;  // this way 'line' remains local
		if ( getline ( file_3d, line ) ) file_ps << line + '\n';
		else break;                                               }

	file_ps << "5 rotxy" << std::endl << "-8 rotyz" << std::endl
          << "0.7 rotxz" << std::endl << std::endl;

	file_ps << "gsave" << std::endl;
	file_ps << translation_x + scale_factor*border << " "
	        << translation_y + scale_factor*border << " translate" << std::endl;
	file_ps << scale_factor << " dup scale" << std::endl << std::endl;

	file_ps << 1. / scale_factor << " setlinewidth" << std::endl;
	file_ps << xmin << " " << ymin << " " << zmin << " " << " proj moveto ";
	file_ps << xmax << " " << ymin << " " << zmin
          << " proj Lineto^ stroke" << std::endl;
	file_ps << xmin << " " << ymin << " " << zmin << " " << " proj moveto ";
	file_ps << xmin << " " << ymax << " " << zmin
          << " proj Lineto^ stroke" << std::endl;
	file_ps << xmin << " " << ymin << " " << zmin << " " << " proj moveto ";
	file_ps << xmin << " " << ymin << " " << zmax
          << " proj Lineto^ stroke" << std::endl << std::endl;
	
	file_ps << "gsave " << 1.5 / scale_factor << " setlinewidth" << std::endl;
	
	{ // just a block for hiding variables
	CellIterator it = this->iterator ( tag::over_segments );
	for ( it.reset() ; it.in_range(); it++ )
	{	Cell seg = *it;
		Cell base = seg.base().reverse();
		Cell tip  = seg.tip();
		file_ps << x(base) << " " << y(base) << " " << z(base)
	          << " proj moveto" << std::endl;
		file_ps << x(tip) << " " << y(tip) << " " << z(tip)
	          << " proj lineto stroke" << std::endl;           }
	} // just a block for hiding variables
																					
	file_ps << "grestore" << std::endl;
	file_ps << "grestore" << std::endl << std::endl;
	file_ps << "showpage" << std::endl;
	file_ps << "%%Trailer" << std::endl;
	file_ps << "%EOF" << std::endl;

	if ( ! file_ps.good() )
	{	std::cerr << "error writing postscript file" << std::endl;
		exit (1);                                                   }

} // end of  Mesh::draw_ps_3d

//----------------------------------------------------------------------------------//


void Mesh::export_msh ( std::string f, Cell::Numbering & ver_numbering )

{	// we use the current manifold
	Manifold space = Manifold::working;
	assert ( space.exists() );
	Function coord = space.coordinates();

	std::ofstream file_msh (f);
	file_msh << "$MeshFormat" << std::endl << "2.2 0 8" << std::endl;
	file_msh << "$EndMeshFormat" << std::endl;
	
	file_msh << "$Nodes" << std::endl << this->number_of(tag::cells_of_dim,0) << std::endl;

	{ // just to make variables local : it, counter, x, y
	CellIterator it = this->iterator ( tag::over_vertices );
	Function x = coord[0], y = coord[1];
	if (coord.nb_of_components() == 2)
	{	for ( it.reset() ; it.in_range(); it++ )
		{	Cell p = *it;
			file_msh << ver_numbering (p) << " "
		           << x(p) << " " << y(p) << " " << 0 << std::endl;  }  }
	else
	{	assert  ( coord.nb_of_components() == 3 );
		Function z = coord[2];
		for ( it.reset() ; it.in_range(); it++ )
		{	Cell p = *it;
			file_msh << ver_numbering (p) << " "
		           << x(p) << " " << y(p) << " " << z(p) << std::endl;  }  }
	file_msh << "$EndNodes" << std::endl;
	} // just to make variables local : it, counter, x, y

	file_msh << "$Elements" << std::endl;
	file_msh << this->number_of ( tag::cells_of_dim, this->dim() ) << std::endl;

	if ( this->dim() == 1 )
	{	CellIterator it = this->iterator ( tag::over_segments );
		size_t counter = 0;
		for ( it.reset() ; it.in_range(); it++)
		{	++counter;
			Cell elem = *it;
			file_msh << counter << " 1 0 ";
			Cell A = elem.base().reverse();
			file_msh << ver_numbering (A) << " ";
			Cell B = elem.tip();
			file_msh << ver_numbering (B) << std::endl;    }  }
	else if ( this->dim() == 2 )
	{	CellIterator it = this->iterator ( tag::over_cells_of_dim, 2 );
		size_t counter = 0;
		for ( it.reset() ; it.in_range(); it++)
		{	++counter;
			Cell elem = *it;
			if ( elem.boundary().number_of ( tag::cells_of_dim, 1 ) == 3 ) // a triangle
				file_msh << counter << " 2 0 ";
			else // a quadrilateral
			{	assert ( elem.boundary().number_of ( tag::cells_of_dim, 1 ) == 4 );
				file_msh << counter << " 3 0 ";                                     }
			CellIterator itt = elem.boundary().iterator ( tag::over_vertices, tag::require_order );
			for ( itt.reset(); itt.in_range(); ++itt )
			{	Cell p = *itt;  file_msh << ver_numbering (p) << " ";   }
			file_msh << std::endl;                                                                 }  }
	else
	{	assert ( this->dim() == 3);
		CellIterator it = this->iterator ( tag::over_cells_of_dim, 3 );
		size_t counter = 0;
		for ( it.reset() ; it.in_range(); it++)
		{	++counter;
			Cell elem = *it;
			size_t n_faces = elem.boundary().number_of ( tag::cells_of_dim, 2 );
			if ( n_faces == 4 ) // a tetrahedron
				file_msh << counter << " 4 0 ";  // to finish !
			else if ( n_faces == 6 )
			{	// 3d parallelogram = 8-node hexahedron = cube
				// see http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
				// and http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
				file_msh << counter << " 5 0 ";
				CellIterator itt = elem.boundary().iterator ( tag::over_cells_of_dim, 2 );
				itt.reset();  Cell back = *itt; // square face behind the cube
				// back is 0321 in gmsh's documentation
				assert ( back.boundary().number_of ( tag::cells_of_dim, 1 ) == 4 );
				CellIterator itv = back.boundary().iterator ( tag::over_vertices, tag::backwards );
				// reverse because we want the vertices ordered as 0, 1, 2, 3
				itv.reset();  Cell ver_0 = *itv;
				Cell seg_03 = back.boundary().cell_in_front_of(ver_0);
				for ( ; itv.in_range(); ++itv )
				{	Cell p = *itv;  file_msh << ver_numbering (p) << " ";   }
				Cell left_wall = elem.boundary().cell_in_front_of(seg_03); // square face on the left
				// left_wall is 0473 in gmsh's documentation
				assert ( left_wall.boundary().number_of ( tag::cells_of_dim, 1 ) == 4 );
				Cell seg_04 = left_wall.boundary().cell_in_front_of(ver_0);
				Cell ver_4 = seg_04.tip();
				Cell seg_47 = left_wall.boundary().cell_in_front_of(ver_4);
				Cell front = elem.boundary().cell_in_front_of(seg_47); // square face in front
				// front is 4567 in gmsh's documentation
				assert ( front.boundary().number_of ( tag::cells_of_dim, 1 ) == 4 );
				CellIterator itvv = front.boundary().iterator ( tag::over_vertices, tag::require_order );
				itvv.reset ( tag::start_at, ver_4 );
				for ( ; itvv.in_range(); ++itvv )
				{	Cell p = *itvv;  file_msh << ver_numbering (p) << " ";   }                }
			else
			{	assert( n_faces == 5 );
				// triangular prism = 6-node prism
				// see http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
				// and http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
				file_msh << counter << " 6 0 ";
				CellIterator itt = elem.boundary().iterator ( tag::over_cells_of_dim, 2 );
				size_t n_tri = 0, n_rect = 0;
				Cell base ( tag::non_existent );  // temporary non-existent cell
				for( itt.reset(); itt.in_range(); ++itt )
				{	Cell face = *itt; // every face elem
			    size_t n_edges = face.boundary().number_of ( tag::cells_of_dim, 1 );
					if ( n_edges == 3 )  { n_tri++;   base = face;              }
					else                 { n_rect++;  assert ( n_edges == 4 );  }         }
				assert ( n_tri == 2 );  assert ( n_rect == 3 );
				// base is 021 in gmsh's documentation
				assert ( base.boundary().number_of ( tag::cells_of_dim, 1 ) == 3 );
				CellIterator itv = base.boundary().iterator ( tag::over_vertices, tag::backwards );
				// we set it reverse because we want the vertices ordered as 0, 1, 2
				itv.reset();  Cell ver_0 = *itv;
				Cell seg_02 = base.boundary().cell_in_front_of(ver_0);
				for (  ; itv.in_range(); ++itv )
				{	Cell p = *itv;  file_msh << ver_numbering (p) << " ";   }
				Cell right_wall = elem.boundary().cell_in_front_of(seg_02);
				// right_wall is 0352 in gmsh's documentation
				assert ( right_wall.boundary().number_of ( tag::cells_of_dim, 1 ) == 4 );
				Cell seg_03 = right_wall.boundary().cell_in_front_of(ver_0);
				Cell ver_3 = seg_03.tip();
				Cell seg_35 = right_wall.boundary().cell_in_front_of(ver_3);
				Cell roof = elem.boundary().cell_in_front_of(seg_35);
				// roof is 345 in gmsh's documentation
				assert ( roof.boundary().number_of ( tag::cells_of_dim, 1 ) == 3 );
				CellIterator itvv = roof.boundary().iterator ( tag::over_vertices, tag::require_order );
				itvv.reset ( tag::start_at, ver_3 );
				for ( ; itvv.in_range(); ++itvv )
				{	Cell p = *itvv;  file_msh << ver_numbering (p) << " ";  }               }
			file_msh << std::endl;                                                                } }
	file_msh << "$EndElements" << std::endl;
	
	if ( ! file_msh.good() )
	{	std::cerr << "error writing msh file" << std::endl;
		exit (1);                                    }

} // end of Mesh::export_msh


void Mesh::export_msh ( std::string f, std::map < Cell, size_t > & numb_map )
	
{	Cell::Numbering::Map numbering ( & numb_map );

	this->export_msh ( f, numbering );

} // end of Mesh::export_msh


void Mesh::export_msh ( std::string f )
	
// the numbering of vertices is produced on-the-fly

{	Cell::Numbering::Map numbering;

	CellIterator it = this->iterator ( tag::over_vertices );
	size_t counter = 0;
	for ( it.reset() ; it.in_range(); it++ )
	{	Cell p = *it;  numbering (p) = counter;  ++counter;  }

	this->export_msh ( f, numbering );

} // end of Mesh::export_msh


size_t & Cell::Numbering::Map::operator() ( const Cell p )  // virtual from Cell::Numbering
{	return (*(this->map))[p];  }
