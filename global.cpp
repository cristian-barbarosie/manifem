
// global.cpp 2021.12.26

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
                   const tag::Winding &, const tag::Util::Action & s )

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
	{	double new_co = coords_q ( B, tag::winding, s );
		coords_Eu ( shadow_of_B ) = new_co;            }
	else
	{	assert ( coords_Eu.nb_of_components() > 1 );
		std::vector < double > new_co = coords_q ( B, tag::winding, s );
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
	// size_t nb_spins = mani_q->winding_nbs.size();
	// for ( size_t i = 0; i < nb_spins; i++ )
	// {	Field::ShortInt & sp = mani_q->winding_nbs[i];
	// 	std::map<Function::ActionGenerator,short int>::const_iterator it =
	// 		s.index_map.lower_bound(mani_q->actions[i]);
	// 	sp.on_cell ( seg.core ) = it->second;                     }
	seg.winding() = s;
	seg.add_to_mesh ( *this, tag::do_not_bother );

	space.set_as_working_manifold();

}  // end of Mesh::build with tag::segment and tag::winding

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
	Mesh::Iterator it = AB.iterator ( tag::over_segments, tag::require_order );
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
	Cell seg_on_CA = CA .cell_in_front_of ( C );
	assert ( seg_on_CA == CA .cell_behind ( Q_CA ) );
	assert ( seg_on_BC == BC .cell_behind ( C ) );
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
	Mesh::Iterator it = seg1.iterator ( tag::over_vertices );
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

	
void build_common
( Mesh msh, const Mesh & AB, const Mesh & BC, const Mesh & CA,
            const Cell & A, const Cell & B, const Cell & C,
  Cell & seg1, Cell & seg2, Manifold::Action & winding_C_from_A,
  bool not_singular                                             )

// return two segments seg1 and seg2 and also 'winding_C_from_A'
	
// C may be singular

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
	
	// we keep a list of horizontal segments (parallel to AB)
	// useful for the next layer of triangles
	std::list < Cell > ground, ceiling;
	{ // just a block of code for hiding 'it'
		Mesh::Iterator it = AB.iterator ( tag::over_segments, tag::require_order );
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

  // we keep winding numbers of vertices relative to A
	Manifold::Action winding_B = 0;
	// C may be singular, so we keep two different winding numbers for it
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = AB.iterator ( tag::over_segments );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell seg = *it;  winding_B += seg.winding();  }
	} // just a block of code for hiding 'it'
	Manifold::Action winding_C_via_B = winding_B;
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = BC.iterator ( tag::over_segments );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell seg = *it;  winding_C_via_B += seg.winding();  }
	} { // just a block of code for hiding names
	winding_C_from_A = 0;
	Mesh::Iterator it = CA.iterator ( tag::over_segments );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell seg = *it;  winding_C_from_A -= seg.winding();  }
	if ( not_singular ) assert ( winding_C_from_A == winding_C_via_B );
	} // just a block of code for hiding 'it'
	
	Manifold::Action winding_Q_AB_ini = 0, winding_P_BC = winding_B, winding_Q_CA = 0;

	// we end before i=N-1, meaning triangles will be built
	// except the last four, near C
	// calling code must build those four triangles
	for ( size_t i = 1; i < N-1; i++ ) // "vertical" movement
	{	// advance one level upwards and slightly right (parallel to CA)
		std::list<Cell>::iterator it_ground = ground.begin();
		Cell ground_seg = *it_ground;
		assert ( ground_seg.base().reverse() == Q_CA );
		Cell ground_ver = ground_seg.tip();
		Manifold::Action winding_ground_ver = winding_Q_CA + ground_seg.winding();
		Cell seg_Q_AB = AB.cell_in_front_of ( Q_AB_ini );
		Q_AB_ini = seg_Q_AB.tip();
		winding_Q_AB_ini += seg_Q_AB.winding();
		Cell Q_AB = Q_AB_ini;
		Manifold::Action winding_Q_AB = winding_Q_AB_ini;
		assert ( seg_P_BC == BC.cell_in_front_of ( P_BC ) );
		P_BC = seg_P_BC.tip();
		winding_P_BC += seg_P_BC.winding();
		Cell seg_Q_CA = CA.cell_behind(Q_CA);
		Cell P_CA = Q_CA = seg_Q_CA.base().reverse();
		winding_Q_CA -= seg_Q_CA.winding();
		Manifold::Action winding_P_CA = winding_Q_CA;
		std::vector < double > v = coords_q ( Q_CA, tag::winding, winding_Q_CA );
		coords_Eu ( shadow_Q_CA ) = v;
		v = coords_q ( P_BC, tag::winding, winding_P_BC );
		coords_Eu ( shadow_P_BC ) = v;
		Cell P_AB = A, Q_BC = C;
		Manifold::Action winding_P_AB = 0, winding_Q_BC = winding_C_via_B;
		// build the first triangle on this layer
		Cell previous_seg ( tag::segment, ground_ver.reverse(), P_CA );
		previous_seg.winding() = - ground_seg.winding() - seg_Q_CA.winding();
		Cell tri ( tag::triangle, ground_seg, previous_seg, seg_Q_CA );
		tri.add_to_mesh ( msh );  // 'msh' is the mesh we are building
		Cell previous_ver = Q_CA;
		Manifold::Action winding_prev_ver = winding_Q_CA;
		ceiling.clear();
		for ( size_t j = i+1; j <= N; j++ ) // "horizontal" movement
		{	// advance one step horizontally (parallel to AB)
			Cell seg_P_AB = AB.cell_in_front_of ( P_AB );
			P_AB = seg_P_AB.tip();
			winding_P_AB += seg_P_AB.winding();
			seg_Q_AB = AB.cell_in_front_of ( Q_AB );
			Q_AB = seg_Q_AB.tip();
			winding_Q_AB += seg_Q_AB.winding();
			Cell seg_Q_BC = BC.cell_behind ( Q_BC );
			Q_BC = seg_Q_BC.base().reverse();
			winding_Q_BC -= seg_Q_BC.winding();
			Cell seg_P_CA = CA.cell_behind ( P_CA );
			P_CA = seg_P_CA.base().reverse();
			winding_P_CA -= seg_P_CA.winding();
			Cell S ( tag::non_existent );  // temporary non-existent cell
			Manifold::Action winding_S = 0;
			if ( j == N )  { S = P_BC; winding_S = winding_P_BC;  }
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
				v = coords_q ( P_AB, tag::winding, winding_P_AB );
				coords_Eu ( shadow_P_AB ) = v;
				v = coords_q ( Q_AB, tag::winding, winding_Q_AB );
				coords_Eu ( shadow_Q_AB ) = v;
				v = coords_q ( Q_BC, tag::winding, winding_Q_BC );
				coords_Eu ( shadow_Q_BC ) = v;
				v = coords_q ( P_CA, tag::winding, winding_P_CA );
				coords_Eu ( shadow_P_CA ) = v;
				mani_Eu.interpolate ( S, frac_AB, shadow_P_AB,  frac_AB, shadow_Q_AB,
													       frac_BC, shadow_P_BC,  frac_BC, shadow_Q_BC,
													       frac_CA, shadow_P_CA,  frac_CA, shadow_Q_CA );  }
			Cell new_seg ( tag::segment, ground_ver.reverse(), S );
			new_seg.winding() = winding_S - winding_ground_ver;
			Cell horizontal_seg ( tag::segment, S.reverse(), previous_ver );
			horizontal_seg.winding() = winding_prev_ver - winding_S;
			assert ( horizontal_seg.winding() + new_seg.winding() - previous_seg.winding() == 0 );
		  Cell tri_1 ( tag::triangle, previous_seg.reverse(), new_seg, horizontal_seg );
		  tri_1.add_to_mesh ( msh );  // 'msh' is the mesh we are building
			it_ground++;  assert ( it_ground != ground.end() );
			ground_seg = *it_ground;
			if ( j == N ) previous_seg = seg_P_BC;
			else
			{	ground_ver = ground_seg.tip();
				previous_seg = Cell ( tag::segment, ground_ver.reverse(), S );
				previous_seg.winding() = winding_S - winding_ground_ver;                }
			assert ( ground_seg.winding() + previous_seg.winding() - new_seg.winding() == 0 );
		  Cell tri_2 ( tag::triangle, ground_seg, previous_seg, new_seg.reverse() );
			tri_2.add_to_mesh ( msh );  // 'msh' is the mesh we are building
			previous_ver = S;  winding_prev_ver = winding_S;
			// add horizontal_seg.reverse() to future ground
			ceiling.push_back ( horizontal_seg.reverse() );                              	  }
		assert ( seg_P_BC.tip() == P_BC );
		seg_P_BC = BC.cell_in_front_of ( P_BC );
		ground = ceiling;                                                                    }
		// improve by moving ceiling to ground, leaving ceiling empty !

	// return two segments (elements of 'ground') and also 'winding_C_from_A'
	
	std::list < Cell > ::iterator it = ground .begin();
	assert ( it != ground .end() );
	seg1 = *it;
	it++;
	assert ( it != ground .end() );
	seg2 = *it;
	it++;
	assert ( it == ground .end() );
	
}  // end of  build_common
	
} // end of anonymous namespace


void Mesh::build ( const tag::Triangle &,
                   const Mesh & AB, const Mesh & BC, const Mesh & CA,
                   const tag::Winding &                              )

// see paragraph 12.4 in the manual
// the tag:::winding tells maniFEM that we are on a quotient manifold
// and that the segments provided (AB, BC, CA) may be winding
	
// beware, sides may be closed loops
	
{	// recover corners from the sides
	// the process is different from the one in 'build' without tag::winding
	// here, sides may be closed loops and then methods 'first_vertex' and 'last_vertex'
	// become meaningless
	Cell A = find_common_vertex ( CA, AB );
	Cell B = find_common_vertex ( AB, BC );
	Cell C = find_common_vertex ( BC, CA );

	// 'build_common' builds all triangles but the last four, near C
	Manifold::Action winding_C_from_A;
	Cell seg1 ( tag::non_existent ), seg2 ( tag::non_existent );
	build_common ( *this, AB, BC, CA, A, B, C, seg1, seg2, winding_C_from_A, true );
	// last argument true means "no singularity", perform all checkings

	// build last four triangles
	// 'build_common' provides two segments (elements of 'ground')
	// and also 'winding_C_from_A' (not used here)
	Cell CE = CA .cell_in_front_of ( C );
	Cell E = CE .tip();
	Cell EF = CA. cell_in_front_of ( E );
	Cell F = EF .tip();
	assert ( F == seg1 .base() .reverse() );
	Cell I = seg1 .tip();
	assert ( I == seg2 .base() .reverse() );
	Cell H = seg2 .tip();
	assert ( H .belongs_to ( BC ) );
	Cell HG = BC .cell_in_front_of ( H );
	Cell G = HG .tip();
	Cell GC = BC .cell_in_front_of ( G );
	assert ( GC .tip() == C );
	Cell IE ( tag::segment, I .reverse(), E );
	IE .winding() = - EF .winding() - seg1 .winding();
	Cell GI ( tag::segment, G .reverse(), I );
	GI .winding() = - seg2 .winding() - HG .winding();
	Cell GE ( tag::segment, G .reverse(), E );
	GE .winding() = GI .winding() + IE .winding();
	assert ( GE .winding() == GC .winding() + CE .winding() );
	Cell FIE ( tag::triangle, seg1, IE, EF );
	FIE .add_to_mesh ( *this );
	Cell GIH ( tag::triangle, GI, seg2, HG );
	GIH .add_to_mesh ( *this );
	Cell IGE ( tag::triangle, GI .reverse(), GE, IE .reverse() );
	IGE .add_to_mesh ( *this );
	Cell CEG ( tag::triangle, CE, GE .reverse(), GC );
	CEG .add_to_mesh ( *this );

} // end of Mesh::build with tag::triangle and tag::winding


void Mesh::build ( const tag::Triangle &,
                   const Mesh & A_B, const Mesh & B_C, const Mesh & C_A,
                   const tag::Winding &, const tag::Singular &, const Cell & O )

// see paragraph 12.4 in the manual
// the tag:::winding tells maniFEM that we are on a quotient manifold
// and that the segments provided (AB, BC, CA) may be winding
	
// beware, sides may be closed loops

// tag::singular means Cell O is special, it is like the vertex of a cone
// of the segments provided as arguments has O as base and other as tip
	
{	Manifold space = Manifold::working;
	assert ( space.exists() );  // we use the current manifold
	Manifold::Quotient * mani_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space.core );
	Function coords_q = space.coordinates();
	Manifold mani_Eu = mani_q->base_space;  // underlying Euclidian manifold
	Function coords_Eu = mani_Eu.coordinates();
	
	// we search for the side which does not contain O and call it AB
	// in this new notation, O will be equal to C
	
	Mesh AB = A_B, BC = B_C, CA = C_A;
	if ( not O .belongs_to ( BC ) )  //  A == O
	{	AB = B_C;  BC = C_A;  CA = A_B;  }
	if ( not O .belongs_to ( CA ) )  //  B == O
	{	AB = C_A;  BC = A_B;  CA = B_C;  }

	Cell A = find_common_vertex ( CA, AB );
	Cell B = find_common_vertex ( AB, BC );	

	// 'build_common' builds all triangles but the last four, near C
	// recall that C == O
	Manifold::Action winding_C_from_A;
	Cell seg1 ( tag::non_existent ), seg2 ( tag::non_existent );
	build_common ( *this, AB, BC, CA, A, B, O, seg1, seg2, winding_C_from_A, false );
	// last argument false means there is a singularity, skip some checkings

	// build last four triangles
	// 'build_common' provides two segments (elements of 'ground')
	// and also 'winding_C_from_A'
	// recall that C == O
	Cell OE = CA .cell_in_front_of ( O );
	Cell E = OE .tip();
	Cell EF = CA. cell_in_front_of ( E );
	Cell F = EF .tip();
	assert ( F == seg1 .base() .reverse() );
	Cell I = seg1 .tip();
	assert ( I == seg2 .base() .reverse() );
	Cell H = seg2 .tip();
	assert ( H .belongs_to ( BC ) );
	Cell HG = BC .cell_in_front_of ( H );
	Cell G = HG .tip();
	Cell GO = BC .cell_in_front_of ( G );
	assert ( GO .tip() == O );
	Cell IE ( tag::segment, I .reverse(), E );
	IE .winding() = - EF .winding() - seg1 .winding();
	Cell GI ( tag::segment, G .reverse(), I );
	GI .winding() = - seg2 .winding() - HG .winding();
	Cell FIE ( tag::triangle, seg1, IE, EF );
	FIE .add_to_mesh ( *this );
	Cell GIH ( tag::triangle, GI, seg2, HG );
	GIH .add_to_mesh ( *this );

	// now the last triangle
	// vertices G and E may be one and the same ( if BC == CA .reverse() )
	// recall that C == O
	if ( E == G )  // we must create one more vertex
	{	Manifold::Action winding_E = winding_C_from_A + OE .winding(),
		                 winding_G = winding_E - IE .winding() - GI .winding();
		Cell shadow_E ( tag::vertex ), shadow_G ( tag::vertex );
		std::vector < double > v = coords_q ( E, tag::winding, winding_E );
		coords_Eu ( shadow_E ) = v;
		v = coords_q ( G, tag::winding, winding_G );
		coords_Eu ( shadow_G ) = v;
		Cell S ( tag::vertex );
		mani_Eu.interpolate ( S, 0.5, shadow_E, 0.5, shadow_G );
		Cell ES ( tag::segment, E .reverse(), S );  // no winding
		Cell IS ( tag::segment, I .reverse(), S );
		IS .winding() = IE .winding();
		Cell ISE ( tag::triangle, IS, ES .reverse(), IE .reverse() );
		ISE .add_to_mesh ( *this );
		Cell SG ( tag::segment, S .reverse(), G );
		SG .winding() = - GI .winding() - IS .winding();
		Cell SIG ( tag::triangle, IS .reverse(), GI .reverse(), SG .reverse() );
		SIG .add_to_mesh ( *this );
		Cell OS ( tag::segment, O .reverse(), S );
		OS .winding() = OE .winding();  // ES .winding() == 0
		Cell OES ( tag::triangle, OE, ES, OS .reverse() );
		OES .add_to_mesh ( *this );
		Cell OSG ( tag::triangle, OS, SG, GO );
		OSG .add_to_mesh ( *this );                                    }
		// triangle OSG does not fulfill the condition of zero winding
	else		
	{	Cell GE ( tag::segment, G .reverse(), E );
		GE .winding() = GI .winding() + IE .winding();
		// assert ( GE .winding() == GC .winding() + CE .winding() );
		// assertion above is usually false
		Cell IGE ( tag::triangle, GI .reverse(), GE, IE .reverse() );
		IGE .add_to_mesh ( *this );
		Cell OEG ( tag::triangle, OE, GE .reverse(), GO );
		OEG .add_to_mesh ( *this );                                   }
		// triangle OEG does not fulfill the condition of zero winding

} // end of Mesh::build with tag::triangle and tag::winding and tag::singular

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
	Mesh::Iterator it = south.iterator ( tag::over_segments, tag::require_order );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell seg = *it;  horizon.push_back ( seg );  }
	} // just a block of code for hiding 'it'

	// start mesh generation
	Mesh::Iterator it_east = east.iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it_west = west.iterator ( tag::over_vertices, tag::backwards );
	Mesh::Iterator it_south = south.iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it_north = north.iterator ( tag::over_vertices, tag::backwards );
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
                   const tag::Winding &                                               )

// see paragraph 12.3 in the manual
// the tag:::winding tells maniFEM that we are on a quotient manifold
// and that the segments provided (south, east, north, west) may be winding

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
	// the process is different from the one in 'build' without tag::winding
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
	Mesh::Iterator it = south.iterator ( tag::over_segments, tag::require_order );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell seg = *it;  horizon.push_back ( seg );  }
	} // just a block of code for hiding 'it'

	// we have to deal with possible winding segments
	// we choose that, at each interpolation operation, i.e., for each new vertex,
	// the new vertex will have winding zero relatively to SW
	// we must keep track of the windings of ver_south, ver_east, ver_north, ver_west, B, D
	// (relatively to SW)
	// we use four shadow vertices for interpolation
	Cell shadow_south ( tag::vertex ), shadow_east ( tag::vertex );
	Cell shadow_north ( tag::vertex ), shadow_west ( tag::vertex );
	
	Manifold::Action winding_NW = 0, winding_SE = 0;
	// winding_SW is zero by our choice, winding_NE is not needed
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = south.iterator ( tag::over_segments );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell seg = *it;  winding_SE += seg.winding();  }
	} { // just a block of code for hiding 'it'
	Mesh::Iterator it = west.iterator ( tag::over_segments );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell seg = *it;  winding_NW -= seg.winding();  }
	} // just a block of code for hiding 'it'

	// start mesh generation
	// sides may be closed loops
	// have SW, SE, NW and NE been correctly defined ?
	Cell seg_east = east.cell_in_front_of ( SE, tag::surely_exists );
	Cell seg_west = west.cell_behind ( SW, tag::surely_exists );
	Manifold::Action winding_ver_west = 0;  // winding_SW is zero by our choice
	Manifold::Action winding_ver_east = winding_SE;
	Manifold::Action winding_B = 0;  // winding_SW is zero by our choice
	for ( size_t i = 1; i < N_vert; ++i )
	{	std::list<Cell>::iterator it = horizon.begin();
		Cell AB = *it;
		Cell A = AB.base().reverse();
	  Cell DA = seg_west;
		Cell D = DA.base().reverse();
		winding_ver_east += seg_east.winding();
		Cell ver_east = seg_east.tip();
		seg_east = east.cell_in_front_of ( ver_east, tag::surely_exists );
		winding_ver_west -= seg_west.winding();
		Cell ver_west = seg_west.base().reverse();
		seg_west = west.cell_behind ( ver_west, tag::surely_exists );
		assert ( ver_west == D );
		double frac_N = double(i) / double(N_vert),  alpha = frac_N * (1-frac_N);
		alpha = alpha*alpha*alpha;
		std::vector < double > v = coords_q ( ver_east, tag::winding, winding_ver_east );
		coords_Eu ( shadow_east ) = v;
		v = coords_q ( ver_west, tag::winding, winding_ver_west );
		coords_Eu ( shadow_west ) = v;
		Cell seg_south = south.cell_in_front_of ( SW, tag::surely_exists );
		Cell seg_north = north.cell_behind ( NW, tag::surely_exists );
		Manifold::Action winding_ver_south = 0;  // winding_SW is zero by our choice
		Manifold::Action winding_ver_north = winding_NW;
		Manifold::Action winding_D = winding_ver_west;
		for ( size_t j = 1; j < N_horiz; j++ )
		{	AB = *it;  // 'it' points into the 'horizon' list of segments
			Cell B = AB.tip();
			winding_B += AB.winding();
			winding_ver_south += seg_south.winding();
			Cell ver_south = seg_south.tip();
			seg_south = south.cell_in_front_of ( ver_south );
			winding_ver_north -= seg_north.winding();
			Cell ver_north = seg_north.base().reverse();
			seg_north = north.cell_behind ( ver_north );
			Cell C ( tag::vertex );  // create a new vertex
			double frac_E = double(j) / double(N_horiz),  beta = frac_E * (1-frac_E);
			beta = beta*beta*beta;
			double sum = alpha + beta, aa = alpha/sum, bb = beta/sum;
			v = coords_q ( ver_south, tag::winding, winding_ver_south );
			coords_Eu ( shadow_south ) = v;
			v = coords_q ( ver_north, tag::winding, winding_ver_north );
			coords_Eu ( shadow_north ) = v;
			mani_Eu.interpolate ( C, bb*(1-frac_N), shadow_south, aa*frac_E,     shadow_east,     
		                           bb*frac_N,     shadow_north, aa*(1-frac_E), shadow_west );
			Cell BC ( tag::segment, B.reverse(), C );  // create a new segment
			Cell CD ( tag::segment, C.reverse(), D );  // create a new segment
			BC.winding() = -winding_B;
			CD.winding() =  winding_D;
		  assert ( AB.winding() + BC.winding() + CD.winding() + DA.winding() == 0 );
		  winding_D = 0;
			if ( cut_rectangles_in_half )
			{	Cell BD ( tag::segment, B.reverse(), D );  // create a new segment
				BD.winding() = winding_D - winding_B;
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
		CD.winding() = - DA.winding() - AB.winding() - BC.winding();
		assert ( winding_D == 0 );
		winding_B = winding_ver_west;
		if ( cut_rectangles_in_half )
		{	Cell BD ( tag::segment, B.reverse(), D );  // create a new segment
			BD.winding() = - AB.winding() - DA.winding();
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
		BC.winding() = - CD.winding() - DA.winding() - AB.winding();
		if ( cut_rectangles_in_half )
		{	Cell BD ( tag::segment, B.reverse(), D );  // create a new segment
			BD.winding() = - DA.winding() - AB.winding();
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
	assert ( AB.winding() + BC.winding() + CD.winding() + DA.winding() == 0 );
	if ( cut_rectangles_in_half )
	{	Cell BD ( tag::segment, B.reverse(), D );  // create a new segment
		BD.winding() = BC.winding() + CD.winding();
		Cell BCD ( tag::triangle, BD.reverse(), BC, CD );
		Cell ABD ( tag::triangle, BD, DA, AB );
		BCD.add_to_mesh (*this);  // 'this' is the mesh we are building
		ABD.add_to_mesh (*this);                          }
	else // with quadrilaterals
	{	Cell Q ( tag::rectangle, AB, BC, CD, DA );
		Q.add_to_mesh (*this);                     }
	it++;  assert ( it == horizon.end() );

} // end of Mesh::build with tag::quadrangle and tag::winding

//----------------------------------------------------------------------------------//


namespace { // anonymous namespace, mimics static linkage

Mesh fold_common ( const Mesh & msh, const std::map < Cell, Cell > & corresp_seg,
                     size_t dim, bool keep_map, std::map < Cell, Cell > & m      )

// new segments have already been built, kept in 'corresp_seg'
// build new cells (polygons)
	
{	assert ( dim == 2 );
	
	Mesh result ( tag::fuzzy, tag::of_dim, 2 );

	Mesh::Iterator it_cll = msh .iterator ( tag::over_cells_of_dim, 2 );
	for ( it_cll .reset(); it_cll .in_range(); it_cll++ )
	{	Cell cll = * it_cll;
		Cell::Positive::HighDim * new_cll_ptr = new Cell::Positive::HighDim
			( tag::whose_boundary_is,
				Mesh ( tag::whose_core_is,
			         new Mesh::Connected::OneDim ( tag::with,
		           cll .boundary() .number_of ( tag::segments ),
	               tag::segments, tag::one_dummy_wrapper       ),
	         tag::freshly_created                                 ),
				tag::one_dummy_wrapper                                     );
		Cell new_cll ( tag::whose_core_is, new_cll_ptr, tag::freshly_created );
		Cell kept_seg ( tag::non_existent );
		Mesh::Iterator it_bdry = cll .boundary() .iterator ( tag::over_segments );
		for ( it_bdry .reset(); it_bdry .in_range(); it_bdry ++ )
		{	Cell seg = * it_bdry;
			if ( seg .is_positive() )
			{	std::map < Cell, Cell > ::const_iterator it = corresp_seg .find ( seg );
				assert ( it != corresp_seg .end() );
				kept_seg = it->second;                                                   }
				// it ->second == corresp_seg [ seg ] );
			else
			{	std::map < Cell, Cell > ::const_iterator it =
					corresp_seg.find ( seg .reverse() );
				assert ( it != corresp_seg .end() );
				kept_seg = it->second .reverse();              }
				// it ->second == corresp_seg [ seg ] );
			assert ( kept_seg .exists() );
			kept_seg .core ->add_to_mesh
				( new_cll .boundary() .core, tag::do_not_bother );                            }
		assert ( kept_seg .exists() );
		new_cll .boundary() .closed_loop ( kept_seg .tip() );
		if ( keep_map )
		{	assert ( m .find ( new_cll ) == m .end() );
			m .insert ( std::pair < Cell, Cell > ( new_cll, cll ) );  }
			// m [ new_cll ] = cll;
		new_cll .add_to_mesh ( result );                                                     }

	return result;

}  // end of  fold_common

//----------------------------------------------------------------------------------//

	
Mesh fold_common_no_sides
( const Mesh & msh,
	const std::map < Cell, std::pair < Cell, Manifold::Action > > & corresp_ver,
  size_t dim, bool keep_map, std::map < Cell, Cell > & m                      )

// build 'corresp_seg' then call fold_common
	
{	if ( msh.dim() == 1 )
	{	Mesh result ( tag::fuzzy, tag::of_dim, 1 );
		Mesh::Iterator it_seg = msh.iterator ( tag::over_segments );
		for ( it_seg.reset(); it_seg.in_range(); it_seg++ )
		{	Cell seg = *it_seg;
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::const_iterator it_base_rev = corresp_ver.find ( seg.base().reverse() );
			assert ( it_base_rev != corresp_ver.end() );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::const_iterator it_tip = corresp_ver.find ( seg.tip() );
			assert ( it_tip != corresp_ver.end() );
			Cell new_seg ( tag::segment,
			               it_base_rev->second.first.reverse(), it_tip->second.first );
			new_seg.winding() = it_tip->second.second - it_base_rev->second.second;
			// new_seg.winding() = corresp_ver [ seg.tip()            ] .second -
			//	                corresp_ver [ seg.base().reverse() ] .second  ;
			new_seg.add_to_mesh ( result );                                            }
			
		return result;                                                                 }

	else
	{	assert ( msh.dim() == 2 );
		// we use a map -- for a faster code, we could use Cell::Core::hook
		std::map < Cell, Cell > corresp_seg;

		Mesh::Iterator it_seg = msh.iterator ( tag::over_segments );
		for ( it_seg.reset(); it_seg.in_range(); it_seg++ )
		{	Cell seg = *it_seg;
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::const_iterator it_base_rev = corresp_ver.find ( seg.base().reverse() );
			assert ( it_base_rev != corresp_ver.end() );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::const_iterator it_tip = corresp_ver.find ( seg.tip() );
			assert ( it_tip != corresp_ver.end() );
			Cell new_seg ( tag::segment,
			               it_base_rev->second.first.reverse(), it_tip->second.first );
			new_seg.winding() = it_tip->second.second - it_base_rev->second.second;
			// new_seg.winding() = corresp_ver [ seg.tip()            ] .second -
			//	                corresp_ver [ seg.base().reverse() ] .second  ;
			// inspired in item 24 of the book : Scott Meyers, Effective STL
			std::map < Cell, Cell > ::iterator it_map = corresp_seg.lower_bound ( seg );
			assert ( ( it_map == corresp_seg.end() ) or
			         ( corresp_seg.key_comp()(seg,it_map->first) ) );
			corresp_seg.emplace_hint ( it_map, std::piecewise_construct,
			  std::forward_as_tuple ( seg ), std::forward_as_tuple ( new_seg ) );     }
		// corresp_seg [ seg ] = new_seg;               

		return fold_common ( msh, corresp_seg, dim, keep_map, m );                       }

}  // end of  fold_common_no_sides

//----------------------------------------------------------------------------------//


Mesh fold_no_sides ( Mesh * that, const tag::BuildNewVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
									size_t dim, bool keep_map, std::map < Cell, Cell > & m )

// take a mesh and fold it around the current working manifold,
// which must be a quotient manifold

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

	// we use a map -- for a faster code, we could use Cell::Core::hook
	std::map < Cell, std::pair < Cell, Manifold::Action > > corresp_ver;

	Mesh::Iterator it_ver = that->iterator ( tag::over_vertices );
	for ( it_ver.reset(); it_ver.in_range(); it_ver++ )
	{	Cell V = *it_ver;
		Cell new_V ( tag::vertex );
		x ( new_V ) = x ( V );   y ( new_V ) = y ( V );
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_map = corresp_ver.lower_bound ( V );
		assert ( ( it_map == corresp_ver.end() ) or
		         ( corresp_ver.key_comp()(V,it_map->first) ) );
		corresp_ver.emplace_hint ( it_map, std::piecewise_construct,
		    std::forward_as_tuple ( V ), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { new_V, 0 } ) );  }
		// corresp_ver [ V ] = { new_V, 0 };

	return fold_common_no_sides ( *that, corresp_ver, dim, keep_map, m );     }
	
//----------------------------------------------------------------------------------//

	
Mesh fold_no_sides ( Mesh * that, const tag::UseExistingVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
									size_t dim, bool keep_map, std::map < Cell, Cell > & m )

// take a mesh and fold it around the current working manifold,
// which must be a quotient manifold

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

	// we use a map -- for a faster code, we could use Cell::Core::hook
	std::map < Cell, std::pair < Cell, Manifold::Action > > corresp_ver;

	Mesh::Iterator it_ver = that->iterator ( tag::over_vertices );
	for ( it_ver.reset(); it_ver.in_range(); it_ver++ )
	{	Cell V = *it_ver;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_map = corresp_ver.lower_bound ( V );
		assert ( ( it_map == corresp_ver.end() ) or
		         ( corresp_ver.key_comp()(V,it_map->first) ) );
		corresp_ver.emplace_hint ( it_map, std::piecewise_construct,
		    std::forward_as_tuple ( V ), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { V, 0 } ) );  }
		// corresp_ver [ V ] = { V, 0 };

	return fold_common_no_sides ( *that, corresp_ver, dim, keep_map, m );               }
	
//----------------------------------------------------------------------------------//

	
Mesh fold_common_two_sides
( const Mesh & msh,
  const std::map < Cell, std::pair < Cell, Manifold::Action > > & corresp_ver,
  const Mesh & side_1, const Mesh & side_2,
  size_t dim, bool keep_map, std::map < Cell, Cell > & m )

// build 'corresp_seg', taking care to identify side_1 with side_2
// then call fold_common
	
{	assert ( msh.dim() == 2 );
	// we use a map -- for a faster code, we could use Cell::Core::hook
	std::map < Cell, Cell > corresp_seg;

	Mesh::Iterator it_seg = msh.iterator ( tag::over_segments );
	for ( it_seg.reset(); it_seg.in_range(); it_seg++ )
	{	Cell seg = *it_seg;
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::const_iterator it_base_rev = corresp_ver.find ( seg.base().reverse() );
		assert ( it_base_rev != corresp_ver.end() );
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::const_iterator it_tip = corresp_ver.find ( seg.tip() );
		assert ( it_tip != corresp_ver.end() );
		if ( seg .belongs_to ( side_2, tag::not_oriented ) ) continue;
			// use corresponding segment on side_1
		Cell new_seg ( tag::segment,
		               it_base_rev->second.first.reverse(), it_tip->second.first );
		new_seg.winding() = it_tip->second.second - it_base_rev->second.second;
		// new_seg.winding() = corresp_ver [ seg.tip()            ] .second -
		//	                corresp_ver [ seg.base().reverse() ] .second  ;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, Cell > ::iterator it_map = corresp_seg.lower_bound ( seg );
		assert ( ( it_map == corresp_seg.end() ) or
		         ( corresp_seg.key_comp()(seg,it_map->first) ) );
		corresp_seg.emplace_hint ( it_map, std::piecewise_construct,
		  std::forward_as_tuple ( seg ), std::forward_as_tuple ( new_seg ) );     }
		// corresp_seg [ seg ] = new_seg;               

	Mesh::Iterator it_seg_1 = side_1.iterator
		( tag::over_segments, tag::require_order, tag::force_positive );
	Mesh::Iterator it_seg_2 = side_2.iterator
		( tag::over_segments, tag::require_order, tag::force_positive );
	for ( it_seg_1.reset(), it_seg_2.reset(); it_seg_1.in_range(); it_seg_1++, it_seg_2++ )
	{	assert ( it_seg_2 .in_range() );
		Cell seg_1 = *it_seg_1;
		Cell seg_2 = *it_seg_2;
		std::map < Cell, Cell > ::iterator it = corresp_seg.find ( seg_1 );
		assert ( it != corresp_seg.end() );
		Cell new_seg_1 = it->second .reverse();
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, Cell > ::iterator it_map = corresp_seg.lower_bound ( seg_2 );
		assert ( ( it_map == corresp_seg.end() ) or
		         ( corresp_seg.key_comp()(seg_2,it_map->first) ) );
		corresp_seg.emplace_hint ( it_map, std::piecewise_construct,
		  std::forward_as_tuple ( seg_2 ), std::forward_as_tuple ( new_seg_1 ) );     }
		// corresp_seg [ seg_2 ] = corresp_seg [ seg_1 ];               
	assert ( not it_seg_2 .in_range() );
		
	return fold_common ( msh, corresp_seg, dim, keep_map, m );

}  // end of  fold_common_two_sides
	
//----------------------------------------------------------------------------------//


Mesh fold_two_sides ( Mesh * that, const tag::Identify &, const Mesh & side_1,
                      const tag::With &, const Mesh & side_2,
                  const tag::BuildNewVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
									size_t dim, bool keep_map, std::map < Cell, Cell > & m )

// take a mesh whose external boundary has two parallel segments
// and identify these two segments

// a quotient manifold will be built with one action generator
	
{	// we use the current (Euclidian) manifold
	Manifold space = Manifold::working;
	assert ( space.exists() );
	Function coord = space.coordinates();
	assert ( coord.nb_of_components() == 2 );
	// we split 'coord' into its components
	Function x = coord[0], y = coord[1];

	// first we need to identify a translation which moves side_1 into side_2

	Cell A = side_1.first_vertex().reverse();
	Cell B = side_1.last_vertex();
	Cell C = side_2.last_vertex();
	Cell D = side_2.first_vertex().reverse();

	double dx = x(D) - x(A), dy = y(D) - y(A);
	double norm = std::sqrt ( dx*dx + dy*dy );
	assert ( std::abs ( dx - ( x(C) - x(B) ) ) < 1.e-4 * norm );
	assert ( std::abs ( dy - ( y(C) - y(B) ) ) < 1.e-4 * norm );

	// the desired translation is ( dx, dy )
	Manifold::Action g ( tag::transforms, coord, tag::into, (x+dx) && ( y+dy) );
	Manifold manif_q = space.quotient ( g );

	// we use a map -- for a faster code, we could use Cell::Core::hook
	std::map < Cell, std::pair < Cell, Manifold::Action > > corresp_ver;

	Mesh::Iterator it_ver = that->iterator ( tag::over_vertices );
	for ( it_ver.reset(); it_ver.in_range(); it_ver++ )
	{	Cell V = *it_ver;
		if ( V.belongs_to ( side_2 ) ) continue;
		Cell new_V ( tag::vertex );
		x ( new_V ) = x ( V );   y ( new_V ) = y ( V );
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_map = corresp_ver.lower_bound ( V );
		assert ( ( it_map == corresp_ver.end() ) or
		         ( corresp_ver.key_comp()(V,it_map->first) ) );
		corresp_ver.emplace_hint ( it_map, std::piecewise_construct,
		    std::forward_as_tuple ( V ), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { new_V, 0 } ) );  }
		// corresp_ver [ V ] = { new_V, 0 };

	assert ( side_1.number_of ( tag::segments ) == side_2.number_of ( tag::segments ) );
	Mesh::Iterator it1 = side_1.iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it2 = side_2.iterator ( tag::over_vertices, tag::require_order );
	for ( it1.reset(), it2.reset(); it1.in_range(); it1++, it2++ )
	{	assert ( it2.in_range() );
		Cell V = *it1;  Cell W = *it2;
		assert ( std::abs ( dx - ( x(W) - x(V) ) ) < 1.e-4 * norm );
		assert ( std::abs ( dy - ( y(W) - y(V) ) ) < 1.e-4 * norm );
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_V = corresp_ver.find ( V );
		assert ( it_V != corresp_ver.end() );
		assert ( it_V->second.second == 0 );
		Cell new_V = it_V->second.first;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_W = corresp_ver.lower_bound ( W );
		assert ( ( it_W == corresp_ver.end() ) or
		         ( corresp_ver.key_comp()(W,it_W->first) ) );
		corresp_ver.emplace_hint ( it_W, std::piecewise_construct,
		    std::forward_as_tuple ( W ), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { new_V, g } ) );  }
		// corresp_ver [ W ] = { new_V, g };
	assert ( not it2.in_range() );
						 
	return fold_common_two_sides ( *that, corresp_ver, side_1, side_2, dim, keep_map, m );

}  // end of fold_two_sides with tag::build_new_vertices

	
//----------------------------------------------------------------------------------//

	
Mesh fold_two_sides ( Mesh * that, const tag::Identify &, const Mesh & side_1,
                      const tag::With &, const Mesh & side_2,
                  const tag::UseExistingVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
									size_t dim, bool keep_map, std::map < Cell, Cell > & m )

// take a mesh whose external boundary has two parallel segments
// and identify these two segments

// a quotient manifold will be built with one action generator
	
{	// we use the current (Euclidian) manifold
	Manifold space = Manifold::working;
	assert ( space.exists() );
	Function coord = space.coordinates();
	assert ( coord.nb_of_components() == 2 );
	// we split 'coord' into its components
	Function x = coord[0], y = coord[1];

	// first we need to identify a translation which moves side_1 into side_2

	Cell A = side_1.first_vertex().reverse();
	Cell B = side_1.last_vertex();
	Cell C = side_2.last_vertex();
	Cell D = side_2.first_vertex().reverse();

	double dx = x(D) - x(A), dy = y(D) - y(A);
	double norm = std::sqrt ( dx*dx + dy*dy );
	assert ( std::abs ( dx - ( x(C) - x(B) ) ) < 1.e-4 * norm );
	assert ( std::abs ( dy - ( y(C) - y(B) ) ) < 1.e-4 * norm );

	// the desired translation is ( dx, dy )
	Manifold::Action g ( tag::transforms, coord, tag::into, (x+dx) && ( y+dy) );
	Manifold manif_q = space.quotient ( g );

	// we use a map -- for a faster code, we could use Cell::Core::hook
	std::map < Cell, std::pair < Cell, Manifold::Action > > corresp_ver;

	Mesh::Iterator it_ver = that->iterator ( tag::over_vertices );
	for ( it_ver.reset(); it_ver.in_range(); it_ver++ )
	{	Cell V = *it_ver;
		if ( V.belongs_to ( side_2 ) ) continue;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_map = corresp_ver.lower_bound ( V );
		assert ( ( it_map == corresp_ver.end() ) or
		         ( corresp_ver.key_comp()(V,it_map->first) ) );
		corresp_ver.emplace_hint ( it_map, std::piecewise_construct,
		    std::forward_as_tuple ( V ), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { V, 0 } ) );  }
		// corresp_ver [ V ] = { V, 0 };

	assert ( side_1.number_of ( tag::segments ) == side_2.number_of ( tag::segments ) );
	Mesh::Iterator it1 = side_1.iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it2 = side_2.iterator ( tag::over_vertices, tag::require_order );
	for ( it1.reset(), it2.reset(); it1.in_range(); it1++, it2++ )
	{	assert ( it2.in_range() );
		Cell V = *it1;  Cell W = *it2;
		assert ( std::abs ( dx - ( x(W) - x(V) ) ) < 1.e-4 * norm );
		assert ( std::abs ( dy - ( y(W) - y(V) ) ) < 1.e-4 * norm );
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_W = corresp_ver.lower_bound ( W );
		assert ( ( it_W == corresp_ver.end() ) or
		         ( corresp_ver.key_comp()(W,it_W->first) ) );
		corresp_ver.emplace_hint ( it_W, std::piecewise_construct,
		    std::forward_as_tuple ( W ), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { V, g } ) );  }
		// corresp_ver [ W ] = { V, g };
	assert ( not it2.in_range() );
						 
	return fold_common_two_sides ( *that, corresp_ver, side_1, side_2, dim, keep_map, m );

}  // end of fold_two_sides with tag::use_existing_vertices

	
//----------------------------------------------------------------------------------//

	
Mesh fold_common_four_sides
( const Mesh & msh,
  const std::map < Cell, std::pair < Cell, Manifold::Action > > & corresp_ver,
  const Mesh & side_1, const Mesh & side_2, const Mesh & side_3, const Mesh & side_4,
  size_t dim, bool keep_map, std::map < Cell, Cell > & m                             )

// build 'corresp_seg', taking care to identify side_1 with side_2 and side_3 with side_4
// then call fold_common
	
{	assert ( msh.dim() == 2 );
	// we use a map -- for a faster code, we could use Cell::Core::hook
	std::map < Cell, Cell > corresp_seg;

	Mesh::Iterator it_seg = msh.iterator ( tag::over_segments );
	for ( it_seg.reset(); it_seg.in_range(); it_seg++ )
	{	Cell seg = *it_seg;
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::const_iterator it_base_rev = corresp_ver.find ( seg.base().reverse() );
		assert ( it_base_rev != corresp_ver.end() );
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::const_iterator it_tip = corresp_ver.find ( seg.tip() );
		assert ( it_tip != corresp_ver.end() );
		if ( seg .belongs_to ( side_2, tag::not_oriented ) ) continue;
			// use corresponding segment on side_1
		if ( seg .belongs_to ( side_4, tag::not_oriented ) ) continue;
			// use corresponding segment on side_3
		Cell new_seg ( tag::segment,
		               it_base_rev->second.first.reverse(), it_tip->second.first );
		new_seg.winding() = it_tip->second.second - it_base_rev->second.second;
		// new_seg.winding() = corresp_ver [ seg.tip()            ] .second -
		//	                   corresp_ver [ seg.base().reverse() ] .second  ;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, Cell > ::iterator it_map = corresp_seg.lower_bound ( seg );
		assert ( ( it_map == corresp_seg.end() ) or
		         ( corresp_seg.key_comp()(seg,it_map->first) ) );
		corresp_seg.emplace_hint ( it_map, std::piecewise_construct,
		  std::forward_as_tuple ( seg ), std::forward_as_tuple ( new_seg ) );     }
		// corresp_seg [ seg ] = new_seg;               

	Mesh::Iterator it_seg_1 = side_1.iterator
		( tag::over_segments, tag::require_order, tag::force_positive );
	Mesh::Iterator it_seg_2 = side_2.iterator
		( tag::over_segments, tag::require_order, tag::force_positive );
	for ( it_seg_1.reset(), it_seg_2.reset(); it_seg_1.in_range(); it_seg_1++, it_seg_2++ )
	{	assert ( it_seg_2 .in_range() );
		Cell seg_1 = *it_seg_1;
		Cell seg_2 = *it_seg_2;
		std::map < Cell, Cell > ::iterator it = corresp_seg.find ( seg_1 );
		assert ( it != corresp_seg.end() );
		Cell new_seg_1 = it->second .reverse();
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, Cell > ::iterator it_map = corresp_seg.lower_bound ( seg_2 );
		assert ( ( it_map == corresp_seg.end() ) or
		         ( corresp_seg.key_comp()(seg_2,it_map->first) ) );
		corresp_seg.emplace_hint ( it_map, std::piecewise_construct,
		  std::forward_as_tuple ( seg_2 ), std::forward_as_tuple ( new_seg_1 ) );     }
		// corresp_seg [ seg_2 ] = corresp_seg [ seg_1 ];               
	assert ( not it_seg_2 .in_range() );
		
	Mesh::Iterator it_seg_3 = side_3.iterator
		( tag::over_segments, tag::require_order, tag::force_positive );
	Mesh::Iterator it_seg_4 = side_4.iterator
		( tag::over_segments, tag::require_order, tag::force_positive );
	for ( it_seg_3.reset(), it_seg_4.reset(); it_seg_3.in_range(); it_seg_3++, it_seg_4++ )
	{	assert ( it_seg_4 .in_range() );
		Cell seg_3 = *it_seg_3;
		Cell seg_4 = *it_seg_4;
		std::map < Cell, Cell > ::iterator it = corresp_seg.find ( seg_3 );
		assert ( it != corresp_seg.end() );
		Cell new_seg_3 = it->second .reverse();
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, Cell > ::iterator it_map = corresp_seg.lower_bound ( seg_4 );
		assert ( ( it_map == corresp_seg.end() ) or
		         ( corresp_seg.key_comp()(seg_4,it_map->first) ) );
		corresp_seg.emplace_hint ( it_map, std::piecewise_construct,
		  std::forward_as_tuple ( seg_4 ), std::forward_as_tuple ( new_seg_3 ) );     }
		// corresp_seg [ seg_4 ] = corresp_seg [ seg_3 ];               
	assert ( not it_seg_4 .in_range() );
		
	return fold_common ( msh, corresp_seg, dim, keep_map, m );

}  // end of  fold_common_four_sides
	
//----------------------------------------------------------------------------------//


Mesh fold_four_sides ( Mesh * that, const tag::Identify &, const Mesh & side_1,
                      const tag::With &, const Mesh & side_2,
                  const tag::Identify &, const Mesh & side_3,
                      const tag::With &, const Mesh & side_4,
                  const tag::BuildNewVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
									size_t dim, bool keep_map, std::map < Cell, Cell > & m )

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

	// we need to identify two translations, one which moves side_1 into side_2
	// and another one which moves side_3 into side_4

	Cell A = side_1.first_vertex().reverse();
	Cell B = side_1.last_vertex();
	Cell C = side_2.last_vertex();
	Cell D = side_2.first_vertex().reverse();

	double dx12 = x(D) - x(A), dy12 = y(D) - y(A);
	double norm12 = std::sqrt ( dx12*dx12 + dy12*dy12 );
	assert ( std::abs ( dx12 - ( x(C) - x(B) ) ) < 1.e-4 * norm12 );
	assert ( std::abs ( dy12 - ( y(C) - y(B) ) ) < 1.e-4 * norm12 );

	A = side_3.first_vertex().reverse();
	B = side_3.last_vertex();
	C = side_4.last_vertex();
	D = side_4.first_vertex().reverse();

	double dx34 = x(D) - x(A), dy34 = y(D) - y(A);
	double norm34 = std::sqrt ( dx34*dx34 + dy34*dy34 );
	assert ( std::abs ( dx34 - ( x(C) - x(B) ) ) < 1.e-4 * norm34 );
	assert ( std::abs ( dy34 - ( y(C) - y(B) ) ) < 1.e-4 * norm34 );

	// the desired translations are ( dx12, dy12 ), ( dx34, dy34 )
	Manifold::Action g12 ( tag::transforms, coord, tag::into, (x+dx12) && ( y+dy12) );
	Manifold::Action g34 ( tag::transforms, coord, tag::into, (x+dx34) && ( y+dy34) );
	Manifold manif_q = space.quotient ( g12, g34 );

	// we use a map -- for a faster code, we could use Cell::Core::hook
	std::map < Cell, std::pair < Cell, Manifold::Action > > corresp_ver;

	Mesh::Iterator it_ver = that->iterator ( tag::over_vertices );
	for ( it_ver.reset(); it_ver.in_range(); it_ver++ )
	{	Cell V = *it_ver;
		if ( V.belongs_to ( side_2 ) ) continue;
		if ( V.belongs_to ( side_4 ) ) continue;
		Cell new_V ( tag::vertex );
		x ( new_V ) = x ( V );   y ( new_V ) = y ( V );
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_map = corresp_ver.lower_bound ( V );
		assert ( ( it_map == corresp_ver.end() ) or
		         ( corresp_ver.key_comp()(V,it_map->first) ) );
		corresp_ver.emplace_hint ( it_map, std::piecewise_construct,
		    std::forward_as_tuple ( V ), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { new_V, 0 } ) );  }
		// corresp_ver [ V ] = { new_V, 0 };

	assert ( side_1.number_of ( tag::segments ) == side_2.number_of ( tag::segments ) );
	Mesh::Iterator it1 = side_1.iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it2 = side_2.iterator ( tag::over_vertices, tag::require_order );
	for ( it1.reset(), it2.reset(); it1.in_range(); it1++, it2++ )
	{	assert ( it2.in_range() );
		Cell V = *it1;  Cell W = *it2;
		assert ( std::abs ( dx12 - ( x(W) - x(V) ) ) < 1.e-4 * norm12 );
		assert ( std::abs ( dy12 - ( y(W) - y(V) ) ) < 1.e-4 * norm12 );
		if ( V.belongs_to ( side_4 ) ) continue;
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_V = corresp_ver.find ( V );
		assert ( it_V != corresp_ver.end() );
		assert ( it_V->second.second == 0 );
		Cell VV = it_V->second.first;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_W = corresp_ver.lower_bound ( W );
		assert ( ( it_W == corresp_ver.end() ) or
		         ( corresp_ver.key_comp()(W,it_W->first) ) );
		corresp_ver.emplace_hint ( it_W, std::piecewise_construct,
		    std::forward_as_tuple ( W ), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { VV, g12 } ) );  }
		// corresp_ver [ W ] = { VV, g12 };
	assert ( not it2.in_range() );

	Cell origin ( tag::non_existent ), corner ( tag::non_existent );
	assert ( side_3.number_of ( tag::segments ) == side_4.number_of ( tag::segments ) );
	Mesh::Iterator it3 = side_3.iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it4 = side_4.iterator ( tag::over_vertices, tag::require_order );
	for ( it3.reset(), it4.reset(); it3.in_range(); it3++, it4++ )
	{	assert ( it4.in_range() );
		Cell V = *it3;  Cell W = *it4;
		assert ( std::abs ( dx34 - ( x(W) - x(V) ) ) < 1.e-4 * norm34 );
		assert ( std::abs ( dy34 - ( y(W) - y(V) ) ) < 1.e-4 * norm34 );
		if ( W .belongs_to ( side_2 ) )
		{	std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_W = corresp_ver.find ( W );
			assert ( it_W == corresp_ver.end() );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V = corresp_ver.find ( V );
			assert ( it_V != corresp_ver.end() );
			assert ( it_V->second .second == g12 );
			origin = it_V->second .first;
			corner = W;
			continue;                                                             }
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_V = corresp_ver.find ( V );
		assert ( it_V != corresp_ver.end() );
		assert ( it_V->second.second == 0 );
		Cell new_V = it_V->second.first;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_W = corresp_ver.lower_bound ( W );
		assert ( ( it_W == corresp_ver.end() ) or
		         ( corresp_ver.key_comp()(W,it_W->first) ) );
		corresp_ver.emplace_hint ( it_W, std::piecewise_construct,
		    std::forward_as_tuple ( W ), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { new_V, g34 } ) );  }
		// corresp_ver [ W ] = { new_V, g34 };
	assert ( not it4.in_range() );

	assert ( origin.exists() );  assert ( corner.exists() );
	// inspired in item 24 of the book : Scott Meyers, Effective STL
	std::map < Cell, std::pair < Cell, Manifold::Action > >
		::iterator it_c = corresp_ver.lower_bound ( corner );
	assert ( ( it_c == corresp_ver.end() ) or
	         ( corresp_ver.key_comp()(corner,it_c->first) ) );
	corresp_ver.emplace_hint ( it_c, std::piecewise_construct,
	    std::forward_as_tuple ( corner ), std::forward_as_tuple
	    ( std::pair < Cell, Manifold::Action > { origin, g12 + g34 } ) );
	// corresp_ver [ corner ] = { origin, g12 + g34 };
	
	return fold_common_four_sides
		( *that, corresp_ver, side_1, side_2, side_3, side_4, dim, keep_map, m );

}  // end of fold_four_sides with tag::build_new_vertices

//----------------------------------------------------------------------------------//

	
Mesh fold_four_sides ( Mesh * that, const tag::Identify &, const Mesh & side_1,
                      const tag::With &, const Mesh & side_2,
                  const tag::Identify &, const Mesh & side_3,
                      const tag::With &, const Mesh & side_4,
                  const tag::UseExistingVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
									size_t dim, bool keep_map, std::map < Cell, Cell > & m )

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

	// we need to identify two translations, one which moves side_1 into side_2
	// and another one which moves side_3 into side_4

	Cell A = side_1.first_vertex().reverse();
	Cell B = side_1.last_vertex();
	Cell C = side_2.last_vertex();
	Cell D = side_2.first_vertex().reverse();

	double dx12 = x(D) - x(A), dy12 = y(D) - y(A);
	double norm12 = std::sqrt ( dx12*dx12 + dy12*dy12 );
	assert ( std::abs ( dx12 - ( x(C) - x(B) ) ) < 1.e-4 * norm12 );
	assert ( std::abs ( dy12 - ( y(C) - y(B) ) ) < 1.e-4 * norm12 );

	A = side_3.first_vertex().reverse();
	B = side_3.last_vertex();
	C = side_4.last_vertex();
	D = side_4.first_vertex().reverse();

	double dx34 = x(D) - x(A), dy34 = y(D) - y(A);
	double norm34 = std::sqrt ( dx34*dx34 + dy34*dy34 );
	assert ( std::abs ( dx34 - ( x(C) - x(B) ) ) < 1.e-4 * norm34 );
	assert ( std::abs ( dy34 - ( y(C) - y(B) ) ) < 1.e-4 * norm34 );

	// the desired translations are ( dx12, dy12 ), ( dx34, dy34 )
	Manifold::Action g12 ( tag::transforms, coord, tag::into, (x+dx12) && ( y+dy12) );
	Manifold::Action g34 ( tag::transforms, coord, tag::into, (x+dx34) && ( y+dy34) );
	Manifold manif_q = space.quotient ( g12, g34 );

	// we use a map -- for a faster code, we could use Cell::Core::hook
	std::map < Cell, std::pair < Cell, Manifold::Action > > corresp_ver;

	Mesh::Iterator it_ver = that->iterator ( tag::over_vertices );
	for ( it_ver.reset(); it_ver.in_range(); it_ver++ )
	{	Cell V = *it_ver;
		if ( V.belongs_to ( side_2 ) ) continue;
		if ( V.belongs_to ( side_4 ) ) continue;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_map = corresp_ver.lower_bound ( V );
		assert ( ( it_map == corresp_ver.end() ) or
		         ( corresp_ver.key_comp()(V,it_map->first) ) );
		corresp_ver.emplace_hint ( it_map, std::piecewise_construct,
		    std::forward_as_tuple ( V ), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { V, 0 } ) );  }
		// corresp_ver [ V ] = { V, 0 };

	assert ( side_1.number_of ( tag::segments ) == side_2.number_of ( tag::segments ) );
	Mesh::Iterator it1 = side_1.iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it2 = side_2.iterator ( tag::over_vertices, tag::require_order );
	for ( it1.reset(), it2.reset(); it1.in_range(); it1++, it2++ )
	{	assert ( it2.in_range() );
		Cell V = *it1;  Cell W = *it2;
		assert ( std::abs ( dx12 - ( x(W) - x(V) ) ) < 1.e-4 * norm12 );
		assert ( std::abs ( dy12 - ( y(W) - y(V) ) ) < 1.e-4 * norm12 );
		if ( V.belongs_to ( side_4 ) ) continue;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_W = corresp_ver.lower_bound ( W );
		assert ( ( it_W == corresp_ver.end() ) or
		         ( corresp_ver.key_comp()(W,it_W->first) ) );
		corresp_ver.emplace_hint ( it_W, std::piecewise_construct,
		    std::forward_as_tuple ( W ), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { V, g12 } ) );  }
		// corresp_ver [ W ] = { V, g12 };
	assert ( not it2.in_range() );

	Cell origin ( tag::non_existent ), corner ( tag::non_existent );
	assert ( side_3.number_of ( tag::segments ) == side_4.number_of ( tag::segments ) );
	Mesh::Iterator it3 = side_3.iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it4 = side_4.iterator ( tag::over_vertices, tag::require_order );
	for ( it3.reset(), it4.reset(); it3.in_range(); it3++, it4++ )
	{	assert ( it4.in_range() );
		Cell V = *it3;  Cell W = *it4;
		assert ( std::abs ( dx34 - ( x(W) - x(V) ) ) < 1.e-4 * norm34 );
		assert ( std::abs ( dy34 - ( y(W) - y(V) ) ) < 1.e-4 * norm34 );
		if ( W .belongs_to ( side_2 ) )
		{	std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_W = corresp_ver.find ( W );
			assert ( it_W == corresp_ver.end() );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V = corresp_ver.find ( V );
			assert ( it_V != corresp_ver.end() );
			assert ( it_V->second .second == g12 );
			origin = it_V->second .first;
			corner = W;
			continue;                                                             }
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_W = corresp_ver.lower_bound ( W );
		assert ( ( it_W == corresp_ver.end() ) or
		         ( corresp_ver.key_comp()(W,it_W->first) ) );
		corresp_ver.emplace_hint ( it_W, std::piecewise_construct,
		    std::forward_as_tuple ( W ), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { V, g34 } ) );  }
		// corresp_ver [ W ] = { V, g34 };
	assert ( not it4.in_range() );

	assert ( origin.exists() );  assert ( corner.exists() );
	// inspired in item 24 of the book : Scott Meyers, Effective STL
	std::map < Cell, std::pair < Cell, Manifold::Action > >
		::iterator it_c = corresp_ver.lower_bound ( corner );
	assert ( ( it_c == corresp_ver.end() ) or
	         ( corresp_ver.key_comp()(corner,it_c->first) ) );
	corresp_ver.emplace_hint ( it_c, std::piecewise_construct,
	    std::forward_as_tuple ( corner ), std::forward_as_tuple
	    ( std::pair < Cell, Manifold::Action > { origin, g12 + g34 } ) );
	// corresp_ver [ corner ] = { origin, g12 + g34 };
	
	return fold_common_four_sides
		( *that, corresp_ver, side_1, side_2, side_3, side_4, dim, keep_map, m );

}  // end of fold_four_sides with tag::use_existing_vertices
	
//----------------------------------------------------------------------------------//

	
Mesh fold_common_six_sides
( const Mesh & msh,
  const std::map < Cell, std::pair < Cell, Manifold::Action > > & corresp_ver,
  const Mesh & side_1, const Mesh & side_2, const Mesh & side_3, const Mesh & side_4,
  const Mesh & side_5, const Mesh & side_6,
  size_t dim, bool keep_map, std::map < Cell, Cell > & m                             )

// build 'corresp_seg', taking care to identify
//   side_1 with side_2, side_3 with side_4 and side_5 with side_6
// then call fold_common
	
{	assert ( msh.dim() == 2 );
	// we use a map -- for a faster code, we could use Cell::Core::hook
	std::map < Cell, Cell > corresp_seg;

	Mesh::Iterator it_seg = msh.iterator ( tag::over_segments );
	for ( it_seg.reset(); it_seg.in_range(); it_seg++ )
	{	Cell seg = *it_seg;
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::const_iterator it_base_rev = corresp_ver.find ( seg.base().reverse() );
		assert ( it_base_rev != corresp_ver.end() );
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::const_iterator it_tip = corresp_ver.find ( seg.tip() );
		assert ( it_tip != corresp_ver.end() );
		if ( seg .belongs_to ( side_2, tag::not_oriented ) ) continue;
			// use corresponding segment on side_1
		if ( seg .belongs_to ( side_4, tag::not_oriented ) ) continue;
			// use corresponding segment on side_3
		if ( seg .belongs_to ( side_6, tag::not_oriented ) ) continue;
			// use corresponding segment on side_5
		Cell new_seg ( tag::segment,
		               it_base_rev->second.first.reverse(), it_tip->second.first );
		new_seg.winding() = it_tip->second.second - it_base_rev->second.second;
		// new_seg.winding() = corresp_ver [ seg.tip()            ] .second -
		//	                corresp_ver [ seg.base().reverse() ] .second  ;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, Cell > ::iterator it_map = corresp_seg.lower_bound ( seg );
		assert ( ( it_map == corresp_seg.end() ) or
		         ( corresp_seg.key_comp()(seg,it_map->first) ) );
		corresp_seg.emplace_hint ( it_map, std::piecewise_construct,
		  std::forward_as_tuple ( seg ), std::forward_as_tuple ( new_seg ) );     }
		// corresp_seg [ seg ] = new_seg;               

	Mesh::Iterator it_seg_1 = side_1.iterator
		( tag::over_segments, tag::require_order, tag::force_positive );
	Mesh::Iterator it_seg_2 = side_2.iterator
		( tag::over_segments, tag::require_order, tag::force_positive );
	for ( it_seg_1.reset(), it_seg_2.reset(); it_seg_1.in_range(); it_seg_1++, it_seg_2++ )
	{	assert ( it_seg_2 .in_range() );
		Cell seg_1 = *it_seg_1;
		Cell seg_2 = *it_seg_2;
		std::map < Cell, Cell > ::iterator it = corresp_seg.find ( seg_1 );
		assert ( it != corresp_seg.end() );
		Cell new_seg_1 = it->second .reverse();
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, Cell > ::iterator it_map = corresp_seg.lower_bound ( seg_2 );
		assert ( ( it_map == corresp_seg.end() ) or
		         ( corresp_seg.key_comp()(seg_2,it_map->first) ) );
		corresp_seg.emplace_hint ( it_map, std::piecewise_construct,
		  std::forward_as_tuple ( seg_2 ), std::forward_as_tuple ( new_seg_1 ) );     }
		// corresp_seg [ seg_2 ] = corresp_seg [ seg_1 ];               
	assert ( not it_seg_2 .in_range() );
		
	Mesh::Iterator it_seg_3 = side_3.iterator
		( tag::over_segments, tag::require_order, tag::force_positive );
	Mesh::Iterator it_seg_4 = side_4.iterator
		( tag::over_segments, tag::require_order, tag::force_positive );
	for ( it_seg_3.reset(), it_seg_4.reset(); it_seg_3.in_range(); it_seg_3++, it_seg_4++ )
	{	assert ( it_seg_4 .in_range() );
		Cell seg_3 = *it_seg_3;
		Cell seg_4 = *it_seg_4;
		std::map < Cell, Cell > ::iterator it = corresp_seg.find ( seg_3 );
		assert ( it != corresp_seg.end() );
		Cell new_seg_3 = it->second .reverse();
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, Cell > ::iterator it_map = corresp_seg.lower_bound ( seg_4 );
		assert ( ( it_map == corresp_seg.end() ) or
		         ( corresp_seg.key_comp()(seg_4,it_map->first) ) );
		corresp_seg.emplace_hint ( it_map, std::piecewise_construct,
		  std::forward_as_tuple ( seg_4 ), std::forward_as_tuple ( new_seg_3 ) );     }
		// corresp_seg [ seg_4 ] = corresp_seg [ seg_3 ];               
	assert ( not it_seg_4 .in_range() );
		
	Mesh::Iterator it_seg_5 = side_5.iterator
		( tag::over_segments, tag::require_order, tag::force_positive );
	Mesh::Iterator it_seg_6 = side_6.iterator
		( tag::over_segments, tag::require_order, tag::force_positive );
	for ( it_seg_5.reset(), it_seg_6.reset(); it_seg_5.in_range(); it_seg_5++, it_seg_6++ )
	{	assert ( it_seg_6 .in_range() );
		Cell seg_5 = *it_seg_5;
		Cell seg_6 = *it_seg_6;
		std::map < Cell, Cell > ::iterator it = corresp_seg.find ( seg_5 );
		assert ( it != corresp_seg.end() );
		Cell new_seg_5 = it->second .reverse();
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, Cell > ::iterator it_map = corresp_seg.lower_bound ( seg_6 );
		assert ( ( it_map == corresp_seg.end() ) or
		         ( corresp_seg.key_comp()(seg_6,it_map->first) ) );
		corresp_seg.emplace_hint ( it_map, std::piecewise_construct,
		  std::forward_as_tuple ( seg_6 ), std::forward_as_tuple ( new_seg_5 ) );     }
		// corresp_seg [ seg_6 ] = corresp_seg [ seg_5 ];               
	assert ( not it_seg_6 .in_range() );
		
	return fold_common ( msh, corresp_seg, dim, keep_map, m );

}  // end of  fold_common_six_sides
	
//----------------------------------------------------------------------------------//

	
Mesh fold_six_sides ( Mesh * that, const tag::Identify &, const Mesh & side_1,
                      const tag::With &, const Mesh & side_2,
                  const tag::Identify &, const Mesh & side3,
                      const tag::With &, const Mesh & side4,
                  const tag::Identify &, const Mesh & side5,
                      const tag::With &, const Mesh & side6,
                  const tag::BuildNewVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
									size_t dim, bool keep_map, std::map < Cell, Cell > & m )

// take a mesh having as external boundary a hexagon
// and identify three pairs of opposite sides

// a quotient manifold with two action generators will be built
	
{	// we use the current manifold
	Manifold space = Manifold::working;
	assert ( space.exists() );
	Function coord = space.coordinates();
	assert ( coord.nb_of_components() == 2 );
	// we split 'coord' into its components
	Function x = coord[0],  y = coord[1];

	// we need to identify two translations, one which moves side_1 into side_2
	// and another one which moves side3 into side4

	Cell A = side_1.first_vertex().reverse();
	Cell B = side_1.last_vertex();
	Cell C = side_2.last_vertex();
	Cell D = side_2.first_vertex().reverse();

	double dx12 = x(D) - x(A), dy12 = y(D) - y(A);
	double norm12 = std::sqrt ( dx12*dx12 + dy12*dy12 );
	assert ( std::abs ( dx12 - ( x(C) - x(B) ) ) < 1.e-4 * norm12 );
	assert ( std::abs ( dy12 - ( y(C) - y(B) ) ) < 1.e-4 * norm12 );

	Mesh side_3 ( tag::non_existent ), side_4 ( tag::non_existent );
	// we want side_1 and side_3 to touch
	if ( A.belongs_to ( side3 ) or B.belongs_to ( side3 ) )
	{	side_3 = side3; side_4 = side4;  }
	else
	{	assert ( A.belongs_to ( side4 ) or B.belongs_to ( side4 ) );
		side_3 = side4; side_4 = side3;                              }

	A = side_3.first_vertex().reverse();
	B = side_3.last_vertex();
	C = side_4.last_vertex();
	D = side_4.first_vertex().reverse();

	double dx34 = x(D) - x(A), dy34 = y(D) - y(A);
	double norm34 = std::sqrt ( dx34*dx34 + dy34*dy34 );
	assert ( std::abs ( dx34 - ( x(C) - x(B) ) ) < 1.e-4 * norm34 );
	assert ( std::abs ( dy34 - ( y(C) - y(B) ) ) < 1.e-4 * norm34 );

	Mesh side_5 ( tag::non_existent ), side_6 ( tag::non_existent );
	// we want side_3 and side_5 to touch
	if ( A.belongs_to ( side5 ) or B.belongs_to ( side5 ) )
	{	side_5 = side5; side_6 = side6;  }
	else
	{	assert ( A.belongs_to ( side6 ) or B.belongs_to ( side6 ) );
		side_5 = side6; side_6 = side5;                              }

	A = side_5.first_vertex().reverse();
	B = side_5.last_vertex();
	C = side_6.last_vertex();
	D = side_6.first_vertex().reverse();

	double dx56 = x(D) - x(A), dy56 = y(D) - y(A);
	double norm56 = std::sqrt ( dx56*dx56 + dy56*dy56 );
	assert ( std::abs ( dx56 - ( x(C) - x(B) ) ) < 1.e-4 * norm56 );
	assert ( std::abs ( dy56 - ( y(C) - y(B) ) ) < 1.e-4 * norm56 );

	// we are confident that side_3 lies between side_1 and side_5
	// orientations may vary, we do not care
	// likewise, side_4 lies between side_2 and side_6

	// d56 is either  d12 + d34  or  d12 - d34  or  -d12 - d34  or  -d12 + d34
	std::vector < std::vector < double > > signs { { 1, 1 }, { 1, -1 }, { -1, -1 }, { -1, 1 } };
	size_t index_min_dif = 0;
	double min_dif = ( std::abs ( dx56 - signs[0][0]*dx12 - signs[0][1]*dx34 ) +
	                   std::abs ( dy56 - signs[0][0]*dy12 - signs[0][1]*dy34 )   );
	for ( size_t i = 1; i < 4; i++ )
	{	double dif = ( std::abs ( dx56 - signs[i][0]*dx12 - signs[i][1]*dx34 ) +
		               std::abs ( dy56 - signs[i][0]*dy12 - signs[i][1]*dy34 )   );
		if ( dif < min_dif )  {  min_dif = dif;  index_min_dif = i;  }              }
	assert ( min_dif < 1.e-4 * norm56 );

	// the desired translations are ( dx12, dy12 ), ( dx34, dy34 )
	Manifold::Action g12 ( tag::transforms, coord, tag::into, (x+dx12) && ( y+dy12) );
	Manifold::Action g34 ( tag::transforms, coord, tag::into, (x+dx34) && ( y+dy34) );
	Manifold manif_q = space.quotient ( g12, g34 );

	// a third translation, not used in the definition of the quotient manifold
	Manifold::Action g56 = signs[index_min_dif][0]*g12 +
	                                     signs[index_min_dif][1]*g34  ;

	// we use a map -- for a faster code, we could use Cell::Core::hook
	std::map < Cell, std::pair < Cell, Manifold::Action > > corresp_ver;

	Mesh::Iterator it_ver = that->iterator ( tag::over_vertices );
	for ( it_ver.reset(); it_ver.in_range(); it_ver++ )
	{	Cell V = *it_ver;
		if ( V.belongs_to ( side_2 ) ) continue;
		if ( V.belongs_to ( side_4 ) ) continue;
		if ( V.belongs_to ( side_6 ) ) continue;
		Cell new_V ( tag::vertex );
		x ( new_V ) = x ( V );   y ( new_V ) = y ( V );
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_map = corresp_ver.lower_bound ( V );
		assert ( ( it_map == corresp_ver.end() ) or
		         ( corresp_ver.key_comp()(V,it_map->first) ) );
		corresp_ver.emplace_hint ( it_map, std::piecewise_construct,
		    std::forward_as_tuple ( V ), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { new_V, 0 } ) );  }
		// corresp_ver [ V ] = { new_V, 0 };

	Cell V16 ( tag::non_existent ), V24 ( tag::non_existent );
	Cell V13 ( tag::non_existent ), V25 ( tag::non_existent );
	assert ( side_1.number_of ( tag::segments ) == side_2.number_of ( tag::segments ) );
	Mesh::Iterator it1 = side_1.iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it2 = side_2.iterator ( tag::over_vertices, tag::require_order );
	for ( it1.reset(), it2.reset(); it1.in_range(); it1++, it2++ )
	{	assert ( it2.in_range() );
		Cell V = *it1;  Cell W = *it2;
		assert ( std::abs ( dx12 - ( x(W) - x(V) ) ) < 1.e-4 * norm12 );
		assert ( std::abs ( dy12 - ( y(W) - y(V) ) ) < 1.e-4 * norm12 );
		if ( V.belongs_to ( side_3 ) )
		{	assert ( W.belongs_to ( side_5 ) );
			V13 = V;  V25 = W;               }
		if ( V.belongs_to ( side_6 ) )
		{	assert ( W.belongs_to ( side_4 ) );
			V16 = V;  V24 = W;  continue;    }
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_V = corresp_ver.find ( V );
		assert ( it_V != corresp_ver.end() );
		assert ( it_V->second.second == 0 );
		Cell VV = it_V->second.first;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_W = corresp_ver.lower_bound ( W );
		assert ( ( it_W == corresp_ver.end() ) or
		         ( corresp_ver.key_comp()(W,it_W->first) ) );
		corresp_ver.emplace_hint ( it_W, std::piecewise_construct,
		    std::forward_as_tuple ( W ), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { VV, g12 } ) );  }
		// corresp_ver [ W ] = { VV, g12 };
	assert ( not it2.in_range() );
	assert ( V16.exists() );  assert ( V24.exists() );
	assert ( V13.exists() );  assert ( V25.exists() );

	Cell V35 ( tag::non_existent ), V46 ( tag::non_existent );
	assert ( side_3.number_of ( tag::segments ) == side_4.number_of ( tag::segments ) );
	Mesh::Iterator it3 = side_3.iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it4 = side_4.iterator ( tag::over_vertices, tag::require_order );
	for ( it3.reset(), it4.reset(); it3.in_range(); it3++, it4++ )
	{	assert ( it4.in_range() );
		Cell V = *it3;  Cell W = *it4;
		assert ( std::abs ( dx34 - ( x(W) - x(V) ) ) < 1.e-4 * norm34 );
		assert ( std::abs ( dy34 - ( y(W) - y(V) ) ) < 1.e-4 * norm34 );
		if ( V .belongs_to ( side_5 ) )
		{	assert ( W.belongs_to ( side_2 ) );
			V35 = V;
			assert ( W == V24 );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_W = corresp_ver.find ( W );
			assert ( it_W == corresp_ver.end() );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V = corresp_ver.find ( V );
			assert ( it_V != corresp_ver.end() );
			assert ( it_V->second .second == 0 );                                }
		if ( V .belongs_to ( side_1 ) )
		{	assert ( W.belongs_to ( side_6 ) );
			V46 = W;
			assert ( V == V13 );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_W = corresp_ver.find ( W );
			assert ( it_W == corresp_ver.end() );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V = corresp_ver.find ( V );
			assert ( it_V != corresp_ver.end() );
			assert ( it_V->second .second == 0 );                                }
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_V = corresp_ver.find ( V );
		assert ( it_V != corresp_ver.end() );
		assert ( it_V->second.second == 0 );
		Cell new_V = it_V->second.first;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_W = corresp_ver.lower_bound ( W );
		assert ( ( it_W == corresp_ver.end() ) or
		         ( corresp_ver.key_comp()(W,it_W->first) ) );
		corresp_ver.emplace_hint ( it_W, std::piecewise_construct,
		    std::forward_as_tuple ( W ), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { new_V, g34 } ) );  }
		// corresp_ver [ W ] = { new_V, g34 };
	assert ( not it4.in_range() );
	assert ( V35.exists() );  assert ( V46.exists() );

	assert ( side_5.number_of ( tag::segments ) == side_6.number_of ( tag::segments ) );
	Mesh::Iterator it5 = side_5.iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it6 = side_6.iterator ( tag::over_vertices, tag::require_order );
	for ( it5.reset(), it6.reset(); it5.in_range(); it5++, it6++ )
	{	assert ( it6.in_range() );
		Cell V = *it5;  Cell W = *it6;
		assert ( std::abs ( dx56 - ( x(W) - x(V) ) ) < 1.e-4 * norm56 );
		assert ( std::abs ( dy56 - ( y(W) - y(V) ) ) < 1.e-4 * norm56 );
		if ( V .belongs_to ( side_2 ) )
		{	assert ( W.belongs_to ( side_4 ) );
			assert ( V == V25 );
			assert ( W == V46 );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_W = corresp_ver.find ( W );
			assert ( it_W != corresp_ver.end() );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V13 = corresp_ver.find ( V13 );
			assert ( it_V13 != corresp_ver.end() );
			assert ( it_V13->second .second == 0 );
			Cell new_V13 = it_V13->second.first;
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V = corresp_ver.find ( V );
			assert ( it_V != corresp_ver.end() );
			assert ( it_V->second .first == new_V13 );
			assert ( it_V->second .second == g12 );
			continue;                                                             }
		if ( V .belongs_to ( side_3 ) )
		{	assert ( W.belongs_to ( side_1 ) );
			assert ( V == V35 );
			assert ( W == V16 );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_W = corresp_ver.find ( W );
			assert ( it_W == corresp_ver.end() );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V35 = corresp_ver.find ( V35 );
			assert ( it_V35 != corresp_ver.end() );
			assert ( it_V35->second .second == 0 );
			Cell new_V35 = it_V35->second.first;
			// inspired in item 24 of the book : Scott Meyers, Effective STL
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V16 = corresp_ver.lower_bound ( V16 );
			assert ( ( it_V16 == corresp_ver.end() ) or
		           ( corresp_ver.key_comp()(V16,it_V16->first) ) );
			corresp_ver.emplace_hint ( it_V16, std::piecewise_construct,
		      std::forward_as_tuple ( V16 ), std::forward_as_tuple
		      ( std::pair < Cell, Manifold::Action > { new_V35, g56 } ) );
			// corresp_ver [ V16 ] = { new_V35, g56 };
			continue;                                                                     }
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_V = corresp_ver.find ( V );
		assert ( it_V != corresp_ver.end() );
		assert ( it_V->second.second == 0 );
		Cell new_V = it_V->second.first;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_W = corresp_ver.lower_bound ( W );
		assert ( ( it_W == corresp_ver.end() ) or
		         ( corresp_ver.key_comp()(W,it_W->first) ) );
		corresp_ver.emplace_hint ( it_W, std::piecewise_construct,
		    std::forward_as_tuple ( W ), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { new_V, g56 } ) );  }
		// corresp_ver [ W ] = { new_V, g56 };
	assert ( not it6.in_range() );

	return fold_common_six_sides
		( *that, corresp_ver, side_1, side_2, side_3, side_4, side_5, side_6, dim, keep_map, m );

}  // end of fold_six_sides with tag::build_new_vertices
	
//----------------------------------------------------------------------------------//

	
Mesh fold_six_sides ( Mesh * that, const tag::Identify &, const Mesh & side_1,
                      const tag::With &, const Mesh & side_2,
                  const tag::Identify &, const Mesh & side3,
                      const tag::With &, const Mesh & side4,
                  const tag::Identify &, const Mesh & side5,
                      const tag::With &, const Mesh & side6,
                  const tag::UseExistingVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
									size_t dim, bool keep_map, std::map < Cell, Cell > & m )

// take a mesh having as external boundary a hexagon
// and identify three pairs of opposite sides

// a quotient manifold with two action generators will be built
		
{	// we use the current manifold
	Manifold space = Manifold::working;
	assert ( space.exists() );
	Function coord = space.coordinates();
	assert ( coord.nb_of_components() == 2 );
	// we split 'coord' into its components
	Function x = coord[0],  y = coord[1];

	// we need to identify two translations, one which moves side_1 into side_2
	// and another one which moves side_3 into side_4

	Cell A = side_1.first_vertex().reverse();
	Cell B = side_1.last_vertex();
	Cell C = side_2.last_vertex();
	Cell D = side_2.first_vertex().reverse();

	double dx12 = x(D) - x(A), dy12 = y(D) - y(A);
	double norm12 = std::sqrt ( dx12*dx12 + dy12*dy12 );
	assert ( std::abs ( dx12 - ( x(C) - x(B) ) ) < 1.e-4 * norm12 );
	assert ( std::abs ( dy12 - ( y(C) - y(B) ) ) < 1.e-4 * norm12 );

	Mesh side_3 ( tag::non_existent ), side_4 ( tag::non_existent );
	// we want side_1 and side_3 to touch
	if ( A.belongs_to ( side3 ) or B.belongs_to ( side3 ) )
	{	side_3 = side3; side_4 = side4;  }
	else
	{	assert ( A.belongs_to ( side4 ) or B.belongs_to ( side4 ) );
		side_3 = side4; side_4 = side3;                                 }

	A = side_3.first_vertex().reverse();
	B = side_3.last_vertex();
	C = side_4.last_vertex();
	D = side_4.first_vertex().reverse();

	double dx34 = x(D) - x(A), dy34 = y(D) - y(A);
	double norm34 = std::sqrt ( dx34*dx34 + dy34*dy34 );
	assert ( std::abs ( dx34 - ( x(C) - x(B) ) ) < 1.e-4 * norm34 );
	assert ( std::abs ( dy34 - ( y(C) - y(B) ) ) < 1.e-4 * norm34 );

	Mesh side_5 ( tag::non_existent ), side_6 ( tag::non_existent );
	// we want side_3 and side_5 to touch
	if ( A.belongs_to ( side5 ) or B.belongs_to ( side5 ) )
	{	side_5 = side5; side_6 = side6;  }
	else
	{	assert ( A.belongs_to ( side6 ) or B.belongs_to ( side6 ) );
		side_5 = side6; side_6 = side5;                                 }

	A = side_5.first_vertex().reverse();
	B = side_5.last_vertex();
	C = side_6.last_vertex();
	D = side_6.first_vertex().reverse();

	double dx56 = x(D) - x(A), dy56 = y(D) - y(A);
	double norm56 = std::sqrt ( dx56*dx56 + dy56*dy56 );
	assert ( std::abs ( dx56 - ( x(C) - x(B) ) ) < 1.e-4 * norm56 );
	assert ( std::abs ( dy56 - ( y(C) - y(B) ) ) < 1.e-4 * norm56 );

	// we are confident that side_3 lies between side_1 and side_5
	// orientations may vary, we do not care
	// likewise, side_4 lies between side_2 and side_6

	// d56 is either  d12 + d34  or  d12 - d34  or  -d12 - d34  or  -d12 + d34
	std::vector < std::vector < double > > signs { { 1, 1 }, { 1, -1 }, { -1, -1 }, { -1, 1 } };
	size_t index_min_dif = 0;
	double min_dif = ( std::abs ( dx56 - signs[0][0]*dx12 - signs[0][1]*dx34 ) +
	                   std::abs ( dy56 - signs[0][0]*dy12 - signs[0][1]*dy34 )   );
	for ( size_t i = 1; i < 4; i++ )
	{	double dif = ( std::abs ( dx56 - signs[i][0]*dx12 - signs[i][1]*dx34 ) +
		               std::abs ( dy56 - signs[i][0]*dy12 - signs[i][1]*dy34 )   );
		if ( dif < min_dif )  {  min_dif = dif;  index_min_dif = i;  }              }
	assert ( min_dif < 1.e-4 * norm56 );

	// the desired translations are ( dx12, dy12 ), ( dx34, dy34 )
	Manifold::Action g12 ( tag::transforms, coord, tag::into, (x+dx12) && ( y+dy12) );
	Manifold::Action g34 ( tag::transforms, coord, tag::into, (x+dx34) && ( y+dy34) );
	Manifold manif_q = space.quotient ( g12, g34 );

	// a third translation, not used in the definition of the quotient manifold
	Manifold::Action g56 = signs[index_min_dif][0]*g12 +
	                                     signs[index_min_dif][1]*g34  ;
	
	// we use a map -- for a faster code, we could use Cell::Core::hook
	std::map < Cell, std::pair < Cell, Manifold::Action > > corresp_ver;

	Mesh::Iterator it_ver = that->iterator ( tag::over_vertices );
	for ( it_ver.reset(); it_ver.in_range(); it_ver++ )
	{	Cell V = *it_ver;
		if ( V.belongs_to ( side_2 ) ) continue;
		if ( V.belongs_to ( side_4 ) ) continue;
		if ( V.belongs_to ( side_6 ) ) continue;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_map = corresp_ver.lower_bound ( V );
		assert ( ( it_map == corresp_ver.end() ) or
		         ( corresp_ver.key_comp()(V,it_map->first) ) );
		corresp_ver.emplace_hint ( it_map, std::piecewise_construct,
		    std::forward_as_tuple ( V ), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { V, 0 } ) );  }
		// corresp_ver [ V ] = { V, 0 };

	Cell V16 ( tag::non_existent ), V24 ( tag::non_existent );
	Cell V13 ( tag::non_existent ), V25 ( tag::non_existent );
	assert ( side_1.number_of ( tag::segments ) == side_2.number_of ( tag::segments ) );
	Mesh::Iterator it1 = side_1.iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it2 = side_2.iterator ( tag::over_vertices, tag::require_order );
	for ( it1.reset(), it2.reset(); it1.in_range(); it1++, it2++ )
	{	assert ( it2.in_range() );
		Cell V = *it1;  Cell W = *it2;
		assert ( std::abs ( dx12 - ( x(W) - x(V) ) ) < 1.e-4 * norm12 );
		assert ( std::abs ( dy12 - ( y(W) - y(V) ) ) < 1.e-4 * norm12 );
		if ( V.belongs_to ( side_3 ) )
		{	assert ( W.belongs_to ( side_5 ) );
			V13 = V;  V25 = W;               }
		if ( V.belongs_to ( side_6 ) )
		{	assert ( W.belongs_to ( side_4 ) );
			V16 = V;  V24 = W;  continue;    }
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_W = corresp_ver.lower_bound ( W );
		assert ( ( it_W == corresp_ver.end() ) or
		         ( corresp_ver.key_comp()(W,it_W->first) ) );
		corresp_ver.emplace_hint ( it_W, std::piecewise_construct,
		    std::forward_as_tuple ( W ), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { V, g12 } ) );  }
		// corresp_ver [ W ] = { V, g12 };
	assert ( not it2.in_range() );
	assert ( V16.exists() );  assert ( V24.exists() );
	assert ( V13.exists() );  assert ( V25.exists() );

	Cell V35 ( tag::non_existent ), V46 ( tag::non_existent );
	assert ( side_3.number_of ( tag::segments ) == side_4.number_of ( tag::segments ) );
	Mesh::Iterator it3 = side_3.iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it4 = side_4.iterator ( tag::over_vertices, tag::require_order );
	for ( it3.reset(), it4.reset(); it3.in_range(); it3++, it4++ )
	{	assert ( it4.in_range() );
		Cell V = *it3;  Cell W = *it4;
		assert ( std::abs ( dx34 - ( x(W) - x(V) ) ) < 1.e-4 * norm34 );
		assert ( std::abs ( dy34 - ( y(W) - y(V) ) ) < 1.e-4 * norm34 );
		if ( V .belongs_to ( side_5 ) )
		{	assert ( W.belongs_to ( side_2 ) );
			V35 = V;
			assert ( W == V24 );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_W = corresp_ver.find ( W );
			assert ( it_W == corresp_ver.end() );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V = corresp_ver.find ( V );
			assert ( it_V != corresp_ver.end() );
			assert ( it_V->second .second == 0 );                                }
		if ( V .belongs_to ( side_1 ) )
		{	assert ( W.belongs_to ( side_6 ) );
			V46 = W;
			assert ( V == V13 );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_W = corresp_ver.find ( W );
			assert ( it_W == corresp_ver.end() );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V = corresp_ver.find ( V );
			assert ( it_V != corresp_ver.end() );
			assert ( it_V->second .second == 0 );                                }
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_W = corresp_ver.lower_bound ( W );
		assert ( ( it_W == corresp_ver.end() ) or
		         ( corresp_ver.key_comp()(W,it_W->first) ) );
		corresp_ver.emplace_hint ( it_W, std::piecewise_construct,
		    std::forward_as_tuple ( W ), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { V, g34 } ) );  }
		// corresp_ver [ W ] = { V, g34 };
	assert ( not it4.in_range() );
	assert ( V35.exists() );  assert ( V46.exists() );

	assert ( side_5.number_of ( tag::segments ) == side_6.number_of ( tag::segments ) );
	Mesh::Iterator it5 = side_5.iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it6 = side_6.iterator ( tag::over_vertices, tag::require_order );
	for ( it5.reset(), it6.reset(); it5.in_range(); it5++, it6++ )
	{	assert ( it6.in_range() );
		Cell V = *it5;  Cell W = *it6;
		assert ( std::abs ( dx56 - ( x(W) - x(V) ) ) < 1.e-4 * norm56 );
		assert ( std::abs ( dy56 - ( y(W) - y(V) ) ) < 1.e-4 * norm56 );
		if ( V .belongs_to ( side_2 ) )
		{	assert ( W.belongs_to ( side_4 ) );
			assert ( V == V25 );
			assert ( W == V46 );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_W = corresp_ver.find ( W );
			assert ( it_W != corresp_ver.end() );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V13 = corresp_ver.find ( V13 );
			assert ( it_V13 != corresp_ver.end() );
			assert ( it_V13->second .first == V13 );
			assert ( it_V13->second .second == 0 );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V = corresp_ver.find ( V );
			assert ( it_V != corresp_ver.end() );
			assert ( it_V->second .first == V13 );
			assert ( it_V->second .second == g12 );
			continue;                                                             }
		if ( V .belongs_to ( side_3 ) )
		{	assert ( W.belongs_to ( side_1 ) );
			assert ( V == V35 );
			assert ( W == V16 );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_W = corresp_ver.find ( W );
			assert ( it_W == corresp_ver.end() );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V35 = corresp_ver.find ( V35 );
			assert ( it_V35 != corresp_ver.end() );
			assert ( it_V35->second .first == V35 );
			assert ( it_V35->second .second == 0 );
			// inspired in item 24 of the book : Scott Meyers, Effective STL
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V16 = corresp_ver.lower_bound ( V16 );
			assert ( ( it_V16 == corresp_ver.end() ) or
		           ( corresp_ver.key_comp()(V16,it_V16->first) ) );
			corresp_ver.emplace_hint ( it_V16, std::piecewise_construct,
		      std::forward_as_tuple ( V16 ), std::forward_as_tuple
		      ( std::pair < Cell, Manifold::Action > { V35, g56 } ) );
			// corresp_ver [ V16 ] = { V35, g56 };
			continue;                                                                     }
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_W = corresp_ver.lower_bound ( W );
		assert ( ( it_W == corresp_ver.end() ) or
		         ( corresp_ver.key_comp()(W,it_W->first) ) );
		corresp_ver.emplace_hint ( it_W, std::piecewise_construct,
		    std::forward_as_tuple ( W ), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { V, g56 } ) );  }
		// corresp_ver [ W ] = { V, g56 };
	assert ( not it6.in_range() );

	return fold_common_six_sides
		( *that, corresp_ver, side_1, side_2, side_3, side_4, side_5, side_6, dim, keep_map, m );

}  // end of fold_six_sides with tag::use_existing_vertices

}  // anonymous namespace

//----------------------------------------------------------------------------------//


// take a mesh and fold it around the current working manifold,
// which must be a quotient manifold

Mesh Mesh::fold ( const tag::BuildNewVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
									size_t dim, std::map < Cell, Cell > & m                 )

{	return fold_no_sides ( this, tag::build_new_vertices, tag::return_map_between,
											tag::cells_of_dim, dim, true, m                     );  }


Mesh Mesh::fold ( const tag::BuildNewVertices & )

{	std::map < Cell, Cell > m;
	return fold_no_sides ( this, tag::build_new_vertices, tag::return_map_between,
											tag::cells_of_dim, 2, false, m                             );  }

	
Mesh Mesh::fold ( const tag::UseExistingVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
									size_t dim, std::map < Cell, Cell > & m                 )

{	return fold_no_sides ( this, tag::use_existing_vertices, tag::return_map_between,
											tag::cells_of_dim, dim, true, m                     );  }
	
	
Mesh Mesh::fold ( const tag::UseExistingVertices & )

{	std::map < Cell, Cell > m;
	return fold_no_sides ( this, tag::use_existing_vertices, tag::return_map_between,
											tag::cells_of_dim, 2, false, m                               );  }

//----------------------------------------------------------------------------------//

	
// take a mesh whose external boundary has two parallel segments
// and identify these two segments

// a quotient manifold will be built with one action generator
	

Mesh Mesh::fold ( const tag::Identify &, const Mesh & side_1,
                      const tag::With &, const Mesh & side_2,
                  const tag::BuildNewVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
									size_t dim, std::map < Cell, Cell > & m                 )

{	return fold_two_sides ( this, tag::identify, side_1, tag::with, side_2,
											tag::build_new_vertices, tag::return_map_between,
											tag::cells_of_dim, dim, true, m                     );  }
	
	
Mesh Mesh::fold ( const tag::Identify &, const Mesh & side_1,
                      const tag::With &, const Mesh & side_2,
                  const tag::BuildNewVertices &           )

{	std::map < Cell, Cell > m;
	return fold_two_sides ( this, tag::identify, side_1, tag::with, side_2,
											tag::build_new_vertices, tag::return_map_between,
											tag::cells_of_dim, 2, false, m                       );  }
	
	
Mesh Mesh::fold ( const tag::Identify &, const Mesh & side_1,
                      const tag::With &, const Mesh & side_2,
                  const tag::UseExistingVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
									size_t dim, std::map < Cell, Cell > & m                 )

{	return fold_two_sides ( this, tag::identify, side_1, tag::with, side_2,
											tag::use_existing_vertices, tag::return_map_between,
											tag::cells_of_dim, dim, true, m                     );  }

	
Mesh Mesh::fold ( const tag::Identify &, const Mesh & side_1,
                      const tag::With &, const Mesh & side_2,
                  const tag::UseExistingVertices &           )

{	std::map < Cell, Cell > m;
	return fold_two_sides ( this, tag::identify, side_1, tag::with, side_2,
											tag::use_existing_vertices, tag::return_map_between,
											tag::cells_of_dim, 2, false, m                       );  }

//----------------------------------------------------------------------------------//

	
// take a mesh having as external boundary a parallelogram
// and identify two pairs of opposite sides

// a quotient manifold with two action generators will be built


Mesh Mesh::fold ( const tag::Identify &, const Mesh & side_1,
                      const tag::With &, const Mesh & side_2,
                  const tag::Identify &, const Mesh & side_3,
                      const tag::With &, const Mesh & side_4,
                  const tag::BuildNewVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
									size_t dim, std::map < Cell, Cell > & m                 )

{	return fold_four_sides ( this, tag::identify, side_1, tag::with, side_2,
											tag::identify, side_3, tag::with, side_4,
											tag::build_new_vertices, tag::return_map_between,
											tag::cells_of_dim, dim, true, m                     );  }
	
	
Mesh Mesh::fold ( const tag::Identify &, const Mesh & side_1,
                      const tag::With &, const Mesh & side_2,
                  const tag::Identify &, const Mesh & side_3,
                      const tag::With &, const Mesh & side_4,
                  const tag::BuildNewVertices &           )

{	std::map < Cell, Cell > m;
	return fold_four_sides ( this, tag::identify, side_1, tag::with, side_2,
											tag::identify, side_3, tag::with, side_4,
											tag::build_new_vertices, tag::return_map_between,
											tag::cells_of_dim, 2, false, m                       );  }

	
Mesh Mesh::fold ( const tag::Identify &, const Mesh & side_1,
                      const tag::With &, const Mesh & side_2,
                  const tag::Identify &, const Mesh & side_3,
                      const tag::With &, const Mesh & side_4,
                  const tag::UseExistingVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
									size_t dim, std::map < Cell, Cell > & m                 )

{	return fold_four_sides ( this, tag::identify, side_1, tag::with, side_2,
											tag::identify, side_3, tag::with, side_4,
											tag::use_existing_vertices, tag::return_map_between,
											tag::cells_of_dim, dim, true, m                     );  }
	
	
Mesh Mesh::fold ( const tag::Identify &, const Mesh & side_1,
                      const tag::With &, const Mesh & side_2,
                  const tag::Identify &, const Mesh & side_3,
                      const tag::With &, const Mesh & side_4,
                  const tag::UseExistingVertices &           )

{	std::map < Cell, Cell > m;
	return fold_four_sides ( this, tag::identify, side_1, tag::with, side_2,
											tag::identify, side_3, tag::with, side_4,
											tag::use_existing_vertices, tag::return_map_between,
											tag::cells_of_dim, 2, false, m                       );  }

//----------------------------------------------------------------------------------//


// take a mesh having as external boundary a hexagon
// and identify three pairs of opposite sides

// a quotient manifold with two action generators will be built
	

Mesh Mesh::fold ( const tag::Identify &, const Mesh & side_1,
                      const tag::With &, const Mesh & side_2,
                  const tag::Identify &, const Mesh & side_3,
                      const tag::With &, const Mesh & side_4,
                  const tag::Identify &, const Mesh & side_5,
                      const tag::With &, const Mesh & side_6,
                  const tag::BuildNewVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
									size_t dim, std::map < Cell, Cell > & m                 )

{	return fold_six_sides ( this, tag::identify, side_1, tag::with, side_2,
											tag::identify, side_3, tag::with, side_4,
											tag::identify, side_5, tag::with, side_6,
											tag::build_new_vertices, tag::return_map_between,
											tag::cells_of_dim, dim, true, m                     );  }

	
Mesh Mesh::fold ( const tag::Identify &, const Mesh & side_1,
                      const tag::With &, const Mesh & side_2,
                  const tag::Identify &, const Mesh & side_3,
                      const tag::With &, const Mesh & side_4,
                  const tag::Identify &, const Mesh & side_5,
                      const tag::With &, const Mesh & side_6,
                  const tag::BuildNewVertices &              )

{	std::map < Cell, Cell > m;
	return fold_six_sides ( this, tag::identify, side_1, tag::with, side_2,
											tag::identify, side_3, tag::with, side_4,
											tag::identify, side_5, tag::with, side_6,
											tag::build_new_vertices, tag::return_map_between,
											tag::cells_of_dim, 2, false, m                       );  }

	
Mesh Mesh::fold ( const tag::Identify &, const Mesh & side_1,
                      const tag::With &, const Mesh & side_2,
                  const tag::Identify &, const Mesh & side_3,
                      const tag::With &, const Mesh & side_4,
                  const tag::Identify &, const Mesh & side_5,
                      const tag::With &, const Mesh & side_6,
                  const tag::UseExistingVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
									size_t dim, std::map < Cell, Cell > & m                 )

{	return fold_six_sides ( this, tag::identify, side_1, tag::with, side_2,
											tag::identify, side_3, tag::with, side_4,
											tag::identify, side_5, tag::with, side_6,
											tag::use_existing_vertices, tag::return_map_between,
											tag::cells_of_dim, dim, true, m                     );  }
  
	
Mesh Mesh::fold ( const tag::Identify &, const Mesh & side_1,
                      const tag::With &, const Mesh & side_2,
                  const tag::Identify &, const Mesh & side_3,
                      const tag::With &, const Mesh & side_4,
                  const tag::Identify &, const Mesh & side_5,
                      const tag::With &, const Mesh & side_6,
                  const tag::UseExistingVertices &              )

{	std::map < Cell, Cell > m;
	return fold_six_sides ( this, tag::identify, side_1, tag::with, side_2,
											tag::identify, side_3, tag::with, side_4,
											tag::identify, side_5, tag::with, side_6,
											tag::use_existing_vertices, tag::return_map_between,
											tag::cells_of_dim, 2, false, m                       );  }

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
	Mesh::Iterator it = this->iterator ( tag::over_vertices );
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
	Mesh::Iterator it = this->iterator ( tag::over_segments );
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
	assert ( n == manif_q->winding_nbs.size() );
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

	// here, the action group has one generator
	assert ( manif_q->actions.size() == 1 );
	assert ( manif_q->winding_nbs.size() == 1 );
	Manifold::Action g = manif_q->actions[0];

	std::ofstream file_ps ( file_name );
	file_ps << "please copy here the preamble from the end of file - after %EOF " << std::endl;
	file_ps	<< "you may also want to define a different contour path" << std::endl;
	file_ps	<< "please erase the preamble from the end of file; erase also these three lines";
	file_ps << std::endl << std::endl;
						 
	double xmin = 1.e8, xmax = -1.e8, ymin = 1.e8, ymax = -1.e8, maxside;

	Cell shadow ( tag::vertex );
	std::vector < double > coords_base, coords_tip;

	Mesh::Iterator it = this->iterator ( tag::over_segments );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell seg = *it;
		Cell base = seg .base() .reverse();
		Cell tip  = seg .tip();
		short int iii = 0;  // is 0 at the first passage, will be 1 at the second passage
		for ( short int dir = -1; dir < 2; dir += 2 )  // dir will be -1 and 1
		{	size_t first_unsuccessful_tries = 1, last_unsuccessful_tries = 0;
			short int ii = iii;
			iii++;  // iii was 0 at the first passage, will be 1 at the second passage
			while ( true )
			{	bool successful_round = false;
				Manifold::Action a = ii*g;
				bool touches_region = false;
				coords_base = coords_q ( base, tag::winding, a );
				coords_Eu ( shadow ) = coords_base;
				touches_region = touches_region or constraints.on_cell ( shadow );
				a += seg.winding();
				coords_tip = coords_q ( tip, tag::winding, a );
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
					if ( yy > ymax ) ymax = yy;                                    }
				if ( successful_round )
				{	first_unsuccessful_tries = 0;
					last_unsuccessful_tries = 0;   }
				else  // either we have not started yet, or we are approaching the end
				{	if ( first_unsuccessful_tries > 0 )  // not started yet
					{	first_unsuccessful_tries++;
						if ( first_unsuccessful_tries > 50 ) break;  }
					else  // approaching end
					{	last_unsuccessful_tries++;
						if ( last_unsuccessful_tries > 10 ) break;   }   }
				ii += dir;                                                          }  }  }
	
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

	file_ps << "gsave" << std::endl;
	// file_ps << "/m{moveto}def" << std::endl;
	// file_ps << "/l{lineto}def" << std::endl;
	// file_ps << "/s{stroke}def" << std::endl;
	file_ps << translation_x + scale_factor*border << " "
	        << translation_y + scale_factor*border << " translate" << std::endl;
	file_ps << scale_factor << " dup scale" << std::endl << std::endl;

  file_ps << "0 setgray contour clip" << std::endl;
  file_ps << 1.5 / scale_factor << " setlinewidth 0 setgray" << std::endl;
	
	if ( ! file_ps.good() )
	{	std::cerr << "error writing postscript file" << std::endl;
		exit (1);                                                   }

} // end of  Mesh::draw_ps with tag::unfold and tag::one_generator

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
	Function x = coords_Eu [0],  y = coords_Eu [1];

	// here, the action group has two generators
	assert ( manif_q->actions.size() == 2 );
	assert ( manif_q->winding_nbs.size() == 2 );
	Manifold::Action g1 = manif_q->actions [0], g2 = manif_q->actions [1];
	
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

	Mesh::Iterator it = this->iterator ( tag::over_segments );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell seg = *it;
		Cell base = seg .base() .reverse();
		Cell tip  = seg .tip();
		// we describe a sort of spiral
		// if the first tries are out of the region, we give up after 50 unsuccsessful rounds
		// at the end, we stop after 10 unsuccessful rounds
		size_t first_unsuccessful_tries = 1, last_unsuccessful_tries = 0;
		size_t size_of_round = 0;
		short int ii = 0, jj = 0;
		while ( true )
		{	size_of_round++;
			bool successful_round = false;
			for ( size_t d = 0; d < 4; d++ )
			{	if ( d == 2 ) size_of_round++;
				for ( size_t i = 0; i < size_of_round; i++ )
				{	Manifold::Action a = ii*g1 + jj*g2;
					bool touches_region = false;
					coords_base = coords_q ( base, tag::winding, a );
					coords_Eu ( shadow ) = coords_base;
					touches_region = touches_region or constraints.on_cell ( shadow );
					a += seg.winding();
					coords_tip = coords_q ( tip, tag::winding, a );
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
	coords_base = coords_q ( shadow, tag::winding, g1 );
	file_ps << coords_base[0] << " " << coords_base[1] << " rlineto ";
	coords_base = coords_q ( shadow, tag::winding, g2 );
	file_ps << coords_base[0] << " " << coords_base[1] << " rlineto ";
	coords_base = coords_q ( shadow, tag::winding, -g1 );
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
         const tag::Unfold &, const std::vector < Manifold::Action > & v,
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

	std::ofstream file_ps ( file_name );
	file_ps << "please copy here the preamble from the end of file - after %EOF " << std::endl;
	file_ps	<< "you can move the shadow by simply changing the initial point" << std::endl;
	file_ps	<< "you may also want to define a different contour path" << std::endl;
	file_ps	<< "please erase the preamble from the end of file; erase also these four lines";
	file_ps << std::endl << std::endl;
						 
	double xmin = 1.e8, xmax = -1.e8, ymin = 1.e8, ymax = -1.e8, maxside;

	Cell shadow ( tag::vertex );
	std::vector < double > coords_base, coords_tip;

	Mesh::Iterator it = this->iterator ( tag::over_segments );
	for ( it.reset() ; it.in_range(); it++ )
	{	Cell seg = *it;
		Cell base = seg.base().reverse();
		Cell tip  = seg.tip();
		for ( std::vector < Manifold::Action > ::const_iterator
						it_v = v.begin(); it_v != v.end(); it_v++                     )
		{	Manifold::Action a = *it_v;
			bool touches_region = false;
			coords_base = coords_q ( base, tag::winding, a );
			coords_Eu ( shadow ) = coords_base;
			touches_region = touches_region or constraints.on_cell ( shadow );
			a += seg.winding();
			coords_tip = coords_q ( tip, tag::winding, a );
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
	coords_base = coords_q ( shadow, tag::winding, 0 );
	file_ps << coords_base[0] << " " << coords_base[1] << " rlineto ";
	coords_base = coords_q ( shadow, tag::winding, 0 );
	file_ps << coords_base[0] << " " << coords_base[1] << " rlineto ";
	coords_base = coords_q ( shadow, tag::winding, 0 );
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
	Mesh::Iterator it = this->iterator ( tag::over_vertices );
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
	Mesh::Iterator it = this->iterator ( tag::over_segments );
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

// 'numb_map' should begin at 0
// we add 1 to each number because gmsh seems to prefer numbers to begin at 1

{	// we use the current manifold
	Manifold space = Manifold::working;
	assert ( space.exists() );
	Function coord = space.coordinates();

	std::ofstream file_msh (f);
	file_msh << "$MeshFormat" << std::endl << "2.2 0 8" << std::endl;
	file_msh << "$EndMeshFormat" << std::endl;
	
	file_msh << "$Nodes" << std::endl << this->number_of(tag::cells_of_dim,0) << std::endl;

	{ // just to make variables local : it, counter, x, y
	Mesh::Iterator it = this->iterator ( tag::over_vertices );
	Function x = coord[0], y = coord[1];
	if (coord.nb_of_components() == 2)
	{	for ( it.reset() ; it.in_range(); it++ )
		{	Cell p = *it;
			file_msh << ver_numbering [p] + 1 << " "
		           << x(p) << " " << y(p) << " " << 0 << std::endl;  }  }
	else
	{	assert  ( coord.nb_of_components() == 3 );
		Function z = coord[2];
		for ( it.reset() ; it.in_range(); it++ )
		{	Cell p = *it;
			file_msh << ver_numbering [p] + 1 << " "
		           << x(p) << " " << y(p) << " " << z(p) << std::endl;  }  }
	file_msh << "$EndNodes" << std::endl;
	} // just to make variables local : it, counter, x, y

	file_msh << "$Elements" << std::endl;
	file_msh << this->number_of ( tag::cells_of_dim, this->dim() ) << std::endl;

	if ( this->dim() == 1 )
	{	Mesh::Iterator it = this->iterator ( tag::over_segments );
		size_t counter = 0;
		for ( it.reset() ; it.in_range(); it++)
		{	++counter;
			Cell elem = *it;
			file_msh << counter << " 1 0 ";
			Cell A = elem.base().reverse();
			file_msh << ver_numbering [A] + 1 << " ";
			Cell B = elem.tip();
			file_msh << ver_numbering [B] + 1 << std::endl;    }  }
	else if ( this->dim() == 2 )
	{	Mesh::Iterator it = this->iterator ( tag::over_cells_of_dim, 2 );
		size_t counter = 0;
		for ( it.reset() ; it.in_range(); it++)
		{	++counter;
			Cell elem = *it;
			if ( elem.boundary().number_of ( tag::cells_of_dim, 1 ) == 3 ) // a triangle
				file_msh << counter << " 2 0 ";
			else // a quadrilateral
			{	assert ( elem.boundary().number_of ( tag::cells_of_dim, 1 ) == 4 );
				file_msh << counter << " 3 0 ";                                     }
			Mesh::Iterator itt = elem.boundary().iterator ( tag::over_vertices, tag::require_order );
			for ( itt.reset(); itt.in_range(); ++itt )
			{	Cell p = *itt;  file_msh << ver_numbering [p] + 1 << " ";   }
			file_msh << std::endl;                                                                 }  }
	else
	{	assert ( this->dim() == 3);
		Mesh::Iterator it = this->iterator ( tag::over_cells_of_dim, 3 );
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
				Mesh::Iterator itt = elem.boundary().iterator ( tag::over_cells_of_dim, 2 );
				itt.reset();  Cell back = *itt; // square face behind the cube
				// back is 0321 in gmsh's documentation
				assert ( back.boundary().number_of ( tag::cells_of_dim, 1 ) == 4 );
				Mesh::Iterator itv = back.boundary().iterator ( tag::over_vertices, tag::backwards );
				// backwards because we want the vertices ordered as 0, 1, 2, 3
				itv.reset();  Cell ver_0 = *itv;
				Cell seg_03 = back.boundary().cell_in_front_of(ver_0);
				for ( ; itv.in_range(); ++itv )
				{	Cell p = *itv;  file_msh << ver_numbering [p] + 1 << " ";   }
				Cell left_wall = elem.boundary().cell_in_front_of(seg_03); // square face on the left
				// left_wall is 0473 in gmsh's documentation
				assert ( left_wall.boundary().number_of ( tag::cells_of_dim, 1 ) == 4 );
				Cell seg_04 = left_wall.boundary().cell_in_front_of(ver_0);
				Cell ver_4 = seg_04.tip();
				Cell seg_47 = left_wall.boundary().cell_in_front_of(ver_4);
				Cell front = elem.boundary().cell_in_front_of(seg_47); // square face in front
				// front is 4567 in gmsh's documentation
				assert ( front.boundary().number_of ( tag::cells_of_dim, 1 ) == 4 );
				Mesh::Iterator itvv = front.boundary().iterator ( tag::over_vertices, tag::require_order );
				itvv.reset ( tag::start_at, ver_4 );
				for ( ; itvv.in_range(); ++itvv )
				{	Cell p = *itvv;  file_msh << ver_numbering [p] + 1 << " ";   }                }
			else
			{	assert( n_faces == 5 );
				// triangular prism = 6-node prism
				// see http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
				// and http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
				file_msh << counter << " 6 0 ";
				Mesh::Iterator itt = elem.boundary().iterator ( tag::over_cells_of_dim, 2 );
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
				Mesh::Iterator itv = base.boundary().iterator ( tag::over_vertices, tag::backwards );
				// backwards because we want the vertices ordered as 0, 1, 2
				itv.reset();  Cell ver_0 = *itv;
				Cell seg_02 = base.boundary().cell_in_front_of(ver_0);
				for (  ; itv.in_range(); ++itv )
				{	Cell p = *itv;  file_msh << ver_numbering [p] + 1 << " ";   }
				Cell right_wall = elem.boundary().cell_in_front_of(seg_02);
				// right_wall is 0352 in gmsh's documentation
				assert ( right_wall.boundary().number_of ( tag::cells_of_dim, 1 ) == 4 );
				Cell seg_03 = right_wall.boundary().cell_in_front_of(ver_0);
				Cell ver_3 = seg_03.tip();
				Cell seg_35 = right_wall.boundary().cell_in_front_of(ver_3);
				Cell roof = elem.boundary().cell_in_front_of(seg_35);
				// roof is 345 in gmsh's documentation
				assert ( roof.boundary().number_of ( tag::cells_of_dim, 1 ) == 3 );
				Mesh::Iterator itvv = roof.boundary().iterator ( tag::over_vertices, tag::require_order );
				itvv.reset ( tag::start_at, ver_3 );
				for ( ; itvv.in_range(); ++itvv )
				{	Cell p = *itvv;  file_msh << ver_numbering [p] + 1 << " ";  }               }
			file_msh << std::endl;                                                                } }
	file_msh << "$EndElements" << std::endl;
	
	if ( ! file_msh.good() )
	{	std::cerr << "error writing msh file" << std::endl;
		exit (1);                                    }

} // end of Mesh::export_msh


void Mesh::export_msh ( std::string f, std::map < Cell, size_t > & numb_map )
	
// 'numb_map' should begin at 0
// later, 'this->export_msh ( f, numbering )' will add 1 to each number

{	Cell::Numbering::Map numbering ( & numb_map );

	this->export_msh ( f, numbering );

} // end of Mesh::export_msh


void Mesh::export_msh ( std::string f )
	
// the numbering of vertices is produced on-the-fly

// we build a 'numbering' map beginning at 0
// later, 'this->export_msh ( f, numbering )' will add 1 to each number
	
{	Cell::Numbering::Map numbering;

	Mesh::Iterator it = this->iterator ( tag::over_vertices );
	size_t counter = 0;
	for ( it.reset() ; it.in_range(); it++ )
	{	Cell p = *it;  numbering [p] = counter;  ++counter;  }

	this->export_msh ( f, numbering );

} // end of Mesh::export_msh

//----------------------------------------------------------------------------------//


namespace {  // anonymous namespace, mimics static linkage
	
inline Mesh unfold_common ( const Mesh & that,
   const std::map < Cell, std::vector < std::pair < Manifold::Action, Cell > > > & built_ver,
                     const tag::OverRegion &, const Function::Inequality::Set & constraints  )
	
{	Mesh result ( tag::fuzzy, tag::of_dim, 2 );

	// build new segments joinig vertices built previously (by the calling code)
	std::map < Cell, std::vector < std::pair < Manifold::Action, Cell > > > built_seg;
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = that .iterator ( tag::over_segments );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell seg = *it;
		std::pair < std::map < Cell, std::vector < std::pair < Manifold::Action, Cell > > >
			::iterator, bool > itbs =
		built_seg .insert ( std::pair < Cell, std::vector < std::pair < Manifold::Action, Cell > > >
		  ( seg, std::vector < std::pair < Manifold::Action, Cell > > () ) );
		assert ( itbs .second );
		Cell base = seg .base() .reverse();
		Cell tip  = seg .tip();
		std::map < Cell, std::vector < std::pair < Manifold::Action, Cell > > >
			::const_iterator it_bvb = built_ver .find ( base );
		assert ( it_bvb != built_ver .end() );
		const std::vector < std::pair < Manifold::Action, Cell > > & bvb = it_bvb->second;
		for ( std::vector < std::pair < Manifold::Action, Cell > > ::const_iterator
		      itb = bvb .begin(); itb != bvb .end(); itb++                          )
		{	Manifold::Action act_base = itb->first;
			std::map < Cell, std::vector < std::pair < Manifold::Action, Cell > > >
				::const_iterator it_bvt = built_ver .find ( tip );
			assert ( it_bvt != built_ver .end() );
			const std::vector < std::pair < Manifold::Action, Cell > > & bvt = it_bvt->second;
			for ( std::vector < std::pair < Manifold::Action, Cell > > ::const_iterator
			      itt = bvt .begin(); itt != bvt .end(); itt++                          )
				if ( itt->first == itb->first + seg .winding() )  // same action
				{	Cell new_seg ( tag::segment, itb->second .reverse(), itt->second );
					itbs .first->second .push_back 
						( std::pair < Manifold::Action, Cell > ( itb->first, new_seg ) ) ;
					break;                                                              }  }  }
	} // just a block of code for hiding 'it'

	// build new cells using segments we just built
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = that .iterator ( tag::over_cells_of_max_dim );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell cll = *it;
		Mesh::Iterator it_seg = cll .boundary() .iterator
			( tag::over_segments, tag::require_order );
		it_seg .reset();  assert ( it_seg .in_range() );
		Cell ini_seg = *it_seg;
		Cell A = ini_seg .base() .reverse();
		std::map < Cell, std::vector < std::pair < Manifold::Action, Cell > > >
			::const_iterator it_bvA = built_ver .find (A);
		assert ( it_bvA != built_ver .end() );
		for ( std::vector < std::pair < Manifold::Action, Cell > > ::const_iterator
		      itA = it_bvA->second .begin(); itA != it_bvA->second .end(); itA++   )
		{	Manifold::Action act_base = itA->first, act_seg = act_base;
			bool build_cell = true;
			for ( it_seg .reset ( tag::start_at, ini_seg ); it_seg .in_range(); it_seg++ )	
			{	Cell other_seg = *it_seg;
				bool segment_found = false;
				if ( other_seg .is_positive() )
				{	std::vector < std::pair < Manifold::Action, Cell > > &
						bss = built_seg [ other_seg ];
					for ( std::vector < std::pair < Manifold::Action, Cell > > ::const_iterator
				        ito = bss .begin(); ito != bss .end(); ito++                          )
						if ( ito->first == act_seg )  {  segment_found = true;  break;  }
					act_seg += other_seg .winding();                                              }
				else  // 'other_seg' is negative
				{	act_seg += other_seg .winding();
					std::vector < std::pair < Manifold::Action, Cell > > &
						bssr = built_seg [ other_seg .reverse() ];
					for ( std::vector < std::pair < Manifold::Action, Cell > > ::const_iterator
				        ito = bssr .begin(); ito != bssr .end(); ito++                        )
						if ( ito->first == act_seg )  {  segment_found = true;  break;  }           }
				if ( not segment_found )  {  build_cell = false;  break;  }                       }
			if ( not build_cell ) continue;
			Cell::Positive::HighDim * new_cll_ptr = new Cell::Positive::HighDim
				( tag::whose_boundary_is,
					Mesh ( tag::whose_core_is,
				         new Mesh::Connected::OneDim ( tag::with,
			           cll .boundary() .number_of ( tag::segments ),
	  	             tag::segments, tag::one_dummy_wrapper       ),
	    	     tag::freshly_created                                 ),
					tag::one_dummy_wrapper                                     );
			Cell new_cll ( tag::whose_core_is, new_cll_ptr, tag::freshly_created );
			Cell kept_ver ( tag::non_existent );
			act_seg = act_base;
			for ( it_seg .reset ( tag::start_at, ini_seg ); it_seg .in_range(); it_seg++ )	
			{	Cell seg = *it_seg;
				bool segment_found = false;
				std::vector < std::pair < Manifold::Action, Cell > > ::iterator ito;
				if ( seg .is_positive() )
				{	std::vector < std::pair < Manifold::Action, Cell > > &
						bss = built_seg [ seg ];
					for ( ito = bss .begin(); ito != bss .end(); ito++ )
						if ( ito->first == act_seg )  {  segment_found = true;  break;  }
					assert ( segment_found );
					assert ( ito->second .exists() );
					kept_ver = ito->second .tip();
					ito->second .core->add_to_mesh
						( new_cll .boundary() .core, tag::do_not_bother );
					act_seg += seg .winding();                                           }
				else  // 'seg' is negative
				{	act_seg += seg .winding();
					std::vector < std::pair < Manifold::Action, Cell > > &
						bssr = built_seg [ seg .reverse() ];
					for ( ito = bssr .begin(); ito != bssr .end(); ito++ )
						if ( ito->first == act_seg )  {  segment_found = true;  break;  }
					assert ( segment_found );
					assert ( ito->second .exists() );
					kept_ver = ito->second .tip();
					ito->second .reverse() .core->add_to_mesh
						( new_cll .boundary() .core, tag::do_not_bother );                   }  }
			assert ( kept_ver .exists() );
			new_cll .boundary() .closed_loop ( kept_ver );
			new_cll .add_to_mesh ( result );
		}	 // end of loop over actions
	}  // end of loop over cells
	} // just a block of code for hiding 'it'

	return result;
	
} // end of unfold_common


inline Mesh unfold_local ( const Mesh & that, const tag::OneGenerator &,
               const tag::OverRegion &, const Function::Inequality::Set & constraints,
               const tag::ReturnMapBetween &, const tag::CellsOfDim &,
               size_t dim, std::map < Cell, std::pair < Cell, Manifold::Action > > & mapping,
               bool fill_mapping                                                             )

// if last argument is true,
// besides the unfolded Mesh, return a mapping giving, for each vertex of the new mesh,
// the corresponding vertex in 'that' mesh together with a winding (action)
	
{	Manifold space = Manifold::working;
	assert ( space.exists() );  // we use the current (quotient) manifold
	Manifold::Quotient * manif_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space.core );
	assert ( manif_q );
	Function coords_q = space.coordinates();
	Manifold mani_Eu = manif_q->base_space;  // underlying Euclidian manifold
	Function coords_Eu = mani_Eu.coordinates();
	assert ( coords_Eu.nb_of_components() == 2 );

	// here, the action group has two generators
	assert ( manif_q->actions.size() == 1 );
	assert ( manif_q->winding_nbs.size() == 1 );
	Manifold::Action g = manif_q->actions[0];

	// we begin by unfolding the vertices
	Cell shadow ( tag::vertex );
	std::vector < double > coords_base, coords_tip;

	// build new vertices within the given region
	std::map < Cell, std::vector < std::pair < Manifold::Action, Cell > > > built_ver;
	assert ( dim == 0 );  // 'mapping' required from vertices to vertices
	Mesh::Iterator it = that .iterator ( tag::over_vertices );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell ver = *it;
		std::pair < std::map < Cell, std::vector < std::pair < Manifold::Action, Cell > > >
			::iterator, bool > itbv =
		built_ver .insert ( std::pair < Cell, std::vector < std::pair < Manifold::Action, Cell > > >
		  ( ver, std::vector < std::pair < Manifold::Action, Cell > > () ) );
		assert ( itbv .second );
		short int iii = 0;  // is 0 at the first passage, will be 1 at the second passage
		for ( short int dir = -1; dir < 2; dir += 2 )  // dir will be -1 and 1
		{	size_t first_unsuccessful_tries = 1, last_unsuccessful_tries = 0;
			short int ii = iii;
			iii++;  // iii was 0 at the first passage, will be 1 at the second passage
			while ( true )
			{	bool successful_round = false;
				Manifold::Action a = ii*g;
				std::vector < double > coords_ver = coords_q ( ver, tag::winding, a );
				coords_Eu ( shadow ) = coords_ver;
				if ( constraints .on_cell ( shadow ) )
				{	successful_round = true;
					Cell new_ver ( tag::vertex );  coords_Eu ( new_ver ) = coords_ver;
					if ( fill_mapping )
					mapping .insert ( std::pair < Cell, std::pair < Cell, Manifold::Action > >
					    ( new_ver, std::pair < Cell, Manifold::Action > ( ver, a ) ) );
					itbv .first->second .push_back
						( std::pair < Manifold::Action, Cell > ( a, new_ver ) );         }
				if ( successful_round )
				{	first_unsuccessful_tries = 0;
					last_unsuccessful_tries = 0;   }
				else  // either we have not started yet, or we are approaching the end
				{	if ( first_unsuccessful_tries > 0 )  // not started yet
					{	first_unsuccessful_tries++;
						if ( first_unsuccessful_tries > 50 ) break;  }
					else  // approaching end
					{	last_unsuccessful_tries++;
						if ( last_unsuccessful_tries > 10 ) break;   }   }
				ii += dir;                                                          }  }  }
	
	// the 'mapping' may contain vertices not actually used
	// (some segments may have never been used) but we leave it like this

	Mesh result = unfold_common ( that, built_ver, tag::over_region, constraints );

	mani_Eu .set_as_working_manifold();

	return result;
	
} // end of  unfold_local  with tag::one_generator


inline Mesh unfold_local ( const Mesh & that, const tag::TwoGenerators &,
               const tag::OverRegion &, const Function::Inequality::Set & constraints,
               const tag::ReturnMapBetween &, const tag::CellsOfDim &,
               size_t dim, std::map < Cell, std::pair < Cell, Manifold::Action > > & mapping,
               bool fill_mapping                                                             )
	
// if last argument is true,
// besides the unfolded Mesh, return a mapping giving, for each vertex of the new mesh,
// the corresponding vertex in 'that' mesh together with a winding (action)
	
{	Manifold space = Manifold::working;
	assert ( space.exists() );  // we use the current (quotient) manifold
	Manifold::Quotient * manif_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space.core );
	assert ( manif_q );
	Function coords_q = space.coordinates();
	Manifold mani_Eu = manif_q->base_space;  // underlying Euclidian manifold
	Function coords_Eu = mani_Eu.coordinates();
	assert ( coords_Eu.nb_of_components() == 2 );

	// here, the action group has two generators
	assert ( manif_q->actions.size() == 2 );
	assert ( manif_q->winding_nbs.size() == 2 );
	Manifold::Action g1 = manif_q->actions[0], g2 = manif_q->actions[1];

	std::vector < std::vector < short int > > directions
		{ { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 } };
	// declare global, here and in progressive.cpp Manifold::Type::Quotient::sq_dist
	// and perhaps in other places

	// we begin by unfolding the vertices
	Cell shadow ( tag::vertex );
	std::vector < double > coords_base, coords_tip;

	// build new vertices within the given region
	std::map < Cell, std::vector < std::pair < Manifold::Action, Cell > > > built_ver;
	assert ( dim == 0 );  // 'mapping' required from vertices to vertices
	Mesh::Iterator it = that .iterator ( tag::over_vertices );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell ver = *it;
		std::pair < std::map < Cell, std::vector < std::pair < Manifold::Action, Cell > > >
			::iterator, bool > itbv =
		built_ver .insert ( std::pair < Cell, std::vector < std::pair < Manifold::Action, Cell > > >
		  ( ver, std::vector < std::pair < Manifold::Action, Cell > > () ) );
		assert ( itbv .second );
		// we describe a sort of spiral
		// if the first tries are out of the region, we give up after 50 unsuccsessful rounds
		// at the end, we stop after 10 unsuccessful rounds
		size_t first_unsuccessful_tries = 1, last_unsuccessful_tries = 0;
		size_t size_of_round = 0;
		short int ii = 0, jj = 0;
		while ( true )
		{	size_of_round++;
			bool successful_round = false;
			for ( size_t d = 0; d < 4; d++ )
			{	if ( d == 2 ) size_of_round++;
				for ( size_t i = 0; i < size_of_round; i++ )
				{	Manifold::Action a = ii*g1 + jj*g2;
					std::vector < double > coords_ver = coords_q ( ver, tag::winding, a );
					coords_Eu ( shadow ) = coords_ver;
					if ( constraints .on_cell ( shadow ) )
					{	successful_round = true;
						Cell new_ver ( tag::vertex );  coords_Eu ( new_ver ) = coords_ver;
						if ( fill_mapping )
						mapping .insert ( std::pair < Cell, std::pair < Cell, Manifold::Action > >
						    ( new_ver, std::pair < Cell, Manifold::Action > ( ver, a ) ) );
						itbv .first->second .push_back
							( std::pair < Manifold::Action, Cell > ( a, new_ver ) );         }  
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
	give_up: ;  // next vertex, please
	}  // end of loop over vertices

	// the 'mapping' may contain vertices not actually used
	// (some segments may have never been used) but we leave it like this

	Mesh result = unfold_common ( that, built_ver, tag::over_region, constraints );

	mani_Eu .set_as_working_manifold();

	return result;
	
} // end of  unfold_local  with tag::two_generators


inline Mesh unfold_local ( const Mesh & that,
               const tag::OverRegion &, const Function::Inequality::Set & constraints,
               const tag::ReturnMapBetween &, const tag::CellsOfDim &,
               size_t dim, std::map < Cell, std::pair < Cell, Manifold::Action > > & mapping,
               bool fill_mapping                                                             )
	
// if last argument is true,
// besides the unfolded Mesh, return a mapping giving, for each vertex of the new mesh,
// the corresponding vertex in 'that' mesh together with a winding (action)
	
{	Manifold space = Manifold::working;
	assert ( space.exists() );  // we use the current (quotient) manifold
	Manifold::Quotient * manif_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space.core );
	assert ( manif_q );

	// the action group may have one or two generators
	size_t n = manif_q->actions .size();
	assert ( n == manif_q->winding_nbs .size() );

	if ( n == 1 )
	return unfold_local ( that, tag::one_generator, tag::over_region, constraints,
	                tag::return_map_between, tag::cells_of_dim, 0, mapping, fill_mapping  );
	
	assert ( n == 2 );
	return unfold_local ( that, tag::two_generators, tag::over_region, constraints,
	                tag::return_map_between, tag::cells_of_dim, 0, mapping, fill_mapping  );

} // end of unfold_local

} // end of anonymous namespace

//----------------------------------------------------------------------------------//


Mesh Mesh::unfold ( const tag::OverRegion &, const Function::Inequality::Set & constraints,
                    const tag::ReturnMapBetween &, const tag::CellsOfDim &,
                    size_t dim, std::map < Cell, std::pair < Cell, Manifold::Action > > & mapping )
const
	
// take a mesh and unfold it over a given region of the plane
// so we can export the resulting mesh in msh format

// besides the unfolded Mesh, return a mapping giving, for each vertex of the new mesh,
// the corresponding vertex in 'this' mesh together with a winding (action)
	
{ assert ( dim == 0 );  // 'mapping' required from vertices to vertices

	return unfold_local ( *this, tag::over_region, constraints,
	           tag::return_map_between, tag::cells_of_dim, 0, mapping, true );
	// last argument true means fill 'mapping'
	
} // end of Mesh::unfold


Mesh Mesh::unfold ( const tag::OverRegion &, const Function::Inequality::Set & constraints )
const
	
// take a mesh and unfold it over a given region of the plane
// so we can export the resulting mesh in msh format

{	std::map < Cell, std::pair < Cell, Manifold::Action > > mapping;

	return unfold_local ( *this, tag::over_region, constraints,
	             tag::return_map_between, tag::cells_of_dim, 0, mapping, false  );
	// last argument false means do not bother with 'mapping'
	
} // end of Mesh::unfold


Mesh Mesh::unfold ( const std::vector < tag::Util::Action > &,
	                  const tag::OverRegion &, const Function::Inequality::Set & constraints,
                    const tag::ReturnMapBetween &, const tag::CellsOfDim &,
                    size_t dim, std::map < Cell, std::pair < Cell, Manifold::Action > > & mapping )
const
	
// take a mesh and unfold it over a given region of the plane
// so we can export the resulting mesh in msh format

// besides the unfolded Mesh, return a mapping giving, for each vertex of the new mesh,
// the corresponding vertex in 'this' mesh together with a winding (action)
	
{ assert ( dim == 0 );  // 'mapping' required from vertices to vertices

	assert ( false );  // not yet implemented
	
	return unfold_local ( *this, tag::over_region, constraints,
	           tag::return_map_between, tag::cells_of_dim, 0, mapping, true );
	// last argument true means fill 'mapping'
	
} // end of Mesh::unfold


Mesh Mesh::unfold ( const std::vector < tag::Util::Action > &,
	                  const tag::OverRegion &, const Function::Inequality::Set & constraints )
const
	
// take a mesh and unfold it over a given region of the plane
// so we can export the resulting mesh in msh format

{	std::map < Cell, std::pair < Cell, Manifold::Action > > mapping;

	assert ( false );  // not yet implemented
	
	return unfold_local ( *this, tag::over_region, constraints,
	             tag::return_map_between, tag::cells_of_dim, 0, mapping, false  );
	// last argument false means do not bother with 'mapping'
	
} // end of Mesh::unfold

//----------------------------------------------------------------------------------//


size_t & Cell::Numbering::Map::operator[] ( const Cell p )  // virtual from Cell::Numbering
{	return (*(this->map))[p];  }
