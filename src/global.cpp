
// global.cpp 2022.02.13

//   This file is part of maniFEM, a C++ library for meshes and finite elements on manifolds.

//   Copyright 2019, 2020, 2021, 2022 Cristian Barbarosie cristian.barbarosie@gmail.com

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

#include <fstream>
#include "maniFEM.h"

using namespace maniFEM;

//----------------------------------------------------------------------------------//

//  https://stackoverflow.com/questions/216823/how-to-trim-a-stdstring

#include <algorithm> 
#include <cctype>
#include <locale>

// trim from start (in place)
static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

// trim from start (copying)
static inline std::string ltrim_copy(std::string s) {
    ltrim(s);
    return s;
}

// trim from end (copying)
static inline std::string rtrim_copy(std::string s) {
    rtrim(s);
    return s;
}

// trim from both ends (copying)
static inline std::string trim_copy(std::string s) {
    trim(s);
    return s;
}

//----------------------------------------------------------------------------------//


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

void build_common                                                   // line 281
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
				// we chose that S will have zero winding relatively to A
				// but we have no control on the winding numbers of the six vertices
				// used for interpolating coordinates, so we use six shadow vertices
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

//----------------------------------------------------------------------------------//


Cell Mesh::common_vertex ( const tag::With &, const Mesh & seg2, const tag::ExactlyOne & ) const

// we look for a common vertex, where 'this' ends and 'seg2' begins (in this order)
// this does not apply to closed loops of course, so we give them a different treatment

{	std::vector < Cell > vec;
	Mesh::Iterator it = this->iterator ( tag::over_vertices );
	it .reset();
	assert ( it .in_range() );
	Cell W = *it;
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;
		if ( V .belongs_to ( seg2 ) ) vec .push_back ( V );  }
	assert ( vec .size() == 1 );
	Cell V = vec [0];
	return V;                                                  }

//----------------------------------------------------------------------------------//


Cell Mesh::common_vertex ( const tag::With &, const Mesh & seg2 ) const

// we look for a common vertex, where 'this' ends and 'seg2' begins (in this order)
// this does not apply to closed loops of course, so we give them a different treatment

{	std::vector < Cell > vec;
	Mesh::Iterator it = this->iterator ( tag::over_vertices );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;
		if ( V .belongs_to ( seg2 ) ) vec .push_back ( V );  }
	assert ( vec .size() > 0 );
	if ( vec .size() == 1 )  // one common vertex only, so we have no doubt
		return vec [0];
	assert ( vec .size() == 2 );  // it makes no sense to have more than two
	// we are looking for the one where 'this' ends and 'seg2' begins
	Cell V = vec [0];
	if ( this->cell_in_front_of ( V, tag::may_not_exist ) .exists() )
	{	// 'this' does not end in V, so this is not the vertex we are looking for
		// perhaps the other one ?
		V = vec [1];
		assert ( not this->cell_in_front_of ( V, tag::may_not_exist ) .exists() );
		// 'this' ends in V
		assert ( not seg2 .cell_behind ( V, tag::may_not_exist ) .exists() );
		// seg2 begins in V
		return V;                                                                 }
	// else : 'this' ends in V
	assert ( not seg2 .cell_behind ( V, tag::may_not_exist ) .exists() );
	// seg2 begins in V
	return V;                                                                      }

//----------------------------------------------------------------------------------//


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
	Cell A = CA .common_vertex ( tag::with, AB );
	Cell B = AB .common_vertex ( tag::with, BC );
	Cell C = BC .common_vertex ( tag::with, CA );

	// 'build_common' builds all triangles but the last four, near C
	Manifold::Action winding_C_from_A;
	Cell seg1 ( tag::non_existent ), seg2 ( tag::non_existent );
	build_common ( *this, AB, BC, CA, A, B, C, seg1, seg2, winding_C_from_A, true );  // line 281
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
// one of the segments provided as arguments has O as base and other as tip
	
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

	Cell A = CA .common_vertex ( tag::with, AB );
	Cell B = AB .common_vertex ( tag::with, BC );	

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
	Cell SW = south .first_vertex() .reverse();
	assert ( SW == west .last_vertex() );
	Cell SE = east .first_vertex() .reverse();
	assert ( SE == south .last_vertex() );
	Cell NE = north .first_vertex() .reverse();
	assert ( NE == east .last_vertex() );
	Cell NW = west .first_vertex() .reverse();
	assert ( NW == north .last_vertex() );
	size_t N_horiz = south .number_of ( tag::segments );
	assert ( N_horiz == north .number_of ( tag::segments ) );
	size_t N_vert = east .number_of ( tag::segments );
	assert ( N_vert == west .number_of ( tag::segments ) );

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
		double frac_N = double(i) / double(N_vert),  alpha = frac_N * (1.-frac_N);
		alpha = alpha*alpha*alpha;
		it_south.reset();  it_south++;
		it_north.reset();  it_north++;
		for ( size_t j = 1; j < N_horiz; j++ )
		{	Cell AB = *it;  // 'it' points into the 'horizon' list of segments
			Cell B = AB.tip();
			Cell C ( tag::vertex );  // create a new vertex
			Cell ver_south = *it_south;
			Cell ver_north = *it_north;
			double frac_E = double(j) / double(N_horiz),  beta = frac_E * (1.-frac_E);
			beta = beta*beta*beta;
			double sum = alpha + beta, aa = alpha/sum, bb = beta/sum;
			space.interpolate ( C, bb*(1.-frac_N), ver_south, aa*frac_E,     ver_east,     
		                         bb*frac_N,     ver_north, aa*(1.-frac_E), ver_west );
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

	#ifndef NDEBUG
	{ // just a block of code for hiding 'spin'
	Manifold::Action spin;
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = south .iterator ( tag::over_cells, tag::of_max_dim );
	for ( it .reset(); it .in_range(); it++ )  spin += (*it) .winding();
	} { // just a block of code for hiding 'it'
	Mesh::Iterator it = east .iterator ( tag::over_cells, tag::of_max_dim );
	for ( it .reset(); it .in_range(); it++ )  spin += (*it) .winding();
	} { // just a block of code for hiding 'it'
	Mesh::Iterator it = north .iterator ( tag::over_cells, tag::of_max_dim );
	for ( it .reset(); it .in_range(); it++ )  spin += (*it) .winding();
	} { // just a block of code for hiding 'it'
	Mesh::Iterator it = west .iterator ( tag::over_cells, tag::of_max_dim );
	for ( it .reset(); it .in_range(); it++ )  spin += (*it) .winding();
	} // just a block of code for hiding 'it'
	assert ( spin == 0 );
	} // just a block of code for hiding 'spin'
	#endif
		
	// recover corners from the sides
	// the process is different from the one in 'build' without tag::winding
	// here, sides may be closed loops and thus methods 'first_vertex' and 'last_vertex'
	// become meaningless
	Cell SW = west  .common_vertex ( tag::with, south );
	Cell SE = south .common_vertex ( tag::with, east );
	Cell NE = east  .common_vertex ( tag::with, north );
	Cell NW = north .common_vertex ( tag::with, west );
	
	size_t N_horiz = south .number_of ( tag::segments );
	assert ( N_horiz == north .number_of ( tag::segments ) );
	size_t N_vert = east .number_of ( tag::segments );
	assert ( N_vert == west .number_of ( tag::segments ) );

	// prepare horizon
	std::list <Cell> horizon;
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = south .iterator ( tag::over_segments, tag::require_order );
	for ( it.reset(); it .in_range(); it++ )
	{	Cell seg = *it;  horizon .push_back ( seg );  }
	} // just a block of code for hiding 'it'

	// we must deal with possible winding segments
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

Mesh Mesh::common_edge ( const tag::With &, const Mesh & m2 ) const

// 'this' and 'm2' are two-dimensional meshes
// the common edge will be built as Connected::OneDim
// in the future, we should introduce a general method 'intersection'

// if 'common_edge' is on the boundary of 'this',
// its orientation will be compatible with 'this' -- thus, opposite to 'm2'

// however, it may happen that the intersection be inner to either of 'this' and 'm2'
// in this case, the orientation is not clearly defined
// this is especially relevant if the working manifold is a quotient manifold
// (if either of 'this' and 'm2' have winding)
	
{	assert ( this->dim() == 2 );
	assert ( m2 .dim() == 2 );

	Mesh result ( tag::fuzzy, tag::of_dim, 1 );
	
	bool on_bdry_of_this = false, on_bdry_of_m2 = false;
	Mesh::Iterator it = this->iterator ( tag::over_segments );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell seg = *it;
		if ( not seg .belongs_to ( m2 ) )  continue;
		if ( ( this->cell_in_front_of ( seg, tag::may_not_exist ) .exists() and
		       not this->cell_behind ( seg, tag::may_not_exist ) .exists()      ) or
		     ( not this->cell_in_front_of ( seg, tag::may_not_exist ) .exists() and
		       this->cell_behind ( seg, tag::may_not_exist ) .exists()              ) )
			on_bdry_of_this = true;
		if ( ( m2 .cell_in_front_of ( seg, tag::may_not_exist ) .exists() and
		       not m2 .cell_behind ( seg, tag::may_not_exist ) .exists()      ) or
		     ( not m2 .cell_in_front_of ( seg, tag::may_not_exist ) .exists() and
		       m2 .cell_behind ( seg, tag::may_not_exist ) .exists()              ) )
			on_bdry_of_m2 = true;
		break;                                                                         }

	if ( not ( on_bdry_of_this or on_bdry_of_m2 ) )
	{	std::cout << "common edge is inner to both meshes" << std::endl;
		std::cout << "maniFEM cannot deal with such a configuration (yet)" << std::endl;
		exit (1);                                                                        }

	if ( on_bdry_of_this )   // 'common_edge' will be oriented accordingly to 'this'
	for ( it .reset(); it .in_range(); it++ )
	{	Cell seg = *it;
		if ( not seg .belongs_to ( m2 ) )  continue;
		if ( this->cell_in_front_of ( seg, tag::may_not_exist ) .exists() )
		{	assert ( not this->cell_behind ( seg, tag::may_not_exist ) .exists() );
			seg .reverse() .add_to_mesh ( result );                                   }
		else
		{	assert ( this->cell_behind ( seg, tag::may_not_exist ) .exists() );
			seg .add_to_mesh ( result );                                          }  }
	else   // 'common_edge' will be oriented opposite to 'm2'
	for ( it .reset(); it .in_range(); it++ )
	{	Cell seg = *it;
		if ( not seg .belongs_to ( m2 ) )  continue;
		if ( m2 .cell_in_front_of ( seg, tag::may_not_exist ) .exists() )
		{	assert ( not m2 .cell_behind ( seg, tag::may_not_exist ) .exists() );
			seg .add_to_mesh ( result );                                         }
		else
		{	assert ( m2 .cell_behind ( seg, tag::may_not_exist ) .exists() );
			seg .reverse() .add_to_mesh ( result );                           }    }
	
	return result .convert_to ( tag::connected, tag::one_dim, tag::surely_exists );
	
}  // end of  Mesh::common_edge

//---------------------------------------------------------------------------------------------//

	
void Mesh::build ( const tag::Hexahedron &, const Mesh & south, const Mesh & north,
	                 Mesh east, Mesh west, const Mesh & up, const Mesh & down        )

// see paragraph 12.5 in the manual
	
// the cube may be in any position, faces and edges may be curved
// we can always rotate it and consider that 'up' is up, 'down' is down
// 'south' is close to the viewer, 'north' is far from the viewer

// even so, we are left with two possibilities :
// either 'east' is towards our left and 'west' is towards our right
// or the other way around

{	Manifold space = Manifold::working;
	assert ( space .exists() );  // we use the current manifold

	Function xyz = space .coordinates();
	Function x = xyz [0], y = xyz [1], z = xyz [2];

	Mesh up_south = up .common_edge ( tag::with, south );  // orientation compatible with 'up'
	Mesh up_north = up .common_edge ( tag::with, north );  // orientation compatible with 'up'
	Mesh up_east = up .common_edge ( tag::with, east );  // orientation compatible with 'up'
	Mesh up_west = up .common_edge ( tag::with, west );  // orientation compatible with 'up'
	Mesh down_south = down .common_edge ( tag::with, south );
		// orientation compatible with 'down'
	Mesh down_north = down .common_edge ( tag::with, north );
		// orientation compatible with 'down'
	Mesh down_east = down .common_edge ( tag::with, east );  // orientation compatible with 'down'
	Mesh down_west = down .common_edge ( tag::with, west );  // orientation compatible with 'down'
	Mesh south_east = south .common_edge ( tag::with, east );
		// orientation compatible with 'south'
	Mesh south_west = south .common_edge ( tag::with, west );
		// orientation compatible with 'south'
	Mesh north_east = north .common_edge ( tag::with, east );
		// orientation compatible with 'north'
	Mesh north_west = north .common_edge ( tag::with, west );
		// orientation compatible with 'north'

	// we let our intuition be guided by drawing in paragraph 12.5 in the manual
	// we are looking at south and assume its four edges are oriented counter-clockwise
	// (the normal points towards us, towards the exterior of the cube)
	if ( down_south .last_vertex() == down_east .first_vertex() .reverse() )
		// 'east' is on our left, we must switch them
	{	Mesh tmp = east;  east = west;  west = tmp;
		tmp = up_east;  up_east = up_west;  up_west = tmp;
		tmp = down_east;  down_east = down_west;  down_west = tmp;
		tmp = south_east;  south_east = south_west;  south_west = tmp;
		tmp = north_east;  north_east = north_west;  north_west = tmp;  }

	Cell ver_down_south_west = down_south .last_vertex();
	assert ( ver_down_south_west == down_west .first_vertex() .reverse() );
	assert ( ver_down_south_west == south_west .last_vertex() );
	Cell ver_down_south_east = down_south .first_vertex() .reverse();
	assert ( ver_down_south_east == down_east .last_vertex() );
	assert ( ver_down_south_east == south_east .first_vertex() .reverse() );
	Cell ver_up_south_east = up_south .last_vertex();
	assert ( ver_up_south_east == up_east .first_vertex() .reverse() );
	assert ( ver_up_south_east == south_east .last_vertex() );
	Cell ver_up_south_west = up_south .first_vertex() .reverse();
	assert ( ver_up_south_west == up_west .last_vertex() );
	assert ( ver_up_south_west == south_west .first_vertex() .reverse() );
	Cell ver_down_north_east = down_north .last_vertex();
	assert ( ver_down_north_east == down_east .first_vertex() .reverse() );
	assert ( ver_down_north_east == north_east .last_vertex() );
	Cell ver_down_north_west = down_north .first_vertex() .reverse();
	assert ( ver_down_north_west == down_west .last_vertex() );
	assert ( ver_down_north_west == north_west .first_vertex() .reverse() );
	Cell ver_up_north_east = up_north .first_vertex() .reverse();
	assert ( ver_up_north_east == up_east .last_vertex() );
	assert ( ver_up_north_east == north_east .first_vertex() .reverse() );

	size_t nb_SN = up_west .number_of ( tag::segments );
	assert ( nb_SN == up_east .number_of ( tag::segments ) );
	assert ( nb_SN == down_west .number_of ( tag::segments ) );
	assert ( nb_SN == down_east .number_of ( tag::segments ) );
	size_t nb_EW = up_south .number_of ( tag::segments );
	assert ( nb_EW == up_north .number_of ( tag::segments ) );
	assert ( nb_EW == down_south .number_of ( tag::segments ) );
	assert ( nb_EW == down_north .number_of ( tag::segments ) );
	size_t nb_ud = south_east .number_of ( tag::segments );
	assert ( nb_ud == south_west .number_of ( tag::segments ) );
	assert ( nb_ud == north_east .number_of ( tag::segments ) );
	assert ( nb_ud == north_west .number_of ( tag::segments ) );

	// we make a copy of half of the boundary and shrink it gradually,
	// until we touch the other half
	Mesh front ( tag::join, down, south, west );

	// we use six vertices, one on each side of the cube, for interpolating coordinates
	// to help us move in the right direction, we keep segments rather than vertices
	assert ( ver_down_south_west == down_west .first_vertex() .reverse() );
	Cell seg_south_west = down_west .cell_in_front_of ( ver_down_south_west, tag::surely_exists );
  Cell seg_north_west = down_west .cell_behind ( ver_down_north_west, tag::surely_exists );
  Cell seg_south_east = down_east .cell_behind ( ver_down_south_east, tag::surely_exists )
		.reverse ( tag::surely_exists );
	// the three segments above are horizontal and point towards north
	
	for ( size_t i_up = 1; i_up < nb_ud; i_up ++ )   // counts jumps upwards
	{	// we move 'seg_south_west' up one square, keeping its base within 'south_west'
		Mesh sq_bdry = west .cell_in_front_of ( seg_south_west, tag::surely_exists ) .boundary();
		seg_south_west = sq_bdry .cell_behind ( seg_south_west .tip(), tag::surely_exists );
		// 'seg_south_west' is now vertical, points down
		seg_south_west =
			sq_bdry .cell_behind ( seg_south_west .base() .reverse(), tag::surely_exists );
		// 'seg_south_west' is horizontal again, points towards north
		// we move 'seg_south_east' up one square, keeping its base within 'south_east'
		sq_bdry = east .cell_behind ( seg_south_east, tag::surely_exists ) .boundary();
		seg_south_east = sq_bdry .cell_in_front_of ( seg_south_east .tip(), tag::surely_exists );
		// 'seg_south_east' is now vertical, points up
		seg_south_east = sq_bdry .cell_in_front_of ( seg_south_east .tip(), tag::surely_exists )
		                 .reverse ( tag::surely_exists );
		// 'seg_south_east' is horizontal again, points towards north
		// we move 'seg_north_west' up one square, keeping its tip within 'north_west'
		sq_bdry = west .cell_in_front_of ( seg_north_west, tag::surely_exists ) .boundary();
		seg_north_west = sq_bdry .cell_behind ( seg_north_west .tip(), tag::surely_exists );
		// 'seg_north_west' is now vertical, points down
		seg_north_west =
			sq_bdry .cell_behind ( seg_north_west .base() .reverse(), tag::surely_exists );
		// 'seg_north_west' is horizontal again, points towards north
		
		Cell seg_west = seg_south_west;
		Cell seg_east = seg_south_east;
		Cell seg_down_west = down_west .cell_in_front_of ( ver_down_south_west, tag::surely_exists );
		Cell seg_up_west = up_west .cell_behind ( ver_up_south_west, tag::surely_exists )
		                   .reverse ( tag::surely_exists );
		// the two segments above are horizontal and point towards north
		double frac_up = double(i_up) / double(nb_ud), alpha = frac_up * (1.-frac_up);
		alpha = alpha * alpha * alpha;

		for ( size_t i_north = 1; i_north < nb_SN; i_north ++ )    // counts jumps toward north
		{	Cell ver_west = seg_west .tip();
			Cell ver_east = seg_east .tip();
			Cell seg_up    = seg_up_west;
			Cell seg_down  = seg_down_west;
			// the two segments above are horizontal and point towards north
			Cell seg_south = south_west .cell_behind  // points down
				( seg_south_west .base() .reverse(), tag::surely_exists );
			sq_bdry = south .cell_behind ( seg_south, tag::surely_exists ) .boundary();
			seg_south = sq_bdry .cell_in_front_of ( seg_south .tip(), tag::surely_exists );
			Cell seg_north = north_west .cell_in_front_of  // points up
				( seg_north_west .tip(), tag::surely_exists );
			sq_bdry = north .cell_behind ( seg_north, tag::surely_exists ) .boundary();
			seg_north = sq_bdry .cell_behind
				( seg_north_west .tip(), tag::surely_exists ) .reverse ( tag::surely_exists );
			// 'seg_south' 'seg_north' point towards east (right)

			// define 'seg1' and 'seg2'
			Cell seg2 = seg_west;  // horizontal, points towards north
			sq_bdry = west .cell_behind ( seg_west, tag::surely_exists ) .boundary();
			Cell seg1 = sq_bdry .cell_in_front_of ( seg2 .tip(), tag::surely_exists );
			// seg1 is vertical, points down
			seg1 = sq_bdry .cell_in_front_of ( seg1 .tip(), tag::surely_exists )
				.reverse ( tag::surely_exists );
			// seg1 is horizontal now, points towards north, it lies one square below seg2
			// these two segments will move together towards east, keeping this relative position
			double frac_N = double(i_north) / double(nb_SN),  beta = frac_N * (1.-frac_N);
			beta = beta * beta * beta;

			for ( size_t i_east = 1; i_east < nb_EW; i_east ++ )    // counts jumps toward east
			{	// move 'seg_up' and 'seg_down' towards east (right)
				sq_bdry = up .cell_in_front_of ( seg_up, tag::surely_exists ) .boundary();
				seg_up  = sq_bdry .cell_behind ( seg_up .tip(), tag::surely_exists );
				// now 'seg_up' points towards west (left)
				seg_up  = sq_bdry .cell_behind ( seg_up .base() .reverse(), tag::surely_exists );
				// 'seg_up' points again towards north
				sq_bdry = down .cell_behind ( seg_down, tag::surely_exists ) .boundary();
				seg_down  = sq_bdry .cell_in_front_of ( seg_down .tip(), tag::surely_exists );
				// now 'seg_down' points towards east (right)
				seg_down  = sq_bdry .cell_in_front_of
					( seg_down .tip(), tag::surely_exists ) .reverse ( tag::surely_exists );
				// 'seg_down' points again towards north
				Cell ver_down  = seg_down  .tip();
				Cell ver_up    = seg_up    .tip();
				Cell ver_south = seg_south .tip();
				Cell ver_north = seg_north .tip();

				// we want to build a new cube
				// three faces exist already
				// we need to build three new faces
				// before that, three new segments
				// and before that, a new vertex
				// looks like a poem, does it not ?
				
				double frac_E = double(i_east) / double(nb_EW),  gamma = frac_E * (1.-frac_E);
				gamma = gamma * gamma * gamma;
				double sum = std::sqrt ( alpha*beta + beta*gamma + alpha*gamma ),
				       aa = alpha / sum, bb = beta / sum, cc = gamma / sum;
				Cell new_ver ( tag::vertex );
				space.interpolate ( new_ver,
	bb*cc*frac_up,      ver_up,   aa*cc*frac_N,      ver_north, aa*bb*frac_E,      ver_east,
	bb*cc*(1.-frac_up), ver_down, aa*cc*(1.-frac_N), ver_south, aa*bb*(1.-frac_E), ver_west );

				// recall that 'seg1' and 'seg2' are two horizontal parallel segments
				// pointing towards north, 'seg1' lying one square below 'seg2'
				Cell face_west = front .cell_behind ( seg2, tag::surely_exists );
				Cell face_down = front .cell_behind ( seg1, tag::surely_exists );
				Cell corner_down_NW = seg1 .tip();
				Cell corner_up_NW = seg2 .tip();
				Cell edge_ud_NW =
					face_west .boundary() .cell_behind ( corner_down_NW, tag::surely_exists );
				Cell edge_WE_down_N =
					face_down .boundary() .cell_in_front_of ( corner_down_NW, tag::surely_exists );
				Cell corner_down_NE = edge_WE_down_N .tip();
				Cell corner_down_SW = seg1 .base() .reverse();
				Cell edge_EW_down_S =
					face_down .boundary() .cell_behind ( corner_down_SW, tag::surely_exists );
				Cell face_south = front .cell_in_front_of ( edge_EW_down_S, tag::surely_exists );
				Cell corner_up_SW = seg2 .base() .reverse();
				Cell edge_EW_up_S =
					face_south .boundary() .cell_behind ( corner_up_SW, tag::surely_exists );
				Cell corner_up_SE = edge_EW_up_S .base() .reverse();
				Cell edge_du_SE = face_south .boundary() .cell_behind ( corner_up_SE );
				Cell corner_down_SE = edge_du_SE .base() .reverse();
				Cell edge_NS_down_E = face_down .boundary()
					.cell_behind ( corner_down_SE, tag::surely_exists );
				Cell edge_WE_up_N ( tag::segment, corner_up_NW .reverse(), new_ver );
				Cell edge_ud_NE ( tag::segment, new_ver .reverse(), corner_down_NE );
				Cell edge_NS_up_E ( tag::segment, new_ver .reverse(), corner_up_SE );
				Cell face_up ( tag::square, seg2 .reverse(),         edge_EW_up_S .reverse(),
				                            edge_NS_up_E .reverse(), edge_WE_up_N .reverse() );
				Cell face_north ( tag::square, edge_ud_NE,            edge_WE_down_N .reverse(),
				                               edge_ud_NW .reverse(), edge_WE_up_N              );
				Cell face_east  ( tag::square, edge_NS_up_E,              edge_du_SE .reverse(),
				                               edge_NS_down_E .reverse(), edge_ud_NE .reverse() );
				Cell new_cube ( tag::cube, face_up, face_down,
				                           face_south, face_north, face_east, face_west );
				new_cube .add_to_mesh ( *this );

				// propagate 'front'
				face_west .remove_from_mesh ( front );
				face_south .remove_from_mesh ( front );
				face_down .remove_from_mesh ( front );
				face_east .reverse() .add_to_mesh ( front );
				face_north .reverse() .add_to_mesh ( front );
				face_up .reverse() .add_to_mesh ( front );

				// move 'seg1' and 'seg2' towards east (right)
				seg1 = edge_NS_down_E .reverse();
				seg2 = edge_NS_up_E .reverse();
	
				// move 'seg_south' and 'seg_north' towards east (right)
				sq_bdry = south .cell_behind ( seg_south, tag::surely_exists ) .boundary();
				// 'sq' is above 'seg_south'
				seg_south  = sq_bdry .cell_in_front_of ( seg_south .tip(), tag::surely_exists );
				// now 'seg_south' points up
				sq_bdry = south .cell_in_front_of ( seg_south, tag::surely_exists ) .boundary();
				// 'sq' is eastern (to the right) to 'seg_south'
				seg_south  =
					sq_bdry .cell_in_front_of ( seg_south .base() .reverse(), tag::surely_exists );
				// 'seg_south' points again towards east (right)
				sq_bdry = north .cell_in_front_of ( seg_north, tag::surely_exists ) .boundary();
				// 'sq' is above 'seg_north'
				seg_north  = sq_bdry .cell_behind ( seg_north .tip(), tag::surely_exists );
				// now 'seg_north' points down
				sq_bdry = north .cell_in_front_of ( seg_north, tag::surely_exists ) .boundary();
				// 'sq' is eastern (to the right) to 'seg_north'
				seg_north  = sq_bdry .cell_behind
					( seg_north .tip(), tag::surely_exists ) .reverse ( tag::surely_exists );
				// 'seg_north' points again towards east (right)
			}  // end of movement  west -> east

			// build the last (eastern) cube of this row

			// recall that 'seg1' and 'seg2' are two horizontal parallel segments
			// pointing towards north, 'seg1' lying one square below 'seg2'
			Cell face_west = front .cell_behind ( seg2, tag::surely_exists );
			Cell face_down = front .cell_behind ( seg1, tag::surely_exists );
			Cell corner_down_NW = seg1 .tip();
			Cell corner_up_NW = seg2 .tip();
			Cell edge_ud_NW =
				face_west .boundary() .cell_behind ( corner_down_NW, tag::surely_exists );
			Cell edge_WE_down_N =
				face_down .boundary() .cell_in_front_of ( corner_down_NW, tag::surely_exists );
			Cell corner_down_NE = edge_WE_down_N .tip();
			Cell corner_down_SW = seg1 .base() .reverse();
			Cell edge_EW_down_S =
				face_down .boundary() .cell_behind ( corner_down_SW, tag::surely_exists );
			Cell face_south = front .cell_in_front_of ( edge_EW_down_S, tag::surely_exists );
			Cell corner_up_SW = seg2 .base() .reverse();
			Cell edge_EW_up_S =
				face_south .boundary() .cell_behind ( corner_up_SW, tag::surely_exists );
			Cell corner_up_SE = edge_EW_up_S .base() .reverse();
			Cell edge_du_SE = face_south .boundary() .cell_behind ( corner_up_SE );
			Cell corner_down_SE = edge_du_SE .base() .reverse();
			Cell edge_NS_down_E = face_down .boundary()
				.cell_behind ( corner_down_SE, tag::surely_exists );
			Cell face_east = east .cell_in_front_of ( edge_NS_down_E, tag::surely_exists );
			Cell edge_du_NE = face_east .boundary()
				.cell_in_front_of ( edge_NS_down_E .base() .reverse(), tag::surely_exists );
			Cell new_ver = edge_du_NE .tip();  // well, not really new ...
			Cell edge_NS_up_E =
				face_east .boundary() .cell_in_front_of ( new_ver, tag::surely_exists );
			Cell edge_WE_up_N ( tag::segment, corner_up_NW .reverse(), new_ver );
			Cell face_up ( tag::square, seg2 .reverse(),         edge_EW_up_S .reverse(),
			                            edge_NS_up_E .reverse(), edge_WE_up_N .reverse() );
			Cell face_north ( tag::square, edge_du_NE .reverse(), edge_WE_down_N .reverse(),
		                                 edge_ud_NW .reverse(), edge_WE_up_N              );
			Cell new_cube ( tag::cube, face_up, face_down,
											face_south, face_north, face_east, face_west );
			new_cube .add_to_mesh ( *this );

			// propagate 'front'
			face_west .remove_from_mesh ( front );
			face_south .remove_from_mesh ( front );
			face_down .remove_from_mesh ( front );
			face_north .reverse() .add_to_mesh ( front );
			face_up .reverse() .add_to_mesh ( front );

			// advance 'seg_west' and 'seg_east' towards north
			sq_bdry = west .cell_in_front_of ( seg_west, tag::surely_exists ) .boundary();
			seg_west  = sq_bdry .cell_behind ( seg_west .tip(), tag::surely_exists );
			// now 'seg_west' points down
			sq_bdry = west .cell_in_front_of ( seg_west, tag::surely_exists ) .boundary();
			seg_west  = sq_bdry .cell_behind ( seg_west .tip(), tag::surely_exists ) .reverse();
			// 'seg_west' points again towards north
			sq_bdry = east .cell_behind ( seg_east, tag::surely_exists ) .boundary();
			seg_east  = sq_bdry .cell_in_front_of ( seg_east .tip(), tag::surely_exists );
			// now 'seg_east' points up
			sq_bdry = east .cell_in_front_of ( seg_east, tag::surely_exists ) .boundary();
			seg_east  = sq_bdry .cell_in_front_of ( seg_east .base() .reverse(), tag::surely_exists );
			// 'seg_east' points again towards north
				
			// move 'seg_down_west' and 'seg_up_west' towards north
			seg_up_west = up_west .cell_behind ( seg_up_west .tip(), tag::surely_exists ) .reverse();
			seg_down_west = down_west .cell_in_front_of ( seg_down_west .tip(), tag::surely_exists );
		}  // end of movement  south -> north

		// build the last (northern) row of this layer

		// define 'seg1' and 'seg2'
		Cell seg2 = seg_west;  // horizontal, points towards north
		sq_bdry = west .cell_behind ( seg_west, tag::surely_exists ) .boundary();
		Cell seg1 = sq_bdry .cell_in_front_of ( seg2 .tip(), tag::surely_exists );
		// seg1 is vertical, points down
		seg1 = sq_bdry .cell_in_front_of ( seg1 .tip(), tag::surely_exists )
			.reverse ( tag::surely_exists );
		// seg1 is horizontal now, points towards north, it lies one square below seg2
		// these two segments will move together towards east, keeping this relative position

		for ( size_t i_east = 1; i_east < nb_EW; i_east ++ )    // counts jumps toward east
		{	// we want to build a new cube
			// four faces exist already
			// we need to build two new faces
			// before that, one new segment (no new vertex)
			
			// recall that 'seg1' and 'seg2' are two horizontal parallel segments
			// pointing towards north, 'seg1' lying one square below 'seg2'
			Cell face_west = front .cell_behind ( seg2, tag::surely_exists );
			Cell face_down = front .cell_behind ( seg1, tag::surely_exists );
			Cell corner_down_NW = seg1 .tip();
			Cell corner_up_NW = seg2 .tip();
			Cell edge_ud_NW =
				face_west .boundary() .cell_behind ( corner_down_NW, tag::surely_exists );
			Cell edge_WE_down_N =
				face_down .boundary() .cell_in_front_of ( corner_down_NW, tag::surely_exists );
			Cell corner_down_NE = edge_WE_down_N .tip();
			Cell corner_down_SW = seg1 .base() .reverse();
			Cell edge_EW_down_S =
				face_down .boundary() .cell_behind ( corner_down_SW, tag::surely_exists );
			Cell face_south = front .cell_in_front_of ( edge_EW_down_S, tag::surely_exists );
			Cell corner_up_SW = seg2 .base() .reverse();
			Cell edge_EW_up_S =
				face_south .boundary() .cell_behind ( corner_up_SW, tag::surely_exists );
			Cell corner_up_SE = edge_EW_up_S .base() .reverse();
			Cell edge_du_SE = face_south .boundary() .cell_behind ( corner_up_SE );
			Cell corner_down_SE = edge_du_SE .base() .reverse();
			Cell edge_NS_down_E = face_down .boundary()
				.cell_behind ( corner_down_SE, tag::surely_exists );
			Cell face_north = north .cell_in_front_of ( edge_ud_NW, tag::surely_exists );
			Cell edge_WE_up_N = face_north .boundary()
				.cell_in_front_of ( edge_ud_NW .base() .reverse(), tag::surely_exists );
			Cell new_ver = edge_WE_up_N .tip();  // well, not really new ...
			Cell edge_ud_NE = face_north .boundary() .cell_in_front_of ( new_ver, tag::surely_exists );
			Cell edge_NS_up_E ( tag::segment, new_ver .reverse(), corner_up_SE );
			Cell face_up ( tag::square, seg2 .reverse(),         edge_EW_up_S .reverse(),
			                            edge_NS_up_E .reverse(), edge_WE_up_N .reverse() );
			Cell face_east  ( tag::square, edge_NS_up_E,              edge_du_SE .reverse(),
			                               edge_NS_down_E .reverse(), edge_ud_NE .reverse() );
			Cell new_cube ( tag::cube, face_up, face_down,
			                           face_south, face_north, face_east, face_west );
			new_cube .add_to_mesh ( *this );

			// propagate 'front'
			face_west .remove_from_mesh ( front );
			face_south .remove_from_mesh ( front );
			face_down .remove_from_mesh ( front );
			face_east .reverse() .add_to_mesh ( front );
			face_up .reverse() .add_to_mesh ( front );

			// move 'seg1' and 'seg2' towards east (right)
			seg1 = edge_NS_down_E .reverse();
			seg2 = edge_NS_up_E .reverse();

		}  // end of movement  west -> east

		// build the last (eastern) cube of this last row of the current layer

		// recall that 'seg1' and 'seg2' are two horizontal parallel segments
		// pointing towards north, 'seg1' lying one square below 'seg2'
		Cell face_west = front .cell_behind ( seg2, tag::surely_exists );
		Cell face_down = front .cell_behind ( seg1, tag::surely_exists );
		Cell corner_down_NW = seg1 .tip();
		Cell corner_up_NW = seg2 .tip();
		Cell edge_ud_NW = face_west .boundary() .cell_behind ( corner_down_NW, tag::surely_exists );
		Cell edge_WE_down_N =
			face_down .boundary() .cell_in_front_of ( corner_down_NW, tag::surely_exists );
		Cell corner_down_NE = edge_WE_down_N .tip();
		Cell corner_down_SW = seg1 .base() .reverse();
		Cell edge_EW_down_S =
			face_down .boundary() .cell_behind ( corner_down_SW, tag::surely_exists );
		Cell face_south = front .cell_in_front_of ( edge_EW_down_S, tag::surely_exists );
		Cell corner_up_SW = seg2 .base() .reverse();
		Cell edge_EW_up_S = face_south .boundary() .cell_behind ( corner_up_SW, tag::surely_exists );
		Cell corner_up_SE = edge_EW_up_S .base() .reverse();
		Cell edge_du_SE = face_south .boundary() .cell_behind ( corner_up_SE );
		Cell corner_down_SE = edge_du_SE .base() .reverse();
		Cell edge_NS_down_E = face_down .boundary()
			.cell_behind ( corner_down_SE, tag::surely_exists );
		Cell face_east = east .cell_in_front_of ( edge_NS_down_E, tag::surely_exists );
		Cell edge_du_NE = face_east .boundary()
			.cell_in_front_of ( edge_NS_down_E .base() .reverse(), tag::surely_exists );
		Cell face_north = north .cell_in_front_of ( edge_ud_NW, tag::surely_exists );
		Cell edge_WE_up_N = face_north .boundary()
			.cell_in_front_of ( edge_ud_NW .base() .reverse(), tag::surely_exists );
		Cell new_ver = edge_WE_up_N .tip();  // well, not really new ...
		Cell edge_NS_up_E = face_east .boundary() .cell_in_front_of ( new_ver, tag::surely_exists );
		Cell face_up ( tag::square, seg2 .reverse(),         edge_EW_up_S .reverse(),
		                            edge_NS_up_E .reverse(), edge_WE_up_N .reverse() );
		Cell new_cube ( tag::cube, face_up, face_down,
										face_south, face_north, face_east, face_west );
		new_cube .add_to_mesh ( *this );

		// propagate 'front'
		face_west .remove_from_mesh ( front );
		face_south .remove_from_mesh ( front );
		face_down .remove_from_mesh ( front );
		face_up .reverse() .add_to_mesh ( front );

	}  // end of movement upwards

	// build the last (upper) layer

	Cell seg_west = up_west .cell_behind ( ver_up_south_west, tag::surely_exists )
		.reverse ( tag::surely_exists );

	for ( size_t i_north = 1; i_north < nb_SN; i_north ++ )    // counts jumps toward north
	{	// define 'seg1' and 'seg2'
		Cell seg2 = seg_west;  // horizontal, points towards north
		Mesh sq_bdry = west .cell_behind ( seg_west, tag::surely_exists ) .boundary();
		Cell seg1 = sq_bdry .cell_in_front_of ( seg2 .tip(), tag::surely_exists );
		// seg1 is vertical, points down
		seg1 = sq_bdry .cell_in_front_of ( seg1 .tip(), tag::surely_exists )
			.reverse ( tag::surely_exists );
		// seg1 is horizontal now, points towards north, it lies one square below seg2
		// these two segments will move together towards east, keeping this relative position

		for ( size_t i_east = 1; i_east < nb_EW; i_east ++ )    // counts jumps toward east
		{	// we want to build a new cube
			// four faces exist already
			// we need to build two new faces
			// one new segment, no new vertex
				
			// recall that 'seg1' and 'seg2' are two horizontal parallel segments
			// pointing towards north, 'seg1' lying one square below 'seg2'
			Cell face_west = front .cell_behind ( seg2, tag::surely_exists );
			Cell face_down = front .cell_behind ( seg1, tag::surely_exists );
			Cell corner_down_NW = seg1 .tip();
			Cell corner_up_NW = seg2 .tip();
			Cell edge_ud_NW =
				face_west .boundary() .cell_behind ( corner_down_NW, tag::surely_exists );
			Cell edge_WE_down_N =
				face_down .boundary() .cell_in_front_of ( corner_down_NW, tag::surely_exists );
			Cell corner_down_NE = edge_WE_down_N .tip();
			Cell corner_down_SW = seg1 .base() .reverse();
			Cell edge_EW_down_S =
				face_down .boundary() .cell_behind ( corner_down_SW, tag::surely_exists );
			Cell face_south = front .cell_in_front_of ( edge_EW_down_S, tag::surely_exists );
			Cell corner_up_SW = seg2 .base() .reverse();
			Cell edge_EW_up_S =
				face_south .boundary() .cell_behind ( corner_up_SW, tag::surely_exists );
			Cell corner_up_SE = edge_EW_up_S .base() .reverse();
			Cell edge_du_SE = face_south .boundary() .cell_behind ( corner_up_SE );
			Cell corner_down_SE = edge_du_SE .base() .reverse();
			Cell edge_NS_down_E = face_down .boundary()
				.cell_behind ( corner_down_SE, tag::surely_exists );
			Cell face_up = up .cell_in_front_of ( seg2, tag::surely_exists );
			Cell edge_WE_up_N = face_up .boundary() .
				cell_behind ( seg2 .tip(), tag::surely_exists ) .reverse ( tag::surely_exists );
			Cell new_ver = edge_WE_up_N .tip();  // well, not really new ...
			Cell edge_NS_up_E = face_up .boundary() .cell_behind
				( new_ver, tag::surely_exists ) .reverse ( tag::surely_exists );
			Cell edge_ud_NE ( tag::segment, new_ver .reverse(), corner_down_NE );
			Cell face_north ( tag::square, edge_ud_NE,            edge_WE_down_N .reverse(),
			                               edge_ud_NW .reverse(), edge_WE_up_N              );
			Cell face_east  ( tag::square, edge_NS_up_E,              edge_du_SE .reverse(),
			                               edge_NS_down_E .reverse(), edge_ud_NE .reverse() );
			Cell new_cube ( tag::cube, face_up, face_down,
			                           face_south, face_north, face_east, face_west );
			new_cube .add_to_mesh ( *this );

			// propagate 'front'
			face_west .remove_from_mesh ( front );
			face_south .remove_from_mesh ( front );
			face_down .remove_from_mesh ( front );
			face_east .reverse() .add_to_mesh ( front );
			face_north .reverse() .add_to_mesh ( front );

			// move 'seg1' and 'seg2' towards east (right)
			seg1 = edge_NS_down_E .reverse();
			seg2 = edge_NS_up_E .reverse();
	
		}  // end of movement  west -> east

		// build the last (eastern) cube of this row (within the last layer)

		// recall that 'seg1' and 'seg2' are two horizontal parallel segments
		// pointing towards north, 'seg1' lying one square below 'seg2'
		Cell face_west = front .cell_behind ( seg2, tag::surely_exists );
		Cell face_down = front .cell_behind ( seg1, tag::surely_exists );
		Cell corner_down_NW = seg1 .tip();
		Cell corner_up_NW = seg2 .tip();
		Cell edge_ud_NW = face_west .boundary() .cell_behind ( corner_down_NW, tag::surely_exists );
		Cell edge_WE_down_N =
			face_down .boundary() .cell_in_front_of ( corner_down_NW, tag::surely_exists );
		Cell corner_down_NE = edge_WE_down_N .tip();
		Cell corner_down_SW = seg1 .base() .reverse();
		Cell edge_EW_down_S =
			face_down .boundary() .cell_behind ( corner_down_SW, tag::surely_exists );
		Cell face_south = front .cell_in_front_of ( edge_EW_down_S, tag::surely_exists );
		Cell corner_up_SW = seg2 .base() .reverse();
		Cell edge_EW_up_S = face_south .boundary() .cell_behind ( corner_up_SW, tag::surely_exists );
		Cell corner_up_SE = edge_EW_up_S .base() .reverse();
		Cell edge_du_SE = face_south .boundary() .cell_behind ( corner_up_SE );
		Cell corner_down_SE = edge_du_SE .base() .reverse();
		Cell edge_NS_down_E = face_down .boundary()
			.cell_behind ( corner_down_SE, tag::surely_exists );
		Cell face_east = east .cell_in_front_of ( edge_NS_down_E, tag::surely_exists );
		Cell edge_du_NE = face_east .boundary()
			.cell_in_front_of ( edge_NS_down_E .base() .reverse(), tag::surely_exists );
		Cell new_ver = edge_du_NE .tip();  // well, not really new ...
		Cell edge_NS_up_E = face_east .boundary() .cell_in_front_of ( new_ver, tag::surely_exists );
		Cell face_up = up .cell_in_front_of ( seg2, tag::surely_exists );
		Cell edge_WE_up_N = face_up .boundary() .
			cell_behind ( seg2 .tip(), tag::surely_exists ) .reverse ( tag::surely_exists );
		assert ( new_ver == edge_WE_up_N .tip() );
		Cell face_north ( tag::square, edge_du_NE .reverse(), edge_WE_down_N .reverse(),
	                                 edge_ud_NW .reverse(), edge_WE_up_N              );
		Cell new_cube ( tag::cube, face_up, face_down,
										face_south, face_north, face_east, face_west );
		new_cube .add_to_mesh ( *this );

		// propagate 'front'
		face_west .remove_from_mesh ( front );
		face_south .remove_from_mesh ( front );
		face_down .remove_from_mesh ( front );
		face_north .reverse() .add_to_mesh ( front );

		seg_west = up_west .cell_behind ( seg_west .tip(), tag::surely_exists )
			.reverse ( tag::surely_exists );
	}  // end of movement  south -> north

	// build the last (northern) row of this last layer

	// define 'seg1' and 'seg2'
	Cell seg2 = seg_west;  // horizontal, points towards north
	Mesh sq_bdry = west .cell_behind ( seg_west, tag::surely_exists ) .boundary();
	Cell seg1 = sq_bdry .cell_in_front_of ( seg2 .tip(), tag::surely_exists );
	// seg1 is vertical, points down
	seg1 = sq_bdry .cell_in_front_of ( seg1 .tip(), tag::surely_exists )
		.reverse ( tag::surely_exists );
	// seg1 is horizontal now, points towards north, it lies one square below seg2
	// these two segments will move together towards east, keeping this relative position

	for ( size_t i_east = 1; i_east < nb_EW; i_east ++ )    // counts jumps toward east
	{	// we want to build a new cube
		// five faces exist already
		// we need to build one new face
		// no new segment, no new vertex
			
		// recall that 'seg1' and 'seg2' are two horizontal parallel segments
		// pointing towards north, 'seg1' lying one square below 'seg2'
		Cell face_west = front .cell_behind ( seg2, tag::surely_exists );
		Cell face_down = front .cell_behind ( seg1, tag::surely_exists );
		Cell corner_down_NW = seg1 .tip();
		Cell corner_up_NW = seg2 .tip();
		Cell edge_ud_NW = face_west .boundary() .cell_behind ( corner_down_NW, tag::surely_exists );
		Cell edge_WE_down_N =
			face_down .boundary() .cell_in_front_of ( corner_down_NW, tag::surely_exists );
		Cell corner_down_NE = edge_WE_down_N .tip();
		Cell corner_down_SW = seg1 .base() .reverse();
		Cell edge_EW_down_S =
			face_down .boundary() .cell_behind ( corner_down_SW, tag::surely_exists );
		Cell face_south = front .cell_in_front_of ( edge_EW_down_S, tag::surely_exists );
		Cell corner_up_SW = seg2 .base() .reverse();
		Cell edge_EW_up_S = face_south .boundary() .cell_behind ( corner_up_SW, tag::surely_exists );
		Cell corner_up_SE = edge_EW_up_S .base() .reverse();
		Cell edge_du_SE = face_south .boundary() .cell_behind ( corner_up_SE );
		Cell corner_down_SE = edge_du_SE .base() .reverse();
		Cell edge_NS_down_E = face_down .boundary()
			.cell_behind ( corner_down_SE, tag::surely_exists );
		Cell face_north = north .cell_in_front_of ( edge_ud_NW, tag::surely_exists );
		Cell edge_WE_up_N = face_north .boundary()
			.cell_in_front_of ( edge_ud_NW .base() .reverse(), tag::surely_exists );
		Cell new_ver = edge_WE_up_N .tip();  // well, not really new ...
		Cell edge_ud_NE = face_north .boundary() .cell_in_front_of ( new_ver, tag::surely_exists );
		Cell face_up = up .cell_in_front_of ( seg2, tag::surely_exists );
		assert ( edge_WE_up_N == face_up .boundary()
		    .cell_behind ( seg2 .tip(), tag::surely_exists ) .reverse ( tag::surely_exists ) );
		Cell edge_NS_up_E = face_up .boundary() .cell_behind
			( new_ver, tag::surely_exists ) .reverse ( tag::surely_exists );
		Cell face_east  ( tag::square, edge_NS_up_E,              edge_du_SE .reverse(),
		                               edge_NS_down_E .reverse(), edge_ud_NE .reverse() );
		Cell new_cube ( tag::cube, face_up, face_down,
		                           face_south, face_north, face_east, face_west );
		new_cube .add_to_mesh ( *this );

		// propagate 'front'
		face_west .remove_from_mesh ( front );
		face_south .remove_from_mesh ( front );
		face_down .remove_from_mesh ( front );
		face_east .reverse() .add_to_mesh ( front );

		// move 'seg1' and 'seg2' towards east (right)
		seg1 = edge_NS_down_E .reverse();
		seg2 = edge_NS_up_E .reverse();

	}  // end of movement  west -> east

	// build the last cube of this last row of the last layer

	// recall that 'seg1' and 'seg2' are two horizontal parallel segments
	// pointing towards north, 'seg1' lying one square below 'seg2'
	Cell face_west = front .cell_behind ( seg2, tag::surely_exists );
	Cell face_down = front .cell_behind ( seg1, tag::surely_exists );
	Cell corner_down_NW = seg1 .tip();
	Cell corner_up_NW = seg2 .tip();
	Cell edge_ud_NW = face_west .boundary() .cell_behind ( corner_down_NW, tag::surely_exists );
	Cell edge_WE_down_N =
		face_down .boundary() .cell_in_front_of ( corner_down_NW, tag::surely_exists );
	Cell corner_down_NE = edge_WE_down_N .tip();
	Cell corner_down_SW = seg1 .base() .reverse();
	Cell edge_EW_down_S =
		face_down .boundary() .cell_behind ( corner_down_SW, tag::surely_exists );
	Cell face_south = front .cell_in_front_of ( edge_EW_down_S, tag::surely_exists );
	Cell corner_up_SW = seg2 .base() .reverse();
	Cell edge_EW_up_S = face_south .boundary() .cell_behind ( corner_up_SW, tag::surely_exists );
	Cell corner_up_SE = edge_EW_up_S .base() .reverse();
	Cell edge_du_SE = face_south .boundary() .cell_behind ( corner_up_SE );
	Cell corner_down_SE = edge_du_SE .base() .reverse();
	Cell edge_NS_down_E = face_down .boundary()
		.cell_behind ( corner_down_SE, tag::surely_exists );
	Cell face_east = east .cell_in_front_of ( edge_NS_down_E, tag::surely_exists );
	Cell face_north = north .cell_in_front_of ( edge_ud_NW, tag::surely_exists );
	Cell face_up = up .cell_in_front_of ( seg2, tag::surely_exists );
	Cell new_cube ( tag::cube, face_up, face_down,
									face_south, face_north, face_east, face_west );
	new_cube .add_to_mesh ( *this );

} // end of Mesh::build with tag::hexahedron

//----------------------------------------------------------------------------------//

	
void Mesh::build ( const tag::Hexahedron &, const Mesh & south, const Mesh & north,
                   Mesh east, Mesh west, const Mesh & up, const Mesh & down,
                   const tag::Winding &                                            )

// the tag:::winding tells maniFEM that we are on a quotient manifold
// and that the faces provided (south, east, north, west, up, down) may have winding edges

// beware, south may be equal to north.reverse, east may be equal to west.reverse
// or they may be not equal but share the same vertices (and segments, reversed)
// beware, the correspondence may be not face-to-face,
// e.g. the first vertex of south_east may show up somewhere in the middle of north_east

// beware, edges may be closed loops, faces may be pieces of cylinder or 2D torus
	
// the cube may be in any position, faces may be curved
// we can always rotate it and consider that 'up' is up, 'down' is down
// 'south' is close to the viewer, 'north' is far from the viewer

// but we are left with two possibilities :
// either 'east' is towards our left and 'west' is towards our right
// or the other way around

{	Manifold space = Manifold::working;
	assert ( space .exists() );  // we use the current (quotient) manifold
	Manifold::Quotient * mani_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space .core );
	Function coords_q = space .coordinates();
	Manifold mani_Eu = mani_q->base_space;  // underlying Euclidian manifold
	Function coords_Eu = mani_Eu .coordinates();

	// Function xyz = space .coordinates();
	// Function x = xyz [0], y = xyz [1], z = xyz [2];

	#ifndef NDEBUG
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = south .iterator ( tag::over_cells, tag::of_max_dim );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell face = *it;
		Manifold::Action s;
		Mesh::Iterator itt = face .boundary() .iterator ( tag::over_cells, tag::of_max_dim );
		for ( itt .reset(); itt .in_range(); itt++ )  s += (*itt) .winding();
		assert ( s == 0 );                                                                    }
	} { // just a block of code for hiding 'it'
	Mesh::Iterator it = north .iterator ( tag::over_cells, tag::of_max_dim );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell face = *it;
		Manifold::Action s;
		Mesh::Iterator itt = face .boundary() .iterator ( tag::over_cells, tag::of_max_dim );
		for ( itt .reset(); itt .in_range(); itt++ )  s += (*itt) .winding();
		assert ( s == 0 );                                                                    }
	} { // just a block of code for hiding 'it'
	Mesh::Iterator it = east .iterator ( tag::over_cells, tag::of_max_dim );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell face = *it;
		Manifold::Action s;
		Mesh::Iterator itt = face .boundary() .iterator ( tag::over_cells, tag::of_max_dim );
		for ( itt .reset(); itt .in_range(); itt++ )  s += (*itt) .winding();
		assert ( s == 0 );                                                                    }
	} { // just a block of code for hiding 'it'
	Mesh::Iterator it = west .iterator ( tag::over_cells, tag::of_max_dim );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell face = *it;
		Manifold::Action s;
		Mesh::Iterator itt = face .boundary() .iterator ( tag::over_cells, tag::of_max_dim );
		for ( itt .reset(); itt .in_range(); itt++ )  s += (*itt) .winding();
		assert ( s == 0 );                                                                    }
	} { // just a block of code for hiding 'it'
	Mesh::Iterator it = up .iterator ( tag::over_cells, tag::of_max_dim );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell face = *it;
		Manifold::Action s;
		Mesh::Iterator itt = face .boundary() .iterator ( tag::over_cells, tag::of_max_dim );
		for ( itt .reset(); itt .in_range(); itt++ )  s += (*itt) .winding();
		assert ( s == 0 );                                                                    }
	} { // just a block of code for hiding 'it'
	Mesh::Iterator it = down .iterator ( tag::over_cells, tag::of_max_dim );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell face = *it;
		Manifold::Action s;
		Mesh::Iterator itt = face .boundary() .iterator ( tag::over_cells, tag::of_max_dim );
		for ( itt .reset(); itt .in_range(); itt++ )  s += (*itt) .winding();
		assert ( s == 0 );                                                                    }
	} // just a block of code for hiding 'it'
	#endif

	Mesh up_south = up .common_edge ( tag::with, south );  // orientation compatible with 'up'
	Mesh up_north = up .common_edge ( tag::with, north );  // orientation compatible with 'up'
	Mesh up_east  = up .common_edge ( tag::with, east );   // orientation compatible with 'up'
	Mesh up_west  = up .common_edge ( tag::with, west );   // orientation compatible with 'up'
	Mesh down_south = down .common_edge ( tag::with, south );
		// orientation compatible with 'down'
	Mesh down_north = down .common_edge ( tag::with, north );
		// orientation compatible with 'down'
	Mesh down_east = down .common_edge ( tag::with, east );  // orientation compatible with 'down'
	Mesh down_west = down .common_edge ( tag::with, west );  // orientation compatible with 'down'
	Mesh south_east = south .common_edge ( tag::with, east );
		// orientation compatible with 'south'
	Mesh south_west = south .common_edge ( tag::with, west );
		// orientation compatible with 'south'
	Mesh north_east = north .common_edge ( tag::with, east );
		// orientation compatible with 'north'
	Mesh north_west = north .common_edge ( tag::with, west );
		// orientation compatible with 'north'

	size_t nb_SN = up_west .number_of ( tag::segments );
	assert ( nb_SN == up_east .number_of ( tag::segments ) );
	assert ( nb_SN == down_west .number_of ( tag::segments ) );
	assert ( nb_SN == down_east .number_of ( tag::segments ) );
	size_t nb_EW = up_south .number_of ( tag::segments );
	assert ( nb_EW == up_north .number_of ( tag::segments ) );
	assert ( nb_EW == down_south .number_of ( tag::segments ) );
	assert ( nb_EW == down_north .number_of ( tag::segments ) );
	size_t nb_ud = south_east .number_of ( tag::segments );
	assert ( nb_ud == south_west .number_of ( tag::segments ) );
	assert ( nb_ud == north_east .number_of ( tag::segments ) );
	assert ( nb_ud == north_west .number_of ( tag::segments ) );

	// it may happen that : up_south   == up_north   .reverse()
	//                      down_south == down_north .reverse()
	//                      up_east    == up_west    .reverse()       and so on ...
	
	// recover corners from the edges
	// the process is different from the one in 'build' without tag::winding
	// here, edges may be closed loops and thus methods 'first_vertex' and 'last_vertex'
	// become meaningless

	// we let our intuition be guided by drawing in paragraph 12.5 in the manual
	// we are looking at south and assume its four edges are oriented counter-clockwise
	// (the normal points towards us, towards the exterior of the cube)

	Cell ver_down_south_west =
		down_south .common_vertex ( tag::with, down_west, tag::exactly_one );
	assert ( ver_down_south_west .belongs_to ( south_west ) );

	// true, 'down_south' is oriented accordingly to 'down', thus points left
	// but is 'west' on our left ?
	if ( not down_south .cell_behind ( ver_down_south_west, tag::may_not_exist ) .exists() or
	     not down_west .cell_in_front_of ( ver_down_south_west, tag::may_not_exist ) .exists() or
	     not south_west .cell_behind ( ver_down_south_west, tag::may_not_exist ) .exists()        )
		// looks like 'east' is on our left, we must switch them
	{	Mesh tmp = east;  east = west;  west = tmp;
		tmp = up_east;  up_east = up_west;  up_west = tmp;
		tmp = down_east;  down_east = down_west;  down_west = tmp;
		tmp = south_east;  south_east = south_west;  south_west = tmp;
		tmp = north_east;  north_east = north_west;  north_west = tmp;  }

	ver_down_south_west = down_south .common_vertex ( tag::with, down_west, tag::exactly_one );
	assert ( ver_down_south_west .belongs_to ( south_west ) );

	// in the case when all edges are closed loops, the above 'if' will never return true
	// (the cells exist always)
	// thus, west may remain on our right and east on our left
	// (opposite to the drawing in the manual)
	// there is nothing we can do about it;  most probably  east == west .reverse()
	// we will simply get a different orientation of the final mesh
		
	// since down_west is oriented accordingly to 'down', it points away from us, towards north :
	assert ( down_west .cell_in_front_of ( ver_down_south_west, tag::may_not_exist ) .exists() );
	// since south_west is oriented accordingly to 'south', it points down :
	assert ( south_west .cell_behind ( ver_down_south_west, tag::may_not_exist ) .exists() );

	Cell ver_down_south_east =
		down_south .common_vertex ( tag::with, down_east, tag::exactly_one );
	assert ( ver_down_south_east .belongs_to ( south_east ) );
	// since down_south is oriented accordingly to 'down', it points left, towards west
	assert ( down_south .cell_in_front_of ( ver_down_south_east, tag::may_not_exist ) .exists() );
	// since down_east is oriented accordingly to 'down', it points towards us, towards south
	assert ( down_east .cell_behind ( ver_down_south_east, tag::may_not_exist ) .exists() );
	// since south_east is oriented accordingly to 'south', it points up
	assert ( south_east .cell_in_front_of ( ver_down_south_east, tag::may_not_exist ) .exists() );
	
	// we must deal with possible winding segments
	// we choose that, at each interpolation operation, i.e., for each new vertex,
	// the new vertex will have winding zero relatively to ver_down_south_west
	// we must keep track of the windings of ver_south, ver_south_west and others
	// (relatively to ver_down_south_west)

	Manifold::Action spin_ver_down_south_east;
	{ // just a block for hiding 'it'
	Mesh::Iterator it = down_south .iterator ( tag::over_segments );
	for ( it .reset(); it .in_range(); it++ )
		spin_ver_down_south_east -= (*it) .winding();
	} // just a block for hiding 'it'
	
	Cell ver_up_south_west = south_west .common_vertex ( tag::with, up, tag::exactly_one );
	// since south_west is oriented accordingly to 'south', it points down
	assert ( south_west .cell_in_front_of ( ver_up_south_west, tag::may_not_exist ) .exists() );
	// since up_south is oriented accordingly to 'up', it points right, towards east
	assert ( up_south .cell_in_front_of ( ver_up_south_west, tag::may_not_exist ) .exists() );
	// since up_west is oriented accordingly to 'up', it points towards us, towards south
	assert ( up_west .cell_behind ( ver_up_south_west, tag::may_not_exist ) .exists() );

	Manifold::Action spin_ver_up_south_west;
	{ // just a block for hiding 'it'
	Mesh::Iterator it = south_west .iterator ( tag::over_segments );
	for ( it .reset(); it .in_range(); it++ )
		spin_ver_up_south_west -= (*it) .winding();
	} // just a block for hiding 'it'
	
	Cell ver_up_south_east = south_east .common_vertex ( tag::with, up, tag::exactly_one );
	// since south_east is oriented accordingly to 'south', it points up
	assert ( south_east .cell_behind ( ver_up_south_east, tag::may_not_exist ) .exists() );
	// since up_south is oriented accordingly to 'up', it points right, towards east
	assert ( up_south .cell_behind ( ver_up_south_east, tag::may_not_exist ) .exists() );
	// since up_east is oriented accordingly to 'up', it points away from us, towards north
	assert ( up_east .cell_in_front_of ( ver_up_south_east, tag::may_not_exist ) .exists() );

	Manifold::Action spin_tmp_1 = spin_ver_up_south_west;
	{ // just a block for hiding 'it'
	Mesh::Iterator it = up_south .iterator ( tag::over_segments );
	for ( it .reset(); it .in_range(); it++ )
		spin_tmp_1 += (*it) .winding();
	} // just a block for hiding 'it'
	#ifndef NDEBUG
	Manifold::Action spin_tmp_2 = spin_ver_down_south_east;
	{ // just a block for hiding 'it'
	Mesh::Iterator it = south_east .iterator ( tag::over_segments );
	for ( it .reset(); it .in_range(); it++ )
		spin_tmp_2 += (*it) .winding();
	} // just a block for hiding 'it'
	assert ( spin_tmp_1 == spin_tmp_2 );
	#endif
	Manifold::Action spin_ver_up_south_east = spin_tmp_1;
	
	Cell ver_down_north_west = north_west .common_vertex ( tag::with, down, tag::exactly_one );
	// since north_west is oriented accordingly to 'north', points up
	assert ( north_west .cell_in_front_of ( ver_down_north_west, tag::may_not_exist ) .exists() );
	// since down_west is oriented accordingly to 'down', points away from us, towards north
	assert ( down_west .cell_behind ( ver_down_north_west, tag::may_not_exist ) .exists() );
	// since down_north is oriented accordingly to 'down', points right, towards east
	assert ( down_north .cell_in_front_of ( ver_down_north_west, tag::may_not_exist ) .exists() );

	Manifold::Action spin_ver_down_north_west;
	{ // just a block for hiding 'it'
	Mesh::Iterator it = down_west .iterator ( tag::over_segments );
	for ( it .reset(); it .in_range(); it++ )
		spin_ver_down_north_west += (*it) .winding();
	} // just a block for hiding 'it'
	
	Cell ver_down_north_east = north_east .common_vertex ( tag::with, down, tag::exactly_one );
	// since north_east is oriented accordingly to 'north', points down
	assert ( north_east .cell_behind ( ver_down_north_east, tag::may_not_exist ) .exists() );
	// since down_east is oriented accordingly to 'down', points towards us, towards south
	assert ( down_east .cell_in_front_of ( ver_down_north_east, tag::may_not_exist ) .exists() );
	// since down_north is oriented accordingly to 'down', points right, towards east
	assert ( down_north .cell_behind ( ver_down_north_east, tag::may_not_exist ) .exists() );

	spin_tmp_1 = spin_ver_down_south_east;
	{ // just a block for hiding 'it'
	Mesh::Iterator it = down_east .iterator ( tag::over_segments );
	for ( it .reset(); it .in_range(); it++ )
		spin_tmp_1 -= (*it) .winding();
	} // just a block for hiding 'it'
	#ifndef NDEBUG
	spin_tmp_2 = spin_ver_down_north_west;
	{ // just a block for hiding 'it'
	Mesh::Iterator it = down_north .iterator ( tag::over_segments );
	for ( it .reset(); it .in_range(); it++ )
		spin_tmp_2 += (*it) .winding();
	} // just a block for hiding 'it'
	assert ( spin_tmp_1 == spin_tmp_2 );
	#endif
	Manifold::Action spin_ver_down_north_east = spin_tmp_1;
		
	// we make a copy of half of the boundary and move it gradually,
	// until we touch the other half
	// Mesh front ( tag::join, down, south, west );
	// unlike when we build a cube without winding,
	// here we cannot join these three faces
	// we use three fronts instead
	Mesh front_down  ( tag::fuzzy, tag::of_dim, 2 ),
	     front_south ( tag::fuzzy, tag::of_dim, 2 ),
       front_west  ( tag::fuzzy, tag::of_dim, 2 );
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = down .iterator ( tag::over_cells_of_max_dim );
	for ( it .reset(); it .in_range(); it++ )  (*it) .add_to_mesh ( front_down );
	} { // just a block of code for hiding 'it'
	Mesh::Iterator it = south .iterator ( tag::over_cells_of_max_dim );
	for ( it .reset(); it .in_range(); it++ )  (*it) .add_to_mesh ( front_south );
	} { // just a block of code for hiding 'it'
	Mesh::Iterator it = west .iterator ( tag::over_cells_of_max_dim );
	for ( it .reset(); it .in_range(); it++ )  (*it) .add_to_mesh ( front_west );
	} // just a block of code for hiding 'it'
	
	// we use six vertices, one on each side of the cube, for interpolating coordinates
	// to help us move in the right direction, we keep segments rather than vertices
	// we keep winding numbers of these vertices, tips of the respective segments
	Cell seg_south_west = down_west .cell_in_front_of ( ver_down_south_west, tag::surely_exists );
	Manifold::Action spin_south_west = seg_south_west .winding();
  Cell seg_north_west = down_west .cell_behind ( ver_down_north_west, tag::surely_exists );
	Manifold::Action spin_north_west = spin_ver_down_north_west;
  Cell seg_south_east = down_east .cell_behind ( ver_down_south_east, tag::surely_exists )
	                      .reverse ( tag::surely_exists );
	Manifold::Action spin_south_east = spin_ver_down_south_east + seg_south_east .winding();
	// the three segments above are horizontal and point towards north

	// we use six shadow vertices for interpolation
	Cell shadow_down  ( tag::vertex ), shadow_up    ( tag::vertex ),
	     shadow_south ( tag::vertex ), shadow_north ( tag::vertex ),
	     shadow_west  ( tag::vertex ), shadow_east  ( tag::vertex );

	for ( size_t i_up = 1; i_up < nb_ud; i_up ++ )   // counts jumps upwards
	{	// we move 'seg_south_west' up one square, keeping its base within 'south_west'
		Mesh sq_bdry = west .cell_in_front_of ( seg_south_west, tag::surely_exists ) .boundary();
		seg_south_west = sq_bdry .cell_behind ( seg_south_west .tip(), tag::surely_exists );
		// 'seg_south_west' is now vertical, points down
		spin_south_west -= seg_south_west .winding();
		seg_south_west =
			sq_bdry .cell_behind ( seg_south_west .base() .reverse(), tag::surely_exists );
		// 'seg_south_west' is horizontal again, points towards north
		// we move 'seg_south_east' up one square, keeping its base within 'south_east'
		sq_bdry = east .cell_behind ( seg_south_east, tag::surely_exists ) .boundary();
		seg_south_east = sq_bdry .cell_in_front_of ( seg_south_east .tip(), tag::surely_exists );
		// 'seg_south_east' is now vertical, points up
		spin_south_east += seg_south_east .winding();
		seg_south_east = sq_bdry .cell_in_front_of ( seg_south_east .tip(), tag::surely_exists )
		                 .reverse ( tag::surely_exists );
		// 'seg_south_east' is horizontal again, points towards north
		// we move 'seg_north_west' up one square, keeping its tip within 'north_west'
		sq_bdry = west .cell_in_front_of ( seg_north_west, tag::surely_exists ) .boundary();
		seg_north_west = sq_bdry .cell_behind ( seg_north_west .tip(), tag::surely_exists );
		// 'seg_north_west' is now vertical, points down
		spin_north_west -= seg_north_west .winding();
		seg_north_west =
			sq_bdry .cell_behind ( seg_north_west .base() .reverse(), tag::surely_exists );
		// 'seg_north_west' is horizontal again, points towards north
		
		Cell seg_west = seg_south_west;
		Manifold::Action spin_west = spin_south_west;
		Cell seg_east = seg_south_east;
		Manifold::Action spin_east = spin_south_east;
		Cell seg_down_west = down_west .cell_in_front_of ( ver_down_south_west, tag::surely_exists );
		Manifold::Action spin_down_west = seg_down_west .winding();
		Cell seg_up_west = up_west .cell_behind ( ver_up_south_west, tag::surely_exists )
		                   .reverse ( tag::surely_exists );
		Manifold::Action spin_up_west = spin_ver_up_south_west + seg_up_west .winding();
		// the two segments above are horizontal and point towards north
		double frac_up = double(i_up) / double(nb_ud), alpha = frac_up * (1.-frac_up);
		alpha = alpha * alpha * alpha;

		for ( size_t i_north = 1; i_north < nb_SN; i_north ++ )    // counts jumps toward north
		{	Cell ver_west = seg_west .tip();
			Cell ver_east = seg_east .tip();
			std::vector < double > v = coords_q ( ver_west, tag::winding, spin_west );
			coords_Eu ( shadow_west ) = v;
			v = coords_q ( ver_east, tag::winding, spin_east );
			coords_Eu ( shadow_east ) = v;
			Cell seg_up    = seg_up_west;
			Manifold::Action spin_up = spin_up_west;
			Cell seg_down  = seg_down_west;
			Manifold::Action spin_down = spin_down_west;
			// the two segments above are horizontal and point towards north
			Cell seg_south = south_west .cell_behind  // points down
				( seg_south_west .base() .reverse(), tag::surely_exists );
			sq_bdry = south .cell_behind ( seg_south, tag::surely_exists ) .boundary();
			seg_south = sq_bdry .cell_in_front_of ( seg_south .tip(), tag::surely_exists );
			Manifold::Action spin_south =
				spin_south_west - seg_south_west .winding() + seg_south .winding();
			Cell seg_north = north_west .cell_in_front_of  // points up
				( seg_north_west .tip(), tag::surely_exists );
			sq_bdry = north .cell_behind ( seg_north, tag::surely_exists ) .boundary();
			seg_north = sq_bdry .cell_behind
				( seg_north_west .tip(), tag::surely_exists ) .reverse ( tag::surely_exists );
			Manifold::Action spin_north = spin_north_west + seg_north .winding();
			// 'seg_south' 'seg_north' point towards east (right)

			// define 'seg1' and 'seg2'
			Cell seg2 = seg_west;  // horizontal, points towards north
			Manifold::Action spin2 = spin_west;
			sq_bdry = west .cell_behind ( seg_west, tag::surely_exists ) .boundary();
			Cell seg1 = sq_bdry .cell_in_front_of ( seg2 .tip(), tag::surely_exists );
			// seg1 is vertical, points down
			Manifold::Action spin1 = spin_west + seg1 .winding();;
			seg1 = sq_bdry .cell_in_front_of ( seg1 .tip(), tag::surely_exists )
			       .reverse ( tag::surely_exists );
			// seg1 is horizontal now, points towards north, it lies one square below seg2
			// these two segments will move together towards east, keeping this relative position
			double frac_N = double(i_north) / double(nb_SN), beta = frac_N * (1.-frac_N);
			beta = beta * beta * beta;

			for ( size_t i_east = 1; i_east < nb_EW; i_east ++ )    // counts jumps toward east
			{	// move 'seg_up' and 'seg_down' towards east (right)
				sq_bdry = up .cell_in_front_of ( seg_up, tag::surely_exists ) .boundary();
				seg_up  = sq_bdry .cell_behind ( seg_up .tip(), tag::surely_exists );
				// now 'seg_up' points towards west (left)
				spin_up -= seg_up .winding();
				seg_up  = sq_bdry .cell_behind ( seg_up .base() .reverse(), tag::surely_exists );
				// 'seg_up' points again towards north
				sq_bdry = down .cell_behind ( seg_down, tag::surely_exists ) .boundary();
				seg_down  = sq_bdry .cell_in_front_of ( seg_down .tip(), tag::surely_exists );
				// now 'seg_down' points towards east (right)
				spin_down += seg_down .winding();
				seg_down  = sq_bdry .cell_in_front_of ( seg_down .tip(), tag::surely_exists )
				            .reverse ( tag::surely_exists );
				// 'seg_down' points again towards north
				Cell ver_down  = seg_down  .tip();
				Cell ver_up    = seg_up    .tip();
				Cell ver_south = seg_south .tip();
				Cell ver_north = seg_north .tip();
				v = coords_q ( ver_down, tag::winding, spin_down );
				coords_Eu ( shadow_down ) = v;
				v = coords_q ( ver_up, tag::winding, spin_up );
				coords_Eu ( shadow_up ) = v;
				v = coords_q ( ver_south, tag::winding, spin_south );
				coords_Eu ( shadow_south ) = v;
				v = coords_q ( ver_north, tag::winding, spin_north );
				coords_Eu ( shadow_north ) = v;

				// we want to build a new cube
				// three faces exist already
				// we need to build three new faces
				// before that, three new segments
				// and before that, a new vertex
				// looks like a poem, does it not ?
				
				double frac_E = double(i_east) / double(nb_EW),  gamma = frac_E * (1.-frac_E);
				gamma = gamma * gamma * gamma;
				double sum = std::sqrt ( alpha*beta + beta*gamma + alpha*gamma ),
				       aa = alpha / sum, bb = beta / sum, cc = gamma / sum;
				Cell new_ver ( tag::vertex );
				// we chose that new_ver will have zero winding relatively to ver_down_south_west
				// but we have no control on the winding numbers of the six vertices
				// used for interpolating coordinates, so we use shadow vertices
				space.interpolate ( new_ver,
				                    bb*cc*frac_up,     shadow_up,    aa*cc*frac_N,       shadow_north,
				                    aa*bb*frac_E,      shadow_east,  bb*cc*(1.-frac_up), shadow_down,
				                    aa*cc*(1.-frac_N), shadow_south, aa*bb*(1.-frac_E),  shadow_west );

				// recall that 'seg1' and 'seg2' are two horizontal parallel segments
				// pointing towards north, 'seg1' lying one square below 'seg2'
				Cell face_west = front_west .cell_behind ( seg2, tag::surely_exists );
				Cell face_down = front_down .cell_behind ( seg1, tag::surely_exists );
				Cell corner_down_NW = seg1 .tip();
				Cell corner_up_NW = seg2 .tip();
				Cell edge_ud_NW =
					face_west .boundary() .cell_behind ( corner_down_NW, tag::surely_exists );
				Cell edge_WE_down_N =
					face_down .boundary() .cell_in_front_of ( corner_down_NW, tag::surely_exists );
				Cell corner_down_NE = edge_WE_down_N .tip();
				Cell corner_down_SW = seg1 .base() .reverse();
				Cell edge_EW_down_S =
					face_down .boundary() .cell_behind ( corner_down_SW, tag::surely_exists );
				Cell face_south = front_south .cell_in_front_of ( edge_EW_down_S, tag::surely_exists );
				Cell corner_up_SW = seg2 .base() .reverse();
				Cell edge_EW_up_S =
					face_south .boundary() .cell_behind ( corner_up_SW, tag::surely_exists );
				Cell corner_up_SE = edge_EW_up_S .base() .reverse();
				Cell edge_du_SE = face_south .boundary() .cell_behind ( corner_up_SE );
				Cell corner_down_SE = edge_du_SE .base() .reverse();
				Cell edge_NS_down_E =
					face_down .boundary() .cell_behind ( corner_down_SE, tag::surely_exists );
				Cell edge_WE_up_N ( tag::segment, corner_up_NW .reverse(), new_ver );
				edge_WE_up_N .winding() = -spin2;
				Cell edge_ud_NE ( tag::segment, new_ver .reverse(), corner_down_NE );
				edge_ud_NE .winding() = spin1 + edge_WE_down_N .winding();
				Cell edge_NS_up_E ( tag::segment, new_ver .reverse(), corner_up_SE );
				edge_NS_up_E .winding() = spin2 - seg2 .winding() - edge_EW_up_S .winding();
				//  seg2  can be seen as  edge_SN_up_W
				Cell face_up ( tag::square, seg2 .reverse(),         edge_EW_up_S .reverse(),
				                            edge_NS_up_E .reverse(), edge_WE_up_N .reverse() );
				Cell face_north ( tag::square, edge_ud_NE,            edge_WE_down_N .reverse(),
				                               edge_ud_NW .reverse(), edge_WE_up_N              );
				Cell face_east  ( tag::square, edge_NS_up_E,              edge_du_SE .reverse(),
				                               edge_NS_down_E .reverse(), edge_ud_NE .reverse() );
				Cell new_cube ( tag::cube, face_up, face_down,
				                           face_south, face_north, face_east, face_west );
				new_cube .add_to_mesh ( *this );

				// propagate 'front'
				face_west .remove_from_mesh ( front_west );
				face_south .remove_from_mesh ( front_south );
				face_down .remove_from_mesh ( front_down );
				face_east .reverse() .add_to_mesh ( front_west );
				face_north .reverse() .add_to_mesh ( front_south );
				face_up .reverse() .add_to_mesh ( front_down );

				// move 'seg1' and 'seg2' towards east (right)
				seg1 = edge_NS_down_E .reverse();
				spin1 += edge_WE_down_N .winding();
				seg2 = edge_NS_up_E .reverse();
				spin2 = 0;
	
				// move 'seg_south' and 'seg_north' towards east (right)
				sq_bdry = south .cell_behind ( seg_south, tag::surely_exists ) .boundary();
				// 'sq' is above 'seg_south'
				seg_south  = sq_bdry .cell_in_front_of ( seg_south .tip(), tag::surely_exists );
				// now 'seg_south' points up
				sq_bdry = south .cell_in_front_of ( seg_south, tag::surely_exists ) .boundary();
				// 'sq' is eastern (to the right) to 'seg_south'
				seg_south  =
					sq_bdry .cell_in_front_of ( seg_south .base() .reverse(), tag::surely_exists );
				spin_south += seg_south .winding();
				// 'seg_south' points again towards east (right)
				sq_bdry = north .cell_in_front_of ( seg_north, tag::surely_exists ) .boundary();
				// 'sq' is above 'seg_north'
				seg_north  = sq_bdry .cell_behind ( seg_north .tip(), tag::surely_exists );
				// now 'seg_north' points down
				sq_bdry = north .cell_in_front_of ( seg_north, tag::surely_exists ) .boundary();
				// 'sq' is eastern (to the right) to 'seg_north'
				seg_north  = sq_bdry .cell_behind
					( seg_north .tip(), tag::surely_exists ) .reverse ( tag::surely_exists );
				spin_north += seg_north .winding();
				// 'seg_north' points again towards east (right)
			}  // end of movement  west -> east

			// build the last (eastern) cube of this row

			// recall that 'seg1' and 'seg2' are two horizontal parallel segments
			// pointing towards north, 'seg1' lying one square below 'seg2'
			Cell face_west = front_west .cell_behind ( seg2, tag::surely_exists );
			Cell face_down = front_down .cell_behind ( seg1, tag::surely_exists );
			Cell corner_down_NW = seg1 .tip();
			Cell corner_up_NW = seg2 .tip();
			Cell edge_ud_NW =
				face_west .boundary() .cell_behind ( corner_down_NW, tag::surely_exists );
			Cell edge_WE_down_N =
				face_down .boundary() .cell_in_front_of ( corner_down_NW, tag::surely_exists );
			Cell corner_down_NE = edge_WE_down_N .tip();
			Cell corner_down_SW = seg1 .base() .reverse();
			Cell edge_EW_down_S =
				face_down .boundary() .cell_behind ( corner_down_SW, tag::surely_exists );
			Cell face_south = front_south .cell_in_front_of ( edge_EW_down_S, tag::surely_exists );
			Cell corner_up_SW = seg2 .base() .reverse();
			Cell edge_EW_up_S =
				face_south .boundary() .cell_behind ( corner_up_SW, tag::surely_exists );
			Cell corner_up_SE = edge_EW_up_S .base() .reverse();
			Cell edge_du_SE = face_south .boundary() .cell_behind ( corner_up_SE );
			Cell corner_down_SE = edge_du_SE .base() .reverse();
			Cell edge_NS_down_E = face_down .boundary()
				.cell_behind ( corner_down_SE, tag::surely_exists );
			Cell face_east = east .cell_in_front_of ( edge_NS_down_E, tag::surely_exists );
			Cell edge_du_NE = face_east .boundary()
				.cell_in_front_of ( edge_NS_down_E .base() .reverse(), tag::surely_exists );
			Cell new_ver = edge_du_NE .tip();  // well, not really new ...
			Cell edge_NS_up_E =
				face_east .boundary() .cell_in_front_of ( new_ver, tag::surely_exists );
			Cell edge_WE_up_N ( tag::segment, corner_up_NW .reverse(), new_ver );
			edge_WE_up_N .winding() =
				edge_ud_NW .winding() + edge_WE_down_N .winding() + edge_du_NE .winding();
			//  seg2  can be seen as  edge_SN_up_W
			assert ( edge_WE_up_N .winding() ==
			         - seg2 .winding() - edge_EW_up_S .winding() - edge_NS_up_E .winding() );
			Cell face_up ( tag::square, seg2 .reverse(),         edge_EW_up_S .reverse(),
			                            edge_NS_up_E .reverse(), edge_WE_up_N .reverse() );
			Cell face_north ( tag::square, edge_du_NE .reverse(), edge_WE_down_N .reverse(),
			                               edge_ud_NW .reverse(), edge_WE_up_N              );
			Cell new_cube ( tag::cube, face_up, face_down,
			                           face_south, face_north, face_east, face_west );
			new_cube .add_to_mesh ( *this );

			// propagate 'front'
			face_west  .remove_from_mesh ( front_west );
			face_south .remove_from_mesh ( front_south );
			face_down  .remove_from_mesh ( front_down );
			face_north .reverse() .add_to_mesh ( front_south );
			face_up    .reverse() .add_to_mesh ( front_down );

			// advance 'seg_west' and 'seg_east' towards north
			sq_bdry = west .cell_in_front_of ( seg_west, tag::surely_exists ) .boundary();
			// square is above seg_west
			seg_west  = sq_bdry .cell_behind ( seg_west .tip(), tag::surely_exists );
			// now 'seg_west' points down
			sq_bdry = west .cell_in_front_of ( seg_west, tag::surely_exists ) .boundary();
			seg_west  = sq_bdry .cell_behind ( seg_west .tip(), tag::surely_exists ) .reverse();
			spin_west += seg_west .winding();
			// 'seg_west' points again towards north
			sq_bdry = east .cell_behind ( seg_east, tag::surely_exists ) .boundary();
			// square is above seg_east
			seg_east  = sq_bdry .cell_in_front_of ( seg_east .tip(), tag::surely_exists );
			// now 'seg_east' points up
			sq_bdry = east .cell_in_front_of ( seg_east, tag::surely_exists ) .boundary();
			seg_east  = sq_bdry .cell_in_front_of ( seg_east .base() .reverse(), tag::surely_exists );
			spin_east += seg_east .winding();
			// 'seg_east' points again towards north
			
			// move 'seg_down_west' and 'seg_up_west' towards north
			seg_up_west = up_west .cell_behind ( seg_up_west .tip(), tag::surely_exists )
			              .reverse ( tag::surely_exists );
			spin_up_west += seg_up_west .winding();
			seg_down_west = down_west .cell_in_front_of ( seg_down_west .tip(), tag::surely_exists );
			spin_down_west += seg_down_west .winding();
		}  // end of movement  south -> north

		// build the last (northern) row of this layer

		// define 'seg1' and 'seg2'
		Cell seg2 = seg_west;  // horizontal, points towards north
		sq_bdry = west .cell_behind ( seg_west, tag::surely_exists ) .boundary();
		Cell seg1 = sq_bdry .cell_in_front_of ( seg2 .tip(), tag::surely_exists );
		// seg1 is vertical, points down
		seg1 = sq_bdry .cell_in_front_of ( seg1 .tip(), tag::surely_exists )
		       .reverse ( tag::surely_exists );
		// seg1 is horizontal now, points towards north, it lies one square below seg2
		// these two segments will move together towards east, keeping this relative position

		for ( size_t i_east = 1; i_east < nb_EW; i_east ++ )    // counts jumps toward east
		{ // we want to build a new cube
			// four faces exist already
			// we need to build two new faces
			// before that, one new segment (no new vertex)
			
			// recall that 'seg1' and 'seg2' are two horizontal parallel segments
			// pointing towards north, 'seg1' lying one square below 'seg2'
			Cell face_west = front_west .cell_behind ( seg2, tag::surely_exists );
			Cell face_down = front_down .cell_behind ( seg1, tag::surely_exists );
			Cell corner_down_NW = seg1 .tip();
			Cell corner_up_NW = seg2 .tip();
			Cell edge_ud_NW =
				face_west .boundary() .cell_behind ( corner_down_NW, tag::surely_exists );
			Cell edge_WE_down_N =
				face_down .boundary() .cell_in_front_of ( corner_down_NW, tag::surely_exists );
			Cell corner_down_NE = edge_WE_down_N .tip();
			Cell corner_down_SW = seg1 .base() .reverse();
			Cell edge_EW_down_S =
				face_down .boundary() .cell_behind ( corner_down_SW, tag::surely_exists );
			Cell face_south = front_south .cell_in_front_of ( edge_EW_down_S, tag::surely_exists );
			Cell corner_up_SW = seg2 .base() .reverse();
			Cell edge_EW_up_S =
				face_south .boundary() .cell_behind ( corner_up_SW, tag::surely_exists );
			Cell corner_up_SE = edge_EW_up_S .base() .reverse();
			Cell edge_du_SE = face_south .boundary() .cell_behind ( corner_up_SE );
			Cell corner_down_SE = edge_du_SE .base() .reverse();
			Cell edge_NS_down_E = face_down .boundary()
				.cell_behind ( corner_down_SE, tag::surely_exists );
			Cell face_north = north .cell_in_front_of ( edge_ud_NW, tag::surely_exists );
			Cell edge_WE_up_N = face_north .boundary()
				.cell_in_front_of ( edge_ud_NW .base() .reverse(), tag::surely_exists );
			Cell new_ver = edge_WE_up_N .tip();  // well, not really new ...
			Cell edge_ud_NE = face_north .boundary() .cell_in_front_of ( new_ver, tag::surely_exists );
			Cell edge_NS_up_E ( tag::segment, new_ver .reverse(), corner_up_SE );
			edge_NS_up_E .winding() =
				edge_ud_NE .winding() + edge_NS_down_E .winding() + edge_du_SE .winding();
			assert ( edge_NS_up_E .winding() ==
							 - edge_WE_up_N .winding() - seg2 .winding() - edge_EW_up_S .winding() );
			Cell face_up ( tag::square, seg2 .reverse(),         edge_EW_up_S .reverse(),
			                            edge_NS_up_E .reverse(), edge_WE_up_N .reverse() );
			Cell face_east  ( tag::square, edge_NS_up_E,              edge_du_SE .reverse(),
			                               edge_NS_down_E .reverse(), edge_ud_NE .reverse() );
			Cell new_cube ( tag::cube, face_up, face_down,
			                           face_south, face_north, face_east, face_west );
			new_cube .add_to_mesh ( *this );

			// propagate 'front'
			face_west  .remove_from_mesh ( front_west );
			face_south .remove_from_mesh ( front_south );
			face_down  .remove_from_mesh ( front_down );
			face_east .reverse() .add_to_mesh ( front_west );
			face_up   .reverse() .add_to_mesh ( front_down );

			// move 'seg1' and 'seg2' towards east (right)
			seg1 = edge_NS_down_E .reverse();
			seg2 = edge_NS_up_E .reverse();

		}  // end of movement  west -> east

		// build the last (eastern) cube of this last row of the current layer

		// recall that 'seg1' and 'seg2' are two horizontal parallel segments
		// pointing towards north, 'seg1' lying one square below 'seg2'
		Cell face_west = front_west .cell_behind ( seg2, tag::surely_exists );
		Cell face_down = front_down .cell_behind ( seg1, tag::surely_exists );
		Cell corner_down_NW = seg1 .tip();
		Cell corner_up_NW = seg2 .tip();
		Cell edge_ud_NW = face_west .boundary() .cell_behind ( corner_down_NW, tag::surely_exists );
		Cell edge_WE_down_N =
			face_down .boundary() .cell_in_front_of ( corner_down_NW, tag::surely_exists );
		Cell corner_down_NE = edge_WE_down_N .tip();
		Cell corner_down_SW = seg1 .base() .reverse();
		Cell edge_EW_down_S =
			face_down .boundary() .cell_behind ( corner_down_SW, tag::surely_exists );
		Cell face_south = front_south .cell_in_front_of ( edge_EW_down_S, tag::surely_exists );
		Cell corner_up_SW = seg2 .base() .reverse();
		Cell edge_EW_up_S = face_south .boundary() .cell_behind ( corner_up_SW, tag::surely_exists );
		Cell corner_up_SE = edge_EW_up_S .base() .reverse();
		Cell edge_du_SE = face_south .boundary() .cell_behind ( corner_up_SE );
		Cell corner_down_SE = edge_du_SE .base() .reverse();
		Cell edge_NS_down_E = face_down .boundary()
			.cell_behind ( corner_down_SE, tag::surely_exists );
		Cell face_east = east .cell_in_front_of ( edge_NS_down_E, tag::surely_exists );
		Cell edge_du_NE = face_east .boundary()
			.cell_in_front_of ( edge_NS_down_E .base() .reverse(), tag::surely_exists );
		Cell face_north = north .cell_in_front_of ( edge_ud_NW, tag::surely_exists );
		Cell edge_WE_up_N = face_north .boundary()
			.cell_in_front_of ( edge_ud_NW .base() .reverse(), tag::surely_exists );
		Cell new_ver = edge_WE_up_N .tip();  // well, not really new ...
		Cell edge_NS_up_E = face_east .boundary() .cell_in_front_of ( new_ver, tag::surely_exists );
		Cell face_up ( tag::square, seg2 .reverse(),         edge_EW_up_S .reverse(),
		                            edge_NS_up_E .reverse(), edge_WE_up_N .reverse() );
		Cell new_cube ( tag::cube, face_up, face_down,
										face_south, face_north, face_east, face_west );
		new_cube .add_to_mesh ( *this );

		// propagate 'front'
		face_west  .remove_from_mesh ( front_west );
		face_south .remove_from_mesh ( front_south );
		face_down  .remove_from_mesh ( front_down );
		face_up .reverse() .add_to_mesh ( front_down );

	}  // end of movement upwards

	// build the last (upper) layer

	Cell seg_west = up_west .cell_behind ( ver_up_south_west, tag::surely_exists )
		.reverse ( tag::surely_exists );

	for ( size_t i_north = 1; i_north < nb_SN; i_north ++ )    // counts jumps toward north
	{	// define 'seg1' and 'seg2'
		Cell seg2 = seg_west;  // horizontal, points towards north
		Mesh sq_bdry = west .cell_behind ( seg_west, tag::surely_exists ) .boundary();
		Cell seg1 = sq_bdry .cell_in_front_of ( seg2 .tip(), tag::surely_exists );
		// seg1 is vertical, points down
		seg1 = sq_bdry .cell_in_front_of ( seg1 .tip(), tag::surely_exists )
			.reverse ( tag::surely_exists );
		// seg1 is horizontal now, points towards north, it lies one square below seg2
		// these two segments will move together towards east, keeping this relative position

		for ( size_t i_east = 1; i_east < nb_EW; i_east ++ )    // counts jumps toward east
		{	// we want to build a new cube
			// four faces exist already
			// we need to build two new faces
			// one new segment, no new vertex
				
			// recall that 'seg1' and 'seg2' are two horizontal parallel segments
			// pointing towards north, 'seg1' lying one square below 'seg2'
			Cell face_west = front_west .cell_behind ( seg2, tag::surely_exists );
			Cell face_down = front_down .cell_behind ( seg1, tag::surely_exists );
			Cell corner_down_NW = seg1 .tip();
			Cell corner_up_NW = seg2 .tip();
			Cell edge_ud_NW =
				face_west .boundary() .cell_behind ( corner_down_NW, tag::surely_exists );
			Cell edge_WE_down_N =
				face_down .boundary() .cell_in_front_of ( corner_down_NW, tag::surely_exists );
			Cell corner_down_NE = edge_WE_down_N .tip();
			Cell corner_down_SW = seg1 .base() .reverse();
			Cell edge_EW_down_S =
				face_down .boundary() .cell_behind ( corner_down_SW, tag::surely_exists );
			Cell face_south = front_south .cell_in_front_of ( edge_EW_down_S, tag::surely_exists );
			Cell corner_up_SW = seg2 .base() .reverse();
			Cell edge_EW_up_S =
				face_south .boundary() .cell_behind ( corner_up_SW, tag::surely_exists );
			Cell corner_up_SE = edge_EW_up_S .base() .reverse();
			Cell edge_du_SE = face_south .boundary() .cell_behind ( corner_up_SE );
			Cell corner_down_SE = edge_du_SE .base() .reverse();
			Cell edge_NS_down_E =
				face_down .boundary() .cell_behind ( corner_down_SE, tag::surely_exists );
			Cell face_up = up .cell_in_front_of ( seg2, tag::surely_exists );
			Cell edge_WE_up_N = face_up .boundary() .cell_behind ( seg2 .tip(), tag::surely_exists )
			                    .reverse ( tag::surely_exists );
			Cell new_ver = edge_WE_up_N .tip();  // well, not really new ...
			Cell edge_NS_up_E = face_up .boundary() .cell_behind ( new_ver, tag::surely_exists )
			                    .reverse ( tag::surely_exists );
			Cell edge_ud_NE ( tag::segment, new_ver .reverse(), corner_down_NE );
			edge_ud_NE .winding() =
				- edge_WE_up_N .winding() + edge_ud_NW .winding() + edge_WE_down_N .winding();
			assert ( edge_ud_NE .winding() ==
							 edge_NS_up_E .winding() - edge_du_SE .winding() - edge_NS_down_E .winding() );
			Cell face_north ( tag::square, edge_ud_NE,            edge_WE_down_N .reverse(),
			                               edge_ud_NW .reverse(), edge_WE_up_N              );
			Cell face_east  ( tag::square, edge_NS_up_E,              edge_du_SE .reverse(),
			                               edge_NS_down_E .reverse(), edge_ud_NE .reverse() );
			Cell new_cube ( tag::cube, face_up, face_down,
			                           face_south, face_north, face_east, face_west );
			new_cube .add_to_mesh ( *this );

			// propagate 'front'
			face_west  .remove_from_mesh ( front_west );
			face_south .remove_from_mesh ( front_south );
			face_down  .remove_from_mesh ( front_down );
			face_east  .reverse() .add_to_mesh ( front_west );
			face_north .reverse() .add_to_mesh ( front_south );

			// move 'seg1' and 'seg2' towards east (right)
			seg1 = edge_NS_down_E .reverse();
			seg2 = edge_NS_up_E .reverse();
	
		}  // end of movement  west -> east

		// build the last cube of this row (within the last layer)

		// recall that 'seg1' and 'seg2' are two horizontal parallel segments
		// pointing towards north, 'seg1' lying one square below 'seg2'
		Cell face_west = front_west .cell_behind ( seg2, tag::surely_exists );
		Cell face_down = front_down .cell_behind ( seg1, tag::surely_exists );
		Cell corner_down_NW = seg1 .tip();
		Cell corner_up_NW = seg2 .tip();
		Cell edge_ud_NW = face_west .boundary() .cell_behind ( corner_down_NW, tag::surely_exists );
		Cell edge_WE_down_N =
			face_down .boundary() .cell_in_front_of ( corner_down_NW, tag::surely_exists );
		Cell corner_down_NE = edge_WE_down_N .tip();
		Cell corner_down_SW = seg1 .base() .reverse();
		Cell edge_EW_down_S =
			face_down .boundary() .cell_behind ( corner_down_SW, tag::surely_exists );
		Cell face_south = front_south .cell_in_front_of ( edge_EW_down_S, tag::surely_exists );
		Cell corner_up_SW = seg2 .base() .reverse();
		Cell edge_EW_up_S = face_south .boundary() .cell_behind ( corner_up_SW, tag::surely_exists );
		Cell corner_up_SE = edge_EW_up_S .base() .reverse();
		Cell edge_du_SE = face_south .boundary() .cell_behind ( corner_up_SE );
		Cell corner_down_SE = edge_du_SE .base() .reverse();
		Cell edge_NS_down_E = face_down .boundary()
			.cell_behind ( corner_down_SE, tag::surely_exists );
		Cell face_east = east .cell_in_front_of ( edge_NS_down_E, tag::surely_exists );
		Cell edge_du_NE = face_east .boundary()
			.cell_in_front_of ( edge_NS_down_E .base() .reverse(), tag::surely_exists );
		Cell new_ver = edge_du_NE .tip();  // well, not really new ...
		Cell edge_NS_up_E = face_east .boundary() .cell_in_front_of ( new_ver, tag::surely_exists );
		Cell face_up = up .cell_in_front_of ( seg2, tag::surely_exists );
		Cell edge_WE_up_N = face_up .boundary() .
			cell_behind ( seg2 .tip(), tag::surely_exists ) .reverse ( tag::surely_exists );
		assert ( new_ver == edge_WE_up_N .tip() );
		Cell face_north ( tag::square, edge_du_NE .reverse(), edge_WE_down_N .reverse(),
	                                 edge_ud_NW .reverse(), edge_WE_up_N              );
		Cell new_cube ( tag::cube, face_up, face_down,
										face_south, face_north, face_east, face_west );
		new_cube .add_to_mesh ( *this );

		// propagate 'front'
		face_west  .remove_from_mesh ( front_west );
		face_south .remove_from_mesh ( front_south );
		face_down  .remove_from_mesh ( front_down );
		face_north .reverse() .add_to_mesh ( front_south );

		seg_west = up_west .cell_behind ( seg_west .tip(), tag::surely_exists )
			.reverse ( tag::surely_exists );
	}  // end of movement  south -> north

	// build the last (northern) row of this last layer

	// define 'seg1' and 'seg2'
	Cell seg2 = seg_west;  // horizontal, points towards north
	Mesh sq_bdry = west .cell_behind ( seg_west, tag::surely_exists ) .boundary();
	Cell seg1 = sq_bdry .cell_in_front_of ( seg2 .tip(), tag::surely_exists );
	// seg1 is vertical, points down
	seg1 = sq_bdry .cell_in_front_of ( seg1 .tip(), tag::surely_exists )
		.reverse ( tag::surely_exists );
	// seg1 is horizontal now, points towards north, it lies one square below seg2
	// these two segments will move together towards east, keeping this relative position

	for ( size_t i_east = 1; i_east < nb_EW; i_east ++ )    // counts jumps toward east
	{	// we want to build a new cube
		// five faces exist already
		// we need to build one new face
		// no new segment, no new vertex
			
		// recall that 'seg1' and 'seg2' are two horizontal parallel segments
		// pointing towards north, 'seg1' lying one square below 'seg2'
		Cell face_west = front_west .cell_behind ( seg2, tag::surely_exists );
		Cell face_down = front_down .cell_behind ( seg1, tag::surely_exists );
		Cell corner_down_NW = seg1 .tip();
		Cell corner_up_NW = seg2 .tip();
		Cell edge_ud_NW = face_west .boundary() .cell_behind ( corner_down_NW, tag::surely_exists );
		Cell edge_WE_down_N =
			face_down .boundary() .cell_in_front_of ( corner_down_NW, tag::surely_exists );
		Cell corner_down_NE = edge_WE_down_N .tip();
		Cell corner_down_SW = seg1 .base() .reverse();
		Cell edge_EW_down_S =
			face_down .boundary() .cell_behind ( corner_down_SW, tag::surely_exists );
		Cell face_south = front_south .cell_in_front_of ( edge_EW_down_S, tag::surely_exists );
		Cell corner_up_SW = seg2 .base() .reverse();
		Cell edge_EW_up_S = face_south .boundary() .cell_behind ( corner_up_SW, tag::surely_exists );
		Cell corner_up_SE = edge_EW_up_S .base() .reverse();
		Cell edge_du_SE = face_south .boundary() .cell_behind ( corner_up_SE );
		Cell corner_down_SE = edge_du_SE .base() .reverse();
		Cell edge_NS_down_E =
			face_down .boundary() .cell_behind ( corner_down_SE, tag::surely_exists );
		Cell face_north = north .cell_in_front_of ( edge_ud_NW, tag::surely_exists );
		Cell edge_WE_up_N = face_north .boundary()
			.cell_in_front_of ( edge_ud_NW .base() .reverse(), tag::surely_exists );
		Cell new_ver = edge_WE_up_N .tip();  // well, not really new ...
		Cell edge_ud_NE = face_north .boundary() .cell_in_front_of ( new_ver, tag::surely_exists );
		Cell face_up = up .cell_in_front_of ( seg2, tag::surely_exists );
		assert ( edge_WE_up_N == face_up .boundary()
		    .cell_behind ( seg2 .tip(), tag::surely_exists ) .reverse ( tag::surely_exists ) );
		Cell edge_NS_up_E = face_up .boundary()
			.cell_behind ( new_ver, tag::surely_exists ) .reverse ( tag::surely_exists );
		Cell face_east  ( tag::square, edge_NS_up_E,              edge_du_SE .reverse(),
		                               edge_NS_down_E .reverse(), edge_ud_NE .reverse() );
		Cell new_cube ( tag::cube, face_up, face_down,
		                           face_south, face_north, face_east, face_west );
		new_cube .add_to_mesh ( *this );

		// propagate 'front'
		face_west  .remove_from_mesh ( front_west );
		face_south .remove_from_mesh ( front_south );
		face_down  .remove_from_mesh ( front_down );
		face_east .reverse() .add_to_mesh ( front_west );

		// move 'seg1' and 'seg2' towards east (right)
		seg1 = edge_NS_down_E .reverse();
		seg2 = edge_NS_up_E .reverse();

	}  // end of movement  west -> east

	// build the last cube of this last row of the last layer

	// recall that 'seg1' and 'seg2' are two horizontal parallel segments
	// pointing towards north, 'seg1' lying one square below 'seg2'
	Cell face_west = front_west .cell_behind ( seg2, tag::surely_exists );
	Cell face_down = front_down .cell_behind ( seg1, tag::surely_exists );
	Cell corner_down_NW = seg1 .tip();
	Cell corner_up_NW = seg2 .tip();
	Cell edge_ud_NW = face_west .boundary() .cell_behind ( corner_down_NW, tag::surely_exists );
	Cell edge_WE_down_N =
		face_down .boundary() .cell_in_front_of ( corner_down_NW, tag::surely_exists );
	Cell corner_down_NE = edge_WE_down_N .tip();
	Cell corner_down_SW = seg1 .base() .reverse();
	Cell edge_EW_down_S =
		face_down .boundary() .cell_behind ( corner_down_SW, tag::surely_exists );
	Cell face_south = front_south .cell_in_front_of ( edge_EW_down_S, tag::surely_exists );
	Cell corner_up_SW = seg2 .base() .reverse();
	Cell edge_EW_up_S = face_south .boundary() .cell_behind ( corner_up_SW, tag::surely_exists );
	Cell corner_up_SE = edge_EW_up_S .base() .reverse();
	Cell edge_du_SE = face_south .boundary() .cell_behind ( corner_up_SE );
	Cell corner_down_SE = edge_du_SE .base() .reverse();
	Cell edge_NS_down_E = face_down .boundary()
		.cell_behind ( corner_down_SE, tag::surely_exists );
	Cell face_east = east .cell_in_front_of ( edge_NS_down_E, tag::surely_exists );
	Cell face_north = north .cell_in_front_of ( edge_ud_NW, tag::surely_exists );
	Cell face_up = up .cell_in_front_of ( seg2, tag::surely_exists );
	Cell new_cube ( tag::cube, face_up, face_down,
									face_south, face_north, face_east, face_west );
	new_cube .add_to_mesh ( *this );

	#ifndef NDEBUG
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = this->iterator ( tag::over_cells_of_dim, 2 );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell face = *it;
		Manifold::Action s;
		Mesh::Iterator itt = face .boundary() .iterator ( tag::over_cells, tag::of_max_dim );
		for ( itt .reset(); itt .in_range(); itt++ )  s += (*itt) .winding();
		assert ( s == 0 );                                                                    }
	} // just a block of code for hiding 'it'
	#endif
	
} // end of Mesh::build with tag::hexahedron and tag::winding

//----------------------------------------------------------------------------------//


void Mesh::export_to_file
( const tag::Msh &, std::string f, Cell::Numbering & ver_numbering ) const

// 'numb_map' should begin at 0
// we add 1 to each number because gmsh seems to prefer numbers to begin at 1

{	// we use the current manifold
	Manifold space = Manifold::working;
	assert ( space .exists() );
	Function coord = space .coordinates();

	std::ofstream file_msh (f);
	file_msh << "$MeshFormat" << std::endl << "4.1 0 4" << std::endl;
	file_msh << "$EndMeshFormat" << std::endl;

	// "Entities" section is optional in gmsh format 4.1, we skip it

	size_t nb_nodes = this->number_of ( tag::cells_of_dim, 0 );
	file_msh << "$Nodes" << std::endl;
	file_msh << "1 " << nb_nodes << " 1 " << nb_nodes << std::endl;

	// only one entity block
	int top_dim = this->dim();
	file_msh << top_dim << " 0 0 " << nb_nodes << std::endl;

	{ // just a block for hiding 'it'
	Mesh::Iterator it = this->iterator ( tag::over_vertices );
	for ( it .reset() ; it .in_range(); it++ )
	{	Cell p = *it;
		file_msh << ver_numbering [p] + 1 << std::endl;  }

	} { // just a block for hiding names : it, x, y
	Mesh::Iterator it = this->iterator ( tag::over_vertices );
	Function x = coord [0], y = coord [1];
	if (coord .nb_of_components() == 2)
	{	for ( it .reset() ; it .in_range(); it++ )
		{	Cell p = *it;
			file_msh << x(p) << " " << y(p) << " " << 0 << std::endl;  }  }
	else
	{	assert ( coord .nb_of_components() == 3 );
		Function z = coord [2];
		for ( it.reset() ; it.in_range(); it++ )
		{	Cell p = *it;
			file_msh << x(p) << " " << y(p) << " " << z(p) << std::endl;  }  }
	file_msh << "$EndNodes" << std::endl;
	} // just a block for hiding names 

	file_msh << "$Elements" << std::endl;

	// we must create one block for each type of cell :
	// segments, triangles, quadrilaterals, tetrahedra, cubes
	size_t nb_seg = 0, nb_tri = 0, nb_quad = 0, nb_tetra = 0, nb_cub = 0, nb_prism = 0;
	{ // just a block for hiding 'it'
	Mesh::Iterator it = this->iterator ( tag::over_cells, tag::of_max_dim );
	if ( top_dim == 1 )
	for ( it .reset() ; it .in_range(); it++ )
	{	Cell cll = *it;
		assert ( cll .dim() == 1 );
		nb_seg ++;                  }
	if ( top_dim == 2 )
	for ( it .reset() ; it .in_range(); it++ )
	{	Cell cll = *it;
		assert ( cll .dim() == 2 );
		size_t n = cll .boundary() .number_of ( tag::cells_of_max_dim );
		if ( n == 3 ) nb_tri ++;
		else  {  assert ( n == 4 );  nb_quad ++;  }                      }
	if ( top_dim == 3 )
	for ( it .reset() ; it .in_range(); it++ )
	{	Cell cll = *it;
		assert ( cll .dim() == 3 );
		size_t n = cll .boundary() .number_of ( tag::cells_of_max_dim );
		if ( n == 4 ) nb_tetra ++;
		else if ( n == 5 ) nb_prism ++;
		else  {  assert ( n == 6 );  nb_cub ++;  }                      }
	} // just a block for hiding 'it'

	// since there is no request to create specific entities,
	// we just create one entity per element type of maximum topological dimension
	size_t nb_blocks = 0;
	if ( top_dim == 1 )   nb_blocks = 1;
	if ( top_dim == 2 )
	{	if ( nb_tri > 0 )   nb_blocks++;
		if ( nb_quad > 0 )  nb_blocks++;  }
	if ( top_dim == 3 )
	{	if ( nb_tetra > 0 ) nb_blocks++;
		if ( nb_prism > 0 ) nb_blocks++;
		if ( nb_cub > 0 )   nb_blocks++;   }
	
	if ( top_dim == 1 )
	{	assert ( nb_tri == 0 );
		assert ( nb_quad == 0 );
		assert ( nb_tetra == 0 );
		assert ( nb_prism == 0 );
		assert ( nb_cub == 0 );
		assert ( nb_seg == this->number_of ( tag::cells_of_max_dim ) );
		assert ( nb_blocks == 1 );
		file_msh << "1 " << nb_seg << " 1 " << nb_seg << std::endl;
		file_msh << "1 0 1 " << nb_seg << std::endl;
		Mesh::Iterator it = this->iterator ( tag::over_segments );
		size_t counter = 0;
		for ( it .reset() ; it .in_range(); it++)
		{	++counter;
			Cell elem = *it;
			file_msh << counter << " ";
			Cell A = elem .base() .reverse();
			file_msh << ver_numbering [A] + 1 << " ";
			Cell B = elem .tip();
			file_msh << ver_numbering [B] + 1 << std::endl;    }          }
	else if ( top_dim == 2 )
	{	assert ( nb_seg == 0 );
		assert ( nb_tetra == 0 );
		assert ( nb_prism == 0 );
		assert ( nb_cub == 0 );
		assert ( nb_tri + nb_quad == this->number_of ( tag::cells_of_max_dim ) );
		assert ( nb_blocks <= 2 );
		Mesh::Iterator it = this->iterator ( tag::over_cells_of_max_dim );
		file_msh << nb_blocks << " " << nb_tri + nb_quad << " 1 " << nb_tri + nb_quad << std::endl;
		size_t counter = 0;
		if ( nb_tri > 0 )
		{	file_msh << "2 0 2 " << nb_tri << std::endl;
			for ( it .reset() ; it .in_range(); it++)
			{	Cell elem = *it;
				if ( elem .boundary() .number_of ( tag::cells_of_max_dim ) == 3 ) // a triangle
				{	++counter;
					file_msh << counter;
					Mesh::Iterator itt = elem .boundary()
						.iterator ( tag::over_vertices, tag::require_order );
					for ( itt .reset(); itt .in_range(); ++itt )
					{	Cell p = *itt;  file_msh << " " << ver_numbering [p] + 1;   }
					file_msh << std::endl;                                            }  }  }
		assert ( counter == nb_tri );
		if ( nb_quad > 0 )
		{	file_msh << "2 0 3 " << nb_quad << std::endl;
			for ( it .reset() ; it .in_range(); it++)
			{	Cell elem = *it;
				if ( elem .boundary() .number_of ( tag::cells_of_max_dim ) == 4 ) // a quadrilateral
				{	++counter;
					file_msh << counter;
					Mesh::Iterator itt = elem .boundary()
						.iterator ( tag::over_vertices, tag::require_order );
					for ( itt .reset(); itt .in_range(); ++itt )
					{	Cell p = *itt;  file_msh << " " << ver_numbering [p] + 1;   }
					file_msh << std::endl;                                            }  }  }
		assert ( counter == nb_tri + nb_quad );                                             }
	else
	{	assert ( top_dim == 3);
		assert ( nb_seg == 0 );
		assert ( nb_tri == 0 );
		assert ( nb_quad == 0 );
		assert ( nb_tetra + nb_cub + nb_prism == this->number_of ( tag::cells_of_max_dim ) );
		assert ( nb_blocks <= 3 );
		file_msh << nb_blocks << " " <<  nb_tetra + nb_prism + nb_cub
		         << " 1 " << nb_tetra + nb_prism + nb_cub << std::endl;
		Mesh::Iterator it = this->iterator ( tag::over_cells_of_max_dim );
		size_t counter = 0;
		if ( nb_tetra > 0 )
		{	file_msh << "3 0 4 " <<  nb_tetra << std::endl;
			for ( it .reset() ; it .in_range(); it++)
			{	++counter;
				Cell elem = *it;
				size_t n_faces = elem .boundary() .number_of ( tag::cells_of_max_dim );
				if ( n_faces == 4 ) // a tetrahedron
				{	file_msh << counter;
					// see http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
					// and http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
					// our orientation is opposite to the one used by gmsh
					// our normals point outwards, in gmsh normals point inwards
					// pick any face :
					Mesh::Iterator itt = elem .boundary() .iterator ( tag::over_cells_of_max_dim );
					itt .reset();  assert ( itt .in_range() );
					Cell ABC = *itt;
					Mesh::Iterator itv = ABC .boundary() .iterator ( tag::over_vertices, tag::backwards );
					for ( itv .reset(); itv .in_range(); ++itv )
					{	Cell p = *itv;  file_msh << " " << ver_numbering [p] + 1;   }
					// now the fourth vertex
					Mesh::Iterator its = ABC .boundary() .iterator ( tag::over_cells_of_max_dim );
					its .reset();  assert ( its .in_range() );
					Cell seg = *its;
					Cell face = elem .boundary() .cell_in_front_of ( seg );
					Cell other_seg = face .boundary() .cell_behind ( seg .tip() );
					Cell D = other_seg .base() .reverse();
					file_msh << " " << ver_numbering [D] + 1 << std::endl;
				}	}	}  // end of if, end of for, end of if
		assert ( counter == nb_tetra );
		if ( nb_cub > 0 )
		{	file_msh << "3 0 5 " <<  nb_cub << std::endl;
			for ( it .reset() ; it .in_range(); it++)
			{	++counter;
				Cell elem = *it;
				size_t n_faces = elem .boundary() .number_of ( tag::cells_of_max_dim );
				if ( n_faces == 6 )
				{	// parallelipiped = 8-node hexahedron = cube
					// see http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
					// and http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
					file_msh << counter;
					Mesh::Iterator itt = elem .boundary() .iterator ( tag::over_cells_of_max_dim );
					itt .reset();  Cell back = *itt; // square face behind the cube
					// back is 0321 in gmsh's documentation
					assert ( back .boundary() .number_of ( tag::cells_of_max_dim ) == 4 );
					Mesh::Iterator itv = back .boundary() .iterator ( tag::over_vertices, tag::backwards );
					// backwards because we want the vertices ordered as 0, 1, 2, 3
					// our orientation is opposite to the one used by gmsh
					// our normals point outwards, in gmsh normals point inwards
					itv .reset();  Cell ver_0 = *itv;
					Cell seg_03 = back .boundary() .cell_in_front_of ( ver_0 );
					for ( ; itv .in_range(); ++itv )
					{	Cell p = *itv;  file_msh << " " << ver_numbering [p] + 1;   }
					Cell left_wall = elem .boundary() .cell_in_front_of ( seg_03 ); // square face on the left
					// left_wall is 0473 in gmsh's documentation
					assert ( left_wall .boundary() .number_of ( tag::cells_of_max_dim ) == 4 );
					Cell seg_04 = left_wall .boundary() .cell_in_front_of ( ver_0 );
					Cell ver_4 = seg_04 .tip();
					Cell seg_47 = left_wall .boundary() .cell_in_front_of ( ver_4 );
					Cell front = elem .boundary() .cell_in_front_of ( seg_47 ); // square face in front
					// front is 4567 in gmsh's documentation
					assert ( front .boundary() .number_of ( tag::cells_of_max_dim ) == 4 );
					Mesh::Iterator itvv = front .boundary() .iterator
						( tag::over_vertices, tag::require_order );
					itvv .reset ( tag::start_at, ver_4 );
					for ( ; itvv .in_range(); ++itvv )
					{	Cell p = *itvv;  file_msh << " " << ver_numbering [p] + 1;   }
					file_msh << std::endl;
				}	}	}  // end of if, end of for, end of it
		assert ( counter == nb_tetra + nb_cub );
		if ( nb_prism > 0 )
		{	file_msh << "3 0 6 " <<  nb_prism << std::endl;
			for ( it .reset() ; it .in_range(); it++)
			{	++counter;
				Cell elem = *it;
				size_t n_faces = elem .boundary() .number_of ( tag::cells_of_max_dim );
				if ( n_faces == 5 )
				{	// triangular prism = 6-node prism
					// see http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
					// and http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
					file_msh << counter;
					Mesh::Iterator itt = elem .boundary() .iterator ( tag::over_cells_of_max_dim );
					size_t n_tri = 0, n_rect = 0;
					Cell floor ( tag::non_existent );  // temporary non-existent cell
					for( itt .reset(); itt .in_range(); ++itt )
					{	Cell face = *itt; // every face elem
						size_t n_edges = face .boundary() .number_of ( tag::cells_of_max_dim );
						if ( n_edges == 3 )  { n_tri++;   floor = face;             }
						else                 { n_rect++;  assert ( n_edges == 4 );  }         }
					assert ( n_tri == 2 );  assert ( n_rect == 3 );
					// 'floor' is 021 in gmsh's documentation
					assert ( floor .boundary() .number_of ( tag::cells_of_max_dim ) == 3 );
					Mesh::Iterator itv = floor .boundary() .iterator
						( tag::over_vertices, tag::backwards );
					// backwards because we want the vertices ordered as 0, 1, 2
					// our orientation is opposite to the one used by gmsh
					// our normals point outwards, in gmsh normals point inwards
					itv .reset();  Cell ver_0 = *itv;
					Cell seg_02 = floor .boundary() .cell_in_front_of ( ver_0 );
					for (  ; itv .in_range(); ++itv )
					{	Cell p = *itv;  file_msh << " " << ver_numbering [p] + 1;   }
					Cell right_wall = elem .boundary() .cell_in_front_of ( seg_02 );
					// right_wall is 0352 in gmsh's documentation
					assert ( right_wall .boundary() .number_of ( tag::cells_of_max_dim ) == 4 );
					Cell seg_03 = right_wall .boundary() .cell_in_front_of ( ver_0 );
					Cell ver_3 = seg_03 .tip();
					Cell seg_35 = right_wall .boundary() .cell_in_front_of ( ver_3 );
					Cell roof = elem .boundary() .cell_in_front_of ( seg_35 );
					// 'roof' is 345 in gmsh's documentation
					assert ( roof .boundary() .number_of ( tag::cells_of_max_dim ) == 3 );
					Mesh::Iterator itvv = roof .boundary() .iterator
						( tag::over_vertices, tag::require_order );
					for ( itvv .reset ( tag::start_at, ver_3 ); itvv .in_range(); ++itvv )
					{	Cell p = *itvv;  file_msh << " " << ver_numbering [p] + 1;  }
				file_msh << std::endl;
				}	}	}  // end of if, end of for, end of it
		assert ( counter == nb_tetra + nb_cub + nb_prism );                                    }
	file_msh << "$EndElements" << std::endl;
	
	if ( ! file_msh .good() )
	{	std::cerr << "error writing msh file" << std::endl;
		exit (1);                                           }

} // end of Mesh::export_to_file



void Mesh::export_to_file
( const tag::Msh &, std::string f, std::map < Cell, size_t > & numb_map ) const
	
// 'numb_map' should begin at 0
// later, 'this->export_to_file ( tag::msh, f, numbering )' will add 1 to each number

{	Cell::Numbering::Map numbering ( & numb_map );

	this->export_to_file ( tag::msh, f, numbering );

} // end of Mesh::export_to_file


void Mesh::export_to_file ( const tag::Msh &, std::string f ) const
	
// the numbering of vertices is produced on-the-fly

// we build a 'numbering' map beginning at 0
// later, 'this->export_to_file ( tag::msh, f, numbering )' will add 1 to each number
	
{	Cell::Numbering::Map numbering;

	Mesh::Iterator it = this->iterator ( tag::over_vertices );
	size_t counter = 0;
	for ( it.reset() ; it.in_range(); it++ )
	{	Cell p = *it;  numbering [p] = counter;  ++counter;  }

	this->export_to_file ( tag::msh, f, numbering );

} // end of Mesh::export_to_file

//----------------------------------------------------------------------------------//


namespace {  // anonymous namespace, mimics static linkage

inline int import_msh_2  // hidden in anonymous namespace
( const std::string filename, int geom_dim,
  std::vector < std::pair < std::pair < int, int >,
    std::vector < std::pair < size_t, std::vector < size_t > > > > > & ver_of_elem,
  std::map < size_t, Cell > & map_ver                                              )

{	assert ( ( geom_dim == 2 ) or ( geom_dim == 3 ) );
	double msh_version;

	std::ifstream ifstr ( filename );
	std::string s;
	std::getline ( ifstr, s );  assert ( trim_copy (s) == "$MeshFormat");
	int binary_file;  ifstr >> msh_version >> binary_file;
	assert ( msh_version < 4.);
	assert ( binary_file == 0 );
	std::getline ( ifstr, s );  // one more number and the endline
	std::getline ( ifstr, s );  assert ( trim_copy (s) == "$EndMeshFormat");

	std::getline ( ifstr, s );  assert ( trim_copy (s) == "$Nodes");

	Manifold RRd ( tag::Euclid, tag::of_dim, geom_dim );
	Function coords = RRd .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	
	size_t nb_nodes;  ifstr >> nb_nodes;
	std::getline ( ifstr, s );  assert ( trim_copy (s) == "");

	for ( size_t i = 0; i < nb_nodes; i++ )
	{	size_t j;  ifstr >> j;  assert ( j == i+1 );
		double xx, yy, zz;  ifstr >> xx >> yy >> zz;
		Cell V ( tag::vertex );
		coords [0] (V) = xx;
		coords [1] (V) = yy;
		if ( geom_dim == 3 )  coords [2] (V) = zz;
		else  assert ( zz == 0. );
		map_ver .insert ( std::pair < size_t, Cell > ( j, V ) );
		std::getline ( ifstr, s );  assert ( trim_copy (s) == "");   }
	std::getline ( ifstr, s );  assert ( trim_copy (s) == "$EndNodes");

	while ( trim_copy (s) != "$Elements")  std::getline ( ifstr, s );
	size_t nb_elem;
	// ver_of_elem == [ { { ent_tag, geom_type }, [ { cll_tag, [ node_tag ] } ] } ]
	// we do not keep the dimension, it is implied by the geometry type
	// we take ent_tag = 0, see assert k == 0 below
	std::vector < int > elem_size { 0, 2, 3, 4, 4, 8, 6, 5 };
	std::vector < int > elem_dim { 0, 1, 2, 2, 3, 3, 3, 3 };
	// 1 segment, 2 triangle, 3 quadrilateral, 4 tetrahedron, 5 cube, 6 prism, 7 pyramid
	int global_top_dim = 0;
	std::vector < int > elem_type_show_up { -1, -1, -1, -1, -1, -1, -1, -1 };
	assert ( elem_size .size() == elem_type_show_up .size() );
	ifstr >> nb_elem;
	std::getline ( ifstr, s );  assert ( trim_copy (s) == "");
	for ( size_t i = 0; i < nb_elem; i++ )
	{	size_t cll_tag;  ifstr >> cll_tag;
		int cll_type;  ifstr >> cll_type;
		if ( elem_dim [ cll_type ] > global_top_dim )  global_top_dim = elem_dim [ cll_type ];
		int k;  ifstr >> k;  ifstr >> k;  ifstr >> k;
		// we just ignore three integers
		assert ( cll_type < int ( elem_size .size() ) );
		if ( elem_type_show_up [ cll_type ] < 0 )
		{	elem_type_show_up [ cll_type ] = ver_of_elem .size();
			ver_of_elem .push_back ( { { 0, cll_type }, { } } );  }
		ver_of_elem [ elem_type_show_up [ cll_type ] ] .second .push_back ( { cll_tag, { } } );
		std::vector < size_t > & local_vec =
			ver_of_elem [ elem_type_show_up [ cll_type ] ] .second .back() .second;
		local_vec .reserve ( elem_size [ cll_type ] );
		for ( int jj = 0; jj < elem_size [ cll_type ]; jj ++ )
		{	size_t j;  ifstr >> j;  local_vec .push_back (j);  }
		std::getline ( ifstr, s );  assert ( trim_copy (s) == "" );                             }
	std::getline ( ifstr, s );  assert ( trim_copy (s) == "$EndElements");

	assert ( ( global_top_dim == 2 ) or ( global_top_dim == 3 ) );
	return global_top_dim;
	
}  // end of  import_msh_2

//----------------------------------------------------------------------------------//


inline int import_msh_4_one_mesh  // hidden in anonymous namespace
( const std::string filename, int geom_dim,
  std::vector < std::pair < std::pair < int, int >,
    std::vector < std::pair < size_t, std::vector < size_t > > > > > & ver_of_elem,
  std::map < size_t, Cell > & map_ver                                              )

{	assert ( ( geom_dim == 2 ) or ( geom_dim == 3 ) );
	double msh_version;

	std::ifstream ifstr ( filename );
	std::string s;
	std::getline ( ifstr, s );  assert ( trim_copy (s) == "$MeshFormat");
	int binary_file;  ifstr >> msh_version >> binary_file;
	assert ( msh_version >= 4.);
	assert ( binary_file == 0 );
	std::getline ( ifstr, s );  // one more number and the endline
	std::getline ( ifstr, s );  assert ( trim_copy (s) == "$EndMeshFormat");

	std::getline ( ifstr, s );
	if ( trim_copy (s) == "$PhysicalNames")
	{	// we skip section PhysicalEntities
		while ( trim_copy (s) != "$EndPhysicalNames")  std::getline ( ifstr, s );
		std::getline ( ifstr, s );                                                   }

	size_t nb_ent_0d = 0, nb_ent_1d = 0, nb_ent_2d = 0, nb_ent_3d = 0;
	if ( trim_copy (s) == "$Entities")
	{	ifstr >> nb_ent_0d >> nb_ent_1d >> nb_ent_2d >> nb_ent_3d;
		// we skip the rest of this section
		while ( trim_copy (s) != "$EndEntities")  std::getline ( ifstr, s );
		std::getline ( ifstr, s );                                           }
	// if there is no $Entities section, all nb_ent remain zero

	if ( trim_copy (s) == "$PartitionedEntities")
	{	// we skip section PartitionedEntities
		while ( trim_copy (s) != "$EndPartitionedEntities")  std::getline ( ifstr, s );
		std::getline ( ifstr, s );                                                      }

	assert ( trim_copy (s) == "$Nodes");

	Manifold RRd ( tag::Euclid, tag::of_dim, geom_dim );
	Function coords = RRd .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	
	size_t numEntityBlocks, minNodeTag, maxNodeTag, nb_nodes;
	ifstr >> numEntityBlocks >> nb_nodes >> minNodeTag >> maxNodeTag;
	std::getline ( ifstr, s );  assert ( trim_copy (s) == "" );

	for ( size_t ent = 0; ent < numEntityBlocks; ent++ )
	{	int top_dim, entity_tag, parametric;
		size_t nbnodes;
		ifstr >> top_dim >> entity_tag >> parametric >> nbnodes;
		std::getline ( ifstr, s );  assert ( trim_copy (s) == "" );
		std::vector < size_t > local_tag;
		for ( size_t i = 0; i < nbnodes; i++ )
		{	size_t j;  ifstr >> j;     // node tag
			std::getline ( ifstr, s ); assert ( trim_copy (s) == "" );
			assert ( j <= maxNodeTag );
			local_tag .push_back (j);                                  }
		for ( size_t i = 0; i < nbnodes; i++ )
		{	double xx, yy, zz;  ifstr >> xx >> yy >> zz;
			std::getline ( ifstr, s );  assert ( trim_copy (s) == "" );
			Cell V ( tag::vertex );
			coords [0] (V) = xx;
			coords [1] (V) = yy;
			if ( geom_dim == 3 ) coords [2] (V) = zz;
			else  assert ( zz == 0. );
			map_ver .insert ( { local_tag [i], V } );                   }  }
	std::getline ( ifstr, s );  assert ( trim_copy (s) == "$EndNodes");

	std::getline ( ifstr, s );
	assert ( trim_copy (s) == "$Elements");
	// ver_of_elem == [ { { ent_tag, geom_type }, [ { cll_tag, [ node_tag ] } ] } ]
	// we do not keep the dimension, it is implied by the geometry type
	// we take ent_tag = 0, see assert k == 0 below
	const std::vector < int > elem_size { 0, 2, 3, 4, 4, 8, 6, 5 };
	const std::vector < int > elem_dim { 0, 1, 2, 2, 3, 3, 3, 3 };
	// 1 segment, 2 triangle, 3 quadrilateral, 4 tetrahedron, 5 cube, 6 prism, 7 pyramid
	int global_top_dim = 0;
	std::vector < int > elem_type_show_up { -1, -1, -1, -1, -1, -1, -1, -1 };
	assert ( elem_size .size() == elem_type_show_up .size() );
	size_t minElemTag, maxElemTag, nb_elem;
	ifstr >> numEntityBlocks >> nb_elem >> minElemTag >> maxElemTag;
	for ( size_t ent = 0; ent < numEntityBlocks; ent++ )
	{	int top_dim;  ifstr >> top_dim;
		if ( global_top_dim < top_dim ) global_top_dim = top_dim;
		size_t j;  ifstr >> j;  //  entity dim, entity tag
		int cll_type;  ifstr >> cll_type >> nb_elem;
		assert ( cll_type < int ( elem_size .size() ) );
		if ( elem_type_show_up [ cll_type ] < 0 )
		{	elem_type_show_up [ cll_type ] = ver_of_elem .size();
			ver_of_elem .push_back ( { { 0, cll_type }, { } } );  }
		for ( size_t i = 0; i < nb_elem; i++ )
		{	size_t cll_tag;  ifstr >> cll_tag;
			ver_of_elem [ elem_type_show_up [ cll_type ] ] .second .push_back ( { cll_tag, { } } );
			std::vector < size_t > & local_vec =
				ver_of_elem [ elem_type_show_up [ cll_type ] ] .second .back() .second;
			local_vec .reserve ( elem_size [ cll_type ] );
			for ( int jj = 0; jj < elem_size [ cll_type ]; jj ++ )
			{	ifstr >> j;  local_vec .push_back (j);  }
			std::getline ( ifstr, s );  assert ( trim_copy (s) == "");                               }  }
	std::getline ( ifstr, s );  assert ( trim_copy (s) == "$EndElements");

	assert ( ( global_top_dim == 2 ) or ( global_top_dim == 3 ) );
	return global_top_dim;
	
}  // end of  import_msh_4_one_mesh

//----------------------------------------------------------------------------------//


inline void import_msh_4_composed  // hidden in anonymous namespace
( const std::string filename, int geom_dim,
  std::vector < std::pair < std::pair < int, int >,
    std::vector < std::pair < size_t, std::vector < size_t > > > > > & ver_of_elem,
  std::map < size_t, Cell > & map_ver                                              )

{	assert ( ( geom_dim == 2 ) or ( geom_dim == 3 ) );
	double msh_version;

	std::ifstream ifstr ( filename );
	std::string s;
	std::getline ( ifstr, s );  assert ( trim_copy (s) == "$MeshFormat");
	int binary_file;  ifstr >> msh_version >> binary_file;
	assert ( msh_version >= 4.);
	assert ( binary_file == 0 );
	std::getline ( ifstr, s );  // one more number and the endline
	std::getline ( ifstr, s );  assert ( trim_copy (s) == "$EndMeshFormat");

	std::getline ( ifstr, s );
	if ( trim_copy (s) == "$PhysicalNames")
	{	// we skip section PhysicalEntities
		while ( trim_copy (s) != "$EndPhysicalNames")  std::getline ( ifstr, s );
		std::getline ( ifstr, s );                                                   }

	size_t nb_ent_0d = 0, nb_ent_1d = 0, nb_ent_2d = 0, nb_ent_3d = 0;
	if ( trim_copy (s) == "$Entities")
	{	ifstr >> nb_ent_0d >> nb_ent_1d >> nb_ent_2d >> nb_ent_3d;
		// we skip the rest of this section
		while ( trim_copy (s) != "$EndEntities")  std::getline ( ifstr, s );
		std::getline ( ifstr, s );                                           }
	// if there is no $Entities section, all nb_ent remain zero

	if ( trim_copy (s) == "$PartitionedEntities")
	{	// we skip section PartitionedEntities
		while ( trim_copy (s) != "$EndPartitionedEntities")  std::getline ( ifstr, s );
		std::getline ( ifstr, s );                                                      }

	assert ( trim_copy (s) == "$Nodes");

	Manifold RRd ( tag::Euclid, tag::of_dim, geom_dim );
	Function coords = RRd .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	
	size_t numEntityBlocks, minNodeTag, maxNodeTag, nb_nodes;
	ifstr >> numEntityBlocks >> nb_nodes >> minNodeTag >> maxNodeTag;
	std::getline ( ifstr, s );  assert ( trim_copy (s) == "" );

	for ( size_t ent = 0; ent < numEntityBlocks; ent++ )
	{	int top_dim, entity_tag, parametric;
		size_t nbnodes;
		ifstr >> top_dim >> entity_tag >> parametric >> nbnodes;
		std::getline ( ifstr, s );  assert ( trim_copy (s) == "" );
		std::vector < size_t > local_tag;
		for ( size_t i = 0; i < nbnodes; i++ )
		{	size_t j;  ifstr >> j;     // node tag
			std::getline ( ifstr, s ); assert ( trim_copy (s) == "" );
			assert ( j <= maxNodeTag );
			local_tag .push_back (j);                                  }
		for ( size_t i = 0; i < nbnodes; i++ )
		{	double xx, yy, zz;  ifstr >> xx >> yy >> zz;
			std::getline ( ifstr, s );  assert ( trim_copy (s) == "" );
			Cell V ( tag::vertex );
			coords [0] (V) = xx;
			coords [1] (V) = yy;
			if ( geom_dim == 3 ) coords [2] (V) = zz;
			else  assert ( zz == 0. );
			map_ver .insert ( { local_tag [i], V } );                   }  }
	std::getline ( ifstr, s );  assert ( trim_copy (s) == "$EndNodes");

	std::getline ( ifstr, s );
	assert ( trim_copy (s) == "$Elements");
	// ver_of_elem == [ { { ent_tag, geom_type }, [ { cll_tag, [ node_tag ] } ] } ]
	// we do not keep the dimension, it is implied by the geometry type
	// we take ent_tag = 0, see assert k == 0 below
	const std::vector < int > elem_size { 1, 2, 3, 4, 4, 8, 6, 5 };
	const std::vector < int > elem_dim { 0, 1, 2, 2, 3, 3, 3, 3 };
	// 0 vertex (changed from 15), 1 segment, 2 triangle, 3 quadrilateral,
	// 4 tetrahedron, 5 cube, 6 prism, 7 pyramid
	size_t minElemTag, maxElemTag, nb_elem;
	ifstr >> numEntityBlocks >> nb_elem >> minElemTag >> maxElemTag;
	for ( size_t ent = 0; ent < numEntityBlocks; ent++ )
	{	int ent_dim, ent_tag;  ifstr >> ent_dim >> ent_tag;
		int cll_type;  ifstr >> cll_type >> nb_elem;
		// std::cout << "global.cpp line 3908, ent " << ent_tag << ", type " << cll_type << std::endl;
		if ( cll_type == 15 )  cll_type = 0;
		assert ( cll_type < int ( elem_size .size() ) );
		ver_of_elem .push_back ( { { ent_tag, cll_type }, { } } );
		for ( size_t i = 0; i < nb_elem; i++ )
		{	size_t cll_tag;  ifstr >> cll_tag;
			ver_of_elem .back() .second .push_back ( { cll_tag, { } } );
			std::vector < size_t > & local_vec =
				ver_of_elem .back() .second .back() .second;
			local_vec .reserve ( elem_size [ cll_type ] );
			for ( int jj = 0; jj < elem_size [ cll_type ]; jj ++ )
			{	size_t j;  ifstr >> j;  local_vec .push_back (j);  }
			std::getline ( ifstr, s );  assert ( trim_copy (s) == "");      }  }
	std::getline ( ifstr, s );  assert ( trim_copy (s) == "$EndElements");

}  // end of  import_msh_4_composed

//----------------------------------------------------------------------------------//


inline Cell retrieve_or_build_and_register  // hidden in anonymous namespace
( const int cll_type, const std::vector < size_t > & nodes,
  std::map < size_t, std::vector < std::pair < std::vector < size_t >, Cell > > > & built_cells,
  const std::map < size_t, Cell > & map_ver                                                     )
	
// since in gmsh the maximum topological dimension is 3,
// we do not keep 3D cells (they are built but not kept)

{	switch ( cll_type )
	{	case 1 :   // segment
		{	assert ( nodes .size() == 2 );
			size_t a = nodes [0], b = nodes [1];
			std::map < size_t, std::vector < std::pair < std::vector < size_t >, Cell > > >
				::const_iterator it = built_cells .find ( a );
			if ( it != built_cells .end() )
			{	std::vector < std::pair < std::vector < size_t >, Cell > >
					::const_iterator itt = it->second .begin();
				for ( ; itt != it->second .end(); itt++ )
				{	const std::vector < size_t > & other_nodes = itt->first;
					if ( other_nodes .size() != 1 ) continue;
					if ( other_nodes [0] == b )  // found it !
						return itt->second;                                    }  }
			it = built_cells .find ( b );
			if ( it != built_cells .end() )
			{	std::vector < std::pair < std::vector < size_t >, Cell > >
					::const_iterator itt = it->second .begin();
				for ( ; itt != it->second .end(); itt++ )
				{	const std::vector < size_t > & other_nodes = itt->first;
					if ( other_nodes .size() != 1 ) continue;
					if ( other_nodes [0] == a )  // found it ! reversed
						return itt->second .reverse();                         }  }
			// we must build a new segment
			// std::cout << "building segment " << a << " " << b << std::endl;
			std::map < size_t, Cell > ::const_iterator it_a = map_ver .find (a),
			                                           it_b = map_ver .find (b);
			assert ( it_a != map_ver .end() );  assert ( it_b != map_ver .end() );
			Cell seg ( tag::segment, it_a->second .reverse(), it_b->second );
			std::vector < std::pair < std::vector < size_t >, Cell > >
				& cells_with_a = built_cells [a];
			// in the above, if key 'a' does not exist yet in 'built_cells',
			// it will be inserted and an empty vector will be associated to it
			cells_with_a .push_back ( { { b }, seg } );
			return seg;                                                                      }
		case 2 :   // triangle
		{	assert ( nodes .size() == 3 );
			size_t a = nodes [0], b = nodes [1], c = nodes [2];
			std::map < size_t, std::vector < std::pair < std::vector < size_t >, Cell > > >
				::const_iterator it = built_cells .find ( a );
			if ( it != built_cells .end() )
			{	std::vector < std::pair < std::vector < size_t >, Cell > >
					::const_iterator itt = it->second .begin();
				for ( ; itt != it->second .end(); itt++ )
				{	const std::vector < size_t > & other_nodes = itt->first;
					if ( other_nodes .size() != 2 ) continue;
					if ( ( other_nodes [0] == b ) and ( other_nodes [1] == c ) )  // found it !
						return itt->second;
					if ( ( other_nodes [0] == c ) and ( other_nodes [1] == b ) )  // found it ! reversed
						return itt->second .reverse();                              }  }
			it = built_cells .find ( b );
			if ( it != built_cells .end() )
			{	std::vector < std::pair < std::vector < size_t >, Cell > >
					::const_iterator itt = it->second .begin();
				for ( ; itt != it->second .end(); itt++ )
				{	const std::vector < size_t > & other_nodes = itt->first;
					if ( other_nodes .size() != 2 ) continue;
					if ( ( other_nodes [0] == c ) and ( other_nodes [1] == a ) )  // found it !
						return itt->second;
					if ( ( other_nodes [0] == a ) and ( other_nodes [1] == c ) )  // found it ! reversed
						return itt->second .reverse();                              }  }
			it = built_cells .find ( c );
			if ( it != built_cells .end() )
			{	std::vector < std::pair < std::vector < size_t >, Cell > >
					::const_iterator itt = it->second .begin();
				for ( ; itt != it->second .end(); itt++ )
				{	const std::vector < size_t > & other_nodes = itt->first;
					if ( other_nodes .size() != 2 ) continue;
					if ( ( other_nodes [0] == a ) and ( other_nodes [1] == b ) )  // found it !
						return itt->second;
					if ( ( other_nodes [0] == b ) and ( other_nodes [1] == a ) )  // found it ! reversed
						return itt->second .reverse();                              }  }
			// we must build a new triangle
			// std::cout << "building triangle " << a << " " << b << " " << c << std::endl;
			Cell AB = retrieve_or_build_and_register ( 1, { a, b }, built_cells, map_ver );
			Cell BC = retrieve_or_build_and_register ( 1, { b, c }, built_cells, map_ver );
			Cell CA = retrieve_or_build_and_register ( 1, { c, a }, built_cells, map_ver );
			Cell tri ( tag::triangle, AB, BC, CA );
			std::vector < std::pair < std::vector < size_t >, Cell > >
				& cells_with_a = built_cells [a];
			// in the above, if key 'a' does not exist yet in 'built_cells',
			// it will be inserted and an empty vector will be associated to it
			cells_with_a .push_back ( { { b, c }, tri } );
			return tri;                                                                        }
		case 3 :   // quadrilateral
		// how do we distinguish a quadrilateral from a tetrahedron ?
		// a mesh should never have a square and a tetrahedron with the same nodes
		// and anyway we do not keep tetrahedra (3D cells are built but not kept)
		{	assert ( nodes .size() == 4 );
			size_t a = nodes [0], b = nodes [1], c = nodes [2], d = nodes [3];
			std::map < size_t, std::vector < std::pair < std::vector < size_t >, Cell > > >
				::const_iterator it = built_cells .find ( a );
			if ( it != built_cells .end() )
			{	std::vector < std::pair < std::vector < size_t >, Cell > >
					::const_iterator itt = it->second .begin();
				for ( ; itt != it->second .end(); itt++ )
				{	const std::vector < size_t > & other_nodes = itt->first;
					if ( other_nodes .size() != 3 ) continue;
					if ( ( other_nodes [0] == b ) and ( other_nodes [1] == c ) and
					     ( other_nodes [2] == d )                                  )  // found it !
						return itt->second;
					if ( ( other_nodes [0] == d ) and ( other_nodes [1] == c ) and
					     ( other_nodes [2] == b )                                  )  // found it ! reversed
						return itt->second .reverse();                                 }  }
			it = built_cells .find ( b );
			if ( it != built_cells .end() )
			{	std::vector < std::pair < std::vector < size_t >, Cell > >
					::const_iterator itt = it->second .begin();
				for ( ; itt != it->second .end(); itt++ )
				{	const std::vector < size_t > & other_nodes = itt->first;
					if ( other_nodes .size() != 3 ) continue;
					if ( ( other_nodes [0] == c ) and ( other_nodes [1] == d ) and
					     ( other_nodes [2] == a )                                  )  // found it !
						return itt->second;
					if ( ( other_nodes [0] == a ) and ( other_nodes [1] == d ) and
					     ( other_nodes [2] == c )                                  )  // found it ! reversed
						return itt->second .reverse();                                 }  }
			it = built_cells .find ( c );
			if ( it != built_cells .end() )
			{	std::vector < std::pair < std::vector < size_t >, Cell > >
					::const_iterator itt = it->second .begin();
				for ( ; itt != it->second .end(); itt++ )
				{	const std::vector < size_t > & other_nodes = itt->first;
					if ( other_nodes .size() != 3 ) continue;
					if ( ( other_nodes [0] == d ) and ( other_nodes [1] == a ) and
					     ( other_nodes [2] == b )                                  )  // found it !
						return itt->second;
					if ( ( other_nodes [0] == b ) and ( other_nodes [1] == a ) and
					     ( other_nodes [2] == d )                                  )  // found it ! reversed
						return itt->second .reverse();                                 }  }
			it = built_cells .find ( d );
			if ( it != built_cells .end() )
			{	std::vector < std::pair < std::vector < size_t >, Cell > >
					::const_iterator itt = it->second .begin();
				for ( ; itt != it->second .end(); itt++ )
				{	const std::vector < size_t > & other_nodes = itt->first;
					if ( other_nodes .size() != 3 ) continue;
					if ( ( other_nodes [0] == a ) and ( other_nodes [1] == b ) and
					     ( other_nodes [2] == c )                                  )  // found it !
						return itt->second;
					if ( ( other_nodes [0] == c ) and ( other_nodes [1] == b ) and
					     ( other_nodes [2] == a )                                  )  // found it ! reversed
						return itt->second .reverse();                                 }  }
			// we must build a new quadrangle
			// std::cout << "building quadrangle " << a << " " << b << " " << c << " " << d << std::endl;
			Cell AB = retrieve_or_build_and_register ( 1, { a, b }, built_cells, map_ver );
			Cell BC = retrieve_or_build_and_register ( 1, { b, c }, built_cells, map_ver );
			Cell CD = retrieve_or_build_and_register ( 1, { c, d }, built_cells, map_ver );
			Cell DA = retrieve_or_build_and_register ( 1, { d, a }, built_cells, map_ver );
			Cell quad ( tag::quadrangle, AB, BC, CD, DA );
			std::vector < std::pair < std::vector < size_t >, Cell > >
				& cells_with_a = built_cells [a];
			// in the above, if key 'a' does not exist yet in 'built_cells',
			// it will be inserted and an empty vector will be associated to it
			cells_with_a .push_back ( { { b, c, d }, quad } );
			return quad;                                                                        }
		case 4 :   // tetrahedron
		// how do we distinguish a quadrilateral from a tetrahedron ?
		// a mesh should never have a square and a tetrahedron with the same nodes
		// and anyway we do not keep tetrahedra (3D cells are built but not kept)
		{	assert ( nodes .size() == 4 );
			size_t a = nodes [0], b = nodes [1], c = nodes [2], d = nodes [3];
			// since in gmsh the maximum topological dimension is 3,
			// there is no use to keep 3D cells in built_cells -- we simply build a new tetrahedron
			// std::cout << "building tetrahedron " << a << " " << b << " " << c << " " << d << std::endl;
			Cell ABC = retrieve_or_build_and_register ( 2, { a, b, c }, built_cells, map_ver );
			Cell ADB = retrieve_or_build_and_register ( 2, { a, d, b }, built_cells, map_ver );
			Cell BDC = retrieve_or_build_and_register ( 2, { b, d, c }, built_cells, map_ver );
			Cell ACD = retrieve_or_build_and_register ( 2, { a, c, d }, built_cells, map_ver );
			Cell tetra ( tag::tetrahedron, ABC, ADB, BDC, ACD );
			return tetra;                                                                        }
		case 5 :   // hexahedron, parallelipiped, cube
		{	assert ( nodes .size() == 8 );
			size_t a = nodes [0], b = nodes [1], c = nodes [2], d = nodes [3],
			       e = nodes [4], f = nodes [5], g = nodes [6], h = nodes [7];
			// since in gmsh the maximum topological dimension is 3,
			// there is no use to keep 3D cells in built_cells -- we simply build a new hexahedron
			// std::cout << "building cube " << a << " " << b << " " << c << " " << d << " "
			//           << e << " " << f << " " << g << " " << h << std::endl;
			Cell ABCD = retrieve_or_build_and_register ( 3, { a, b, c, d }, built_cells, map_ver );
			Cell HGFE = retrieve_or_build_and_register ( 3, { h, g, f, e }, built_cells, map_ver );
			Cell BAEF = retrieve_or_build_and_register ( 3, { b, a, e, f }, built_cells, map_ver );
			Cell ADHE = retrieve_or_build_and_register ( 3, { a, d, h, e }, built_cells, map_ver );
			Cell DCGH = retrieve_or_build_and_register ( 3, { d, c, g, h }, built_cells, map_ver );
			Cell CBFG = retrieve_or_build_and_register ( 3, { c, b, f, g }, built_cells, map_ver );
			Cell cube ( tag::hexahedron, ABCD, HGFE, BAEF, ADHE, DCGH, CBFG );
			return cube;                                                                            }
		default :  assert ( false );
	}  // end of switch ( cll_type );

	assert ( false );
	return Cell ( tag::non_existent );  // just to avoid compilation errors
		
}	 // end of retrieve_or_build_and_register
	
//----------------------------------------------------------------------------------//


inline void import_msh_common  // hidden in anonymous namespace
( Mesh * that, std::vector < std::pair < std::pair < int, int >,
     std::vector < std::pair < size_t, std::vector < size_t > > > > > & ver_of_elem,
  std::map < size_t, Cell > & map_ver, size_t global_top_dim                        )

// ver_of_elem == [ { { ent_tag, geom_type }, [ { cll_tag, [ node_tag ] } ] } ]
// we do not keep the dimension, it is implied by the geometry type

{	// vertices have already been built, kept in map_ver
	// here we build cells of dimension >= 1
	std::map < size_t, std::vector < std::pair < std::vector < size_t >, Cell > > > built_cells;
	// if a cell, say abcd, has already been built, is gets registered in the above
	// under built_cells [ tag_a, { { tag_b, tag_c, tag_d }, ABCD } ]

	const std::vector < size_t > elem_dim { 0, 1, 2, 2, 3, 3, 3, 3 };

	Mesh result ( tag::fuzzy, tag::of_dim, global_top_dim );

	for ( size_t ent = 0; ent < ver_of_elem .size(); ent ++ )
	{	size_t geom_type = ver_of_elem [ent] .first .second;
		if ( geom_type == 15 )  continue;  // vertex
		assert ( ( 1 <= geom_type ) and ( geom_type <= 7 ) );
		std::vector < std::pair < size_t, std::vector < size_t > > > & elem =
			ver_of_elem [ent] .second;
		for ( size_t el = 0; el < elem .size(); el++ )
		{	std::vector < size_t > & nodes = elem [el] .second;
			Cell cll = retrieve_or_build_and_register ( geom_type, nodes, built_cells, map_ver );
			assert ( cll .dim() == elem_dim [ geom_type ] );
			if ( cll .dim() == global_top_dim )  cll .add_to_mesh ( result );                     }  }

	// after all cells around a vertex or segment or face have been built,
	// that vertex or segment or face can be eliminated from the database
	// thus keeping its size down

	*that = result;

}  // end of  import_msh_common

//----------------------------------------------------------------------------------//


inline void import_msh ( Mesh * that, const std::string filename )
// hidden in anonymous namespace

{	double msh_version;
	int geom_dim;

	// we perform a preliminary passage and find the version and
	// the geometric dimension
	{ // just a block of code for hiding names
	std::ifstream ifstr ( filename );
	std::string s;
	std::getline ( ifstr, s );  assert ( trim_copy (s) == "$MeshFormat");
	ifstr >> msh_version;
	int binary_file;  ifstr >> binary_file;
	if ( binary_file != 0 )
	{	std::cout << "binary files not accepted" <<std::endl;
		exit (1);                                             }
	std::getline ( ifstr, s );  // one more number and the endline
	std::getline ( ifstr, s );  assert ( trim_copy (s) == "$EndMeshFormat");

	while ( trim_copy (s) != "$Nodes")  std::getline ( ifstr, s );

	size_t numEntityBlocks, minNodeTag, maxNodeTag;
	bool all_z_zero = true;
	size_t nb_nodes;
	if ( msh_version < 4. )
	{	ifstr >> nb_nodes;  maxNodeTag = nb_nodes;
		std::getline ( ifstr, s );  assert ( trim_copy (s) == "");  }
	else  // msh_version >= 4.
	{	ifstr >> numEntityBlocks >> nb_nodes >> minNodeTag >> maxNodeTag;
		std::getline ( ifstr, s );  assert ( trim_copy (s) == "" );       }

	if ( msh_version < 4. )
	{	for ( size_t i = 0; i < nb_nodes; i++ )
		{	size_t j;  ifstr >> j;  assert ( j == i+1 );
			double xx, yy, zz;  ifstr >> xx >> yy >> zz;
			if ( zz != 0. )  all_z_zero = false;
			std::getline ( ifstr, s );  assert ( trim_copy (s) == "");   }
		std::getline ( ifstr, s );  assert ( trim_copy (s) == "$EndNodes");  }
	else  // msh_version >= 4.
	{	int top_dim, entity_tag, parametric;
		size_t nbnodes;
		for ( size_t ent = 0; ent < numEntityBlocks; ent++ )
		{	ifstr >> top_dim >> entity_tag >> parametric >> nbnodes;
			std::getline ( ifstr, s );  assert ( trim_copy (s) == "" );
			for ( size_t i = 0; i < nbnodes; i++ )
			{	size_t j;  ifstr >> j;     // node tag
				std::getline ( ifstr, s ); assert ( trim_copy (s) == "" );
				assert ( j <= maxNodeTag );                                 }
			for ( size_t i = 0; i < nbnodes; i++ )
			{	double xx, yy, zz;  ifstr >> xx >> yy >> zz;
				std::getline ( ifstr, s );  assert ( trim_copy (s) == "" );
				if ( zz != 0. ) all_z_zero = false;                                }  }
		std::getline ( ifstr, s );  assert ( trim_copy (s) == "$EndNodes");         }

	if ( all_z_zero )  geom_dim = 2;
	else  geom_dim = 3;
	} // just a block of code for hiding 'ifstr'
	
	// std::cout << "msh version " << msh_version << ", geom dim " << geom_dim << std::endl;
	
  std::vector < std::pair < std::pair < int, int >,
	    std::vector < std::pair < size_t, std::vector < size_t > > > > > ver_of_elem;
	// ver_of_elem == [ { { ent_tag, geom_type }, [ { cll_tag, [ node_tag ] } ] } ]
	// we do not keep the dimension, it is implied by the geometry type

	std::map < size_t, Cell > map_ver;

	size_t global_top_dim;
	if ( msh_version < 4. )
		global_top_dim = import_msh_2 ( filename, geom_dim, ver_of_elem, map_ver );
	else  global_top_dim = import_msh_4_one_mesh ( filename, geom_dim, ver_of_elem, map_ver );

	import_msh_common ( that, ver_of_elem, map_ver, global_top_dim );

}  // end of  import_msh

//----------------------------------------------------------------------------------//


inline std::map < std::pair < int, int >, Mesh > import_msh_common  // hidden in anonymous namespace
( std::vector < std::pair < std::pair < int, int >,
     std::vector < std::pair < size_t, std::vector < size_t > > > > > & ver_of_elem,
  std::map < size_t, Cell > & map_ver, size_t global_top_dim                        )

// ver_of_elem == [ { { ent_tag, geom_type }, [ { cll_tag, [ node_tag ] } ] } ]
// we do not keep the dimension, it is implied by the geometry type

{	// vertices have already been built, kept in map_ver
	// here we build cells of dimension >= 1
	std::map < size_t, std::vector < std::pair < std::vector < size_t >, Cell > > > built_cells;
	// if a cell, say abcd, has already been built, is gets registered in the above
	// under built_cells [ tag_a, { { tag_b, tag_c, tag_d }, ABCD } ]

	const std::vector < size_t > elem_dim { 0, 1, 2, 2, 3, 3, 3, 3 };

	std::map < std::pair < int, int >, Mesh > meshes;
	
	for ( size_t ent = 0; ent < ver_of_elem .size(); ent ++ )
	{	int ent_tag = ver_of_elem [ent] .first .first;
		int geom_type = ver_of_elem [ent] .first .second;
		if ( geom_type == 0 )  continue;  // vertex (changed from 15 to 0)
		assert ( ( 1 <= geom_type ) and ( geom_type <= 7 ) );
		if ( global_top_dim > 0 )  assert ( elem_dim [ geom_type ] == global_top_dim );
		Mesh result ( tag::fuzzy, tag::of_dim, elem_dim [ geom_type ] );
		std::vector < std::pair < size_t, std::vector < size_t > > > & elem =
			ver_of_elem [ent] .second;
		for ( size_t el = 0; el < elem .size(); el++ )
		{	std::vector < size_t > & nodes = elem [el] .second;
			Cell cll = retrieve_or_build_and_register ( geom_type, nodes, built_cells, map_ver );
			assert ( cll .dim() == elem_dim [ geom_type ] );
		  cll .add_to_mesh ( result );                                                         }
		// if elem_dim [ geom_type ] == 1  we should try to build a Connected Mesh instead
		std::cout << "global.cpp line 4287 " << ent << " " << ver_of_elem [ent] .first .first
							<< " " << ver_of_elem [ent] .first .second << std::endl;
		std::pair < std::map < std::pair < int, int >, Mesh > ::iterator, bool > p =
			meshes .insert ( { { elem_dim [ geom_type ], ent_tag }, result } );
		assert ( p .second );
	}

	// after all cells around a vertex or segment or face have been built,
	// that vertex or segment or face can be eliminated from the database
	// thus keeping its size down

	return meshes;

}  // end of  import_msh_common

} // end of anonymous namespace

//----------------------------------------------------------------------------------//


std::map < std::pair < int, int >, Mesh > import_msh ( const std::string filename )
// hidden in anonymous namespace

{	double msh_version;
	int geom_dim;

	// we perform a preliminary passage and find the version and
	// the geometric dimension
	{ // just a block of code for hiding names
	std::ifstream ifstr ( filename );
	std::string s;
	std::getline ( ifstr, s );  assert ( trim_copy (s) == "$MeshFormat");
	ifstr >> msh_version;
	int binary_file;  ifstr >> binary_file;
	if ( binary_file != 0 )
	{	std::cout << "binary files not accepted" <<std::endl;
		exit (1);                                             }
	std::getline ( ifstr, s );  // one more number and the endline
	std::getline ( ifstr, s );  assert ( trim_copy (s) == "$EndMeshFormat");

	while ( trim_copy (s) != "$Nodes")  std::getline ( ifstr, s );

	size_t numEntityBlocks, minNodeTag, maxNodeTag;
	bool all_z_zero = true;
	size_t nb_nodes;
	if ( msh_version < 4. )
	{	ifstr >> nb_nodes;  maxNodeTag = nb_nodes;
		std::getline ( ifstr, s );  assert ( trim_copy (s) == "");  }
	else  // msh_version >= 4.
	{	ifstr >> numEntityBlocks >> nb_nodes >> minNodeTag >> maxNodeTag;
		std::getline ( ifstr, s );  assert ( trim_copy (s) == "" );       }

	if ( msh_version < 4. )
	{	for ( size_t i = 0; i < nb_nodes; i++ )
		{	size_t j;  ifstr >> j;  assert ( j == i+1 );
			double xx, yy, zz;  ifstr >> xx >> yy >> zz;
			if ( zz != 0. )  all_z_zero = false;
			std::getline ( ifstr, s );  assert ( trim_copy (s) == "");   }
		std::getline ( ifstr, s );  assert ( trim_copy (s) == "$EndNodes");  }
	else  // msh_version >= 4.
	{	int top_dim, entity_tag, parametric;
		size_t nbnodes;
		for ( size_t ent = 0; ent < numEntityBlocks; ent++ )
		{	ifstr >> top_dim >> entity_tag >> parametric >> nbnodes;
			std::getline ( ifstr, s );  assert ( trim_copy (s) == "" );
			for ( size_t i = 0; i < nbnodes; i++ )
			{	size_t j;  ifstr >> j;     // node tag
				std::getline ( ifstr, s ); assert ( trim_copy (s) == "" );
				assert ( j <= maxNodeTag );                                 }
			for ( size_t i = 0; i < nbnodes; i++ )
			{	double xx, yy, zz;  ifstr >> xx >> yy >> zz;
				std::getline ( ifstr, s );  assert ( trim_copy (s) == "" );
				if ( zz != 0. ) all_z_zero = false;                                }  }
		std::getline ( ifstr, s );  assert ( trim_copy (s) == "$EndNodes");         }

	if ( all_z_zero )  geom_dim = 2;
	else  geom_dim = 3;
	} // just a block of code for hiding 'ifstr'
	
	// std::cout << "msh version " << msh_version << ", geom dim " << geom_dim << std::endl;
	
  std::vector < std::pair < std::pair < int, int >,
	    std::vector < std::pair < size_t, std::vector < size_t > > > > > ver_of_elem;
	// ver_of_elem == [ { { ent_tag, geom_type }, [ { cll_tag, [ node_tag ] } ] } ]
	// we do not keep the dimension, it is implied by the geometry type

	std::map < size_t, Cell > map_ver;

	size_t global_top_dim = 0;
	if ( msh_version < 4. )
		global_top_dim = import_msh_2 ( filename, geom_dim, ver_of_elem, map_ver );
	else  import_msh_4_composed ( filename, geom_dim, ver_of_elem, map_ver );

	return import_msh_common ( ver_of_elem, map_ver, global_top_dim );

}  // end of  import_msh

//----------------------------------------------------------------------------------//



Mesh::Mesh ( const tag::Import &, const tag::Msh &, const std::string filename )
:	Mesh ( tag::non_existent )
{	import_msh ( this, filename );  }

//----------------------------------------------------------------------------------//


namespace { // anonymous namespace, mimics static linkage

Mesh fold_common ( const Mesh & msh, const std::map < Cell, Cell > & corresp_seg,
                     size_t dim, bool keep_map, std::map < Cell, Cell > & m      )
// hidden in anonymous namespace

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
	               tag::segments, tag::one_dummy_wrapper      ),
	         tag::freshly_created                               ),
				tag::one_dummy_wrapper                                  );
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

	
Mesh fold_common_no_sides   // hidden in anonymous namespace
( const Mesh & msh,
	const std::map < Cell, std::pair < Cell, Manifold::Action > > & corresp_ver,
  size_t dim, bool keep_map, std::map < Cell, Cell > & m                      )

// build 'corresp_seg' then call fold_common
	
{	if ( msh .dim() == 1 )
	{	Mesh result ( tag::fuzzy, tag::of_dim, 1 );
		Mesh::Iterator it_seg = msh .iterator ( tag::over_segments );
		for ( it_seg .reset(); it_seg .in_range(); it_seg++ )
		{	Cell seg = *it_seg;
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::const_iterator it_base_rev = corresp_ver.find ( seg .base() .reverse() );
			assert ( it_base_rev != corresp_ver .end() );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::const_iterator it_tip = corresp_ver .find ( seg .tip() );
			assert ( it_tip != corresp_ver .end() );
			Cell new_seg ( tag::segment,
			               it_base_rev->second .first .reverse(), it_tip->second .first );
			new_seg .winding() = it_tip->second .second - it_base_rev->second .second;
			// new_seg.winding() = corresp_ver [ seg.tip()            ] .second -
			//	                corresp_ver [ seg.base().reverse() ] .second  ;
			new_seg.add_to_mesh ( result );                                            }
			
		return result;                                                                 }

	else
	{	assert ( msh.dim() == 2 );
		// we use a map -- for a faster code, we could use Cell::Core::hook
		std::map < Cell, Cell > corresp_seg;

		Mesh::Iterator it_seg = msh .iterator ( tag::over_segments );
		for ( it_seg .reset(); it_seg .in_range(); it_seg++ )
		{	Cell seg = *it_seg;
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::const_iterator it_base_rev = corresp_ver .find ( seg .base() .reverse() );
			assert ( it_base_rev != corresp_ver .end() );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::const_iterator it_tip = corresp_ver .find ( seg .tip() );
			assert ( it_tip != corresp_ver .end() );
			Cell new_seg ( tag::segment,
			               it_base_rev->second .first .reverse(), it_tip->second .first );
			new_seg.winding() = it_tip->second .second - it_base_rev->second .second;
			// new_seg.winding() = corresp_ver [ seg.tip()            ] .second -
			//	                corresp_ver [ seg.base().reverse() ] .second  ;
			// inspired in item 24 of the book : Scott Meyers, Effective STL
			std::map < Cell, Cell > ::iterator it_map = corresp_seg .lower_bound ( seg );
			assert ( ( it_map == corresp_seg .end() ) or
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
// hidden in anonymous namespace

// take a mesh and fold it around the current working manifold,
// which must be a quotient manifold

{	Manifold space = Manifold::working;
	assert ( space.exists() );  // we use the current (quotient) manifold
	Manifold::Quotient * manif_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space .core );
	assert ( manif_q );
	Function coords_q = space .coordinates();
	Manifold mani_Eu = manif_q->base_space;  // underlying Euclidian manifold
	Function coords_Eu = mani_Eu .coordinates();

	// we use a map -- for a faster code, we could use Cell::Core::hook
	std::map < Cell, std::pair < Cell, Manifold::Action > > corresp_ver;

	Mesh::Iterator it_ver = that->iterator ( tag::over_vertices );
	for ( it_ver .reset(); it_ver .in_range(); it_ver++ )
	{	Cell V = *it_ver;
		Cell new_V ( tag::vertex );
		coords_Eu ( new_V ) = coords_Eu (V);
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_map = corresp_ver .lower_bound ( V );
		assert ( ( it_map == corresp_ver .end() ) or
		         ( corresp_ver.key_comp()(V,it_map->first) ) );
		corresp_ver .emplace_hint ( it_map, std::piecewise_construct,
		    std::forward_as_tuple ( V ), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { new_V, 0 } ) );  }
		// corresp_ver [ V ] = { new_V, 0 };

	return fold_common_no_sides ( *that, corresp_ver, dim, keep_map, m );     }
	
//----------------------------------------------------------------------------------//

	
Mesh fold_no_sides ( Mesh * that, const tag::UseExistingVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
									size_t dim, bool keep_map, std::map < Cell, Cell > & m )
// hidden in anonymous namespace

// take a mesh and fold it around the current working manifold,
// which must be a quotient manifold

{	Manifold space = Manifold::working;
	assert ( space.exists() );  // we use the current (quotient) manifold
	Manifold::Quotient * manif_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space .core );
	assert ( manif_q );
	Function coords_q = space .coordinates();
	Manifold mani_Eu = manif_q->base_space;  // underlying Euclidian manifold

	// we use a map -- for a faster code, we could use Cell::Core::hook
	std::map < Cell, std::pair < Cell, Manifold::Action > > corresp_ver;

	Mesh::Iterator it_ver = that->iterator ( tag::over_vertices );
	for ( it_ver .reset(); it_ver .in_range(); it_ver++ )
	{	Cell V = *it_ver;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_map = corresp_ver .lower_bound ( V );
		assert ( ( it_map == corresp_ver .end() ) or
		         ( corresp_ver.key_comp()(V,it_map->first) ) );
		corresp_ver.emplace_hint ( it_map, std::piecewise_construct,
		    std::forward_as_tuple ( V ), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { V, 0 } ) );  }
		// corresp_ver [ V ] = { V, 0 };

	return fold_common_no_sides ( *that, corresp_ver, dim, keep_map, m );               }
	
//----------------------------------------------------------------------------------//

	
Mesh fold_common_two_sides   // hidden in anonymous namespace
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
	for ( it_seg .reset(); it_seg .in_range(); it_seg++ )
	{	Cell seg = *it_seg;
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::const_iterator it_base_rev = corresp_ver .find ( seg .base() .reverse() );
		assert ( it_base_rev != corresp_ver .end() );
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::const_iterator it_tip = corresp_ver .find ( seg .tip() );
		assert ( it_tip != corresp_ver .end() );
		if ( seg .belongs_to ( side_2, tag::not_oriented ) ) continue;
			// use corresponding segment on side_1
		Cell new_seg ( tag::segment,
		               it_base_rev->second .first .reverse(), it_tip->second .first );
		new_seg .winding() = it_tip->second .second - it_base_rev->second .second;
		// new_seg .winding() = corresp_ver [ seg .tip()             ] .second -
		//	                    corresp_ver [ seg .base() .reverse() ] .second  ;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, Cell > ::iterator it_map = corresp_seg.lower_bound ( seg );
		assert ( ( it_map == corresp_seg .end() ) or
		         ( corresp_seg.key_comp()(seg,it_map->first) ) );
		corresp_seg.emplace_hint ( it_map, std::piecewise_construct,
		  std::forward_as_tuple ( seg ), std::forward_as_tuple ( new_seg ) );     }
		// corresp_seg [ seg ] = new_seg;               

	Mesh::Iterator it_seg_1 = side_1 .iterator
		( tag::over_segments, tag::require_order, tag::force_positive );
	Mesh::Iterator it_seg_2 = side_2 .iterator
		( tag::over_segments, tag::require_order, tag::force_positive );
	for ( it_seg_1 .reset(), it_seg_2 .reset(); it_seg_1 .in_range(); it_seg_1++, it_seg_2++ )
	{	assert ( it_seg_2 .in_range() );
		Cell seg_1 = *it_seg_1;
		Cell seg_2 = *it_seg_2;
		std::map < Cell, Cell > ::iterator it = corresp_seg .find ( seg_1 );
		assert ( it != corresp_seg .end() );
		Cell new_seg_1 = it->second .reverse();
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, Cell > ::iterator it_map = corresp_seg .lower_bound ( seg_2 );
		assert ( ( it_map == corresp_seg .end() ) or
		         ( corresp_seg.key_comp()(seg_2,it_map->first) ) );
		corresp_seg.emplace_hint ( it_map, std::piecewise_construct,
		  std::forward_as_tuple ( seg_2 ), std::forward_as_tuple ( new_seg_1 ) );     }
		// corresp_seg [ seg_2 ] = corresp_seg [ seg_1 ];               
	assert ( not it_seg_2 .in_range() );
		
	return fold_common ( msh, corresp_seg, dim, keep_map, m );

}  // end of  fold_common_two_sides
	
//----------------------------------------------------------------------------------//


Mesh fold_two_sides    // hidden in anonymous namespace    // line 1283
( Mesh * that, const tag::Identify &, const Mesh & side_1,
               const tag::With &,     const Mesh & side_2,
  const tag::BuildNewVertices &,
  const tag::ReturnMapBetween &, const tag::CellsOfDim &,
  size_t dim, bool keep_map, std::map < Cell, Cell > & m   )

// take a mesh whose external boundary has two parallel segments
// and identify these two segments

// a quotient manifold will be built with one action generator
	
{	// we use the current (Euclidian) manifold
	Manifold space = Manifold::working;
	assert ( space .exists() );
	Function coord = space .coordinates();
	assert ( coord .nb_of_components() == 2 );
	// we split 'coord' into its components
	Function x = coord [0], y = coord [1];

	// first we need to identify a translation which moves side_1 into side_2
	// if these two 'sides' are Connected::OneDim, we may use first and last vertices
	// but they may be (one-dim or higher-dim) Fuzzy meshes ... too complicated !
	// almost impossible to find a mapping between vertices of, say, two opposite square faces

	Cell A = side_1 .first_vertex() .reverse();
	Cell B = side_1 .last_vertex();
	Cell C = side_2 .last_vertex();
	Cell D = side_2 .first_vertex() .reverse();

	double dx = x(D) - x(A), dy = y(D) - y(A);
	double norm = std::sqrt ( dx*dx + dy*dy );
	assert ( std::abs ( dx - ( x(C) - x(B) ) ) < 1.e-4 * norm );
	assert ( std::abs ( dy - ( y(C) - y(B) ) ) < 1.e-4 * norm );

	// the desired translation is ( dx, dy )
	Manifold::Action g ( tag::transforms, coord, tag::into, (x+dx) && ( y+dy) );
	Manifold manif_q = space .quotient ( g );

	// we use a map -- for a faster code, we could use Cell::Core::hook
	std::map < Cell, std::pair < Cell, Manifold::Action > > corresp_ver;

	Mesh::Iterator it_ver = that->iterator ( tag::over_vertices );
	for ( it_ver .reset(); it_ver .in_range(); it_ver++ )
	{	Cell V = *it_ver;
		if ( V.belongs_to ( side_2 ) ) continue;
		Cell new_V ( tag::vertex );
		x ( new_V ) = x (V);   y ( new_V ) = y (V);
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_map = corresp_ver .lower_bound ( V );
		assert ( ( it_map == corresp_ver .end() ) or
		         ( corresp_ver.key_comp()(V,it_map->first) ) );
		corresp_ver.emplace_hint ( it_map, std::piecewise_construct,
		    std::forward_as_tuple ( V ), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { new_V, 0 } ) );  }
		// corresp_ver [ V ] = { new_V, 0 };

	assert ( side_1 .number_of ( tag::segments ) == side_2 .number_of ( tag::segments ) );
	Mesh::Iterator it1 = side_1 .iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it2 = side_2 .iterator ( tag::over_vertices, tag::require_order );
	for ( it1 .reset(), it2 .reset(); it1 .in_range(); it1++, it2++ )
	{	assert ( it2 .in_range() );
		Cell V = *it1;  Cell W = *it2;
		assert ( std::abs ( dx - ( x(W) - x(V) ) ) < 1.e-4 * norm );
		assert ( std::abs ( dy - ( y(W) - y(V) ) ) < 1.e-4 * norm );
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_V = corresp_ver .find ( V );
		assert ( it_V != corresp_ver .end() );
		assert ( it_V->second .second == 0 );
		Cell new_V = it_V->second .first;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_W = corresp_ver .lower_bound ( W );
		assert ( ( it_W == corresp_ver .end() ) or
		         ( corresp_ver.key_comp()(W,it_W->first) ) );
		corresp_ver.emplace_hint ( it_W, std::piecewise_construct,
		    std::forward_as_tuple (W), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { new_V, g } ) );  }
		// corresp_ver [ W ] = { new_V, g };
	assert ( not it2 .in_range() );
						 
	return fold_common_two_sides ( *that, corresp_ver, side_1, side_2, dim, keep_map, m );

}  // end of fold_two_sides with tag::build_new_vertices

	
//----------------------------------------------------------------------------------//

	
Mesh fold_two_sides    // hidden in anonymous namespace    // line 1373
( Mesh * that, const tag::Identify &, const Mesh & side_1,
               const tag::With &,     const Mesh & side_2,
  const tag::UseExistingVertices &,
  const tag::ReturnMapBetween &, const tag::CellsOfDim &,
  size_t dim, bool keep_map, std::map < Cell, Cell > & m  )

// take a mesh whose external boundary has two parallel faces
// (which may be segments or higher-dimensional meshes)
// and identify these two faces
// they do not touch each other

// a quotient manifold will be built with one action generator
	
{	// we use the current (Euclidian) manifold
	Manifold space = Manifold::working;
	assert ( space .exists() );
	Function coord = space .coordinates();
	assert ( coord .nb_of_components() == 2 );
	// we split 'coord' into its components
	Function x = coord [0], y = coord [1];

	// first we need to identify a translation which moves side_1 into side_2
	// if these two 'sides' are Connected::OneDim, we may use first and last vertices
	// but they may be (one-dim or higher-dim) Fuzzy meshes ... too complicated !
	// almost impossible to find a mapping between vertices of, say, two opposite square faces
	
	Cell A = side_1 .first_vertex() .reverse();
	Cell B = side_1 .last_vertex();
	Cell C = side_2 .last_vertex();
	Cell D = side_2 .first_vertex() .reverse();

	double dx = x(D) - x(A), dy = y(D) - y(A);
	double norm = std::sqrt ( dx*dx + dy*dy );
	assert ( std::abs ( dx - ( x(C) - x(B) ) ) < 1.e-4 * norm );
	assert ( std::abs ( dy - ( y(C) - y(B) ) ) < 1.e-4 * norm );

	// the desired translation is ( dx, dy )
	Manifold::Action g ( tag::transforms, coord, tag::into, (x+dx) && ( y+dy) );
	Manifold manif_q = space .quotient ( g );

	// we use a map -- for a faster code, we could use Cell::Core::hook
	std::map < Cell, std::pair < Cell, Manifold::Action > > corresp_ver;

	Mesh::Iterator it_ver = that->iterator ( tag::over_vertices );
	for ( it_ver .reset(); it_ver .in_range(); it_ver++ )
	{	Cell V = *it_ver;
		if ( V.belongs_to ( side_2 ) ) continue;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_map = corresp_ver .lower_bound (V);
		assert ( ( it_map == corresp_ver .end() ) or
		         ( corresp_ver.key_comp()(V,it_map->first) ) );
		corresp_ver.emplace_hint ( it_map, std::piecewise_construct,
		    std::forward_as_tuple (V), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { V, 0 } ) );  }
		// corresp_ver [V] = { V, 0 };

	assert ( side_1 .number_of ( tag::segments ) == side_2 .number_of ( tag::segments ) );
	Mesh::Iterator it1 = side_1 .iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it2 = side_2 .iterator ( tag::over_vertices, tag::require_order );
	for ( it1 .reset(), it2 .reset(); it1 .in_range(); it1++, it2++ )
	{	assert ( it2 .in_range() );
		Cell V = *it1;  Cell W = *it2;
		assert ( std::abs ( dx - ( x(W) - x(V) ) ) < 1.e-4 * norm );
		assert ( std::abs ( dy - ( y(W) - y(V) ) ) < 1.e-4 * norm );
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_W = corresp_ver .lower_bound (W);
		assert ( ( it_W == corresp_ver .end() ) or
		         ( corresp_ver.key_comp()(W,it_W->first) ) );
		corresp_ver.emplace_hint ( it_W, std::piecewise_construct,
		    std::forward_as_tuple (W), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { V, g } ) );  }
		// corresp_ver [W] = { V, g };
	assert ( not it2 .in_range() );
						 
	return fold_common_two_sides ( *that, corresp_ver, side_1, side_2, dim, keep_map, m );

}  // end of fold_two_sides with tag::use_existing_vertices

	
//----------------------------------------------------------------------------------//

	
Mesh fold_common_four_sides   // hidden in anonymous namespace
( const Mesh & msh,
  const std::map < Cell, std::pair < Cell, Manifold::Action > > & corresp_ver,
  const Mesh & side_1, const Mesh & side_2, const Mesh & side_3, const Mesh & side_4,
  size_t dim, bool keep_map, std::map < Cell, Cell > & m                             )

// build 'corresp_seg', taking care to identify side_1 with side_2 and side_3 with side_4
// then call fold_common
	
{	assert ( msh .dim() == 2 );
	// we use a map -- for a faster code, we could use Cell::Core::hook
	std::map < Cell, Cell > corresp_seg;

	Mesh::Iterator it_seg = msh .iterator ( tag::over_segments );
	for ( it_seg .reset(); it_seg .in_range(); it_seg++ )
	{	Cell seg = *it_seg;
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::const_iterator it_base_rev = corresp_ver .find ( seg .base() .reverse() );
		assert ( it_base_rev != corresp_ver .end() );
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::const_iterator it_tip = corresp_ver .find ( seg .tip() );
		assert ( it_tip != corresp_ver .end() );
		if ( seg .belongs_to ( side_2, tag::not_oriented ) ) continue;
			// use corresponding segment on side_1
		if ( seg .belongs_to ( side_4, tag::not_oriented ) ) continue;
			// use corresponding segment on side_3
		Cell new_seg ( tag::segment,
		               it_base_rev->second .first .reverse(), it_tip->second .first );
		new_seg .winding() = it_tip->second .second - it_base_rev->second .second;
		// new_seg.winding() = corresp_ver [ seg .tip()             ] .second -
		//	                   corresp_ver [ seg .base() .reverse() ] .second  ;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, Cell > ::iterator it_map = corresp_seg.lower_bound ( seg );
		assert ( ( it_map == corresp_seg .end() ) or
		         ( corresp_seg.key_comp()(seg,it_map->first) ) );
		corresp_seg.emplace_hint ( it_map, std::piecewise_construct,
		  std::forward_as_tuple ( seg ), std::forward_as_tuple ( new_seg ) );     }
		// corresp_seg [ seg ] = new_seg;               

	Mesh::Iterator it_seg_1 = side_1 .iterator
		( tag::over_segments, tag::require_order, tag::force_positive );
	Mesh::Iterator it_seg_2 = side_2 .iterator
		( tag::over_segments, tag::require_order, tag::force_positive );
	for ( it_seg_1 .reset(), it_seg_2 .reset(); it_seg_1 .in_range(); it_seg_1++, it_seg_2++ )
	{	assert ( it_seg_2 .in_range() );
		Cell seg_1 = *it_seg_1;
		Cell seg_2 = *it_seg_2;
		std::map < Cell, Cell > ::iterator it = corresp_seg .find ( seg_1 );
		assert ( it != corresp_seg .end() );
		Cell new_seg_1 = it->second .reverse();
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, Cell > ::iterator it_map = corresp_seg .lower_bound ( seg_2 );
		assert ( ( it_map == corresp_seg .end() ) or
		         ( corresp_seg.key_comp()(seg_2,it_map->first) ) );
		corresp_seg .emplace_hint ( it_map, std::piecewise_construct,
		  std::forward_as_tuple ( seg_2 ), std::forward_as_tuple ( new_seg_1 ) );     }
		// corresp_seg [ seg_2 ] = corresp_seg [ seg_1 ];               
	assert ( not it_seg_2 .in_range() );
		
	Mesh::Iterator it_seg_3 = side_3 .iterator
		( tag::over_segments, tag::require_order, tag::force_positive );
	Mesh::Iterator it_seg_4 = side_4 .iterator
		( tag::over_segments, tag::require_order, tag::force_positive );
	for ( it_seg_3 .reset(), it_seg_4 .reset(); it_seg_3 .in_range(); it_seg_3++, it_seg_4++ )
	{	assert ( it_seg_4 .in_range() );
		Cell seg_3 = *it_seg_3;
		Cell seg_4 = *it_seg_4;
		std::map < Cell, Cell > ::iterator it = corresp_seg .find ( seg_3 );
		assert ( it != corresp_seg .end() );
		Cell new_seg_3 = it->second .reverse();
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, Cell > ::iterator it_map = corresp_seg .lower_bound ( seg_4 );
		assert ( ( it_map == corresp_seg .end() ) or
		         ( corresp_seg.key_comp()(seg_4,it_map->first) ) );
		corresp_seg .emplace_hint ( it_map, std::piecewise_construct,
		  std::forward_as_tuple ( seg_4 ), std::forward_as_tuple ( new_seg_3 ) );     }
		// corresp_seg [ seg_4 ] = corresp_seg [ seg_3 ];               
	assert ( not it_seg_4 .in_range() );
		
	return fold_common ( msh, corresp_seg, dim, keep_map, m );

}  // end of  fold_common_four_sides
	
//----------------------------------------------------------------------------------//


Mesh fold_four_sides   // hidden in anonymous namespace    // line 1544
( Mesh * that, const tag::Identify &, const Mesh & side_1,
               const tag::With &,     const Mesh & side_2,
               const tag::Identify &, const Mesh & side_3,
               const tag::With &,     const Mesh & side_4,
  const tag::BuildNewVertices &,
  const tag::ReturnMapBetween &, const tag::CellsOfDim &,
  size_t dim, bool keep_map, std::map < Cell, Cell > & m  )

// take a mesh having as external boundary a parallelogram
// and identify two pairs of opposite sides

// a quotient manifold with two action generators will be built
	
{	// we use the current manifold
	Manifold space = Manifold::working;
	assert ( space.exists() );
	Function coord = space.coordinates();
	assert ( coord.nb_of_components() == 2 );
	// we split 'coord' into its components
	Function x = coord [0],  y = coord [1];

	// we need to identify two translations, one which moves side_1 into side_2
	// and another one which moves side_3 into side_4

	Cell A = side_1 .first_vertex() .reverse();
	Cell B = side_1 .last_vertex();
	Cell C = side_2 .last_vertex();
	Cell D = side_2 .first_vertex() .reverse();

	double dx12 = x(D) - x(A), dy12 = y(D) - y(A);
	double norm12 = std::sqrt ( dx12*dx12 + dy12*dy12 );
	assert ( std::abs ( dx12 - ( x(C) - x(B) ) ) < 1.e-4 * norm12 );
	assert ( std::abs ( dy12 - ( y(C) - y(B) ) ) < 1.e-4 * norm12 );

	A = side_3 .first_vertex() .reverse();
	B = side_3 .last_vertex();
	C = side_4 .last_vertex();
	D = side_4 .first_vertex() .reverse();

	double dx34 = x(D) - x(A), dy34 = y(D) - y(A);
	double norm34 = std::sqrt ( dx34*dx34 + dy34*dy34 );
	assert ( std::abs ( dx34 - ( x(C) - x(B) ) ) < 1.e-4 * norm34 );
	assert ( std::abs ( dy34 - ( y(C) - y(B) ) ) < 1.e-4 * norm34 );

	// the desired translations are ( dx12, dy12 ), ( dx34, dy34 )
	Manifold::Action g12 ( tag::transforms, coord, tag::into, (x+dx12) && ( y+dy12) );
	Manifold::Action g34 ( tag::transforms, coord, tag::into, (x+dx34) && ( y+dy34) );
	Manifold manif_q = space .quotient ( g12, g34 );

	// we use a map -- for a faster code, we could use Cell::Core::hook
	std::map < Cell, std::pair < Cell, Manifold::Action > > corresp_ver;

	Mesh::Iterator it_ver = that->iterator ( tag::over_vertices );
	for ( it_ver .reset(); it_ver .in_range(); it_ver++ )
	{	Cell V = *it_ver;
		if ( V .belongs_to ( side_2 ) ) continue;
		if ( V .belongs_to ( side_4 ) ) continue;
		Cell new_V ( tag::vertex );
		x ( new_V ) = x (V);   y ( new_V ) = y (V);
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_map = corresp_ver .lower_bound (V);
		assert ( ( it_map == corresp_ver .end() ) or
		         ( corresp_ver.key_comp()(V,it_map->first) ) );
		corresp_ver .emplace_hint ( it_map, std::piecewise_construct,
		    std::forward_as_tuple (V), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { new_V, 0 } ) );  }
		// corresp_ver [ V ] = { new_V, 0 };

	assert ( side_1 .number_of ( tag::segments ) == side_2 .number_of ( tag::segments ) );
	Mesh::Iterator it1 = side_1 .iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it2 = side_2 .iterator ( tag::over_vertices, tag::require_order );
	for ( it1 .reset(), it2 .reset(); it1 .in_range(); it1++, it2++ )
	{	assert ( it2 .in_range() );
		Cell V = *it1;  Cell W = *it2;
		assert ( std::abs ( dx12 - ( x(W) - x(V) ) ) < 1.e-4 * norm12 );
		assert ( std::abs ( dy12 - ( y(W) - y(V) ) ) < 1.e-4 * norm12 );
		if ( V .belongs_to ( side_4 ) ) continue;
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_V = corresp_ver .find (V);
		assert ( it_V != corresp_ver .end() );
		assert ( it_V->second .second == 0 );
		Cell VV = it_V->second .first;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_W = corresp_ver .lower_bound (W);
		assert ( ( it_W == corresp_ver .end() ) or
		         ( corresp_ver.key_comp()(W,it_W->first) ) );
		corresp_ver .emplace_hint ( it_W, std::piecewise_construct,
		    std::forward_as_tuple (W), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { VV, g12 } ) );  }
		// corresp_ver [W] = { VV, g12 };
	assert ( not it2 .in_range() );

	Cell origin ( tag::non_existent ), corner ( tag::non_existent );
	assert ( side_3 .number_of ( tag::segments ) == side_4 .number_of ( tag::segments ) );
	Mesh::Iterator it3 = side_3 .iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it4 = side_4 .iterator ( tag::over_vertices, tag::require_order );
	for ( it3 .reset(), it4 .reset(); it3 .in_range(); it3++, it4++ )
	{	assert ( it4 .in_range() );
		Cell V = *it3;  Cell W = *it4;
		assert ( std::abs ( dx34 - ( x(W) - x(V) ) ) < 1.e-4 * norm34 );
		assert ( std::abs ( dy34 - ( y(W) - y(V) ) ) < 1.e-4 * norm34 );
		if ( W .belongs_to ( side_2 ) )
		{	std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_W = corresp_ver .find (W);
			assert ( it_W == corresp_ver .end() );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V = corresp_ver .find (V);
			assert ( it_V != corresp_ver .end() );
			assert ( it_V->second .second == g12 );
			origin = it_V->second .first;
			corner = W;
			continue;                                                             }
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_V = corresp_ver .find (V);
		assert ( it_V != corresp_ver .end() );
		assert ( it_V->second .second == 0 );
		Cell new_V = it_V->second .first;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_W = corresp_ver .lower_bound (W);
		assert ( ( it_W == corresp_ver .end() ) or
		         ( corresp_ver.key_comp()(W,it_W->first) ) );
		corresp_ver.emplace_hint ( it_W, std::piecewise_construct,
		    std::forward_as_tuple (W), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { new_V, g34 } ) );  }
		// corresp_ver [W] = { new_V, g34 };
	assert ( not it4 .in_range() );

	assert ( origin .exists() );  assert ( corner .exists() );
	// inspired in item 24 of the book : Scott Meyers, Effective STL
	std::map < Cell, std::pair < Cell, Manifold::Action > >
		::iterator it_c = corresp_ver .lower_bound ( corner );
	assert ( ( it_c == corresp_ver .end() ) or
	         ( corresp_ver.key_comp()(corner,it_c->first) ) );
	corresp_ver .emplace_hint ( it_c, std::piecewise_construct,
	    std::forward_as_tuple ( corner ), std::forward_as_tuple
	    ( std::pair < Cell, Manifold::Action > { origin, g12 + g34 } ) );
	// corresp_ver [ corner ] = { origin, g12 + g34 };
	
	return fold_common_four_sides
		( *that, corresp_ver, side_1, side_2, side_3, side_4, dim, keep_map, m );

}  // end of fold_four_sides with tag::build_new_vertices

//----------------------------------------------------------------------------------//

	
Mesh fold_four_sides    // hidden in anonymous namespace      // line 1694
( Mesh * that, const tag::Identify &, const Mesh & side_1,
               const tag::With &, const Mesh & side_2,
               const tag::Identify &, const Mesh & side_3,
               const tag::With &, const Mesh & side_4,
  const tag::UseExistingVertices &,
  const tag::ReturnMapBetween &, const tag::CellsOfDim &,
  size_t dim, bool keep_map, std::map < Cell, Cell > & m  )

// take a mesh having as external boundary a parallelogram
// and identify two pairs of opposite sides

// a quotient manifold with two action generators will be built
	
{	// we use the current manifold
	Manifold space = Manifold::working;
	assert ( space .exists() );
	Function coord = space .coordinates();
	assert ( coord .nb_of_components() == 2 );
	// we split 'coord' into its components
	Function x = coord [0], y = coord [1];

	// we need to identify two translations, one which moves side_1 into side_2
	// and another one which moves side_3 into side_4

	Cell A = side_1 .first_vertex() .reverse();
	Cell B = side_1 .last_vertex();
	Cell C = side_2 .last_vertex();
	Cell D = side_2 .first_vertex() .reverse();

	double dx12 = x(D) - x(A), dy12 = y(D) - y(A);
	double norm12 = std::sqrt ( dx12*dx12 + dy12*dy12 );
	assert ( std::abs ( dx12 - ( x(C) - x(B) ) ) < 1.e-4 * norm12 );
	assert ( std::abs ( dy12 - ( y(C) - y(B) ) ) < 1.e-4 * norm12 );

	A = side_3 .first_vertex() .reverse();
	B = side_3 .last_vertex();
	C = side_4 .last_vertex();
	D = side_4 .first_vertex() .reverse();

	double dx34 = x(D) - x(A), dy34 = y(D) - y(A);
	double norm34 = std::sqrt ( dx34*dx34 + dy34*dy34 );
	assert ( std::abs ( dx34 - ( x(C) - x(B) ) ) < 1.e-4 * norm34 );
	assert ( std::abs ( dy34 - ( y(C) - y(B) ) ) < 1.e-4 * norm34 );

	// the desired translations are ( dx12, dy12 ), ( dx34, dy34 )
	Manifold::Action g12 ( tag::transforms, coord, tag::into, (x+dx12) && ( y+dy12) );
	Manifold::Action g34 ( tag::transforms, coord, tag::into, (x+dx34) && ( y+dy34) );
	Manifold manif_q = space .quotient ( g12, g34 );

	// we use a map -- for a faster code, we could use Cell::Core::hook
	std::map < Cell, std::pair < Cell, Manifold::Action > > corresp_ver;

	Mesh::Iterator it_ver = that->iterator ( tag::over_vertices );
	for ( it_ver .reset(); it_ver .in_range(); it_ver++ )
	{	Cell V = *it_ver;
		if ( V .belongs_to ( side_2 ) ) continue;
		if ( V .belongs_to ( side_4 ) ) continue;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_map = corresp_ver .lower_bound (V);
		assert ( ( it_map == corresp_ver .end() ) or
		         ( corresp_ver.key_comp()(V,it_map->first) ) );
		corresp_ver .emplace_hint ( it_map, std::piecewise_construct,
		    std::forward_as_tuple (V), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { V, 0 } ) );  }
		// corresp_ver [V] = { V, 0 };

	assert ( side_1 .number_of ( tag::segments ) == side_2 .number_of ( tag::segments ) );
	Mesh::Iterator it1 = side_1 .iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it2 = side_2 .iterator ( tag::over_vertices, tag::require_order );
	for ( it1 .reset(), it2 .reset(); it1 .in_range(); it1++, it2++ )
	{	assert ( it2 .in_range() );
		Cell V = *it1;  Cell W = *it2;
		assert ( std::abs ( dx12 - ( x(W) - x(V) ) ) < 1.e-4 * norm12 );
		assert ( std::abs ( dy12 - ( y(W) - y(V) ) ) < 1.e-4 * norm12 );
		if ( V.belongs_to ( side_4 ) ) continue;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_W = corresp_ver .lower_bound (W);
		assert ( ( it_W == corresp_ver .end() ) or
		         ( corresp_ver.key_comp()(W,it_W->first) ) );
		corresp_ver .emplace_hint ( it_W, std::piecewise_construct,
		    std::forward_as_tuple (W), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { V, g12 } ) );  }
		// corresp_ver [W] = { V, g12 };
	assert ( not it2 .in_range() );

	Cell origin ( tag::non_existent ), corner ( tag::non_existent );
	assert ( side_3 .number_of ( tag::segments ) == side_4 .number_of ( tag::segments ) );
	Mesh::Iterator it3 = side_3 .iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it4 = side_4 .iterator ( tag::over_vertices, tag::require_order );
	for ( it3 .reset(), it4 .reset(); it3 .in_range(); it3++, it4++ )
	{	assert ( it4 .in_range() );
		Cell V = *it3;  Cell W = *it4;
		assert ( std::abs ( dx34 - ( x(W) - x(V) ) ) < 1.e-4 * norm34 );
		assert ( std::abs ( dy34 - ( y(W) - y(V) ) ) < 1.e-4 * norm34 );
		if ( W .belongs_to ( side_2 ) )
		{	std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_W = corresp_ver .find (W);
			assert ( it_W == corresp_ver .end() );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V = corresp_ver .find ( V );
			assert ( it_V != corresp_ver .end() );
			assert ( it_V->second .second == g12 );
			origin = it_V->second .first;
			corner = W;
			continue;                                                             }
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_W = corresp_ver .lower_bound (W);
		assert ( ( it_W == corresp_ver .end() ) or
		         ( corresp_ver.key_comp()(W,it_W->first) ) );
		corresp_ver .emplace_hint ( it_W, std::piecewise_construct,
		    std::forward_as_tuple (W), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { V, g34 } ) );  }
		// corresp_ver [W] = { V, g34 };
	assert ( not it4 .in_range() );

	assert ( origin .exists() );  assert ( corner .exists() );
	// inspired in item 24 of the book : Scott Meyers, Effective STL
	std::map < Cell, std::pair < Cell, Manifold::Action > >
		::iterator it_c = corresp_ver .lower_bound ( corner );
	assert ( ( it_c == corresp_ver .end() ) or
	         ( corresp_ver.key_comp()(corner,it_c->first) ) );
	corresp_ver .emplace_hint ( it_c, std::piecewise_construct,
	    std::forward_as_tuple ( corner ), std::forward_as_tuple
	    ( std::pair < Cell, Manifold::Action > { origin, g12 + g34 } ) );
	// corresp_ver [ corner ] = { origin, g12 + g34 };
	
	return fold_common_four_sides
		( *that, corresp_ver, side_1, side_2, side_3, side_4, dim, keep_map, m );

}  // end of fold_four_sides with tag::use_existing_vertices
	
//----------------------------------------------------------------------------------//

	
Mesh fold_common_six_sides    // hidden in anonymous namespace
( const Mesh & msh,
  const std::map < Cell, std::pair < Cell, Manifold::Action > > & corresp_ver,
  const Mesh & side_1, const Mesh & side_2, const Mesh & side_3, const Mesh & side_4,
  const Mesh & side_5, const Mesh & side_6,
  size_t dim, bool keep_map, std::map < Cell, Cell > & m                             )

// build 'corresp_seg', taking care to identify
//   side_1 with side_2, side_3 with side_4 and side_5 with side_6
// then call fold_common
	
{	assert ( msh .dim() == 2 );
	// we use a map -- for a faster code, we could use Cell::Core::hook
	std::map < Cell, Cell > corresp_seg;

	Mesh::Iterator it_seg = msh .iterator ( tag::over_segments );
	for ( it_seg .reset(); it_seg .in_range(); it_seg++ )
	{	Cell seg = *it_seg;
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::const_iterator it_base_rev = corresp_ver .find ( seg .base() .reverse() );
		assert ( it_base_rev != corresp_ver .end() );
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::const_iterator it_tip = corresp_ver .find ( seg .tip() );
		assert ( it_tip != corresp_ver .end() );
		if ( seg .belongs_to ( side_2, tag::not_oriented ) ) continue;
			// use corresponding segment on side_1
		if ( seg .belongs_to ( side_4, tag::not_oriented ) ) continue;
			// use corresponding segment on side_3
		if ( seg .belongs_to ( side_6, tag::not_oriented ) ) continue;
			// use corresponding segment on side_5
		Cell new_seg ( tag::segment,
		               it_base_rev->second .first .reverse(), it_tip->second .first );
		new_seg.winding() = it_tip->second .second - it_base_rev->second .second;
		// new_seg.winding() = corresp_ver [ seg .tip()             ] .second -
		//	                   corresp_ver [ seg .base() .reverse() ] .second  ;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, Cell > ::iterator it_map = corresp_seg .lower_bound ( seg );
		assert ( ( it_map == corresp_seg .end() ) or
		         ( corresp_seg.key_comp()(seg,it_map->first) ) );
		corresp_seg .emplace_hint ( it_map, std::piecewise_construct,
		  std::forward_as_tuple ( seg ), std::forward_as_tuple ( new_seg ) );     }
		// corresp_seg [ seg ] = new_seg;               

	Mesh::Iterator it_seg_1 = side_1 .iterator
		( tag::over_segments, tag::require_order, tag::force_positive );
	Mesh::Iterator it_seg_2 = side_2 .iterator
		( tag::over_segments, tag::require_order, tag::force_positive );
	for ( it_seg_1 .reset(), it_seg_2 .reset(); it_seg_1 .in_range(); it_seg_1++, it_seg_2++ )
	{	assert ( it_seg_2 .in_range() );
		Cell seg_1 = *it_seg_1;
		Cell seg_2 = *it_seg_2;
		std::map < Cell, Cell > ::iterator it = corresp_seg.find ( seg_1 );
		assert ( it != corresp_seg .end() );
		Cell new_seg_1 = it->second .reverse();
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, Cell > ::iterator it_map = corresp_seg .lower_bound ( seg_2 );
		assert ( ( it_map == corresp_seg .end() ) or
		         ( corresp_seg.key_comp()(seg_2,it_map->first) ) );
		corresp_seg.emplace_hint ( it_map, std::piecewise_construct,
		  std::forward_as_tuple ( seg_2 ), std::forward_as_tuple ( new_seg_1 ) );     }
		// corresp_seg [ seg_2 ] = corresp_seg [ seg_1 ];               
	assert ( not it_seg_2 .in_range() );
		
	Mesh::Iterator it_seg_3 = side_3 .iterator
		( tag::over_segments, tag::require_order, tag::force_positive );
	Mesh::Iterator it_seg_4 = side_4 .iterator
		( tag::over_segments, tag::require_order, tag::force_positive );
	for ( it_seg_3 .reset(), it_seg_4 .reset(); it_seg_3 .in_range(); it_seg_3++, it_seg_4++ )
	{	assert ( it_seg_4 .in_range() );
		Cell seg_3 = *it_seg_3;
		Cell seg_4 = *it_seg_4;
		std::map < Cell, Cell > ::iterator it = corresp_seg .find ( seg_3 );
		assert ( it != corresp_seg .end() );
		Cell new_seg_3 = it->second .reverse();
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, Cell > ::iterator it_map = corresp_seg .lower_bound ( seg_4 );
		assert ( ( it_map == corresp_seg .end() ) or
		         ( corresp_seg.key_comp()(seg_4,it_map->first) ) );
		corresp_seg .emplace_hint ( it_map, std::piecewise_construct,
		  std::forward_as_tuple ( seg_4 ), std::forward_as_tuple ( new_seg_3 ) );     }
		// corresp_seg [ seg_4 ] = corresp_seg [ seg_3 ];               
	assert ( not it_seg_4 .in_range() );
		
	Mesh::Iterator it_seg_5 = side_5 .iterator
		( tag::over_segments, tag::require_order, tag::force_positive );
	Mesh::Iterator it_seg_6 = side_6 .iterator
		( tag::over_segments, tag::require_order, tag::force_positive );
	for ( it_seg_5 .reset(), it_seg_6 .reset(); it_seg_5 .in_range(); it_seg_5++, it_seg_6++ )
	{	assert ( it_seg_6 .in_range() );
		Cell seg_5 = *it_seg_5;
		Cell seg_6 = *it_seg_6;
		std::map < Cell, Cell > ::iterator it = corresp_seg .find ( seg_5 );
		assert ( it != corresp_seg .end() );
		Cell new_seg_5 = it->second .reverse();
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, Cell > ::iterator it_map = corresp_seg .lower_bound ( seg_6 );
		assert ( ( it_map == corresp_seg .end() ) or
		         ( corresp_seg.key_comp()(seg_6,it_map->first) ) );
		corresp_seg .emplace_hint ( it_map, std::piecewise_construct,
		  std::forward_as_tuple ( seg_6 ), std::forward_as_tuple ( new_seg_5 ) );     }
		// corresp_seg [ seg_6 ] = corresp_seg [ seg_5 ];               
	assert ( not it_seg_6 .in_range() );
		
	return fold_common ( msh, corresp_seg, dim, keep_map, m );

}  // end of  fold_common_six_sides
	
//----------------------------------------------------------------------------------//

	
Mesh fold_six_sides       // hidden in anonymous namespace   // line 1942
( Mesh * that, const tag::Identify &, const Mesh & side_1,
               const tag::With &,     const Mesh & side_2,
               const tag::Identify &, const Mesh & side3,
               const tag::With &,     const Mesh & side4,
               const tag::Identify &, const Mesh & side5,
               const tag::With &,     const Mesh & side6,
  const tag::BuildNewVertices &,
  const tag::ReturnMapBetween &, const tag::CellsOfDim &,
  size_t dim, bool keep_map, std::map < Cell, Cell > & m  )

// take a mesh having as external boundary a hexagon
// and identify three pairs of opposite sides

// a quotient manifold with two action generators will be built
	
{	// we use the current manifold
	Manifold space = Manifold::working;
	assert ( space .exists() );
	Function coord = space .coordinates();
	assert ( coord .nb_of_components() == 2 );
	// we split 'coord' into its components
	Function x = coord [0], y = coord [1];

	// we need to identify two translations, one which moves side_1 into side_2
	// and another one which moves side3 into side4

	Cell A = side_1 .first_vertex() .reverse();
	Cell B = side_1 .last_vertex();
	Cell C = side_2 .last_vertex();
	Cell D = side_2 .first_vertex() .reverse();

	double dx12 = x(D) - x(A), dy12 = y(D) - y(A);
	double norm12 = std::sqrt ( dx12*dx12 + dy12*dy12 );
	assert ( std::abs ( dx12 - ( x(C) - x(B) ) ) < 1.e-4 * norm12 );
	assert ( std::abs ( dy12 - ( y(C) - y(B) ) ) < 1.e-4 * norm12 );

	Mesh side_3 ( tag::non_existent ), side_4 ( tag::non_existent );
	// we want side_1 and side_3 to touch
	if ( A .belongs_to ( side3 ) or B .belongs_to ( side3 ) )
	{	side_3 = side3; side_4 = side4;  }
	else
	{	assert ( A .belongs_to ( side4 ) or B .belongs_to ( side4 ) );
		side_3 = side4; side_4 = side3;                               }

	A = side_3 .first_vertex() .reverse();
	B = side_3 .last_vertex();
	C = side_4 .last_vertex();
	D = side_4 .first_vertex() .reverse();

	double dx34 = x(D) - x(A), dy34 = y(D) - y(A);
	double norm34 = std::sqrt ( dx34*dx34 + dy34*dy34 );
	assert ( std::abs ( dx34 - ( x(C) - x(B) ) ) < 1.e-4 * norm34 );
	assert ( std::abs ( dy34 - ( y(C) - y(B) ) ) < 1.e-4 * norm34 );

	Mesh side_5 ( tag::non_existent ), side_6 ( tag::non_existent );
	// we want side_3 and side_5 to touch
	if ( A .belongs_to ( side5 ) or B .belongs_to ( side5 ) )
	{	side_5 = side5; side_6 = side6;  }
	else
	{	assert ( A .belongs_to ( side6 ) or B .belongs_to ( side6 ) );
		side_5 = side6; side_6 = side5;                               }

	A = side_5 .first_vertex() .reverse();
	B = side_5 .last_vertex();
	C = side_6 .last_vertex();
	D = side_6 .first_vertex() .reverse();

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
	for ( it_ver .reset(); it_ver .in_range(); it_ver++ )
	{	Cell V = *it_ver;
		if ( V .belongs_to ( side_2 ) ) continue;
		if ( V .belongs_to ( side_4 ) ) continue;
		if ( V .belongs_to ( side_6 ) ) continue;
		Cell new_V ( tag::vertex );
		x ( new_V ) = x (V);   y ( new_V ) = y (V);
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_map = corresp_ver .lower_bound (V);
		assert ( ( it_map == corresp_ver .end() ) or
		         ( corresp_ver.key_comp()(V,it_map->first) ) );
		corresp_ver .emplace_hint ( it_map, std::piecewise_construct,
		    std::forward_as_tuple ( V ), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { new_V, 0 } ) );  }
		// corresp_ver [V] = { new_V, 0 };

	Cell V16 ( tag::non_existent ), V24 ( tag::non_existent );
	Cell V13 ( tag::non_existent ), V25 ( tag::non_existent );
	assert ( side_1 .number_of ( tag::segments ) == side_2 .number_of ( tag::segments ) );
	Mesh::Iterator it1 = side_1 .iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it2 = side_2 .iterator ( tag::over_vertices, tag::require_order );
	for ( it1 .reset(), it2 .reset(); it1 .in_range(); it1++, it2++ )
	{	assert ( it2 .in_range() );
		Cell V = *it1;  Cell W = *it2;
		assert ( std::abs ( dx12 - ( x(W) - x(V) ) ) < 1.e-4 * norm12 );
		assert ( std::abs ( dy12 - ( y(W) - y(V) ) ) < 1.e-4 * norm12 );
		if ( V .belongs_to ( side_3 ) )
		{	assert ( W .belongs_to ( side_5 ) );
			V13 = V;  V25 = W;               }
		if ( V.belongs_to ( side_6 ) )
		{	assert ( W .belongs_to ( side_4 ) );
			V16 = V;  V24 = W;  continue;    }
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_V = corresp_ver .find ( V );
		assert ( it_V != corresp_ver .end() );
		assert ( it_V->second .second == 0 );
		Cell VV = it_V->second .first;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_W = corresp_ver .lower_bound (W);
		assert ( ( it_W == corresp_ver .end() ) or
		         ( corresp_ver.key_comp()(W,it_W->first) ) );
		corresp_ver.emplace_hint ( it_W, std::piecewise_construct,
		    std::forward_as_tuple (W), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { VV, g12 } ) );  }
		// corresp_ver [W] = { VV, g12 };
	assert ( not it2.in_range() );
	assert ( V16 .exists() );  assert ( V24 .exists() );
	assert ( V13 .exists() );  assert ( V25 .exists() );

	Cell V35 ( tag::non_existent ), V46 ( tag::non_existent );
	assert ( side_3 .number_of ( tag::segments ) == side_4 .number_of ( tag::segments ) );
	Mesh::Iterator it3 = side_3 .iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it4 = side_4 .iterator ( tag::over_vertices, tag::require_order );
	for ( it3 .reset(), it4 .reset(); it3 .in_range(); it3++, it4++ )
	{	assert ( it4 .in_range() );
		Cell V = *it3;  Cell W = *it4;
		assert ( std::abs ( dx34 - ( x(W) - x(V) ) ) < 1.e-4 * norm34 );
		assert ( std::abs ( dy34 - ( y(W) - y(V) ) ) < 1.e-4 * norm34 );
		if ( V .belongs_to ( side_5 ) )
		{	assert ( W .belongs_to ( side_2 ) );
			V35 = V;
			assert ( W == V24 );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_W = corresp_ver .find (W);
			assert ( it_W == corresp_ver .end() );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V = corresp_ver .find (V);
			assert ( it_V != corresp_ver .end() );
			assert ( it_V->second .second == 0 );                                }
		if ( V .belongs_to ( side_1 ) )
		{	assert ( W .belongs_to ( side_6 ) );
			V46 = W;
			assert ( V == V13 );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_W = corresp_ver .find (W);
			assert ( it_W == corresp_ver .end() );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V = corresp_ver .find (V);
			assert ( it_V != corresp_ver .end() );
			assert ( it_V->second .second == 0 );                                }
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_V = corresp_ver .find (V);
		assert ( it_V != corresp_ver .end() );
		assert ( it_V->second .second == 0 );
		Cell new_V = it_V->second .first;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_W = corresp_ver .lower_bound (W);
		assert ( ( it_W == corresp_ver .end() ) or
		         ( corresp_ver.key_comp()(W,it_W->first) ) );
		corresp_ver .emplace_hint ( it_W, std::piecewise_construct,
		    std::forward_as_tuple (W), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { new_V, g34 } ) );  }
		// corresp_ver [W] = { new_V, g34 };
	assert ( not it4 .in_range() );
	assert ( V35 .exists() );  assert ( V46 .exists() );

	assert ( side_5 .number_of ( tag::segments ) == side_6 .number_of ( tag::segments ) );
	Mesh::Iterator it5 = side_5 .iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it6 = side_6 .iterator ( tag::over_vertices, tag::require_order );
	for ( it5 .reset(), it6 .reset(); it5 .in_range(); it5++, it6++ )
	{	assert ( it6 .in_range() );
		Cell V = *it5;  Cell W = *it6;
		assert ( std::abs ( dx56 - ( x(W) - x(V) ) ) < 1.e-4 * norm56 );
		assert ( std::abs ( dy56 - ( y(W) - y(V) ) ) < 1.e-4 * norm56 );
		if ( V .belongs_to ( side_2 ) )
		{	assert ( W .belongs_to ( side_4 ) );
			assert ( V == V25 );
			assert ( W == V46 );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_W = corresp_ver .find (W);
			assert ( it_W != corresp_ver .end() );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V13 = corresp_ver .find ( V13 );
			assert ( it_V13 != corresp_ver .end() );
			assert ( it_V13->second .second == 0 );
			Cell new_V13 = it_V13->second .first;
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V = corresp_ver .find (V);
			assert ( it_V != corresp_ver .end() );
			assert ( it_V->second .first == new_V13 );
			assert ( it_V->second .second == g12 );
			continue;                                                             }
		if ( V .belongs_to ( side_3 ) )
		{	assert ( W .belongs_to ( side_1 ) );
			assert ( V == V35 );
			assert ( W == V16 );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_W = corresp_ver .find (W);
			assert ( it_W == corresp_ver .end() );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V35 = corresp_ver .find ( V35 );
			assert ( it_V35 != corresp_ver .end() );
			assert ( it_V35->second .second == 0 );
			Cell new_V35 = it_V35->second .first;
			// inspired in item 24 of the book : Scott Meyers, Effective STL
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V16 = corresp_ver .lower_bound ( V16 );
			assert ( ( it_V16 == corresp_ver .end() ) or
		           ( corresp_ver.key_comp()(V16,it_V16->first) ) );
			corresp_ver .emplace_hint ( it_V16, std::piecewise_construct,
		      std::forward_as_tuple ( V16 ), std::forward_as_tuple
		      ( std::pair < Cell, Manifold::Action > { new_V35, g56 } ) );
			// corresp_ver [ V16 ] = { new_V35, g56 };
			continue;                                                                     }
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_V = corresp_ver .find ( V );
		assert ( it_V != corresp_ver .end() );
		assert ( it_V->second.second == 0 );
		Cell new_V = it_V->second.first;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_W = corresp_ver .lower_bound (W);
		assert ( ( it_W == corresp_ver .end() ) or
		         ( corresp_ver.key_comp()(W,it_W->first) ) );
		corresp_ver .emplace_hint ( it_W, std::piecewise_construct,
		    std::forward_as_tuple (W), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { new_V, g56 } ) );  }
		// corresp_ver [W] = { new_V, g56 };
	assert ( not it6 .in_range() );

	return fold_common_six_sides
		( *that, corresp_ver, side_1, side_2, side_3, side_4, side_5, side_6, dim, keep_map, m );

}  // end of fold_six_sides with tag::build_new_vertices
	
//----------------------------------------------------------------------------------//

	
Mesh fold_six_sides    // hidden in anonymous namespace      // line 2214
( Mesh * that, const tag::Identify &, const Mesh & side_1,
               const tag::With &,     const Mesh & side_2,
               const tag::Identify &, const Mesh & side3,
               const tag::With &, const Mesh & side4,
               const tag::Identify &, const Mesh & side5,
               const tag::With &, const Mesh & side6,
  const tag::UseExistingVertices &,
  const tag::ReturnMapBetween &, const tag::CellsOfDim &,
  size_t dim, bool keep_map, std::map < Cell, Cell > & m  )

// take a mesh having as external boundary a hexagon
// and identify three pairs of opposite sides

// a quotient manifold with two action generators will be built
		
{	// we use the current manifold
	Manifold space = Manifold::working;
	assert ( space .exists() );
	Function coord = space .coordinates();
	assert ( coord .nb_of_components() == 2 );
	// we split 'coord' into its components
	Function x = coord [0], y = coord [1];

	// we need to identify two translations, one which moves side_1 into side_2
	// and another one which moves side_3 into side_4

	Cell A = side_1 .first_vertex() .reverse();
	Cell B = side_1 .last_vertex();
	Cell C = side_2 .last_vertex();
	Cell D = side_2 .first_vertex() .reverse();

	double dx12 = x(D) - x(A), dy12 = y(D) - y(A);
	double norm12 = std::sqrt ( dx12*dx12 + dy12*dy12 );
	assert ( std::abs ( dx12 - ( x(C) - x(B) ) ) < 1.e-4 * norm12 );
	assert ( std::abs ( dy12 - ( y(C) - y(B) ) ) < 1.e-4 * norm12 );

	Mesh side_3 ( tag::non_existent ), side_4 ( tag::non_existent );
	// we want side_1 and side_3 to touch
	if ( A .belongs_to ( side3 ) or B .belongs_to ( side3 ) )
	{	side_3 = side3; side_4 = side4;  }
	else
	{	assert ( A .belongs_to ( side4 ) or B .belongs_to ( side4 ) );
		side_3 = side4; side_4 = side3;                                 }

	A = side_3 .first_vertex() .reverse();
	B = side_3 .last_vertex();
	C = side_4 .last_vertex();
	D = side_4 .first_vertex() .reverse();

	double dx34 = x(D) - x(A), dy34 = y(D) - y(A);
	double norm34 = std::sqrt ( dx34*dx34 + dy34*dy34 );
	assert ( std::abs ( dx34 - ( x(C) - x(B) ) ) < 1.e-4 * norm34 );
	assert ( std::abs ( dy34 - ( y(C) - y(B) ) ) < 1.e-4 * norm34 );

	Mesh side_5 ( tag::non_existent ), side_6 ( tag::non_existent );
	// we want side_3 and side_5 to touch
	if ( A .belongs_to ( side5 ) or B .belongs_to ( side5 ) )
	{	side_5 = side5; side_6 = side6;  }
	else
	{	assert ( A .belongs_to ( side6 ) or B .belongs_to ( side6 ) );
		side_5 = side6; side_6 = side5;                                 }

	A = side_5 .first_vertex() .reverse();
	B = side_5 .last_vertex();
	C = side_6 .last_vertex();
	D = side_6 .first_vertex() .reverse();

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
	Manifold manif_q = space .quotient ( g12, g34 );

	// a third translation, not used in the definition of the quotient manifold
	Manifold::Action g56 = signs[index_min_dif][0]*g12 + signs[index_min_dif][1]*g34;
	
	// we use a map -- for a faster code, we could use Cell::Core::hook
	std::map < Cell, std::pair < Cell, Manifold::Action > > corresp_ver;

	Mesh::Iterator it_ver = that->iterator ( tag::over_vertices );
	for ( it_ver .reset(); it_ver .in_range(); it_ver++ )
	{	Cell V = *it_ver;
		if ( V .belongs_to ( side_2 ) ) continue;
		if ( V .belongs_to ( side_4 ) ) continue;
		if ( V .belongs_to ( side_6 ) ) continue;
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_map = corresp_ver .lower_bound (V);
		assert ( ( it_map == corresp_ver .end() ) or
		         ( corresp_ver.key_comp()(V,it_map->first) ) );
		corresp_ver .emplace_hint ( it_map, std::piecewise_construct,
		    std::forward_as_tuple (V), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { V, 0 } ) );  }
		// corresp_ver [V] = { V, 0 };

	Cell V16 ( tag::non_existent ), V24 ( tag::non_existent );
	Cell V13 ( tag::non_existent ), V25 ( tag::non_existent );
	assert ( side_1 .number_of ( tag::segments ) == side_2 .number_of ( tag::segments ) );
	Mesh::Iterator it1 = side_1 .iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it2 = side_2 .iterator ( tag::over_vertices, tag::require_order );
	for ( it1 .reset(), it2 .reset(); it1 .in_range(); it1++, it2++ )
	{	assert ( it2 .in_range() );
		Cell V = *it1;  Cell W = *it2;
		assert ( std::abs ( dx12 - ( x(W) - x(V) ) ) < 1.e-4 * norm12 );
		assert ( std::abs ( dy12 - ( y(W) - y(V) ) ) < 1.e-4 * norm12 );
		if ( V .belongs_to ( side_3 ) )
		{	assert ( W .belongs_to ( side_5 ) );
			V13 = V;  V25 = W;                  }
		if ( V .belongs_to ( side_6 ) )
		{	assert ( W .belongs_to ( side_4 ) );
			V16 = V;  V24 = W;  continue;       }
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_W = corresp_ver .lower_bound (W);
		assert ( ( it_W == corresp_ver .end() ) or
		         ( corresp_ver.key_comp()(W,it_W->first) ) );
		corresp_ver .emplace_hint ( it_W, std::piecewise_construct,
		    std::forward_as_tuple (W), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { V, g12 } ) );  }
		// corresp_ver [W] = { V, g12 };
	assert ( not it2 .in_range() );
	assert ( V16 .exists() );  assert ( V24 .exists() );
	assert ( V13 .exists() );  assert ( V25 .exists() );

	Cell V35 ( tag::non_existent ), V46 ( tag::non_existent );
	assert ( side_3 .number_of ( tag::segments ) == side_4 .number_of ( tag::segments ) );
	Mesh::Iterator it3 = side_3 .iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it4 = side_4 .iterator ( tag::over_vertices, tag::require_order );
	for ( it3 .reset(), it4 .reset(); it3 .in_range(); it3++, it4++ )
	{	assert ( it4 .in_range() );
		Cell V = *it3;  Cell W = *it4;
		assert ( std::abs ( dx34 - ( x(W) - x(V) ) ) < 1.e-4 * norm34 );
		assert ( std::abs ( dy34 - ( y(W) - y(V) ) ) < 1.e-4 * norm34 );
		if ( V .belongs_to ( side_5 ) )
		{	assert ( W .belongs_to ( side_2 ) );
			V35 = V;
			assert ( W == V24 );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_W = corresp_ver .find (W);
			assert ( it_W == corresp_ver .end() );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V = corresp_ver .find (V);
			assert ( it_V != corresp_ver .end() );
			assert ( it_V->second .second == 0 );                                }
		if ( V .belongs_to ( side_1 ) )
		{	assert ( W .belongs_to ( side_6 ) );
			V46 = W;
			assert ( V == V13 );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_W = corresp_ver .find (W);
			assert ( it_W == corresp_ver .end() );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V = corresp_ver .find (V);
			assert ( it_V != corresp_ver .end() );
			assert ( it_V->second .second == 0 );                                }
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_W = corresp_ver .lower_bound (W);
		assert ( ( it_W == corresp_ver .end() ) or
		         ( corresp_ver.key_comp()(W,it_W->first) ) );
		corresp_ver .emplace_hint ( it_W, std::piecewise_construct,
		    std::forward_as_tuple (W), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { V, g34 } ) );  }
		// corresp_ver [W] = { V, g34 };
	assert ( not it4 .in_range() );
	assert ( V35 .exists() );  assert ( V46 .exists() );

	assert ( side_5 .number_of ( tag::segments ) == side_6 .number_of ( tag::segments ) );
	Mesh::Iterator it5 = side_5 .iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it6 = side_6 .iterator ( tag::over_vertices, tag::require_order );
	for ( it5 .reset(), it6 .reset(); it5 .in_range(); it5++, it6++ )
	{	assert ( it6 .in_range() );
		Cell V = *it5;  Cell W = *it6;
		assert ( std::abs ( dx56 - ( x(W) - x(V) ) ) < 1.e-4 * norm56 );
		assert ( std::abs ( dy56 - ( y(W) - y(V) ) ) < 1.e-4 * norm56 );
		if ( V .belongs_to ( side_2 ) )
		{	assert ( W .belongs_to ( side_4 ) );
			assert ( V == V25 );
			assert ( W == V46 );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_W = corresp_ver .find (W);
			assert ( it_W != corresp_ver .end() );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V13 = corresp_ver .find ( V13 );
			assert ( it_V13 != corresp_ver .end() );
			assert ( it_V13->second .first == V13 );
			assert ( it_V13->second .second == 0 );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V = corresp_ver .find (V);
			assert ( it_V != corresp_ver .end() );
			assert ( it_V->second .first == V13 );
			assert ( it_V->second .second == g12 );
			continue;                                                             }
		if ( V .belongs_to ( side_3 ) )
		{	assert ( W .belongs_to ( side_1 ) );
			assert ( V == V35 );
			assert ( W == V16 );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_W = corresp_ver .find (W);
			assert ( it_W == corresp_ver .end() );
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V35 = corresp_ver .find ( V35 );
			assert ( it_V35 != corresp_ver .end() );
			assert ( it_V35->second .first == V35 );
			assert ( it_V35->second .second == 0 );
			// inspired in item 24 of the book : Scott Meyers, Effective STL
			std::map < Cell, std::pair < Cell, Manifold::Action > >
				::iterator it_V16 = corresp_ver .lower_bound ( V16 );
			assert ( ( it_V16 == corresp_ver .end() ) or
		           ( corresp_ver.key_comp()(V16,it_V16->first) ) );
			corresp_ver .emplace_hint ( it_V16, std::piecewise_construct,
		      std::forward_as_tuple ( V16 ), std::forward_as_tuple
		      ( std::pair < Cell, Manifold::Action > { V35, g56 } ) );
			// corresp_ver [ V16 ] = { V35, g56 };
			continue;                                                                     }
		// inspired in item 24 of the book : Scott Meyers, Effective STL
		std::map < Cell, std::pair < Cell, Manifold::Action > >
			::iterator it_W = corresp_ver .lower_bound (W);
		assert ( ( it_W == corresp_ver .end() ) or
		         ( corresp_ver.key_comp()(W,it_W->first) ) );
		corresp_ver .emplace_hint ( it_W, std::piecewise_construct,
		    std::forward_as_tuple (W), std::forward_as_tuple
		    ( std::pair < Cell, Manifold::Action > { V, g56 } ) );  }
		// corresp_ver [W] = { V, g56 };
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
	                       tag::cells_of_dim, dim, true, m                        );  }


Mesh Mesh::fold ( const tag::BuildNewVertices & )

{	std::map < Cell, Cell > m;
	return fold_no_sides ( this, tag::build_new_vertices, tag::return_map_between,
	                       tag::cells_of_dim, 2, false, m                         );  }

	
Mesh Mesh::fold ( const tag::UseExistingVertices &,
                  const tag::ReturnMapBetween &, const tag::CellsOfDim &,
                  size_t dim, std::map < Cell, Cell > & m                )

{	return fold_no_sides ( this, tag::use_existing_vertices, tag::return_map_between,
	                       tag::cells_of_dim, dim, true, m                           );  }
	
	
Mesh Mesh::fold ( const tag::UseExistingVertices & )

{	std::map < Cell, Cell > m;
	return fold_no_sides ( this, tag::use_existing_vertices, tag::return_map_between,
	                       tag::cells_of_dim, 2, false, m                            );  }

//----------------------------------------------------------------------------------//

	
// take a mesh whose external boundary has two parallel sides and identify these two sides
// they do not touch each other

// a quotient manifold will be built with one action generator
	

Mesh Mesh::fold ( const tag::Identify &, const Mesh & side_1,
                  const tag::With &,     const Mesh & side_2,
                  const tag::BuildNewVertices &,
                  const tag::ReturnMapBetween &, const tag::CellsOfDim &,
                  size_t dim, std::map < Cell, Cell > & m                )

{	return fold_two_sides ( this, tag::identify, side_1, tag::with, side_2,  // line 1283
	                        tag::build_new_vertices, tag::return_map_between,
	                        tag::cells_of_dim, dim, true, m                  );  }
	
	
Mesh Mesh::fold ( const tag::Identify &, const Mesh & side_1,
                  const tag::With &,     const Mesh & side_2,
                  const tag::BuildNewVertices &              )

{	std::map < Cell, Cell > m;
	return fold_two_sides ( this, tag::identify, side_1, tag::with, side_2,  // line 1283
	                        tag::build_new_vertices, tag::return_map_between,
	                        tag::cells_of_dim, 2, false, m                   );  }
	
	
Mesh Mesh::fold ( const tag::Identify &, const Mesh & side_1,
                  const tag::With &,     const Mesh & side_2,
                  const tag::UseExistingVertices &,
                  const tag::ReturnMapBetween &, const tag::CellsOfDim &,
                  size_t dim, std::map < Cell, Cell > & m                )

{	return fold_two_sides ( this, tag::identify, side_1, tag::with, side_2,  // line 1373
	                        tag::use_existing_vertices, tag::return_map_between,
	                        tag::cells_of_dim, dim, true, m                     );  }

	
Mesh Mesh::fold ( const tag::Identify &, const Mesh & side_1,
                  const tag::With &,     const Mesh & side_2,
                  const tag::UseExistingVertices &           )

{	std::map < Cell, Cell > m;
	return fold_two_sides ( this, tag::identify, side_1, tag::with, side_2,  // line 1373
	                        tag::use_existing_vertices, tag::return_map_between,
	                        tag::cells_of_dim, 2, false, m                      );  }

//----------------------------------------------------------------------------------//

	
// take a mesh having as external boundary a parallelogram
// and identify two pairs of opposite sides

// a quotient manifold with two action generators will be built


Mesh Mesh::fold ( const tag::Identify &, const Mesh & side_1,
                  const tag::With &,     const Mesh & side_2,
                  const tag::Identify &, const Mesh & side_3,
                  const tag::With &,     const Mesh & side_4,
                  const tag::BuildNewVertices &,
                  const tag::ReturnMapBetween &, const tag::CellsOfDim &,
                  size_t dim, std::map < Cell, Cell > & m                )

{	return fold_four_sides ( this, tag::identify, side_1, tag::with, side_2,   // line 1544
	                         tag::identify, side_3, tag::with, side_4,
	                         tag::build_new_vertices, tag::return_map_between,
	                         tag::cells_of_dim, dim, true, m                  );  }
	
	
Mesh Mesh::fold ( const tag::Identify &, const Mesh & side_1,
                  const tag::With &,     const Mesh & side_2,
                  const tag::Identify &, const Mesh & side_3,
                  const tag::With &,     const Mesh & side_4,
                  const tag::BuildNewVertices &              )

{	std::map < Cell, Cell > m;
	return fold_four_sides ( this, tag::identify, side_1, tag::with, side_2,   // line 1544
	                         tag::identify, side_3, tag::with, side_4,
	                         tag::build_new_vertices, tag::return_map_between,
	                         tag::cells_of_dim, 2, false, m                   );  }

	
Mesh Mesh::fold ( const tag::Identify &, const Mesh & side_1,
                  const tag::With &,     const Mesh & side_2,
                  const tag::Identify &, const Mesh & side_3,
                  const tag::With &,     const Mesh & side_4,
                  const tag::UseExistingVertices &,
                  const tag::ReturnMapBetween &, const tag::CellsOfDim &,
                  size_t dim, std::map < Cell, Cell > & m                )

{	return fold_four_sides ( this, tag::identify, side_1, tag::with, side_2,  // line 1694
	                         tag::identify, side_3, tag::with, side_4,
	                         tag::use_existing_vertices, tag::return_map_between,
	                         tag::cells_of_dim, dim, true, m                     );  }
	
	
Mesh Mesh::fold ( const tag::Identify &, const Mesh & side_1,
                  const tag::With &,     const Mesh & side_2,
                  const tag::Identify &, const Mesh & side_3,
                  const tag::With &,     const Mesh & side_4,
                  const tag::UseExistingVertices &           )

{	std::map < Cell, Cell > m;
	return fold_four_sides ( this, tag::identify, side_1, tag::with, side_2,  // line 1694
	                         tag::identify, side_3, tag::with, side_4,
	                         tag::use_existing_vertices, tag::return_map_between,
	                         tag::cells_of_dim, 2, false, m                      );  }

//----------------------------------------------------------------------------------//


// take a mesh having as external boundary a hexagon
// and identify three pairs of opposite sides

// a quotient manifold with two action generators will be built
	

Mesh Mesh::fold ( const tag::Identify &, const Mesh & side_1,
                  const tag::With &,     const Mesh & side_2,
                  const tag::Identify &, const Mesh & side_3,
                  const tag::With &,     const Mesh & side_4,
                  const tag::Identify &, const Mesh & side_5,
                  const tag::With &,     const Mesh & side_6,
                  const tag::BuildNewVertices &,
                  const tag::ReturnMapBetween &, const tag::CellsOfDim &,
                  size_t dim, std::map < Cell, Cell > & m                )

{	return fold_six_sides ( this, tag::identify, side_1, tag::with, side_2,  // line 1942
	                        tag::identify, side_3, tag::with, side_4,
	                        tag::identify, side_5, tag::with, side_6,
	                        tag::build_new_vertices, tag::return_map_between,
	                        tag::cells_of_dim, dim, true, m                  );  }

	
Mesh Mesh::fold ( const tag::Identify &, const Mesh & side_1,
                  const tag::With &,     const Mesh & side_2,
                  const tag::Identify &, const Mesh & side_3,
                  const tag::With &,     const Mesh & side_4,
                  const tag::Identify &, const Mesh & side_5,
                  const tag::With &,     const Mesh & side_6,
                  const tag::BuildNewVertices &              )

{	std::map < Cell, Cell > m;
	return fold_six_sides ( this, tag::identify, side_1, tag::with, side_2,  // line 1942
	                        tag::identify, side_3, tag::with, side_4,
	                        tag::identify, side_5, tag::with, side_6,
	                        tag::build_new_vertices, tag::return_map_between,
	                        tag::cells_of_dim, 2, false, m                   );  }

	
Mesh Mesh::fold ( const tag::Identify &, const Mesh & side_1,
                  const tag::With &,     const Mesh & side_2,
                  const tag::Identify &, const Mesh & side_3,
                  const tag::With &,     const Mesh & side_4,
                  const tag::Identify &, const Mesh & side_5,
                  const tag::With &,     const Mesh & side_6,
                  const tag::UseExistingVertices &,
                  const tag::ReturnMapBetween &, const tag::CellsOfDim &,
                  size_t dim, std::map < Cell, Cell > & m                )

{	return fold_six_sides ( this, tag::identify, side_1, tag::with, side_2,   // line 2214
	                        tag::identify, side_3, tag::with, side_4,
	                        tag::identify, side_5, tag::with, side_6,
	                        tag::use_existing_vertices, tag::return_map_between,
	                        tag::cells_of_dim, dim, true, m                     );  }
  
	
Mesh Mesh::fold ( const tag::Identify &, const Mesh & side_1,
                  const tag::With &,     const Mesh & side_2,
                  const tag::Identify &, const Mesh & side_3,
                  const tag::With &,     const Mesh & side_4,
                  const tag::Identify &, const Mesh & side_5,
                  const tag::With &,     const Mesh & side_6,
                  const tag::UseExistingVertices &           )

{	std::map < Cell, Cell > m;
	return fold_six_sides ( this, tag::identify, side_1, tag::with, side_2,   // line 2214
	                         tag::identify, side_3, tag::with, side_4,
	                         tag::identify, side_5, tag::with, side_6,
	                         tag::use_existing_vertices, tag::return_map_between,
	                         tag::cells_of_dim, 2, false, m                      );  }

//----------------------------------------------------------------------------------//


// gs -q -dNOPAUSE -dBATCH -sDEVICE=png16 -sOUTPUTFILE=out.png in.ps


void Mesh::draw_ps ( std::string file_name ) const
	
{	// we use the current manifold
	Manifold space = Manifold::working;
	assert ( space .exists() );
	Function coord = space .coordinates();
	assert ( coord .nb_of_components() == 2 );
	// we split 'coord' into its components
	Function x = coord [0], y = coord [1];

	double xmin, xmax, ymin, ymax, maxside;

	{ // just a block for hiding variables
	Mesh::Iterator it = this->iterator ( tag::over_vertices );
	it .reset();
	assert( it .in_range() );
	Cell Vfirst = *it;
	xmin = xmax = x ( Vfirst );
	ymin = ymax = y ( Vfirst );	
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
	for ( it .reset() ; it .in_range(); it++ )
	{	Cell seg = *it;
		Cell base = seg .base() .reverse();
		Cell tip  = seg .tip();
		file_ps << x ( base ) << " " << y ( base ) << " moveto" << std::endl;
		file_ps << x ( tip ) << " " << y ( tip ) << " lineto stroke" << std::endl;  }
	} // just a block for hiding variables
																					
	file_ps << "grestore" << std::endl;
	file_ps << "grestore" << std::endl << std::endl;
	file_ps << "showpage" << std::endl;
	file_ps << "%%Trailer" << std::endl;
	file_ps << "%EOF" << std::endl;

	if ( not file_ps .good() )
	{	std::cerr << "error writing postscript file" << std::endl;
		exit (1);                                                   }

} // end of  Mesh::draw_ps

//----------------------------------------------------------------------------------//


void Mesh::draw_ps ( std::string file_name, const tag::Unfold &,
         const tag::OverRegion &, const Function::Inequality::Set & constraints ) const

// we draw several translations (or, more generally, transformations)
// of each segment, within the region described by the constraints
	
{	Manifold space = Manifold::working;
	assert ( space .exists() );  // we use the current (quotient) manifold
	Manifold::Quotient * manif_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space .core );
	assert ( manif_q );

	// the action group may have one or two generators
	size_t n = manif_q->actions .size();
	assert ( n == manif_q->winding_nbs .size() );
	if ( n == 1 ) this->draw_ps ( file_name, tag::unfold, tag::one_generator,
	                              tag::over_region, constraints               );
	else
	{	assert ( n == 2 );
		this->draw_ps ( file_name, tag::unfold, tag::two_generators,
		                tag::over_region, constraints               );  }

} // end of  Mesh::draw_ps with tag::unfold
	
//----------------------------------------------------------------------------------//


void Mesh::draw_ps ( std::string file_name, const tag::Unfold &, const tag::OneGenerator &,
               const tag::OverRegion &, const Function::Inequality::Set & constraints ) const

// we draw several translations (or, more generally, transformations)
// of each segment, within the region described by the constraints
	
{	Manifold space = Manifold::working;
	assert ( space .exists() );  // we use the current (quotient) manifold
	Manifold::Quotient * manif_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space .core );
	assert ( manif_q );
	Function coords_q = space .coordinates();
	Manifold mani_Eu = manif_q->base_space;  // underlying Euclidian manifold
	Function coords_Eu = mani_Eu .coordinates();
	assert ( coords_Eu .nb_of_components() == 2 );
	Function x = coords_Eu [0], y = coords_Eu [1];

	// here, the action group has one generator
	assert ( manif_q->actions .size() == 1 );
	assert ( manif_q->winding_nbs .size() == 1 );
	Manifold::Action g = manif_q->actions [0];

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
				touches_region = touches_region or constraints .on_cell ( shadow );
				a += seg .winding();
				coords_tip = coords_q ( tip, tag::winding, a );
				coords_Eu ( shadow ) = coords_tip;
				touches_region = touches_region or constraints .on_cell ( shadow );
				if ( touches_region )
				{	successful_round = true;
					double xx = coords_base [0], yy = coords_base [1];
					file_ps << xx << " " << yy << " moveto" << std::endl;
					if ( xx < xmin ) xmin = xx;
					if ( xx > xmax ) xmax = xx;
					if ( yy < ymin ) ymin = yy;
					if ( yy > ymax ) ymax = yy;
					xx = coords_tip [0];  yy = coords_tip [1];
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
	
	if ( not file_ps .good() )
	{	std::cerr << "error writing postscript file" << std::endl;
		exit (1);                                                   }

} // end of  Mesh::draw_ps with tag::unfold and tag::one_generator

//----------------------------------------------------------------------------------//


void Mesh::draw_ps ( std::string file_name, const tag::Unfold &, const tag::TwoGenerators &,
              const tag::OverRegion &, const Function::Inequality::Set & constraints ) const

// we draw several translations (or, more generally, transformations)
// of each segment, within the region described by the constraints
	
{	Manifold space = Manifold::working;
	assert ( space .exists() );  // we use the current (quotient) manifold
	Manifold::Quotient * manif_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space .core );
	assert ( manif_q );
	Function coords_q = space .coordinates();
	Manifold mani_Eu = manif_q->base_space;  // underlying Euclidian manifold
	Function coords_Eu = mani_Eu .coordinates();
	assert ( coords_Eu .nb_of_components() == 2 );
	Function x = coords_Eu [0], y = coords_Eu [1];

	// here, the action group has two generators
	assert ( manif_q->actions .size() == 2 );
	assert ( manif_q->winding_nbs .size() == 2 );
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
					touches_region = touches_region or constraints .on_cell ( shadow );
					a += seg .winding();
					coords_tip = coords_q ( tip, tag::winding, a );
					coords_Eu ( shadow ) = coords_tip;
					touches_region = touches_region or constraints .on_cell ( shadow );
					if ( touches_region )
					{	successful_round = true;
						double xx = coords_base [0], yy = coords_base [1];
						file_ps << xx << " " << yy << " moveto" << std::endl;
						if ( xx < xmin ) xmin = xx;
						if ( xx > xmax ) xmax = xx;
						if ( yy < ymin ) ymin = yy;
						if ( yy > ymax ) ymax = yy;
						xx = coords_tip [0];  yy = coords_tip [1];
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
	
	if ( not file_ps .good() )
	{	std::cerr << "error writing postscript file" << std::endl;
		exit (1);                                                   }

} // end of  Mesh::draw_ps with tag::unfold and tag::two_generators

//----------------------------------------------------------------------------------//


void Mesh::draw_ps ( std::string file_name,
         const tag::Unfold &, const std::vector < Manifold::Action > & v,
         const tag::OverRegion &, const Function::Inequality::Set & constraints ) const

// we draw several translations (or, more generally, transformations)
// of each segment, within the region described by the constraints
	
{	Manifold space = Manifold::working;
	assert ( space .exists() );  // we use the current (quotient) manifold
	Manifold::Quotient * manif_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space .core );
	Function coords_q = space .coordinates();
	Manifold mani_Eu = manif_q->base_space;  // underlying Euclidian manifold
	Function coords_Eu = mani_Eu .coordinates();
	assert ( coords_Eu .nb_of_components() == 2 );
	Function x = coords_Eu [0], y = coords_Eu [1];

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
	for ( it .reset() ; it .in_range(); it++ )
	{	Cell seg = *it;
		Cell base = seg .base() .reverse();
		Cell tip  = seg .tip();
		for ( std::vector < Manifold::Action > ::const_iterator
						it_v = v .begin(); it_v != v .end(); it_v++    )
		{	Manifold::Action a = *it_v;
			bool touches_region = false;
			coords_base = coords_q ( base, tag::winding, a );
			coords_Eu ( shadow ) = coords_base;
			touches_region = touches_region or constraints .on_cell ( shadow );
			a += seg .winding();
			coords_tip = coords_q ( tip, tag::winding, a );
			coords_Eu ( shadow ) = coords_tip;
			touches_region = touches_region or constraints .on_cell ( shadow );
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
	
	if ( not file_ps .good() )
	{	std::cerr << "error writing postscript file" << std::endl;
		exit (1);                                                   }

} // end of  Mesh::draw_ps with tag::unfold and tag::two_generators

//----------------------------------------------------------------------------------//


//  method below relies on some postscript macros for (very) rudimentary 3d drawings
//  available at https://github.com/cristian-barbarosie/manifem, file 3d.ps

void Mesh::draw_ps_3d ( std::string file_name ) const
	
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

	if ( not file_ps .good() )
	{	std::cerr << "error writing postscript file" << std::endl;
		exit (1);                                                   }

} // end of  Mesh::draw_ps_3d

//----------------------------------------------------------------------------------//


namespace {  // anonymous namespace, mimics static linkage
	
inline Mesh unfold_common ( const Mesh & that,
   const std::map < Cell, std::vector < std::pair < Manifold::Action, Cell > > > & built_ver,
                     const tag::OverRegion &, const Function::Inequality::Set & constraints  )

// 'constraints' not used
	
{	// build new segments joinig vertices built previously (by the calling code)
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
	std::map < Cell, std::pair < Cell,
		std::vector < std::pair < Manifold::Action, Cell > > > > built_faces;
	// similar to built_ver, but keep also a "base" vertex
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = that .iterator ( tag::over_cells_of_max_dim );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell cll = *it;
		Mesh::Iterator it_seg =
			cll .boundary() .iterator ( tag::over_segments, tag::require_order );
		it_seg .reset();  assert ( it_seg .in_range() );
		Cell ini_seg = *it_seg;
		Cell A = ini_seg .base() .reverse();
		std::map < Cell, std::vector < std::pair < Manifold::Action, Cell > > >
			::const_iterator it_bvA = built_ver .find (A);
		assert ( it_bvA != built_ver .end() );  // or rather 'continue' if not found ...
		std::pair < std::map < Cell, std::pair < Cell,
			std::vector < std::pair < Manifold::Action, Cell > > > > ::iterator, bool >
			p = built_faces .insert ( std::pair < Cell, std::pair < Cell,
			std::vector < std::pair < Manifold::Action, Cell > > > > ( cll, { A, { } } ) );
		assert ( p .second );
		std::vector < std::pair < Manifold::Action, Cell > > & vec = p .first->second .second;
		// 'vec' is empty for now
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
	  		           tag::segments, tag::one_dummy_wrapper        ),
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
			vec .push_back ( { act_base, new_cll } );
		}	 // end of loop over actions
	}  // end of loop over cells
	} // just a block of code for hiding 'it'

	if ( true )  // dim == 2
	{	Mesh result ( tag::fuzzy, tag::of_dim, 2 );
		Mesh::Iterator it = that .iterator ( tag::over_cells_of_max_dim );
		for ( it .reset(); it .in_range(); it++ )
		{	Cell cll = *it;
			std::map < Cell, std::pair < Cell,
				std::vector < std::pair < Manifold::Action, Cell > > > > ::iterator
				it_f = built_faces .find ( cll );
			assert ( it_f != built_faces .end() );
			Cell base_vertex = it_f->first;
			std::vector < std::pair < Manifold::Action, Cell > > & vec = it_f->second .second;
			for ( std::vector < std::pair < Manifold::Action, Cell > > ::const_iterator
			      itA = vec .begin(); itA != vec .end(); itA++                         )
			{ Cell new_cll = itA->second;;
				new_cll .add_to_mesh ( result );  }                                               }
		return result;                                                                           }

	assert ( false );
	Mesh result ( tag::fuzzy, tag::of_dim, 3 );
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
               const tag::ReturnMapBetween &, const tag::CellsOfDim &, size_t dim,
               std::map < Cell, std::pair < Cell, Manifold::Action > > & mapping,
               bool fill_mapping                                                      )
	
// if last argument is true,
// besides the unfolded Mesh, return a mapping giving, for each vertex of the new mesh,
// the corresponding vertex in 'that' mesh together with a winding (action)
	
{	Manifold space = Manifold::working;
	assert ( space.exists() );  // we use the current (quotient) manifold
	Manifold::Quotient * manif_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space.core );
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
							( std::pair < Manifold::Action, Cell > ( a, new_ver ) );                 }
					ii += directions[d][0];
					jj += directions[d][1];                                                        }  }
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
               const tag::ReturnMapBetween &, const tag::CellsOfDim &, size_t dim,
               std::map < Cell, std::pair < Cell, Manifold::Action > > & mapping,
               bool fill_mapping                                                      )
	
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


inline Mesh unfold_local ( const Mesh & that, const std::vector < tag::Util::Action > & aa,
               const tag::OverRegion &, const Function::Inequality::Set & constraints,
               const tag::ReturnMapBetween &, const tag::CellsOfDim &, size_t dim,
               std::map < Cell, std::pair < Cell, Manifold::Action > > & mapping,
               bool fill_mapping                                                           )
	
// if last argument is true,
// besides the unfolded Mesh, return a mapping giving, for each vertex of the new mesh,
// the corresponding vertex in 'that' mesh together with a winding (action)
	
{	Manifold space = Manifold::working;
	assert ( space.exists() );  // we use the current (quotient) manifold
	Manifold::Quotient * manif_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space.core );
	Function coords_q = space .coordinates();
	Manifold mani_Eu = manif_q->base_space;  // underlying Euclidian manifold
	Function coords_Eu = mani_Eu .coordinates();
	std::cout << "global.cpp line 6479" << std::endl;

	Cell shadow ( tag::vertex );
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
		for ( std::vector < Manifold::Action > ::const_iterator
						it_aa = aa .begin(); it_aa != aa .end(); it_aa ++ )
		{	Manifold::Action a = *it_aa;
			std::vector < double > coords_ver = coords_q ( ver, tag::winding, a );
			coords_Eu ( shadow ) = coords_ver;
			if ( constraints .on_cell ( shadow ) )
			{	Cell new_ver ( tag::vertex );  coords_Eu ( new_ver ) = coords_ver;
				if ( fill_mapping )
				mapping .insert ( std::pair < Cell, std::pair < Cell, Manifold::Action > >
													( new_ver, std::pair < Cell, Manifold::Action > ( ver, a ) ) );
				itbv .first->second .push_back
					( std::pair < Manifold::Action, Cell > ( a, new_ver ) );                        }  }  }
	
	Mesh result = unfold_common ( that, built_ver, tag::over_region, constraints );

	mani_Eu .set_as_working_manifold();

	return result;

} // end of unfold_local

} // end of anonymous namespace

//----------------------------------------------------------------------------------//


Mesh Mesh::unfold ( const tag::OverRegion &, const Function::Inequality::Set & constraints,
                    const tag::ReturnMapBetween &, const tag::CellsOfDim &, size_t dim,
                    std::map < Cell, std::pair < Cell, Manifold::Action > > & mapping      )
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


Mesh Mesh::unfold ( const std::vector < tag::Util::Action > & aa,
	                  const tag::OverRegion &, const Function::Inequality::Set & constraints,
                    const tag::ReturnMapBetween &, const tag::CellsOfDim &, size_t dim,
                    std::map < Cell, std::pair < Cell, Manifold::Action > > & mapping      )
const
	
// take a mesh and unfold it over a given region of the plane
// so we can export the resulting mesh in msh format

// besides the unfolded Mesh, return a mapping giving, for each vertex of the new mesh,
// the corresponding vertex in 'this' mesh together with a winding (action)
	
{ assert ( dim == 0 );  // 'mapping' required from vertices to vertices

	return unfold_local ( *this, aa, tag::over_region, constraints,
	           tag::return_map_between, tag::cells_of_dim, 0, mapping, true );
	// last argument true means fill 'mapping'
	
} // end of Mesh::unfold


Mesh Mesh::unfold ( const std::vector < tag::Util::Action > & aa,
	                  const tag::OverRegion &, const Function::Inequality::Set & constraints )
const
	
// take a mesh and unfold it over a given region of the plane
// so we can export the resulting mesh in msh format

{	std::map < Cell, std::pair < Cell, Manifold::Action > > mapping;

	return unfold_local ( *this, aa, tag::over_region, constraints,
	                      tag::return_map_between, tag::cells_of_dim, 0, mapping, false  );
	// last argument false means do not bother with 'mapping'
	
} // end of Mesh::unfold

//----------------------------------------------------------------------------------//


size_t & Cell::Numbering::Map::operator[] ( const Cell p )  // virtual from Cell::Numbering
{	return (*(this->map))[p];  }
