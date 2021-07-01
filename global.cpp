
// global.cpp 2021.04.06

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


void Mesh::pretty_constructor
( const tag::Segment &, const Cell & A, const Cell & B, const tag::DividedIn &, size_t n )
	
// see paragraph 11.2 in the manual
	
{	// we use the current manifold
	Manifold space = Manifold::working;
	assert ( space.exists() );
	
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
		seg.add_to (*this);
		prev_point = P.reverse();                          }
	Cell seg ( tag::segment, prev_point, B );
	seg.add_to (*this);                                        }

// code below does the same as code above
// above is pretty, slightly slower - below is ugly, slightly faster


Mesh::OneDim::Positive::Positive
( const tag::Segment &, Cell::Negative::Vertex * A, Cell::Positive::Vertex * B,
  const tag::DividedIn &, size_t n                                              )

:	Positive()

{	// we use the current manifold
	Manifold::Core * space = Manifold::working.core;
	assert ( space );
	
	assert ( not A->is_positive() );
	assert ( A->reverse_p );
	Cell::Positive::Vertex * posA = ( Cell::Positive::Vertex * ) A->reverse_p;
	assert ( posA->is_positive() );
	assert ( B->is_positive() );
	assert ( A->get_dim() == 0 );
	assert ( B->get_dim() == 0 );
	
	Cell::Negative::Vertex * prev_point = A;
	for ( size_t i=1; i < n; ++i )
	{	Cell::Positive::Vertex * P = new Cell::Positive::Vertex;
		double frac = double(i)/double(n);
		space->interpolate ( P, 1.-frac, posA, frac, B );
		Cell::Positive::Segment * seg = new Cell::Positive::Segment ( prev_point, P );
		seg->add_to (this);
		prev_point = ( Cell::Negative::Vertex * ) P->build_reverse();                  }
	Cell::Positive::Segment * seg = new Cell::Positive::Segment ( prev_point, B );
	seg->add_to (this);                                                                 }

//----------------------------------------------------------------------------------//


void Mesh::pretty_constructor
( const tag::Triangle &, const Mesh & AB, const Mesh & BC, const Mesh & CA )

// see paragraph 11.4 in the manual
	
{	// we use the current manifold
	Manifold space = Manifold::working;
	assert ( space.exists() );

	// sides must be split in the same number of segments :
	size_t N = AB.number_of ( tag::cells_of_dim, 1 );
	assert ( N == BC.number_of ( tag::cells_of_dim, 1 ) );
	assert ( N == CA.number_of ( tag::cells_of_dim, 1 ) );

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
	CellIterator it = AB.iter_over ( tag::segments );
	for ( it.reset(); it.in_range(); it++ ) ground.push_back ( *it );
	} // just a block of code for hiding 'it'

	// we shall use six points, two on AB, two on BC, two on CA
	// like shadows of the point currently buing built
	Cell Q_AB = A, P_BC = B, Q_CA = A;
	Cell seg_on_BC = BC.cell_in_front_of(B);

	for ( size_t i = 1; i < N; i++ ) // "vertical" movement
	{	// advance one level upwards and slightly right (parallel to CA)
		Q_AB = AB.cell_in_front_of(Q_AB).tip();
		P_BC = seg_on_BC.tip();
		Cell TA = CA.cell_behind(Q_CA);
		Cell P_CA = Q_CA = TA.base().reverse();
		Cell P_AB = A;  Cell Q_BC = C;
		Cell Q_AB_c = Q_AB;  // keep a copy of Q_AB
		std::list<Cell>::iterator it_ground = ground.begin();
		// build the first triangle on this layer
		Cell AP = *it_ground;
		Cell ground_ver = AP.tip();
		Cell previous_seg ( tag::segment, ground_ver.reverse(), P_CA );
		Cell tri ( tag::triangle, AP, previous_seg, TA );
		tri.add_to ( *this );  // 'this' is the mesh we are building
		Cell previous_ver = Q_CA;
		ceiling.clear();
		for ( size_t j = i+1; j <= N; j++ ) // "horizontal" movement
		{	// advance one step horizontally (parallel to AB)
			P_AB = AB.cell_in_front_of(P_AB).tip();
			Q_AB_c = AB.cell_in_front_of(Q_AB_c).tip();
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
				double s = frac_AB + frac_BC + frac_CA;  s *= 2.;
				frac_AB /= s;  frac_BC /= s;  frac_CA /= s;
				S = Cell ( tag::vertex );
				space.interpolate ( S, frac_AB, P_AB,  frac_AB, Q_AB_c,
														   frac_BC, P_BC,  frac_BC, Q_BC,
														   frac_CA, P_CA,  frac_CA, Q_CA  );  }
			Cell new_seg ( tag::segment, ground_ver.reverse(), S );
			Cell horizontal_seg ( tag::segment, S.reverse(), previous_ver );
		  Cell tri_1 ( tag::triangle, previous_seg.reverse(), new_seg, horizontal_seg );
		  tri_1.add_to ( *this );
			it_ground++; assert ( it_ground != ground.end() );
			Cell ground_seg = *it_ground;
			if ( j == N ) previous_seg = seg_on_BC;
			else
			{	ground_ver = ground_seg.tip();
				previous_seg = Cell ( tag::segment, ground_ver.reverse(), S );  }
		  Cell tri_2 ( tag::triangle, ground_seg, previous_seg, new_seg.reverse() );
			tri_2.add_to ( *this );
			previous_ver = S;
			// add horizontal_seg.reverse() to future ground
			ceiling.push_back ( horizontal_seg.reverse() );                              	  }
		seg_on_BC = BC.cell_in_front_of(P_BC);
		ground = ceiling;                                                                    }
		// improve by moving ceiling to ground, leaving ceiling empty !

	// last triangle
	Cell TA = CA.cell_behind(Q_CA);
	assert ( not ground.empty() );
	Cell AP = *(ground.begin());
	Cell tri ( tag::triangle, seg_on_BC, TA, AP );
	tri.add_to ( *this );
		
} // end of Mesh::pretty_constructor ( tag::triangle,... )

	
// code below does the same as code above - inline !


Mesh::Positive::Positive
( const tag::Triangle &, const Mesh & AB, const Mesh & BC, const Mesh & CA )

:	Mesh::Core ( tag::of_dimension, 3, tag::minus_one )
	
{	Mesh wrapper ( tag::whose_core_is, this );
	wrapper.pretty_constructor ( tag::triangle, AB, BC, CA );  }


//----------------------------------------------------------------------------------//


void Mesh::pretty_constructor ( const tag::Quadrangle &, const Mesh & south,
  const Mesh & east, const Mesh & north, const Mesh & west, const tag::WithTriangles & wt )

// 'wt' defaults to 'tag::not_with_triangles',
// which means 'cut_rectangles_in_half' defaults to 'false'

// see paragraph 11.3 in the manual
	
{	bool cut_rectangles_in_half = ( wt == tag::with_triangles );

	// we use the current manifold
	Manifold space = Manifold::working;
	assert ( space.exists() );

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
	CellIterator it = south.iter_over ( tag::segments );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell seg = *it;  horizon.push_back ( seg );  }
	} // just a block of code for hiding 'it'

	// start mesh generation
	CellIterator it_east = east.iter_over ( tag::vertices );
	CellIterator it_west = west.iter_over ( tag::vertices, tag::reverse );
	CellIterator it_south = south.iter_over ( tag::vertices );
	CellIterator it_north = north.iter_over ( tag::vertices, tag::reverse );
	it_east.reset(); it_east++;
	it_west.reset(); it_west++;
	for ( size_t i = 1; i < N_vert; ++i )
	{	std::list<Cell>::iterator it = horizon.begin();
		Cell seg = *it;
		Cell A = seg.base().reverse();
		Cell s4 = west.cell_behind ( A, tag::surely_exists );
		Cell D = s4.base().reverse();
		Cell ver_east = *it_east;
		Cell ver_west = *it_west;
		double frac_N = double(i) / double(N_vert),  alpha = frac_N * (1-frac_N);
		alpha = alpha*alpha*alpha;
		it_south.reset(); it_south++;
		it_north.reset(); it_north++;
		for ( size_t j = 1; j < N_horiz; j++ )
		{	Cell s1 = *it;  // 'it' points into the 'horizon' list
			Cell B = s1.tip();
			Cell C ( tag::vertex );  // create a new vertex
			Cell ver_south = *it_south;
			Cell ver_north = *it_north;
			double frac_E = double(j) / double(N_horiz),  beta = frac_E * (1-frac_E);
			beta = beta*beta*beta;
			double sum = alpha + beta,  aa = alpha/sum,  bb = beta/sum;
			space.interpolate ( C, bb*(1-frac_N), ver_south, aa*frac_E,     ver_east,     
		                         bb*frac_N,     ver_north, aa*(1-frac_E), ver_west );
			Cell s2 ( tag::segment, B.reverse(), C );  // create a new segment
			Cell s3 ( tag::segment, C.reverse(), D );  // create a new segment
			if ( cut_rectangles_in_half )
			{	Cell BD ( tag::segment, B.reverse(), D );  // create a new segment
				Cell T1 ( tag::triangle, BD.reverse(), s2, s3 );  // create a new triangle
				Cell T2 ( tag::triangle, BD, s4, s1 );  // create a new triangle
				T1.add_to (*this);  // 'this' is the mesh we are building
				T2.add_to (*this);                                 }
			else // with quadrilaterals
			{	Cell Q ( tag::rectangle, s1, s2, s3, s4 );  // create a new rectangle
				Q.add_to (*this);                                 }
			// 's3' is on the ceiling, we keep it in the 'horizon' list
			// it will be on the ground when we build the next layer of cells
			*it = s3.reverse(); // 'it' points into the 'horizon' list
			it++;
			D = C;
			s4 = s2.reverse();
			it_south++;  it_north++;
		} // end of for j
		it_south++;  it_north++;
		assert ( not it_south.in_range() );
		assert ( not it_north.in_range() );
		// last rectangle of this row, east side already exists
		Cell s1 = *it;
		Cell B = s1.tip();
		Cell s2 = east.cell_in_front_of ( B, tag::surely_exists );
		Cell C = s2.tip();
		Cell s3 ( tag::segment, C.reverse(), D );  // create a new segment
		if ( cut_rectangles_in_half )
		{	Cell BD ( tag::segment, B.reverse(), D );  // create a new segment
			Cell T1 ( tag::triangle, BD.reverse(), s2, s3 );  // create a new triangle
			Cell T2 ( tag::triangle, BD, s4, s1 );  // create a new triangle
			T1.add_to (*this);
			T2.add_to (*this);                                 }
		else // with quadrilaterals
		{	Cell Q ( tag::rectangle, s1, s2, s3, s4 );  // create a new rectangle
			Q.add_to (*this);                                 }
		*it = s3.reverse();
		it_east++;  it_west++;
		it++;  assert ( it == horizon.end() );
	} // end of for i
	it_east++;  it_west++;
	assert ( not it_east.in_range() );
	assert ( not it_west.in_range() );
	// last row of rectangles is different, north sides already exist
	std::list<Cell>::iterator it = horizon.begin();
	Cell s4 = west.cell_in_front_of ( NW, tag::surely_exists );
	Cell D = NW;
	for (size_t j=1; j < N_horiz; j++)
	{	Cell s1 = *it;
		// s1 == south.cell_in_front_of (A);
		Cell B = s1.tip();
		Cell s3 = north.cell_behind ( D );
		Cell C = s3.base().reverse();
		Cell s2 ( tag::segment, B.reverse(), C );  // create a new segment
		if ( cut_rectangles_in_half )
		{	Cell BD ( tag::segment, B.reverse(), D );  // create a new segment
			Cell T1 ( tag::triangle, BD.reverse(), s2, s3 );  // create a new triangle
			Cell T2 ( tag::triangle, BD, s4, s1 );  // create a new triangle
			T1.add_to (*this);
			T2.add_to (*this);                                 }
		else // with quadrilaterals
		{	Cell Q ( tag::rectangle, s1, s2, s3, s4 );  // create a new rectangle
			Q.add_to (*this);                           }
		it++;
		D = C;
		s4 = s2.reverse();                                          }
	// and the last rectangle of the last row
	Cell s1 = *it;
	Cell B = s1.tip();
	Cell s2 = east.cell_in_front_of (B);
	Cell C = s2.tip();
	assert ( C == NE );
	Cell s3 = north.cell_behind ( D );
	if ( cut_rectangles_in_half )
	{	Cell BD ( tag::segment, B.reverse(), D );
		Cell T1 ( tag::triangle, BD.reverse(), s2, s3 );
		Cell T2 ( tag::triangle, BD, s4, s1 );
		T1.add_to (*this);
		T2.add_to (*this);                                 }
	else // with quadrilaterals
	{	Cell Q ( tag::rectangle, s1, s2, s3, s4 );
		Q.add_to (*this);                                 }
	it++;  assert ( it == horizon.end() );

} // end of Mesh::pretty_constructor ( tag::Quadrangle...)

// code below does the same as code above
// above is pretty, slightly slower - below is ugly, slightly faster


void Mesh::Positive::build_rectangle ( const Mesh & south, const Mesh & east,
  const Mesh & north, const Mesh & west, bool cut_rectangles_in_half )

{	// we use the current manifold
	Manifold::Core * space = Manifold::working.core;
	assert ( space );
	
	// recover corners from the sides
	Cell::Positive::Vertex * SW = (Cell::Positive::Vertex*) south.first_vertex().core->reverse_p;
	assert ( SW == west.last_vertex().core );
	Cell::Positive::Vertex * SE = (Cell::Positive::Vertex*) east.first_vertex().core->reverse_p;
	assert ( SE == south.last_vertex().core );
	Cell::Positive::Vertex * NE = (Cell::Positive::Vertex*) north.first_vertex().core->reverse_p;
	assert ( NE == east.last_vertex().core );
	Cell::Positive::Vertex * NW = (Cell::Positive::Vertex*) west.first_vertex().core->reverse_p;
	assert ( NW == north.last_vertex().core );
	size_t N_horiz = south.number_of ( tag::segments );
	assert ( N_horiz == north.number_of ( tag::segments ) );
	size_t N_vert = east.number_of ( tag::segments );
	assert ( N_vert == west.number_of ( tag::segments ) );

	// since the boundary is closed, all vertices must have reverse
	
	// prepare horizon
	std::list < Cell::Core * > horizon;
	{ // just a block of code for hiding 'it'
	CellIterator it = south.iter_over ( tag::segments );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * seg = (*it).core;  horizon.push_back ( seg );  }
	} // just a block of code for hiding 'it'

	// start mesh generation
	CellIterator it_east = east.iter_over ( tag::vertices );
	CellIterator it_west = west.iter_over ( tag::vertices, tag::reverse );
	CellIterator it_south = south.iter_over ( tag::vertices );
	CellIterator it_north = north.iter_over ( tag::vertices, tag::reverse );
	it_east.reset(); it_east++;
	it_west.reset(); it_west++;
	for ( size_t i = 1; i < N_vert; ++i )
	{	std::list<Cell::Core*>::iterator it = horizon.begin();
		Cell::Core * seg = *it;
		Cell::Positive::Vertex * A = (Cell::Positive::Vertex*) seg->base()->reverse_p;
		Cell::Core * s4 = west.cell_behind ( A, tag::surely_exists );
		Cell::Positive::Vertex * D = (Cell::Positive::Vertex*) s4->base()->reverse_p;
		Cell::Positive::Vertex * ver_east = (Cell::Positive::Vertex*) (*it_east).core;
		Cell::Positive::Vertex * ver_west = (Cell::Positive::Vertex*) (*it_west).core;
		double frac_N = double(i) / double(N_vert),  alpha = frac_N * (1-frac_N);
		alpha = alpha*alpha*alpha;
		it_south.reset(); it_south++;
		it_north.reset(); it_north++;
		for ( size_t j = 1; j < N_horiz; j++ )
		{	Cell::Core * s1 = *it;  // 'it' points into the 'horizon' list
			Cell::Positive::Vertex * B = (Cell::Positive::Vertex*) s1->tip();
			Cell::Positive::Vertex * C = new Cell::Positive::Vertex;
			C->reverse_p = new Cell::Negative::Vertex ( tag::reverse_of, C );
			// we could use C->reverse ( tag::build_if_not_exists )
			// but here we know for a fact that the reverse does not exist
			Cell::Positive::Vertex * ver_south = (Cell::Positive::Vertex*) (*it_south).core;
			Cell::Positive::Vertex * ver_north = (Cell::Positive::Vertex*) (*it_north).core;
			double frac_E = double(j) / double(N_horiz),  beta = frac_E * (1-frac_E);
			beta = beta*beta*beta;
			double sum = alpha + beta,  aa = alpha/sum,  bb = beta/sum;
			space->interpolate ( C, bb*(1-frac_N), ver_south, aa*frac_E,     ver_east,     
		                         bb*frac_N,     ver_north, aa*(1-frac_E), ver_west );
			Cell::Positive::Segment * s2 = new Cell::Positive::Segment
				( (Cell::Negative::Vertex*) B->reverse_p, C );
			s2->reverse_p = new Cell::Negative::Segment ( tag::reverse_of, s2 );
			Cell::Positive::Segment * s3 = new Cell::Positive::Segment
				( (Cell::Negative::Vertex*) C->reverse_p, D );
			s3->reverse_p = new Cell::Negative::Segment ( tag::reverse_of, s3 );
			if ( cut_rectangles_in_half )
			{	Cell::Positive::Segment * BD = new Cell::Positive::Segment
					( (Cell::Negative::Vertex*) B->reverse_p, D );
				BD->reverse_p = new Cell::Negative::Segment ( tag::reverse_of, BD );
				Cell::Positive * T1 = new Cell::Positive ( tag::triangle, BD->reverse_p, s2, s3 );
				Cell::Positive * T2 = new Cell::Positive ( tag::triangle, BD, s4, s1 );
				T1->add_to (this);  // 'this' is the mesh we are building
				T2->add_to (this);                                                                  }
			else // with quadrilaterals
			{	Cell::Positive * Q = new Cell::Positive ( tag::rectangle, s1, s2, s3, s4 );
				Q->add_to (this);                                                             }
			// 's3' is on the ceiling, we keep it in the 'horizon' list
			// it will be on the ground when we build the next layer of cells
			*it = s3->reverse_p; // 'it' points into the 'horizon' list
			it++;
			D = C;
			s4 = s2->reverse_p;
			it_south++;  it_north++;
		} // end of for j
		it_south++;  it_north++;
		assert ( not it_south.in_range() );
		assert ( not it_north.in_range() );
		// last rectangle of this row, east side already exists
		Cell::Core * s1 = *it;
		Cell::Positive::Vertex * B = (Cell::Positive::Vertex*) s1->tip();
		Cell::Core * s2 = east.cell_in_front_of ( B, tag::surely_exists );
		Cell::Positive::Vertex * C = (Cell::Positive::Vertex*) s2->tip();
		Cell::Positive::Segment * s3 = new Cell::Positive::Segment
			( (Cell::Negative::Vertex*) C->reverse_p, D );
		s3->reverse_p = new Cell::Negative::Segment ( tag::reverse_of, s3 );
		if ( cut_rectangles_in_half )
		{	Cell::Positive::Segment * BD = new Cell::Positive::Segment
				( (Cell::Negative::Vertex*) B->reverse_p, D );
			BD->reverse_p = new Cell::Negative::Segment ( tag::reverse_of, BD );
			Cell::Positive * T1 = new Cell::Positive ( tag::triangle, BD->reverse_p, s2, s3 );
			Cell::Positive * T2 = new Cell::Positive ( tag::triangle, BD, s4, s1 );
			T1->add_to (this);
			T2->add_to (this);                                                                  }
		else // with quadrilaterals
		{	Cell::Positive * Q = new Cell::Positive ( tag::rectangle, s1, s2, s3, s4 );
			Q->add_to (this);                                                            }
		*it = s3->reverse_p;
		it_east++;  it_west++;
		it++;  assert ( it == horizon.end() );
	} // end of for i
	it_east++;  it_west++;
	assert ( not it_east.in_range() );
	assert ( not it_west.in_range() );
	// last row of rectangles is different, north sides already exist
	std::list<Cell::Core*>::iterator it = horizon.begin();
	Cell::Core * s4 = west.cell_in_front_of ( NW, tag::surely_exists );
	Cell::Positive::Vertex * D = NW;
	for (size_t j=1; j < N_horiz; j++)
	{	Cell::Core * s1 = *it;
		// s1 == south.cell_in_front_of (A);
	  Cell::Positive::Vertex * B = (Cell::Positive::Vertex*) s1->tip();
		Cell::Core * s3 = north.cell_behind ( D );
		Cell::Positive::Vertex * C = (Cell::Positive::Vertex*) s3->base()->reverse_p;
		Cell::Positive::Segment * s2 = new Cell::Positive::Segment
			( (Cell::Negative::Vertex*) B->reverse_p, C );
		s2->reverse_p = new Cell::Negative::Segment ( tag::reverse_of, s2 );
		if ( cut_rectangles_in_half )
		{	Cell::Positive::Segment * BD = new Cell::Positive::Segment
				( (Cell::Negative::Vertex*) B->reverse_p, D );
			BD->reverse_p = new Cell::Negative::Segment ( tag::reverse_of, BD );
			Cell::Positive * T1 = new Cell::Positive ( tag::triangle, BD->reverse_p, s2, s3 );
			Cell::Positive * T2 = new Cell::Positive ( tag::triangle, BD, s4, s1 );
			T1->add_to (this);
			T2->add_to (this);                                                                  }
		else // with quadrilaterals
		{	Cell::Positive * Q = new Cell::Positive ( tag::rectangle, s1, s2, s3, s4 );
			Q->add_to (this);                                                             }
		it++;
		D = C;
		s4 = s2->reverse_p;                             }
	// and the last rectangle of the last row
	Cell::Core * s1 = *it;
	Cell::Positive::Vertex * B = (Cell::Positive::Vertex*) s1->tip();
	Cell::Core * s2 = east.cell_in_front_of (B);
	Cell::Positive::Vertex * C = (Cell::Positive::Vertex*) s2->tip();
	assert ( C == NE );
	Cell::Core * s3 = north.cell_behind ( D );
	if ( cut_rectangles_in_half )
	{	Cell::Positive::Segment * BD = new Cell::Positive::Segment
			( (Cell::Negative::Vertex*) B->reverse_p, D );
		BD->reverse_p = new Cell::Negative::Segment ( tag::reverse_of, BD );
		Cell::Positive * T1 = new Cell::Positive ( tag::triangle, BD->reverse_p, s2, s3 );
		Cell::Positive * T2 = new Cell::Positive ( tag::triangle, BD, s4, s1 );
		T1->add_to (this);
		T2->add_to (this);                                                                   }
	else // with quadrilaterals
	{	Cell::Positive * Q = new Cell::Positive ( tag::rectangle, s1, s2, s3, s4 );
		Q->add_to (this);                                                             }
	it++;  assert ( it == horizon.end() );

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
	CellIterator it = this->iter_over ( tag::cells_of_dim, 0 );
	it.reset();  assert( it.in_range() );
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
	CellIterator it = this->iter_over ( tag::cells_of_dim, 1 );
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
	CellIterator it = this->iter_over ( tag::cells_of_dim, 0 );
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
	CellIterator it = this->iter_over ( tag::cells_of_dim, 1 );
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


void Mesh::export_msh ( std::string f, std::map<Cell::Core*,size_t> & ver_numbering )

{	// we use the current manifold
	Manifold space = Manifold::working;
	assert ( space.exists() );
	Function coord = space.coordinates();

	std::ofstream file_msh (f);
	file_msh << "$MeshFormat" << std::endl << "2.2 0 8" << std::endl;
	file_msh << "$EndMeshFormat" << std::endl;
	
	file_msh << "$Nodes" << std::endl << this->number_of(tag::cells_of_dim,0) << std::endl;

	{ // just to make variables local : it, counter, x, y
	CellIterator it = this->iter_over ( tag::cells_of_dim, 0 );
	Function x = coord[0], y = coord[1];
	if (coord.nb_of_components() == 2)
	{	for ( it.reset() ; it.in_range(); it++ )
		{	Cell p = *it;
			file_msh << ver_numbering [p.core] << " "
		           << x(p) << " " << y(p) << " " << 0 << std::endl;  }  }
	else
	{	assert  ( coord.nb_of_components() == 3 );
		Function z = coord[2];
		for ( it.reset() ; it.in_range(); it++ )
		{	Cell p = *it;
			file_msh << ver_numbering [p.core] << " "
		           << x(p) << " " << y(p) << " " << z(p) << std::endl;  }  }
	file_msh << "$EndNodes" << std::endl;
	} // just to make variables local : it, counter, x, y

	file_msh << "$Elements" << std::endl;
	file_msh << this->number_of ( tag::cells_of_dim, this->dim() ) << std::endl;

	if ( this->dim() == 1 )
	{	CellIterator it = this->iter_over ( tag::cells_of_dim, 1 );
		size_t counter = 0;
		for ( it.reset() ; it.in_range(); it++)
		{	++counter;
			Cell elem = *it;
			file_msh << counter << " 1 0 ";
			Cell A = elem.base().reverse();
			file_msh << ver_numbering [A.core] << " ";
			Cell B = elem.tip();
			file_msh << ver_numbering [B.core] << std::endl;    }  }
	else if ( this->dim() == 2 )
	{	CellIterator it = this->iter_over ( tag::cells_of_dim, 2 );
		size_t counter = 0;
		for ( it.reset() ; it.in_range(); it++)
		{	++counter;
			Cell elem = *it;
			if ( elem.boundary().number_of ( tag::cells_of_dim, 1 ) == 3 ) // a triangle
				file_msh << counter << " 2 0 ";
			else // a quadrilateral
			{	assert ( elem.boundary().number_of ( tag::cells_of_dim, 1 ) == 4 );
				file_msh << counter << " 3 0 ";                            }
			CellIterator itt = elem.boundary().iter_over ( tag::vertices );
			for ( itt.reset(); itt.in_range(); ++itt )
			{	Cell p = *itt;
				file_msh << ver_numbering [p.core] << " ";   }
			file_msh << std::endl;                                                   }  }
	else
	{	assert ( this->dim() == 3);
		CellIterator it = this->iter_over ( tag::cells_of_dim, 3 );
		size_t counter = 0;
		for ( it.reset() ; it.in_range(); it++)
		{	++counter;
			Cell elem = *it;
			size_t n_faces = elem.boundary().number_of ( tag::cells_of_dim, 2 );
			if ( n_faces == 4 ) // a tetrahedron
				file_msh << counter << " 4 0 ";  // to finish !
			else if ( n_faces == 6 )
			{	// 3d parallelogram = 8-node hexaedron = cube
				// see http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
				// and http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
				file_msh << counter << " 5 0 ";
				CellIterator itt = elem.boundary().iter_over ( tag::cells_of_dim, 2 );
				itt.reset();  Cell back = *itt; // square face behind the cube
				// back is 0321 in gmsh's documentation
				assert ( back.boundary().number_of ( tag::cells_of_dim, 1 ) == 4 );
				CellIterator itv = back.boundary().iter_over ( tag::vertices, tag::reverse );
				// reverse because we want the vertices ordered as 0, 1, 2, 3
				itv.reset();  Cell ver_0 = *itv;
				Cell seg_03 = back.boundary().cell_in_front_of(ver_0);
				for ( ; itv.in_range(); ++itv )
				{	Cell p = *itv;  file_msh << ver_numbering [p.core] << " ";   }
				Cell left_wall = elem.boundary().cell_in_front_of(seg_03); // square face on the left
				// left_wall is 0473 in gmsh's documentation
				assert ( left_wall.boundary().number_of ( tag::cells_of_dim, 1 ) == 4 );
				Cell seg_04 = left_wall.boundary().cell_in_front_of(ver_0);
				Cell ver_4 = seg_04.tip();
				Cell seg_47 = left_wall.boundary().cell_in_front_of(ver_4);
				Cell front = elem.boundary().cell_in_front_of(seg_47); // square face in front
				// front is 4567 in gmsh's documentation
				assert ( front.boundary().number_of ( tag::cells_of_dim, 1 ) == 4 );
				CellIterator itvv = front.boundary().iter_over ( tag::vertices );
				itvv.reset(ver_4);
				for ( ; itvv.in_range(); ++itvv )
				{	Cell p = *itvv;  file_msh << ver_numbering [p.core] << " ";   }                }
			else
			{	assert( n_faces == 5 );
				// triangular prism = 6-node prism
				// see http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
				// and http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
				file_msh << counter << " 6 0 ";
				CellIterator itt = elem.boundary().iter_over ( tag::cells_of_dim, 2 );
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
				CellIterator itv = base.boundary().iter_over ( tag::vertices, tag::reverse );
				// we set it reverse because we want the vertices ordered as 0, 1, 2
				itv.reset();  Cell ver_0 = *itv;
				Cell seg_02 = base.boundary().cell_in_front_of(ver_0);
				for (  ; itv.in_range(); ++itv )
				{	Cell p = *itv;
					file_msh << ver_numbering [p.core] << " ";   }
				Cell right_wall = elem.boundary().cell_in_front_of(seg_02);
				// right_wall is 0352 in gmsh's documentation
				assert ( right_wall.boundary().number_of ( tag::cells_of_dim, 1 ) == 4 );
				Cell seg_03 = right_wall.boundary().cell_in_front_of(ver_0);
				Cell ver_3 = seg_03.tip();
				Cell seg_35 = right_wall.boundary().cell_in_front_of(ver_3);
				Cell roof = elem.boundary().cell_in_front_of(seg_35);
				// roof is 345 in gmsh's documentation
				assert ( roof.boundary().number_of ( tag::cells_of_dim, 1 ) == 3 );
				CellIterator itvv = roof.boundary().iter_over ( tag::vertices );
				itvv.reset(ver_3);
				for ( ; itvv.in_range(); ++itvv )
				{	Cell p = *itvv;  file_msh << ver_numbering [p.core] << " ";  }               }
			file_msh << std::endl;                                                                } }
	file_msh << "$EndElements" << std::endl;
	
	if ( ! file_msh.good() )
	{	std::cerr << "error writing msh file" << std::endl;
		exit (1);                                    }

} // end of Mesh::export_msh


void Mesh::export_msh ( std::string f )
	
// the numbering of vertices is produced on-the-fly

{	std::map < Cell::Core *, size_t > numbering;

	CellIterator it = this->iter_over ( tag::cells_of_dim, 0 );
	size_t counter = 0;
	for ( it.reset() ; it.in_range(); it++ )
	{	++counter;  Cell p = *it;  numbering [p.core] = counter;  }

	this->export_msh ( f, numbering );

} // end of Mesh::export_msh

//----------------------------------------------------------------------------------//


void Mesh::join_list ( const std::list<Mesh> & l )

{	// check the dimensions
	#ifndef NDEBUG
	std::list<Mesh>::const_iterator it0 = l.begin();
	assert ( it0 != l.end() );
	for ( ; it0 != l.end(); it0++ ) assert ( this->dim() == it0->dim() );
	#endif  // DEBUG
	// sweep the list of meshes
	std::list<Mesh>::const_iterator it = l.begin();
	for ( ; it != l.end(); it++ )
	{	Mesh m = *it;  // sweep all cells of m
		// we should implement iterator with tag::cells_of_max_dim, which use cells.back()
		CellIterator itt = m.iter_over ( tag::cells_of_dim, this->dim() );
		for ( itt.reset(); itt.in_range(); itt++ )
		{	Cell cll = *itt;
			cll.add_to ( *this );  }                                           }  }


