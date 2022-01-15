
// example presented in paragraph 7.9 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// a skew torus made of two triangular patches

#include "maniFEM.h"
using namespace maniFEM;


int main ( )

{	// begin with the usual two-dimensional space
	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy [0], y = xy [1];

	// define two actions on RR2 (translations)
	Manifold::Action g1 ( tag::transforms, xy, tag::into, (x+1.) && y ),
	                 g2 ( tag::transforms, xy, tag::into, (x+0.5) && (y+0.866) );

	// and divide RR2 by the group of translations generated by {g1,g2}
	Manifold torus_manif = RR2 .quotient ( g1, g2 );

	// one vertex is enough to start the process
	Cell A ( tag::vertex );  x (A) = 0. ;  y (A) = 0.;

	// with this vertex, we build three segments
	Mesh seg_horiz ( tag::segment, A .reverse(), A, tag::divided_in, 10, tag::winding, g1 ),
	     seg1      ( tag::segment, A .reverse(), A, tag::divided_in, 10, tag::winding, g2 ),
	     seg2      ( tag::segment, A .reverse(), A, tag::divided_in, 10, tag::winding, g2 - g1 );

	Mesh tri1 ( tag::triangle, seg_horiz, seg2, seg1 .reverse(), tag::winding ),
	     tri2 ( tag::triangle, seg_horiz .reverse(), seg2 .reverse(), seg1, tag::winding );

	Mesh torus ( tag::join, tri1, tri2 );

	std::vector < Cell > vec;
	Mesh::Iterator it = torus .iterator ( tag::over_cells, tag::of_dim, 2, tag::around, A );
	for ( it .reset(); it .in_range(); it++ )  vec .push_back ( *it );
	std::vector < Cell > ::iterator itv;
	for ( itv = vec .begin(); itv != vec .end(); itv++ )
	{	Cell sq = *itv;  sq .remove_from_mesh ( torus );  }

	// it makes no sense to export 'torus' in msh format
	// but we can draw an unfolded mesh
	torus .draw_ps ("torus.eps", tag::unfold,
	                 tag::over_region, -1.6 < x < 1.7, -1.3 < y < 1.3 );

	std::cout << "produced file torus.eps - please edit before viewing" << std::endl;

}  // end of main
