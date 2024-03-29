
// example presented in paragraph 7.7 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// a skew flat torus

#include "maniFEM.h"
using namespace maniFEM;


int main ( )

{	// begin with the usual two-dimensional space
	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0], y = xy[1];

	// define two actions on RR2 (translations)
	Manifold::Action g1 ( tag::transforms, xy, tag::into, (x+1.) && (y+0.1) ),
	                 g2 ( tag::transforms, xy, tag::into, (x+0.1) && (y+1.) );

	// and divide RR2 by the group of translations generated by {g1,g2}
	Manifold torus_manif = RR2 .quotient ( g1, g2 );

	// one vertex is enough to start the process
	Cell A ( tag::vertex );  x(A) = 0.02;  y(A) = 0.02;

	// with this vertex, we build two segments
	Mesh seg_horiz ( tag::segment, A.reverse(), A, tag::divided_in, 10, tag::winding, g1 ),
	     seg_vert  ( tag::segment, A.reverse(), A, tag::divided_in, 10, tag::winding, g2 );

	// two segments are enough to define a rectangle
	Mesh torus ( tag::rectangle,
               seg_horiz, seg_vert, seg_horiz.reverse(), seg_vert.reverse(),
	             tag::winding                                                 );

	std::vector < Cell > vec;
	Mesh::Iterator it = torus .iterator ( tag::over_cells, tag::of_dim, 2, tag::around, A );
	for ( it .reset(); it .in_range(); it++ )  vec .push_back ( *it );
	std::vector < Cell > ::iterator itv;
	for ( itv = vec .begin(); itv != vec .end(); itv++ )
	{	Cell sq = *itv;  sq .remove_from_mesh ( torus );  }

	// we can draw an unfolded mesh
	torus .draw_ps ("torus.eps", tag::unfold,
	                 tag::over_region, -1.2 < x < 1.7, -0.3 < y < 1.5 );

	std::cout << "produced file torus.eps - please edit before viewing" << std::endl;

}  // end of main
