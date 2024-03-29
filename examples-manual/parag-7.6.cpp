
// example presented in paragraph 7.6 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// drill a hole in the flat torus RR2/ZZ2

#include "maniFEM.h"
using namespace maniFEM;


int main ( )

{	// begin with the usual two-dimensional space
	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0], y = xy[1];

	// define two actions on RR2 (translations)
	Manifold::Action g_horiz ( tag::transforms, xy, tag::into, (x+1.) && y ),
	                 g_vert  ( tag::transforms, xy, tag::into, x && (y+1.) );

	// and divide RR2 by the group of translations generated by { g_horiz, g_vert }
	Manifold torus_manif = RR2 .quotient ( g_horiz, g_vert );

	// one vertex is enough to start the process
	Cell A ( tag::vertex );  x(A) = 0.02;  y(A) = 0.02;

	// with this vertex, we build two segments
	Mesh seg_horiz ( tag::segment, A .reverse(), A, tag::divided_in, 10, tag::winding, g_horiz ),
	     seg_vert  ( tag::segment, A .reverse(), A, tag::divided_in, 10, tag::winding, g_vert );

	// two segments are enough to define a rectangle
	Mesh torus ( tag::rectangle,
               seg_horiz, seg_vert, seg_horiz.reverse(), seg_vert.reverse(),
	             tag::winding                                                    );
	// the tag::winding provides no specific information,
	// it just warns maniFEM that we are on a quotient manifold
	// and that it must take winding segments into account
	// specific information about winding numbers is included in the two segments

	// we drill a small hole around A
	std::vector < Cell > vec;
	Mesh::Iterator it = torus .iterator ( tag::over_cells, tag::of_dim, 2, tag::around, A );
	for ( it .reset(); it .in_range(); it++ ) vec .push_back ( *it );
	std::vector < Cell > ::iterator itv;
	for ( itv = vec .begin(); itv != vec .end(); itv++ )
	{	Cell sq = *itv;  sq .remove_from_mesh ( torus );  }

	// we can draw an unfolded mesh
	torus .draw_ps ("torus.eps", tag::unfold,
	                 tag::over_region, x*x + 2.*y*y < 3.5 );
	std::cout << "produced file torus.eps - please edit before viewing" << std::endl;

	// or we can build a new, unfolded mesh and subsequently export as msh file
	Mesh unfolded = torus .unfold ( tag::over_region, x*x + 2.*y*y < 3.5 );
	unfolded .export_to_file ( tag::msh, "torus.msh");
	std::cout << "produced file torus.msh" << std::endl;

}  // end of main
