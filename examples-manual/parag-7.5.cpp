
// example presented in paragraph 7.5 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// mesh a doughnut

#include "maniFEM.h"
using namespace maniFEM;


int main ( )

{	// begin with the usual two-dimensional space
	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function ab = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function alpha = ab [0], beta = ab [1];
	const double pi = 3.1415926536;

	// define two actions on RR2 (translations)
	Manifold::Action g1 ( tag::transforms, ab, tag::into, (alpha+2.*pi) && beta ),
	                 g2 ( tag::transforms, ab, tag::into, alpha && (beta+2.*pi) );

	// and divide RR2 by the group of translations generated by { g1, g2 }
	Manifold torus_manif = RR2 .quotient ( g1, g2 );

	// one vertex is enough to start the process
	Cell A ( tag::vertex );  alpha (A) = 0.;  beta (A) = 0.;

	// with this vertex, we build two segments
	Mesh seg_horiz ( tag::segment, A .reverse(), A, tag::divided_in, 40, tag::winding, g1 ),
	     seg_vert  ( tag::segment, A .reverse(), A, tag::divided_in, 20, tag::winding, g2 );

	// two segments are enough to define a rectangle
	Mesh torus ( tag::rectangle,
               seg_horiz, seg_vert, seg_horiz.reverse(), seg_vert.reverse(),
	             tag::winding                                                    );
	// the tag::winding provides no specific information,
	// it just warns maniFEM that we are on a quotient manifold
	// and that it must take winding segments into account
	// specific information about winding numbers is included in the two segments

	// parametrize the doughnut :
	const double big_radius = 3., small_radius = 1.;
	// define x, y and z as functions of alpha and beta
	Function x = ( big_radius + small_radius*cos(beta) ) * cos(alpha),
	         y = ( big_radius + small_radius*cos(beta) ) * sin(alpha),
	         z = small_radius*sin(beta);

	// forget about alpha and beta :
	torus_manif .set_coordinates ( x && y && z );
	// in future statements (e.g. for graphical representation)
	// x, y and z will be used, not alpha nor beta :
	torus .export_to_file ( tag::msh, "torus.msh");

	std::cout << "produced file torus.msh" << std::endl;

}  // end of main
