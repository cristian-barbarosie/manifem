
// example presented in paragraph 7.17 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// a sector of an annulus, repeated

#include "maniFEM.h"
using namespace maniFEM;


int main ( )

{	// begin with the usual two-dimensional space
	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy [0], y = xy [1];

	const double theta = 1.,  // radians
		cos_theta = std::cos ( theta ), sin_theta = std::sin ( theta );
	
	// define an action on RR2 (a linear map)
	Manifold::Action g ( tag::transforms, xy, tag::into,
		( 0.8 * cos_theta * x - 0.8 * sin_theta * y ) && ( 0.8 * sin_theta * x + 0.8 * cos_theta * y ) );

	// and divide RR2 by the group generated by g
	Manifold torus_manif = RR2 .quotient ( g );

	Cell A ( tag::vertex );  x (A) = 0.;  y (A) = 1.;
	Cell B ( tag::vertex );  x (B) = 0.;  y (B) = 2.;

	Mesh AA ( tag::segment, A .reverse(), A, tag::divided_in, 10, tag::winding, -g );
	Mesh AB ( tag::segment, A .reverse(), B, tag::divided_in, 7 );  // no winding
	Mesh BB ( tag::segment, B .reverse(), B, tag::divided_in, 10, tag::winding, g );

	Mesh sector ( tag::rectangle, AA, AB, BB, AB.reverse(), tag::winding );

	// we can draw an unfolded mesh
	sector .draw_ps ("sector.eps", tag::unfold,
	                  { -2*g, -g, 0, g, 2*g, 3*g, 4*g, 5*g, 6*g,
	                    7*g, 8*g, 9*g, 10*g, 11*g, 12*g, 13*g,
	                    14*g, 15*g, 16*g, 17*g, 18*g, 19*g      },
	                  tag::over_region, -2.1 < x < 2.3, -0.6 < y < 2.1 );
	std::cout << "produced file sector.eps - please edit before viewing" << std::endl;

	// or we can build a new, unfolded mesh and subsequently export as msh file
	Mesh unf = sector .unfold ( { -2*g, -g, 0, g, 2*g, 3*g, 4*g, 5*g, 6*g,
	                              7*g, 8*g, 9*g, 10*g, 11*g, 12*g, 13*g,
	                              14*g, 15*g, 16*g, 17*g, 18*g, 19*g      },
	                            tag::over_region, -2.1 < x < 2.3, -0.6 < y < 2.1 );
	unf .export_to_file ( tag::msh, "sector.msh");
	std::cout << "produced file sector.msh" << std::endl;

}  // end of main
