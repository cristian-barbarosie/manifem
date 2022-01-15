
// example presented in paragraph 7.13 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// a cone-like shape

#include "maniFEM.h"
using namespace maniFEM;


int main ( )

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy [0], y = xy [1];

	const double pi = 3.1415926536;
	const double theta = pi/3.,  // radians
		cos_theta = std::cos ( theta ), sin_theta = std::sin ( theta );
	
	// define an action on RR2 (a rotation)
	Manifold::Action rot ( tag::transforms, xy, tag::into,
		( cos_theta * x - sin_theta * y ) && ( sin_theta * x + cos_theta * y ) );

	// and divide RR2 by the group generated by rot
	Manifold torus_manif = RR2 .quotient ( rot );

	Cell O ( tag::vertex );  x (O) = 0.;  y (O) = 0.;
	Cell A ( tag::vertex );  x (A) = 1.;  y (A) = 0.;

	Mesh OA ( tag::segment, O .reverse(), A, tag::divided_in, 10 );  // no winding
	Mesh AA ( tag::segment, A .reverse(), A, tag::divided_in, 10, tag::winding, rot );

	Mesh sector ( tag::triangle, OA, AA, OA .reverse(), tag::winding, tag::singular, O );
	// on the last triangle, the one having O as vertex, the sum of windings is not zero

	sector .draw_ps ("sector.eps", tag::unfold, { -rot, 0, rot },
	                  tag::over_region, -2.1 < x < 2.1, -0.3 < y < 2.1 );

	std::cout << "produced file sector.eps - please edit before viewing" << std::endl;

}  // end of main
