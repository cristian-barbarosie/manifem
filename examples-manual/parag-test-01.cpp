

#include "maniFEM.h"
using namespace maniFEM;


int main ( )

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0], y = xy[1];

	const double pi = 3.1415926536;
	const double theta = pi/2.,  // radians
		cos_theta = std::cos ( theta ), sin_theta = std::sin ( theta );
	
	// define an action on RR2 (a rotation)
	Manifold::Action rot ( tag::transforms, xy, tag::into,
		( cos_theta * x - sin_theta * y ) && ( sin_theta * x + cos_theta * y ) );

	// and divide RR2 by the group generated by rot
	Manifold torus_manif = RR2.quotient ( rot );

	Cell O ( tag::vertex );  x(O) = 0.;  y(O) = 0.;
	Cell A ( tag::vertex );  x(A) = 1.;  y(A) = 0.;
	Cell B ( tag::vertex );  x(B) = 1.;  y(B) = 1.;

	Mesh OA ( tag::segment, O.reverse(), A, tag::divided_in, 10 );  // no winding
	Mesh OB ( tag::segment, O.reverse(), B, tag::divided_in, 10 );  // no winding
	Mesh AB ( tag::segment, A.reverse(), B, tag::divided_in, 10 );  // no winding
	Mesh BA ( tag::segment, B.reverse(), A, tag::divided_in, 10, tag::winding, rot );

	Mesh tri_1 ( tag::triangle, OA, AB, OB.reverse() );  // no winding
	Mesh tri_2 ( tag::triangle, OB, BA, OA.reverse(), tag::winding, tag::singular, O );
	// on the last triangle, the one having O as vertex, the sum of windings is not zero

	Mesh sector ( tag::join, tri_1, tri_2 );
	
	sector.draw_ps ( "sector.eps", tag::unfold, { -rot, 0, rot },
	                 tag::over_region, -2.1 < x < 2.1, -0.3 < y < 2.1 );

	std::cout << "produced file sector.eps - please edit before viewing" << std::endl;
}
