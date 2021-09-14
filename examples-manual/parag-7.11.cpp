
// example presented in paragraph 7.11 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// a sector of an annulus, repeated

#include "maniFEM.h"
using namespace maniFEM;


int main ( )

{	// begin with the usual two-dimensional space
	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0], y = xy[1];

	const double theta = 1.,  // radians
		cos_theta = std::cos ( theta ), sin_theta = std::sin ( theta );
	
	// define an action on RR2 (a rotation)
	Function::Action g ( tag::transforms, xy, tag::into,
		( cos_theta * x - sin_theta *y ) && ( sin_theta * x + cos_theta * y ) );

	// and divide RR2 by the group generated by g
	Manifold torus_manif = RR2.quotient ( g );

	Cell A ( tag::vertex );  x(A) = 0.;  y(A) = 1.;
	Cell B ( tag::vertex );  x(B) = 0.;  y(B) = 2.;

	Mesh AA ( tag::segment, A.reverse(), A, tag::divided_in, 10, tag::spin, -g );
	Mesh AB ( tag::segment, A.reverse(), B, tag::divided_in, 10 );  // no jump
	Mesh BB ( tag::segment, B.reverse(), B, tag::divided_in, 10, tag::spin, g );

	Mesh cylinder ( tag::rectangle, AA, AB, BB, AB.reverse(), tag::spin );

	// it makes no sense to export 'torus' in msh format
	// but we can draw an unfolded mesh
	torus.draw_ps ( "sector.eps", tag::unfold, { -2*g, -g, 0, g },
	                tag::over_region, -2.1 < x < 2.1, -0.3 < y < 2.1 );

	std::cout << "produced file torus.eps - please edit before viewing" << std::endl;
}
