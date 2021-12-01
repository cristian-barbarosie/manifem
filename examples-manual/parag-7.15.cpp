
// example presented in paragraph 7.15 of the manual
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
	
	// define an action on RR2 (a rotation)
	Manifold::Action rot ( tag::transforms, xy, tag::into,
		( cos_theta * x - sin_theta * y ) && ( sin_theta * x + cos_theta * y ) );
	// and a zoom
	Manifold::Action zoom ( tag::transforms, xy, tag::into, ( 2. * x ) && ( 2. * y ) );

	// and divide RR2 by the group generated by { rot, zoom }
	Manifold torus_manif = RR2 .quotient ( rot, zoom );

	Cell A ( tag::vertex );  x (A) = 0.;  y (A) = 1.;

	Mesh AA ( tag::segment, A .reverse(), A, tag::divided_in, 10, tag::winding, -rot );
	Mesh AA_vert ( tag::segment, A .reverse(), A, tag::divided_in, 7, tag::winding, zoom );

	Mesh sector ( tag::rectangle, AA, AA_vert, AA .reverse(), AA_vert .reverse(), tag::winding );

	// it makes no sense to export 'torus' in msh format
	// but we can draw an unfolded mesh
	sector .draw_ps ("sector.eps", tag::unfold,
	                  { -2*rot, -rot, 0, rot, 2*rot, zoom-2*rot, zoom-rot, zoom, zoom+rot,
	                    zoom+2*rot, -zoom-2*rot, -zoom-rot, -zoom, -zoom+rot, -zoom+2*rot,
	                    -2*zoom-2*rot, -2*zoom-rot, -2*zoom, -2*zoom+rot, -2*zoom+2*rot,
	                    -3*zoom-2*rot, -3*zoom-rot, -3*zoom, -3*zoom+rot, -3*zoom+2*rot   },
	                  tag::over_region, -2.1 < x < 2.1, -0.3 < y < 2.1                      );

	std::cout << "produced file sector.eps - please edit before viewing" << std::endl;
}
