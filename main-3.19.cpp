
// example presented in paragraph 3.19 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// a sphere with a cylinder-shaped tunnel

#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;
using namespace std;


int main ()

{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz[0],  y = xyz[1],  z = xyz[2];

	const double rs = 1.;    // radius of the sphere
	const double rc = 0.45;  // radius of the cylinder
	const double seg_size = 0.1;
	
	Manifold cylinder = RR3.implicit ( y*y + (z-0.5)*(z-0.5) == rc*rc );
	Manifold intersection = cylinder.implicit ( x*x + y*y + z*z == rs*rs );
	Cell start1 ( tag::vertex );  x(start1) = 1.;  y(start1) = 0.;  z(start1) = 0.5 - rc;
	intersection.project ( start1 );
	Mesh circle_1 ( tag::progressive, tag::start_at, start1, tag::towards, { 0., 1., 0. },
	                tag::desired_length, seg_size                                          );

	Cell start2 ( tag::vertex );  x(start2) = -1.;  y(start2) = 0.;  z(start2) = 0.5 - rc;
	intersection.project ( start2 );
	Mesh circle_2 ( tag::progressive, tag::start_at, start2, tag::towards, { 0., 1., 0. },
                  tag::desired_length, seg_size                                          );

	Mesh two_circles ( tag::join, circle_1.reverse(), circle_2 );
	cylinder.set_as_working_manifold();
	Mesh cyl ( tag::progressive, tag::boundary, two_circles, tag::start_at, start1,
						 tag::towards, { -1., 0., 0. }, tag::desired_length, seg_size         );

	Mesh two_circles_rev ( tag::join, circle_1, circle_2.reverse() );
	RR3.implicit ( x*x + y*y + z*z == rs*rs );
	Mesh sph ( tag::progressive, tag::boundary, two_circles_rev, tag::start_at, start1,
             tag::towards, { 0., 0., -1. }, tag::desired_length, seg_size             );

	Mesh all ( tag::join, cyl, sph );
	all.export_msh ("sphere-tunnel.msh");
	cout << "produced file sphere-tunnel.msh" << endl;
}
