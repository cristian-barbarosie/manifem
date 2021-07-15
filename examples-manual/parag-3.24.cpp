
// example presented in paragraph 3.24 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// a disk with an eccentric hole, anisotropic mesh
// the code shown in the manual does not work (yet)
// we fake the desired result by meshing a surface in 3D

#include "maniFEM.h"
using namespace maniFEM;


int main ( )

{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz[0],  y = xyz[1],  z = xyz[2];

	Manifold parab_surf = RR3.implicit ( z == 2. * smooth_max ( 3.*y - x - 3.2 , 0., tag::threshold, 1. ) );

	double d = 0.095;
	Manifold circle = parab_surf.implicit ( x*x + y*y == 1. );
	Cell P ( tag::vertex );  x(P) = 1.;  y(P) = 0.;  z(P) = 0.;
	circle.project ( P );
	Mesh outer ( tag::progressive, tag::start_at, P, tag::desired_length, d, tag::random_orientation );

	double y0 = 0.37;
	Manifold ellipse = parab_surf.implicit ( x*x + (y-y0)*(y-y0) + 0.3*x*y == 0.25 );
	Mesh inner ( tag::progressive, tag::desired_length, d, tag::random_orientation );

	Mesh circles ( tag::join, outer, inner.reverse() );

	parab_surf.set_as_working_manifold();
	Mesh disk ( tag::progressive, tag::boundary, circles,
              tag::start_at, P, tag::towards, { -1., 0., 0. },
              tag::desired_length, d                           );

	disk.export_msh ("disk.msh");
	RR3.set_as_working_manifold();
	RR3.set_coordinates ( x && y );
	disk.draw_ps ("disk.eps");

	std::cout << "produced files disk.eps and disk.msh" << std::endl;
}
