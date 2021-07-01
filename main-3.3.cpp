
// example presented in paragraph 3.3 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// builds a circle with an eccentric hole

#include "maniFEM.h"
using namespace maniFEM;


int main ( )

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0],  y = xy[1];

	double d = 0.065;

	Manifold circle = RR2.implicit ( x*x + y*y == 1. );
	Mesh outer ( tag::progressive, tag::entire_manifold, circle, tag::desired_length, d );

	double y0 = 0.36;
	Manifold ellipse = RR2.implicit ( x*x + (y-y0)*(y-y0) + 0.3*x*y == 0.25 );
	Mesh inner ( tag::progressive, tag::entire_manifold, ellipse, tag::desired_length, d );

	Mesh circles ( tag::join, outer, inner.reverse() );

	RR2.set_as_working_manifold();
	Mesh disk ( tag::progressive, tag::boundary, circles, tag::desired_length, d );

	disk.export_msh ("disk.msh");
	disk.draw_ps ("disk.eps");

	std::cout << "produced files disk.msh and disk.eps" << std::endl;
}
