

#include "maniFEM.h"

using namespace maniFEM;

int main () 

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0],  y = xy[1];
	Function d = 0.1/(x*x+y*y+0.5);

	Manifold circle_manif = RR2.implicit ( x*x + 2*y*y == 1. );
	Mesh circle ( tag::progressive, tag::entire_manifold, circle_manif, tag::desired_length, d );
	// Mesh circle ( tag::progressive, tag::desired_length, 0.2 );

	RR2.set_as_working_manifold();
	Mesh disk ( tag::progressive, tag::boundary, circle, tag::desired_length, d );

	disk.export_msh ("disk.msh");
	disk.draw_ps ("disk.eps");
	std::cout << "produced files disk.eps and disk.msh" << std::endl;
}
