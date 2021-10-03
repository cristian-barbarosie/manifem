

#include "maniFEM.h"

using namespace maniFEM;


	
int main ( )

{	// begin with the usual two-dimensional space
	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0], y = xy[1];

	const double e = 1.5;
	Manifold curve = RR2.implicit 
		( smooth_min ( 100.*power((x+y)*(x+y),e) + power((x-y-1.)*(x-y-1.),e),
		               100.*power((x-y)*(x-y),e) + power((x+y+1.)*(x+y+1.),e),
	                 tag::threshold, 2.                     )  == 1. );

	Mesh loop ( tag::progressive, tag::desired_length, 0.04 );

	// loop.export_msh ("loop.msh");
	loop.draw_ps ("loop.eps");
}

