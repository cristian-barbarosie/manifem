

#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;

int main ()

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0],  y = xy[1];

	// we introduce two equivalence relations on RR2
	Function::Action g1, g2;
	g1 ( xy ) = (x+1.) && y;
	g2 ( xy ) = x && (y+1.);

	// create multi-functions xm, ym (for testing purposes)

	// and divide RR2 by these equivalence relations
	// Manifold torus_manif = RR2.quotient ( g1, g2 );

}
