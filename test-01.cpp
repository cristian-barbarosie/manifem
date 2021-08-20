

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

	// and divide RR2 by these equivalence relations
	// Manifold torus_manif = RR2.quotient ( g1, g2 );

	// create a multi-valued function xm (for testing purposes)
	// Function xm ( tag::Lagrange, tag::of_degree, 1, tag::multivalued );
	// describe the action of g1 and g2 on xm :
	// g1 ( xm ) = xm + 1.5;
	// g2 ( xm ) = xm - 0.2;
	// xm.property ( tag::through, g1, tag::becomes, xm + 1.5 );
	// xm.property ( tag::through, g2, tag::becomes, xm - 0.2 );

	Cell V ( tag::vertex );  x(V) = 1.1;  y(V) = 22.2;

	// std::cout << xm(V) << std::endl << std::flush;
	// std::cout << "0 0 " << std::endl << std::flush;
	// std::cout << xm ( V, tag::spin, {0,0} ) << std::endl << std::flush;
	// std::cout << "1 0 " << std::endl << std::flush;
	// std::cout << xm ( V, tag::spin, {1,0} ) << std::endl << std::flush;
	// std::cout << "1 1 " << std::endl << std::flush;
	// std::cout << xm ( V, tag::spin, {1,1} ) << std::endl << std::flush;
	// std::cout << "1 2 " << std::endl << std::flush;
	// std::cout << xm ( V, tag::spin, {1,2} ) << std::endl << std::flush;
	
}
