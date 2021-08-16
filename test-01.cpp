

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
	Function xm = x.multivalued ( tag::through, g1, tag::becomes, x+1.5,
	                              tag::through, g2, tag::becomes, x-0.2 );

	Cell V ( tag::vertex );  x(V) = 1.1;  y(V) = 22.2;

	std::cout << xm(V) << std::endl;
	std::cout << "00 " << xm ( V, tag::spin, {0,0} ) << std::endl;
	std::cout << "01 " << xm ( V, tag::spin, {0,1} ) << std::endl;
	std::cout << "10 " << xm ( V, tag::spin, {1,0} ) << std::endl;
	std::cout << "11 " << xm ( V, tag::spin, {1,1} ) << std::endl;
	std::cout << "12 " << xm ( V, tag::spin, {1,2} ) << std::endl;
	
	// and divide RR2 by these equivalence relations
	// Manifold torus_manif = RR2.quotient ( g1, g2 );

}
