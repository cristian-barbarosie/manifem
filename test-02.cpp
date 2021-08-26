

#include "maniFEM.h"

using namespace maniFEM;

int main ()

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0],  y = xy[1];

	// we introduce two equivalence relations on RR2
	Function::Action g1 ( tag::transforms, xy, tag::into, (x+3.) && (y+0.1) );
	Function::Action g2 ( tag::transforms, xy, tag::into, (x-0.1) && (y+1.) );

	Manifold torus = RR2.quotient ( g1, g2 );
	
	Cell A ( tag::vertex ), B ( tag::vertex );
	Cell AB ( tag::segment, A.reverse(), B );
	AB.spin() = g1;
	std::cout << "AB.spin() == 0  ";
	if ( AB.spin() == 0 ) std::cout << "true" << std::endl;
	else std::cout << " false" << std::endl;
	std::cout << "AB.spin() == g1  ";
	if ( AB.spin() == g1 ) std::cout << "true" << std::endl;
	else std::cout << " false" << std::endl;
	std::cout << "AB.spin() == g2  ";
	if ( AB.spin() == g2 ) std::cout << "true" << std::endl;
	else std::cout << " false" << std::endl;
	std::cout << "AB.spin() == g1+g2  ";
	if ( AB.spin() == g1+g2 ) std::cout << "true" << std::endl;
	else std::cout << " false" << std::endl;


	std::cout << "end of main" << std::endl << std::flush;
}
