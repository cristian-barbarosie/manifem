

#include "maniFEM.h"

using namespace maniFEM;

inline void print_spin ( tag::Util::CompositionOfActions a )

{	std::cout << "{";
	auto it = a.index_map.begin();
	for ( ; it != a.index_map.end(); it++ )
		std::cout << it->first.id << ":" << it->second << ",";
	std::cout << "}" << std::endl;                           }


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

	std::cout << "AB.spin()   ";
	print_spin ( AB.spin() );

	AB.spin() = g1+g2;
	std::cout << "g1+g2   ";
	print_spin ( AB.spin() );
	std::cout << "g1+2g2   ";
	print_spin ( AB.spin() + g2 );

	AB.spin() += g1;
	std::cout << "2g1+g2   ";
	print_spin ( AB.spin() );

	AB.spin() -= AB.spin();
	std::cout << "zero   ";
	print_spin ( AB.spin() );

	std::cout << "end of main" << std::endl << std::flush;
}
