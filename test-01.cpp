

#include "maniFEM.h"

using namespace maniFEM;

int main () 

{	Manifold spiral ( tag::Euclid, tag::of_dim, 1 );
	Function t = spiral.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	const double pi = 4.*atan(1.);

	Cell::Numbering::Field num ( tag::vertices );
	
	Cell A ( tag::vertex );  t(A) = pi/2.;
	Cell B ( tag::vertex );  t(B) = 5.*pi;
	Mesh arc_of_spiral ( tag::segment, A.reverse(), B, tag::divided_in, 50 );
	 
	CellIterator it = arc_of_spiral.iterator ( tag::over_vertices );
//	for ( it.reset(); it.in_range(); it++ )
//	{	Cell P = *it;
//		std::cout << "vertex number " << num(P) << " has coords " << x(P) << " " << y(P) << std::endl; }

	std::cout << "so we have " << num.size() << " labels, but only "
						<< arc_of_spiral.number_of ( tag::vertices ) << " vertices" << std::endl;
}
