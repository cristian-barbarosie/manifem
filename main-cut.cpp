
#include "maniFEM.h"

using namespace maniFEM;

int main ()
	
{	Manifold RR2 ( tag::Euclid, tag::of_dim, 3 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0],  y = xy[1];

	
