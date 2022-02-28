
#include "maniFEM.h"

using namespace maniFEM;

int main ()

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy [0], y = xy [1];

	Mesh::STSI * msh = new Mesh::STSI ( tag::of_dim, 3, tag::minus_one, tag::one_dummy_wrapper );

}
