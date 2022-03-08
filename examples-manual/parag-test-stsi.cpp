
#include "maniFEM.h"

using namespace maniFEM;

int main ()

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy [0], y = xy [1];

	Cell O ( tag::vertex, tag::of_coords, { 0., 0. } );
	Cell A ( tag::vertex, tag::of_coords, { 1., 0. } );
	Cell B ( tag::vertex, tag::of_coords, { 1., 1. } );
	Cell C ( tag::vertex, tag::of_coords, { 0., 1. } );
	Cell D ( tag::vertex, tag::of_coords, { 0., -1. } );
	Cell E ( tag::vertex, tag::of_coords, { -1., -1. } );
	Cell F ( tag::vertex, tag::of_coords, { -1., 0. } );

	Cell OA ( tag::segment, O .reverse(), A );
	Cell AB ( tag::segment, A .reverse(), B );
	Cell BC ( tag::segment, B .reverse(), C );
	Cell CO ( tag::segment, C .reverse(), O );
	Cell OD ( tag::segment, O .reverse(), D );
	Cell DE ( tag::segment, D .reverse(), E );
	Cell EF ( tag::segment, E .reverse(), F );
	Cell FO ( tag::segment, F .reverse(), O );

	Mesh chain ( tag::stsi, tag::of_dim, 1 );

	OA .add_to_mesh ( chain );
	AB .add_to_mesh ( chain );
	BC .add_to_mesh ( chain );
	CO .add_to_mesh ( chain );
	OD .reverse() .add_to_mesh ( chain );
	DE .reverse() .add_to_mesh ( chain );
	FO .reverse() .add_to_mesh ( chain );

	assert ( chain .cell_in_front_of ( A, tag::surely_exists ) == AB );
	assert ( chain .cell_in_front_of ( B, tag::surely_exists ) == BC );
	assert ( chain .cell_in_front_of ( O, tag::seen_from, CO, tag::surely_exists ) == OA );
	assert ( chain .cell_behind ( O, tag::seen_from, OD, tag::surely_exists ) == FO );
	
}
