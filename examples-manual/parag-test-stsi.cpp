
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

	std::cout << "main line 31" << std::endl;
	OA .add_to_mesh ( chain );
	std::cout << "main line 33" << std::endl;
	AB .add_to_mesh ( chain );
	std::cout << "main line 35" << std::endl;
	BC .add_to_mesh ( chain );
	std::cout << "main line 37" << std::endl;
	CO .add_to_mesh ( chain );
	std::cout << "main line 39" << std::endl;
	OD .reverse() .add_to_mesh ( chain );
	std::cout << "main line 51" << std::endl;
	DE .reverse() .add_to_mesh ( chain );
	std::cout << "main line 43" << std::endl;
	FO .reverse() .add_to_mesh ( chain );

	std::cout << "main line 46" << std::endl;
	assert ( chain .cell_in_front_of ( A, tag::surely_exists ) == AB );
	std::cout << "main line 48" << std::endl;
	assert ( chain .cell_in_front_of ( B, tag::surely_exists ) == BC );
	std::cout << "main line 50" << std::endl;
	assert ( chain .cell_in_front_of ( O, tag::seen_from, CO, tag::surely_exists ) == OA );
	std::cout << "main line 52" << std::endl;
	assert ( chain .cell_behind ( O, tag::seen_from, FO .reverse(), tag::surely_exists ) == OD .reverse() );
	std::cout << "main line 54" << std::endl;
	assert ( not chain .cell_behind ( E, tag::seen_from, DE .reverse(), tag::may_not_exist ) .exists() );
	std::cout << "main line 56" << std::endl;
	
}
