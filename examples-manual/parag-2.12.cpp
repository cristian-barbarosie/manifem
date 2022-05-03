
// example presented in paragraph 2.12 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// using zero-dimensional manifolds to compute the intersection of two parabolas

#include "maniFEM.h"

using namespace maniFEM;

int main ()

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0], y = xy[1];

	Manifold parab_vert  = RR2 .implicit ( y == x * x - 2. );
	Manifold parab_horiz = RR2 .implicit ( x == y * y - 2. );

	Manifold four_points ( tag::intersect, parab_vert, parab_horiz );
	
	Cell A ( tag::vertex, tag::of_coords, {-1.,-1.}, tag::project );
	Cell B ( tag::vertex, tag::of_coords, { 1.,-1.}, tag::project );
	Cell C ( tag::vertex, tag::of_coords, { 1., 1.}, tag::project );
	Cell D ( tag::vertex, tag::of_coords, {-1., 1.}, tag::project );

	// four arcs of parabola :
	parab_vert .set_as_working_manifold();
	Mesh BC ( tag::segment, B .reverse(), C, tag::divided_in, 12 );
	Mesh DA ( tag::segment, D .reverse(), A, tag::divided_in, 12 );
	parab_horiz .set_as_working_manifold();
	Mesh AB ( tag::segment, A .reverse(), B, tag::divided_in, 10 );
	Mesh CD ( tag::segment, C .reverse(), D, tag::divided_in, 10 );

	RR2 .set_as_working_manifold();

	Mesh ABCD ( tag::rectangle, AB, BC, CD, DA );

	ABCD .export_to_file ( tag::msh, "rectangle.msh");
	
	std::cout << "produced file rectangle.msh" << std::endl;
	
}  // end of main
