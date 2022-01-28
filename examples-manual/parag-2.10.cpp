
// example presented in paragraph 2.10 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// half of a disk

#include "maniFEM.h"

using namespace maniFEM;

int main ()

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy [0], y = xy [1];

	Cell O ( tag::vertex );  x (O) =  0.;   y (O) = 0.;
	Cell A ( tag::vertex );  x (A) =  1.;   y (A) = 0.;
	Cell B ( tag::vertex );  x (B) =  0.5;  y (B) = 0.8;
	Cell C ( tag::vertex );  x (C) = -0.5;  y (C) = 0.8;
	Cell D ( tag::vertex );  x (D) = -1.;   y (D) = 0.;

	Manifold circle = RR2 .implicit ( x*x + y*y == 1. );
	circle .project (B);  circle .project (C);

	// three arcs of circle :
	Mesh AB ( tag::segment, A .reverse(), B, tag::divided_in, 10 );
	Mesh BC ( tag::segment, B .reverse(), C, tag::divided_in, 10 );
	Mesh CD ( tag::segment, C .reverse(), D, tag::divided_in, 10 );

	RR2 .set_as_working_manifold();

	Mesh OA ( tag::segment, O .reverse(), A, tag::divided_in, 10 );
	Mesh OB ( tag::segment, O .reverse(), B, tag::divided_in, 10 );
	Mesh OC ( tag::segment, O .reverse(), C, tag::divided_in, 10 );
	Mesh OD ( tag::segment, O .reverse(), D, tag::divided_in, 10 );

	Mesh OAB ( tag::triangle, OA, AB, OB .reverse() );
	Mesh OBC ( tag::triangle, OB, BC, OC .reverse() );
	Mesh OCD ( tag::triangle, OC, CD, OD .reverse() );
	Mesh half_disk ( tag::join, OAB, OBC, OCD );

	half_disk .export_to_file ( tag::msh, "half-disk.msh");
	
	std::cout << "produced file half-disk.msh" << std::endl;
	
}  // end of main
