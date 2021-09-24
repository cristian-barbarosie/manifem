
// example presented in paragraph 7.14 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// folds a square

#include "maniFEM.h"
using namespace maniFEM;

int main ( )

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0],  y = xy[1];

	size_t n = 30;
	double d = 2. / double(n);

	Cell A ( tag::vertex );  x(A) = -1.2;  y(A) = -1.3;
	Cell B ( tag::vertex );  x(B) =  1.2;  y(B) = -1.3;
	Cell C ( tag::vertex );  x(C) =  1.2;  y(C) =  1.3;
	Cell D ( tag::vertex );  x(D) = -1.2;  y(D) =  1.3;

	Mesh AB ( tag::segment, A.reverse(), B, tag::divided_in, n );
	Mesh BC ( tag::segment, B.reverse(), C, tag::divided_in, n );
	Mesh CD ( tag::segment, C.reverse(), D, tag::divided_in, n );
	Mesh DA ( tag::segment, D.reverse(), A, tag::divided_in, n );

	Manifold circle = RR2.implicit ( x*x + y*y == 1. );
	Mesh inner ( tag::progressive, tag::entire_manifold, circle, tag::desired_length, d );

	Mesh bdry ( tag::join, AB, BC, CD, DA, inner.reverse() );

	RR2.set_as_working_manifold();
	Mesh square ( tag::progressive, tag::boundary, bdry, tag::desired_length, d );

	Mesh cyl = square.fold ( tag::identify, AB, tag::with, CD, tag::build_new_vertices );
}
