
// example presented in paragraph 7.14 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// folds a square

#include "maniFEM.h"
using namespace maniFEM;

int main ( )

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0], y = xy[1];

	size_t n = 20;
	double d = 2.6 / double(n);

	Cell A ( tag::vertex );  x(A) = -1.3;  y(A) = -1.3;
	Cell B ( tag::vertex );  x(B) =  1.3;  y(B) = -1.3;
	Cell C ( tag::vertex );  x(C) =  1.3;  y(C) =  1.3;
	Cell D ( tag::vertex );  x(D) = -1.3;  y(D) =  1.3;

	Mesh AB ( tag::segment, A.reverse(), B, tag::divided_in, n );
	Mesh BC ( tag::segment, B.reverse(), C, tag::divided_in, n );
	Mesh CD ( tag::segment, C.reverse(), D, tag::divided_in, n );
	Mesh DA ( tag::segment, D.reverse(), A, tag::divided_in, n );

	Manifold circle = RR2.implicit ( x*x + y*y == 1. );
	Mesh inner ( tag::progressive, tag::entire_manifold, circle, tag::desired_length, d );

	Mesh bdry ( tag::join, AB, BC, CD, DA, inner.reverse() );

	RR2.set_as_working_manifold();
	Mesh square ( tag::progressive, tag::boundary, bdry, tag::desired_length, d );

	Mesh cyl = square.fold
		( tag::identify, BC, tag::with, DA.reverse(), tag::use_existing_vertices );

	cyl.draw_ps ( "cylinder.eps", tag::unfold,
                tag::over_region, -2.1 < x < 4.3, -3.6 < y < 2.1 );

	CellIterator it = cyl.iterator ( tag::over_vertices );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell P = *it;
		if ( P.belongs_to ( inner ) ) continue;
		if ( P.belongs_to ( AB ) ) continue;
		if ( P.belongs_to ( CD ) ) continue;
		cyl.baricenter ( P, tag::spin );        }

	cyl.draw_ps ( "cylinder-smooth.eps", tag::unfold,
               tag::over_region, -2.1 < x < 4.3, -3.6 < y < 2.1 );

	std::cout << "produced files cylinder.eps and cylinder-smooth.eps - please edit before viewing" << std::endl;	
}
