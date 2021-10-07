
// example presented in paragraph 7.16 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// folds a six-sided figure

#include "maniFEM.h"
using namespace maniFEM;

int main ( )

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0], y = xy[1];

	size_t n = 10;
	double d = 1.3 / double(n);

	Cell A ( tag::vertex );  x(A) = -1.3;  y(A) = -1.3;
	Cell B ( tag::vertex );  x(B) =  0. ;  y(B) = -1.3;
	Cell C ( tag::vertex );  x(C) =  1.3;  y(C) = -1.3;
	Cell D ( tag::vertex );  x(D) =  1.3;  y(D) =  1.3;
	Cell E ( tag::vertex );  x(E) =  0. ;  y(E) =  1.3;
	Cell F ( tag::vertex );  x(F) = -1.3;  y(F) =  1.3;

	Mesh AB ( tag::segment, A.reverse(), B, tag::divided_in,   n );
	Mesh BC ( tag::segment, B.reverse(), C, tag::divided_in,   n );
	Mesh CD ( tag::segment, C.reverse(), D, tag::divided_in, 2*n );
	Mesh DE ( tag::segment, D.reverse(), E, tag::divided_in,   n );
	Mesh EF ( tag::segment, E.reverse(), F, tag::divided_in,   n );
	Mesh FA ( tag::segment, F.reverse(), A, tag::divided_in, 2*n );

	Manifold circle = RR2.implicit ( x*x + y*y == 1. );
	Mesh inner ( tag::progressive, tag::entire_manifold, circle, tag::desired_length, d );

	Mesh bdry ( tag::join, { AB, BC, CD, DE, EF, FA, inner.reverse() } );

	RR2.set_as_working_manifold();
	Mesh square ( tag::progressive, tag::boundary, bdry, tag::desired_length, d );

	Mesh torus = square.fold ( tag::identify, CD, tag::with, FA.reverse(),
	                           tag::identify, BC, tag::with, EF.reverse(),
	                           tag::identify, AB, tag::with, DE.reverse(),
	                           tag::use_existing_vertices                 );

	std::cout << "produced folded mesh, now drawing, please wait" << std::endl << std::flush;
	
	torus.draw_ps ( "torus.eps", tag::unfold,
                  tag::over_region, -2.1 < x < 4.3, -2.6 < y < 2.1 );

	std::cout << "produced file torus.eps - please edit before viewing" << std::endl;	
}
