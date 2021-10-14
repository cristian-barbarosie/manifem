
// solve a celullar problem


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

	Manifold circle = RR2.implicit ( x*x + y*y == 0.7 );
	Mesh inner ( tag::progressive, tag::entire_manifold, circle, tag::desired_length, d );

	Mesh bdry ( tag::join, AB, BC, CD, DA, inner.reverse() );

	RR2.set_as_working_manifold();
	Mesh square ( tag::progressive, tag::boundary, bdry, tag::desired_length, d );

	Mesh torus = square.fold ( tag::identify, AB, tag::with, CD.reverse(),
	                           tag::identify, BC, tag::with, DA.reverse(),
	                           tag::use_existing_vertices                 );

	std::cout << "produced folded mesh, now drawing, please wait" << std::endl << std::flush;
	
	torus.draw_ps ( "torus.eps", tag::unfold,
                  tag::over_region, -2.1 < x < 4.3, -3.6 < y < 2.1 );

	std::cout << "now smoothening ... " << std::flush;

	CellIterator it = torus.iterator ( tag::over_vertices );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell P = *it;
		if ( P.belongs_to ( inner ) ) continue;
		torus.baricenter ( P, tag::spin );        }

	std::cout << "and drawing again, please wait" << std::endl << std::flush;

	torus.draw_ps ( "torus-smooth.eps", tag::unfold,
                  tag::over_region, -2.1 < x < 4.3, -3.6 < y < 2.1 );

	std::cout << "produced files torus.eps and torus-smooth.eps - please edit before viewing" << std::endl;	
}
