

#include "maniFEM.h"
using namespace maniFEM;


int main ( )

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0], y = xy[1];

	size_t n = 10;
	double l = 2.6;

	Cell A ( tag::vertex );  x(A) = -l/2.;  y(A) = -l/2.;
	Cell B ( tag::vertex );  x(B) =  l/2.;  y(B) = -l/2.;
	Cell C ( tag::vertex );  x(C) =  l/2.;  y(C) =  l/2.;
	Cell D ( tag::vertex );  x(D) = -l/2.;  y(D) =  l/2.;

	Mesh AB ( tag::segment, A.reverse(), B, tag::divided_in, n );
	Mesh BC ( tag::segment, B.reverse(), C, tag::divided_in, n );
	Mesh CD ( tag::segment, C.reverse(), D, tag::divided_in, n );
	Mesh DA ( tag::segment, D.reverse(), A, tag::divided_in, n );

	Mesh square ( tag::rectangle, AB, BC, CD, DA );

	Mesh torus = square.fold ( tag::identify, AB, tag::with, CD.reverse(),
	                           tag::identify, BC, tag::with, DA.reverse(),
	                           tag::use_existing_vertices                 );

	std::list < Cell > li;
	CellIterator it = torus .iterator ( tag::over_cells_of_dim, 2,
																			tag::around, B );
	for ( it .reset(); it .in_range(); it ++ ) li .push_back ( * it );

	for ( std::list < Cell > :: iterator itt = li .begin();
				itt != li .end(); itt ++ )
		( * itt ) .remove_from_mesh ( torus );	
	
	torus.draw_ps ( "torus.eps", tag::unfold,
	                tag::over_region, -3. < x < 3.7, -3.4 < y < 3.5 );

	std::cout << "produced file torus.eps - please edit before viewing" << std::endl;

	return 0;
}
