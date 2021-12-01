
// example presented in paragraph 10.2 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// in a mesh of squares, cut one square in two triangles

#include "maniFEM.h"
#include <fstream>
#include <iostream>

using namespace std;
using namespace maniFEM;


int main ( )

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy [0], y = xy [1];
	
	Cell SW ( tag::vertex );  x (SW) = 0.;  y (SW) = 0.;
	Cell SE ( tag::vertex );  x (SE) = 1.;  y (SE) = 0.;
	Cell NE ( tag::vertex );  x (NE) = 1.;  y (NE) = 1.;
	Cell NW ( tag::vertex );  x (NW) = 0.;  y (NW) = 1.;
	
	Mesh south ( tag::segment, SW .reverse(), SE, tag::divided_in, 20 );
	Mesh east  ( tag::segment, SE .reverse(), NE, tag::divided_in, 20 );
	Mesh north ( tag::segment, NE .reverse(), NW, tag::divided_in, 20 );
	Mesh west  ( tag::segment, NW .reverse(), SW, tag::divided_in, 20 );

	Mesh msh ( tag::rectangle, south, east, north, west );

	Cell A = SW;
	Cell AB = south .cell_in_front_of ( A );
	Cell B = AB .tip();
	Cell ABCD = msh .cell_behind ( AB );
	Cell BC = ABCD .boundary() .cell_in_front_of ( B );
	Cell C = BC.tip();
	Cell CD = ABCD .boundary() .cell_in_front_of ( C );
	Cell D = CD .tip();
	Cell DA = ABCD .boundary() .cell_in_front_of ( D );
	assert ( DA .tip() == A );

	CD .cut_from_bdry_of ( ABCD );
	DA .cut_from_bdry_of ( ABCD );
	// at this point, ABCD is a cell whose boundary is incomplete
	// it has only two sides and an opening
	// however, it still is part of 'msh'
	Cell AC ( tag::segment, A .reverse(), C );
	AC .reverse() .glue_on_bdry_of ( ABCD );
	// at this point, and in spite of its name, ABCD is no longer a square
	// it is a triangle, still part of 'msh'
	Cell CDA ( tag::triangle, CD, DA, AC );
	CDA .add_to_mesh ( msh );

	msh .export_msh ( "cut-square.msh" );
	std::cout << "produced file cut-square.msh" << std::endl;
	
}  // end of  main
