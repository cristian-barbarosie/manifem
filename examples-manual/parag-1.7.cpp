
// example presented in paragraph 1.7 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// meshing a cube

#include "maniFEM.h"
#include <fstream>

using namespace maniFEM;

int main ()

{	// we choose our (geometric) space dimension :
	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	
	// xyz is a map defined on our future mesh with values in RR2 :
	Function xyz = RR3 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );

	// we can extract components of xyz using the [] operator :
	Function x = xyz [0], y = xyz [1], z = xyz [2];

	Cell A ( tag::vertex );  x (A) = 0.;  y (A) = 0.;  z (A) = 0.;
	Cell B ( tag::vertex );  x (B) = 1.;  y (B) = 0.;  z (B) = 0.;
	Cell C ( tag::vertex );  x (C) = 1.;  y (C) = 0.;  z (C) = 1.;
	Cell D ( tag::vertex );  x (D) = 0.;  y (D) = 0.;  z (D) = 1.;
	Cell E ( tag::vertex );  x (E) = 0.;  y (E) = 1.;  z (E) = 0.;
	Cell F ( tag::vertex );  x (F) = 1.;  y (F) = 1.;  z (F) = 0.;
	Cell G ( tag::vertex );  x (G) = 1.;  y (G) = 1.;  z (G) = 1.;
	Cell H ( tag::vertex );  x (H) = 0.;  y (H) = 1.;  z (H) = 1.;
	Mesh AB ( tag::segment, A .reverse(), B, tag::divided_in, 10 );
	Mesh BC ( tag::segment, B .reverse(), C, tag::divided_in, 10 );
	Mesh CD ( tag::segment, C .reverse(), D, tag::divided_in, 10 );
	Mesh DA ( tag::segment, D .reverse(), A, tag::divided_in, 10 );
	Mesh EF ( tag::segment, E .reverse(), F, tag::divided_in, 10 );
	Mesh FG ( tag::segment, F .reverse(), G, tag::divided_in, 10 );
	Mesh GH ( tag::segment, G .reverse(), H, tag::divided_in, 10 );
	Mesh HE ( tag::segment, H .reverse(), E, tag::divided_in, 10 );
	Mesh AE ( tag::segment, A .reverse(), E, tag::divided_in, 10 );
	Mesh BF ( tag::segment, B .reverse(), F, tag::divided_in, 10 );
	Mesh CG ( tag::segment, C .reverse(), G, tag::divided_in, 10 );
	Mesh DH ( tag::segment, D .reverse(), H, tag::divided_in, 10 );

	Mesh ABCD ( tag::rectangle, AB, BC, CD, DA );
	Mesh EFGH ( tag::rectangle, EF, FG, GH, HE );
	Mesh AEFB ( tag::rectangle, AE, EF, BF .reverse(), AB .reverse() );
	Mesh DHEA ( tag::rectangle, DH, HE, AE .reverse(), DA .reverse() );
	Mesh CGHD ( tag::rectangle, CG, GH, DH .reverse(), CD .reverse() );
	Mesh BFGC ( tag::rectangle, BF, FG, CG .reverse(), BC .reverse() );

	Mesh cube ( tag::cube, ABCD, EFGH .reverse(), BFGC, DHEA, CGHD, AEFB );

	cube .export_to_file ( tag::msh, "cube.msh");
	
	std::cout << "produced file cube.msh" << std::endl;

}  // end of main


