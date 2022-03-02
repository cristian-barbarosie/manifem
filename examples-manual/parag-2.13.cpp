
// example presented in paragraph 2.13 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// a 3D star-like shape

#include "maniFEM.h"

using namespace maniFEM;

int main ()

{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz [0], y = xyz [1], z = xyz [2];

	Cell A ( tag::vertex );  x (A) = -1.;  y (A) = -1.;  z (A) = -1.;
	Cell B ( tag::vertex );  x (B) =  1.;  y (B) = -1.;  z (B) = -1.;
	Cell C ( tag::vertex );  x (C) =  1.;  y (C) = -1.;  z (C) =  1.;
	Cell D ( tag::vertex );  x (D) = -1.;  y (D) = -1.;  z (D) =  1.;
	Cell E ( tag::vertex );  x (E) = -1.;  y (E) =  1.;  z (E) = -1.;
	Cell F ( tag::vertex );  x (F) =  1.;  y (F) =  1.;  z (F) = -1.;
	Cell G ( tag::vertex );  x (G) =  1.;  y (G) =  1.;  z (G) =  1.;
	Cell H ( tag::vertex );  x (H) = -1.;  y (H) =  1.;  z (H) =  1.;

	const double p = 3.3;  // recommended values p > 3.
	const double q = (1.-p) / 4.;

	// although this restriction is not explicitly imposed,
	// segments AB and GH will stay within the plane y == z
	// because A, B, G and H do and because the equation given below
	// is symmetric in y and z

	// in paragraph 2.15, we impose explicitly a second equation
	
	RR3 .implicit ( x*x + q*(y+z)*(y+z) == 2.-p );
	Mesh AB ( tag::segment, A .reverse(), B, tag::divided_in, 15 );
	Mesh GH ( tag::segment, G .reverse(), H, tag::divided_in, 15 );

	// likewise, segments CD and EF will stay within the plane y == -z
	RR3 .implicit ( x*x + q*(y-z)*(y-z) == 2.-p );
	Mesh CD ( tag::segment, C .reverse(), D, tag::divided_in, 15 );
	Mesh EF ( tag::segment, E .reverse(), F, tag::divided_in, 15 );

	// segments AE and CG will stay within the plane x == z, and so on
	RR3 .implicit ( y*y + q*(x+z)*(x+z) == 2.-p );
	Mesh AE ( tag::segment, A .reverse(), E, tag::divided_in, 15 );
	Mesh CG ( tag::segment, C .reverse(), G, tag::divided_in, 15 );

	RR3 .implicit ( y*y + q*(x-z)*(x-z) == 2.-p );
	Mesh BF ( tag::segment, B .reverse(), F, tag::divided_in, 15 );
	Mesh DH ( tag::segment, D .reverse(), H, tag::divided_in, 15 );

	RR3 .implicit ( z*z + q*(x+y)*(x+y) == 2.-p );
	Mesh DA ( tag::segment, D .reverse(), A, tag::divided_in, 15 );
	Mesh FG ( tag::segment, F .reverse(), G, tag::divided_in, 15 );

	RR3 .implicit ( z*z + q*(x-y)*(x-y) == 2.-p );
	Mesh BC ( tag::segment, B .reverse(), C, tag::divided_in, 15 );
	Mesh HE ( tag::segment, H .reverse(), E, tag::divided_in, 15 );

	RR3 .implicit ( x*x + y*y - p*z*z == 2.-p );
	Mesh AEFB ( tag::rectangle, AE, EF, BF .reverse(), AB .reverse() );
	Mesh CGHD ( tag::rectangle, CG, GH, DH .reverse(), CD .reverse() );
	
	RR3 .implicit ( x*x - p*y*y + z*z == 2.-p );
	Mesh ABCD ( tag::rectangle, AB, BC, CD, DA );
	Mesh EFGH ( tag::rectangle, EF, FG, GH, HE );

	RR3 .implicit ( - p*x*x + y*y + z*z == 2.-p );
	Mesh DHEA ( tag::rectangle, DH, HE, AE .reverse(), DA .reverse() );
	Mesh BFGC ( tag::rectangle, BF, FG, CG .reverse(), BC .reverse() );

	// back to the initial 3D Euclidian space :
	
	RR3 .set_as_working_manifold();
	Mesh cube ( tag::cube, ABCD, EFGH .reverse(), BFGC, DHEA, CGHD, AEFB );

	cube .export_to_file ( tag::msh, "cube.msh");
	
	std::cout << "produced file cube.msh" << std::endl;

}  // end of main


