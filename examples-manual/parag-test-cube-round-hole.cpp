
#include "maniFEM.h"

using namespace maniFEM;

int main ()

{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz [0], y = xyz [1], z = xyz [2];

	const double r2 = 0.2;
	Manifold sphere = RR3 .implicit ( x*x + y*y + z*z == r2 );

	Cell A ( tag::vertex );  x (A) = -1.;  y (A) =  1.;  z (A) = -1.;
	Cell B ( tag::vertex );  x (B) =  1.;  y (B) =  1.;  z (B) = -1.;
	Cell C ( tag::vertex );  x (C) =  1.;  y (C) =  1.;  z (C) =  1.;
	Cell D ( tag::vertex );  x (D) = -1.;  y (D) =  1.;  z (D) =  1.;
	Cell E ( tag::vertex );  x (E) = -1.;  y (E) =  1.;  z (E) = -1.;
	Cell F ( tag::vertex );  x (F) =  1.;  y (F) =  1.;  z (F) = -1.;
	Cell G ( tag::vertex );  x (G) =  1.;  y (G) =  1.;  z (G) =  1.;
	Cell H ( tag::vertex );  x (H) = -1.;  y (H) =  1.;  z (H) =  1.;
	sphere .project (E);  sphere .project (F);  sphere .project (G);  sphere .project (H);  

	Mesh EF ( tag::segment, E .reverse(), F, tag::divided_in, 15 );
	Mesh GH ( tag::segment, G .reverse(), H, tag::divided_in, 15 );
	Mesh FG ( tag::segment, F .reverse(), G, tag::divided_in, 15 );
	Mesh HE ( tag::segment, H .reverse(), E, tag::divided_in, 15 );

	Mesh EFGH ( tag::rectangle, EF, FG, GH, HE );

	RR3 .set_as_working_manifold();

	Mesh AB ( tag::segment, A .reverse(), B, tag::divided_in, 15 );
	Mesh CD ( tag::segment, C .reverse(), D, tag::divided_in, 15 );

	Mesh AE ( tag::segment, A .reverse(), E, tag::divided_in, 8 );
	Mesh BF ( tag::segment, B .reverse(), F, tag::divided_in, 8 );
	Mesh CG ( tag::segment, C .reverse(), G, tag::divided_in, 8 );
	Mesh DH ( tag::segment, D .reverse(), H, tag::divided_in, 8 );

	Mesh BC ( tag::segment, B .reverse(), C, tag::divided_in, 15 );
	Mesh DA ( tag::segment, D .reverse(), A, tag::divided_in, 15 );

	Mesh AEFB ( tag::rectangle, AE, EF, BF .reverse(), AB .reverse() );
	Mesh CGHD ( tag::rectangle, CG, GH, DH .reverse(), CD .reverse() );
	Mesh ABCD ( tag::rectangle, AB, BC, CD, DA );
	Mesh DHEA ( tag::rectangle, DH, HE, AE .reverse(), DA .reverse() );
	Mesh BFGC ( tag::rectangle, BF, FG, CG .reverse(), BC .reverse() );

	RR3 .set_as_working_manifold();
	Mesh cube ( tag::cube, ABCD, EFGH .reverse(), BFGC, DHEA, CGHD, AEFB );

	cube .export_to_file ( tag::msh, "cube.msh");
	
	std::cout << "produced file cube.msh" << std::endl;

}  // end of main

