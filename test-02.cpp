

#include "maniFEM.h"

using namespace maniFEM;


	
int main ( )

{	// begin with the usual two-dimensional space
	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0], y = xy[1];

	Cell A ( tag::vertex );  x(A) =  0.  ;  y(A) = -0.35;
	Cell B ( tag::vertex );  x(B) =  0.  ;  y(B) = -1.05;
	Cell C ( tag::vertex );  x(C) =  1.15;  y(C) = -2.15;

	A.core->name = "A";
	B.core->name = "B";
	C.core->name = "C";

	Cell BA ( tag::segment, B.reverse(), A );
	Cell BC ( tag::segment, B.reverse(), C );
	Cell CA ( tag::segment, C.reverse(), A );

	BA.core->name = "BA";
	BC.core->name = "BC";
	CA.core->name = "CA";
	
	Cell ABC ( tag::triangle, BA.reverse(), BC, CA );
	ABC.core->name = "ABC";

	Mesh msh ( tag::whose_core_is,
     new Mesh::Fuzzy ( tag::of_dimension, 3, tag::minus_one,
											 tag::one_dummy_wrapper ),
						 tag::freshly_created, tag::is_positive           );
	msh.core->name = "msh";

	ABC .add_to_mesh ( msh );
	
  BA .reverse() .cut_from_bdry_of ( ABC );

	std::cout << "reached end of main" << std::endl;
}

