

#include "maniFEM.h"

using namespace maniFEM;

int main () {

	// we choose our (geometric) space dimension :
	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	
	// xyz is a map defined on our future mesh with values in RR3 :
	Function xyz = RR3.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );

	// we can extract components of xyz using the [] operator :
	Function x = xyz[0],  y = xyz[1],  z = xyz[2];

	// Let's build a rectangular mesh. First, the four corners :
	Cell SW ( tag::vertex );  x(SW) =  0.;  y(SW) =  0.;  z(SW) =  0.;
	Cell SE ( tag::vertex );  x(SE) =  1.;  y(SE) =  0.;  z(SE) =  0.;
	Cell NE ( tag::vertex );  x(NE) =  1.;  y(NE) =  1.;  z(NE) =  0.;
	Cell NW ( tag::vertex );  x(NW) =  0.;  y(NW) =  1.;  z(NW) =  1.;

	Cell s ( tag::segment, SW.reverse(), SE );
	Mesh south ( tag::segment, SW.reverse(), SE, tag::divided_in, 10 );
	
	std::cout << "main, before exit" << std::endl;  exit(0);

	// now build the four sides of the rectangle :
	Mesh east  ( tag::segment, SE.reverse(), NE, tag::divided_in, 10 );
	Mesh north ( tag::segment, NE.reverse(), NW, tag::divided_in, 10 );
	Mesh west  ( tag::segment, NW.reverse(), SW, tag::divided_in, 10 );

	// and now the rectangle :
	Mesh rect_mesh ( tag::rectangle, south, east, north, west );

	CellIterator it1 = rect_mesh.iterator ( tag::over_vertices );
	for ( it1.reset(); it1.in_range(); it1++ )
	{	Cell P = *it1;
		CellIterator it2 = rect_mesh.iterator ( tag::over_segments, tag::around, P );
		CellIterator it3 = rect_mesh.iterator ( tag::over_vertices, tag::around, P );
		for ( it2.reset(), it3.reset(); it2.in_range(); it2++, it3++ )
		{	assert ( it3.in_range() );
			assert ( (*it2).tip() == P );
			assert ( (*it2).base().reverse() == *it3 );  }
		assert ( not it3.in_range() );                                                }

	std::cout << "end of main" << std::endl << std::flush;
}
