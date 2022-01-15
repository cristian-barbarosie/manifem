
#include "maniFEM.h"

using namespace maniFEM;


//-----------------------------------------------------------------------------------//


int main ()
	
{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0],  y = xy[1];

	Cell X ( tag::vertex );  x(X) = 0.;  y(X) = 0.;
	Cell Y ( tag::vertex );  x(Y) = 1.;  y(Y) = 0.;
	Cell Z ( tag::vertex );  x(Z) = 1.;  y(Z) = 1.;
	Cell T ( tag::vertex );  x(T) = 0.;  y(T) = 1.;

	Mesh XY ( tag::segment, X.reverse(), Y, tag::divided_in, 2 );
	Mesh YZ ( tag::segment, Y.reverse(), Z, tag::divided_in, 2 );
	Mesh ZT ( tag::segment, Z.reverse(), T, tag::divided_in, 2 );
	Mesh TX ( tag::segment, T.reverse(), X, tag::divided_in, 2 );

	Mesh rect_mesh ( tag::rectangle, XY, YZ, ZT, TX );

	// reverse each square cell
	std::vector < Cell > vec_of_cells;
	vec_of_cells .reserve (4);
	{  // just a block of code for hiding 'it'
	Mesh::Iterator it = rect_mesh .iterator ( tag::over_cells_of_max_dim );
	for ( it .reset(); it .in_range(); it++ )
		vec_of_cells .push_back ( *it );
	} {  // just a block of code for hiding 'it'
	std::vector < Cell > ::iterator it;
	for ( it = vec_of_cells .begin(); it != vec_of_cells .end(); it++ )
	{	Cell sq = *it;  sq .remove_from_mesh ( rect_mesh );  }
	for ( it = vec_of_cells .begin(); it != vec_of_cells .end(); it++ )
	{	Cell sq = *it;  sq .reverse() .add_to_mesh ( rect_mesh );  }
	} {  // just a block of code for hiding 'it'

	// try to use an iterator around a vertex
	Mesh::Iterator it = rect_mesh .iterator ( tag::over_vertices );
	it .reset();  assert ( it .in_range() );
	for ( size_t count = 0; count < 2; count++ )
		{ it++;  assert ( it .in_range() );  }
	Cell ver = * it;
	std::cout << x(ver) << " " << y(ver) << std::endl;
	Mesh::Iterator itt = rect_mesh .iterator ( tag::over_vertices, tag::around, ver );
	itt .reset();  assert ( itt .in_range() );
	}  // just a block of code for hiding 'it'
	
}  // end of main	


		
