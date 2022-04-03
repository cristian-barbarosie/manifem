

#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;

int main ()

{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz [0], y = xyz [1], z = xyz[2];

	size_t n = 10;

	Cell A ( tag::vertex, tag::of_coords, { -0.5, 0., 0. } );
	Cell B ( tag::vertex, tag::of_coords, {  0.5, 0., 0. } );

	Mesh AB ( tag::segment, A .reverse(), B, tag::divided_in, n );

	double height = 0.8;
	Manifold circle = RR3 .implicit ( z == height, x*x + y*y == 1. );
	Cell C1 ( tag::vertex, tag::of_coords, {  1. ,  0.  , height }, tag::project );
	Cell C2 ( tag::vertex, tag::of_coords, {  0.5,  0.86, height }, tag::project );
	Cell C3 ( tag::vertex, tag::of_coords, { -0.5,  0.86, height }, tag::project );	
	Cell C4 ( tag::vertex, tag::of_coords, { -1. ,  0.  , height }, tag::project );
	Cell C5 ( tag::vertex, tag::of_coords, { -0.5, -0.86, height }, tag::project );
	Cell C6 ( tag::vertex, tag::of_coords, {  0.5, -0.86, height }, tag::project );
	Mesh C1C2 ( tag::segment, C1 .reverse(), C2, tag::divided_in, n );
	Mesh C2C3 ( tag::segment, C2 .reverse(), C3, tag::divided_in, n );
	Mesh C3C4 ( tag::segment, C3 .reverse(), C4, tag::divided_in, n );
	Mesh C4C5 ( tag::segment, C4 .reverse(), C5, tag::divided_in, n );
	Mesh C5C6 ( tag::segment, C5 .reverse(), C6, tag::divided_in, n );
	Mesh C6C1 ( tag::segment, C6 .reverse(), C1, tag::divided_in, n );

	RR3 .set_as_working_manifold();
	Mesh BC1 ( tag::segment, B .reverse(), C1, tag::divided_in, n );
	Mesh BC2 ( tag::segment, B .reverse(), C2, tag::divided_in, n );
	Mesh AC3 ( tag::segment, A .reverse(), C3, tag::divided_in, n );
	Mesh AC4 ( tag::segment, A .reverse(), C4, tag::divided_in, n );
	Mesh AC5 ( tag::segment, A .reverse(), C5, tag::divided_in, n );
	Mesh BC6 ( tag::segment, B .reverse(), C6, tag::divided_in, n );

	Mesh ABC2C3 ( tag::rectangle, AB, BC2, C2C3, AC3 .reverse(), tag::with_triangles );
	Mesh BAC5C6 ( tag::rectangle, AB .reverse(), AC5, C5C6, BC6 .reverse(), tag::with_triangles );
	Mesh AC3C4 ( tag::triangle, AC3, C3C4, AC4 .reverse() );
	Mesh AC4C5 ( tag::triangle, AC4, C4C5, AC5 .reverse() );
	Mesh BC6C1 ( tag::triangle, BC6, C6C1, BC1 .reverse() );
	Mesh BC1C2 ( tag::triangle, BC1, C1C2, BC2 .reverse() );

	Mesh filter ( tag::join, { ABC2C3, BAC5C6, AC3C4, AC4C5, BC6C1, BC1C2 } );

	Mesh::Iterator it = filter .iterator ( tag::over_vertices );
	for ( int i = 1; i < 20; i++ )
	for ( it .reset(); it .in_range(); it++ )
	{	Cell ver = *it;
		if ( ver .belongs_to ( AB ) ) continue;
		if ( ver .is_inner_to ( filter ) ) filter .baricenter ( ver );  }

	filter .export_to_file ( tag::msh, "filter.msh");

}
	
