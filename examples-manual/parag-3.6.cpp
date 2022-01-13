
// example presented in paragraph 3.6 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// meshes a sphere progressively

#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;

int main ()

{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz [0], y = xyz [1], z = xyz [2];

	RR3 .implicit ( x*x + y*y + z*z == 1. );
	Mesh sphere ( tag::progressive, tag::desired_length, 0.11 );

	// should work just the same with     tag::orientation, tag::inherent
	// should work with                   tag::orientation, tag::random
	// should produce error message with  tag::orientation, tag::intrinsic

	// uncomment the three lines below to slightly improve the quality of the mesh
	
	// Mesh::Iterator it = sphere .iterator ( tag::over_vertices );
	// for ( it .reset(); it .in_range(); it++ )
	// {	Cell ver = *it;  sphere .baricenter ( ver );  }
	
	sphere .export_msh ("sphere.msh");

	std::cout << "produced file sphere.msh" << std::endl;
}
