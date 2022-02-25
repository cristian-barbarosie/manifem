
// example presented in paragraph 3.9 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// a "bumpy" hemisphere meshed progressively

#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;


int main ()

{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz [0], y = xyz [1], z = xyz [2];

	Manifold nut = RR3 .implicit ( x*x + y*y + z*z + 1.5*x*y*z == 1. );

	// build the base (a closed curve)
	nut .implicit ( x*x + 3.*z == 0. );

	Mesh circle ( tag::frontal, tag::desired_length, 0.1, tag::orientation, tag::random );

	// should produce error message if no orientation is provided
	// should produce error message  with tag::orientation, tag::intrinsic
	// should produce error message  with tag::orientation, tag::inherent

	nut .set_as_working_manifold();
	Mesh bumpy ( tag::frontal, tag::boundary, circle, tag::desired_length, 0.098,
							 tag::orientation, tag::random                                   );

	Mesh::Iterator it = bumpy .iterator ( tag::over_vertices );
	for ( it .reset();  it .in_range(); it++ )
	{	Cell ver = *it;
		if ( ver .is_inner_to ( bumpy ) )  bumpy .baricenter ( ver );  }
	
	bumpy .export_to_file ( tag::msh, "bumpy.msh");
	
	std::cout << "produced file bumpy.msh" << std::endl;

}  // end of main
