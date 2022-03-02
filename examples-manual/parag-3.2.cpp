
// example presented in paragraph 3.2 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// builds a circle then fills the disk, all frontal

#include "maniFEM.h"
using namespace maniFEM;


int main ( )

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy [0], y = xy [1];

	Manifold circle_manif = RR2 .implicit ( x*x + y*y == 1. );

	// Mesh circle ( tag::frontal, tag::entire_manifold, circle_manif, tag::desired_length, 0.2 );
	// Mesh circle ( tag::frontal, tag::desired_length, 0.2 );

	Cell A ( tag::vertex );  x (A) = 0.;  y (A) = 1.;
	Mesh circle ( tag::frontal, tag::start_at, A, tag::desired_length, 0.2 );

	// should work just the same with     tag::orientation, tag::inherent
	// should work with                   tag::orientation, tag::random,
	//        but constructor  Mesh disk  below may enter in an endless process
	//        of meshing the exterior of the disk
	// should produce error message with  tag::orientation, tag::intrinsic

	RR2 .set_as_working_manifold();
	Mesh disk ( tag::frontal, tag::boundary, circle, tag::desired_length, 0.2 );

	// should work just the same with     tag::orientation, tag::intrinsic
	// should produce error message with  tag::orientation, tag::random
	// should produce error message with  tag::orientation, tag::inherent

	disk .export_to_file ( tag::msh, "disk.msh");
	std::cout << "produced file disk.msh" << std::endl;

}  // end of main
