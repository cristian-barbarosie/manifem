
// example presented in paragraph 3.2 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// builds a circle then fills the disk, all progressive

#include "maniFEM.h"
using namespace maniFEM;


int main ( )

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0],  y = xy[1];

	Manifold circle_manif = RR2.implicit ( x*x + y*y == 1. );
	Mesh circle ( tag::progressive, tag::entire_manifold, circle_manif, tag::desired_length, 0.2 );
	// Mesh circle ( tag::progressive, tag::desired_length, 0.2 );

	RR2.set_as_working_manifold();
	Mesh disk ( tag::progressive, tag::boundary, circle, tag::desired_length, 0.2 );

	disk.export_msh ("disk.msh");
	std::cout << "produced file disk.msh" << std::endl;
}
