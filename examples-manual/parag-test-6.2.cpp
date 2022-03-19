
// example presented in paragraph 6.2 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// integrates functions on a sphere

#include "maniFEM.h"

using namespace maniFEM;

int main ()

{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz [0], y = xyz [1], z = xyz [2];

	RR3 .implicit ( x*x + y*y + z*z == 1. );
	Mesh sphere ( tag::frontal, tag::desired_length, 0.1 );

	sphere .export_to_file ( tag::msh, "sphere.msh");
	
	Integrator integr ( tag::Gauss, tag::tri_6 );

	std::cout << "integral of 1 = " << integr ( 1., tag::on, sphere ) << std::endl;
	std::cout << "integral of x = " << integr ( x, tag::on, sphere ) << std::endl;
	std::cout << "integral of x*x = " << integr ( x*x, tag::on, sphere ) << std::endl;
	std::cout << "integral of x*y = " << integr ( x*y, tag::on, sphere ) << std::endl;

}  // end of main
