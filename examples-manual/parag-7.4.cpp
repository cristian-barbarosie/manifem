
// example presented in paragraph 7.4 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// build a circle

#include "maniFEM.h"
using namespace maniFEM;


int main ( )

{	Manifold RR ( tag::Euclid, tag::of_dim, 1 );
	Function theta = RR.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	const double pi = 3.1415926536;
	Function::Action g ( tag::transforms, theta, tag::into, theta + 2*pi );
	Manifold circle = RR.quotient ( g );
	
	Cell A ( tag::vertex );  theta ( A ) = 0.;
	Mesh seg ( tag::segment, A.reverse(), A, tag::divided_in, 30, tag::spin, g );

	// define new coordinates x and y as arithmetic expressions of theta
	Function x = cos(theta), y = sin(theta);
	// forget about theta; in future statements, x and y will be used
	circle.set_coordinates ( x && y );
	seg.export_msh ("circle.msh");
	
	std::cout << "produced file circle.msh" << std::endl;
}