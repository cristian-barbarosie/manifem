
// example presented in paragraph 7.3 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// mesh a cylinder

#include "maniFEM.h"
using namespace maniFEM;


int main ( )

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function theta_z = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function theta = theta_z [0], z = theta_z [1];
	const double pi = 3.1415926536;
	Manifold::Action g ( tag::transforms, theta_z, tag::into, (theta+2*pi) && z );
	Manifold cylinder_manif = RR2 .quotient ( g );
	
	Cell A ( tag::vertex );  theta (A) = 0.;  z (A) = -1.;
	Cell B ( tag::vertex );  theta (B) = 0.,  z (B) =  1.;

	Mesh AA ( tag::segment, A .reverse(), A, tag::divided_in, 30, tag::winding, g );
	Mesh AB ( tag::segment, A .reverse(), B, tag::divided_in, 15 );  // no winding
	Mesh BB ( tag::segment, B .reverse(), B, tag::divided_in, 30, tag::winding, -g );

	Mesh cylinder ( tag::rectangle, AA, AB, BB, AB .reverse(), tag::winding );
	// the tag::winding provides no specific information,
	// it just warns maniFEM that we are on a quotient manifold
	// and that it must take winding segments into account
	// specific information about winding numbers is included in the three segments AA, AB, BB

	// define new coordinates x and y as arithmetic expressions of theta
	Function x = cos ( theta ), y = sin ( theta );
	// forget about theta; in future statements, x, y and z will be used
	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );	
	RR3 .set_coordinates ( x && y && z );
	cylinder .export_to_file ( tag::msh, "cylinder.msh");
	
	std::cout << "produced file cylinder.msh" << std::endl;

}  // end of main

