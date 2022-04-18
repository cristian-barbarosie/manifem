
// example presented in paragraph 7.27 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// a rotating vector field on a circle, defined by  coef * theta

#include <fstream>

#include "maniFEM.h"
using namespace maniFEM;


int main ( )

{	Manifold RR ( tag::Euclid, tag::of_dim, 1 );
	Function theta = RR .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	const double pi = 3.1415926536;
	Manifold::Action g ( tag::transforms, theta, tag::into, theta + 2.*pi );
	Manifold circle_manif = RR .quotient (g);

	Cell A ( tag::vertex );  theta (A) = 0.;
	Mesh circle ( tag::segment, A .reverse(), A, tag::divided_in, 20, tag::winding, g );

	// it should be noted that, at this stage, theta is not a multifunction
	// theta is the only coordinate on RR, which is a Euclidian manifold
	// the coordinate on circle_manif is a multi-function :
	Function theta_mv = circle_manif .coordinates();

	const double coef = 0.3333;
	Function u = cos ( coef * theta );
	Function v = sin ( coef * theta );
	Function uv = u && v;

	const double c = std::cos ( coef * 2.* pi ), s = std::sin ( coef * 2.* pi );
	Function uv_mv = uv .make_multivalued
		( tag::through, g, tag::becomes, ( c*u - s*v ) && ( s*u + c*v ) );

	// it would be comfortable to extract components of uv_mv :
	// Function u_mv = uv_mv [0], v_mv = uv_mv [1]
	// but maniFEM does not allow extracting components of multi-valued vector fields
	// so we just use uv_mv and extract components after evaluation at a vertex

	{ // just a block of code for hiding 'it'
	// first, we compute finite diferences without taking windings into account
	// and we note a concentration
	Mesh::Iterator it = circle .iterator ( tag::over_segments, tag::require_order );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell seg = *it;
		Cell AA = seg .base() .reverse(), BB = seg .tip();
		double du = u (BB) - u (AA),
		       dv = v (BB) - v (AA),
		       dt = theta_mv ( BB, tag::winding, seg .winding() ) - theta_mv (AA);
		std::cout << ( du*du + dv*dv ) / dt / dt << std::endl;                    }
	std::cout << std::endl;
	// now we take windings into account and notice that all derivatives have the same magnitude
	for ( it .reset(); it .in_range(); it++ )
	{	Cell seg = *it;
		Cell AA = seg .base() .reverse(), BB = seg .tip();
		std::vector < double > tmp_BB = uv_mv ( BB, tag::winding, seg .winding() ),
		                       tmp_AA = uv_mv (AA);
		double du = tmp_BB [0] - tmp_AA [0],
		       dv = tmp_BB [1] - tmp_AA [1],
		       dt = theta_mv ( BB, tag::winding, seg .winding() ) - theta_mv (AA);
		std::cout << ( du*du + dv*dv ) / dt / dt << std::endl;                     }
	} // just a block of code
	
	// numbering is explained in paragraph 6.3 of the manual
	std::map < Cell, size_t > numbering;
	{ // just a block of code for hiding 'it' and 'counter'
	Mesh::Iterator it = circle .iterator ( tag::over_vertices );
	size_t counter = 0;
	for ( it .reset() ; it .in_range(); it ++ )
	{	Cell V = *it;  numbering [V] = counter;  ++ counter;  }
	assert ( counter == numbering .size() );
	} // just a block of code

	// define new coordinates x and y as arithmetic expressions of theta
	Function x = cos ( theta ), y = sin ( theta );

	// forget about theta; in future statements, x and y will be used
	Manifold RR2 ( tag::Euclid, tag::of_dimension, 2 );
	RR2 .set_coordinates ( x && y );
	circle .export_to_file ( tag::msh, "circle.msh", numbering );

	{ // just a block of code for hiding variables
	std::ofstream solution_file ("circle.msh", std::fstream::app );
	solution_file << "$NodeData" << std::endl;
	solution_file << "1" << std::endl;   // one string follows
	solution_file << "\"eigenvector\"" << std::endl;
	solution_file << "1" << std::endl;   //  one real follows
	solution_file << "0.0" << std::endl;  // time [??]
	solution_file << "3" << std::endl;   // three integers follow
	solution_file << "0" << std::endl;   // time step [??]
	solution_file << "3" << std::endl;  // scalar values of u
	solution_file << circle .number_of ( tag::vertices ) << std::endl;
	// number of values listed below
	Mesh::Iterator it = circle .iterator ( tag::over_vertices );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell P = *it;
		std::vector < double > tmp = uv_mv (P);
		solution_file << numbering [P] + 1 << " " << tmp [0] << " " << tmp [1] << " 0. "<< std::endl;  }
	} // just a block of code

	std::cout << "produced file circle.msh" << std::endl;

	return 0;
	
}  // end of main

