
#include <fstream>

#include "maniFEM.h"
using namespace maniFEM;


int main1 ( )

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy [0], y = xy [1];

	Manifold circle_manif = RR2 .implicit ( x*x + y*y == 1. );

	// Mesh circle ( tag::frontal, tag::entire_manifold, circle_manif, tag::desired_length, 0.2 );
	// Mesh circle ( tag::frontal, tag::desired_length, 0.2 );

	Mesh circle ( tag::frontal, tag::desired_length, 0.2 );

	// numbering is explained in paragraph 6.3 of the manual
	std::map < Cell, size_t > numbering;
	{ // just a block of code for hiding 'it' and 'counter'
	Mesh::Iterator it = circle .iterator ( tag::over_vertices );
	size_t counter = 0;
	for ( it .reset() ; it .in_range(); it ++ )
	{	Cell V = *it;  numbering [V] = counter;  ++ counter;  }
	assert ( counter == numbering .size() );
	} // just a block of code

	Function r = sqrt ( x*x + y*y );  // equal to 1. actually

	// this is one way of defining the vector field (u,v) :
	// Function u = sqrt ( (r+x) / 2. / r );
	// Function v = y / sqrt ( 2. * r * (r+x) );
	// the above has the disadvantage of a possible division by zero, when x == -r

	// below we define the same vector field using the sign of y
	// Function u = sqrt ( (r+x) / 2. / r );
	// Function v = sign (y) * sqrt ( (r-x) / 2. / r );
	// the above still has a disandvantage : the functions are not differentiable
	// because of the sign on one hand
	// most importantly, because of the square root in the neighbourhood of zero

	// below we define the same vector field by different formulas in two regions :
	// for negative x we only use r-x which is far from 0.
	// for positive x we only use r+x which is far from 0.
	// so functions stay smooth except for the cuts
	Function u ( tag::piecewise,
	             y / sqrt ( 2. * r * (r-x) ), tag::iff, x, tag::less_than, 0.,
	             sqrt ( (r+x) / 2. / r ),     tag::otherwise                  );
	Function v ( tag::piecewise,
	             sqrt ( (r-x) / 2. / r ),     tag::iff, x, tag::less_than, 0.,
	             y / sqrt ( 2. * r * (r+x) ), tag::otherwise                  );
	
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
		solution_file << numbering [P] + 1 << " " << u(P) << " " << v(P) << " 0. "<< std::endl;  }
	} // just a block of code
	
	return 0;
	
}  // end of main1


int main2 ( )

{	Manifold RR ( tag::Euclid, tag::of_dim, 1 );
	Function theta = RR .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	const double pi = 3.1415926536;
	Manifold::Action g ( tag::transforms, theta, tag::into, theta + 2*pi );
	Manifold circle_manif = RR .quotient (g);

	Cell A ( tag::vertex );  theta (A) = 0.01;
	// if we set  theta (A) = 0.  a division by zero will occur in  v
	Mesh circle ( tag::segment, A .reverse(), A, tag::divided_in, 20, tag::winding, g );

	// define new coordinates x and y as arithmetic expressions of theta
	Function x = cos ( theta ), y = sin ( theta );

	// numbering is explained in paragraph 6.3 of the manual
	std::map < Cell, size_t > numbering;
	{ // just a block of code for hiding 'it' and 'counter'
	Mesh::Iterator it = circle .iterator ( tag::over_vertices );
	size_t counter = 0;
	for ( it .reset() ; it .in_range(); it ++ )
	{	Cell V = *it;  numbering [V] = counter;  ++ counter;  }
	assert ( counter == numbering .size() );
	} // just a block of code

	Function u = cos ( theta / 2. );
	Function v = sin ( theta / 2. );

	Function u_mv = u .make_multivalued ( tag::through, g, tag::becomes, -u );
	Function v_mv = v .make_multivalued ( tag::through, g, tag::becomes, -v );

	Function theta_mv = circle_manif .coordinates();

	{ // just a block of code for hiding variables
	Mesh::Iterator it = circle .iterator ( tag::over_segments, tag::require_order );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell seg = *it;
		Cell AA = seg .base() .reverse(), BB = seg .tip();
		double du = u (BB) - u (AA),
		       dv = v (BB) - v (AA),
		       dt = theta_mv ( BB, tag::winding, seg .winding() ) - theta_mv (AA);
		std::cout << ( du*du + dv*dv ) / dt / dt << std::endl;                    }
	std::cout << std::endl;
	for ( it .reset(); it .in_range(); it++ )
	{	Cell seg = *it;
		Cell AA = seg .base() .reverse(), BB = seg .tip();
		double du = u_mv ( BB, tag::winding, seg .winding() ) - u_mv (AA),
		       dv = v_mv ( BB, tag::winding, seg .winding() ) - v_mv (AA),
		       dt = theta_mv ( BB, tag::winding, seg .winding() ) - theta_mv (AA);
		std::cout << ( du*du + dv*dv ) / dt / dt << std::endl;                    }
	} // just a block of code for hiding variables
	
	// forget about theta; in future statements, x and y will be used
	circle_manif .set_coordinates ( x && y );
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
		solution_file << numbering [P] + 1 << " " << u_mv (P) << " " << v_mv (P) << " 0. "<< std::endl;  }
	} // just a block of code

	return 0;
	
}  // end of main2


int main () { main2();  }
