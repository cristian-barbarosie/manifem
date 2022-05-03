
// example presented in paragraph 7.24 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// a rotating vector field on a circle

#include <fstream>

#include "maniFEM.h"
using namespace maniFEM;


int main ( )

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0], y = xy[1];

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

	std::cout << "produced file circle.msh" << std::endl;
	
	return 0;
	
}  // end of main
