
// example presented in paragraph 7.26 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// a rotating vector field on a circle, defined as eigenvector

#include <fstream>

#include "maniFEM.h"
using namespace maniFEM;

#include <Eigen/Eigenvalues>


int main ( )

{	Manifold RR ( tag::Euclid, tag::of_dim, 1 );
	Function theta = RR .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	const double pi = 3.1415926536;
	Manifold::Action g ( tag::transforms, theta, tag::into, theta + 2.*pi );
	Manifold circle_manif = RR .quotient (g);

	// uv  will be a vector field defined at vertices
	Function uv ( tag::lives_on, tag::vertices, tag::has_size, 2 );
	Function u = uv[0], v = uv[1];

	Cell A ( tag::vertex );  theta (A) = 0.;
	Mesh circle ( tag::segment, A .reverse(), A, tag::divided_in, 20, tag::winding, g );

	// it should be noted that, at this stage, theta is not a multifunction
	// theta is the only coordinate on RR, which is a Euclidian manifold
	// the coordinate on circle_manif is a multi-function :
	Function theta_mv = circle_manif .coordinates();

	// first, we compute finite diferences without taking windings into account
	// and we note a concentration

	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = circle .iterator ( tag::over_segments, tag::require_order );
	it .reset();  assert ( it .in_range() );
	Cell first_seg = *it;
	Cell BB = first_seg .tip();
	Eigen::Matrix2d M;
	M ( 0, 0 ) =  std::cos ( theta (BB) );
	M ( 0, 1 ) =  M ( 1, 0 ) = std::sin ( theta (BB) );
	M ( 1, 1 ) = - std::cos ( theta (BB) );
	Eigen::SelfAdjointEigenSolver < Eigen::Matrix2d > es;
	es .compute (M);
	Eigen::Matrix2d eigenvec = es .eigenvectors();
	u (BB) = eigenvec .col (0) (0);
	v (BB) = eigenvec .col (0) (1);
	for ( it++; it .in_range(); it++ )
	{	Cell seg = *it;
		Cell AA = seg .base() .reverse();
		BB = seg .tip();
		Cell V = *it;
		M ( 0, 0 ) = std::cos ( theta (BB) );
		M ( 0, 1 ) =  M ( 1, 0 ) = std::sin ( theta (BB) );
		M ( 1, 1 ) = - std::cos ( theta (BB) );
		es .compute (M);
		eigenvec = es .eigenvectors();
		// among four candidates, we choose the one which is closest to previous eigenvec
		short int index = -1, sign = 0;
		double dist_min = 100.;
		Eigen::Vector2d previous_eigenvec ( u (AA), v (AA) );
		for ( short int i = 0; i < 2; i++ )
		for ( short int s = -1; s < 2; s += 2 )
		{	double d = ( eigenvec .col(i) - s * previous_eigenvec ) .norm();
			if ( d < dist_min )
			{	dist_min = d;  index = i;  sign = s;  }                        }
		assert ( index >= 0 );  assert ( sign != 0 );
		u (BB) = sign * eigenvec .col (index) (0);
		v (BB) = sign * eigenvec .col (index) (1);
		double du = u (BB) - u (AA),
		       dv = v (BB) - v (AA),
		       dt = theta_mv ( BB, tag::winding, seg .winding() ) - theta_mv (AA);
		std::cout << ( du*du + dv*dv ) / dt / dt << std::endl;                    }
	// the above works fine, but on the last segment u and v have a jump :
	Cell AA = first_seg .base() .reverse();
	BB = first_seg .tip();
	double du = u (BB) - u (AA),
	       dv = v (BB) - v (AA),
	       dt = theta_mv ( BB, tag::winding, first_seg .winding() ) - theta_mv (AA);
	std::cout << ( du*du + dv*dv ) / dt / dt << std::endl;	
	std::cout << std::endl;
	} // just a block of code

	// by using multi-functions, everything goes into place :

	Function u_mv = u .make_multivalued ( tag::through, g, tag::becomes, -u );
	Function v_mv = v .make_multivalued ( tag::through, g, tag::becomes, -v );

	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = circle .iterator ( tag::over_segments, tag::require_order );
	it .reset();  assert ( it .in_range() );
	Cell first_seg = *it;
	Cell BB = first_seg .tip();
	Eigen::Matrix2d M;
	M ( 0, 0 ) =  std::cos ( theta (BB) );
	M ( 0, 1 ) =  M ( 1, 0 ) = std::sin ( theta (BB) );
	M ( 1, 1 ) = - std::cos ( theta (BB) );
	Eigen::SelfAdjointEigenSolver < Eigen::Matrix2d > es;
	es .compute (M);
	Eigen::Matrix2d eigenvec = es .eigenvectors();
	u (BB) = eigenvec .col (0) (0);
	v (BB) = eigenvec .col (0) (1);
	for ( it++; it .in_range(); it++ )
	{	Cell seg = *it;
		Cell AA = seg .base() .reverse();
		BB = seg .tip();
		Cell V = *it;
		M ( 0, 0 ) = std::cos ( theta (BB) );
		M ( 0, 1 ) =  M ( 1, 0 ) = std::sin ( theta (BB) );
		M ( 1, 1 ) = - std::cos ( theta (BB) );
		es .compute (M);
		eigenvec = es .eigenvectors();
		// among four candidates, we choose the one which is closest to previous eigenvec
		short int index = -1, sign = 0;
		double dist_min = 100.;
		// we want to match the current vector with the previous vector
		// since we are looking back, we change the sign of the winding number
		Eigen::Vector2d previous_eigenvec ( u_mv ( AA, tag::winding, - seg .winding() ),
		                                    v_mv ( AA, tag::winding, - seg .winding() ) );
		for ( short int i = 0; i < 2; i++ )
		for ( short int s = -1; s < 2; s += 2 )
		{	double d = ( eigenvec .col(i) - s * previous_eigenvec ) .norm();
			if ( d < dist_min )
			{	dist_min = d;  index = i;  sign = s;  }                        }
		assert ( index >= 0 );  assert ( sign != 0 );
		u (BB) = sign * eigenvec .col (index) (0);
		v (BB) = sign * eigenvec .col (index) (1);
		double du = u_mv ( BB, tag::winding, seg .winding() ) - u_mv (AA),
		       dv = v_mv ( BB, tag::winding, seg .winding() ) - v_mv (AA),
		       dt = theta_mv ( BB, tag::winding, seg .winding() ) - theta_mv (AA);
		std::cout << ( du*du + dv*dv ) / dt / dt << std::endl;                    }
	// let us check now the first segment
	Cell AA = first_seg .base() .reverse();
	BB = first_seg .tip();
	double du = u_mv ( BB, tag::winding, first_seg .winding() ) - u_mv (AA),
	       dv = v_mv ( BB, tag::winding, first_seg .winding() ) - v_mv (AA),
	       dt = theta_mv ( BB, tag::winding, first_seg .winding() ) - theta_mv (AA);
	std::cout << ( du*du + dv*dv ) / dt / dt << std::endl;
	std::cout << std::endl;
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
		solution_file << numbering [P] + 1 << " " << u_mv (P) << " " << v_mv (P) << " 0. "<< std::endl;  }
	} // just a block of code

	std::cout << "produced file circle.msh" << std::endl;

	return 0;
	
}  // end of main

