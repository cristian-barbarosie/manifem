
// solve a celullar problem


#include "maniFEM.h"
#include <Eigen/Sparse>
#include <Eigen/OrderingMethods>
using namespace maniFEM;

int main ( )

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0], y = xy[1];

	Cell A ( tag::vertex );  x(A) = 0.;  y(A) = 0.;
	Cell B ( tag::vertex );  x(B) = 1.;  y(B) = 0.;
	Cell C ( tag::vertex );  x(C) = 1.;  y(C) = 1.;
	Cell D ( tag::vertex );  x(D) = 0.;  y(D) = 1.;

	Mesh AB ( tag::segment, A.reverse(), B, tag::divided_in, 10 );
	Mesh BC ( tag::segment, B.reverse(), C, tag::divided_in, 10 );
	Mesh CD ( tag::segment, C.reverse(), D, tag::divided_in, 10 );
	Mesh DA ( tag::segment, D.reverse(), A, tag::divided_in, 10 );

	Mesh square ( tag::rectangle, AB, BC, CD, DA );

	// declare the type of finite element
	FiniteElement fe ( tag::with_master, tag::quadrangle, tag::Lagrange, tag::of_degree, 1 );
	Integrator integ = fe.set_integrator ( tag::Gauss, tag::quad_9 );

	std::map < Cell, size_t > numbering;
	{ // just a block of code for hiding 'it' and 'counter'
	CellIterator it = square.iterator ( tag::over_vertices );
	size_t counter = 0;
	for ( it.reset() ; it.in_range(); it++ )
	{	Cell V = *it;  numbering[V] = counter;  ++counter;  }
	assert ( counter == numbering.size() );
	} // just a block of code

	size_t size_matrix = numbering.size();
	std::cout << "global matrix " << size_matrix + 1 << "x" << size_matrix << std::endl;
	Eigen::SparseMatrix < double > matrix_A ( size_matrix + 1, size_matrix );
	
	matrix_A.reserve ( Eigen::VectorXi::Constant ( size_matrix, 10 ) );

	Eigen::VectorXd vector_b ( size_matrix + 1 ), vector_sol ( size_matrix );
	vector_b.setZero();

	// run over all square cells composing 'square'
	{ // just a block of code for hiding 'it'
	CellIterator it = square.iterator ( tag::over_cells_of_max_dim );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell small_tri = *it;
		std::cout << "******************" << std::endl;
		fe.dock_on ( small_tri );
		// run twice over the four vertices of 'small_tri'
		CellIterator it1 = small_tri.boundary().iterator ( tag::over_vertices );
		CellIterator it2 = small_tri.boundary().iterator ( tag::over_vertices );
		for ( it1.reset(); it1.in_range(); it1++ )
		for ( it2.reset(); it2.in_range(); it2++ )
		{	Cell V = *it1, W = *it2;  // V may be the same as W, no problem about that
			// std::cout << "vertices V=(" << x(V) << "," << y(V) << ") " << numbering[V] << ", W=("
			// 					<< x(W) << "," << y(W) << ") " << numbering[W]) << std::endl;
			Function psi_V = fe .basis_function ( V ),
			         psi_W = fe .basis_function ( W ),
			         d_psi_V_dx = psi_V .deriv ( x ),
			         d_psi_V_dy = psi_V .deriv ( y ),
			         d_psi_W_dx = psi_W .deriv ( x ),
			         d_psi_W_dy = psi_W .deriv ( y );
			// 'fe' is already docked on 'small_tri' so this will be the domain of integration
			matrix_A.coeffRef ( numbering[V], numbering[W] ) +=
				fe.integrate ( d_psi_V_dx * d_psi_W_dx + d_psi_V_dy * d_psi_W_dy );
			std::cout << numbering[V] << " " << numbering[W] << " "
								<< matrix_A.coeffRef ( numbering[V], numbering[W] ) << std::endl;
		}  }
	} // just a block of code 
	
}
