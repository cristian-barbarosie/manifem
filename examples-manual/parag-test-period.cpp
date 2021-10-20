
// solve a celullar problem
// square periodicity, triangular elements, circular hole


#include "maniFEM.h"
#include <Eigen/Sparse>
#include <Eigen/OrderingMethods>
using namespace maniFEM;

int main ( )

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0], y = xy[1];

	size_t n = 20;
	double d = 2.6 / double(n);

	Cell A ( tag::vertex );  x(A) = -1.3;  y(A) = -1.3;
	Cell B ( tag::vertex );  x(B) =  1.3;  y(B) = -1.3;
	Cell C ( tag::vertex );  x(C) =  1.3;  y(C) =  1.3;
	Cell D ( tag::vertex );  x(D) = -1.3;  y(D) =  1.3;

	Mesh AB ( tag::segment, A.reverse(), B, tag::divided_in, n );
	Mesh BC ( tag::segment, B.reverse(), C, tag::divided_in, n );
	Mesh CD ( tag::segment, C.reverse(), D, tag::divided_in, n );
	Mesh DA ( tag::segment, D.reverse(), A, tag::divided_in, n );

	Manifold circle = RR2.implicit ( x*x + y*y == 0.7 );
	Mesh inner ( tag::progressive, tag::entire_manifold, circle, tag::desired_length, d );

	Mesh bdry ( tag::join, AB, BC, CD, DA, inner.reverse() );

	RR2.set_as_working_manifold();
	Mesh square ( tag::progressive, tag::boundary, bdry, tag::desired_length, d );

	Mesh torus = square.fold ( tag::identify, AB, tag::with, CD.reverse(),
	                           tag::identify, BC, tag::with, DA.reverse(),
	                           tag::use_existing_vertices                 );

	// declare the type of finite element
	FiniteElement fe ( tag::with_master, tag::triangle, tag::Lagrange, tag::of_degree, 1 );
	Integrator integ = fe.set_integrator ( tag::Gauss, tag::tri_4 );

	std::map < Cell, size_t > numbering;
	{ // just a block of code for hiding 'it' and 'counter'
	CellIterator it = torus.iterator ( tag::over_vertices );
	size_t counter = 0;
	for ( it.reset() ; it.in_range(); it++ )
	{	Cell V = *it;  numbering[V] = counter;  ++counter;  }
	assert ( counter == numbering.size() );
	} // just a block of code

	size_t size_matrix = numbering.size();
	std::cout << "global matrix " << size_matrix + 1 << "x" << size_matrix << std::endl;
	Eigen::SparseMatrix < double > matrix_A ( size_matrix + 1, size_matrix );
	
	matrix_A.reserve ( Eigen::VectorXi::Constant ( size_matrix, 8 ) );
	// since we will be working with a mesh of triangles,
	// there will be, in average, eight non-zero elements per column
	// the diagonal entry plus six neighbour vertices plus the last equation

	Eigen::VectorXd vector_b ( size_matrix + 1 ), vector_sol ( size_matrix );
	vector_b.setZero();

	// run over all square cells composing 'torus'
	{ // just a block of code for hiding 'it'
	CellIterator it = torus.iterator ( tag::over_cells_of_max_dim );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell small_tri = *it;
		std::cout << "*****************************" << std::endl;
		fe.dock_on ( small_tri, tag::spin );
		// run twice over the four vertices of 'small_tri'
		CellIterator it1 = small_tri.boundary().iterator ( tag::over_vertices );
		CellIterator it2 = small_tri.boundary().iterator ( tag::over_vertices );
		for ( it1.reset(); it1.in_range(); it1++ )
		for ( it2.reset(); it2.in_range(); it2++ )
		{	Cell V = *it1, W = *it2;  // V may be the same as W, no problem about that
			// std::cout << "vertices V=(" << x(V) << "," << y(V) << ") " << numbering[V] << ", W=("
			// 					<< x(W) << "," << y(W) << ") " << numbering[W]) << std::endl;
			Function psiV = fe.basis_function(V),
			         psiW = fe.basis_function(W),
			         d_psiV_dx = psiV.deriv(x),
			         d_psiV_dy = psiV.deriv(y),
			         d_psiW_dx = psiW.deriv(x),
			         d_psiW_dy = psiW.deriv(y);
			// 'fe' is already docked on 'small_tri' so this will be the domain of integration
			matrix_A.coeffRef ( numbering[V], numbering[W] ) +=
				fe.integrate ( d_psiV_dx * d_psiW_dx + d_psiV_dy * d_psiW_dy );
		}  }
	} // just a block of code 

}
