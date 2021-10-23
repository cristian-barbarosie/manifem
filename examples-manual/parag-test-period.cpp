
// solve a celullar problem
// square periodicity, triangular elements, circular hole
// identification of opposite sides based order (when exporting msh)


#include "maniFEM.h"
#include <fstream>
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

	// we number all nodes in 'square', not only those belonging to 'torus'
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
	
	matrix_A.reserve ( Eigen::VectorXi::Constant ( size_matrix, 8 ) );
	// since we will be working with a mesh of triangles,
	// there will be, in average, eight non-zero elements per column
	// the diagonal entry plus six neighbour vertices plus the last equation

	// we fill the main diagonal with ones
	// then we put zero for vertices belonging to 'torus'
	{ // just a block of code for hiding 'it'
	for ( size_t i = 0; i < size_matrix; i++ ) matrix_A.coeffRef ( i, i ) = 1.;
	CellIterator it = torus.iterator ( tag::over_vertices );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell V = *it;
		matrix_A.coeffRef ( numbering[V], numbering[V] ) = 0.;  }
	} // just a block of code for hiding 'it'
	
	Eigen::VectorXd vector_b ( size_matrix + 1 ), vector_sol ( size_matrix );
	vector_b.setZero();

	xy = Manifold::working.coordinates();
	x = xy[0];  y = xy[1];

	// macroscopic temperature gradient
	std::vector < double > macro_grad { 1., 0.5 };
	Function::Jump jump_of_solution = macro_grad[0] * x.jump() + macro_grad[1] * y.jump();

	// run over all square cells composing 'torus'
	{ // just a block of code for hiding 'it'
	CellIterator it = torus.iterator ( tag::over_cells_of_max_dim );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell small_tri = *it;
		fe.dock_on ( small_tri, tag::spin );
		// run twice over the four vertices of 'small_tri'
		CellIterator it_V = small_tri.boundary().iterator ( tag::over_vertices );
		for ( it_V.reset(); it_V.in_range(); it_V++ )
		{	Cell V = *it_V;
			// perhaps implement an interator returning a vertex and a segment
			Cell seg = small_tri .boundary(). cell_in_front_of ( V );
			Cell W = V;
			double jump_V_W = 0.;
			while ( true )
			{	assert ( W == seg.base().reverse() );
				// V may be the same as W, no problem about that
				Function psi_V = fe .basis_function ( V ),
				         psi_W = fe .basis_function ( W ),
				         d_psi_V_dx = psi_V .deriv ( x ),
				         d_psi_V_dy = psi_V .deriv ( y ),
				         d_psi_W_dx = psi_W .deriv ( x ),
				         d_psi_W_dy = psi_W .deriv ( y );
				// 'fe' is already docked on 'small_tri' so this will be the domain of integration
				double integral = fe.integrate ( d_psi_V_dx * d_psi_W_dx + d_psi_V_dy * d_psi_W_dy );
				matrix_A.coeffRef ( numbering[V], numbering[W] ) += integral;
				vector_b ( numbering[V] ) -= jump_V_W * integral;
				jump_V_W += jump_of_solution ( seg.spin() );
				W = seg.tip();
				if ( V == W ) break;
				seg = small_tri .boundary() .cell_in_front_of ( seg.tip() );                          }
			// here  jump_V_W  should be zero again
			// but we do not assert that, rounding errors may mess up things
		}  }
	} // just a block of code for hiding 'it'

	// we add, as last equation, the condition of zero average
	// actually, a more rudimentary condition : a line of ones
	for ( size_t i = 0; i < size_matrix; i++ ) matrix_A.coeffRef ( size_matrix, i ) = 1.;

	matrix_A .makeCompressed();

	Eigen::SparseQR < Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;

	solver.compute ( matrix_A );
	if ( solver.info() != Eigen::Success )
	{	std::cout << "Eigen solver.compute failed" << std::endl;
		exit ( 0 );                                  }

	vector_sol = solver.solve ( vector_b );
	if ( solver.info() != Eigen::Success )
	{	std::cout << "Eigen solver.solve failed" << std::endl;
		exit ( 0 );                                  }

	RR2 .set_as_working_manifold();
	xy = Manifold::working.coordinates();
	x = xy[0];  y = xy[1];
	square.export_msh ("cell.msh", numbering );
	
	{ // just a block of code for hiding variables
	std::ofstream solution_file ("cell.msh", std::fstream::app );
	solution_file << "$NodeData" << std::endl;
	solution_file << "1" << std::endl;   // one string follows
	solution_file << "\"temperature\"" << std::endl;
	solution_file << "1" << std::endl;   //  one real follows
	solution_file << "0.0" << std::endl;  // time [??]
	solution_file << "3" << std::endl;   // three integers follow
	solution_file << "0" << std::endl;   // time step [??]
	solution_file << "1" << std::endl;  // scalar values of u
	solution_file << square .number_of ( tag::vertices ) << std::endl;
  // number of values listed below
	CellIterator it = torus.iterator ( tag::over_vertices );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell P = *it;
		size_t j = numbering [ P ];
		solution_file << j << " " << vector_sol[j] << std::endl;  }
	size_t j = numbering [ B ];
	assert ( not A .belongs_to ( torus ) );
	solution_file << numbering[A] << " "
		<< vector_sol[j] + macro_grad[0] * ( x ( A ) - x ( B ) )
		                 + macro_grad[1] * ( y ( A ) - y ( B ) ) << std::endl;
	assert ( not C .belongs_to ( torus ) );
	solution_file << numbering[C] << " "
		<< vector_sol[j] + macro_grad[0] * ( x ( C ) - x ( B ) )
		                 + macro_grad[1] * ( y ( C ) - y ( B ) ) << std::endl;
	assert ( not D .belongs_to ( torus ) );
	solution_file << numbering[D] << " "
		<< vector_sol[j] + macro_grad[0] * ( x ( D ) - x ( B ) )
		                 + macro_grad[1] * ( y ( D ) - y ( B ) ) << std::endl;
	CellIterator it_AB = AB .iterator ( tag::over_vertices, tag::require_order );
	CellIterator it_CD = CD .iterator ( tag::over_vertices, tag::backwards );
	it_AB .reset();  assert ( it_AB .in_range() );
	assert ( *it_AB == A );  it_AB ++;
	it_CD .reset();  assert ( it_CD .in_range() );
	assert ( *it_CD == D );  it_CD ++;
	for ( ; ; it_AB ++, it_CD ++ )
	{	assert ( it_AB .in_range() );  assert ( it_CD .in_range() );
		Cell V = *it_AB, W = *it_CD;
		if ( V == B )  {  assert ( W == C );  break;  }
		assert ( not W .belongs_to ( torus ) );
		j = numbering [ V ];
		solution_file << numbering[W] << " "
			<< vector_sol[j] + macro_grad[0] * ( x ( W ) - x ( V ) )
			                 + macro_grad[1] * ( y ( W ) - y ( V ) ) << std::endl;  }
	CellIterator it_BC = BC .iterator ( tag::over_vertices, tag::require_order );
	CellIterator it_DA = DA .iterator ( tag::over_vertices, tag::backwards );
	it_BC .reset();  assert ( it_BC .in_range() );
	assert ( *it_BC == B );  it_BC ++;
	it_DA .reset();  assert ( it_DA .in_range() );
	assert ( *it_DA == A );  it_DA ++;
	for ( ; ; it_BC ++, it_DA ++ )
	{	assert ( it_BC .in_range() );  assert ( it_DA .in_range() );
		Cell V = *it_BC, W = *it_DA;
		if ( V == C )  {  assert ( W == D );  break;  }
		assert ( not W .belongs_to ( torus ) );
		j = numbering [ V ];
		solution_file << numbering[W] << " "
			<< vector_sol[j] + macro_grad[0] * ( x ( W ) - x ( V ) )
			                 + macro_grad[1] * ( y ( W ) - y ( V ) ) << std::endl;  }
	} // just a block of code

	std::cout << "produced file cell.msh" << std::endl;

	return 0;

}
