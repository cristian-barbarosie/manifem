
// example presented in paragraph 6.3 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// rudimentary use of finite elements, 1D integrals

#include "maniFEM.h"
#include "math.h"

#include <fstream>
#include <Eigen/Sparse>

using namespace maniFEM;
using namespace std;

void impose_value_of_unknown
(	Eigen::SparseMatrix < double > & matrix_A, Eigen::VectorXd & vector_b,
	size_t i, double val                                                 )

// in a system of linear equations, destroy equation 'i' and impose u(i) = val
// change also column 'i' of the matrix, just to preserve symmetry

// used for imposing Dirichlet boundary conditions

{	size_t size_matrix = matrix_A .innerSize();
	vector_b (i) = val;
	for ( size_t j = 0; j < size_matrix; j++ )
		matrix_A .coeffRef ( i, j ) = 0.;
	matrix_A .coeffRef ( i, i ) = 1.;
	for ( size_t j = 0; j < size_matrix; j++ )
	{	if ( i == j ) continue;
		vector_b(j) -= matrix_A .coeffRef ( j, i ) * val;
		matrix_A .coeffRef ( j, i ) = 0.;                 }  }

	
int main ()

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy [0], y = xy [1];

	// declare the type of finite element
	FiniteElement fe ( tag::with_master, tag::quadrangle, tag::Lagrange, tag::of_degree, 1 );
	Integrator integ = fe .set_integrator ( tag::Gauss, tag::quad_4 );
	FiniteElement fe_bdry ( tag::with_master, tag::segment, tag::Lagrange, tag::of_degree, 1 );
	Integrator integ_bdry = fe_bdry .set_integrator ( tag::Gauss, tag::seg_3 );

	// build a 10x10 square mesh
	Cell A ( tag::vertex );  x (A) = 0.;   y (A) = 0.;
	Cell B ( tag::vertex );  x (B) = 1.;   y (B) = 0.;
	Cell C ( tag::vertex );  x (C) = 1.;   y (C) = 1.;
	Cell D ( tag::vertex );  x (D) = 0.;   y (D) = 1.;
	Mesh AB ( tag::segment, A .reverse(), B, tag::divided_in, 10 );
	Mesh BC ( tag::segment, B .reverse(), C, tag::divided_in, 12 );
	Mesh CD ( tag::segment, C .reverse(), D, tag::divided_in, 10 );
	Mesh DA ( tag::segment, D .reverse(), A, tag::divided_in, 12 );
	Mesh ABCD ( tag::rectangle, AB, BC, CD, DA );

	// below we use a Cell::Numbering::Map
	// which is essentially an  std::map < Cell, size_t >  disguised
	// the only difference is the syntax for accessing the index of a vertex
	// here we use operator()
	// in order to keep the presentation simple, we do not show
	// this different syntax in paragraph 6.3 of the manual
	// a more efficient numbering is shown in paragraph 6.4 of the manual
	Cell::Numbering::Map numbering;
	{ // just a block of code for hiding 'it' and 'counter'
	Mesh::Iterator it = ABCD .iterator ( tag::over_vertices );
	size_t counter = 0;
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;  numbering (V) = counter;  ++counter;  }
	assert ( counter == numbering .size() );
	} // just a block of code

	size_t size_matrix = ABCD.number_of ( tag::vertices );
	assert ( size_matrix == numbering.size() );
	std::cout << "global matrix " << size_matrix << "x" << size_matrix << std::endl;
	Eigen::SparseMatrix < double > matrix_A ( size_matrix, size_matrix );
	Eigen::VectorXd vector_b ( size_matrix );
	vector_b .setZero();

	// run over all square cells composing ABCD
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = ABCD .iterator ( tag::over_cells_of_max_dim );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell small_square = *it;
		fe .dock_on ( small_square );
		// run twice over the four vertices of 'small_square'
		Mesh::Iterator it1 = small_square .boundary() .iterator ( tag::over_vertices );
		Mesh::Iterator it2 = small_square .boundary() .iterator ( tag::over_vertices );
		for ( it1 .reset(); it1 .in_range(); it1++ )
		for ( it2 .reset(); it2 .in_range(); it2++ )
		{	Cell V = *it1, W = *it2;  // V may be the same as W, no problem about that
			Function psiV = fe .basis_function (V),
			         psiW = fe .basis_function (W),
			         d_psiV_dx = psiV .deriv (x),
			         d_psiV_dy = psiV .deriv (y),
			         d_psiW_dx = psiW .deriv (x),
			         d_psiW_dy = psiW .deriv (y);
			// 'fe' is already docked on 'small_square' so this will be the domain of integration
			matrix_A .coeffRef ( numbering(V), numbering(W) ) +=
				fe .integrate ( d_psiV_dx * d_psiW_dx + d_psiV_dy * d_psiW_dy );  }  }
	} // just a block of code 

	Function heat_source = y*y;
	// impose Neumann boundary conditions du/dn = 1.
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = DA .iterator ( tag::over_segments );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell seg = *it;
		fe_bdry .dock_on ( seg );
		Cell V = seg .base() .reverse();
		assert ( V .is_positive() );
		size_t i = numbering (V);
		Function psiV = fe_bdry .basis_function (V);
		vector_b [i] += fe_bdry .integrate ( heat_source * psiV );
		Cell W = seg .tip();
		assert ( W .is_positive() );
		size_t j = numbering (W);
		Function psiW = fe_bdry .basis_function (W);
		vector_b [j] += fe_bdry .integrate ( heat_source * psiW );   }
	} // just a block of code 
	
	// impose Dirichlet boundary conditions  u = 0
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = AB .iterator ( tag::over_vertices );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell P = *it;
		size_t i = numbering (P);
		impose_value_of_unknown ( matrix_A, vector_b, i, 0. );  }
	} { // just a block of code for hiding 'it' 
	Mesh::Iterator it = BC .iterator ( tag::over_vertices );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell P = *it;
		size_t i = numbering (P);
		impose_value_of_unknown ( matrix_A, vector_b, i, 0. );  }
	} { // just a block of code for hiding 'it'
	Mesh::Iterator it = CD .iterator ( tag::over_vertices );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell P = *it;
		size_t i = numbering (P);
		impose_value_of_unknown ( matrix_A, vector_b, i, 0. );  }
	} // just a block of code 

	// solve the system of linear equations
	Eigen::ConjugateGradient < Eigen::SparseMatrix < double >,
	                           Eigen::Lower|Eigen::Upper    > cg;
	cg .compute ( matrix_A );

	Eigen::VectorXd vector_sol = cg .solve ( vector_b );
	if ( cg .info() != Eigen::Success )
	{	std::cout << "Eigen solver.solve failed" << std::endl;
		exit ( 0 );                                            }
		
	ABCD .export_msh ("square-Neumann.msh", numbering );
	{ // just a block of code for hiding variables
	ofstream solution_file ("square-Neumann.msh", fstream::app );
	solution_file << "$NodeData" << endl;
	solution_file << "1" << endl;   // one string follows
	solution_file << "\"temperature\"" << endl;
	solution_file << "1" << endl;   //  one real follows
	solution_file << "0.0" << endl;  // time [??]
	solution_file << "3" << endl;   // three integers follow
	solution_file << "0" << endl;   // time step [??]
	solution_file << "1" << endl;  // scalar values of u
	solution_file << ABCD .number_of ( tag::vertices ) << endl;  // number of values listed below
	Mesh::Iterator it = ABCD .iterator ( tag::over_vertices );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell P = *it;
		size_t i = numbering (P);
		solution_file << i + 1 << " " << vector_sol[i] << std::endl;   }
	} // just a block of code

	std::cout << "produced file square-Neumann.msh" << std::endl;

	return 0;
}

