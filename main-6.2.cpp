
// example presented in paragraph 6.2 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// rudimentary use of finite elements

#include "maniFEM.h"
#include "math.h"

#include <fstream>
#include <Eigen/Sparse>

using namespace maniFEM;
using namespace std;

void impose_value_of_unknown
(	Eigen::SparseMatrix <double> & matrix_A, Eigen::VectorXd & vector_b,
	size_t i, double val                                                 )

// in a system of linear equations, destroy equation 'i' and impose u(i) = val
// change also column 'i' of the matrix, just to preserve symmetry of the matrix

// used for imposing Dirichlet boundary conditions

{	size_t size_matrix = matrix_A.innerSize();
	vector_b(i) = val;
	for ( size_t j = 0; j < size_matrix; j++ )
		matrix_A.coeffRef ( i, j ) = 0.;
	matrix_A.coeffRef ( i, i ) = 1.;
	for ( size_t j = 0; j < size_matrix; j++ )
	{	if ( i == j ) continue;
		vector_b(j) -= matrix_A.coeffRef ( j, i ) * val;
		matrix_A.coeffRef ( j, i ) = 0.;                  }  }

	
int main ()

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0],  y = xy[1];

	// build a 10x10 square mesh
	Cell A ( tag::vertex );  x(A) = 0.;   y(A) = 0.;
	Cell B ( tag::vertex );  x(B) = 1.;   y(B) = 0.;
	Cell C ( tag::vertex );  x(C) = 1.;   y(C) = 1.;
	Cell D ( tag::vertex );  x(D) = 0.;   y(D) = 1.;
	Mesh AB ( tag::segment, A.reverse(), B, tag::divided_in, 10 );
	Mesh BC ( tag::segment, B.reverse(), C, tag::divided_in, 12 );
	Mesh CD ( tag::segment, C.reverse(), D, tag::divided_in, 10 );
	Mesh DA ( tag::segment, D.reverse(), A, tag::divided_in, 12 );
	Mesh ABCD ( tag::rectangle, AB, BC, CD, DA );

	// declare the type of finite element
	FiniteElement fe ( tag::with_master, tag::quadrangle, tag::Lagrange, tag::of_degree, 1 );
	Integrator integ = fe.set_integrator ( tag::Gauss, tag::quad_4 );

	// there will be a more elegant and efficient way of producing the numbering
	std::map < Cell::Core *, size_t > numbering;
	{ // just a block of code for hiding 'it' and 'counter'
	CellIterator it = ABCD.iter_over ( tag::vertices );
	size_t counter = 0;
	for ( it.reset() ; it.in_range(); it++ )
	{	Cell p = *it;  ++counter;  numbering [p.core] = counter;  }
	} // just a block of code 
	
	size_t size_matrix = ABCD.number_of ( tag::vertices );
	assert ( size_matrix == numbering.size() );
	std::cout << "global matrix " << size_matrix << "x" << size_matrix << std::endl;
	Eigen::SparseMatrix <double> matrix_A ( size_matrix, size_matrix );
	Eigen::VectorXd vector_b ( size_matrix ), vector_sol ( size_matrix );
	vector_b.setZero();

	Function::name[x.core] = "x";
	Function::name[y.core] = "y";

	// run over all square cells composing ABCD
	{ // just a block of code for hiding 'it'
	int counter = 0;	
	CellIterator it = ABCD.iter_over ( tag::cells_of_dim, 2 );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell small_square = *it;
		fe.dock_on ( small_square );
		// run twice over the four vertices of 'small_square'
		CellIterator it1 = small_square.boundary().iter_over ( tag::vertices );
		CellIterator it2 = small_square.boundary().iter_over ( tag::vertices );
		for ( it1.reset(); it1.in_range(); it1++ )
		for ( it2.reset(); it2.in_range(); it2++ )
		{	Cell V = *it1, W = *it2;  // V may be the same as W, no problem about that
			// std::cout << "vertices V=(" << x(V) << "," << y(V) << ") " << numbering[V.core] << ", W=("
			// 					<< x(W) << "," << y(W) << ") " << numbering[W.core] << std::endl;
			Function psiV = fe.basis_function(V),
			         psiW = fe.basis_function(W),
			         d_psiV_dx = psiV.deriv(x),
			         d_psiV_dy = psiV.deriv(y),
			         d_psiW_dx = psiW.deriv(x),
			         d_psiW_dy = psiW.deriv(y);
			// 'fe' is already docked on 'small_square' so this will be the domain of integration
			matrix_A.coeffRef ( numbering[V.core]-1, numbering[W.core]-1 ) +=
				fe.integrate ( d_psiV_dx * d_psiW_dx + d_psiV_dy * d_psiW_dy );
		}  }
	} // just a block of code 

	// impose Dirichlet boundary conditions  u = xy
	{ // just a block of code for hiding 'it'
	CellIterator it = AB.iter_over ( tag::vertices );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell P = *it;
		size_t i = numbering[P.core]-1;
		impose_value_of_unknown ( matrix_A, vector_b, i, 0. );  }
	} { // just a block of code for hiding 'it' 
	CellIterator it = BC.iter_over ( tag::vertices );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell P = *it;
		size_t i = numbering[P.core]-1;
		impose_value_of_unknown ( matrix_A, vector_b, i, y(P) );  }
	} { // just a block of code for hiding 'it'
	CellIterator it = CD.iter_over ( tag::vertices );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell P = *it;
		size_t i = numbering[P.core]-1;
		impose_value_of_unknown ( matrix_A, vector_b, i, x(P) );  }
	} { // just a block of code for hiding 'it'
	CellIterator it = DA.iter_over ( tag::vertices );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell P = *it;
		size_t i = numbering[P.core]-1;
		impose_value_of_unknown ( matrix_A, vector_b, i, 0. );  }
	} // just a block of code 
	
//	for ( size_t i = 0; i < size_matrix; i++ )
//	{	std::cout << i << " ";
//		for ( size_t j = 0; j < size_matrix; j++ )
//			std::cout << matrix_A.coeffRef (i,j) << " ";
//		std::cout << std::endl;                           }

//	std::cout << std::endl;
//	for ( size_t j = 0; j < size_matrix; j++ )
//		std::cout << vector_b (j) << " ";
//	std::cout << std::endl;
//	for ( size_t i = 0; i < size_matrix; i++ )
//		std::cout << vector_sol (i) << " ";
//	std::cout << std::endl;

	// solve the system of linear equations
	Eigen::ConjugateGradient < Eigen::SparseMatrix<double>,
	                           Eigen::Lower|Eigen::Upper    > cg;
	cg.compute ( matrix_A );
	vector_sol = cg.solve ( vector_b );

	ABCD.export_msh ("square-Dirichlet.msh", numbering );
	{ // just a block of code for hiding variables
	ofstream solution_file ("square-Dirichlet.msh", fstream::app );
	solution_file << "$NodeData" << endl;
	solution_file << "1" << endl;   // one string follows
	solution_file << "\"temperature\"" << endl;
	solution_file << "1" << endl;   //  one real follows
	solution_file << "0.0" << endl;  // time [??]
	solution_file << "3" << endl;   // three integers follow
	solution_file << "0" << endl;   // time step [??]
	solution_file << "1" << endl;  // scalar values of u
	solution_file << ABCD.number_of ( tag::vertices ) << endl;  // number of values listed below
	CellIterator it = ABCD.iter_over ( tag::vertices );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell P = *it;
		size_t i = numbering[P.core];
		solution_file << i << " " << vector_sol[i-1] << std::endl;   }
	} // just a block of code

	std::cout << "produced file square-Dirichlet.msh" << std::endl;

	return 0;
}

/*

		FiniteElement::WithMaster * fe_core =
			dynamic_cast < FiniteElement::WithMaster * > ( fe.core );
		assert ( fe_core );
		Function::Diffeomorphism * tran = Function::core_to_diffeom ( fe_core->transf.core );
		Function::name[tran->master_coords[0].core] = "xi";
		Function::name[tran->master_coords[1].core] = "eta";

			std::cout << "psiV : " << psiV.repr() << std::endl;
			std::cout << "psiW : " << psiW.repr() << std::endl;

	std::cout << "d psiV / dx : " << d_psiV_dx.repr() << std::endl;
	std::cout << "d psiV / dy : " << d_psiV_dy.repr() << std::endl;
	std::cout << "d psiW / dx : " << d_psiW_dx.repr() << std::endl;
	std::cout << "d psiW / dy : " << d_psiW_dy.repr() << std::endl;


*/	

