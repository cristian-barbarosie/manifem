
// example presented in paragraph 6.5 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// comparison between Lagrange P1 and P2 finite elements, triangular
// solve the Laplace operator on a disk with Dirichlet boundary conditions

#include "maniFEM.h"
#include "math.h"

#include <fstream>
#include <Eigen/Sparse>

using namespace maniFEM;

	
void impose_value_of_unknown
(	Eigen::SparseMatrix < double > & matrix_A, Eigen::VectorXd & vector_b,
	size_t i, double val                                                   )

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

//---------------------------------------------------------------------------------------------

	
int main ()

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy [0], y = xy [1];
	Function r2 = x*x + y*y;

	Manifold circle_manif = RR2 .implicit ( x*x + y*y == 1. );

	Cell A ( tag::vertex );  x (A) = 1.;  y (A) = 0.;
	Mesh circle ( tag::frontal, tag::start_at, A, tag::desired_length, 0.1 );

	RR2 .set_as_working_manifold();
	Mesh disk ( tag::frontal, tag::boundary, circle, tag::desired_length, 0.1 );

	// we solve the elliptic equation
	//   - lapl u + 8 (x^2+y^2) / (1+x^2+y^2)^2 u = 4 / (1+x^2+y^2)^2
	// whose exact solution is  1 / (1+x^2+y^2)
	// whose integral is  pi ln 2
	
	const double exact = 4. * std::atan (1.) * std::log (2.);
	std::cout << "integral, exact " << exact << std::endl;

	Function alpha = 8. * r2 / (1+r2) / (1+r2),
	         beta  = 4. / (1+r2) / (1+r2);

	// declare the type of finite element
	FiniteElement fe1 ( tag::with_master, tag::triangle, tag::Lagrange, tag::of_degree, 1 );
	fe1 .set_integrator ( tag::Gauss, tag::tri_6 );

	// a different solution for numbering vertices is shown in paragraph 6.4 of the manual
	// below we use an std::map<Cell,size_t>
	std::map < Cell, size_t > numbering;
	{ // just a block of code for hiding 'it' and 'counter'
	size_t counter = 0;
	Mesh::Iterator it_ver = disk .iterator ( tag::over_vertices );
	for ( it_ver .reset() ; it_ver .in_range(); it_ver ++ )
	{	Cell V = *it_ver;  numbering [V] = counter;  ++ counter;  }
	assert ( counter == numbering .size() );
	} // just a block of code

	size_t size_matrix = numbering .size();
	std::cout << "global matrix " << size_matrix << "x" << size_matrix << std::endl;
	Eigen::SparseMatrix <double> matrix_A1 ( size_matrix, size_matrix );
	
	matrix_A1 .reserve ( Eigen::VectorXi::Constant ( size_matrix, 7 ) );
	// since we will be working with a mesh of triangles,
	// there will be about 7 non-zero elements per column
	// the diagonal entry plus six neighbour vertices

	Eigen::VectorXd vector_b1 ( size_matrix );
	vector_b1 .setZero();

	// run over all triangular cells composing disk
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = disk .iterator ( tag::over_cells_of_max_dim );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell small_tri = *it;
		fe1 .dock_on ( small_tri );

		// run twice over the three vertices of 'small_tri'
		Mesh::Iterator it_V = small_tri .boundary() .iterator ( tag::over_vertices );
		for ( it_V .reset(); it_V .in_range(); it_V ++ )
		{	Cell V = * it_V;
			// perhaps implement an iterator returning a vertex and a segment
			Function psi_V = fe1 .basis_function ( V ),
			         d_psi_V_dx = psi_V .deriv ( x ),
			         d_psi_V_dy = psi_V .deriv ( y );

			Cell W = V;
			while ( true )  // V may be the same as W, no problem about that
			{	Function psi_W = fe1 .basis_function ( W ),
				         d_psi_W_dx = psi_W .deriv ( x ),
				         d_psi_W_dy = psi_W .deriv ( y );

				// 'fe' is already docked on 'small_tri' so this will be the domain of integration
				matrix_A1 .coeffRef ( numbering[V], numbering[W] ) +=
					fe1 .integrate ( d_psi_V_dx * d_psi_W_dx + d_psi_V_dy * d_psi_W_dy ) +
					fe1 .integrate ( alpha * psi_V * psi_W );

				Cell seg = small_tri .boundary() .cell_in_front_of ( W );
				W = seg.tip();
				if ( V == W ) break;
			}  // end of while
			
			vector_b1 ( numbering[V] ) += fe1 .integrate ( beta * psi_V );
			
	}	}  // end of two for loops
	} // just a block of code 

	// impose Dirichlet boundary conditions  u = 0.5
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = circle .iterator ( tag::over_vertices );
	for ( it .reset(); it .in_range(); it++ )
	// perhaps implement an iterator returning a vertex and a segment
	{	Cell P = *it;
		size_t i = numbering [P];
		impose_value_of_unknown ( matrix_A1, vector_b1, i, 0.5 );          }
	} // just a block of code 
	
	// solve the system of linear equations
	Eigen::ConjugateGradient < Eigen::SparseMatrix < double >,
	                           Eigen::Lower | Eigen::Upper    > cg;
	cg .compute ( matrix_A1 );

	Eigen::VectorXd vector_sol = cg .solve ( vector_b1 );
	if ( cg .info() != Eigen::Success )
	{	std::cout << "Eigen solver.solve failed" << std::endl;
		exit ( 0 );                                            }

	// now compute the integral of u
	double integral = 0.;
	// run over all triangular cells composing disk
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = disk .iterator ( tag::over_cells_of_max_dim );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell small_tri = *it;
		fe1 .dock_on ( small_tri );
		// run twice over the three vertices of 'small_tri'
		Mesh::Iterator it_V = small_tri .boundary() .iterator ( tag::over_vertices );
		for ( it_V .reset(); it_V .in_range(); it_V ++ )
		{	Cell V = * it_V;
			Function psi_V = fe1 .basis_function ( V );
			size_t i = numbering [V];
			integral += vector_sol (i) * fe1 .integrate ( psi_V );			
	}	}  // end of two for loops
	} // just a block of code

	std::cout << "integral computed with degree 1 finite elements " << integral
						<< ", error " << std::abs ( integral - exact ) << std::endl;

	// declare the type of finite element
	FiniteElement fe2 ( tag::with_master, tag::triangle,
	                    tag::Lagrange, tag::of_degree, 2, tag::straight );
	fe2 .set_integrator ( tag::Gauss, tag::tri_6 );

	{ // just a block of code for hiding 'it' and 'counter'
	size_t counter = numbering .size();
	Mesh::Iterator it_seg = disk .iterator ( tag::over_segments );  // positive segments
	for ( it_seg .reset() ; it_seg .in_range(); it_seg ++ )
	{	Cell seg = *it_seg;  numbering [ seg ] = counter;  ++ counter;  }
	assert ( counter == numbering .size() );
	} // just a block of code

	size_matrix = numbering .size();
	std::cout << "global matrix " << size_matrix << "x" << size_matrix << std::endl;
	Eigen::SparseMatrix <double> matrix_A2 ( size_matrix, size_matrix );
	
	matrix_A2 .reserve ( Eigen::VectorXi::Constant ( size_matrix, 12 ) );
	// since we will be working with a mesh of triangles,
	// there will be about 12 non-zero elements per column
	// a vertex has six neighbour vertices plus twelve neighbour segments
	// a segment has four neighbour vertices plus four neighbour segments
	// there are three times more segments than vertices

	Eigen::VectorXd vector_b2 ( size_matrix );
	vector_b2 .setZero();

	// run over all triangular cells composing disk
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = disk .iterator ( tag::over_cells_of_max_dim );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell small_tri = *it;
		fe2 .dock_on ( small_tri );

		// run twice over the three vertices of 'small_tri'
		Mesh::Iterator it_V = small_tri .boundary() .iterator ( tag::over_vertices );
		for ( it_V .reset(); it_V .in_range(); it_V ++ )
		{	Cell V = * it_V;
			// perhaps implement an iterator returning a vertex and a segment
			Function psi_V = fe2 .basis_function ( V ),
			         d_psi_V_dx = psi_V .deriv ( x ),
			         d_psi_V_dy = psi_V .deriv ( y );
			Cell seg = small_tri .boundary() .cell_in_front_of ( V );
			Function psi_seg_V = fe2 .basis_function ( seg ),
			         d_psi_seg_V_dx = psi_seg_V .deriv ( x ),
			         d_psi_seg_V_dy = psi_seg_V .deriv ( y );
			Cell sV = seg .get_positive();
			Cell W = V;
			while ( true )  // V may be the same as W, no problem about that
			{	assert ( W == seg .base() .reverse() );
				Function psi_W = fe2 .basis_function ( W ),
				         d_psi_W_dx = psi_W .deriv ( x ),
				         d_psi_W_dy = psi_W .deriv ( y );
				Function psi_seg_W = fe2 .basis_function ( seg ),
				         d_psi_seg_W_dx = psi_seg_W .deriv ( x ),
				         d_psi_seg_W_dy = psi_seg_W .deriv ( y );
				Cell sW = seg .get_positive();

				// 'fe' is already docked on 'small_tri' so this will be the domain of integration
				matrix_A2 .coeffRef ( numbering[V], numbering[W] ) +=
					fe2 .integrate ( d_psi_V_dx * d_psi_W_dx + d_psi_V_dy * d_psi_W_dy ) +
					fe2 .integrate ( alpha * psi_V * psi_W );
				matrix_A2 .coeffRef ( numbering[V], numbering[sW] ) +=
					fe2 .integrate ( d_psi_V_dx * d_psi_seg_W_dx + d_psi_V_dy * d_psi_seg_W_dy ) +
					fe2 .integrate ( alpha * psi_V * psi_seg_W );
				matrix_A2 .coeffRef ( numbering[sV], numbering[W] ) +=
					fe2 .integrate ( d_psi_seg_V_dx * d_psi_W_dx + d_psi_seg_V_dy * d_psi_W_dy ) +
					fe2 .integrate ( alpha * psi_seg_V * psi_W );
				matrix_A2 .coeffRef ( numbering[sV], numbering[sW] ) +=
					fe2 .integrate ( d_psi_seg_V_dx * d_psi_seg_W_dx + d_psi_seg_V_dy * d_psi_seg_W_dy ) +
					fe2 .integrate ( alpha * psi_seg_V * psi_seg_W );

				W = seg.tip();
				if ( V == W ) break;
				seg = small_tri .boundary() .cell_in_front_of ( W );
			}  // end of while
			
			vector_b2 ( numbering[V] ) += fe2 .integrate ( beta * psi_V );
			vector_b2 ( numbering[sV] ) += fe2 .integrate ( beta * psi_seg_V );
			
	}	}  // end of two for loops
	} // just a block of code 

	// impose Dirichlet boundary conditions  u = 0.5
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = circle .iterator ( tag::over_vertices );
	for ( it .reset(); it .in_range(); it++ )
	// perhaps implement an iterator returning a vertex and a segment
	{	Cell P = *it;
		size_t i = numbering [P];
		impose_value_of_unknown ( matrix_A2, vector_b2, i, 0.5 );
		Cell seg = circle .cell_in_front_of ( P, tag::surely_exists );
		i = numbering [ seg .get_positive() ];
		impose_value_of_unknown ( matrix_A2, vector_b2, i, 0.5 );        }
	} // just a block of code 
	
	// solve the system of linear equations
	cg .compute ( matrix_A2 );

	vector_sol = cg .solve ( vector_b2 );
	if ( cg .info() != Eigen::Success )
	{	std::cout << "Eigen solver.solve failed" << std::endl;
		exit ( 0 );                                            }

	// now compute the integral of u
	integral = 0.;
	// run over all triangular cells composing disk
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = disk .iterator ( tag::over_cells_of_max_dim );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell small_tri = *it;
		fe2 .dock_on ( small_tri );
		// run twice over the three vertices of 'small_tri'
		Mesh::Iterator it_V = small_tri .boundary() .iterator ( tag::over_vertices );
		for ( it_V .reset(); it_V .in_range(); it_V ++ )
		{	Cell V = * it_V;
			// perhaps implement an iterator returning a vertex and a segment
			Function psi_V = fe2 .basis_function ( V );
			size_t i = numbering [V];
			integral += vector_sol (i) * fe2 .integrate ( psi_V );
			Cell seg = small_tri .boundary() .cell_in_front_of (V);
			i = numbering [ seg .get_positive() ];
			Function psi_seg = fe2 .basis_function ( seg );
			integral += vector_sol (i) * fe2 .integrate ( psi_seg );
	}	}  // end of two for loops
	} // just a block of code

	std::cout << "integral computed with degree 2 finite elements " << integral
						<< ", error " << std::abs ( integral - exact ) << std::endl;

}  // end of main

//---------------------------------------------------------------------------------------------


