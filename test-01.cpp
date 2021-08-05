

#include "maniFEM.h"
#include "math.h"

#include <fstream>

using namespace maniFEM;

	
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

	std::cout << "total number of vertices " << ABCD.number_of ( tag::vertices );
	std::cout << "we have " << fe.numbers[0]->size() << " numbered vertices" <<  std::endl;

	// there will be a more elegant and efficient way of producing the numbering
	std::map < Cell::Core *, size_t > numbering;
	{ // just a block of code for hiding 'it' and 'counter'
	CellIterator it = ABCD.iterator ( tag::over_vertices );
	size_t counter = 0;
	for ( it.reset() ; it.in_range(); it++ )
	{	Cell p = *it;  ++counter;  numbering [p.core] = counter;  }
	} // just a block of code 
	
	// run over all square cells composing ABCD
	{ // just a block of code for hiding 'it'
	int counter = 0;	
	CellIterator it = ABCD.iterator ( tag::over_cells_of_max_dim );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell small_square = *it;
		fe.dock_on ( small_square );
		// run twice over the four vertices of 'small_square'
		CellIterator it1 = small_square.boundary().iterator ( tag::over_vertices );
		CellIterator it2 = small_square.boundary().iterator ( tag::over_vertices );
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
			//matrix_A.coeffRef ( numbering[V.core]-1, numbering[W.core]-1 ) +=
				fe.integrate ( d_psiV_dx * d_psiW_dx + d_psiV_dy * d_psiW_dy );
		}  }
	} // just a block of code 

	ABCD.export_msh ("square-Dirichlet.msh", * fe.numbers[0] );
	std::cout << "produced file square-Dirichlet.msh" << std::endl;

	return 0;
}
