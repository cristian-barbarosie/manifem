


#include "maniFEM.h"

using namespace maniFEM;

		
int main ()

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy [0], y = xy [1];

	FiniteElement fe ( tag::with_master, tag::triangle,
	                   tag::Lagrange, tag::of_degree, 1, tag::enumerate_cells );
	Integrator integ = fe .set_integrator ( tag::Gauss, tag::tri_6 );

	// build a 10x10 square mesh
	Cell A ( tag::vertex );  x (A) = 0.;   y (A) = 0.;
	Cell B ( tag::vertex );  x (B) = 1.;   y (B) = 0.1;
	Cell C ( tag::vertex );  x (C) = 1.;   y (C) = 1.2;

	Cell AB ( tag::segment, A .reverse(), B );
	Cell BC ( tag::segment, B .reverse(), C );
	Cell CA ( tag::segment, C .reverse(), A );
	Cell ABC ( tag::triangle, AB, BC, CA );

	fe .dock_on ( ABC );
	Function psi_A = fe .basis_function (A);
	Function psi_B = fe .basis_function (B);
	std::cout << fe .integrate ( psi_A .deriv (x) ) << " "
            << fe .integrate ( psi_A .deriv (y) ) << std::endl;
	std::cout << fe .integrate ( psi_B .deriv (x) ) << " "
            << fe .integrate ( psi_B .deriv (y) ) << std::endl;

	return 0;
}
