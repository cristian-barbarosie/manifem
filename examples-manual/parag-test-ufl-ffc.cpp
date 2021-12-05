


#include "maniFEM.h"

using namespace maniFEM;

	
	
int main ()

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy [0], y = xy [1];

	// declare the type of finite element
	FiniteElement fe ( tag::triangle, tag::Lagrange, tag::of_degree, 1 );

	Integrator integ = fe .set_integrator ( tag::ufl_ffc );
	// this type of integrator benefits from an early declaration of
	// the integrals we intend to compute later (after docking 'fe' on a cell)
	{ // just a block of code for hiding 'f' and 'vec_fun'
	Function bf ( tag::basis_function, tag::within, fe );
	integ .pre_compute ( tag::for_a_given, tag::basis_function, bf,
                       tag::integral_of, { bf .deriv (x), bf .deriv (y) } );
	} // just a block of code

	Cell A ( tag::vertex );  x (A) = 0.;   y (A) = 0.;
	Cell B ( tag::vertex );  x (B) = 1.;   y (B) = 0.;
	Cell C ( tag::vertex );  x (C) = 1.;   y (C) = 1.2;

	Cell AB ( tag::segment, A .reverse(), B );
	Cell BC ( tag::segment, B .reverse(), C );
	Cell CA ( tag::segment, C .reverse(), A );
	Cell ABC ( tag::triangle, AB, BC, CA );

	fe .dock_on ( ABC );
	Function psi_A = fe .basis_function (A);
	Function psi_B = fe .basis_function (B);
	std::vector < double > result = fe .integrate ( tag::pre_computed, psi_A );
	assert ( result .size() == 2 );
	std::cout << result [0] << " " << result [1] << std::endl;
	result = fe .integrate ( tag::pre_computed, psi_B );
	assert ( result .size() == 2 );
	std::cout << result [0] << " " << result [1] << std::endl;

	return 0;
}
