

#include "maniFEM.h"

using namespace maniFEM;

	
	
int main ()

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy [0], y = xy [1];

	// declare the type of finite element
	FiniteElement fe_ufl_ffc ( tag::triangle, tag::Lagrange, tag::of_degree, 1 );
	FiniteElement fe_gauss ( tag::with_master, tag::triangle, tag::Lagrange, tag::of_degree, 1 );

	Integrator integ_ufl_ffc = fe_ufl_ffc .set_integrator ( tag::ufl_ffc );
	Integrator integ_gauss = fe_gauss .set_integrator ( tag::Gauss, tag::tri_3 );

	// UFL-FFC integrators require an early declaration of
	// the integrals we intend to compute later (after docking 'fe' on a cell)
	{ // just a block of code for hiding 'bf'
	Function bf ( tag::basis_function, tag::within, fe_ufl_ffc );
	integ_ufl_ffc .pre_compute ( tag::for_a_given, tag::basis_function, bf,
                               tag::integral_of, { bf .deriv(y) * bf .deriv(y) } );
	} // just a block of code

	Cell A ( tag::vertex );  x (A) = -0.13;   y (A) = 0.191;
	Cell B ( tag::vertex );  x (B) =  1.;     y (B) = 0.;
	Cell C ( tag::vertex );  x (C) =  1.01;   y (C) = 1.2;

	Cell AB ( tag::segment, A .reverse(), B );
	Cell BC ( tag::segment, B .reverse(), C );
	Cell CA ( tag::segment, C .reverse(), A );
	Cell ABC ( tag::triangle, AB, BC, CA );

	{ // just a block of code for hiding names
	fe_gauss .dock_on ( ABC );
	Function psi_A = fe_gauss .basis_function (A);
	Function psi_B = fe_gauss .basis_function (B);
	Function psi_C = fe_gauss .basis_function (C);
	std::cout << fe_gauss .integrate ( psi_A .deriv(y) * psi_B .deriv(y) ) << " || ";
	std::cout << fe_gauss .integrate ( psi_A .deriv(y) * psi_C .deriv(y) ) << " || ";
	std::cout << fe_gauss .integrate ( psi_B .deriv(y) * psi_C .deriv(y) ) << " || ";
	std::cout << fe_gauss .integrate ( psi_A .deriv(y) * psi_A .deriv(y) ) << " || ";
	std::cout << fe_gauss .integrate ( psi_B .deriv(y) * psi_A .deriv(y) ) << " || ";
	std::cout << fe_gauss .integrate ( psi_C .deriv(y) * psi_C .deriv(y) ) << std::endl;
	} { // just a block of code for hiding names
	fe_ufl_ffc .dock_on ( ABC );
	Function psi_A = fe_ufl_ffc .basis_function (A);
	Function psi_B = fe_ufl_ffc .basis_function (B);
	Function psi_C = fe_ufl_ffc .basis_function (C);
	Integrator::Result result =
		fe_ufl_ffc .integrate ( tag::pre_computed, psi_A, psi_B );
	assert ( result .size() == 1 );
	std::cout << result [0] << " || ";
	result = fe_ufl_ffc .integrate ( tag::pre_computed, psi_A, psi_C );
	assert ( result .size() == 1 );
	std::cout << result [0] << " || ";
	result = fe_ufl_ffc .integrate ( tag::pre_computed, psi_B, psi_C );
	assert ( result .size() == 1 );
	std::cout << result [0] << " || ";
	result = fe_ufl_ffc .integrate ( tag::pre_computed, psi_A, psi_A );
	assert ( result .size() == 1 );
	std::cout << result [0] << " || ";
	result = fe_ufl_ffc .integrate ( tag::pre_computed, psi_B, psi_A );
	assert ( result .size() == 1 );
	std::cout << result [0] << " || ";
	result = fe_ufl_ffc .integrate ( tag::pre_computed, psi_C, psi_C );
	assert ( result .size() == 1 );
	std::cout << result [0] << std::endl;
	} // just a block of code for hiding names

	return 0;
}
