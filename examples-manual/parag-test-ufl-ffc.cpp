

#include "maniFEM.h"

using namespace maniFEM;

	
	
int main ()

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy [0], y = xy [1];

	// declare the type of finite element
	FiniteElement fe_ufl_ffc ( tag::quadrangle, tag::Lagrange, tag::of_degree, 1 );
	FiniteElement fe_gauss ( tag::with_master, tag::quadrangle, tag::Lagrange, tag::of_degree, 1 );

	Integrator integ_ufl_ffc = fe_ufl_ffc .set_integrator ( tag::ufl_ffc );
	fe_gauss .set_integrator ( tag::Gauss, tag::quad_9 );

	// UFL-FFC integrators require an early declaration of
	// the integrals we intend to compute later (after docking 'fe' on a cell)
	Function bf1 ( tag::basis_function, tag::within, fe_ufl_ffc ),
	         bf2 ( tag::basis_function, tag::within, fe_ufl_ffc );
	integ_ufl_ffc .pre_compute ( tag::for_given, tag::basis_functions, bf1, bf2,
															 tag::integral_of, { bf1 .deriv(x) * bf2 .deriv(x) } );

	Cell A ( tag::vertex );  x (A) = 0.;   y (A) = 0.;
	Cell B ( tag::vertex );  x (B) = 1.;   y (B) = 0.;
	Cell C ( tag::vertex );  x (C) = 1.;   y (C) = 1.2;
	Cell D ( tag::vertex );  x (D) = 0.;   y (D) = 1.2;
	x (A) = 0.;    y (A) = 0.06;
	x (B) = 1.;    y (B) = 0.;
	x (C) = 1.2;   y (C) = 1.;
	x (D) = 0.;    y (D) = 1.13;

	Cell AB ( tag::segment, A .reverse(), B );
	Cell BC ( tag::segment, B .reverse(), C );
	Cell CD ( tag::segment, C .reverse(), D );
	Cell DA ( tag::segment, D .reverse(), A );
	Cell ABCD ( tag::rectangle, BC, CD, DA, AB );

	{ // just a block of code for hiding names
	fe_gauss .dock_on ( ABCD );
	Function psi_A = fe_gauss .basis_function (A),
	         psi_B = fe_gauss .basis_function (B),
	         psi_C = fe_gauss .basis_function (C),
	         psi_D = fe_gauss .basis_function (D);
	std::cout << "       ";
	std::cout << fe_gauss .integrate ( psi_A .deriv(x) * psi_A .deriv(x) ) << " || ";
	std::cout << fe_gauss .integrate ( psi_A .deriv(x) * psi_B .deriv(x) ) << " || ";
	std::cout << fe_gauss .integrate ( psi_A .deriv(x) * psi_C .deriv(x) ) << " || ";
	std::cout << fe_gauss .integrate ( psi_A .deriv(x) * psi_D .deriv(x) ) << std::endl;
	} { // just a block of code for hiding names
	fe_ufl_ffc .dock_on ( ABCD );
	Function psi_A = fe_ufl_ffc .basis_function (A),
	         psi_B = fe_ufl_ffc .basis_function (B),
	         psi_C = fe_ufl_ffc .basis_function (C),
	         psi_D = fe_ufl_ffc .basis_function (D);
	std::vector < double > result =
		fe_ufl_ffc .integrate ( tag::pre_computed, tag::replace, bf1, tag::by, psi_A,
	                                             tag::replace, bf2, tag::by, psi_A );
	assert ( result .size() == 1 );
	std::cout << result [0] << " || ";
	result = fe_ufl_ffc .integrate
		( tag::pre_computed, tag::replace, bf1, tag::by, psi_A,
	                       tag::replace, bf2, tag::by, psi_B );
	assert ( result .size() == 1 );
	std::cout << result [0] << " || ";
	result = fe_ufl_ffc .integrate
		( tag::pre_computed, tag::replace, bf1, tag::by, psi_A,
	                       tag::replace, bf2, tag::by, psi_C );
	assert ( result .size() == 1 );
	std::cout << result [0] << " || ";
	result = fe_ufl_ffc .integrate
		( tag::pre_computed, tag::replace, bf1, tag::by, psi_A,
	                       tag::replace, bf2, tag::by, psi_D );
	assert ( result .size() == 1 );
	std::cout << result [0] << std::endl;
	} // just a block of code for hiding names

	return 0;
}
