

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
	fe_gauss .set_integrator ( tag::Gauss, tag::tri_3 );

	// UFL-FFC integrators require an early declaration of
	// the integrals we intend to compute later (after docking 'fe' on a cell)
	Function bf1 ( tag::basis_function, tag::within, fe_ufl_ffc ),
	         bf2 ( tag::basis_function, tag::within, fe_ufl_ffc );
	integ_ufl_ffc .pre_compute ( tag::for_given, tag::basis_functions, bf1, bf2,
	    tag::integral_of, { bf1 .deriv(x) * bf2,
											    bf1 .deriv(x) * bf2 .deriv(x) + bf1 .deriv(y) * bf2 .deriv(y),
                          bf1 * bf2, bf2, bf2 .deriv(x), bf1 .deriv(x) * bf2 .deriv(x)  } );

	Cell A ( tag::vertex );  x (A) = -0.13;   y (A) = 0.191;
	Cell B ( tag::vertex );  x (B) =  1.;     y (B) = 0.;
	Cell C ( tag::vertex );  x (C) =  1.01;   y (C) = 1.2;

	Cell AB ( tag::segment, A .reverse(), B );
	Cell BC ( tag::segment, B .reverse(), C );
	Cell CA ( tag::segment, C .reverse(), A );
	Cell ABC ( tag::triangle, AB, BC, CA );

	{ // just a block of code for hiding names
	fe_gauss .dock_on ( ABC );
	Function psi_A = fe_gauss .basis_function (A),
	         psi_B = fe_gauss .basis_function (B),
	         psi_C = fe_gauss .basis_function (C);
	std::cout << "        ";
	std::cout << fe_gauss .integrate ( psi_A .deriv(x) * psi_A ) << " "
						<< fe_gauss .integrate ( psi_A .deriv(x) * psi_A .deriv(x) +
																		 psi_A .deriv(y) * psi_A .deriv(y) ) << " "
						<< fe_gauss .integrate ( psi_A * psi_A ) << " "
						<< fe_gauss .integrate ( psi_A ) << " "
						<< fe_gauss .integrate ( psi_A .deriv(x) ) << " || ";
	std::cout << fe_gauss .integrate ( psi_A .deriv(x) * psi_B ) << " "
						<< fe_gauss .integrate ( psi_A .deriv(x) * psi_B .deriv(x) +
																		 psi_A .deriv(y) * psi_B .deriv(y) ) << " "
						<< fe_gauss .integrate ( psi_A * psi_B ) << " "
						<< fe_gauss .integrate ( psi_B ) << " "
						<< fe_gauss .integrate ( psi_B .deriv(x) ) << " || ";
	std::cout << fe_gauss .integrate ( psi_A .deriv(x) * psi_C ) << " "
						<< fe_gauss .integrate ( psi_A .deriv(x) * psi_C .deriv(x) +
																		 psi_A .deriv(y) * psi_C .deriv(y) ) << " "
						<< fe_gauss .integrate ( psi_A * psi_C ) << " "
						<< fe_gauss .integrate ( psi_C ) << " "
						<< fe_gauss .integrate ( psi_C .deriv(x) ) << std::endl;
	} { // just a block of code for hiding names
		fe_ufl_ffc .dock_on ( ABC, tag::winding );
	Function psi_A = fe_ufl_ffc .basis_function (A),
	         psi_B = fe_ufl_ffc .basis_function (B),
	         psi_C = fe_ufl_ffc .basis_function (C);
	std::vector < double > result =
		fe_ufl_ffc .integrate ( tag::pre_computed, tag::replace, bf1, tag::by, psi_A,
	                                             tag::replace, bf2, tag::by, psi_A );
	assert ( result .size() == 6 );
	std::cout << result [0] << " " << result [1] << " " << result [2] << " " << result [3] << " " << result [4] << " || ";
	result = fe_ufl_ffc .integrate
		( tag::pre_computed, tag::replace, bf1, tag::by, psi_A,
		                     tag::replace, bf2, tag::by, psi_B );
	assert ( result .size() == 6 );
	std::cout << result [0] << " " << result [1] << " " << result [2] << " " << result [3] << " " << result [4] << " || ";
	result = fe_ufl_ffc .integrate
		( tag::pre_computed, tag::replace, bf1, tag::by, psi_A,
		                     tag::replace, bf2, tag::by, psi_C );
	assert ( result .size() == 6 );
	std::cout << result [0] << " " << result [1] << " " << result [2] << " " << result [3] << " " << result [4] << std::endl;
	} // just a block of code for hiding names

	std::cout << "A=2=R B=0=P C=1=Q" << std::endl;
	
	return 0;
}
