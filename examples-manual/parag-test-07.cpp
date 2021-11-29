



#include "maniFEM.h"
#include "math.h"
using namespace maniFEM;


int main ( )

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0], y = xy[1];

	Cell A ( tag::vertex);  x (A) = 0.;  y (A) = 0.;
	Cell B ( tag::vertex);  x (B) = 1.;  y (B) = 0.;
	Cell C ( tag::vertex);  x (C) = 0.5;  y (C) = std::sqrt ( 0.75 );

	Cell AB ( tag::segment, A .reverse(), B );  AB .core->name = "AB";
	Cell BC ( tag::segment, B .reverse(), C );  BC .core->name = "BC";
	Cell CA ( tag::segment, C .reverse(), A );  CA .core->name = "CA";

	Cell ABC ( tag::triangle, AB, BC, CA );

	FiniteElement fe ( tag::with_master, tag::triangle, tag::Lagrange, tag::of_degree, 1 );
	fe .set_integrator ( tag::Gauss, tag::tri_4 );

	fe .dock_on ( ABC );
	Function psi = fe .basis_function (B);
	std::cout << fe .integrate ( psi .deriv (x) ) << std::endl;
	std::cout << fe .integrate ( psi .deriv (y) ) << std::endl;
	
}
