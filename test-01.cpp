

#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;

int main ()

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0],  y = xy[1];

	// we introduce two equivalence relations on RR2
	Function::Action g1 ( tag::transforms, xy, tag::into, (x+3.) && (y+0.1) );
	Function::Action g2 ( tag::transforms, xy, tag::into, (x-0.1) && (y+1.) );

	// and divide RR2 by these equivalence relations
	Manifold torus_manif = RR2.quotient ( g1, g2 );

	Cell V ( tag::vertex );  x(V) = 1.1;  y(V) = 22.2;

	Function xym = torus_manif.coordinates();
	// Function xm = xym[0], ym = xym[1];
	std::vector < double > xyV = xym(V);
	std::cout << "g1(3.,0.1) g2(-0.1,1.)" << std::endl;
	xyV = xym ( V, tag::spin, {0,0} );
	std::cout << " 0 0 : " << xyV[0] << " " << xyV[1] << std::endl;
	xyV = xym ( V, tag::spin, {0,1} );
	std::cout << " 0 1 : " << xyV[0] << " " << xyV[1] << std::endl;
	xyV = xym ( V, tag::spin, {1,0} );
	std::cout << " 1 0 : " << xyV[0] << " " << xyV[1] << std::endl;
	xyV = xym ( V, tag::spin, {5,0} );
	std::cout << " 5 0 : " << xyV[0] << " " << xyV[1] << std::endl;
	xyV = xym ( V, tag::spin, {-1,0} );
	std::cout << "-1 0 : " << xyV[0] << " " << xyV[1] << std::endl;
	xyV = xym ( V, tag::spin, {-3,7} );
	std::cout << "-3 7 : " << xyV[0] << " " << xyV[1] << std::endl;

	// std::vector < double > vv =
	// 	Function::Vector::MultiValued::JumpIsSum::analyse_linear_expression
	// 	( (4.1+x) && (y-2), xy );
	// for ( size_t i = 0; i < vv.size(); i++ )
	// 	std::cout << vv[i] << " ";
	// std::cout << std::endl;
		
	// create a multi-valued function test (for testing purposes)
	// Function test ( tag::Lagrange, tag::of_degree, 1, tag::multivalued );
	// describe the action of g1 and g2 on test :
	// g1 * test = test + 1.5;
	// g2 * test = test - 0.2;
	// test.property ( tag::through, g1, tag::becomes, test + 1.5 );
	// test.property ( tag::through, g2, tag::becomes, test - 0.2 );

	// std::cout << test(V) << std::endl << std::flush;
	// std::cout << "0 0 " << std::endl << std::flush;
	// std::cout << test ( V, tag::spin, {0,0} ) << std::endl << std::flush;
	// std::cout << "1 0 " << std::endl << std::flush;
	// std::cout << test ( V, tag::spin, {1,0} ) << std::endl << std::flush;
	// std::cout << "1 1 " << std::endl << std::flush;
	// std::cout << test ( V, tag::spin, {1,1} ) << std::endl << std::flush;
	// std::cout << "1 2 " << std::endl << std::flush;
	// std::cout << test ( V, tag::spin, {1,2} ) << std::endl << std::flush;
	
}
