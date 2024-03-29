
// example presented in paragraph 3.7 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// a complicated surface, a kind of convolution between two tori

#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;
using namespace std;

int main ()
	
{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz [0], y = xyz [1], z = xyz [2];

	cout << "this example takes time" << endl;
	Function f1 = max ( x*x + y*y, 0.1 );  // we cut with 0.1 to avoid singularities
	Function f2 = 1. - power ( f1, -0.5 );
	Function d1 = f1 * f2 * f2 + z*z;
	Function f3 = max ( (x-0.4)*(x-0.4) + z*z, 0.1 );
	Function f4 = 1. - power ( f3, -0.5 );
	Function d2 = y*y + f3 * f4 * f4;
	RR3 .implicit ( smooth_min ( d1, d2, tag::threshold, 0.2 ) == 0.15 );

	Mesh two_tori ( tag::frontal, tag::desired_length, 0.09 );
	// tag::desired_length, 0.03 + 0.05 * abs ( d2-d1 )

	two_tori .export_to_file ( tag::msh, "two-tori.msh");
	std::cout << "produced file two-tori.msh" << std::endl;

}  // end of main


void main_1 ()
	
{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz [0], y = xyz [1], z = xyz [2];

	cout << "this example takes time" << endl;
	Function f1 = x*x + y*y + 0.1;
	Function f2 = 1. - power ( f1, -0.5 );
	Function d1 = f1 * f2 * f2 + z*z;
	Function f3 = x*x + z*z + 0.1;
	Function f4 = 1. - power ( f3, -0.5 );
	Function d2 = y*y + f3 * f4 * f4;
	RR3 .implicit ( smooth_min ( d1, d2, tag::threshold, 0.1 ) == 0.15 );

	Mesh two_tori ( tag::frontal, tag::desired_length, 0.09 );

	two_tori .export_to_file ( tag::msh, "two-tori.msh");
	std::cout << "produced file two-tori.msh" << std::endl;

}  // end of main_1

