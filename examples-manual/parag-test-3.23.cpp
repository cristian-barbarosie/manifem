
// example presented in paragraph 3.23 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// two intersecting tori

#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;
using namespace std;

int main ()
	
{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz [0], y = xyz [1], z = xyz [2];

	// here, we don't need to add 0.1 to f1 and f3,
	//       the shape of each torus avoids the singularities
	// we change slightly f3 (replace 0.4 by 0.5) to highlight the intersection
	Function f1 = x*x + y*y;
	Function f2 = 1. - power ( f1, -0.5 );
	Function d1 = f1 * f2 * f2 + z*z;  // squared distance to a circle in the xy plane
	Function f3 = (x-0.5)*(x-0.5) + z*z;
	Function f4 = 1. - power ( f3, -0.5 );
	Function d2 = y*y + f3 * f4 * f4;  // squared distance to a circle in the xz plane

	const double seg_size = 0.1;
	// Function seg_size = 0.03 + 0.05 * abs ( d2-d1 );
	// Function seg_size = 0.1 * power ( 0.15 + abs ( d2-d1 ), 0.4 );

	Manifold intersection = RR3 .implicit ( d1 == 0.15, d2 == 0.15 );

	Cell start_1 ( tag::vertex );
	x ( start_1 ) = 1.3;  y ( start_1 ) = -0.4;  z ( start_1 ) = 0.;
	intersection .project ( start_1 );
	Mesh circle_1 ( tag::frontal, tag::start_at, start_1, tag::towards, { 0., 0., 1. },
	                tag::desired_length, seg_size                                      );

	Cell start_2 ( tag::vertex );
	x ( start_2 ) = -0.9;  y ( start_2 ) = 0.;  z ( start_2 ) = 0.4;
	intersection .project ( start_2 );
	Mesh circle_2 ( tag::frontal, tag::start_at, start_2, tag::towards, { 0., 1., 0. },
	                tag::desired_length, seg_size                                      );

	Mesh two_circles ( tag::join, circle_1, circle_2 );

	RR3 .implicit ( d1 == 0.15 );
	Mesh torus_1 ( tag::frontal, tag::boundary, two_circles,
	               tag::start_at, start_1, tag::towards, { -0.2, -1., 0. },
	               tag::desired_length, seg_size                           );

	RR3 .implicit ( d2 == 0.15 );
  Mesh torus_2 ( tag::frontal, tag::boundary, two_circles .reverse(),
	               tag::start_at, start_2, tag::towards, { 0.2, 0., 1. },
	               tag::desired_length, seg_size                           );

	Mesh two_tori ( tag::join, torus_1, torus_2 );

	two_tori .export_to_file ( tag::msh, "two-tori.msh");
	std::cout << "produced file two-tori.msh" << std::endl;

}  // end of main



