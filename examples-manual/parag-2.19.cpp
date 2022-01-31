
// example presented in paragraph 2.19 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// builds an imperfect donut (parametric, starts with high dimension)

#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;
using namespace std;

int main ()

{	Manifold RR5 ( tag::Euclid, tag::of_dim, 5 );
	Function xyzab = RR5 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyzab [0], y = xyzab [1], z = xyzab [2], alpha = xyzab [3], beta = xyzab [4];

	const double big_radius = 3, small_radius = 1;
	Manifold torus = RR5 .parametric
		( x == ( big_radius + small_radius * cos(beta) ) * cos(alpha),
	    y == ( big_radius + small_radius * cos(beta) ) * sin(alpha),
	    z == small_radius * sin(beta)                               );
	
	const double pi = 4. * std::atan(1.);
	Cell A ( tag::vertex );  alpha (A) = 0.;       beta (A) = 0.;      torus .project (A);
	Cell B ( tag::vertex );  alpha (B) = 0.;       beta (B) = 1.9*pi;  torus .project (B);
	Cell C ( tag::vertex );  alpha (C) = 1.95*pi;  beta (C) = 1.9*pi;  torus .project (C);
	Cell D ( tag::vertex );  alpha (D) = 1.95*pi;  beta (D) = 0.;      torus .project (D);

	Mesh AB ( tag::segment, A .reverse(), B, tag::divided_in, 19 );
	Mesh BC ( tag::segment, B .reverse(), C, tag::divided_in, 39 );
	Mesh CD ( tag::segment, C .reverse(), D, tag::divided_in, 19 );
	Mesh DA ( tag::segment, D .reverse(), A, tag::divided_in, 39 );

	Mesh ABCD ( tag::rectangle, AB, BC, CD, DA );

	// forget about alpha and beta
	Manifold RR3 ( tag::Euclid, tag::of_dimension, 3 );
	RR3 .set_coordinates ( x && y && z );

	ABCD .export_to_file ( tag::msh, "torus.msh");

	cout << "produced file torus.msh" << endl;

}  // end of main
