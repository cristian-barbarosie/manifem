
// example presented in paragraph 2.15 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// builds an imperfect donut (parametric)

#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;
using namespace std;

int main ()

{	Manifold torus ( tag::Euclid, tag::of_dim, 2 );
	Function alpha_beta = torus.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function alpha = alpha_beta[0], beta = alpha_beta[1];

	const double big_radius = 3, small_radius = 1;
	Function x = ( big_radius + small_radius * cos(beta) ) * cos(alpha),
	         y = ( big_radius + small_radius * cos(beta) ) * sin(alpha),
	         z = small_radius * sin(beta);

	const double pi = 4.*atan(1.);
	Cell A ( tag::vertex );  alpha(A) = 0.;       beta(A) = 0.;
	Cell B ( tag::vertex );  alpha(B) = 0.;       beta(B) = 1.9*pi;
	Cell C ( tag::vertex );  alpha(C) = 1.95*pi;  beta(C) = 1.9*pi;
	Cell D ( tag::vertex );  alpha(D) = 1.95*pi;  beta(D) = 0.;

	Mesh AB ( tag::segment, A.reverse(), B, tag::divided_in, 19 );
	Mesh BC ( tag::segment, B.reverse(), C, tag::divided_in, 39 );
	Mesh CD ( tag::segment, C.reverse(), D, tag::divided_in, 19 );
	Mesh DA ( tag::segment, D.reverse(), A, tag::divided_in, 39 );

	torus.set_as_working_manifold();
	Mesh ABCD ( tag::rectangle, AB, BC, CD, DA );

	torus.set_coordinates ( x && y && z );

	ABCD.export_msh ("torus.msh");

	cout << "produced file torus.msh" << endl;
}
