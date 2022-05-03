
// example presented in paragraph 2.18 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// builds a spiral mesh (parametric)

#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;
using namespace std;

int main ()

{	Manifold spiral ( tag::Euclid, tag::of_dim, 1 );
	Function t = spiral .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	const double pi = 4. * std::atan(1.);
	
	Cell A ( tag::vertex );  t(A) = pi/2.;
	Cell B ( tag::vertex );  t(B) = 5.*pi;
	Mesh arc_of_spiral ( tag::segment, A .reverse(), B, tag::divided_in, 50 );

	// forget about t, in future statements x and y will be used
	Function x = t * cos(t), y = t * sin(t);
	Manifold RR2 ( tag::Euclid, tag::of_dimension, 2 );
	RR2 .set_coordinates ( x && y );

	arc_of_spiral .draw_ps ("spiral.eps");
	arc_of_spiral .export_to_file ( tag::msh, "spiral.msh");

	cout << "produced files spiral.eps and spiral.msh" << endl;

}  // end of main
