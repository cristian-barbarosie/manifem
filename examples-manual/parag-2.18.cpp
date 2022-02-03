
// example presented in paragraph 2.18 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// closes a circle in a cumbersome manner

#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;
using namespace std;

int main ()

{	Manifold circle_manif ( tag::Euclid, tag::of_dim, 1 );
	Function t = circle_manif .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	const double pi = 4. * std::atan(1.);
	
	Cell A ( tag::vertex );  t (A) = 0.;
	Cell B ( tag::vertex );  t (B) = 1.9*pi;
	Mesh incomplete_circle ( tag::segment, A .reverse(), B, tag::divided_in, 19 );

	Mesh small_piece ( tag::segment, B .reverse(), A, tag::divided_in, 1 );
	Mesh circle ( tag::join, incomplete_circle, small_piece );

	// forget about t, in future statements x and y will be used
	Function x = cos(t), y = sin(t);
	Manifold RR2 ( tag::Euclid, tag::of_dimension, 2 );
	RR2 .set_coordinates ( x && y );

	circle .draw_ps ("circle.eps");
	circle .export_to_file ( tag::msh, "circle.msh");

	cout << "produced files circle.eps and circle.msh" << endl;

}  // end of main
