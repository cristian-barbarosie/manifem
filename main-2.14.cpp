
// example presented in paragraph 2.14 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// closes a circle in a cumbersome manner

#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;
using namespace std;

int main ()

{	Manifold circle_manif ( tag::Euclid, tag::of_dim, 1 );
	Function t = circle_manif.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	const double pi = 4.*atan(1.);
	
	Cell A ( tag::vertex );  t(A) = 0.;
	Cell B ( tag::vertex );  t(B) = 1.9*pi;
	Mesh circle ( tag::segment, A.reverse(), B, tag::divided_in, 19 );

	Cell BA ( tag::segment, B.reverse(), A );
	BA.add_to ( circle );

	Function x = cos(t), y = sin(t);
	circle_manif.set_coordinates ( x && y );

	circle.draw_ps ("circle.eps");
	circle.export_msh ("circle.msh");

	cout << "produced files circle.eps and circle.msh" << endl;
}
