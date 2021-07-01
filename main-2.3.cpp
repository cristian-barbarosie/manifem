
// example presented in paragraph 2.3 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// an arc of hiperbola on an implicit manifold

#include "maniFEM.h"

using namespace maniFEM;
using namespace std;

int main () {

	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0],  y = xy[1];

	RR2.implicit ( x*y == 1. );

	Cell A ( tag::vertex );  x(A) =  0.5;   y(A) =  2.;
	Cell B ( tag::vertex );  x(B) =  3.;    y(B) =  0.333333333333;

	Mesh arc_of_hiperbola ( tag::segment, A.reverse(), B, tag::divided_in, 7 );

	arc_of_hiperbola.draw_ps ("hiperbola.eps");
	arc_of_hiperbola.export_msh ("hiperbola.msh");
	
	cout << "produced files hiperbola.eps and hiperbola.msh" << endl;
}
