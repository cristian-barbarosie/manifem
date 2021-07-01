
// example presented in paragraph 2.4 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// builds a circle by joining four segments on an implicit manifold

#include "maniFEM.h"

using namespace maniFEM;
using namespace std;

int main () {

	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0],  y = xy[1];

	RR2.implicit ( x*x + y*y == 1. );

	Cell N ( tag::vertex );  x(N) =  0.;  y(N) =  1.;
	Cell W ( tag::vertex );  x(W) = -1.;  y(W) =  0.;
	Cell S ( tag::vertex );  x(S) =  0.;  y(S) = -1.;
	Cell E ( tag::vertex );  x(E) =  1.;  y(E) =  0.;

	Mesh NW ( tag::segment, N.reverse(), W, tag::divided_in, 5 );
	Mesh WS ( tag::segment, W.reverse(), S, tag::divided_in, 5 );
	Mesh SE ( tag::segment, S.reverse(), E, tag::divided_in, 5 );
	Mesh EN ( tag::segment, E.reverse(), N, tag::divided_in, 5 );
	Mesh circle ( tag::join, NW, WS, SE, EN );

	circle.draw_ps ("circle.eps");
	circle.export_msh ("circle.msh");
	
	cout << "produced files circle.eps and circle.msh" << endl;
}
