
// example presented in paragraph 2.9 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// a diamond-shape domain in RR2 (alternating between manifolds)

#include "maniFEM.h"

using namespace maniFEM;
using namespace std;

int main () {

	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0],  y = xy[1];

	Cell N ( tag::vertex );  x(N) =  0.;  y(N) =  1.;
	Cell W ( tag::vertex );  x(W) = -1.;  y(W) =  0.;
	Cell S ( tag::vertex );  x(S) =  0.;  y(S) = -1.;
	Cell E ( tag::vertex );  x(E) =  1.;  y(E) =  0.;

	RR2.implicit ( x*y + x - y == -1. );
	Mesh NW ( tag::segment, N.reverse(), W, tag::divided_in, 10 );
	RR2.implicit ( x*y - x - y ==  1. );
	Mesh WS ( tag::segment, W.reverse(), S, tag::divided_in, 10 );
	RR2.implicit ( x*y - x + y == -1. );
	Mesh SE ( tag::segment, S.reverse(), E, tag::divided_in, 10 );
	RR2.implicit ( x*y + x + y ==  1. );
	Mesh EN ( tag::segment, E.reverse(), N, tag::divided_in, 10 );

	RR2.set_as_working_manifold();

	Mesh diamond ( tag::rectangle, NW, WS, SE, EN );

	diamond.draw_ps ("diamond.eps");
	diamond.export_msh ("diamond.msh");
	
	cout << "produced files diamond.eps and diamond.msh" << endl;
}
