
// example presented in paragraph 2.11 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// builds a circle in 3d (implicit manifold with two equations)

#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;
using namespace std;

int main ()

{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz[0],  y = xyz[1],  z = xyz[2];

	RR3.implicit ( x*x + y*y == 1., x*y == 4.*z );
	Cell S ( tag::vertex );  x(S) =  0.;  y(S) = -1.;  z(S) = 0.;
	Cell E ( tag::vertex );  x(E) =  1.;  y(E) =  0.;  z(E) = 0.;
	Cell N ( tag::vertex );  x(N) =  0.;  y(N) =  1.;  z(N) = 0.;
	Cell W ( tag::vertex );  x(W) = -1.;  y(W) =  0.;  z(W) = 0.;
	Mesh SE ( tag::segment, S.reverse(), E, tag::divided_in, 5 );
	Mesh EN ( tag::segment, E.reverse(), N, tag::divided_in, 5 );
	Mesh NW ( tag::segment, N.reverse(), W, tag::divided_in, 5 );
	Mesh WS ( tag::segment, W.reverse(), S, tag::divided_in, 5 );
	Mesh circle ( tag::join, SE, EN, NW, WS );

	circle.draw_ps_3d ("circle-3d.eps");
	circle.export_msh ("circle-3d.msh");

	cout << "produced files circle-3d.eps and circle-3d.msh" << endl;
}
