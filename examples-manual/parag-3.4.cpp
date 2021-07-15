
// example presented in paragraph 3.4 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// builds a circle in 3D (implicit manifold with two equations) progressively

#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;
using namespace std;

int main ()

{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz[0],  y = xyz[1],  z = xyz[2];

	RR3.implicit ( x*x + y*y == 1., x*y == 4.*z );
	Mesh circle ( tag::progressive, tag::desired_length, 0.1, tag::random_orientation );

	circle.draw_ps_3d ("circle-3d.eps");
	circle.export_msh ("circle-3d.msh");

	cout << "produced files circle-3d.eps and circle-3d.msh" << endl;
}
