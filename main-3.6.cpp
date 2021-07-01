
// example presented in paragraph 3.6 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// meshes a sphere progressively

#include "maniFEM.h"
#include <iomanip>
#include "math.h"

using namespace maniFEM;
using namespace std;

int main ()

{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function  x = xyz[0],  y = xyz[1],  z = xyz[2];

	cout << "this example takes some time" << endl;
	// cout << setprecision(10);

	RR3.implicit ( x*x + y*y + z*z == 1. );
	Mesh sphere ( tag::progressive, tag::desired_length, 0.06 );

	sphere.export_msh ("sphere.msh");

	cout << "produced file sphere.msh" << endl;
}
