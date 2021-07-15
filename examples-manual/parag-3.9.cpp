
// example presented in paragraph 3.9 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// a "bumpy" hemisphere meshed progressively

#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;
using namespace std;


int main ()

{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz[0],  y = xyz[1],  z = xyz[2];

	Manifold nut = RR3.implicit ( x*x + y*y + z*z + 1.5*x*y*z == 1. );

	// build the base (a closed curve)
	nut.implicit ( x*x + 3.*z == 0. );

	Mesh circle ( tag::progressive, tag::desired_length, 0.1, tag::random_orientation );

	nut.set_as_working_manifold();
	Mesh bumpy ( tag::progressive, tag::boundary, circle, tag::desired_length, 0.098 );
	bumpy.export_msh ("bumpy.msh");
	
	cout << "produced file bumpy.msh" << endl;
}
