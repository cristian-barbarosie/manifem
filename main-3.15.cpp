
// example presented in paragraph 3.15 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// a "bumpy" hemisphere built progressively

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

	Cell S ( tag::vertex );  x(S) =  0.;  y(S) = -1.;  z(S) =  0.;
	std::vector < double > tau { 1., 0., 0. };
	Mesh circle ( tag::progressive, tag::start_at, S, tag::towards, tau,
                tag::desired_length, 0.1                               );

	nut.set_as_working_manifold();
	std::vector < double > N { 0., 0., -1. };
	Mesh bumpy ( tag::progressive, tag::boundary, circle,
               tag::start_at, S, tag::towards, N,
               tag::desired_length, 0.1                 );
	bumpy.export_msh ("bumpy.msh");
	
	cout << "produced file bumpy.msh" << endl;
}
