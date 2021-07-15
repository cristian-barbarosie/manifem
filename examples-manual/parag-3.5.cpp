
// example presented in paragraph 3.5 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// builds a spiral in RR2

#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;
using namespace std;

int main ()

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0],  y = xy[1];

	cout << "this example takes some time" << endl;

	Function r = power ( x*x + y*y, 0.25 );
	RR2.implicit ( x*sin(r) == y*cos(r) );

	Cell A ( tag::vertex );  x(A) = std::pow(  3.14159,2.) ;  y(A) =  0.  ;
	Cell B ( tag::vertex );  x(B) = std::pow(9*3.14159,2.) ;  y(B) =  0.  ;
	Mesh spiral ( tag::progressive, tag::start_at, A, tag::stop_at, B, tag::desired_length, 1. );

	spiral.export_msh ("spiral.msh");

	cout << "produced file spiral.msh" << endl;
}
