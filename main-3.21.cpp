
// example presented in paragraph 3.21 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// a singular point at the tangency point between two manifolds
// the code shown in the manual does not work (yet)
// we fake the result by building by hand the triangles around the singular point


#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;
using namespace std;


int main ()

{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz[0],  y = xyz[1],  z = xyz[2];
	double seg_size = 0.1;
	// std::cout << "segment size : ";
	// std::cin >> seg_size;

	Manifold cyl_manif = RR3.implicit ( y*y + (z-0.5)*(z-0.5) == 0.25 );

	Cell V ( tag::vertex );  x(V) = 0.;  y(V) = seg_size;  z(V) = 1.;
	cyl_manif.project(V);
	Cell W ( tag::vertex );  x(W) = 0.;  y(W) = -seg_size;  z(W) = 1.;
	cyl_manif.project(W);
	
	Manifold intersection = cyl_manif.implicit ( x*x + y*y + z*z == 1. );
	Cell O ( tag::vertex );  x(O) = 0.;  y(O) = 0.;  z(O) = 1.;
	Cell A ( tag::vertex );  x(A) = 0.7*seg_size;  y(A) = 0.7*seg_size;  z(A) = 1.;
	intersection.project(A);
	Cell B ( tag::vertex );  x(B) = -0.7*seg_size;  y(B) = 0.7*seg_size;  z(B) = 1.;
	intersection.project(B);
	Cell C ( tag::vertex );  x(C) = -0.7*seg_size;  y(C) = -0.7*seg_size;  z(C) = 1.;
	intersection.project(C);
	Cell D ( tag::vertex );  x(D) = 0.7*seg_size;  y(D) = -0.7*seg_size;  z(D) = 1.;
	intersection.project(D);

	Cell OA ( tag::segment, O.reverse(), A );
	Cell OB ( tag::segment, O.reverse(), B );
	Cell OV ( tag::segment, O.reverse(), V );
	Cell AV ( tag::segment, A.reverse(), V );
	Cell VB ( tag::segment, V.reverse(), B );
	Cell OAV ( tag::triangle, OA, AV, OV.reverse() );
	Cell OVB ( tag::triangle, OV, VB, OB.reverse() );
	Cell OC ( tag::segment, O.reverse(), C );
	Cell OD ( tag::segment, O.reverse(), D );
	Cell OW ( tag::segment, O.reverse(), W );
	Cell CW ( tag::segment, C.reverse(), W );
	Cell WD ( tag::segment, W.reverse(), D );
	Cell OCW ( tag::triangle, OC, CW, OW.reverse() );
	Cell OWD ( tag::triangle, OW, WD, OD.reverse() );
	
	std::vector < double > tau { 1., 1., 0. };
	Mesh circle_1 ( tag::progressive, tag::start_at, A, tag::towards, tau,
	                tag::stop_at, D, tag::desired_length, seg_size         );
	tau = { -1., 1., 0. };
	Mesh circle_2 ( tag::progressive, tag::start_at, B, tag::towards, tau,
	                tag::stop_at, C, tag::desired_length, seg_size         );

	Mesh circles ( tag::join, circle_1, circle_2.reverse() );
	AV.reverse().add_to ( circles );  VB.reverse().add_to ( circles );
	CW.reverse().add_to ( circles );  WD.reverse().add_to ( circles );

	cyl_manif.set_as_working_manifold();
	tau = { 0., 1., 0. };
	Mesh cone( tag::progressive, tag::boundary, circles,
	           tag::start_at, V, tag::towards, tau,
	           tag::desired_length, seg_size );

	OAV.add_to ( cone );  OVB.add_to ( cone );
	OCW.add_to ( cone );  OWD.add_to ( cone );

	cone.export_msh ("napkin.msh");
	cout << "produced file napkin.msh" << endl;
}
