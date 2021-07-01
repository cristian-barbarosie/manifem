
// example presented in paragraph 3.20 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// a sharp cone
// the code shown in the manual does not work (yet)
// we fake the result by building by hand the triangles around the vertex

#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;
using namespace std;


int main ()

{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz[0],  y = xyz[1],  z = xyz[2];
	double seg_size = 0.13;
	// std::cout << "segment size : ";
	// std::cin >> seg_size;

	Manifold cone_manif = RR3.implicit ( x*x + y*y == z*z );
	
	Cell O ( tag::vertex );  x(O) = 0.;  y(O) = 0.;  z(O) = 0.;
	Cell A ( tag::vertex );  x(A) = 0.7*seg_size;  y(A) = 0.;  z(A) = 0.7*seg_size;
	cone_manif.project(A);
	Cell B ( tag::vertex );  x(B) = 0.2*seg_size;  y(B) = 0.6*seg_size;  z(B) = 0.7*seg_size;
	cone_manif.project(B);
	Cell C ( tag::vertex );  x(C) = -0.5*seg_size;  y(C) = 0.4*seg_size;  z(C) = 0.7*seg_size;
	cone_manif.project(C);
	Cell D ( tag::vertex );  x(D) = -0.5*seg_size;  y(D) = -0.4*seg_size;  z(D) = 0.7*seg_size;
	cone_manif.project(D);
	Cell E ( tag::vertex );  x(E) = 0.2*seg_size;  y(E) = -0.6*seg_size;  z(E) = 0.7*seg_size;
	cone_manif.project(E);
	Cell OA ( tag::segment, O.reverse(), A );
	Cell OB ( tag::segment, O.reverse(), B );
	Cell OC ( tag::segment, O.reverse(), C );
	Cell OD ( tag::segment, O.reverse(), D );
	Cell OE ( tag::segment, O.reverse(), E );
	Cell AB ( tag::segment, A.reverse(), B );
	Cell BC ( tag::segment, B.reverse(), C );
	Cell CD ( tag::segment, C.reverse(), D );
	Cell DE ( tag::segment, D.reverse(), E );
	Cell EA ( tag::segment, E.reverse(), A );
	Cell OAB ( tag::triangle, OA, AB, OB.reverse() );
	Cell OBC ( tag::triangle, OB, BC, OC.reverse() );
	Cell OCD ( tag::triangle, OC, CD, OD.reverse() );
	Cell ODE ( tag::triangle, OD, DE, OE.reverse() );
	Cell OEA ( tag::triangle, OE, EA, OA.reverse() );
	Mesh ponta ( tag::of_dimension, 2, tag::greater_than_one );
	OAB.add_to ( ponta );
	OBC.add_to ( ponta );
	OCD.add_to ( ponta );
	ODE.add_to ( ponta );
	OEA.add_to ( ponta );

	Mesh small_circle ( tag::of_dimension_one );
	AB.add_to ( small_circle );
	BC.add_to ( small_circle );
	CD.add_to ( small_circle );
	DE.add_to ( small_circle );
	EA.add_to ( small_circle );

	Cell V ( tag::vertex );  x(V) = 1.;  y(V) = 0.;  z(V) = 1.;
	cone_manif.implicit ( z == 1. );
	std::vector < double > tau = { 0., 1., 0. };
	Mesh big_circle ( tag::progressive, tag::start_at, V, tag::towards, tau,
	                  tag::desired_length, seg_size );
	Mesh two_circles ( tag::join, small_circle.reverse(), big_circle );

	cone_manif.set_as_working_manifold();
	tau = { -1., 0., -1. };
	Mesh cone_up ( tag::progressive, tag::boundary, two_circles,
	               tag::start_at, V, tag::towards, tau,
	               tag::desired_length, seg_size );

	Mesh cone ( tag::join, ponta, cone_up );
	cone.export_msh ("cone.msh");
	cout << "produced file cone.msh" << endl;
}
