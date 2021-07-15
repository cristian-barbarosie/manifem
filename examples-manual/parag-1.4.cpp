
// example presented in paragraph 1.4 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// builds a trapezoidal mesh made of triangles

#include "maniFEM.h"

using namespace maniFEM;
using namespace std;

int main () {

	// we choose our (geometric) space dimension :
	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	
	// xy is a map defined on our future mesh with values in RR2 :
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );

	// we can extract components of xyz using the [] operator :
	Function x = xy[0],  y = xy[1];

	Cell A ( tag::vertex );  x(A) = -1. ;  y(A) = 0.;
	Cell B ( tag::vertex );  x(B) =  0. ;  y(B) = 0.;
	Cell C ( tag::vertex );  x(C) =  1. ;  y(C) = 0.;
	Cell D ( tag::vertex );  x(D) = -0.5;  y(D) = 1.;
	Cell E ( tag::vertex );  x(E) =  0.5;  y(E) = 1.;
	Mesh AB ( tag::segment, A.reverse(), B, tag::divided_in, 8 );
	Mesh BC ( tag::segment, B.reverse(), C, tag::divided_in, 8 );
	Mesh AD ( tag::segment, A.reverse(), D, tag::divided_in, 8 );
	Mesh BD ( tag::segment, B.reverse(), D, tag::divided_in, 8 );
	Mesh BE ( tag::segment, B.reverse(), E, tag::divided_in, 8 );
	Mesh CE ( tag::segment, C.reverse(), E, tag::divided_in, 8 );
	Mesh ED ( tag::segment, E.reverse(), D, tag::divided_in, 8 );
	Mesh ABD ( tag::triangle, AB, BD, AD.reverse() );
	Mesh BCE ( tag::triangle, BC, CE, BE.reverse() );
	Mesh BED ( tag::triangle, BE, ED, BD.reverse() );
	Mesh three_tri ( tag::join, ABD, BCE, BED );

	three_tri.export_msh ("three-tri.msh");
	
	cout << "produced file three-tri.msh" << endl;
}
