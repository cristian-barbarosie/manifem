
// example presented in paragraph 1.5 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// builds a mesh which mixes triangles and rectangles

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

	Cell A ( tag::vertex );  x(A) = -1. ;  y(A) = 0. ;
	Cell B ( tag::vertex );  x(B) =  0. ;  y(B) = 0. ;
	Cell C ( tag::vertex );  x(C) =  1. ;  y(C) = 0. ;
	Cell D ( tag::vertex );  x(D) = -0.5;  y(D) = 0.8;
	Cell E ( tag::vertex );  x(E) =  0.5;  y(E) = 0.8;
	Cell F ( tag::vertex );  x(F) =  0. ;  y(F) = 1. ;
	Mesh AB ( tag::segment, A.reverse(), B, tag::divided_in, 8 );
	Mesh BC ( tag::segment, B.reverse(), C, tag::divided_in, 8 );
	Mesh AD ( tag::segment, A.reverse(), D, tag::divided_in, 8 );
	Mesh BD ( tag::segment, B.reverse(), D, tag::divided_in, 8 );
	Mesh BE ( tag::segment, B.reverse(), E, tag::divided_in, 8 );
	Mesh CE ( tag::segment, C.reverse(), E, tag::divided_in, 8 );
	Mesh EF ( tag::segment, E.reverse(), F, tag::divided_in, 8 );
	Mesh FD ( tag::segment, F.reverse(), D, tag::divided_in, 8 );
	Mesh ABD ( tag::triangle, AB, BD, AD.reverse() );
	Mesh BCE ( tag::triangle, BC, CE, BE.reverse() );
	Mesh BEFD ( tag::rectangle, BE, EF, FD, BD.reverse() );
	Mesh two_tri_one_rect ( tag::join, ABD, BEFD, BCE );

	two_tri_one_rect.export_msh ("two-tri-one-rect.msh");
	
	cout << "produced file two-tri-one-rect.msh" << endl;
}
