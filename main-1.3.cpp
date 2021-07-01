
// example presented in paragraph 1.3 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// builds an L-shaped mesh by joining rectangles

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

	Cell A ( tag::vertex );  x(A) = -1.;  y(A) = 0.;
	Cell B ( tag::vertex );  x(B) =  0.;  y(B) = 0.;
	Cell C ( tag::vertex );  x(C) =  0.;  y(C) = 0.5;
	Cell D ( tag::vertex );  x(D) = -1.;  y(D) = 0.5;
	Cell E ( tag::vertex );  x(E) =  0.;  y(E) = 1.;
	Cell F ( tag::vertex );  x(F) = -1.;  y(F) = 1.;
	Cell G ( tag::vertex );  x(G) =  1.;  y(G) = 0.;
	Cell H ( tag::vertex );  x(H) =  1.;  y(H) = 0.5;
	Mesh AB ( tag::segment, A.reverse(), B, tag::divided_in, 10 );
	Mesh BC ( tag::segment, B.reverse(), C, tag::divided_in, 8 );
	Mesh CD ( tag::segment, C.reverse(), D, tag::divided_in, 10 );
	Mesh DA ( tag::segment, D.reverse(), A, tag::divided_in, 8 );
	Mesh CE ( tag::segment, C.reverse(), E, tag::divided_in, 7 );
	Mesh EF ( tag::segment, E.reverse(), F, tag::divided_in, 10 );
	Mesh FD ( tag::segment, F.reverse(), D, tag::divided_in, 7 );
	Mesh BG ( tag::segment, B.reverse(), G, tag::divided_in, 12 );
	Mesh GH ( tag::segment, G.reverse(), H, tag::divided_in, 8 );
	Mesh HC ( tag::segment, H.reverse(), C, tag::divided_in, 12 );

	Mesh ABCD ( tag::rectangle, AB, BC, CD, DA );
	Mesh CEFD ( tag::rectangle, CE, EF, FD, CD.reverse() );
	Mesh BGHC ( tag::rectangle, GH, HC, BC.reverse(), BG );
	Mesh L_shaped ( tag::join, ABCD, CEFD, BGHC );

	L_shaped.export_msh ("L-shaped.msh");
	L_shaped.draw_ps ( "L-shaped.eps");
	
	cout << "produced files L-shaped.msh and L-shaped.eps" << endl;
}
