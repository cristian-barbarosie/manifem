
// example presented in paragraph 2.1 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// builds an L-shaped mesh by joining two segments and then two rectangles

#include "maniFEM.h"

using namespace maniFEM;
using namespace std;

int main () {

	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0],  y = xy[1];

	Cell A ( tag::vertex );  x(A) = -1.;  y(A) =  0.;
	Cell C ( tag::vertex );  x(C) =  0.;  y(C) =  0.5;
	Cell D ( tag::vertex );  x(D) = -1.;  y(D) =  0.5;
	Cell E ( tag::vertex );  x(E) =  0.;  y(E) =  1.;
	Cell F ( tag::vertex );  x(F) = -1.;  y(F) =  1.;
	Cell G ( tag::vertex );  x(G) =  1.;  y(G) =  0.;
	Cell H ( tag::vertex );  x(H) =  1.;  y(H) =  0.5;
	Mesh AG ( tag::segment, A.reverse(), G, tag::divided_in, 22 );
	Mesh GH ( tag::segment, G.reverse(), H, tag::divided_in,  8 );
	Mesh HC ( tag::segment, H.reverse(), C, tag::divided_in, 12 );
	Mesh CD ( tag::segment, C.reverse(), D, tag::divided_in, 10 );
	Mesh HD ( tag::join, HC, CD );
	Mesh DA ( tag::segment, D.reverse(), A, tag::divided_in,  8 );
	Mesh CE ( tag::segment, C.reverse(), E, tag::divided_in,  7 );
	Mesh EF ( tag::segment, E.reverse(), F, tag::divided_in, 10 );
	Mesh FD ( tag::segment, F.reverse(), D, tag::divided_in,  7 );
	Mesh AGHD ( tag::rectangle, AG, GH, HD, DA );
	Mesh CEFD ( tag::rectangle, CE, EF, FD, CD.reverse() );
	Mesh L_shaped ( tag::join, AGHD, CEFD );

	L_shaped.draw_ps ("L-shaped.eps");
	L_shaped.export_msh ("L-shaped.msh");

	cout << "produced files L-shaped.eps and L-shaped.msh" << endl;
}
