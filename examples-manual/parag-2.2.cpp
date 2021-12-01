
// example presented in paragraph 2.2 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// builds an arrow-shaped mesh

#include "maniFEM.h"

using namespace maniFEM;
using namespace std;

int main () {

	// we choose our (geometric) space dimension :
	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	
	// xy is a map defined on our future mesh with values in RR2 :
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );

	// we can extract components of xy using the [] operator :
	Function x = xy [0], y = xy[1];

	Cell A ( tag::vertex );  x (A) = 0.;  y (A) =  0.  ;
	Cell B ( tag::vertex );  x (B) = 1.;  y (B) = -0.5 ;
	Cell C ( tag::vertex );  x (C) = 1.;  y (C) = -0.25;
	Cell D ( tag::vertex );  x (D) = 3.;  y (D) = -0.25;
	Cell E ( tag::vertex );  x (E) = 3.;  y (E) =  0.25;
	Cell F ( tag::vertex );  x (F) = 1.;  y (F) =  0.25;
	Cell G ( tag::vertex );  x (G) = 1.;  y (G) =  0.5 ;
	Mesh AB ( tag::segment, A .reverse(), B, tag::divided_in, 15 );
	Mesh BC ( tag::segment, B .reverse(), C, tag::divided_in,  4 );
	Mesh CF ( tag::segment, C .reverse(), F, tag::divided_in,  7 );
	Mesh CD ( tag::segment, C .reverse(), D, tag::divided_in, 25 );
	Mesh DE ( tag::segment, D .reverse(), E, tag::divided_in,  7 );
	Mesh EF ( tag::segment, E .reverse(), F, tag::divided_in, 25 );
	Mesh FG ( tag::segment, F .reverse(), G, tag::divided_in,  4 );
	Mesh GA ( tag::segment, G .reverse(), A, tag::divided_in, 15 );

	Mesh BG ( tag::join, BC, CF, FG );
	
	Mesh ABG ( tag::triangle, AB, BG, GA );
	Mesh CDEF ( tag::rectangle, CD, DE, EF, CF .reverse() );
	Mesh arrow ( tag::join, ABG, CDEF );

	arrow .draw_ps ( "arrow.eps");
	arrow .export_msh ("arrow.msh");
	
	cout << "produced files arrow.eps and arrow.msh" << endl;
}
