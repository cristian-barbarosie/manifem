
// example presented in paragraph 7.1 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// mesh the one-dimensional quotient manifold RR/ZZ

#include "maniFEM.h"
using namespace maniFEM;


int main ( )

{	Manifold RR ( tag::Euclid, tag::of_dim, 1 );
	Function x = RR .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );

	// define an action on RR (a translation)
	Manifold::Action g ( tag::transforms, x, tag::into, x+1. );

	// and divide RR by this equivalence relation
	Manifold circle = RR .quotient ( g );

	// one vertex is enough to start the process
	Cell A ( tag::vertex );  x (A) = 0.02;

	// with this vertex, we build a segment
	Mesh seg ( tag::segment, A .reverse(), A, tag::divided_in, 10, tag::winding, g );

}  // end of main
