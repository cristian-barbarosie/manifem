
// example presented in paragraph 10.3 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// in a mesh of squares, cut some of the squares in two triangles

#include "maniFEM.h"
#include <fstream>
#include <iostream>

using namespace std;
using namespace maniFEM;


int main ( )

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy [0], y = xy [1];
	
	Cell A ( tag::vertex );  x (A) = 0.;  y (A) = 0.;
	Cell B ( tag::vertex );  x (B) = 1.;  y (B) = 0.;
	Cell C ( tag::vertex );  x (C) = 1.;  y (C) = 1.;
	Cell D ( tag::vertex );  x (D) = 0.;  y (D) = 1.;
	
	Mesh AB ( tag::segment, A .reverse(), B, tag::divided_in, 20 );
	Mesh BC ( tag::segment, B .reverse(), C, tag::divided_in, 20 );
	Mesh CD ( tag::segment, C .reverse(), D, tag::divided_in, 20 );
	Mesh DA ( tag::segment, D .reverse(), A, tag::divided_in, 20 );
	
	Mesh ABCD ( tag::rectangle, AB, BC, CD, DA) ;

	// we choose a function psi
	// and cut a square in two triangles if the values of psi
	// at two opposite vertices are close to zero
	
	Function psi = (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) - 0.09 ;
	double tol = 1.01e-2;

	std::list < Cell > list_of_squares;

	Mesh::Iterator it = ABCD .iterator ( tag::over_cells_of_dim, 2 );
	for ( it .reset() ; it .in_range(); it++ )
	{	Cell square = *it;  
		Mesh::Iterator itsq = square .boundary() .iterator ( tag::over_vertices );
		itsq .reset();  assert ( itsq .in_range() );
		Cell P = * itsq;  // first vertex of the square
		itsq ++;  assert ( itsq .in_range() );
		Cell Q = * itsq;  // second vertex of the square
		itsq ++;  assert ( itsq .in_range() );
		Cell R = * itsq;  // third vertex of the square
		itsq ++;  assert ( itsq .in_range() );
		Cell S = * itsq;  // fourth vertex of the square
		itsq ++;  assert ( not itsq .in_range() );
		if ( ( ( std::abs ( psi (P) ) < tol ) and ( std::abs ( psi (R) ) < tol ) ) or
		     ( ( std::abs ( psi (Q) ) < tol ) and ( std::abs ( psi (S) ) < tol ) )    )
			list_of_squares .push_back ( square );                                            }

	std::list < Cell > ::iterator it_list;
	for ( it_list = list_of_squares .begin(); it_list != list_of_squares .end(); it_list ++ )
	{	Cell square = *it_list;
		Mesh::Iterator itsq = square .boundary() .iterator ( tag::over_vertices );
		itsq .reset();  assert ( itsq .in_range() );
		Cell P = * itsq;  // first vertex of the square
		itsq ++;  assert ( itsq .in_range() );
		Cell Q = * itsq;  // second vertex of the square
		itsq ++;  assert ( itsq .in_range() );
		Cell R = * itsq;  // third vertex of the square
		itsq ++;  assert ( itsq .in_range() );
		Cell S = * itsq;  // fourth vertex of the square
		itsq ++;  assert ( not itsq .in_range() );
		Cell PQ = square .boundary() .cell_in_front_of ( P );
		assert ( PQ .tip() == Q );
		Cell QR = square .boundary() .cell_in_front_of ( Q );
		assert ( QR .tip() == R );
		Cell RS = square .boundary() .cell_in_front_of ( R );
		assert ( RS .tip() == S );
		Cell SP = square .boundary() .cell_in_front_of ( S );
		assert ( SP .tip() == P );
		if ( ( std::abs ( psi (P) ) < tol ) and ( std::abs ( psi (R) ) < tol ) )
		{	// we cut along diagonal PR
			square .remove_from_mesh ( ABCD );
			Cell PR ( tag::segment, P .reverse(), R );
			Cell PQR ( tag::triangle, PQ, QR, PR .reverse() );
			PQR .add_to_mesh ( ABCD );
			Cell RSP ( tag::triangle, RS, SP, PR );
			RSP .add_to_mesh ( ABCD );                               }
		else if ( ( std::abs ( psi (Q) ) < tol ) and ( std::abs ( psi (S) ) < tol ) )
		{	// we cut along diagonal QS
			square .remove_from_mesh ( ABCD );
			Cell QS ( tag::segment, Q .reverse(), S );
			Cell PQS ( tag::triangle, PQ, QS, SP );
			PQS .add_to_mesh ( ABCD );
			Cell QRS ( tag::triangle, QR, RS, QS .reverse() );
			QRS .add_to_mesh ( ABCD );                               }
		else assert ( false );                                                              }
	
	ABCD .export_to_file ( tag::msh, "cut-squares.msh");
	std::cout << "produced file cut-squares.msh" << std::endl;

}  // end of main
