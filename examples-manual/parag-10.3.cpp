
// example presented in paragraph 10.3 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// in a mesh of squares, cuts some of the squares in two triangles


#include "maniFEM.h"
#include <fstream>
#include <iostream>

using namespace std;
using namespace maniFEM;


int main ( )

{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz[0], y = xyz[1], z = xyz[2] ;
	
	Cell A ( tag::vertex ); x(A) = 0. ; y(A) = 0.; z(A) = 0. ;
	Cell B ( tag::vertex ); x(B) = 1. ; y(B) = 0.; z(B) = 0. ;
	Cell C ( tag::vertex ); x(C) = 1. ; y(C) = 1.; z(C) = 0. ;
	Cell D ( tag::vertex ); x(D) = 0. ; y(D) = 1.; z(D) = 0. ;
	
	Mesh AB ( tag::segment, A.reverse(), B, tag::divided_in, 20);
	Mesh BC ( tag::segment, B.reverse(), C, tag::divided_in, 20);
	Mesh CD ( tag::segment, C.reverse(), D, tag::divided_in, 20);
	Mesh DA ( tag::segment, D.reverse(), A, tag::divided_in, 20);
	
	Mesh ABCD ( tag::rectangle, AB, BC, CD, DA) ;
	
	Function Psi = (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) - 0.09 ;

	//	std::map < Cell::Core *, size_t > numbering;
	//	{ // just a block of code for hiding 'it' and 'counter'
	//	CellIterator it = ABCD.iterator ( tag::over_vertices );
	//	size_t counter = 0;
	//	for ( it.reset() ; it.in_range(); it++ )
	//	{	Cell p = *it;  ++counter;  numbering [p.core] = counter;  }
	//	} // just a block of code 

	//	{ // just a block of code for hiding 'it' and 'counter'
	//	CellIterator it = ABCD.iterator ( tag::over_vertices );
	//	for ( it.reset() ; it.in_range(); it++ )
	//	{	Cell p = *it;  z(p) = Psi(p) ;  }
	//	} // just a block of code 
	
	//	ABCD.export_msh ("psi.msh", numbering);
	//	std::cout << "produced file psi.msh" << std::endl;
	
//	{ // just a block of code for hiding variables
//	ofstream solution_file ("square.msh", fstream::app );
//	solution_file << "$NodeData" << endl;
//	solution_file << "1" << endl;   // one string follows
//	solution_file << "\"temperature\"" << endl;
//	solution_file << "1" << endl;   //  one real follows
//	solution_file << "0.0" << endl;  // time [??]
//	solution_file << "3" << endl;   // three integers follow
//	solution_file << "0" << endl;   // time step [??]
//	solution_file << "1" << endl;  // scalar values of u
//	solution_file << ABCD.number_of ( tag::vertices ) << endl;  // number of values listed below
//	CellIterator it = ABCD.iterator ( tag::over_vertices );
//	for ( it.reset(); it.in_range(); it++ )
//	{	Cell P = *it;
//		size_t i = numbering[P.core];
//		solution_file << i << " " << Psi(P) << std::endl;   }
//	} // just a block of code
	

//////////////////////////////////////////////////////////////
	list<Cell> list_of_squares, list_of_vertices;
//////////////////////////////////////////////////////////////


	{ // just a block of code for hiding 'it' and 'counter'
	CellIterator it = ABCD.iterator ( tag::over_cells_of_dim, 2 );
	for ( it.reset() ; it.in_range(); it++ )
	{	Cell square = *it;  
		CellIterator itsq = square.boundary().iterator ( tag::over_vertices );
		int contador = 0; 
		for ( itsq.reset() ; itsq.in_range(); itsq++ )
		{	Cell p = *itsq;  
			if (Psi(p)<0.) contador++ ;  }
//		cout << contador << endl;
		if ((contador==4) or (contador == 0)) continue ;
		CellIterator itside = square.boundary().iterator ( tag::over_segments );
		int conta = 0 , acerta = 0 ;
		double alpha[2], beta[2] ;
		int guarda_aresta[2] ;
		for ( itside.reset() ; itside.in_range(); itside++ )
		{	conta++ ;
			Cell seg = *itside;  Cell P = seg.base().reverse() ; Cell Q = seg.tip() ;  	
			if ((Psi(P) <= 0. and Psi(Q) > 0.) or (Psi(P) >= 0. and Psi(Q) < 0.))
			{	assert(acerta < 2);
				double lambda = - Psi(P)/Psi(Q) ; 
				beta[acerta] = 1./(1. + lambda) ;
				alpha[acerta] = beta[acerta] * lambda ;
				guarda_aresta[acerta] = conta ;
				acerta++ ;
			}
		}
		assert (acerta == 2);
		double porcent;
		if (guarda_aresta[1]==guarda_aresta[0]+2) // opposite sides
		{
			porcent= (alpha[0]+beta[1])/2;
		}
		else 
		{	if (guarda_aresta[1] == guarda_aresta[0]+1) porcent = beta[0]*alpha[1]/2.;
			else{
				assert( guarda_aresta[0] == 1 );
				assert( guarda_aresta[1] == 4 );
				porcent = beta[1]*alpha[0]/2.;
			}
		}
		if ( porcent > 0.65 or porcent < 0.35 ) continue;
		map <int,double> abs_Psi;
		contador = 0;
		for ( itsq.reset() ; itsq.in_range(); itsq++ )
    {	Cell p = *itsq;  
			abs_Psi[contador] = abs(Psi(p)) ;
		contador++ ; }
		double max =abs_Psi[0];
		for (contador =1; contador <4; contador++)
		{	if (abs_Psi[contador] > max)
			{	acerta = contador;
				max = abs_Psi[contador];
			} 
		}
		abs_Psi.erase(acerta);
		assert ( abs_Psi.size() == 3);
		max = -1;
		for ( map<int,double>::iterator itt = abs_Psi.begin(); itt != abs_Psi.end(); itt++)
		{	double val = itt -> second;
			if( val > max )
			{	max = val;
				acerta = itt -> first;
			}
		}
		abs_Psi.erase ( acerta );
		assert ( abs_Psi.size() == 2);
		//		cout << "this square will be cut in two halves"<< " " ;
		//		for ( map<int,double>::iterator itt = abs_Psi.begin(); itt != abs_Psi.end(); itt++)
		//			cout << itt->first << " " ;
		//		cout << endl;
		map<int,double>::iterator itt = abs_Psi.begin(); 


//////////////////////////////////////////////////////////////
		itsq.reset() ; assert(itsq.in_range());
		Cell P=*itsq; itsq++; assert(itsq.in_range());
		Cell PQ = square.boundary().cell_in_front_of(P);
		Cell Q = PQ.tip();
		if ( itt->first == 0 ) 
		{	// cut along diagonal PR
			list_of_squares.push_back(square);
			list_of_vertices.push_back(P);
		}
		else
		{	assert ( itt-> first == 1 );
			// cut along diagonal QS
			list_of_squares.push_back(square);
			list_of_vertices.push_back(Q);
		}
	}
	//	cout << "two lists, of lengths " << list_of_squares.size() << " "
	//	     << list_of_vertices.size() << endl;
	for ( list<Cell>::iterator it1 = list_of_squares.begin(),
        it2 = list_of_vertices.begin(); it1 != list_of_squares.end(); it1++, it2++ )
	{	assert ( it2 != list_of_vertices.end() );
		Cell square = *it1;
		Cell P = *it2;
		Cell PQ = square.boundary().cell_in_front_of(P);
		Cell Q = PQ.tip();
		Cell QR = square.boundary().cell_in_front_of(Q);
		Cell R = QR.tip();
		Cell RS = square.boundary().cell_in_front_of(R);
		Cell S = RS.tip();
		Cell SP = square.boundary().cell_in_front_of(S);
		assert ( SP.tip() == P );
		square.remove_from_mesh ( ABCD );
		Cell PR ( tag::segment, P.reverse(), R );
		Cell PQR ( tag::triangle, PQ, QR, PR.reverse() );
		PQR.add_to_mesh ( ABCD );
		Cell RSP ( tag::triangle, RS, SP, PR );
		RSP.add_to_mesh ( ABCD );                            }
	} // just a block of code 
////////////////////////////////////////////////////////////////

	
	ABCD.export_msh ( "cut-squares.msh" );
	std::cout << "produced file cut-squares.msh" << std::endl;
}
