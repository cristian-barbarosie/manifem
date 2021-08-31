

#include "maniFEM.h"
using namespace maniFEM;


void print_segment ( Cell seg )

{	Manifold::Quotient * manif_q = dynamic_cast
		< Manifold::Quotient* > ( Manifold::working.core );
	assert ( manif_q );
	Function xy = manif_q->base_space.coordinates();
	Function x = xy[0], y = xy[1];
	size_t n = manif_q->actions.size();
	assert ( n == manif_q->spins.size() );
	Cell V = seg.base().reverse();
	Cell W = seg.tip();
	std::cout << "[(" << x(V) << "," << y(V) << "),(" << x(W) << "," <<y(W) << ")] ";
	Function::CompositionOfActions a = seg.spin();
	std::cout << "(";
	for ( size_t i = 0; i < n; i++ )
	{	Function::Action & g = manif_q->actions[i];
		std::map<Function::Action,short int>::const_iterator itt = a.index_map.find ( g );
		if ( itt == a.index_map.end() )
		{	std::cout << "0,"; continue;  }
		short int exp = itt->second;
		assert ( exp != 0 );
		std::cout << exp << ",";                                                            }
	std::cout << ")" << std::endl;                                                           }
	

void print_spin ( Function::CompositionOfActions a )

{	Manifold::Quotient * manif_q = dynamic_cast
		< Manifold::Quotient* > ( Manifold::working.core );
	assert ( manif_q );
	Function xy = manif_q->base_space.coordinates();
	Function x = xy[0], y = xy[1];
	size_t n = manif_q->actions.size();
	assert ( n == manif_q->spins.size() );
	std::cout << "(";
	for ( size_t i = 0; i < n; i++ )
	{	Function::Action & g = manif_q->actions[i];
		std::map<Function::Action,short int>::const_iterator itt = a.index_map.find ( g );
		if ( itt == a.index_map.end() )
		{	std::cout << "0,"; continue;  }
		short int exp = itt->second;
		assert ( exp != 0 );
		std::cout << exp << ",";                                                            }
	std::cout << ")" << std::endl;                                                           }


int main ( )

{	// begin with the usual two-dimensional space
	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0], y = xy[1];

	// define two actions on RR2 (translations)
	Function::Action g_horiz ( tag::transforms, xy, tag::into, (x+1.) && y ),
	                 g_vert ( tag::transforms, xy, tag::into, x && (y+1.) );

	// and divide RR2 by the group of translations generated by {g1,g2}
	Manifold torus_manif = RR2.quotient ( g_horiz, g_vert );

	// try with four rectangles 
	
	Cell A ( tag::vertex );  x(A) = 0.;   y(A) = 0.;
	Cell B ( tag::vertex );  x(B) = 0.5;  y(B) = 0.;
	Cell C ( tag::vertex );  x(C) = 0.5;  y(C) = 0.5;
	Cell D ( tag::vertex );  x(D) = 0.;   y(D) = 0.5;

	Mesh AB ( tag::segment, A.reverse(), B, tag::divided_in, 5 ),
	     BC ( tag::segment, B.reverse(), C, tag::divided_in, 5 ),
	     CD ( tag::segment, C.reverse(), D, tag::divided_in, 5 ),
	     DA ( tag::segment, D.reverse(), A, tag::divided_in, 5 ),
	     BA1 ( tag::segment, B.reverse(), A, tag::divided_in, 5, tag::spin, g_horiz ),
	     CD1 ( tag::segment, C.reverse(), D, tag::divided_in, 5, tag::spin, g_horiz ),
	     CB2 ( tag::segment, C.reverse(), B, tag::divided_in, 5, tag::spin, g_vert ),
	     DA2 ( tag::segment, D.reverse(), A, tag::divided_in, 5, tag::spin, g_vert );

	Mesh sq1 ( tag::rectangle, AB, BC, CD, DA ),
	     sq2 ( tag::rectangle, BA1, DA.reverse(), CD1.reverse(), BC.reverse(), tag::spin ),
	     sq3 ( tag::rectangle, DA2.reverse(), CD.reverse(), CB2, AB.reverse(), tag::spin ),
	     sq4 ( tag::rectangle, CD1, DA2, BA1.reverse(), CB2.reverse(), tag::spin );
	
	exit ( 0 );

	Mesh torus ( tag::join, sq1, sq3 );

	std::cout << "sq3" << std::endl;
	CellIterator it = sq3.iterator ( tag::over_segments );
	for ( it.reset(); it.in_range(); it++ ) print_segment ( *it );
	exit ( 0 );
	
	sq3.draw_ps ( "torus.eps", tag::unfold,
	                tag::over_region, -0.45 < x < 1.55, -0.45 < y < 1.55 );


	std::cout << "produced file torus.eps - please edit it before viewing" << std::endl;

}
