
// example presented in paragraph 7.21 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// folds a hexagon with a star-shaped hole

#include "maniFEM.h"
using namespace maniFEM;


int main ( )

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0], y = xy[1];

	const double coef = 0.7;

	Function d_AA = ( x - 0.75*coef ) * ( x - 0.75*coef ) + ( y - 0.43*coef) * ( y - 0.43*coef);
	Function d_BB = x * x + ( y - 0.86*coef) * ( y - 0.86*coef);
	Function d_CC = ( x + 0.75*coef ) * ( x + 0.75*coef ) + ( y - 0.43*coef) * ( y - 0.43*coef);
	Function d_DD = ( x + 0.75*coef ) * ( x + 0.75*coef ) + ( y + 0.43*coef) * ( y + 0.43*coef);
	Function d_EE = x * x + ( y + 0.86*coef) * ( y + 0.86*coef);
	Function d_FF = ( x - 0.75*coef ) * ( x - 0.75*coef ) + ( y + 0.43*coef) * ( y + 0.43*coef);
	// Function d = 0.05 + 0.055 * power ( d_AA * d_BB * d_CC * d_DD * d_EE * d_FF, 0.15 );
	Function dmin = smooth_min ( d_AA, d_BB, d_CC, d_DD, d_EE, d_FF, tag::threshold, 0.01 );
	Function d = 0.15 * power ( 0.001 + dmin, 0.25 ) - 0.002;

	Cell AA ( tag::vertex );  x(AA) =  0.75 * coef;  y(AA) =  0.43 * coef;
	Cell BB ( tag::vertex );  x(BB) =  0.         ;  y(BB) =  0.86 * coef;
	Cell CC ( tag::vertex );  x(CC) = -0.75 * coef;  y(CC) =  0.43 * coef;
	Cell DD ( tag::vertex );  x(DD) = -0.75 * coef;  y(DD) = -0.43 * coef;
	Cell EE ( tag::vertex );  x(EE) =  0.         ;  y(EE) = -0.86 * coef;
	Cell FF ( tag::vertex );  x(FF) =  0.75 * coef;  y(FF) = -0.43 * coef;

	const double r2 = (1.-0.75*coef)*(1.-0.75*coef) + 0.43*0.43*coef*coef;
	RR2 .implicit ( (x-1.)*(x-1.) + y*y == r2 );
	Mesh AF ( tag::frontal, tag::start_at, AA, tag::stop_at, FF,
	          tag::desired_length, d, tag::shortest_path        );
	RR2 .implicit ( (x-0.5)*(x-0.5) + (y+0.86)*(y+0.86) == r2 );
	Mesh FE ( tag::frontal, tag::start_at, FF, tag::stop_at, EE,
	          tag::desired_length, d, tag::shortest_path        );
	RR2 .implicit ( (x+0.5)*(x+0.5) + (y+0.86)*(y+0.86) == r2 );
	Mesh ED ( tag::frontal, tag::start_at, EE, tag::stop_at, DD,
	          tag::desired_length, d, tag::shortest_path        );
	RR2 .implicit ( (x+1.)*(x+1.) + y*y == r2 );
	Mesh DC ( tag::frontal, tag::start_at, DD, tag::stop_at, CC,
	          tag::desired_length, d, tag::shortest_path        );
	RR2 .implicit ( (x+0.5)*(x+0.5) + (y-0.86)*(y-0.86) == r2 );
	Mesh CB ( tag::frontal, tag::start_at, CC, tag::stop_at, BB,
	          tag::desired_length, d, tag::shortest_path        );
	RR2 .implicit ( (x-0.5)*(x-0.5) + (y-0.86)*(y-0.86) == r2 );
	Mesh BA ( tag::frontal, tag::start_at, BB, tag::stop_at, AA,
	          tag::desired_length, d, tag::shortest_path        );

	Cell A ( tag::vertex );  x (A) =  1. ;  y (A) =  0.  ;
	Cell B ( tag::vertex );  x (B) =  0.5;  y (B) =  0.86;
	Cell C ( tag::vertex );  x (C) = -0.5;  y (C) =  0.86;
	Cell D ( tag::vertex );  x (D) = -1. ;  y (D) =  0.  ;
	Cell E ( tag::vertex );  x (E) = -0.5;  y (E) = -0.86;
	Cell F ( tag::vertex );  x (F) =  0.5;  y (F) = -0.86;

	RR2 .implicit ( 2.*x + y/0.86 == 2. );
	Mesh AB ( tag::frontal, tag::start_at, A, tag::stop_at, B,
	          tag::desired_length, d, tag::shortest_path      );
	RR2 .implicit ( y == 0.86 );
	Mesh BC ( tag::frontal, tag::start_at, B, tag::stop_at, C,
	          tag::desired_length, d, tag::shortest_path      );
	RR2 .implicit ( 2.*x - y/0.86 == -2. );
	Mesh CD ( tag::frontal, tag::start_at, C, tag::stop_at, D,
	          tag::desired_length, d, tag::shortest_path      );
	RR2 .implicit ( 2.*x + y/0.86 == -2. );
	Mesh DE ( tag::frontal, tag::start_at, D, tag::stop_at, E,
	          tag::desired_length, d, tag::shortest_path      );
	RR2 .implicit ( y == -0.86 );
	Mesh EF ( tag::frontal, tag::start_at, E, tag::stop_at, F,
	          tag::desired_length, d, tag::shortest_path      );
	RR2 .implicit ( 2.*x - y/0.86 == 2. );
	Mesh FA ( tag::frontal, tag::start_at, F, tag::stop_at, A,
	          tag::desired_length, d, tag::shortest_path      );

	// we must ensure opposite faces have the same number of elements
	// and the corresponding vertices have the same positions
	// otherwise method Mesh::fold will not work

	size_t n_AB = AB .number_of ( tag::segments );
	assert ( n_AB == DE .number_of ( tag::segments ) );
	size_t n_BC = BC .number_of ( tag::segments );
	assert ( n_BC == EF .number_of ( tag::segments ) );
	size_t n_CD = CD .number_of ( tag::segments );
	assert ( n_CD == FA .number_of ( tag::segments ) );
	{ // just a block for hiding it1, it2
	Mesh::Iterator it1 = AB .iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it2 = DE .iterator ( tag::over_vertices, tag::backwards );
	it1 .reset();  assert ( it1 .in_range() );
	it2 .reset();  assert ( it2 .in_range() );
	for ( size_t i = 1; i < n_AB; i++ )
	{	it1++;  assert ( it1 .in_range() );
		it2++;  assert ( it2 .in_range() );
		Cell V = *it1, W = *it2;
		double xx = x(V), yy = y(V);
		xx = ( xx - 1.5 + x(W) ) / 2.;
		yy = ( yy - 0.86 + y(W) ) / 2.;
		x(V) = xx + 1.5;  y(V) = yy + 0.86;
		x(W) = xx;  y(W) = yy;              }
	it1++;  assert ( it1 .in_range() );
	it2++;  assert ( it2 .in_range() );
	it1++;  assert ( not it1 .in_range() );
	it2++;  assert ( not it2 .in_range() );
	} { // just a block for hiding it1, it2
	Mesh::Iterator it1 = BC .iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it2 = EF .iterator ( tag::over_vertices, tag::backwards );
	it1 .reset();  assert ( it1 .in_range() );
	it2 .reset();  assert ( it2 .in_range() );
	for ( size_t i = 1; i < n_BC; i++ )
	{	it1++;  assert ( it1 .in_range() );
		it2++;  assert ( it2 .in_range() );
		Cell V = *it1, W = *it2;
		double xx = ( x(V) + x(W) ) /2.;
		x(V) = xx;  x(W) = xx;           }
	it1++;  assert ( it1 .in_range() );
	it2++;  assert ( it2 .in_range() );
	it1++;  assert ( not it1 .in_range() );
	it2++;  assert ( not it2 .in_range() );
	} { // just a block for hiding it1, it2
	Mesh::Iterator it1 = CD .iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it2 = FA .iterator ( tag::over_vertices, tag::backwards );
	it1 .reset();  assert ( it1 .in_range() );
	it2 .reset();  assert ( it2 .in_range() );
	for ( size_t i = 1; i < n_CD; i++ )
	{	it1++;  assert ( it1 .in_range() );
		it2++;  assert ( it2 .in_range() );
		Cell V = *it1, W = *it2;
		double xx = x(V), yy = y(V);
		xx = ( xx + 1.5 + x(W) ) / 2.;
		yy = ( yy - 0.86 + y(W) ) / 2.;
		x(V) = xx - 1.5;  y(V) = yy + 0.86;
		x(W) = xx;  y(W) = yy;              }
	it1++;  assert ( it1 .in_range() );
	it2++;  assert ( it2 .in_range() );
	it1++;  assert ( not it1 .in_range() );
	it2++;  assert ( not it2 .in_range() );
	} // just a block of code

	Mesh bdry ( tag::join, { AB, BC, CD, DE, EF, FA, AF, FE, ED, DC, CB, BA } );

	RR2 .set_as_working_manifold();
	Mesh hexa ( tag::frontal, tag::boundary, bdry, tag::desired_length, d );

	Mesh torus = hexa .fold ( tag::identify, CD, tag::with, FA .reverse(),
	                          tag::identify, BC, tag::with, EF .reverse(),
	                          tag::identify, AB, tag::with, DE .reverse(),
	                          tag::use_existing_vertices                  );

	{ // just a block for hiding 'it'
	Mesh::Iterator it = torus .iterator ( tag::over_vertices );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;
		if ( V .is_inner_to ( torus ) ) torus .baricenter ( V, tag::winding );  }
	} // just a block for hiding 'it'
	
	std::cout << "produced folded mesh, now drawing, please wait" << std::endl << std::flush;
	
	torus .draw_ps ("torus.eps", tag::unfold,
	                 tag::over_region, -1.5 < x < 3.4, -2.8 < y < 0.4 );

	std::cout << "produced file torus.eps - please edit before viewing" << std::endl;	

}  // end of main


//-----------------------------------------------------------------------------------------


