

#include "maniFEM.h"

using namespace maniFEM;


	
int main ( )

{	// begin with the usual two-dimensional space
	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0], y = xy[1];

	Cell A ( tag::vertex );  x(A) =  0.  ;  y(A) = -0.35;
	Cell B ( tag::vertex );  x(B) =  0.  ;  y(B) = -1.05;
	Cell C ( tag::vertex );  x(C) =  1.15;  y(C) = -2.15;
	Cell D ( tag::vertex );  x(D) =  2.3 ;  y(D) = -1.05;
	Cell E ( tag::vertex );  x(E) =  2.3 ;  y(E) = -0.35;
	Cell F ( tag::vertex );  x(F) =  1.15;  y(F) = -1.45;
	Cell G ( tag::vertex );  x(G) =  1.15;  y(G) = -0.75;
	Cell H ( tag::vertex );  x(H) =  0.  ;  y(H) =  0.35;
	Cell I ( tag::vertex );  x(I) = -1.15;  y(I) = -0.75;
	Cell J ( tag::vertex );  x(J) = -1.15;  y(J) = -1.45;

	Mesh AB ( tag::segment, A.reverse(), B, tag::divided_in, 14 );
	Mesh BC ( tag::segment, B.reverse(), C, tag::divided_in, 32 );
	Mesh CD ( tag::segment, C.reverse(), D, tag::divided_in, 32 );
	Mesh DE ( tag::segment, D.reverse(), E, tag::divided_in, 14 );
	Mesh EF ( tag::segment, E.reverse(), F, tag::divided_in, 32 );
	Mesh FG ( tag::segment, F.reverse(), G, tag::divided_in, 14 );
	Mesh GH ( tag::segment, G.reverse(), H, tag::divided_in, 32 );
	Mesh HI ( tag::segment, H.reverse(), I, tag::divided_in, 32 );
	Mesh IJ ( tag::segment, I.reverse(), J, tag::divided_in, 14 );
	Mesh JA ( tag::segment, J.reverse(), A, tag::divided_in, 32 );

	Mesh zig_zag ( tag::join, { AB, BC, CD, DE, EF, FG, GH, HI, IJ, JA } );

	const double e = 1.5;
	Manifold curve1 = RR2.implicit 
		( smooth_min ( 300.*power((x+y)*(x+y),e) + power((x-y-1.)*(x-y-1.),e),
		               300.*power((x-y)*(x-y),e) + power((x+y+1.)*(x+y+1.),e),
	                 tag::threshold, 1.5                     )  == 1. );

	Mesh loop1 ( tag::progressive, tag::desired_length, 0.05 );

	double a = 0.15, b = 0.82;
	Manifold curve2 = RR2.implicit 
		( smooth_min ( 300.*power((x-a+y+b)*(x-a+y+b),e) + power((x-a-y-b-1.)*(x-a-y-b-1.),e),
		               300.*power((x-2-a-y-b)*(x-2-a-y-b),e) + power((x-2-a+y+b+1.)*(x-2-a+y+b+1.),e),
	                 tag::threshold, 1.5                     )  == 1. );

	Mesh loop2 ( tag::progressive, tag::desired_length, 0.05 );

	b = 1.4;
	Manifold curve3 = RR2.implicit 
		( smooth_min ( 300.*power((x+y+b)*(x+y+b),e) + power((x-y-b-1.)*(x-y-b-1.),e),
		               300.*power((x-y-b)*(x-y-b),e) + power((x+y+b+1.)*(x+y+b+1.),e),
	                 tag::threshold, 1.5                     )  == 1. );

	Mesh loop3 ( tag::progressive, tag::desired_length, 0.05 );

	b = 0.82-1.4;
	Manifold curve4 = RR2.implicit 
		( smooth_min ( 300.*power((x-a+y+b)*(x-a+y+b),e) + power((x-a-y-b-1.)*(x-a-y-b-1.),e),
		               300.*power((x-2-a-y-b)*(x-2-a-y-b),e) + power((x-2-a+y+b+1.)*(x-2-a+y+b+1.),e),
	                 tag::threshold, 1.5                     )  == 1. );

	Mesh loop4 ( tag::progressive, tag::desired_length, 0.05 );

	a -= 2.3;  b = 0.82;
	Manifold curve5 = RR2.implicit 
		( smooth_min ( 300.*power((x-a+y+b)*(x-a+y+b),e) + power((x-a-y-b-1.)*(x-a-y-b-1.),e),
		               300.*power((x-2-a-y-b)*(x-2-a-y-b),e) + power((x-2-a+y+b+1.)*(x-2-a+y+b+1.),e),
	                 tag::threshold, 1.5                     )  == 1. );

	Mesh loop5 ( tag::progressive, tag::desired_length, 0.05 );

	Mesh loop ( tag::join, loop1.reverse(), loop2.reverse(), zig_zag );

	RR2.set_as_working_manifold();
	Mesh domain ( tag::progressive, tag::boundary, loop, tag::desired_length, 0.05 );

	Mesh IJAB ( tag::join, IJ, JA, AB );
	Mesh GFED ( tag::join, DE.reverse(), EF.reverse(), FG.reverse() );

	Mesh VV = domain.fold ( tag::identify, IJAB, tag::with, GFED,
                          tag::identify, BC, tag::with, GH.reverse(),
													tag::identify, CD, tag::with, HI.reverse(),
													tag::use_existing_vertices                 );

	VV.draw_ps ( "VV.eps", tag::unfold, tag::over_region, -1.5 < x < 2.5, -3 < y < 0.5 );
}

