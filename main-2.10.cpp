
// example presented in paragraph 2.10 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// a physalis-like surface

#include "maniFEM.h"

using namespace maniFEM;


int main ( )

{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz[0], y = xyz[1], z = xyz[2];

	Function r2 = x*x + y*y + z*z;
	const double pi = 3.1415926536;
	Manifold apple = RR3.implicit ( power(r2,0.5) * sin(r2-pi/6.) == z );

	Cell A ( tag::vertex );  x(A) =  0.;  y(A) =  0.;  z(A) = std::sqrt ( 2.*pi/3. );
	Cell B1 ( tag::vertex );  x(B1) =  1.;  y(B1) =  0.;  z(B1) =  1.;
	Cell C1 ( tag::vertex );  x(C1) =  1.;  y(C1) =  0.;  z(C1) =  0.;
	apple.project (B1);  apple.project (C1);
	Cell D ( tag::vertex );  x(D) =  0.;  y(D) =  0.;  z(D) =  0.;
	Mesh AB1 ( tag::segment, A.reverse(), B1, tag::divided_in, 10 );
	Mesh B1C1 ( tag::segment, B1.reverse(), C1, tag::divided_in, 10 );
	Mesh C1D ( tag::segment, C1.reverse(), D, tag::divided_in, 10 );
	Cell B2 ( tag::vertex );  x(B2) =  0.707;  y(B2) =  0.707;  z(B2) =  1.;
	Cell C2 ( tag::vertex );  x(C2) =  0.707;  y(C2) =  0.707;  z(C2) =  0.;
	apple.project (B2);  apple.project (C2);
	Mesh AB2 ( tag::segment, A.reverse(), B2, tag::divided_in, 10 );
	Mesh B2C2 ( tag::segment, B2.reverse(), C2, tag::divided_in, 10 );
	Mesh C2D ( tag::segment, C2.reverse(), D, tag::divided_in, 10 );
	Cell B3 ( tag::vertex );  x(B3) =  0.;  y(B3) =  1.;  z(B3) =  1.;
	Cell C3 ( tag::vertex );  x(C3) =  0.;  y(C3) =  1.;  z(C3) =  0.;
	apple.project (B3);  apple.project (C3);
	Mesh AB3 ( tag::segment, A.reverse(), B3, tag::divided_in, 10 );
	Mesh B3C3 ( tag::segment, B3.reverse(), C3, tag::divided_in, 10 );
	Mesh C3D ( tag::segment, C3.reverse(), D, tag::divided_in, 10 );
	Cell B4 ( tag::vertex );  x(B4) = -0.707;  y(B4) =  0.707;  z(B4) =  1.;
	Cell C4 ( tag::vertex );  x(C4) = -0.707;  y(C4) =  0.707;  z(C4) =  0.;
	apple.project (B4);  apple.project (C4);
	Mesh AB4 ( tag::segment, A.reverse(), B4, tag::divided_in, 10 );
	Mesh B4C4 ( tag::segment, B4.reverse(), C4, tag::divided_in, 10 );
	Mesh C4D ( tag::segment, C4.reverse(), D, tag::divided_in, 10 );
	Cell B5 ( tag::vertex );  x(B5) =  -1.;  y(B5) =  0.;  z(B5) =  1.;
	Cell C5 ( tag::vertex );  x(C5) =  -1.;  y(C5) =  0.;  z(C5) =  0.;
	apple.project (B5);  apple.project (C5);
	Mesh AB5 ( tag::segment, A.reverse(), B5, tag::divided_in, 10 );
	Mesh B5C5 ( tag::segment, B5.reverse(), C5, tag::divided_in, 10 );
	Mesh C5D ( tag::segment, C5.reverse(), D, tag::divided_in, 10 );
	Cell B6 ( tag::vertex );  x(B6) =  -0.707;  y(B6) =  -0.707;  z(B6) =  1.;
	Cell C6 ( tag::vertex );  x(C6) =  -0.707;  y(C6) =  -0.707;  z(C6) =  0.;
	apple.project (B6);  apple.project (C6);
	Mesh AB6 ( tag::segment, A.reverse(), B6, tag::divided_in, 10 );
	Mesh B6C6 ( tag::segment, B6.reverse(), C6, tag::divided_in, 10 );
	Mesh C6D ( tag::segment, C6.reverse(), D, tag::divided_in, 10 );
	Cell B7 ( tag::vertex );  x(B7) =  0.;  y(B7) =  -1.;  z(B7) =  1.;
	Cell C7 ( tag::vertex );  x(C7) =  0.;  y(C7) =  -1.;  z(C7) =  0.;
	apple.project (B7);  apple.project (C7);
	Mesh AB7 ( tag::segment, A.reverse(), B7, tag::divided_in, 10 );
	Mesh B7C7 ( tag::segment, B7.reverse(), C7, tag::divided_in, 10 );
	Mesh C7D ( tag::segment, C7.reverse(), D, tag::divided_in, 10 );
	Cell B8 ( tag::vertex );  x(B8) =  0.707;  y(B8) =  -0.707;  z(B8) =  1.;
	Cell C8 ( tag::vertex );  x(C8) =  0.707;  y(C8) =  -0.707;  z(C8) =  0.;
	apple.project (B8);  apple.project (C8);
	Mesh AB8 ( tag::segment, A.reverse(), B8, tag::divided_in, 10 );
	Mesh B8C8 ( tag::segment, B8.reverse(), C8, tag::divided_in, 10 );
	Mesh C8D ( tag::segment, C8.reverse(), D, tag::divided_in, 10 );

	RR3.set_as_working_manifold();
	Mesh B1B2 ( tag::segment, B1.reverse(), B2, tag::divided_in, 10 );
	Mesh B2B3 ( tag::segment, B2.reverse(), B3, tag::divided_in, 10 );
	Mesh B3B4 ( tag::segment, B3.reverse(), B4, tag::divided_in, 10 );
	Mesh B4B5 ( tag::segment, B4.reverse(), B5, tag::divided_in, 10 );
	Mesh B5B6 ( tag::segment, B5.reverse(), B6, tag::divided_in, 10 );
	Mesh B6B7 ( tag::segment, B6.reverse(), B7, tag::divided_in, 10 );
	Mesh B7B8 ( tag::segment, B7.reverse(), B8, tag::divided_in, 10 );
	Mesh B8B1 ( tag::segment, B8.reverse(), B1, tag::divided_in, 10 );
	Mesh C1C2 ( tag::segment, C1.reverse(), C2, tag::divided_in, 10 );
	Mesh C2C3 ( tag::segment, C2.reverse(), C3, tag::divided_in, 10 );
	Mesh C3C4 ( tag::segment, C3.reverse(), C4, tag::divided_in, 10 );
	Mesh C4C5 ( tag::segment, C4.reverse(), C5, tag::divided_in, 10 );
	Mesh C5C6 ( tag::segment, C5.reverse(), C6, tag::divided_in, 10 );
	Mesh C6C7 ( tag::segment, C6.reverse(), C7, tag::divided_in, 10 );
	Mesh C7C8 ( tag::segment, C7.reverse(), C8, tag::divided_in, 10 );
	Mesh C8C1 ( tag::segment, C8.reverse(), C1, tag::divided_in, 10 );
	
	Mesh AB1B2 ( tag::triangle, AB1, B1B2, AB2.reverse() );
	Mesh AB2B3 ( tag::triangle, AB2, B2B3, AB3.reverse() );
	Mesh AB3B4 ( tag::triangle, AB3, B3B4, AB4.reverse() );
	Mesh AB4B5 ( tag::triangle, AB4, B4B5, AB5.reverse() );
	Mesh AB5B6 ( tag::triangle, AB5, B5B6, AB6.reverse() );
	Mesh AB6B7 ( tag::triangle, AB6, B6B7, AB7.reverse() );
	Mesh AB7B8 ( tag::triangle, AB7, B7B8, AB8.reverse() );
	Mesh AB8B1 ( tag::triangle, AB8, B8B1, AB1.reverse() );

	Mesh B1C1C2B2 ( tag::quadrangle, B1C1, C1C2, B2C2.reverse(), B1B2.reverse(), tag::with_triangles );
	Mesh B2C2C3B3 ( tag::quadrangle, B2C2, C2C3, B3C3.reverse(), B2B3.reverse(), tag::with_triangles );
	Mesh B3C3C4B4 ( tag::quadrangle, B3C3, C3C4, B4C4.reverse(), B3B4.reverse(), tag::with_triangles );
	Mesh B4C4C5B5 ( tag::quadrangle, B4C4, C4C5, B5C5.reverse(), B4B5.reverse(), tag::with_triangles );
	Mesh B5C5C6B6 ( tag::quadrangle, B5C5, C5C6, B6C6.reverse(), B5B6.reverse(), tag::with_triangles );
	Mesh B6C6C7B7 ( tag::quadrangle, B6C6, C6C7, B7C7.reverse(), B6B7.reverse(), tag::with_triangles );
	Mesh B7C7C8B8 ( tag::quadrangle, B7C7, C7C8, B8C8.reverse(), B7B8.reverse(), tag::with_triangles );
	Mesh B8C8C1B1 ( tag::quadrangle, B8C8, C8C1, B1C1.reverse(), B8B1.reverse(), tag::with_triangles );

	Mesh C1DC2 ( tag::triangle, C1D, C2D.reverse(), C1C2.reverse() );
	Mesh C2DC3 ( tag::triangle, C2D, C3D.reverse(), C2C3.reverse() );
	Mesh C3DC4 ( tag::triangle, C3D, C4D.reverse(), C3C4.reverse() );
	Mesh C4DC5 ( tag::triangle, C4D, C5D.reverse(), C4C5.reverse() );
	Mesh C5DC6 ( tag::triangle, C5D, C6D.reverse(), C5C6.reverse() );
	Mesh C6DC7 ( tag::triangle, C6D, C7D.reverse(), C6C7.reverse() );
	Mesh C7DC8 ( tag::triangle, C7D, C8D.reverse(), C7C8.reverse() );
	Mesh C8DC1 ( tag::triangle, C8D, C1D.reverse(), C8C1.reverse() );

	Mesh sect1 ( tag::join, AB1B2, B1C1C2B2, C1DC2 );
	Mesh sect2 ( tag::join, AB2B3, B2C2C3B3, C2DC3 );
	Mesh sect3 ( tag::join, AB3B4, B3C3C4B4, C3DC4 );
	Mesh sect4 ( tag::join, AB4B5, B4C4C5B5, C4DC5 );
	Mesh sect5 ( tag::join, AB5B6, B5C5C6B6, C5DC6 );
	Mesh sect6 ( tag::join, AB6B7, B6C6C7B7, C6DC7 );
	Mesh sect7 ( tag::join, AB7B8, B7C7C8B8, C7DC8 );
	Mesh sect8 ( tag::join, AB8B1, B8C8C1B1, C8DC1 );

	std::list < Mesh > lm { sect1, sect2, sect3, sect4, sect5, sect6, sect7, sect8 };
	Mesh fisalis ( tag::join, lm ); 
	
	CellIterator it = fisalis.iter_over ( tag::vertices );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell P = *it;
		x(P) *= 0.8;  y(P) *= 0.8;
		if ( z(P) > 1.3 )
		{	x(P) = x(P) / ( 1. + 300. * std::pow ( z(P) - 1.3, 3. ) );
			y(P) = y(P) / ( 1. + 300. * std::pow ( z(P) - 1.3, 3. ) );
			z(P) = z(P) * ( 1. + 10. * ( z(P) - 1.3 ) * ( z(P) - 1.3 ) );  }
		if ( z(P) > 0. ) z(P) *= 0.8;                                       }

	std::list<Mesh>::iterator it1;
	for ( it1 = lm.begin(); it1 != lm.end(); it1++ )
	{	Mesh sect = *it1;
		CellIterator it2 = sect.iter_over ( tag::cells_of_dim, 1 );
		for ( int i = 1; i < 20; i++ )
		for ( it2.reset(); it2.in_range(); it2++ )
		{	Cell seg = *it2;
			sect.baricenter ( seg.tip(), seg );
			seg = seg.reverse();
			sect.baricenter ( seg.tip(), seg );  }                     }

	// does the 'baricenter' method really need a segment as second argument ?
	// defined in manifold.h
	
	fisalis.export_msh ("physalis.msh");

	std::cout << "produced file physalis.msh" << std::endl;
}
