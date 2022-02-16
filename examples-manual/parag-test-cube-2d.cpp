
#include "maniFEM.h"

using namespace maniFEM;

int main ()

{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz [0], y = xyz [1], z = xyz [2];

	// define three actions on RR3 (translations)
	Manifold::Action gx ( tag::transforms, xyz, tag::into, (x+2.) && y && z ),
	                 gy ( tag::transforms, xyz, tag::into, x && (y+2.) && z ),
	                 gz ( tag::transforms, xyz, tag::into, x && y && (z+2.) );
	
	// and divide RR3 by the group of translations generated by { gx, gy }
	Manifold torus_manif = RR3 .quotient ( gx, gy );

	const double r2 = 0.2;
	Manifold sphere = RR3 .implicit ( x*x + y*y + z*z == r2 );

	Cell AA ( tag::vertex );  x (AA) = -0.5;  y (AA) = -0.5;  z (AA) = -0.5;
	Cell BB ( tag::vertex );  x (BB) = -0.5;  y (BB) =  0.5;  z (BB) = -0.5;
	Cell CC ( tag::vertex );  x (CC) =  0.5;  y (CC) =  0.5;  z (CC) = -0.5;
	Cell DD ( tag::vertex );  x (DD) =  0.5;  y (DD) = -0.5;  z (DD) = -0.5;
	Cell EE ( tag::vertex );  x (EE) = -0.5;  y (EE) = -0.5;  z (EE) =  0.5;
	Cell FF ( tag::vertex );  x (FF) = -0.5;  y (FF) =  0.5;  z (FF) =  0.5;
	Cell GG ( tag::vertex );  x (GG) =  0.5;  y (GG) =  0.5;  z (GG) =  0.5;
	Cell HH ( tag::vertex );  x (HH) =  0.5;  y (HH) = -0.5;  z (HH) =  0.5;
	sphere .project (AA);  sphere .project (BB);  sphere .project (CC);  sphere .project (DD);  
	sphere .project (EE);  sphere .project (FF);  sphere .project (GG);  sphere .project (HH);  

	Mesh AABB ( tag::segment, AA .reverse(), BB, tag::divided_in, 15 );
	Mesh BBCC ( tag::segment, BB .reverse(), CC, tag::divided_in, 15 );
	Mesh CCDD ( tag::segment, CC .reverse(), DD, tag::divided_in, 15 );
	Mesh DDAA ( tag::segment, DD .reverse(), AA, tag::divided_in, 15 );
	Mesh AAEE ( tag::segment, AA .reverse(), EE, tag::divided_in, 15 );
	Mesh BBFF ( tag::segment, BB .reverse(), FF, tag::divided_in, 15 );
	Mesh CCGG ( tag::segment, CC .reverse(), GG, tag::divided_in, 15 );
	Mesh DDHH ( tag::segment, DD .reverse(), HH, tag::divided_in, 15 );
	Mesh EEFF ( tag::segment, EE .reverse(), FF, tag::divided_in, 15 );
	Mesh FFGG ( tag::segment, FF .reverse(), GG, tag::divided_in, 15 );
	Mesh GGHH ( tag::segment, GG .reverse(), HH, tag::divided_in, 15 );
	Mesh HHEE ( tag::segment, HH .reverse(), EE, tag::divided_in, 15 );

	Mesh ABCD_s ( tag::quadrangle, AABB, BBCC, CCDD, DDAA );
	Mesh EFGH_s ( tag::quadrangle, EEFF, FFGG, GGHH, HHEE );
	Mesh ADHE_s ( tag::quadrangle, DDAA .reverse(), DDHH, HHEE, AAEE .reverse() );
	Mesh BAEF_s ( tag::quadrangle, AABB .reverse(), AAEE, EEFF, BBFF .reverse() );
	Mesh CBFG_s ( tag::quadrangle, BBCC .reverse(), BBFF, FFGG, CCGG .reverse() );
	Mesh DCGH_s ( tag::quadrangle, CCDD .reverse(), CCGG, GGHH, DDHH .reverse() );

	// back to the torus (forget about the sphere manifold)
	torus_manif .set_as_working_manifold();
	
	Cell A ( tag::vertex );  x (A) = -1.;  y (A) = -1.;  z (A) = -1.;
	Cell B = A, C = A, D = A;
	Cell E ( tag::vertex );  x (E) = -1.;  y (E) = -1.;  z (E) =  1.;
	Cell F = E, G = E, H = E;

	Mesh AB ( tag::segment, A .reverse(), B, tag::divided_in, 15, tag::winding, gy );
	Mesh BC ( tag::segment, B .reverse(), C, tag::divided_in, 15, tag::winding, gx );
	Mesh CD = AB .reverse(), DA = BC .reverse();
	Mesh AE ( tag::segment, A .reverse(), E, tag::divided_in, 15 );
	Mesh BF = AE, CG = AE, DH = AE;
	Mesh EF ( tag::segment, E .reverse(), F, tag::divided_in, 15, tag::winding, gy );
	Mesh FG ( tag::segment, F .reverse(), G, tag::divided_in, 15, tag::winding, gx );
	Mesh GH = EF .reverse(), HE = FG .reverse();

	Mesh ABCD ( tag::quadrangle, AB, BC, CD, DA, tag::winding );
	Mesh EFGH ( tag::quadrangle, EF, FG, GH, HE, tag::winding );
	Mesh ADHE ( tag::quadrangle, DA .reverse(), DH, HE, AE .reverse(), tag::winding );
	Mesh BAEF ( tag::quadrangle, AB .reverse(), AE, EF, BF .reverse(), tag::winding );
	Mesh CBFG ( tag::quadrangle, BC .reverse(), BF, FG, CG .reverse(), tag::winding );
	Mesh DCGH ( tag::quadrangle, CD .reverse(), CG, GH, DH .reverse(), tag::winding );

	Mesh AAA ( tag::segment, A .reverse(), AA, tag::divided_in, 8 );
	Mesh BBB ( tag::segment, B .reverse(), BB, tag::divided_in, 8, tag::winding, -gy );
	Mesh CCC ( tag::segment, C .reverse(), CC, tag::divided_in, 8, tag::winding, -gx-gy );
	Mesh DDD ( tag::segment, D .reverse(), DD, tag::divided_in, 8, tag::winding, -gx );
	Mesh EEE ( tag::segment, E .reverse(), EE, tag::divided_in, 8 );
	Mesh FFF ( tag::segment, F .reverse(), FF, tag::divided_in, 8, tag::winding, -gy );
	Mesh GGG ( tag::segment, G .reverse(), GG, tag::divided_in, 8, tag::winding, -gx-gy );
	Mesh HHH ( tag::segment, H .reverse(), HH, tag::divided_in, 8, tag::winding, -gx );

	Mesh AEEA ( tag::quadrangle, AE, EEE, AAEE .reverse(), AAA .reverse(), tag::winding );
	Mesh BFFB ( tag::quadrangle, BF, FFF, BBFF .reverse(), BBB .reverse(), tag::winding );
	Mesh CGGC ( tag::quadrangle, CG, GGG, CCGG .reverse(), CCC .reverse(), tag::winding );
	Mesh DHHD ( tag::quadrangle, DH, HHH, DDHH .reverse(), DDD .reverse(), tag::winding );
	Mesh ABBA ( tag::quadrangle, AB, BBB, AABB .reverse(), AAA .reverse(), tag::winding );
	Mesh BCCB ( tag::quadrangle, BC, CCC, BBCC .reverse(), BBB .reverse(), tag::winding );
	Mesh CDDC ( tag::quadrangle, CD, DDD, CCDD .reverse(), CCC .reverse(), tag::winding );
	Mesh DAAD ( tag::quadrangle, DA, AAA, DDAA .reverse(), DDD .reverse(), tag::winding );
	Mesh EFFE ( tag::quadrangle, EF, FFF, EEFF .reverse(), EEE .reverse(), tag::winding );
	Mesh FGGF ( tag::quadrangle, FG, GGG, FFGG .reverse(), FFF .reverse(), tag::winding );
	Mesh GHHG ( tag::quadrangle, GH, HHH, GGHH .reverse(), GGG .reverse(), tag::winding );
	Mesh HEEH ( tag::quadrangle, HE, EEE, HHEE .reverse(), HHH .reverse(), tag::winding );

	Mesh cube_ABCD ( tag::cube, ABCD, ABCD_s .reverse(), ABBA .reverse(), CDDC .reverse(),
	                            DAAD .reverse(), BCCB .reverse(), tag::winding            );
	Mesh cube_EFGH ( tag::cube, EFGH, EFGH_s .reverse(), EFFE .reverse(), GHHG .reverse(),
	                            FGGF .reverse(), HEEH .reverse(), tag::winding            );
	Mesh cube_BAEF ( tag::cube, BAEF, BAEF_s .reverse(), ABBA, EFFE .reverse(),
	                            AEEA .reverse(), BFFB, tag::winding            );
	Mesh cube_CBFG ( tag::cube, CBFG, CBFG_s .reverse(), BCCB, FGGF .reverse(),
	                            BFFB .reverse(), CGGC, tag::winding            );
	Mesh cube_DCGH ( tag::cube, DCGH, DCGH_s .reverse(), CDDC, GHHG .reverse(),
	                            CGGC .reverse(), DHHD, tag::winding            );
	Mesh cube_ADHE ( tag::cube, ADHE, ADHE_s .reverse(), DAAD, HEEH .reverse(),
	                            DHHD .reverse(), AEEA, tag::winding            );

	std::cout << "built six 3D meshes, now joining" << std::endl;
	Mesh torus ( tag::join, { cube_ABCD, cube_EFGH .reverse(), cube_BAEF,
	                          cube_CBFG, cube_DCGH, cube_ADHE            } );

	std::cout << "built 3D torus with round hole, total " << torus .number_of ( tag::cells_of_max_dim )
	          << " cubes, now unfolding" << std::endl;
	Mesh torus_unfolded = torus .unfold ( tag::over_region, x*x + 2*y*y + 3*z*z < 10 );

	std::cout << "built unfolded mesh with " << torus_unfolded .number_of ( tag::cells_of_max_dim )
	          << " cubes, now exporting" << std::endl;
	torus_unfolded .export_to_file ( tag::msh, "torus_of_cubes.msh");
	
	std::cout << "produced file torus_of_cubes.msh" << std::endl;

}  // end of main


