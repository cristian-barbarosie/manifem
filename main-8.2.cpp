
// example presented in paragraph 8.2 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// builds a ring-shaped mesh

#include "maniFEM.h"
using namespace maniFEM;


int main ( )

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0],  y = xy[1];

	short int n_sectors = 15;
	double step_theta = 8*atan(1.)/n_sectors;
	short int radial_divisions = 10;
	short int rot_divisions = 5;

	// start the process by building a segment
	Cell ini_A ( tag::vertex );  x(ini_A) = 1.;  y(ini_A) = 0.;
	Cell ini_B ( tag::vertex );  x(ini_B) = 2.;  y(ini_B) = 0.;
	Mesh ini_seg ( tag::segment, ini_A.reverse(), ini_B, tag::divided_in, radial_divisions );
	Mesh prev_seg = ini_seg;
	Cell  A = ini_A,  B = ini_B;
	std::list < Mesh > sectors;

	for ( short int i = 1; i < n_sectors; i++ )
	{	double theta = i * step_theta;
		// we build two new points
		Cell C ( tag::vertex );  x(C) =    cos(theta);  y(C) =    sin(theta);
		Cell D ( tag::vertex );  x(D) = 2.*cos(theta);  y(D) = 2.*sin(theta);
		// and three new segments
		Mesh BD ( tag::segment, B.reverse(), D, tag::divided_in, rot_divisions );
		Mesh DC ( tag::segment, D.reverse(), C, tag::divided_in, radial_divisions );
		Mesh CA ( tag::segment, C.reverse(), A, tag::divided_in, rot_divisions );
		// and a quadrangle
		Mesh quadr ( tag::quadrangle, prev_seg, BD, DC, CA );
		sectors.push_back ( quadr );
		prev_seg = DC.reverse();
		A = C;  B = D;                                                                }

	// we now build the last sector, thus closing the ring
	// prev_seg, A and B have rotated during the construction process
	// but ini_seg, ini_A and ini_B are the same, initial, ones
	Mesh outer ( tag::segment, B.reverse(), ini_B, tag::divided_in, rot_divisions );
	Mesh inner ( tag::segment, ini_A.reverse(), A, tag::divided_in, rot_divisions );
	Mesh quadr ( tag::quadrangle, outer, ini_seg.reverse(), inner, prev_seg );
	sectors.push_back ( quadr );

	Mesh ring ( tag::join, sectors );

	ring.export_msh ("ring.msh");
	std::cout << "produced file ring.msh" << std::endl;
}
