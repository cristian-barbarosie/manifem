
#include "maniFEM.h"
#include <set>
#include <fstream>

using namespace maniFEM;


bool opposite_signs ( double a, double b )
{	return ( ( a >= 0. ) and ( b < 0. ) ) or ( ( a < 0. ) and ( b >= 0. ) );  }

void special_draw ( Mesh & square, Mesh & cut, std::string f );

Mesh build_interface ( Mesh ambient, Function psi );



int main ()
	
{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0],  y = xy[1];

	Cell A ( tag::vertex );  x(A) = -1.;  y(A) = 0.;
	Cell B ( tag::vertex );  x(B) =  1.;  y(B) = 0.;
	Cell C ( tag::vertex );  x(C) =  1.;  y(C) = 1.;
	Cell D ( tag::vertex );  x(D) = -1.;  y(D) = 1.;

	Mesh AB ( tag::segment, A.reverse(), B, tag::divided_in, 24 );
	Mesh BC ( tag::segment, B.reverse(), C, tag::divided_in, 12 );
	Mesh CD ( tag::segment, C.reverse(), D, tag::divided_in, 24 );
	Mesh DA ( tag::segment, D.reverse(), A, tag::divided_in, 12 );

	Mesh rect_mesh ( tag::rectangle, AB, BC, CD, DA );

	double radius = 0.35;
	Function psi = 0.5 * ( ( x*x + (y-0.2)*(y-0.2) ) / radius - radius );

	// 0.5 * ( ( x*x + (y-0.2)*(y-0.2) ) / radius - radius )  circulo

	// x*x + (y-0.2)*(y-0.2) - 0.3 + 0.2*x*y - 1.35*x*x*y*y
	// x*x + (y-0.2)*(y-0.2) - 0.3 + 0.2*x*y - 1.5*x*x*y*y

	Mesh cut = build_interface ( rect_mesh, psi );
	// just finds the segments where the level line  psi == 0  passes

	// improve_angles ( rect_mesh, psi, cut );
	// splits some rectangles in two triangles
	
	special_draw ( rect_mesh, cut, "square-cut.eps" );

	std::cout << "reached end" << std::endl;
	
} // end of main
	
//-----------------------------------------------------------------------------------//


Mesh build_interface ( Mesh ambient, Function psi )

{	// we use the current manifold
	Manifold space = Manifold::working;
	assert ( space.exists() );
	Function coord = space.coordinates();
	assert ( coord.nb_of_components() == 2 );
	// we split 'coord' into its components
	Function x = coord[0],  y = coord[1];

	Function psi_x = psi.deriv(x), psi_y = psi.deriv(y);

	Mesh interf ( tag::fuzzy, tag::of_dimension, 1 );
	// empty mesh, it will be returned after adding segments to it

	// we run over all segments of 'ambient' mesh
	// looking for those where psi changes sign

	{ // just a block of code for hiding 'it'
	CellIterator it = ambient.iterator ( tag::over_segments );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell seg = *it;
		Cell A = seg.base().reverse();
		Cell B = seg.tip();
		if ( not opposite_signs ( psi(A), psi(B) ) ) continue;
		// now psi takes opposite signs on the extremities of 'seg'
		// we want to add 'seg' to 'interf'
		// but we must choose between 'seg' and 'seg.reverse()'
		// we choose the orientation such that the gradient of psi
		// points to the right hand side
		double dx = x(B) - x(A);
		double dy = y(B) - y(A);
		double gx = psi_x(A);
		double gy = psi_y(A);
		if ( gx*dy > gy*dx ) seg.add_to_mesh ( interf );
		else seg.reverse().add_to_mesh ( interf );                    }
	} // just a block of code for hiding 'it'

	// at this moment we have a mesh 'interf' but it is highly disconnected
	// if you want to see it, uncomment next line
	// return interf;
	
	// in order to close the curve,
	// we must add some other segments, although psi does not change sign there

	// we look for segments where 'interf' is disconnected
	// we keep a list (set) of problematic squares

	std::set < Cell > set_of_squares;
	{ // just a block of code for hiding 'it'
	CellIterator it = interf.iterator ( tag::over_segments );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell seg = *it;
		Cell A = seg.base().reverse();
		Cell B = seg.tip();
		if ( interf.cell_behind(A,tag::may_not_exist).exists()
	       and interf.cell_in_front_of(B,tag::may_not_exist).exists() )
			continue;  // not problematic
		// now, 'interf' is disconnected at this 'seg' - we keep both neighbour squares
		Cell sq = ambient.cell_behind ( seg, tag::may_not_exist );
		if ( sq.exists() ) set_of_squares.insert ( sq );
		sq = ambient.cell_in_front_of ( seg, tag::may_not_exist );
		if ( sq.exists() ) set_of_squares.insert ( sq );                                }
	} // just a block of code for hiding 'it'

	std::cout << "we found " << set_of_squares.size() << " squares, ";
	
	while ( true )
	{	std::forward_list<Cell> to_eliminate;
		std::set<Cell>::iterator it_set;
		for ( it_set = set_of_squares.begin(); it_set != set_of_squares.end(); it_set++ )
		{	Cell sq = *it_set;
			bool elim = false;
			CellIterator it = sq.boundary().iterator ( tag::over_vertices );
			for ( it.reset(); it.in_range(); it++ )
			{	Cell P = *it;
				elim = elim or P.belongs_to ( interf ); }
			if ( elim )
			{	to_eliminate.push_front ( sq );  continue;  }                           }
		std::forward_list<Cell>::iterator it_list;
		//  can we use a list of iterators over set_of_squares ?
		for ( it_list = to_eliminate.begin(); it_list != to eliminate.end(); it_list++ )
		{	Cell sq = *it_set;  set_of_squares.erase ( sq );  }
		break;                                                                           }

	std::cout << "but only " << set_of_squares.size() << " are problematic" << std::endl;
	
	return interf;                                                       }

//-----------------------------------------------------------------------------------//

	
void special_draw ( Mesh & square, Mesh & cut, std::string f )

{	// we use the current manifold
	Manifold space = Manifold::working;
	assert ( space.exists() );
	Function coord = space.coordinates();
	assert ( coord.nb_of_components() == 2 );
	// we split 'coord' into its components
	Function x = coord[0],  y = coord[1];

	double xmin, xmax, ymin, ymax, maxside;
	
	{ // just a block for hiding variables
	CellIterator it = square.iterator ( tag::over_vertices );
	it.reset(); assert ( it.in_range() );
	Cell Vfirst = *it;
	xmin = xmax = x(Vfirst);
	ymin = ymax = y(Vfirst);	
	for ( it++; it.in_range(); it++ )
	{ Cell V = *it; 
		double xV = x(V), yV = y(V);
		if ( xV < xmin ) xmin = xV;
	  if ( xV > xmax ) xmax = xV;
	  if ( yV < ymin ) ymin = yV;
	  if ( yV > ymax ) ymax = yV;      }
	} // just a block for hiding variables 
	if ( xmax-xmin < ymax-ymin ) maxside = ymax-ymin;
	else maxside = xmax-xmin; 
	double border = 0.02*maxside;
	double scale_factor = 500/maxside;
	double translation_x = -scale_factor*xmin;
	double translation_y = -scale_factor*ymin;
	
	std::ofstream file_ps (f);
	file_ps << "%!PS-Adobe-3.0 EPSF-3.0" << std::endl;
	file_ps << "%%Title:                     malha" << std::endl;
	file_ps << "%%BoundingBox:  " << -border*scale_factor << " " << -border*scale_factor << " "
					<< scale_factor*(xmax+border) + translation_x << "   "
					<< scale_factor*(ymax+border) + translation_y << std::endl;
	file_ps << "%%EndComments" << std::endl << std::endl;

	file_ps << "gsave" << std::endl;
	// file_ps << "/m{moveto}def" << std::endl;
	// file_ps << "/l{lineto}def" << std::endl;
	// file_ps << "/s{stroke}def" << std::endl;
	file_ps << translation_x << " " << translation_y << " translate" << std::endl;
	file_ps << scale_factor << " " << scale_factor << " scale" << std::endl << std::endl;
	
	{ // just a block of code for hiding 'it'
	file_ps << "gsave 0.004 setlinewidth 0.5 setgray" << std::endl;
	CellIterator it = square.iterator ( tag::over_segments );
	for ( it.reset(); it.in_range(); it++)
	{	Cell seg = *it;
		Cell base = seg.base().reverse();
		Cell tip  = seg.tip();
		file_ps << x(base) << " " << y(base) << " moveto" << std::endl;
		file_ps << x(tip) << " " << y(tip) << " lineto stroke" << std::endl;  }
	} // just a block of code for hiding 'it'
	
#ifndef NDEBUG
	{ // just a block of code for hiding 'it'
	CellIterator it = square.iterator ( tag::over_vertices );
	file_ps << "/Courier findfont 0.2 scalefont setfont" << std::endl;
	for ( it.reset(); it.in_range(); it++ )
	{	Cell p = *it;
		if ( p.core->name == "" ) continue;
		file_ps << x(p) << " " << y(p) << " moveto (" << p.core->name << ") show" << std::endl;  }
	} // just a block of code for hiding 'it'
#endif

	std::ifstream file_arrows ("arrows.ps");
	while ( true )
	{	std::string line;  // this way 'line' remains local
		if ( getline ( file_arrows, line ) ) file_ps << line + '\n';
		else break;                                               }

	{ // just a block of code for hiding 'it'
	file_ps << "0.8 0 0 setrgbcolor 0.007 setlinewidth" << std::endl;
	CellIterator it = cut.iterator ( tag::over_segments );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell seg = *it;
		Cell base = seg.base().reverse();
		Cell tip  = seg.tip();
		file_ps << x(base) << " " << y(base) << " moveto" << std::endl;
		file_ps << x(tip) << " " << y(tip) << " Lineto^ stroke" << std::endl;  }		
	} // just a block of code for hiding 'it'
	
	file_ps << "grestore" << std::endl;
	file_ps << "grestore" << std::endl;
	file_ps << "showpage" << std::endl;
	file_ps << "%%Trailer" << std::endl;
	file_ps << "%EOF" << std::endl;

	if ( ! file_ps.good() )
	{	std::cerr << "error writing postscript file" << std::endl;
		exit (1);                                                  }

} // end of special_draw



