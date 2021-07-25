
#include "maniFEM.h"
#include <set>
#include <fstream>

using namespace maniFEM;

int global_counter = 0;

inline bool opposite_signs ( double a, double b )
{	return ( ( a >= 0. ) and ( b <= 0. ) ) or ( ( a <= 0. ) and ( b >= 0. ) );  }

void special_draw ( Mesh & square, Mesh & cut, std::string f );

Mesh build_interface ( Mesh ambient, Function psi );

bool on_boundary ( Cell A, Mesh msh );

//-----------------------------------------------------------------------------------//


int main ()
	
{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0],  y = xy[1];

	Cell A ( tag::vertex );  x(A) = -1.;  y(A) = 0.;
	Cell B ( tag::vertex );  x(B) =  1.;  y(B) = 0.;
	Cell C ( tag::vertex );  x(C) =  1.;  y(C) = 1.;
	Cell D ( tag::vertex );  x(D) = -1.;  y(D) = 1.;

	Mesh AB ( tag::segment, A.reverse(), B, tag::divided_in, 40 );
	Mesh BC ( tag::segment, B.reverse(), C, tag::divided_in, 20 );
	Mesh CD ( tag::segment, C.reverse(), D, tag::divided_in, 40 );
	Mesh DA ( tag::segment, D.reverse(), A, tag::divided_in, 20 );

	Mesh rect_mesh ( tag::rectangle, AB, BC, CD, DA );

//	double radius = 0.4;
	// std::cout << "radius = ";  std::cin >> radius;
	Function psi = x*x + (y-0.2)*(y-0.2) - 0.3 + 0.2*x*y - 1.5*x*x*y*y;

	// 0.5 * ( ( x*x + (y-0.2)*(y-0.2) ) / radius - radius )  circulo

	// x*x + (y-0.2)*(y-0.2) - 0.3 + 0.2*x*y - 1.35*x*x*y*y
	// x*x + (y-0.2)*(y-0.2) - 0.3 + 0.2*x*y - 1.5*x*x*y*y

	Mesh cut = build_interface ( rect_mesh, psi );
	// builds a set of vertices next to the level line  psi == 0
	// then builds a one-dimensional mesh

	Manifold level = RR2.implicit ( psi == 0. );

	CellIterator it=cut.iterator(tag::over_vertices);
	for(it.reset(); it.in_range(); it++) {
		Cell V=*it;
		level.project(V);
		if (V.belongs_to(AB)) y(V)=0.;
		if (V.belongs_to(BC)) x(V)=1.;
		if (V.belongs_to(CD)) y(V)=1.;
		if (V.belongs_to(DA)) x(V)=-1.;
	}

	// 'level' is working manifold
	
	for ( it.reset(); it.in_range(); it++ )  // 'it' runs over vertices of cut 
	{	Cell W = *it; cut.baricenter ( W );  }

	RR2.set_as_working_manifold();

	for ( it.reset(); it.in_range(); it++)  // 'it' runs over vertices of cut 
	{	Cell V = *it;
		CellIterator itv = rect_mesh.iterator ( tag::over_vertices, tag::around, V ); 
		for ( itv.reset(); itv.in_range(); itv++ )
		{	Cell W = *itv;
			if ( W.belongs_to(cut) ) continue;
			if ( W.belongs_to(AB) )  { AB.baricenter(W); continue; }
			if ( W.belongs_to(BC) )  { BC.baricenter(W); continue; }
			if ( W.belongs_to(CD) )  { CD.baricenter(W); continue; }
			if ( W.belongs_to(DA) )  { DA.baricenter(W); continue; }
			rect_mesh.baricenter ( W );                              }                 }

	special_draw ( rect_mesh, cut, "square-cut.eps" );
	rect_mesh.export_msh("background.msh");
	std::cout << "produced files square-cut.eps and background.msh" << std::endl;
	exit ( 0 );
	
} // end of main
	
//-----------------------------------------------------------------------------------//


class compare_values_of

{ public :
	Function f;
	inline compare_values_of ( const Function ff ) : f { ff } { }
	inline bool operator() ( Cell A, Cell B ) const
	{	return this->f ( A ) > this->f ( B );  }                     };

//-----------------------------------------------------------------------------------//


bool on_boundary ( Cell A, Mesh msh )

{	assert ( A.dim() == 0 );
	assert ( A.is_positive() );

	CellIterator it = msh.iterator ( tag::over_segments, tag::around, A );
	it.reset();
	assert ( it.in_range() );
	return not msh.cell_behind ( *it, tag::may_not_exist ) .exists();         }
	
//-----------------------------------------------------------------------------------//

void verify_tri ( Cell tri )

{	assert ( tri.dim() == 2 );
	assert ( tri.boundary().number_of ( tag::vertices ) == 3 );

	CellIterator it = tri.boundary().iterator ( tag::over_vertices );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell V = *it;
		Cell::Positive * cll = tag::Util::assert_cast
			< Cell::Core*, Cell::Positive* > ( V.core );
		typedef std::map < Mesh::Core*, Cell::field_to_meshes > maptype;
		// typedef std::map < Mesh::Core*, Cell::field_to_meshes_same_dim > maptype_sd;
		maptype & cm2 = cll->meshes[1];
		typename maptype::iterator itt = cm2.find ( tri.boundary().core );
		assert ( itt != cm2.end() );
		Cell::field_to_meshes f = itt->second;
		std::cout << "f, neg " << f.counter_neg << ", pos " << f.counter_pos << std::endl;
	}		
}		

//-----------------------------------------------------------------------------------//

inline Cell find_segment ( Cell & A, Cell & B, Mesh & ambient )

{	assert ( A.is_positive() );
	assert ( B.is_positive() );
	assert ( A.dim() == 0 );
	assert ( B.dim() == 0 );
	CellIterator it_seg = ambient.iterator ( tag::over_segments, tag::around, A );
	// how about tag::around, A.reverse() ?
	for ( it_seg.reset(); it_seg.in_range(); it_seg++ )
	{	Cell seg = *it_seg;
		assert ( seg.tip() == A );
		if ( seg.base().reverse() == B ) return seg.reverse();      }
	// seg not found, A and B must be diagonally opposed in some square
	CellIterator it_sq = ambient.iterator ( tag::over_cells_of_dim, 2, tag::around, A );
	for ( it_sq.reset(); it_sq.in_range(); it_sq++ )
	{	Cell sq = *it_sq;
		assert ( A.belongs_to ( sq.boundary() ) );
		if ( not B.belongs_to ( sq.boundary() ) ) continue;
		Cell AC = sq.boundary().cell_in_front_of ( A );
		Cell C = AC.tip();
		Cell CB = sq.boundary().cell_in_front_of ( C );
		assert ( CB.tip() == B );
		Cell BD = sq.boundary().cell_in_front_of ( B );
		Cell D = BD.tip();
		Cell DA = sq.boundary().cell_in_front_of ( D );
		assert ( DA.tip() == A );
		sq.remove_from_mesh ( ambient );
		Cell AB ( tag::segment, A.reverse(), B );
		Cell ACB ( tag::triangle, AC, CB, AB.reverse() );
		Cell BDA ( tag::triangle, BD, DA, AB );
		// verify_tri ( ACB );  verify_tri ( BDA );
		ACB.add_to_mesh ( ambient );
		BDA.add_to_mesh ( ambient );
		return AB;                                            }
	assert ( false );

}  // end of  find_segment

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

	// we run over all segments of 'ambient' mesh
	// looking for those where psi changes sign
	// and build a set of vertices, ordered by the absolute values of psi

	compare_values_of comp_psi ( abs(psi) );
	std::multiset < Cell, compare_values_of > ms ( comp_psi );
	
	{  // just a block of code for hiding names
	CellIterator it_seg = ambient.iterator ( tag::over_segments );
	for ( it_seg.reset(); it_seg.in_range(); it_seg++ )
	{	Cell seg = *it_seg;
		Cell A = seg.base().reverse();
		Cell B = seg.tip();
		if ( not opposite_signs ( psi(A), psi(B) ) ) continue;
		ms.insert ( A );
		ms.insert ( B );                                        }
	}  // just a block of code for hiding names

	// at this moment we have a cloud of points, possibly repeated

	Mesh interf ( tag::fuzzy, tag::of_dimension, 1 );
	// empty mesh, it will be returned after adding segments to it

	{  // just a block of code for hiding names
	std::set < Cell > set_useful, set_needed;
	std::forward_list < Cell > list_useful;
	std::multiset<Cell,compare_values_of>::iterator it_ms;

	for ( it_ms = ms.begin(); it_ms != ms.end(); it_ms++ )
	{	Cell A = * it_ms;
		if ( set_useful.find ( A ) == set_useful.end() )
		{	set_useful.insert ( A );
			list_useful.push_front ( A );                   }  }

	std::forward_list<Cell>::iterator it_list;
	for ( it_list = list_useful.begin(); it_list != list_useful.end(); it_list++ )
	
	{	// this process begins with values of psi close to zero
		Cell A = *it_list;
		if ( set_useful.find ( A ) == set_useful.end() ) continue;
		
		double xA = x(A), yA = y(A);
		double gx = psi_x(A), gy = psi_y(A);

		CellIterator it_around = ambient.iterator
			( tag::over_cells_of_dim, 2, tag::around, A );
		std::set < Cell > set_of_neigh_1, set_of_neigh_2;
		for ( it_around.reset(); it_around.in_range(); it_around++ )
		{	Cell cll = *it_around;
			CellIterator it_ver = cll.boundary().iterator ( tag::over_vertices );
			for ( it_ver.reset(); it_ver.in_range(); it_ver++ )
			{	Cell B = *it_ver;
				if ( set_useful.find(B) == set_useful.end() ) continue;
				if ( B == A ) continue;
				double dx = x(B) - xA, dy = y(B) - yA;
				if ( dx*gy > dy*gx ) set_of_neigh_1.insert ( B );
			  else set_of_neigh_2.insert ( B );                 }     }
		std::set<Cell>::iterator it_neigh;
		int found = 0;
		Cell B_found ( tag::non_existent );
		for ( it_neigh = set_of_neigh_1.begin(); it_neigh != set_of_neigh_1.end(); it_neigh++ )
		{	Cell B = *it_neigh;
			if ( set_needed.find ( B ) != set_needed.end() )
			{	found++; B_found = B;  }                        }
		assert ( ( found == 0 ) or ( found == 1 ) );
		if ( found )
		{	for ( it_neigh = set_of_neigh_1.begin(); it_neigh != set_of_neigh_1.end(); it_neigh++ )
			{	Cell B = *it_neigh;
				if ( B == B_found ) continue;
				set_useful.erase ( B );        }  }
		else  // not found
		{	it_neigh = set_of_neigh_1.begin();
			if ( it_neigh != set_of_neigh_1.end() )
			{	Cell Bmin = *it_neigh;
				double psi_min = std::abs ( psi ( Bmin ) );
				for ( it_neigh++; it_neigh != set_of_neigh_1.end(); it_neigh++ )
				{	Cell B = *it_neigh;
					double g = std::abs ( psi ( B ) );
					if ( g < psi_min ) { psi_min = g; Bmin = B;  }     }
				for ( it_neigh = set_of_neigh_1.begin(); it_neigh != set_of_neigh_1.end(); it_neigh++ )
				{	Cell B = *it_neigh;
					if ( B != Bmin )
					if ( B != Bmin ) set_useful.erase ( B );  }
				if ( on_boundary ( A, ambient ) and on_boundary ( Bmin, ambient ) )
					set_useful.erase ( Bmin );
				else
				{	Cell AB ( tag::segment, A.reverse(), Bmin );
					set_needed.insert ( Bmin );                  }           }  }
		found = 0;
		B_found = Cell ( tag::non_existent );
		for ( it_neigh = set_of_neigh_2.begin(); it_neigh != set_of_neigh_2.end(); it_neigh++ )
		{	Cell B = *it_neigh;
			if ( set_needed.find ( B ) != set_needed.end() )
			{	found++; B_found = B;  }                        }
		assert ( ( found == 0 ) or ( found == 1 ) );
		if ( found )
		{	for ( it_neigh = set_of_neigh_2.begin(); it_neigh != set_of_neigh_2.end(); it_neigh++ )
			{	Cell B = *it_neigh;
				if ( B == B_found ) continue;
				set_useful.erase ( B );        }  }
		else  // not found
		{	it_neigh = set_of_neigh_2.begin();
			if ( it_neigh != set_of_neigh_2.end() )
			{	Cell Bmin = *it_neigh;
				double psi_min = std::abs ( psi ( Bmin ) );
				for ( it_neigh++; it_neigh != set_of_neigh_2.end(); it_neigh++ )
				{	Cell B = *it_neigh;
					double g = std::abs ( psi ( B ) );
					if ( g < psi_min ) { psi_min = g; Bmin = B;  }     }
				for ( it_neigh = set_of_neigh_2.begin(); it_neigh != set_of_neigh_2.end(); it_neigh++ )
				{	Cell B = *it_neigh;
					if ( B!= Bmin )
					if ( B != Bmin ) set_useful.erase ( B );  }
				if ( on_boundary ( A, ambient ) and on_boundary ( Bmin, ambient ) )
					set_useful.erase ( Bmin );
				else
				{	Cell BA ( tag::segment, Bmin.reverse(), A );
					set_needed.insert ( Bmin );                  }           }   }

	}  // end of for it_list
		
	std::set<Cell>::iterator it_set;
	for ( it_set = set_useful.begin(); it_set != set_useful.end(); it_set++ )

	{	Cell A = *it_set;

		bool front_exists = interf.cell_in_front_of ( A, tag::may_not_exist ) .exists();
		bool back_exists = interf.cell_behind ( A, tag::may_not_exist ) .exists();

		if ( front_exists and back_exists ) continue;
		
		double xA = x(A), yA = y(A);
		double gx = psi_x(A), gy = psi_y(A);

		CellIterator it_around = ambient.iterator
			( tag::over_cells_of_dim, 2, tag::around, A );
		std::set < Cell > set_of_neigh;
		for ( it_around.reset(); it_around.in_range(); it_around++ )
		{	Cell cll = *it_around;
			CellIterator it_ver = cll.boundary().iterator ( tag::over_vertices );
			for ( it_ver.reset(); it_ver.in_range(); it_ver++ )
			{	Cell B = *it_ver;
				if ( B == A ) continue;
				if ( set_useful.find(B) == set_useful.end() ) continue;
				set_of_neigh.insert ( B );                       }                   }
		std::set<Cell>::iterator it_neigh;
		for ( it_neigh = set_of_neigh.begin(); it_neigh != set_of_neigh.end(); it_neigh++ )
		{	Cell B = *it_neigh;
			double dx = x(B) - xA, dy = y(B) - yA;
			// we add a new segment to 'interf'
			if ( ( not front_exists ) and ( dx*gy < dy*gx ) )
			{	Cell seg = find_segment ( A, B, ambient );
				seg.add_to_mesh ( interf );                       }
			if ( ( not back_exists ) and ( dx*gy > dy*gx ) )
			{	Cell seg = find_segment ( B, A, ambient );
				seg.add_to_mesh ( interf );                       }  }

	}  // end of for it_set in set_useful
	}  // just a block of code for hiding names
                  
	// we cut in halves some more squares, at places like :      |
	{  // just a block of code for hiding names                 /
	size_t counter = 0;                //                      |
	CellIterator it_seg = interf.iterator ( tag::over_segments );
	for ( it_seg.reset(); it_seg.in_range(); it_seg++ )
	{	Cell seg = *it_seg;
		Cell A = seg.base().reverse(), B = seg.tip();
		Cell cll1 = ambient.cell_in_front_of ( seg, tag::surely_exists );
		if ( cll1.boundary().number_of ( tag::segments ) != 3 ) continue;
		Cell cll2 = ambient.cell_behind ( seg, tag::surely_exists );
		assert ( cll2.boundary().number_of ( tag::segments ) == 3 );
		Cell front = interf.cell_in_front_of ( B, tag::may_not_exist );
		if ( not front.exists() ) continue;
		Cell back = interf.cell_behind ( A, tag::may_not_exist );
		if ( not back.exists() ) continue;
		Cell cll3 = ambient.cell_in_front_of ( front, tag::surely_exists );
		if ( cll3.boundary().number_of ( tag::segments ) != 4 ) continue;
		Cell cll4 = ambient.cell_behind ( front, tag::surely_exists );
		assert ( cll4.boundary().number_of ( tag::segments ) == 4 );
		Cell cll5 = ambient.cell_in_front_of ( back, tag::surely_exists );
		if ( cll5.boundary().number_of ( tag::segments ) != 4 ) continue;
		Cell cll6 = ambient.cell_behind ( back, tag::surely_exists );
		assert ( cll6.boundary().number_of ( tag::segments ) == 4 );
		size_t neigh_1 = 0, neigh_2 = 0, neigh_3 = 0, neigh_4 = 0;
		CellIterator it1 = ambient.iterator ( tag::over_segments, tag::around, A.reverse() );
		std::cout << "aha 1" << std::endl;
		it1.reset ( tag::start_at, back.reverse() );
		std::cout << "aha 2" << std::endl;
		for ( it1.reset ( tag::start_at, back.reverse() ); it1.in_range(); it1++ )
		{	if ( *it1 == seg ) break;
			neigh_1++;                }
		assert ( ( neigh_1 == 2 ) or ( neigh_1 == 3 ) );
		for ( it1.reset ( tag::start_at, seg ); it1.in_range(); it1++ )
		{	if ( *it1 == back.reverse() ) break;
			neigh_2++;                           }
		assert ( ( neigh_2 == 2 ) or ( neigh_2 == 3 ) );
		if ( neigh_1 + neigh_2 == 6 ) continue;
		assert ( neigh_1 + neigh_2 == 5 );
		CellIterator it2 = ambient.iterator ( tag::over_segments, tag::around, B.reverse() );
		for ( it2.reset ( tag::start_at, seg.reverse() ); it2.in_range(); it2++ )
		{	if ( *it2 == front ) break;
			neigh_3++;                  }
		assert ( ( neigh_3 == 2 ) or ( neigh_3 == 3 ) );
		for ( it2.reset ( tag::start_at, front ); it2.in_range(); it2++ )
		{	if ( *it2 == seg.reverse() ) break;
			neigh_4++;                          }
		assert ( ( neigh_4 == 2 ) or ( neigh_4 == 3 ) );
		if ( neigh_3 + neigh_4 == 6 ) continue;
		assert ( neigh_3 + neigh_4 == 5 );
		
		counter++;                                                                    }
	std::cout << "found " << counter << " interesting configurations" << std::endl;
	}  // just a block of code for hiding names

	return interf;

}  // end of  build_interface

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
	file_ps << "%%Title:                " << f << std::endl;
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
	file_ps << "0.8 0 0 setrgbcolor 0.005 setlinewidth" << std::endl;
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



