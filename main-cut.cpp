
#include "maniFEM.h"
#include <set>
#include <fstream>

using namespace maniFEM;


bool opposite_signs ( double a, double b )
{	return ( ( a >= 0. ) and ( b < 0. ) ) or ( ( a < 0. ) and ( b >= 0. ) );  }

void special_draw ( Mesh & square, Mesh & cut, std::string f );

Mesh build_interface ( Mesh ambient, Function psi );


//-----------------------------------------------------------------------------------//


class compare_values_of

{ public :
	Function f;
	inline compare_values_of ( const Function ff ) : f { ff } { }
	inline bool operator() ( Cell A, Cell B ) const
	{	return this->f(A) < this->f(B);  }                          };

//-----------------------------------------------------------------------------------//


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

	double radius = 0.41;
	// std::cout << "radius = ";  std::cin >> radius;
	Function psi = 0.5*x*x + (y-0.35)*(y-0.35) - 0.15 ;

	// 0.5 * ( ( x*x + (y-0.2)*(y-0.2) ) / radius - radius )  circulo

	// x*x + (y-0.2)*(y-0.2) - 0.3 + 0.2*x*y - 1.35*x*x*y*y
	// x*x + (y-0.2)*(y-0.2) - 0.3 + 0.2*x*y - 1.5*x*x*y*y

	compare_values_of comp_psi ( -abs(psi) );
	std::multiset < Cell, compare_values_of > m ( comp_psi );

//	CellIterator it = rect_mesh.iterator ( tag::over_vertices );
//	for ( it.reset(); it.in_range(); it++ ) m.insert ( *it );
//	std::multiset<Cell,compare_values_of>::iterator itt;
//	for ( itt = m.begin(); itt != m.end(); itt++ )
//		std::cout << psi(*itt) << std::endl;

	Mesh cut = build_interface ( rect_mesh, psi );
	// just finds the segments where the level line  psi == 0  passes

	// improve_angles ( rect_mesh, psi, cut );
	// splits some rectangles in two triangles
	
	special_draw ( rect_mesh, cut, "square-cut.eps" );

	std::cout << "produced file square-cut.eps" << std::endl;
	
} // end of main
	
//-----------------------------------------------------------------------------------//


void join_two_opposite_segs_pos
( Cell square, Cell seg1, Cell seg2, Mesh & ambient, Mesh & interf,
  const Function psi, const Function dpsi_dx, const Function dpsi_dy,
  std::set<Cell> & set_of_squares                                    )

{	assert ( seg1.belongs_to ( square.boundary(), tag::oriented ) );
	assert ( seg2.belongs_to ( square.boundary(), tag::oriented ) );
	assert ( seg1.belongs_to ( interf, tag::oriented ) );
	assert ( seg2.belongs_to ( interf, tag::oriented ) );
	
	bool tip_of_seg1_free = not interf.cell_in_front_of
	       ( seg1.tip(), tag::may_not_exist ) .exists(),
	     base_of_seg1_free = not interf.cell_behind
	       ( seg1.base().reverse(), tag::may_not_exist ) .exists(),
	     tip_of_seg2_free = not interf.cell_in_front_of
	       ( seg2.tip(), tag::may_not_exist ) .exists(),
	     base_of_seg2_free = not interf.cell_behind
	       ( seg2.base().reverse(), tag::may_not_exist ) .exists();

	Cell other = ambient.cell_in_front_of ( seg1, tag::may_not_exist );
	if ( other.exists() ) set_of_squares.insert ( other );
	other = ambient.cell_in_front_of ( seg2, tag::may_not_exist );
	if ( other.exists() ) set_of_squares.insert ( other );

	if ( tip_of_seg1_free and tip_of_seg2_free and
	     base_of_seg1_free and base_of_seg2_free ) // much freedom
	// but we can choose according to the orientation of grad psi
	{	// we use the current manifold
		Manifold space = Manifold::working;
		assert ( space.exists() );
		Function coord = space.coordinates();
		assert ( coord.nb_of_components() == 2 );
		// we split 'coord' into its components
		Function x = coord[0],  y = coord[1];
		Cell seg = square.boundary().cell_in_front_of ( seg1.tip() );
		double dx = x ( seg.tip() ) - x ( seg.base().reverse() ),
		       dy = y ( seg.tip() ) - y ( seg.base().reverse() );
		double gx = dpsi_dx(seg.tip()), gy = dpsi_dy(seg.tip());
		if ( gx*dy > gy*dx ) seg.add_to_mesh ( interf );
		else
		{	seg = square.boundary().cell_in_front_of ( seg2.tip() );
			seg.add_to_mesh ( interf );                              }
		return;                                                       }

	if ( tip_of_seg1_free and base_of_seg2_free )
	{	Cell seg = square.boundary().cell_in_front_of ( seg1.tip() );
		seg.add_to_mesh ( interf );
		return;                                                      }

	if ( tip_of_seg2_free and base_of_seg1_free )
	{	Cell seg = square.boundary().cell_in_front_of ( seg2.tip() );
		seg.add_to_mesh ( interf );
		return;                                                      }

	assert ( false );

}  // end of  join_two_opposite_segs_pos
		
//-----------------------------------------------------------------------------------//


void join_two_opposite_segs_neg
( Cell square, Cell seg1, Cell seg2, Mesh & ambient, Mesh & interf,
  const Function psi, const Function dpsi_dx, const Function dpsi_dy,
  std::set<Cell> & set_of_squares                                    )

{	assert ( seg1.reverse().belongs_to ( square.boundary(), tag::oriented ) );
	assert ( seg2.reverse().belongs_to ( square.boundary(), tag::oriented ) );
	assert ( seg1.belongs_to ( interf, tag::oriented ) );
	assert ( seg2.belongs_to ( interf, tag::oriented ) );
	
	bool tip_of_seg1_free = not interf.cell_in_front_of
	       ( seg1.tip(), tag::may_not_exist ) .exists(),
	     base_of_seg1_free = not interf.cell_behind
	       ( seg1.base().reverse(), tag::may_not_exist ) .exists(),
	     tip_of_seg2_free = not interf.cell_in_front_of
	       ( seg2.tip(), tag::may_not_exist ) .exists(),
	     base_of_seg2_free = not interf.cell_behind
	       ( seg2.base().reverse(), tag::may_not_exist ) .exists();

	Cell other = ambient.cell_behind ( seg1, tag::may_not_exist );
	if ( other.exists() ) set_of_squares.insert ( other );
	other = ambient.cell_behind ( seg2, tag::may_not_exist );
	if ( other.exists() ) set_of_squares.insert ( other );

	if ( tip_of_seg1_free and tip_of_seg2_free and
	     base_of_seg1_free and base_of_seg2_free ) // much freedom
	// but we can choose according to the orientation of grad psi
	{	// we use the current manifold
		Manifold space = Manifold::working;
		assert ( space.exists() );
		Function coord = space.coordinates();
		assert ( coord.nb_of_components() == 2 );
		// we split 'coord' into its components
		Function x = coord[0],  y = coord[1];
		Cell seg = square.boundary().cell_behind ( seg1.tip() ) .reverse();
		double dx = x ( seg.tip() ) - x ( seg.base().reverse() ),
		       dy = y ( seg.tip() ) - y ( seg.base().reverse() );
		double gx = dpsi_dx(seg.tip()), gy = dpsi_dy(seg.tip());
		if ( gx*dy > gy*dx ) seg.add_to_mesh ( interf );
		else
		{	seg = square.boundary().cell_behind ( seg2.tip() ) .reverse();
			seg.add_to_mesh ( interf );                                   }
		return;                                                           }

	if ( tip_of_seg1_free and base_of_seg2_free )
	{	Cell seg = square.boundary().cell_behind ( seg1.tip() ) .reverse();
		seg.add_to_mesh ( interf );
		return;                                                            }

	if ( tip_of_seg2_free and base_of_seg1_free )
	{	Cell seg = square.boundary().cell_behind ( seg2.tip() ) .reverse();
		seg.add_to_mesh ( interf );
		return;                                                            }

	assert ( false );

}  // end of  join_two_opposite_segs_neg
		
//-----------------------------------------------------------------------------------//


void join_two_parallel_segs
( Cell square, Cell seg1, Cell seg2, Mesh & ambient, Mesh & interf,
  const Function psi, const Function dpsi_dx, const Function dpsi_dy,
  std::set<Cell> & set_of_squares                                    )

// tries to "solve" the problem of a square having two disconnected segments of 'interf'
// returns true is square "solved" (will be eliminated by the calling functions)
	
{	// we use the current manifold
	Manifold space = Manifold::working;
	assert ( space.exists() );
	Function coord = space.coordinates();
	assert ( coord.nb_of_components() == 2 );
	// we split 'coord' into its components
	Function x = coord[0],  y = coord[1];

	assert ( seg1.belongs_to ( square.boundary(), tag::oriented ) );
	assert ( seg2.reverse().belongs_to ( square.boundary(), tag::oriented ) );
	assert ( seg1.belongs_to ( interf, tag::oriented ) );
	assert ( seg2.belongs_to ( interf, tag::oriented ) );
	
	bool tip_of_seg1_free = not interf.cell_in_front_of
	       ( seg1.tip(), tag::may_not_exist ) .exists(),
	     base_of_seg1_free = not interf.cell_behind
	       ( seg1.base().reverse(), tag::may_not_exist ) .exists(),
	     tip_of_seg2_free = not interf.cell_in_front_of
	       ( seg2.tip(), tag::may_not_exist ) .exists(),
	     base_of_seg2_free = not interf.cell_behind
	       ( seg2.base().reverse(), tag::may_not_exist ) .exists(),
	     seg1_on_bdry = not ambient.cell_in_front_of ( seg1, tag::may_not_exist ) .exists(),
	     seg2_on_bdry = not ambient.cell_behind ( seg2, tag::may_not_exist ) .exists();

	if ( tip_of_seg1_free and tip_of_seg2_free and
	     base_of_seg1_free and base_of_seg2_free and
	     not seg1_on_bdry and not seg2_on_bdry       ) // too much freedom
		return;

	if ( ( not base_of_seg1_free ) and tip_of_seg1_free and base_of_seg2_free )
	{	seg1.remove_from_mesh ( interf );
		Cell seg = square.boundary().cell_in_front_of ( seg2.base().reverse() );
		seg.reverse().add_to_mesh ( interf );
		if ( not seg2_on_bdry )
			set_of_squares.insert ( ambient.cell_behind ( seg2 ) );
		return;                                                                 }

	if ( seg1_on_bdry and tip_of_seg1_free and base_of_seg2_free )
	{	Cell seg = square.boundary().cell_in_front_of ( seg2.base().reverse() );
		Cell B = seg.base().reverse(), A = seg.tip();
		double dx = x(B) - x(A);
		double dy = y(B) - y(A);
		double gx = dpsi_dx(A);
		double gy = dpsi_dy(A);
		if ( gx*dy > gy*dx )
		{	seg1.remove_from_mesh ( interf );
			seg.reverse().add_to_mesh ( interf );
			if ( not seg2_on_bdry )
				set_of_squares.insert ( ambient.cell_behind ( seg2 ) );
			return;                                                   }             }

	if ( ( not tip_of_seg1_free ) and base_of_seg1_free and tip_of_seg2_free )
	{	seg1.remove_from_mesh ( interf );
		Cell seg = square.boundary().cell_behind ( seg2.tip() );
		seg.reverse().add_to_mesh ( interf );
		if ( not seg2_on_bdry )
			set_of_squares.insert ( ambient.cell_behind ( seg2 ) );
		return;                                                  }

	if ( seg1_on_bdry and base_of_seg1_free and tip_of_seg2_free )
	{	Cell seg = square.boundary().cell_behind ( seg2.tip() );
		Cell B = seg.base().reverse(), A = seg.tip();
		double dx = x(B) - x(A);
		double dy = y(B) - y(A);
		double gx = dpsi_dx(A);
		double gy = dpsi_dy(A);
		if ( gx*dy > gy*dx )
		{	seg1.remove_from_mesh ( interf );
			seg.reverse().add_to_mesh ( interf );
			if ( not seg2_on_bdry )
				set_of_squares.insert ( ambient.cell_behind ( seg2 ) );
			return;                                                  }     }

	if ( ( not base_of_seg2_free ) and tip_of_seg2_free and base_of_seg1_free )
	{	seg2.remove_from_mesh ( interf );
		Cell seg = square.boundary().cell_in_front_of ( seg2.base().reverse() );
		seg.add_to_mesh ( interf );
		if ( not seg1_on_bdry )
			set_of_squares.insert ( ambient.cell_in_front_of ( seg1 ) );
		return;                                                                 }

	if ( seg2_on_bdry and tip_of_seg2_free and base_of_seg1_free )
	{	Cell seg = square.boundary().cell_in_front_of ( seg2.base().reverse() );
		Cell A = seg.base().reverse(), B = seg.tip();
		double dx = x(B) - x(A);
		double dy = y(B) - y(A);
		double gx = dpsi_dx(A);
		double gy = dpsi_dy(A);
		if ( gx*dy > gy*dx )
		{	seg2.remove_from_mesh ( interf );
			seg.add_to_mesh ( interf );
			if ( not seg1_on_bdry )
				set_of_squares.insert ( ambient.cell_in_front_of ( seg1 ) );
			return;                                                        }         }

	if ( ( not tip_of_seg2_free ) and base_of_seg2_free and tip_of_seg1_free )
	{	seg2.remove_from_mesh ( interf );
		Cell seg = square.boundary().cell_in_front_of ( seg1.tip() );
		seg.add_to_mesh ( interf );
		if ( not seg1_on_bdry )
			set_of_squares.insert ( ambient.cell_in_front_of ( seg1 ) );
		return;                                                       }

	if ( seg2_on_bdry and base_of_seg2_free and tip_of_seg1_free )
	{	Cell seg = square.boundary().cell_in_front_of ( seg1.tip() );
		Cell A = seg.base().reverse(), B = seg.tip();
		double dx = x(B) - x(A);
		double dy = y(B) - y(A);
		double gx = dpsi_dx(A);
		double gy = dpsi_dy(A);
		if ( gx*dy > gy*dx )
		{	seg2.remove_from_mesh ( interf );
			seg.add_to_mesh ( interf );
			if ( not seg1_on_bdry )
				set_of_squares.insert ( ambient.cell_in_front_of ( seg1 ) );
			return;                                                       }  }

}  // end of  join_two_parallel_segs

//-----------------------------------------------------------------------------------//


void try_to_solve_square ( Cell sq, Mesh & ambient, Mesh & interf,
													 const Function psi, const Function psi_x, const Function psi_y,
													 std::set<Cell> & set_of_squares                                 )
			
{ CellIterator it_seg = sq.boundary().iterator ( tag::over_segments );
	int found = 0;
	Cell seg1 ( tag::non_existent );
	for ( it_seg.reset(); it_seg.in_range(); it_seg++ )
	{	seg1 = *it_seg;
		if ( seg1.belongs_to ( interf, tag::oriented ) )
		{	found = 1;  break;  }
		if ( seg1.reverse().belongs_to ( interf, tag::oriented ) )
		{	found = -1;  break;  }                                    }
	assert ( found != 0 );
	Cell seg2 = sq.boundary().cell_in_front_of ( seg1.tip() );
	if ( seg2.belongs_to ( interf, tag::oriented ) ) return;
	if ( seg2.reverse().belongs_to ( interf, tag::oriented ) ) return;
	seg2 = sq.boundary().cell_in_front_of ( seg2.tip() );
	assert ( ( seg2.belongs_to ( interf, tag::oriented ) ) or
	         ( seg2.reverse().belongs_to ( interf, tag::oriented ) ) );
	Cell seg3 = sq.boundary().cell_in_front_of ( seg2.tip() );
	if ( seg3.belongs_to ( interf, tag::oriented ) ) return;
	if ( seg3.reverse().belongs_to ( interf, tag::oriented ) ) return;

	if ( seg1.reverse().belongs_to ( interf, tag::oriented ) and
			 seg2.reverse().belongs_to ( interf, tag::oriented )     )
		if ( found == 1 ) join_two_opposite_segs_pos
			( sq, seg1.reverse(), seg2.reverse(), ambient, interf, psi, psi_x, psi_y, set_of_squares );
		else join_two_opposite_segs_neg
			( sq, seg1.reverse(), seg2.reverse(), ambient, interf, psi, psi_x, psi_y, set_of_squares );
	else if ( seg1.belongs_to ( interf, tag::oriented ) and
	          seg2.belongs_to ( interf, tag::oriented )     )
		if ( found == 1 ) join_two_opposite_segs_pos
			( sq, seg1, seg2, ambient, interf, psi, psi_x, psi_y, set_of_squares );
		else join_two_opposite_segs_neg
			( sq, seg1, seg2, ambient, interf, psi, psi_x, psi_y, set_of_squares );
	else if ( seg1.belongs_to ( interf, tag::oriented ) and
            seg2.reverse().belongs_to ( interf, tag::oriented ) )
		join_two_parallel_segs
			( sq, seg1, seg2.reverse(), ambient, interf, psi, psi_x, psi_y, set_of_squares );
	else if ( seg1.reverse().belongs_to ( interf, tag::oriented ) and
            seg2.belongs_to ( interf, tag::oriented )               )
		join_two_parallel_segs
			( sq, seg2, seg1.reverse(), ambient, interf, psi, psi_x, psi_y, set_of_squares );

}  // end of try_to_solve_square

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

	while ( set_of_squares.size() )  // while there are squares to "solve"

	{	std::forward_list<Cell> list_of_squares;
		std::set<Cell>::iterator it_set;
		// we keep in list_of_squares only those whose four vertices belong to 'interf'
		for ( it_set = set_of_squares.begin(); it_set != set_of_squares.end(); it_set++ )
		{	Cell sq = *it_set;
			bool keep = true;
			CellIterator it_ver = sq.boundary().iterator ( tag::over_vertices );
			for ( it_ver.reset(); it_ver.in_range(); it_ver++ )
			{	Cell P = *it_ver;
				keep = keep and P.belongs_to ( interf ); }
			if ( keep )
				list_of_squares.push_front ( sq );                                 }

		set_of_squares.clear();
		std::forward_list<Cell>::iterator it_list;
		for ( it_list = list_of_squares.begin(); it_list != list_of_squares.end(); it_list++ )
			try_to_solve_square ( *it_list, ambient, interf, psi, psi_x, psi_y, set_of_squares );  }
		// if not solved, *it will somehow appear later
		// if solved, one or two new squares will be added to set_of_squares
		
	// eliminate segments of 'interf' which are on the boundary of 'ambient' :
	{ // just a block of code for hiding 'it'
	CellIterator it = interf.iterator ( tag::over_segments );
	std::forward_list<Cell> list_of_segs;
	for ( it.reset(); it.in_range(); it++ )
	{	Cell seg = *it;
		if ( ( not ambient.cell_in_front_of ( seg, tag::may_not_exist ) .exists() )
		     or ( not ambient.cell_behind ( seg, tag::may_not_exist ) .exists()     ) )
			list_of_segs.push_front ( seg );                                                }
	std::forward_list<Cell>::iterator it_list;
	for ( it_list = list_of_segs.begin(); it_list != list_of_segs.end(); it_list++ )
		(*it_list).remove_from_mesh ( interf );
	} // just a block of code for hiding 'it'
		
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



