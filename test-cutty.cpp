
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

	double radius = 0.4;
	// std::cout << "radius = ";  std::cin >> radius;
	Function psi = x*x + (y-0.2)*(y-0.2) - 0.3 + 0.2*x*y - 1.5*x*x*y*y;

	// 0.5 * ( ( x*x + (y-0.2)*(y-0.2) ) / radius - radius )  circulo

	// x*x + (y-0.2)*(y-0.2) - 0.3 + 0.2*x*y - 1.35*x*x*y*y
	// x*x + (y-0.2)*(y-0.2) - 0.3 + 0.2*x*y - 1.5*x*x*y*y

	Mesh cut = build_interface ( rect_mesh, psi );
	// builds a set of vertices next to the level line  psi == 0
	// then builds a one-dimensional mesh

	// improve_interf_90 ( rect_mesh, cut, psi );
	// reverts some of the angles

	// improve_interf_135 ( rect_mesh, cut, psi );
	// splits some rectangles in two triangles
	
	special_draw ( rect_mesh, cut, "square-cut.eps" );

	std::cout << "produced file square-cut.eps" << std::endl;
	
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
	
	CellIterator it_seg = ambient.iterator ( tag::over_segments );
	for ( it_seg.reset(); it_seg.in_range(); it_seg++ )
	{	Cell seg = *it_seg;
		Cell A = seg.base().reverse();
		Cell B = seg.tip();
		if ( not opposite_signs ( psi(A), psi(B) ) ) continue;
		ms.insert ( A );
		ms.insert ( B );                                        }

	// at this moment we have a cloud of points, possibly repeated

	Mesh interf ( tag::fuzzy, tag::of_dimension, 1 );
	// empty mesh, it will be returned after adding segments to it

	std::set < Cell > set_useful, set_needed;
	std::forward_list < Cell > list_useful;
	std::multiset<Cell,compare_values_of>::iterator it_ms;

	for ( it_ms = ms.begin(); it_ms != ms.end(); it_ms++ )
	{	Cell A = * it_ms;
		if ( set_useful.find ( A ) == set_useful.end() )
		{	set_useful.insert ( A );
			list_useful.push_front ( A );                   }  }

	std::cout << "at the beginning, we have " << set_useful.size() << " vertices" << std::endl;
	int counter = 0;
	
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
		
	std::cout << "in the end, we have " << set_useful.size() << " vertices" << std::endl << std::flush;
	
	for ( it_list = list_useful.begin(); it_list != list_useful.end(); it_list++ )
	{	Cell A = *it_list;
		if ( set_useful.find ( A ) == set_useful.end() ) continue;
		std::cout << x(A) << " " << y(A) << " moveto ";
		std::cout << x(A) << " " << y(A) << " 0.016 0 360 arc fill" << std::endl;  }

	exit ( 0 );
	return interf;
	
	for ( it_list = list_useful.begin(); it_list != list_useful.end(); it_list++ )

	{	Cell A = *it_list;
		
		double xA = x(A), yA = y(A);
		double gx = psi_x(A), gy = psi_y(A);

		CellIterator it_around = ambient.iterator
			( tag::over_cells_of_dim, 2, tag::around, A );
		int found_pos = 0, found_neg = 0;
		std::set < Cell > set_of_neigh;
		for ( it_around.reset(); it_around.in_range(); it_around++ )
		{	Cell cll = *it_around;
			CellIterator it_ver = cll.boundary().iterator ( tag::over_vertices );
			for ( it_ver.reset(); it_ver.in_range(); it_ver++ )
			{	Cell B = *it_ver;
				if ( set_useful.find(B) == set_useful.end() ) continue;
				set_of_neigh.insert ( B );                       }                   }
		std::set<Cell>::iterator it_neigh;
		for ( it_neigh = set_of_neigh.begin(); it_neigh != set_of_neigh.end(); it_neigh++ )
		{	Cell B = *it_neigh;
			double dx = x(B) - xA, dy = y(B) - yA;
			if ( dx*gy > dy*gx ) found_pos++;
			if ( dx*gy < dy*gx ) found_neg++;      }
		std::cout << found_pos << " " << found_neg << std::endl << std::flush;

	}  // end of for it_useful

	return interf;

}  // end of  build_interface

//-----------------------------------------------------------------------------------//



void improve_interf_90 ( Mesh ambient, Mesh interf, Function psi )

{	std::cout << "improve_interf_90" << std::endl;
	
	// we build a multi-set of vertices, ordered according to the values of psi
	// more precisely, the first ones will be those where the absolute value of psi is large

	compare_values_of comp_psi ( abs(psi) );
	std::multiset < Cell, compare_values_of > ms ( comp_psi );

	CellIterator it_ver = interf.iterator ( tag::over_vertices );
	for ( it_ver.reset(); it_ver.in_range(); it_ver++ ) ms.insert ( *it_ver );

	std::multiset<Cell,compare_values_of>::iterator it_ms;
	for ( it_ms = ms.begin(); it_ms != ms.end(); it_ms++ )

	{	Cell A = *it_ms;

		// suppose the interface has moved away from A in the meanwhile :
		if ( not A.belongs_to ( interf ) ) continue;

		Cell prev_seg = interf.cell_behind ( A, tag::may_not_exist );
		if ( not prev_seg.exists() )  // we are at the boundary
		{	std::cout << "we are at the boundary, not prev_seg.exists() " << std::endl;
			Cell AB = interf.cell_in_front_of ( A, tag::surely_exists );
			Cell ABDC = ambient.cell_behind ( AB );
			Cell CA = ABDC.boundary().cell_behind ( A );
			assert ( not ambient.cell_in_front_of ( CA, tag::may_not_exist ) .exists() );
			Cell C = CA.base().reverse();
			Cell B = AB.tip();
 			Cell BD = ABDC.boundary().cell_in_front_of ( B );
			if ( BD.belongs_to ( interf, tag::oriented ) )
			if ( ( abs ( psi ( C ) ) < abs ( psi ( A ) ) ) and
			     ( abs ( psi ( C ) ) < abs ( psi ( B ) ) )     )
			// we give up A and B and insert C
			{	AB.remove_from_mesh ( interf );
				BD.remove_from_mesh ( interf );
				Cell D = BD.tip();
				Cell DC = ABDC.boundary().cell_in_front_of ( D );
				assert ( DC.tip() == C );
				DC.reverse().add_to_mesh ( interf );
				ms.insert ( C );
				continue;                                           }
			Cell AEFB = ambient.cell_in_front_of ( AB );
			Cell FB = AEFB.boundary().cell_behind ( B );
			Cell F = FB.base().reverse();
			Cell EF = AEFB.boundary().cell_behind ( F );
			assert ( EF.tip() == F );
			Cell E = EF.base().reverse();
			if ( FB.reverse().belongs_to ( interf, tag::oriented ) )
			if ( ( abs ( psi ( E ) ) < abs ( psi ( A ) ) ) and
			     ( abs ( psi ( E ) ) < abs ( psi ( B ) ) )     )
			// we give up A and B and insert E
			{	AB.remove_from_mesh ( interf );
				FB.reverse().remove_from_mesh ( interf );
				EF.add_to_mesh ( interf );
				ms.insert ( E );
				continue;                                                }
			continue;                                                    }
		// end of if

		Cell next_seg = interf.cell_in_front_of ( A, tag::may_not_exist );
		if ( not next_seg.exists() )  // we are at the boundary
		{	std::cout << "we are at the boundary, not next_seg.exists()" << std::endl;
			Cell BA = interf.cell_behind ( A, tag::surely_exists );
			Cell B = BA.base().reverse();
			Cell BACD = ambient.cell_behind ( BA );
			Cell DB = BACD.boundary().cell_behind ( B );
			Cell D = DB.base().reverse();
			Cell CD = BACD.boundary().cell_behind ( D );
			assert ( CD.tip() == D );
			Cell C = CD.base().reverse();
			if ( DB.belongs_to ( interf, tag::oriented ) )
			if ( ( abs ( psi ( C ) ) < abs ( psi ( A ) ) ) and
			     ( abs ( psi ( C ) ) < abs ( psi ( B ) ) )     )
			// we give up A and B and insert C
			{	BA.remove_from_mesh ( interf );
				DB.remove_from_mesh ( interf );
				CD.reverse().add_to_mesh ( interf );
				ms.insert ( C );
				continue;                                           }
			Cell ABFE = ambient.cell_in_front_of ( BA );
			Cell BF = ABFE.boundary().cell_in_front_of ( B );
			Cell F = BF.tip();
			Cell FE = ABFE.boundary().cell_in_front_of ( F );
			assert ( FE.base() == F.reverse() );
			Cell E = FE.tip();
			if ( BF.reverse().belongs_to ( interf, tag::oriented ) )
			if ( ( abs ( psi ( E ) ) < abs ( psi ( A ) ) ) and
			     ( abs ( psi ( E ) ) < abs ( psi ( B ) ) )     )
			// we give up A and B and insert E
			{	BA.remove_from_mesh ( interf );
				BF.reverse().remove_from_mesh ( interf );
				FE.add_to_mesh ( interf );
				ms.insert ( E );
				continue;                                                }
			continue;                                                    }
		// end of if

		Cell sq1 = ambient.cell_in_front_of ( next_seg ),
		     sq2 = ambient.cell_behind ( next_seg ),
		     sq3 = ambient.cell_in_front_of ( prev_seg ),
		     sq4 = ambient.cell_behind ( prev_seg );
	  assert ( ( sq1 != sq2 ) and ( sq1 != sq4 ) and ( sq2 != sq3 ) and ( sq3 != sq4 ) );
		if ( sq1 == sq3 )  // 'interf' turns right at A
		{	std::cout << "interf turns right" << std::endl;
			Cell B = next_seg.tip();
			Cell CB = sq1.boundary().cell_behind ( B );
			Cell C = CB.base().reverse();
			Cell DC = sq1.boundary().cell_behind ( C );
			Cell D = DC.base().reverse();
			assert ( D == prev_seg.base().reverse() );
			if ( abs ( psi ( C ) ) < abs ( psi ( A ) ) )
			{	prev_seg.remove_from_mesh ( interf );
				next_seg.remove_from_mesh ( interf );
				bool bc_ok = true, cd_ok = true;
				if ( CB.reverse().belongs_to ( interf, tag::oriented ) )
				{	CB.reverse().remove_from_mesh ( interf );  bc_ok = false;  }
				else if ( not interf.cell_in_front_of ( B, tag::may_not_exist ) .exists() )
				{	assert ( not ambient.cell_in_front_of ( CB, tag::may_not_exist ) .exists() );
					bc_ok = false;  }
				if ( DC.reverse().belongs_to ( interf, tag::oriented ) )
				{	DC.reverse().remove_from_mesh ( interf );  cd_ok = false;  }
				else if ( not interf.cell_behind ( D, tag::may_not_exist ) .exists() )
				{	assert ( not ambient.cell_in_front_of ( DC, tag::may_not_exist ) .exists() );
					cd_ok = false;  }
				if ( bc_ok ) CB.add_to_mesh ( interf );
				if ( cd_ok ) DC.add_to_mesh ( interf );
				ms.insert ( C );                                         }
			continue;                                                        }
		// end of if  sq1 == sq3

		if ( sq2 == sq4 )  // 'interf' turns left at A
		{	std::cout << "interf turns left" << std::endl;
			Cell B = next_seg.tip();
			Cell BC = sq2.boundary().cell_in_front_of ( B );
			Cell C = BC.tip();
			Cell CD = sq2.boundary().cell_in_front_of ( C );
			Cell D = CD.tip();
			assert ( D == prev_seg.base().reverse() );
			if ( abs ( psi ( C ) ) < abs ( psi ( A ) ) )
			{	prev_seg.remove_from_mesh ( interf );
				next_seg.remove_from_mesh ( interf );
				bool bc_ok = true, cd_ok = true;
				if ( BC.belongs_to ( interf, tag::oriented ) )
				{	BC.remove_from_mesh ( interf );  bc_ok = false;  }
				else if ( not interf.cell_in_front_of ( B, tag::may_not_exist ) .exists() )
				{	assert ( not ambient.cell_in_front_of ( BC, tag::may_not_exist ) .exists() );
					bc_ok = false;  }
				if ( CD.belongs_to ( interf, tag::oriented ) )
				{	CD.remove_from_mesh ( interf );  cd_ok = false;  }
				else if ( not interf.cell_behind ( D, tag::may_not_exist ) .exists() )
				{	assert ( not ambient.cell_in_front_of ( CD, tag::may_not_exist ) .exists() );
					cd_ok = false;  }
				if ( bc_ok ) BC.reverse().add_to_mesh ( interf );
				if ( cd_ok ) CD.reverse().add_to_mesh ( interf );
				ms.insert ( C );
				continue;                                                }
			continue;                                                        }
		// end of if  sq2 == sq4
	} // end of for it_ms

}  // end of  improve_interf_90

//-----------------------------------------------------------------------------------//


void improve_interf_135 ( Mesh ambient, Mesh interf, Function psi )

{	std::cout << "improve_interf_135" << std::endl;
	
	// we build a multi-set of vertices, ordered according to the values of psi
	// more precisely, the first ones will be those where the absolute value of psi is large

	compare_values_of comp_psi ( abs(psi) );
	std::multiset < Cell, compare_values_of > ms ( comp_psi );

	CellIterator it_ver = interf.iterator ( tag::over_vertices );
	for ( it_ver.reset(); it_ver.in_range(); it_ver++ ) ms.insert ( *it_ver );

	std::multiset<Cell,compare_values_of>::iterator it_ms;
	for ( it_ms = ms.begin(); it_ms != ms.end(); it_ms++ )

	{	Cell A = *it_ms;
		
		Cell prev_seg = interf.cell_behind ( A, tag::may_not_exist );
		if ( not prev_seg.exists() )  // we are at the boundary
		{	std::cout << "we are at the boundary, not prev_seg.exists(), ";
			Cell AB = interf.cell_in_front_of ( A, tag::surely_exists );
			if ( psi ( A ) > 0. )
			{	Cell ABDC = ambient.cell_behind ( AB );
				std::cout << "psi ( A ) > 0." << std::endl;
				Cell CA = ABDC.boundary().cell_behind ( A );
				assert ( not ambient.cell_in_front_of ( CA, tag::may_not_exist ) .exists() );
				Cell C = CA.base().reverse();
				if ( - psi ( C ) > psi ( A ) ) continue;
				Cell B = AB.tip();
				Cell BD = ABDC.boundary().cell_in_front_of ( B );
				Cell D = BD.tip();
				Cell DC = ABDC.boundary().cell_in_front_of ( D );
				assert ( DC.tip() == C );
				ABDC.remove_from_mesh ( ambient );
				AB.remove_from_mesh ( interf );
				Cell BC ( tag::segment, B.reverse(), C );
				Cell ABC ( tag::triangle, AB, BC, CA );
				Cell BDC ( tag::triangle, BD, DC, BC.reverse() );
				ABC.add_to_mesh ( ambient );
				BDC.add_to_mesh ( ambient );
				BC.reverse().add_to_mesh ( interf );          }
			else  // psi ( A ) <= 0
			{ Cell ACDB = ambient.cell_in_front_of ( AB );
				std::cout << "psi ( A ) <= 0." << std::endl;
				Cell AC = ACDB.boundary().cell_in_front_of ( A );
				assert ( not ambient.cell_in_front_of ( AC, tag::may_not_exist ) .exists() );
				Cell C = AC.tip();
				if ( psi ( C ) > - psi ( A ) ) continue;
				Cell B = AB.tip();
				Cell DB = ACDB.boundary().cell_behind ( B );
				Cell D = DB.base().reverse();
				Cell CD = ACDB.boundary().cell_in_front_of ( C );
				assert ( CD.tip() == D );
				ACDB.remove_from_mesh ( ambient );
				AB.remove_from_mesh ( interf );
				Cell BC ( tag::segment, B.reverse(), C );
				Cell ACB ( tag::triangle, AC, BC.reverse(), AB.reverse() );
				Cell BCD ( tag::triangle, BC, CD, DB );
				ACB.add_to_mesh ( ambient );
				BCD.add_to_mesh ( ambient );
				BC.reverse().add_to_mesh ( interf );                           }
		}

	}  // end of for it_ms
}  // end of  improve_interf_135
	
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



