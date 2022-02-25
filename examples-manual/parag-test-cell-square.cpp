
// square periodicity, circular hole, triangular elements, formulation in strain

#include "maniFEM.h"
#include <fstream>
#include <set>
#include <Eigen/Sparse>
#include <Eigen/OrderingMethods>
using namespace maniFEM;
using namespace std;

template <typename T>
class myTensor
// constructors initialize to zero
{ // data:
	public:
	vector<T> elements;
	list<size_t> dimensions, cumulative_dims;
	size_t total_dim;
	// constructors:
	myTensor () {};
	myTensor (const list<size_t> dims)
	{	dimensions = dims;
		allocate_space (); };
	myTensor (const list<char> dims)
	{	dimensions = str2list(dims);
		allocate_space ();           };
	myTensor (const char dims[])
	{	dimensions = str2list(dims);
		allocate_space ();           };
	myTensor (size_t i, size_t j,
	          size_t k, size_t l)
	{	dimensions.push_back(i);
		dimensions.push_back(j);
		dimensions.push_back(k);
		dimensions.push_back(l);
		allocate_space ();       };
	myTensor ( size_t i, size_t j, size_t k )
	{	dimensions.push_back(i);
		dimensions.push_back(j);
		dimensions.push_back(k);
		allocate_space ();       };
  myTensor (size_t i, size_t j)
	{	dimensions.push_back(i);
		dimensions.push_back(j);
		allocate_space ();       };
  myTensor (int i)
	{	dimensions.push_back(i);
		allocate_space ();       };
	~myTensor () { };
	// methods:
	list<size_t> str2list (const list<char> lc)
	{	const size_t izero = int('0');
		list<char>::iterator i;
		for ( i=lc.begin(); i!=lc.end(); i++ )
		{	assert (*i >= '0');
			assert (*i <= '9'); }
		list<size_t> li;
		for ( i=lc.begin(); i!=lc.end(); i++ )
			li.push_back(int(*i)-izero);
		return li;                             }
	list<size_t> str2list (const char lc[])
	{	const size_t izero = int('0');
		for (int i=0;i<strlen(lc);i++)
		{	assert (lc[i] >= '0');
			assert (lc[i] <= '9'); }
		list<size_t> li;
		for (int i=0;i<strlen(lc);i++)
			li.push_back(int(lc[i])-izero);
		return li;                        }
	void allocate_space ()
	{	total_dim = 1;
		list<size_t>::iterator k;
		for ( k=dimensions.begin(); k!=dimensions.end(); k++ )
		{	cumulative_dims.push_back(total_dim);
			total_dim *= *k;                       }
		// elements.reserve (total_dim);
		elements.insert (elements.end(),total_dim,0.0);
		assert (elements.size() == total_dim);                 };
	T& operator()(list<size_t> index)
	{	assert (index.size() == dimensions.size());
		size_t pointer = 0;
		list<size_t>::iterator i,d,cd;
		for ( i=index.begin(), d=dimensions.begin(),
	        cd=cumulative_dims.begin();
	        i!=index.end(); i++,d++,cd++          )
		{	assert (*i >= 0);
			assert (*i < *d);
			pointer += (*i)*(*cd); }
		return (elements[pointer]);                    }
	T& operator()(const char index[])
	{	return operator()(str2list(index)); }
	T& operator() ( size_t i, size_t j,
	                size_t k, size_t l )
	{	assert ( dimensions.size() == 4 );
		list<size_t> index;
		index.push_back(i);
		index.push_back(j);
		index.push_back(k);
		index.push_back(l);
		return operator()(index);         }
	T& operator() ( size_t i, size_t j, size_t k )
	{	assert ( dimensions.size() == 3 );
		list<size_t> index;
		index.push_back(i);
		index.push_back(j);
		index.push_back(k);
		return operator()(index);         }
	T& operator() ( size_t i, size_t j )
	{	assert ( dimensions.size() == 2 );
		list<size_t> index;
		index.push_back(i);
		index.push_back(j);
		return operator()(index);          }
  T& operator()(int i)
	{	assert ( dimensions.size() == 1 );
		list<size_t> index(1,i);
		return operator()(index);         }
};  // end of  class myTensor


void limit_number_of_neighbours ( Mesh msh );
	
void remove_short_segments ( Mesh & msh, double threshold );

void flip_split_long_segments ( Mesh & msh, double threshold );

void baricenters ( Mesh & msh );



int main ( )

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0], y = xy[1];

	size_t n = 10;
	double l = 2.6;
	double d = l / double(n);
	double areaY = l*l ;

	Cell A ( tag::vertex );  x(A) = -l/2.;  y(A) = -l/2.;
	Cell B ( tag::vertex );  x(B) =  l/2.;  y(B) = -l/2.;
	Cell C ( tag::vertex );  x(C) =  l/2.;  y(C) =  l/2.;
	Cell D ( tag::vertex );  x(D) = -l/2.;  y(D) =  l/2.;

	Mesh AB ( tag::segment, A.reverse(), B, tag::divided_in, n );
	Mesh BC ( tag::segment, B.reverse(), C, tag::divided_in, n );
	Mesh CD ( tag::segment, C.reverse(), D, tag::divided_in, n );
	Mesh DA ( tag::segment, D.reverse(), A, tag::divided_in, n );

	Manifold circle = RR2.implicit ( x*x + y*y == 0.7 );
	Mesh inner ( tag::frontal, tag::entire_manifold, circle, tag::desired_length, d );

	Mesh bdry ( tag::join, AB, BC, CD, DA, inner.reverse() );

	RR2.set_as_working_manifold();
	Mesh square ( tag::frontal, tag::boundary, bdry, tag::desired_length, d );

	Mesh torus = square.fold ( tag::identify, AB, tag::with, CD.reverse(),
	                           tag::identify, BC, tag::with, DA.reverse(),
	                           tag::use_existing_vertices                 );

	// std::cout << "produced folded mesh, now drawing, please wait" << std::endl << std::flush;
	
	// torus.draw_ps ( "torus.eps", tag::unfold,
  //                 tag::over_region, -2.1 < x < 4.3, -3.6 < y < 2.1 );

	// Hooke's Law,  E = 1.,  nu = 0.3
	double lambda = 0.576923, mu = 0.38461538;

	myTensor <double> Hooke (2,2,2,2);
	Hooke (0,0,0,0) = 2*mu + lambda;
	Hooke (0,0,0,1) = 0.;
	Hooke (0,0,1,0) = 0.;
	Hooke (0,0,1,1) = lambda;
	Hooke (0,1,0,0) = 0.;
	Hooke (0,1,0,1) = mu;
	Hooke (0,1,1,0) = mu;
	Hooke (0,1,1,1) = 0.;
	Hooke (1,0,0,0) = 0.;
	Hooke (1,0,0,1) = mu;
	Hooke (1,0,1,0) = mu;
	Hooke (1,0,1,1) = 0.;
	Hooke (1,1,0,0) = lambda;
	Hooke (1,1,0,1) = 0.;
	Hooke (1,1,1,0) = 0.;
	Hooke (1,1,1,1) = 2*mu + lambda;

	// declare the type of finite element
	FiniteElement fe ( tag::with_master, tag::triangle, tag::Lagrange, tag::of_degree, 1 );
	fe.set_integrator ( tag::Gauss, tag::tri_4 );

	std::map < Cell, size_t > numbering;
	{ // just a block of code for hiding 'it' and 'counter'
	Mesh::Iterator it = square .iterator ( tag::over_vertices );
	size_t counter = 0;
	for ( it .reset() ; it .in_range(); it ++ )
	{	Cell V = *it;  numbering [V] = counter;  ++ counter;  }
	assert ( counter == numbering .size() );
	} // just a block of code

	size_t number_dofs = numbering.size();
	std::cout << "global matrix " << 2*number_dofs << "x" << 2*number_dofs << std::endl;
	Eigen::SparseMatrix < double > matrix_A ( 2*number_dofs, 2*number_dofs );
	
	matrix_A.reserve ( Eigen::VectorXi::Constant ( 2*number_dofs, 14 ) );
	// since we will be working with a mesh of triangles,
	// there will be, in average, 14 = 2*(6+1) non-zero elements per column
	// the diagonal entry plus six neighbour vertices

	// we fill the main diagonal with ones
	// then we put zero for vertices belonging to 'torus'
	{ // just a block of code for hiding 'it'
	for ( size_t i = 0; i < 2*number_dofs; i++ ) matrix_A.coeffRef ( i, i ) = 1.;
	Mesh::Iterator it = torus.iterator ( tag::over_vertices );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell V = *it;
		matrix_A.coeffRef ( 2*numbering[V], 2*numbering[V] ) = 0.;
		matrix_A.coeffRef ( 2*numbering[V]+1, 2*numbering[V]+1 ) = 0.;		}
	} // just a block of code for hiding 'it'
	
	Eigen::VectorXd vector_b ( 2*number_dofs ), vector_sol ( 2*number_dofs );
	vector_b.setZero();

	xy = Manifold::working.coordinates();
	x = xy[0];  y = xy[1];

	// macroscopic temperature gradient
	myTensor < double > macro_strain (2,2);
	macro_strain (0,0) = 1.;
	macro_strain (0,1) = 0.;
	macro_strain (1,0) = 0.;
	macro_strain (1,1) = 1.;

	Function::Jump jump_of_u_x = macro_strain(0,0) * x.jump() + macro_strain(0,1) * y.jump();
	Function::Jump jump_of_u_y = macro_strain(1,0) * x.jump() + macro_strain(1,1) * y.jump();

	// run over all triangular cells composing 'torus'
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = torus .iterator ( tag::over_cells_of_max_dim );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell small_tri = *it;
		fe .dock_on ( small_tri, tag::winding );
		// run twice over the three vertices of 'small_tri'
		Mesh::Iterator it_V = small_tri .boundary() .iterator ( tag::over_vertices );
		for ( it_V .reset(); it_V .in_range(); it_V ++ )
		{	Cell V = * it_V;
			// perhaps implement an iterator returning a vertex and a segment
			Cell seg = small_tri .boundary() .cell_in_front_of ( V );
			Cell W = V;
			Function psi_V = fe .basis_function ( V ),
			         d_psiV_dx = psi_V .deriv ( x ),
			         d_psiV_dy = psi_V .deriv ( y );
			double jump_V_W_u_x = 0.;
			double jump_V_W_u_y = 0.;
			while ( true )
			{	assert ( W == seg .base() .reverse() );
				// V may be the same as W, no problem about that
				Function psi_W = fe .basis_function ( W ),
				         d_psiW_dx = psi_W .deriv ( x ),
				         d_psiW_dy = psi_W .deriv ( y );
				// 'fe' is already docked on 'small_tri' so this will be the domain of integration
				double int_d_psiW_dx_d_psiV_dx = fe .integrate ( d_psiW_dx * d_psiV_dx );
				double int_d_psiW_dx_d_psiV_dy = fe .integrate ( d_psiW_dx * d_psiV_dy );
				double int_d_psiW_dy_d_psiV_dx = fe .integrate ( d_psiW_dy * d_psiV_dx );
				double int_d_psiW_dy_d_psiV_dy = fe .integrate ( d_psiW_dy * d_psiV_dy );
				myTensor < double > energy ( 2, 2 );
				for ( size_t i = 0; i < 2; i ++ )
				for ( size_t k = 0; k < 2; k ++ )
					energy (i,k) = Hooke (i,0,k,0) * int_d_psiW_dx_d_psiV_dx +
					               Hooke (i,0,k,1) * int_d_psiW_dy_d_psiV_dx +
					               Hooke (i,1,k,0) * int_d_psiW_dx_d_psiV_dy +
				                 Hooke (i,1,k,1) * int_d_psiW_dy_d_psiV_dy ;
				matrix_A .coeffRef ( 2*numbering[V], 2*numbering[W] ) += energy(0,0);
				matrix_A .coeffRef ( 2*numbering[V], 2*numbering[W]+1 ) += energy(0,1);
				matrix_A .coeffRef ( 2*numbering[V]+1, 2*numbering[W] ) += energy(1,0);
				matrix_A .coeffRef ( 2*numbering[V]+1, 2*numbering[W]+1 ) += energy(1,1);
				vector_b ( 2*numbering[V] ) -=
					jump_V_W_u_x * energy(0,0) + jump_V_W_u_y * energy(0,1);
				vector_b ( 2*numbering[V]+1 ) -=
					jump_V_W_u_x * energy(1,0) + jump_V_W_u_y * energy(1,1);
				jump_V_W_u_x += jump_of_u_x ( seg .winding() );
				jump_V_W_u_y += jump_of_u_y ( seg .winding() );
				W = seg.tip();
				if ( V == W ) break;
				seg = small_tri .boundary() .cell_in_front_of ( W );                        }
			// end of loop in W
			// here  jump_V_W_u_x and jump_V_W_u_y  should be zero again
			// but we do not assert that, rounding errors may mess up things
		}  }  // end of loop in V, end of loop in small_tri
	} // just a block of code for hiding 'it'

	std::cout << "now solving the system of linear equations" << std::endl;
	
	matrix_A .makeCompressed();

	Eigen::SparseQR < Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;

	solver.compute ( matrix_A );
	if ( solver.info() != Eigen::Success )
	{	std::cout << "Eigen solver.compute failed" << std::endl;
		exit ( 0 );                                              }

	vector_sol = solver.solve ( vector_b );
	if ( solver.info() != Eigen::Success )
	{	std::cout << "Eigen solver.solve failed" << std::endl;
		exit ( 0 );                                            }
		
	Eigen::VectorXd new_b = matrix_A * vector_sol;

	myTensor<double> macro_stress(2,2);
	for ( size_t i=0; i<2; i++ )
	for ( size_t j=0; j<2; j++ )  macro_stress (i,j) = 0.;

	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = torus .iterator ( tag::over_cells_of_max_dim );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell small_tri = *it;
		fe .dock_on ( small_tri, tag::winding );
		// run twice over the four vertices of 'small_tri'
		Mesh::Iterator it_V = small_tri .boundary() .iterator ( tag::over_vertices );
		double jump_V_u_x = 0., jump_V_u_y = 0.; 
		for ( it_V .reset(); it_V .in_range(); it_V++ )
		{	Cell V = *it_V;
			Function psi_V = fe .basis_function ( V ),
			         d_psi_V_dx = psi_V .deriv ( x ),
			         d_psi_V_dy = psi_V .deriv ( y );
			for ( size_t i=0; i<2; i++ )
			for ( size_t j=0; j<2; j++ )
				macro_stress (i,j) +=
					( vector_sol ( 2*numbering[V] ) + jump_V_u_x ) *
					  ( Hooke (i,j,0,0) * fe .integrate ( d_psi_V_dx ) +
					    Hooke (i,j,0,1) * fe .integrate ( d_psi_V_dy )  ) +
					( vector_sol ( 2*numbering[V]+1 ) + jump_V_u_y ) *
					  ( Hooke (i,j,1,0) * fe .integrate ( d_psi_V_dx ) +
					    Hooke (i,j,1,1) * fe .integrate ( d_psi_V_dy )  )  ;
			Cell seg = small_tri .boundary() .cell_in_front_of (V);
			jump_V_u_x += jump_of_u_x ( seg .winding() );
			jump_V_u_y += jump_of_u_y ( seg .winding() );                      }          }
			// end of loop in V, end of loop in small_tri
			// at the end of each loop in V,  jump_V_u_x and jump_V_u_y  should be zero again
			// but we do not assert that, rounding errors may mess up things
	} // just a block of code for hiding 'it'
	
	for ( size_t i=0; i<2; i++ )
	for ( size_t j=0; j<2; j++ )  macro_stress (i,j) /= areaY;
	
	cout << "macro strain " << macro_strain(0,0) << " " << macro_strain(0,1) << " "
	                        << macro_strain(1,0) << " " << macro_strain(1,1) << " " << endl;
	cout << "macro stress " << macro_stress(0,0) << " " << macro_stress(0,1) << " "
	                        << macro_stress(1,0) << " " << macro_stress(1,1) << " " << endl;

	RR2 .set_as_working_manifold();
	xy = Manifold::working .coordinates();
	x = xy[0];  y = xy[1];
	
	// we define the solution in all vertices of 'square'	(those not belonging to 'torus')
	// it is easier to impose the zero average condition on 'square'
	// it is also easier to export_to_file
	
	{ // just a block of code for hiding variables
	size_t j = numbering [ B ];
	assert ( not A .belongs_to ( torus ) );
	vector_sol [ 2 * numbering[A] ] = vector_sol [ 2*j ]
		+ macro_strain(0,0) * ( x ( A ) - x ( B ) )
		+ macro_strain(0,1) * ( y ( A ) - y ( B ) );
	vector_sol [ 2 * numbering[A] + 1 ] =  vector_sol [ 2*j+1 ]
		+ macro_strain(1,0) * ( x ( A ) - x ( B ) )
		+ macro_strain(1,1) * ( y ( A ) - y ( B ) );
	assert ( not C .belongs_to ( torus ) );
	vector_sol [ 2 * numbering[C] ] = vector_sol [ 2*j ]
		+ macro_strain(0,0) * ( x ( C ) - x ( B ) )
		+ macro_strain(0,1) * ( y ( C ) - y ( B ) );
	vector_sol [ 2 * numbering[C] + 1 ] =  vector_sol [ 2*j+1 ]
		+ macro_strain(1,0) * ( x ( C ) - x ( B ) )
		+ macro_strain(1,1) * ( y ( C ) - y ( B ) );
	assert ( not D .belongs_to ( torus ) );
	vector_sol [ 2 * numbering[D] ] = vector_sol [ 2*j ]
		+ macro_strain(0,0) * ( x ( D ) - x ( B ) )
		+ macro_strain(0,1) * ( y ( D ) - y ( B ) );
	vector_sol [ 2 * numbering[D] + 1 ] =  vector_sol [ 2*j+1 ]
		+ macro_strain(1,0) * ( x ( D ) - x ( B ) )
		+ macro_strain(1,1) * ( y ( D ) - y ( B ) );
	Mesh::Iterator it_AB = AB .iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it_CD = CD .iterator ( tag::over_vertices, tag::backwards );
	it_AB .reset();  assert ( it_AB .in_range() );
	assert ( *it_AB == A );  it_AB ++;
	it_CD .reset();  assert ( it_CD .in_range() );
	assert ( *it_CD == D );  it_CD ++;
	for ( ; ; it_AB ++, it_CD ++ )
	{	assert ( it_AB .in_range() );  assert ( it_CD .in_range() );
		Cell V = *it_AB, W = *it_CD;
		if ( V == B )  {  assert ( W == C );  break;  }
		assert ( not W .belongs_to ( torus ) );
		j = numbering [ V ];
		vector_sol [ 2 * numbering[W] ] = vector_sol [ 2*j ]
			+ macro_strain(0,0) * ( x ( W ) - x ( V ) )
			+ macro_strain(0,1) * ( y ( W ) - y ( V ) );
		vector_sol [ 2 * numbering[W] + 1 ] =  vector_sol [ 2*j+1 ]
			+ macro_strain(1,0) * ( x ( W ) - x ( V ) )
			+ macro_strain(1,1) * ( y ( W ) - y ( V ) );                 }
	Mesh::Iterator it_BC = BC .iterator ( tag::over_vertices, tag::require_order );
	Mesh::Iterator it_DA = DA .iterator ( tag::over_vertices, tag::backwards );
	it_BC .reset();  assert ( it_BC .in_range() );
	assert ( *it_BC == B );  it_BC ++;
	it_DA .reset();  assert ( it_DA .in_range() );
	assert ( *it_DA == A );  it_DA ++;
	for ( ; ; it_BC ++, it_DA ++ )
	{	assert ( it_BC .in_range() );  assert ( it_DA .in_range() );
		Cell V = *it_BC, W = *it_DA;
		if ( V == C )  {  assert ( W == D );  break;  }
		assert ( not W .belongs_to ( torus ) );
		j = numbering [ V ];
		vector_sol [ 2 * numbering[W] ] = vector_sol [ 2*j ]
			+ macro_strain(0,0) * ( x ( W ) - x ( V ) )
			+ macro_strain(0,1) * ( y ( W ) - y ( V ) );
		vector_sol [ 2 * numbering[W] + 1 ] =  vector_sol [ 2*j+1 ]
			+ macro_strain(1,0) * ( x ( W ) - x ( V ) )
			+ macro_strain(1,1) * ( y ( W ) - y ( V ) );                 }
	} // just a block of code

	// impose sum u_x = 0, sum u_y = 0
	{ // just a block of code for hiding variables
	double sum_x = 0., sum_y = 0.;
	for ( size_t i = 0; i < number_dofs; i ++ )
	{	sum_x += vector_sol [ 2*i ];
		sum_y += vector_sol [ 2*i+1 ];  }
	sum_x /= number_dofs;
	sum_y /= number_dofs;
	for ( size_t i = 0; i < number_dofs; i ++ )
	{	vector_sol [ 2*i ] -= sum_x;
		vector_sol [ 2*i+1 ] -= sum_y;  }	
	} // just a block of code
	
	square.export_to_file ( tag::msh, "cell-elast-strain-cell-sq.msh", numbering );

	{ // just a block of code for hiding variables
	std::ofstream solution_file ("cell-elast-strain-cell-sq.msh", std::fstream::app );
	solution_file << "$NodeData" << std::endl;
	solution_file << "1" << std::endl;   // one string follows
	solution_file << "\"elastic displacement\"" << std::endl;
	solution_file << "1" << std::endl;   //  one real follows
	solution_file << "0.0" << std::endl;  // time [??]
	solution_file << "3" << std::endl;   // three integers follow
	solution_file << "0" << std::endl;   // time step [??]
	solution_file << "3" << std::endl;  // scalar values of u
	solution_file << square .number_of ( tag::vertices ) << std::endl;
	// number of values listed below
	Mesh::Iterator it = square .iterator ( tag::over_vertices );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell P = *it;
		size_t j = numbering [ P ];
		solution_file << j + 1 << " " << vector_sol [ 2*j ] << " "
	                            << vector_sol [ 2*j+1 ] << " 0. "<< std::endl;  }
	} // just a block of code

	std::cout << "produced file cell-elast-strain-cell-sq.msh" << std::endl;	
							
}  // end of main


//-----------------------------------------------------------------------------------------


inline bool flip_segment ( Mesh & msh, Cell & seg )

// flip 'seg' if it is inner to 'msh'
// equilibrates (baricenter) the four neighbour vertices

// return true if the segment has been flipped, false if not

// assumes there are only triangular cells
	
// assumes there is no higher-dimensional mesh "above" 'msh'

// assumes that the current manifold is a quotient manifold (manipulates windings)

{	Cell tri2 = msh.cell_in_front_of ( seg, tag::may_not_exist );
	if ( not tri2.exists() ) return false; 
	Cell tri1 = msh.cell_behind ( seg, tag::may_not_exist );
	if ( not tri1.exists() ) return false;
	// or, equivalently :  if ( not seg.is_inner_to ( msh ) ) return false

	Cell A = seg.base().reverse();
	Cell B = seg.tip();

	Cell BC = tri1 .boundary() .cell_in_front_of ( B, tag::surely_exists );
	Cell CA = tri1 .boundary() .cell_behind ( A, tag::surely_exists );
	Cell AD = tri2 .boundary() .cell_in_front_of ( A, tag::surely_exists );
	Cell DB = tri2 .boundary() .cell_behind ( B, tag::surely_exists );
	Cell C = BC.tip();
	assert ( CA.base().reverse() == C );
	Cell D = AD.tip();
	assert ( DB.base().reverse() == D );
	
	assert ( seg.winding() == - BC.winding() - CA.winding() );
	assert ( seg.winding() == AD.winding() + DB.winding() );
	seg.winding() = CA.winding() + AD.winding();
	assert ( seg.winding() == - BC.winding() - DB.winding() );

	B .cut_from_bdry_of ( seg, tag::do_not_bother );
	A .reverse() .cut_from_bdry_of ( seg, tag::do_not_bother );
	CA .cut_from_bdry_of ( tri1, tag::do_not_bother );
	DB .cut_from_bdry_of ( tri2, tag::do_not_bother );
	C .reverse() .glue_on_bdry_of ( seg, tag::do_not_bother );
	D .glue_on_bdry_of ( seg, tag::do_not_bother );
	DB .glue_on_bdry_of ( tri1, tag::do_not_bother );
	CA .glue_on_bdry_of ( tri2, tag::do_not_bother );

	tri1 .boundary() .closed_loop ( B );
	tri2 .boundary() .closed_loop ( A );

	if ( A .is_inner_to ( msh ) ) msh .baricenter ( A, tag::winding );
	if ( B .is_inner_to ( msh ) ) msh .baricenter ( B, tag::winding );
	if ( C .is_inner_to ( msh ) ) msh .baricenter ( C, tag::winding );
	if ( D .is_inner_to ( msh ) ) msh .baricenter ( D, tag::winding );

	return true;

}  // end of  flip_segment
	
//-----------------------------------------------------------------------------------------


double length_square ( Cell AB )
	
{	Manifold space = Manifold::working;
	assert ( space.exists() );  // we use the current (quotient) manifold
	Manifold::Quotient * mani_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space.core );
	Function coords_q = space.coordinates();
	Manifold mani_Eu = mani_q->base_space;  // underlying Euclidian manifold
	Function coords_Eu = mani_Eu.coordinates();
	Cell A = AB.base().reverse();
	Cell B = AB.tip();
	Manifold::Action s = AB.winding();
	std::vector < double > A_co = coords_q ( A );  // same as coords_Eu ( A )
	std::vector < double > B_co = coords_q ( B, tag::winding, s );
	size_t n = A_co.size();
	assert ( n == B_co .size() );
	double len_AB_2 = 0.;
	for ( size_t i = 0; i < n; i++ )
		{	double v = B_co[i] - A_co[i];  len_AB_2 += v*v;  }
	return len_AB_2;                                              }
	

class compare_lenghts_of_segs

{ public :
	inline bool operator() ( Cell AB, Cell CD ) const
	{	return length_square ( AB ) > length_square ( CD );  }
};

//-----------------------------------------------------------------------------------//


inline bool split_segment ( Mesh & msh, Cell & seg )

// splits 'seg' in two if it is inner to 'msh'
// splits also the two neighbour triangles
// equilibrates (baricenter) the four neighbour vertices
// among the four segments "around" the newly created vertex,
// flips the longest

// returns true if the segment has been split, false if not
	
// assumes there is no higher-dimensional mesh "above" 'msh'

// assumes also that the current manifold is not implicit
// (for an implicit manifold, a projection operation should be added)

// assumes that the current manifold is a quotient manifold (manipulates windings)

{	Manifold space = Manifold::working;
	assert ( space.exists() );  // we use the current (quotient) manifold
	Manifold::Quotient * mani_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space.core );
	Function coords_q = space.coordinates();
	Manifold mani_Eu = mani_q->base_space;  // underlying Euclidian manifold
	Function coords_Eu = mani_Eu.coordinates();

	if ( not seg .belongs_to ( msh ) ) return false;

	Cell tri2 = msh.cell_in_front_of ( seg, tag::may_not_exist );
	if ( not tri2.exists() ) return false;
	Cell tri1 = msh.cell_behind ( seg, tag::may_not_exist );
	if ( not tri1.exists() ) return false;
	// or, equivalently :  if ( not seg.is_inner_to ( msh ) ) return false
	// refuse to split segments on the boundary

	Cell A = seg .base() .reverse();
	Cell B = seg .tip();
	Manifold::Action s = seg.winding();
	std::vector < double > A_co = coords_q ( A );  // same as coords_Eu ( A )
	std::vector < double > B_co = coords_q ( B, tag::winding, s );
	size_t n = A_co.size();

	Cell BC = tri1 .boundary() .cell_in_front_of ( B, tag::surely_exists );
	Cell CA = tri1 .boundary() .cell_behind ( A, tag::surely_exists );
	Cell AD = tri2 .boundary() .cell_in_front_of ( A, tag::surely_exists );
	Cell DB = tri2 .boundary() .cell_behind ( B, tag::surely_exists );
	Cell C = BC.tip();
	assert ( CA.base().reverse() == C );
	Cell D = AD.tip();
	assert ( DB.base().reverse() == D );

	assert ( seg.winding() == - BC.winding() - CA.winding() );
	assert ( seg.winding() == AD.winding() + DB.winding() );

	Cell E ( tag::vertex );   // put 'E' at the middle of 'seg' (aka AB)
	for ( size_t i = 0; i < n; i++ )
		coords_Eu[i] ( E ) = ( A_co[i] + B_co[i] ) / 2.;

	A .reverse() .cut_from_bdry_of ( seg, tag::do_not_bother );
	E .reverse() .glue_on_bdry_of ( seg, tag::do_not_bother );
	CA .cut_from_bdry_of ( tri1, tag::do_not_bother );
	AD .cut_from_bdry_of ( tri2, tag::do_not_bother );

	Cell AE ( tag::segment, A.reverse(), E );  // no winding
	Cell CE ( tag::segment, C.reverse(), E );
	CE .winding() = CA .winding();
	Cell DE ( tag::segment, D.reverse(), E );
	DE .winding() = - AD .winding();

	CE .glue_on_bdry_of ( tri1, tag::do_not_bother );
	DE .reverse() .glue_on_bdry_of ( tri2, tag::do_not_bother );

	tri1 .boundary() .closed_loop ( B );
	tri2 .boundary() .closed_loop ( B );

	Cell AEC ( tag::triangle, AE, CE.reverse(), CA );
	AEC .add_to_mesh ( msh );
	Cell ADE ( tag::triangle, AD, DE, AE.reverse() );
	ADE .add_to_mesh ( msh );

	// we want to sweep over the four segments "around" E : CA, AD, DB, BC
	// we want to flip the longest one (could be the two longest)
	// if the longest one is stuck (is on the boundary) we flip the next one

	// we use a map to order them
	compare_lenghts_of_segs comp_len;
	std::multiset < Cell, compare_lenghts_of_segs > ms ( comp_len );
	ms .insert ( CA );
	ms .insert ( AD );
	ms .insert ( DB );
	ms .insert ( BC );

	for ( std::multiset < Cell, compare_lenghts_of_segs > ::iterator it = ms .begin();
				it != ms .end(); it ++                                                      )
	{	Cell sseg = * it;
		if ( flip_segment ( msh, sseg ) )  break;  }
		// we choose to to flip only one segment

	if ( A .is_inner_to ( msh ) ) msh .baricenter ( A, tag::winding );
	if ( B .is_inner_to ( msh ) ) msh .baricenter ( B, tag::winding );
	if ( C .is_inner_to ( msh ) ) msh .baricenter ( C, tag::winding );
	if ( D .is_inner_to ( msh ) ) msh .baricenter ( D, tag::winding );
	if ( E .is_inner_to ( msh ) ) msh .baricenter ( E, tag::winding );

	return true;

}  // end of  split_segment


//-----------------------------------------------------------------------------------------

void limit_number_of_neighbours ( Mesh msh )

// only applies to meshes of triangular cells
	
// assumes there is no higher-dimensional mesh "above" 'msh'

// calls 'flip_segment' which assumes that the current manifold is a quotient manifold
// (manipulates windings)

{	std::forward_list < Cell > has_few_neighbours, has_many_neighbours;

	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = msh .iterator ( tag::over_vertices );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell P = * it;
		if ( not P .is_inner_to ( msh ) ) continue;
		// how many neighbours does P have ?
		size_t counter = 0;
		Mesh::Iterator it_around_P = msh .iterator ( tag::over_segments, tag::around, P );
		for ( it_around_P .reset(); it_around_P .in_range(); it_around_P ++ )
			counter ++;
		if ( counter < 5 ) has_few_neighbours .push_front ( P );
		if ( counter > 7 ) has_many_neighbours .push_front ( P );                         }
	} // just a block of code for hiding 'it'

	// we use a map to order neighbour segments, we flip the longest one
	compare_lenghts_of_segs comp_len;
	for ( std::forward_list < Cell > ::iterator it = has_few_neighbours .begin();
				it != has_few_neighbours .end(); it ++                                 )
	{	Cell P = * it;
		// we count again the neighbours, configuration may have changed in the meanwhile
		size_t counter = 0;
		std::multiset < Cell, compare_lenghts_of_segs > ms ( comp_len );
		Mesh::Iterator it_around_P =
			msh .iterator ( tag::over_cells_of_dim, 2, tag::around, P );
		for ( it_around_P .reset(); it_around_P .in_range(); it_around_P ++ )
		{	Cell tri = * it_around_P;
			counter ++;
			assert ( tri.dim() == 2 );
			Cell SP = tri .boundary() .cell_behind ( P, tag::surely_exists );
			Cell PT = tri .boundary() .cell_in_front_of ( P, tag::surely_exists );
			Cell T = PT .tip();
			Cell TS = tri .boundary() .cell_in_front_of ( T, tag::surely_exists );
			assert ( TS .tip() == SP .base() .reverse() );
			ms .insert ( TS );                                                      }
		if ( counter > 4 ) continue;
		for ( std::multiset < Cell, compare_lenghts_of_segs > ::iterator itt = ms .begin();
					itt != ms .end(); itt ++                                                      )
		{	Cell seg = * itt;
			if ( flip_segment ( msh, seg ) )  break;  }                                         }
			// we choose to to flip only one segment

	for ( std::forward_list < Cell > ::iterator it = has_many_neighbours .begin();
				it != has_many_neighbours .end(); it ++                                 )
	{	Cell P = * it;
		// we count again the neighbours, configuration may have changed in the meanwhile
		size_t counter = 0;
		std::multiset < Cell, compare_lenghts_of_segs > ms ( comp_len );
		Mesh::Iterator it_around_P =
			msh .iterator ( tag::over_segments, tag::around, P );
		for ( it_around_P .reset(); it_around_P .in_range(); it_around_P ++ )
		{	Cell seg = * it_around_P;
			assert ( seg.tip() == P );
			counter ++;
			ms .insert ( seg );        }
		if ( counter < 8 ) continue;
		for ( std::multiset < Cell, compare_lenghts_of_segs > ::iterator itt = ms .begin();
					itt != ms .end(); itt ++                                                      )
		{	Cell seg = * itt;
			if ( flip_segment ( msh, seg ) )  break;  }                                         }
			// we choose to to flip only one segment

}  // end of  limit_number_of_neighbours
	
//-----------------------------------------------------------------------------------------

	
void flip_split_long_segments ( Mesh & msh, double threshold )

// flip long segments, make baricenters on the four neighbour vertices
// if a segment is still long after flip, split it

// assumes there is no higher-dimensional mesh "above" 'msh'

// assumes also that the current manifold is not implicit
// (for an implicit manifold, a projection operation should be added)

// calls 'flip_segment' and 'split_segment' which assume that the current manifold
// is a quotient manifold (manipulate windings)

{
	double thr_sq = threshold * threshold;

	std::list < Cell > list_of_segments;

	Mesh::Iterator it = msh.iterator ( tag::over_segments );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell seg = *it;

	  Cell tri1 = msh.cell_in_front_of ( seg, tag::may_not_exist );
		if ( not tri1.exists() ) continue;
	  Cell tri2 = msh.cell_behind ( seg, tag::may_not_exist );
		if ( not tri2.exists() ) continue;
		// or, equivalently :  if ( not seg.is_inner_to ( msh ) ) continue
		// segments on the boundary cannot be flipped

		if ( length_square ( seg ) > thr_sq )  list_of_segments .push_back ( seg );   }

	for ( size_t ii = 0; ii < 2; ii++ )
	for ( std::list < Cell > ::iterator itt = list_of_segments .begin();
        itt != list_of_segments .end();                                )
	{	Cell seg = * itt;

		// 'seg' may have been eliminated from 'msh' in the meanwhile
		if ( not seg .belongs_to ( msh ) )
			{	itt = list_of_segments .erase ( itt );  continue;  }

		// the configuration around 'seg' may have changed in the meanwhile
		Cell tri2 = msh.cell_in_front_of ( seg, tag::may_not_exist );
		if ( not tri2.exists() ) 
		{	itt = list_of_segments .erase ( itt );  continue;  }
		Cell tri1 = msh.cell_behind ( seg, tag::may_not_exist );
		if ( not tri1.exists() )
		{	itt = list_of_segments .erase ( itt );  continue;  }
		// or, equivalently :  if ( not seg.is_inner_to ( msh ) ) continue
		// segments on the boundary cannot be flipped

		// length of 'seg' may have changed in the meanwhile
		if ( length_square ( seg ) < thr_sq ) 
		{	itt = list_of_segments .erase ( itt );  continue;  }

		flip_segment ( msh, seg );
		itt ++;                                                                   }

	// segments still on the list will be split

	for ( std::list < Cell > ::iterator itt = list_of_segments .begin();
        itt != list_of_segments .end(); itt ++                         )
	{	Cell seg = * itt;
		// the configuration around 'seg' may have changed in the meanwhile
	  Cell tri1 = msh.cell_in_front_of ( seg, tag::may_not_exist );
		if ( not tri1.exists() ) continue;
	  Cell tri2 = msh.cell_behind ( seg, tag::may_not_exist );
		if ( not tri2.exists() ) continue;
		// or, equivalently :  if ( not seg.is_inner_to ( msh ) ) continue
		if ( length_square ( seg ) > thr_sq )  split_segment ( msh, seg );  }

}

//-----------------------------------------------------------------------------------------


void remove_short_segments ( Mesh & msh, double threshold )

// remove segments shorter than the given threshold

// assumes there is no higher-dimensional mesh "above" 'msh'

// assumes also that the current manifold is not implicit
// (for an implicit manifold, a projection operation should be added)

// assumes that the current manifold is a quotient manifold
// (it manipulates windings)

{	Manifold space = Manifold::working;
	assert ( space.exists() );  // we use the current (quotient) manifold
	Manifold::Quotient * mani_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space.core );
	Function coords_q = space.coordinates();
	Manifold mani_Eu = mani_q->base_space;  // underlying Euclidian manifold
	Function coords_Eu = mani_Eu.coordinates();

	double thr_sq = threshold * threshold;

	std::forward_list < Cell > list_of_segments;

	Mesh::Iterator it = msh.iterator ( tag::over_segments );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell seg = *it;

	  Cell tri1 = msh.cell_in_front_of ( seg, tag::may_not_exist );
		if ( not tri1.exists() ) continue;
	  Cell tri2 = msh.cell_behind ( seg, tag::may_not_exist );
		if ( not tri2.exists() ) continue;
		// or, equivalently :  if ( not seg.is_inner_to ( msh ) ) continue
		// segments on the boundary are left unchanged

		Cell A = seg.base().reverse();
		Cell B = seg.tip();
		if ( length_square ( seg ) > thr_sq ) continue;
		// long segments are left unchanged

		// method 'is_inner_to' says "vertex is not on the boundary of mesh"
		bool wing_A_free = A .is_inner_to (msh );
		bool wing_B_free = B .is_inner_to (msh );

		if ( wing_A_free or wing_B_free ) list_of_segments .push_front ( seg );  }

	for ( std::forward_list < Cell > ::iterator	itt = list_of_segments .begin();
	      itt != list_of_segments .end(); itt ++                                 )
	{	Cell seg = *itt;
		assert ( seg.exists() );

		// 'seg' may have been already eliminated
		if ( not seg .belongs_to ( msh ) ) continue;
		
	  Cell tri1 = msh.cell_in_front_of ( seg, tag::may_not_exist );
		if ( not tri1.exists() ) continue;
	  Cell tri2 = msh.cell_behind ( seg, tag::may_not_exist );
		if ( not tri2.exists() ) continue;
		// or, equivalently :  if ( not seg.is_inner_to ( msh ) ) continue
		// segments on the boundary are left unchanged

		Cell A = seg.base().reverse();
		Cell B = seg.tip();
		Manifold::Action s = seg.winding();
		std::vector < double > A_co = coords_q ( A );
		std::vector < double > B_co = coords_q ( B, tag::winding, s );
		size_t n = A_co.size();

		Cell CB = tri1 .boundary() .cell_behind ( B, tag::surely_exists );
		Cell AC = tri1 .boundary() .cell_in_front_of ( A, tag::surely_exists );
		Cell DA = tri2 .boundary() .cell_behind ( A, tag::surely_exists );
		Cell BD = tri2 .boundary() .cell_in_front_of ( B, tag::surely_exists );
		assert ( CB.base().reverse() == AC.tip() );
		assert ( DA.base().reverse() == BD.tip() );
		Cell tri3 = msh .cell_in_front_of ( CB, tag::may_not_exist );
		Cell tri4 = msh .cell_in_front_of ( BD, tag::may_not_exist );
		Cell tri5 = msh .cell_in_front_of ( AC, tag::may_not_exist );
		Cell tri6 = msh .cell_in_front_of ( DA, tag::may_not_exist );

		// method 'is_inner_to' says "vertex is not on the boundary of mesh"
		bool wing_A_free = A .is_inner_to (msh );
		bool wing_B_free = B .is_inner_to (msh );

		// due to previous changes, wings may be not free any longer
		if ( ( not wing_A_free ) and ( not wing_B_free ) ) continue;

		// beware, if msh happens to be the boundary of some higher-dimensional cell,
		// and that higher-dim cell has neighbours in some higher-dim mesh,
		// statements below do not work as expected

		if ( wing_B_free and wing_A_free )
			// both wings free, we collapse one of them at random
			// we collapse wing_B, vertex B disappears,
			// A will occupy middle of segment seg

		{	// make a list of segments pointing towards B
			std::forward_list < Cell > list_of_segs;
			Mesh::Iterator it_around_B = msh.iterator ( tag::over_segments, tag::around, B );
			it_around_B .reset ( tag::start_at, BD.reverse() );
			assert ( it_around_B .in_range() );
			assert ( *it_around_B == BD.reverse() );
			it_around_B++;
			assert ( it_around_B .in_range() );
			assert ( *it_around_B == seg );
			it_around_B++;
			assert ( it_around_B .in_range() );
			assert ( *it_around_B == CB );
			// seg, BD and CB will disappear
			for ( it_around_B++; it_around_B .in_range(); it_around_B++ )
				list_of_segs .push_front ( *it_around_B );
			
			// change winding number of these segments
			if ( s != 0 )
				for ( std::forward_list < Cell > ::iterator it_list = list_of_segs .begin();
			        it_list != list_of_segs .end(); it_list ++                            )
				{	Cell segm = *it_list;      // 'segm' points towards B
				  segm.winding() -= s;     }

			tri1 .remove_from_mesh ( msh );
			tri2 .remove_from_mesh ( msh );
			assert ( tri3.exists() );
			assert ( tri4.exists() );
			assert ( tri5.exists() );
			assert ( tri6.exists() );
			CB .reverse() .cut_from_bdry_of ( tri3, tag::do_not_bother );
			BD .reverse() .cut_from_bdry_of ( tri4, tag::do_not_bother );
			AC .glue_on_bdry_of ( tri3, tag::do_not_bother );
			DA .glue_on_bdry_of ( tri4, tag::do_not_bother );

			// replace B by A in segments neighbour to B
			for ( std::forward_list < Cell > ::iterator it_list = list_of_segs .begin();
			      it_list != list_of_segs .end(); it_list ++                            )
			{	Cell segm = *it_list;
				assert ( segm.tip() == B );
				B .cut_from_bdry_of ( segm, tag::do_not_bother );
				A .glue_on_bdry_of ( segm, tag::do_not_bother );  }

			tri3 .boundary() .closed_loop ( A );
			tri4 .boundary() .closed_loop ( A );

			// change coordinates of A
			// if we are on an implicit manifold, we should project A
			for ( size_t i = 0; i < n; i++ ) A_co[i] = ( A_co[i] + B_co[i] ) / 2.;
			coords_Eu ( A ) = A_co;                                                        }

		else if ( wing_B_free )
			// we collapse wing_B, vertex B disappears, vertex A stays where it is

		{	// make a list of segments pointing towards B
			std::forward_list < Cell > list_of_segs;
			Mesh::Iterator it_around_B = msh.iterator ( tag::over_segments, tag::around, B );
			it_around_B .reset ( tag::start_at, BD.reverse() );
			assert ( it_around_B .in_range() );
			assert ( *it_around_B == BD.reverse() );
			it_around_B++;
			assert ( it_around_B .in_range() );
			assert ( *it_around_B == seg );
			it_around_B++;
			assert ( it_around_B .in_range() );
			assert ( *it_around_B == CB );
			// seg, BD and CB will disappear
			for ( it_around_B++; it_around_B .in_range(); it_around_B++ )
				list_of_segs .push_front ( *it_around_B );
			
			// change winding number of these segments
			if ( s != 0 )
				for ( std::forward_list < Cell > ::iterator it_list = list_of_segs .begin();
			        it_list != list_of_segs .end(); it_list ++                            )
				{	Cell segm = *it_list;      // 'segm' points towards B
				  segm.winding() -= s;     }

			tri1 .remove_from_mesh ( msh );
			tri2 .remove_from_mesh ( msh );
			assert ( tri3.exists() );
			assert ( tri4.exists() );
			CB .reverse() .cut_from_bdry_of ( tri3, tag::do_not_bother );
			BD .reverse() .cut_from_bdry_of ( tri4, tag::do_not_bother );
			AC .glue_on_bdry_of ( tri3, tag::do_not_bother );
			DA .glue_on_bdry_of ( tri4, tag::do_not_bother );

			// replace B by A in segments neighbour to B
			for ( std::forward_list < Cell > ::iterator it_list = list_of_segs .begin();
			      it_list != list_of_segs .end(); it_list ++                            )
			{	Cell segm = *it_list;
				assert ( segm.tip() == B );
				assert ( B.exists() );  assert ( segm.exists() );
				B .cut_from_bdry_of ( segm, tag::do_not_bother );
				A .glue_on_bdry_of ( segm, tag::do_not_bother );  }

			tri3 .boundary() .closed_loop ( A );
			tri4 .boundary() .closed_loop ( A );                                           }
		
		else if ( wing_A_free )
			// we collapse wing_A, vertex A disappears,, vertex B stays where it is

		{	// make a list of segments pointing towards A
			std::forward_list < Cell > list_of_segs;
			Mesh::Iterator it_around_A = msh.iterator ( tag::over_segments, tag::around, A );
			it_around_A .reset ( tag::start_at, AC.reverse() );
			assert ( it_around_A .in_range() );
			assert ( *it_around_A == AC.reverse() );
			it_around_A++;
			assert ( it_around_A .in_range() );
			assert ( *it_around_A == seg.reverse() );
			it_around_A++;
			assert ( it_around_A .in_range() );
			assert ( *it_around_A == DA );
			// seg, AC and DA will disappear
			for ( it_around_A++; it_around_A .in_range(); it_around_A++ )
				list_of_segs .push_front ( *it_around_A );
			
			// change winding number of these segments
			if ( s != 0 )
				for ( std::forward_list < Cell > ::iterator it_list = list_of_segs .begin();
			        it_list != list_of_segs .end(); it_list ++                            )
				{	Cell segm = *it_list;      // 'segm' points towards A
				  segm.winding() += s;     }

			tri1 .remove_from_mesh ( msh );
			tri2 .remove_from_mesh ( msh );
			assert ( tri5.exists() );
			assert ( tri6.exists() );
			AC .reverse() .cut_from_bdry_of ( tri5, tag::do_not_bother );
			DA .reverse() .cut_from_bdry_of ( tri6, tag::do_not_bother );
			CB .glue_on_bdry_of ( tri5, tag::do_not_bother );
			BD .glue_on_bdry_of ( tri6, tag::do_not_bother );

			// replace A by B in segments neighbour to A
			for ( std::forward_list < Cell > ::iterator it_list = list_of_segs .begin();
			      it_list != list_of_segs .end(); it_list ++                            )
			{	Cell segm = *it_list;
				assert ( segm.tip() == A );
				A .cut_from_bdry_of ( segm, tag::do_not_bother );
				B .glue_on_bdry_of ( segm, tag::do_not_bother );  }

			tri5 .boundary() .closed_loop ( B );
			tri6 .boundary() .closed_loop ( B );                                          }

		else assert ( false );    // at least one wing should be free
		
	}  // end of  for  over segments of msh
	
}  // end of  remove_short_segments
	

void baricenters ( Mesh & msh )

// segments on the boundary are not moved -- we use method  is_inner_to

{	Mesh::Iterator it = msh.iterator ( tag::over_vertices );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell P = *it;
		if ( P.is_inner_to ( msh ) ) msh.baricenter ( P, tag::winding );  }  }
	
