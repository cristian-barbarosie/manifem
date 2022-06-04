
// 2D cantilever

#include "maniFEM.h"
#include <fstream>
#include <set>
#include <Eigen/Sparse>
#include <Eigen/OrderingMethods>
using namespace maniFEM;


void impose_value_of_unknown
(	Eigen::SparseMatrix < double > & matrix_A, Eigen::VectorXd & vector_b,
	size_t i, double val                                                   )

// in a system of linear equations, destroy equation 'i' and impose u(i) = val
// change also column 'i' of the matrix, just to preserve symmetry

// used for imposing Dirichlet boundary conditions

{	size_t size_matrix = matrix_A .innerSize();
	vector_b (i) = val;
	for ( size_t j = 0; j < size_matrix; j++ )
		matrix_A .coeffRef ( i, j ) = 0.;
	matrix_A .coeffRef ( i, i ) = 1.;
	for ( size_t j = 0; j < size_matrix; j++ )
	{	if ( i == j ) continue;
		vector_b(j) -= matrix_A .coeff ( j, i ) * val;
		matrix_A .coeffRef ( j, i ) = 0.;                 }  }


inline std::vector < double > tangent ( const Cell & seg )
{	Function xy = Manifold::working .coordinates();
	Function x = xy [0], y = xy [1];
	Cell A = seg .base() .reverse();
	Cell B = seg .tip();
	std::vector < double > t (2);
	assert ( t .size() == 2 );
	t[0] = x(B) - x(A);
	t[1] = y(B) - y(A);
	return t;                                         }


inline double sq_length ( const Cell & seg )
{	std::vector < double > t = tangent ( seg );
	double l = 0.;
	for ( size_t i = 0; i < 2; i++ ) l += t[i]*t[i];
	return l;                                        }

		
double area ( Cell tri )

{	assert ( tri .dim() == 2 );
	Function xy = Manifold::working .coordinates();
	Function x = xy [0], y = xy [1];
	Mesh::Iterator it = tri .boundary() .iterator ( tag::over_vertices );
	it .reset();  assert ( it .in_range() );
	Cell A = *it;
	it ++;  assert ( it .in_range() );
	Cell B = *it;
	it ++;  assert ( it .in_range() );
	Cell C = *it;
	it ++;  assert ( not it .in_range() );
	double ABx = x(B) - x(A), ABy = y(B)-y(A), ACx = x(C) - x(A), ACy = y(C) - y(A);
	double a = ABx * ACy - ABy * ACx;
	assert ( a > 0. );
	return a;                                                                         }

void split_elim_segs ( Mesh & chain, const double lower_lim, const double upper_lim );
// defined after main

void smoothen ( Mesh & chain );
// defined after main


int main1 ( )  // produce initial shape

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy [0], y = xy [1];

	Cell A ( tag::vertex, tag::of_coords, { 0., 0. } );
	Cell B ( tag::vertex, tag::of_coords, { 0.4, 0. } );
	Cell C ( tag::vertex, tag::of_coords, { 0.4, 1. } );
	Cell D ( tag::vertex, tag::of_coords, { 0., 1. } );
	Cell E ( tag::vertex, tag::of_coords, { 0.4, 0.4 } );
	Cell F ( tag::vertex, tag::of_coords, { 0.4, 0.6 } );

	const double seg_len = 0.03;
	const size_t n = 1. / seg_len;

	Mesh BE ( tag::segment, B.reverse(), E, tag::divided_in, int(0.4*n) );
	Mesh EF ( tag::segment, E.reverse(), F, tag::divided_in, int(0.2*n) );
	Mesh FC ( tag::segment, F.reverse(), C, tag::divided_in, int(0.4*n) );

	// Manifold arc_inf = RR2 .implicit ( y == 0.15 * x * (x-2.) );
	Mesh AB_old ( tag::frontal, tag::start_at, A, tag::towards, { 1., 0. },
	                        tag::stop_at,  B, tag::desired_length, seg_len );
	Mesh AB ( tag::segment, A .reverse(), B, tag::divided_in, int(0.3*n) );

	// Manifold arc_sup = RR2 .implicit ( y == 1. - 0.15 * x * (x-2.) );
	Mesh CD_old ( tag::frontal, tag::start_at, C, tag::towards, { -1., 0. },
	                        tag::stop_at,  D, tag::desired_length, seg_len );
	Mesh CD ( tag::segment, C .reverse(), D, tag::divided_in, int(0.3*n) );

	Manifold clamped = RR2 .implicit ( x == 0. );
	Manifold circle_manif_clamped = RR2 .implicit ( x*x + (y-0.5)*(y-0.5) == 0.03 );
	Manifold two_points ( tag::intersect, circle_manif_clamped, clamped );
	Cell G ( tag::vertex, tag::of_coords, { 0., 0.9 }, tag::project );
	Cell H ( tag::vertex, tag::of_coords, { 0., 0.1 }, tag::project );

	clamped .set_as_working_manifold();
	Mesh DG ( tag::frontal, tag::start_at, D, tag::towards, { 0., -1. },
	                        tag::stop_at,  G, tag::desired_length, seg_len );
	Mesh HA ( tag::frontal, tag::start_at, H, tag::towards, { 0., -1. },
	                        tag::stop_at,  A, tag::desired_length, seg_len );

	circle_manif_clamped .set_as_working_manifold();
	Mesh semi_circle ( tag::frontal, tag::start_at, G, tag::towards, { 1., 0. },
	                                 tag::stop_at,  H, tag::desired_length, seg_len );
	
	Manifold circle_manif_1 = RR2 .implicit ( (x-0.4)*(x-0.4) + (y-0.75)*(y-0.75) == 0.025 );
	Mesh circle_1 ( tag::frontal, tag::desired_length, seg_len );

	Manifold circle_manif_2 = RR2 .implicit ( (x-0.4)*(x-0.4) + (y-0.25)*(y-0.25) == 0.025 );
	Mesh circle_2 ( tag::frontal, tag::desired_length, seg_len );

	Manifold circle_manif_3 = RR2 .implicit ( (x-0.9)*(x-0.9) + (y-0.95)*(y-0.95) == 0.02 );
	Mesh circle_3 ( tag::frontal, tag::desired_length, seg_len );

	Manifold circle_manif_4 = RR2 .implicit ( (x-0.9)*(x-0.9) + (y-0.5)*(y-0.5) == 0.025 );
	Mesh circle_4 ( tag::frontal, tag::desired_length, seg_len );

	Manifold circle_manif_5 = RR2 .implicit ( (x-0.9)*(x-0.9) + (y-0.05)*(y-0.05) == 0.02 );
	Mesh circle_5 ( tag::frontal, tag::desired_length, seg_len );

	Manifold circle_manif_6 = RR2 .implicit ( (x-1.3)*(x-1.3) + (y-0.75)*(y-0.75) == 0.02 );
	Mesh circle_6 ( tag::frontal, tag::desired_length, seg_len );

	Manifold circle_manif_7 = RR2 .implicit ( (x-1.3)*(x-1.3) + (y-0.25)*(y-0.25) == 0.02 );
	Mesh circle_7 ( tag::frontal, tag::desired_length, seg_len );

	Manifold circle_manif_8 = RR2 .implicit ( (x-1.7)*(x-1.7) + (y-0.5)*(y-0.5) == 0.025 );
	Mesh circle_8 ( tag::frontal, tag::desired_length, seg_len );

	std::list < Mesh > inner_list_old
		{ circle_1 .reverse(), circle_2 .reverse(), circle_4 .reverse(),
		  circle_6 .reverse(), circle_7 .reverse(), circle_8 .reverse() };
	std::list < Mesh > inner_list;

	std::list < Mesh > bdry_list = inner_list;
	bdry_list .push_back ( AB );
	bdry_list .push_back ( BE );
	bdry_list .push_back ( EF );
	bdry_list .push_back ( FC );
	bdry_list .push_back ( CD );
	bdry_list .push_back ( DG );
	bdry_list .push_back ( semi_circle );
	bdry_list .push_back ( HA );

	Mesh bdry ( tag::join, bdry_list );
	RR2 .set_as_working_manifold();
	// Mesh cantilever ( tag::frontal, tag::boundary, bdry, tag::desired_length, seg_len );
	// cantilever .draw_ps ("cantilever.eps");
	// exit(0);

	// write to file coordinates of boundary vertices
	std::ofstream f ("bdry.dat");

	f << seg_len << std::endl;

	// begin at lower left corner
	{ // just a block of code for hiding 'it'
	f << AB .number_of ( tag::segments ) << std::endl;
	Mesh::Iterator it = AB .iterator ( tag::over_vertices, tag::require_order );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;
		f << x(V) << " " << y(V) << std::endl;  }
	} // just a block of code

	// three pieces on the right hand side
	{ // just a block of code for hiding 'it'
	f << BE .number_of ( tag::segments ) << std::endl;
	Mesh::Iterator it = BE .iterator ( tag::over_vertices, tag::require_order );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;
		f << x(V) << " " << y(V) << std::endl;  }
	} { // just a block of code for hiding 'it'
	f << EF .number_of ( tag::segments ) << std::endl;
	Mesh::Iterator it = EF .iterator ( tag::over_vertices, tag::require_order );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;
		f << x(V) << " " << y(V) << std::endl;  }
	} { // just a block of code for hiding 'it'
	f << FC .number_of ( tag::segments ) << std::endl;
	Mesh::Iterator it = FC .iterator ( tag::over_vertices, tag::require_order );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;
		f << x(V) << " " << y(V) << std::endl;  }
	} // just a block of code
	
	// upper side
	{ // just a block of code for hiding 'it'
	f << CD .number_of ( tag::segments ) << std::endl;
	Mesh::Iterator it = CD .iterator ( tag::over_vertices, tag::require_order );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;
		f << x(V) << " " << y(V) << std::endl;  }
	} // just a block of code

	// three pieces on the left hand side
	{ // just a block of code for hiding 'it'
	f << DG .number_of ( tag::segments ) << std::endl;
	Mesh::Iterator it = DG .iterator ( tag::over_vertices, tag::require_order );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;
		f << x(V) << " " << y(V) << std::endl;  }
	} { // just a block of code for hiding 'it'
	f << semi_circle .number_of ( tag::segments ) << std::endl;
	Mesh::Iterator it = semi_circle .iterator ( tag::over_vertices, tag::require_order );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;
		f << x(V) << " " << y(V) << std::endl;  }
	} { // just a block of code for hiding 'it'
	f << HA .number_of ( tag::segments ) << std::endl;
	Mesh::Iterator it = HA .iterator ( tag::over_vertices, tag::require_order );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;
		f << x(V) << " " << y(V) << std::endl;  }
	} // just a block of code
	
	// several inner curves
	f << inner_list .size() << std::endl;
	for ( std::list < Mesh > ::iterator il = inner_list .begin();
				il != inner_list .end(); il++ )
	{	Mesh circle = *il;
		f << circle .number_of ( tag::segments ) << std::endl;
		Mesh::Iterator it = circle .iterator ( tag::over_vertices, tag::require_order );
		for ( it .reset(); it .in_range(); it++ )
		{	Cell V = *it;
			f << x(V) << " " << y(V) << std::endl;  }                                    }

	return 0;
	
}  // end of  main1
	
//-----------------------------------------------------------------------------------------------------//	


int main ()
	
{	const double eta = 0.00005;
	const double lagrange = 100.;
	const size_t iter_max = 5;

	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy [0], y = xy [1];

	// declare two types of finite elements
	// FiniteElement fe_1d ( tag::segment, tag::Lagrange, tag::of_degree, 1 );
	// Integrator integ_1d = fe_1d .set_integrator ( tag::hand_coded );
	// no 1d finite hand coded finite element for now - we will simply apply forces at vertices
	
	FiniteElement fe_2d ( tag::triangle, tag::Lagrange, tag::of_degree, 1 );
	Integrator integ_2d = fe_2d .set_integrator ( tag::hand_coded );

	// read from file coordinates of boundary vertices
	// file will be overwritten at the end !
	std::ifstream ff ("bdry.dat");
	size_t n;
	double xx, yy;
	
	double seg_len;  ff >> seg_len;  // seg_len *= 0.9;

	// begin at lower left corner
	ff >> n;
	Cell A ( tag::vertex );
	Cell B ( tag::vertex );
	Mesh AB ( tag::segment, A .reverse(), B, tag::divided_in, n );
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = AB .iterator ( tag::over_vertices, tag::require_order );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;
		ff >> xx >> yy;  x(V) = xx;  y(V) = yy;   }
	} // just a block of code

	// three pieces on the right hand side
	ff >> n;
	ff >> xx >> yy;  // coords of B, previously read
	Cell E ( tag::vertex );
	Mesh BE ( tag::segment, B .reverse(), E, tag::divided_in, n );
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = BE .iterator ( tag::over_vertices, tag::require_order );
	it .reset();  assert ( it .in_range() );
	for ( it++; it .in_range(); it++ )
	{	Cell V = *it;
		ff >> xx >> yy;  x(V) = xx;  y(V) = yy;   }
	} // just a block of code
	ff >> n;
	ff >> xx >> yy;  // coords of E, previously read
	Cell F ( tag::vertex );
	Mesh EF ( tag::segment, E .reverse(), F, tag::divided_in, n );
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = EF .iterator ( tag::over_vertices, tag::require_order );
	it .reset();  assert ( it .in_range() );
	for ( it++; it .in_range(); it++ )
	{	Cell V = *it;
		ff >> xx >> yy;  x(V) = xx;  y(V) = yy;   }
	} // just a block of code
	ff >> n;
	ff >> xx >> yy;  // coords of F, previously read
	Cell C ( tag::vertex );
	Mesh FC ( tag::segment, F .reverse(), C, tag::divided_in, n );
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = FC .iterator ( tag::over_vertices, tag::require_order );
	it .reset();  assert ( it .in_range() );
	for ( it++; it .in_range(); it++ )
	{	Cell V = *it;
		ff >> xx >> yy;  x(V) = xx;  y(V) = yy;   }
	} // just a block of code

	// upper side
	ff >> n;
	ff >> xx >> yy;  // coords of C, previously read
	Cell D ( tag::vertex );
	Mesh CD ( tag::segment, C .reverse(), D, tag::divided_in, n );
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = CD .iterator ( tag::over_vertices, tag::require_order );
	it .reset();  assert ( it .in_range() );
	for ( it++; it .in_range(); it++ )
	{	Cell V = *it;
		ff >> xx >> yy;  x(V) = xx;  y(V) = yy;   }
	} // just a block of code

	// three pieces on the left hand side
	ff >> n;
	ff >> xx >> yy;  // coords of D, previously read
	Cell G ( tag::vertex );
	Mesh DG ( tag::segment, D .reverse(), G, tag::divided_in, n );
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = DG .iterator ( tag::over_vertices, tag::require_order );
	it .reset();  assert ( it .in_range() );
	for ( it++; it .in_range(); it++ )
	{	Cell V = *it;
		ff >> xx >> yy;  x(V) = xx;  y(V) = yy;   }
	} // just a block of code
	ff >> n;
	ff >> xx >> yy;  // coords of G, previously read
	Cell H ( tag::vertex );
	Mesh semi_circle ( tag::segment, G .reverse(), H, tag::divided_in, n );
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = semi_circle .iterator ( tag::over_vertices, tag::require_order );
	it .reset();  assert ( it .in_range() );
	for ( it++; it .in_range(); it++ )
	{	Cell V = *it;
		ff >> xx >> yy;  x(V) = xx;  y(V) = yy;   }
	} // just a block of code
	ff >> n;
	ff >> xx >> yy;  // coords of H, previously read
	Mesh HA ( tag::segment, H .reverse(), A, tag::divided_in, n );
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = HA .iterator ( tag::over_vertices, tag::require_order );
	it .reset();  assert ( it .in_range() );
	for ( it++; it .in_range(); it++ )
	{	Cell V = *it;
		ff >> xx >> yy;  x(V) = xx;  y(V) = yy;   }
	} // just a block of code

	// several inner circles
	std::list < Mesh > inner;
	{ // just a block of code for hiding names
	size_t m;  ff >> m;
	Manifold circle_manif = RR2 .implicit ( x*x + y*y == 1. );
	for ( size_t i = 0; i < m; i++ )
	{	ff >> n;
		Cell V ( tag::vertex ), W ( tag::vertex );
		Cell WV ( tag::segment, W .reverse(), V );
		Mesh circle ( tag::segment, V .reverse(), W, tag::divided_in, n-1 );
		WV .add_to_mesh ( circle );
		assert ( n == circle .number_of ( tag::segments ) );
		Mesh::Iterator it = circle .iterator ( tag::over_vertices, tag::require_order );
		for ( it .reset(); it .in_range(); it++ )
		{	Cell P = *it;
			ff >> xx >> yy;  x(P) = xx;  y(P) = yy;   }
		inner .push_back ( circle );                                                      }
	assert ( m == inner .size() );
	} // just a block of code

	std::list < Mesh > bdry_pieces { AB, BE, EF, FC, CD, DG, semi_circle, HA };
	std::list < Mesh > free_bdry_pieces { AB, BE, FC, CD, semi_circle };
	std::list < Mesh > Dirichlet_bdry_pieces { DG, HA };

	for ( std::list < Mesh > ::iterator il = inner .begin();
				il != inner .end(); il++ )
	{	Mesh circle = *il;
		bdry_pieces .push_back ( circle );
		free_bdry_pieces .push_back ( circle );  }
		
	Mesh bdry ( tag::join, bdry_pieces );
	Mesh free_bdry ( tag::join, free_bdry_pieces );

	RR2 .set_as_working_manifold();
	Mesh cantilever ( tag::frontal, tag::boundary, bdry, tag::desired_length, seg_len );

	// Hooke's Law
	// sigma = 2 mu e + lambda tr e I
	double lambda = 1., mu = 3.;

	tag::Util::Tensor <double> Hooke(2,2,2,2);
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
	
	// hand-coded integrators require an early declaration of
	// the integrals we intend to compute later (after docking 'fe_2d' on a cell)
	Function bf1 ( tag::basis_function, tag::within, fe_2d ),
	         bf2 ( tag::basis_function, tag::within, fe_2d );
	integ_2d .pre_compute ( tag::for_given, tag::basis_functions, bf1, bf2,
	     tag::integral_of, { bf1 .deriv(x) * bf2 .deriv(x), bf1 .deriv(x) * bf2 .deriv(y),
	                         bf1 .deriv(y) * bf2 .deriv(x), bf1 .deriv(y) * bf2 .deriv(y) } );

	for ( size_t iter = 0; iter < iter_max; iter++ )
	{ std::cout << "iter " << iter << ", ";
		
	std::map < Cell::Core *, size_t > numbering;
	{ // just a block of code for hiding 'it' and 'counter'
	Mesh::Iterator it = cantilever .iterator ( tag::over_vertices );
	size_t counter = 0;
	for ( it.reset() ; it.in_range(); it++ )
	{	Cell p = *it;  ++counter;  numbering [p.core] = counter;  }
	} // just a block of code 
	// in each node P with number i>=0 ,
	// u_x(P)=vector_sol(i) , u_y(P)=vector_sol(i+N) where N is the number of nodes
	//  i=numbering(P)-1	
	size_t N = cantilever .number_of ( tag::vertices );
	assert ( N == numbering .size() );
	size_t size_matrix = 2*N;
	std::cout << "global matrix " << size_matrix << "x" << size_matrix << std::endl;
	Eigen::SparseMatrix <double> matrix_A ( size_matrix, size_matrix );
	matrix_A .reserve ( Eigen::VectorXi::Constant ( size_matrix, 14 ) );
	// since we will be working with a mesh of triangles,
	// there will be, in average, 14 = 2*(6+1) non-zero elements per column
	// the diagonal entry plus six neighbour vertices
	Eigen::VectorXd vector_b ( size_matrix );
	vector_b .setZero();

	// run over all triangular cells composing 'cantilever'
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = cantilever .iterator ( tag::over, tag::cells_of_max_dim );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell small_tri = *it;
		fe_2d .dock_on ( small_tri );
		// run twice over the four vertices of 'small_tri'
		Mesh::Iterator it1 = small_tri .boundary().iterator ( tag::over_vertices );
		Mesh::Iterator it2 = small_tri .boundary().iterator ( tag::over_vertices );
		for ( it1 .reset(); it1 .in_range(); it1++ )
		{	Cell V = *it1 ;
			Function psiV = fe_2d .basis_function(V) ;
			for ( it2 .reset(); it2 .in_range(); it2++ )
			{	Cell W = *it2;  // V may be the same as W, no problem about that
				Function psiW = fe_2d .basis_function (W);
				std::vector < double > result = fe_2d .integrate 
					( tag::pre_computed, tag::replace, bf1, tag::by, psiW,
				                       tag::replace, bf2, tag::by, psiV );
				assert ( result .size() == 4 );
				double int_d_psiW_dx_d_psiV_dx = result[0];
				double int_d_psiW_dx_d_psiV_dy = result[1];
				double int_d_psiW_dy_d_psiV_dx = result[2];
				double int_d_psiW_dy_d_psiV_dy = result[3];
				// 'fe' is already docked on 'small_tri' so this will be the domain of integration
				matrix_A .coeffRef ( numbering[V.core]-1, numbering[W.core]-1 ) += 
				   Hooke (0,0,0,0) * int_d_psiW_dx_d_psiV_dx +
				   Hooke (0,0,0,1) * int_d_psiW_dy_d_psiV_dx +
				   Hooke (0,1,0,0) * int_d_psiW_dx_d_psiV_dy +
				   Hooke (0,1,0,1) * int_d_psiW_dy_d_psiV_dy ;
				matrix_A .coeffRef ( numbering[V.core]-1, numbering[W.core]-1+N ) += 
				   Hooke (0,0,1,0) * int_d_psiW_dx_d_psiV_dx +
				   Hooke (0,0,1,1) * int_d_psiW_dy_d_psiV_dx +
				   Hooke (0,1,1,0) * int_d_psiW_dx_d_psiV_dy +
				   Hooke (0,1,1,1) * int_d_psiW_dy_d_psiV_dy ;
				matrix_A .coeffRef ( numbering[V.core]-1+N, numbering[W.core]-1 ) += 
				   Hooke (1,0,0,0) * int_d_psiW_dx_d_psiV_dx +
				   Hooke (1,0,0,1) * int_d_psiW_dy_d_psiV_dx +
				   Hooke (1,1,0,0) * int_d_psiW_dx_d_psiV_dy +
				   Hooke (1,1,0,1) * int_d_psiW_dy_d_psiV_dy ;
				matrix_A .coeffRef ( numbering[V.core]-1+N, numbering[W.core]-1+N ) += 
				   Hooke (1,0,1,0) * int_d_psiW_dx_d_psiV_dx +
				   Hooke (1,0,1,1) * int_d_psiW_dy_d_psiV_dx +
				   Hooke (1,1,1,0) * int_d_psiW_dx_d_psiV_dy +
				   Hooke (1,1,1,1) * int_d_psiW_dy_d_psiV_dy ;                          }  }  }  
	} // just a block of code 
				 
	// impose Dirichlet boundary condition (clamping on left side)
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = DG .iterator ( tag::over_vertices );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;
		impose_value_of_unknown ( matrix_A, vector_b, numbering[V.core]-1, 0. );
		impose_value_of_unknown ( matrix_A, vector_b, numbering[V.core]-1+N, 0. );  }
	} { // just a block of code for hiding 'it'
	Mesh::Iterator it = HA .iterator ( tag::over_vertices );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;
		impose_value_of_unknown ( matrix_A, vector_b, numbering[V.core]-1, 0. );
		impose_value_of_unknown ( matrix_A, vector_b, numbering[V.core]-1+N, 0. );  }
	} // just a block of code 

	// apply vertical force along EF
	{ // just a block of code for hiding 'it'
	const double force_x =  0.;
	const double force_y =  1.;
	Mesh::Iterator it = EF .iterator ( tag::over_segments );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell seg = *it;
		Cell V = seg .base() .reverse(), W = seg .tip();
		double len_seg = y(W) - y(V);
		assert ( len_seg > 0. );
		vector_b [ numbering[V.core]-1   ] += len_seg / .2 * force_x;
		vector_b [ numbering[V.core]-1+N ] += len_seg / .2 * force_y;
		vector_b [ numbering[W.core]-1   ] += len_seg / .2 * force_x;
		vector_b [ numbering[W.core]-1+N ] += len_seg / .2 * force_y;  }
	} // just a block of code 

	// solve the system of linear equations
	Eigen::ConjugateGradient < Eigen::SparseMatrix<double>,
	                           Eigen::Lower | Eigen::Upper > cg;
	cg .compute ( matrix_A );
	Eigen::VectorXd vector_sol = cg .solve ( vector_b );
	// std::cout << "error " << cg .error() << std::endl;
	
	// elastic deformation of the mesh
	// run over all vertices of 'cantilever'
	if ( false )
	{ // just a block of code for hiding 'it'
	const double coef = 0.1;
	Mesh::Iterator it = cantilever .iterator ( tag::over_vertices );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;
		x(V) += coef * vector_sol [ numbering [V.core]-1   ];
		y(V) += coef * vector_sol [ numbering [V.core]-1+N ];  }
	} // just a block of code

	// compute the shape derivative = elastic energy
	Eigen::VectorXd optim_deform ( size_matrix );
	optim_deform .setZero();

	// we temporarily use the field optim_deform [ numbering[T.core]-1 ]
	// for keeping the elastic energy density
	// run over all segments of free_bdry
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = free_bdry .iterator ( tag::over_segments );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell seg = *it;
		Cell S = seg .base() .reverse(), T = seg .tip();
		Cell small_tri = cantilever .cell_behind ( seg, tag::surely_exists );
		// we compute the elastic energy in this small_tri
		double energy = 0.;
		fe_2d .dock_on ( small_tri );
		Mesh::Iterator it1 = small_tri .boundary().iterator ( tag::over_vertices );
		Mesh::Iterator it2 = small_tri .boundary().iterator ( tag::over_vertices );
		for ( it1 .reset(); it1 .in_range(); it1++ )
		{	Cell V = *it1 ;
			Function psiV = fe_2d .basis_function(V) ;
			for ( it2 .reset(); it2 .in_range(); it2++ )
			{	Cell W = *it2;  // V may be the same as W, no problem about that
				Function psiW = fe_2d .basis_function (W);
				std::vector < double > result = fe_2d .integrate 
					( tag::pre_computed, tag::replace, bf1, tag::by, psiW,
				                       tag::replace, bf2, tag::by, psiV );
				assert ( result .size() == 4 );
				double int_d_psiW_dx_d_psiV_dx = result[0];
				double int_d_psiW_dx_d_psiV_dy = result[1];
				double int_d_psiW_dy_d_psiV_dx = result[2];
				double int_d_psiW_dy_d_psiV_dy = result[3];
				size_t nV = numbering[V.core]-1, nW = numbering[W.core]-1;
				energy += vector_sol [ nV ] * vector_sol [ nW ] *
				  ( Hooke (0,0,0,0) * int_d_psiW_dx_d_psiV_dx +
				    Hooke (0,1,0,0) * int_d_psiW_dy_d_psiV_dx +
				    Hooke (0,0,0,1) * int_d_psiW_dx_d_psiV_dy +
				    Hooke (0,1,0,1) * int_d_psiW_dy_d_psiV_dy  );
				energy += vector_sol [ nV ] * vector_sol [ nW+N ] *
				  ( Hooke (1,0,0,0) * int_d_psiW_dx_d_psiV_dx +
				    Hooke (1,1,0,0) * int_d_psiW_dy_d_psiV_dx +
				    Hooke (1,0,0,1) * int_d_psiW_dx_d_psiV_dy +
				    Hooke (1,1,0,1) * int_d_psiW_dy_d_psiV_dy  );
				energy += vector_sol [ nV+N ] * vector_sol [ nW ] *
				  ( Hooke (0,0,1,0) * int_d_psiW_dx_d_psiV_dx +
				    Hooke (0,1,1,0) * int_d_psiW_dy_d_psiV_dx +
				    Hooke (0,0,1,1) * int_d_psiW_dx_d_psiV_dy +
				    Hooke (0,1,1,1) * int_d_psiW_dy_d_psiV_dy  );
				energy += vector_sol [ nV+N ] * vector_sol [ nW+N ] *
				  ( Hooke (1,0,1,0) * int_d_psiW_dx_d_psiV_dx +
				    Hooke (1,1,1,0) * int_d_psiW_dy_d_psiV_dx +
				    Hooke (1,0,1,1) * int_d_psiW_dx_d_psiV_dy +
				    Hooke (1,1,1,1) * int_d_psiW_dy_d_psiV_dy  );             }  }
		// we divide by the area of the triangle
		energy /= area ( small_tri );  // well, it's twice the area ...
		if ( not ( ( S == E ) or ( S == F ) ) )
		{	assert ( not S .belongs_to ( EF, tag::not_oriented ) );
			// assert ( not S .belongs_to ( DA, tag::not_oriented ) );
			optim_deform [ numbering[S.core]-1 ] += energy;       }
		if ( not ( ( T == E ) or ( T == F ) ) )
		{	assert ( not T .belongs_to ( EF, tag::not_oriented ) );
			// assert ( not T .belongs_to ( DA, tag::not_oriented ) );
			optim_deform [ numbering[T.core]-1 ] += energy;       }              }
	} // just a block of code

	// multiply by normal vector
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = free_bdry .iterator ( tag::over_vertices );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;
		Cell s1 = free_bdry .cell_in_front_of ( V, tag::may_not_exist );
		Cell s2 = free_bdry .cell_behind      ( V, tag::may_not_exist );
		if ( not ( s1 .exists() and s2 .exists() ) )  continue;
		std::vector < double > t1 = tangent ( s1 );
		double norm = std::sqrt ( t1[0] * t1[0] + t1[1] * t1[1] );
		for ( size_t i = 0; i < 2; i++ )  t1[i] /= norm;
		std::vector < double > t2 = tangent ( s2 );
		norm = std::sqrt ( t2[0] * t2[0] + t2[1] * t2[1] );
		for ( size_t i = 0; i < 2; i++ )  t2[i] /= norm;
		double energy = optim_deform [ numbering[V.core]-1 ];
		double coef = 1.;
		for ( size_t i = 0; i < 2; i++ )  coef += t1[i]*t2[i];
		coef = ( energy - lagrange ) / coef;
		optim_deform [ numbering[V.core]-1   ] =   ( t1[1] + t2[1] ) * coef;
		optim_deform [ numbering[V.core]-1+N ] = - ( t1[0] + t2[0] ) * coef;  }
	} // just a block of code

	// Dirichlet boundary CG & HA deserves a special treatment
	{ // just a block of code for hiding names
	// we multiply by two because we are at an extremity of the free boundary
	double energy_D = 2. * optim_deform [ numbering[D.core]-1 ],
	       energy_G = 2. * optim_deform [ numbering[G.core]-1 ];
	if ( ( energy_D > lagrange ) and y(D) >= 1. )
	{	energy_D = lagrange;  y(D) = 1.;  }
	if ( ( energy_G < lagrange ) and y(G) >= 0.85 )
	{	energy_G = lagrange;  y(G) = 0.85;  }
	if ( ( energy_G > lagrange ) and y(G) <= 0.673205 )
	{	energy_G = lagrange;  y(G) = 0.673205;  }
	// speed at D is  energy_D - lagrange
	// speed at G is  lagrange - energy_G
	size_t nn = DG .number_of ( tag::segments );
	double step = - ( energy_G + energy_D - 2.*lagrange ) / nn;
	Mesh::Iterator it = DG .iterator ( tag::over_vertices, tag::require_order );
	it .reset();
	for ( size_t i = 0; i <= nn; i++ )
	{	assert ( it .in_range() );
		Cell V = *it;
		it++;
		optim_deform [ numbering[V.core]-1   ] = 0.;
		optim_deform [ numbering[V.core]-1+N ] = energy_D - lagrange + i * step;   }
  assert ( not it .in_range() );
	} { // just a block of code for hiding names
	// we multiply by two because we are at an extremity of the free boundary
	double energy_A = 2. * optim_deform [ numbering[A.core]-1 ],
	       energy_H = 2. * optim_deform [ numbering[H.core]-1 ];
	if ( ( energy_A > lagrange ) and y(A) <= 0. )
	{	energy_A = lagrange;  y(A) = 0.;  }
	if ( ( energy_H < lagrange ) and y(H) <= 0.15 )
	{	energy_H = lagrange;  y(H) = 0.15;  }
	if ( ( energy_H > lagrange ) and y(H) >= 0.326795 )
	{	energy_H = lagrange;  y(H) = 0.326795;  }
	// speed at A is  lagrange - energy_A
	// speed at H is  energy_H - lagrange
	size_t nn = HA .number_of ( tag::segments );
	double step = - ( energy_A + energy_H - 2.*lagrange ) / nn;
	Mesh::Iterator it = HA .iterator ( tag::over_vertices, tag::require_order );
	it .reset();
	for ( size_t i = 0; i <= nn; i++ )
	{	assert ( it .in_range() );
		Cell V = *it;
		it++;
		optim_deform [ numbering[V.core]-1   ] = 0.;
		optim_deform [ numbering[V.core]-1+N ] = energy_H - lagrange + i * step;   }
  assert ( not it .in_range() );
	} // just a block of code

	// propagate the deformation on the entire mesh
	if ( false )
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = cantilever .iterator ( tag::over_vertices );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;
		if ( V .belongs_to ( bdry ) ) continue;
		size_t nb_neigh = 0;
		double def_x = 0., def_y = 0.;
		Mesh::Iterator it1 = cantilever .iterator ( tag::over_vertices, tag::around, V );
		for ( it1 .reset(); it1 .in_range(); it1++ )
		{	nb_neigh ++;
			Cell W = *it1;
			def_x += optim_deform [ numbering[W.core]-1   ];
			def_y += optim_deform [ numbering[W.core]-1+N ];  }
		optim_deform [ numbering[V.core]-1   ] = def_x / nb_neigh;
		optim_deform [ numbering[V.core]-1+N ] = def_y / nb_neigh;                         }
	} // just a block of code
			
	// move all vertices on the boundary
	{ // just a block of code for hiding 'it'
	Mesh::Iterator it = bdry .iterator ( tag::over_vertices );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;
		x(V) += eta * optim_deform [ numbering[V.core]-1   ];
		y(V) += eta * optim_deform [ numbering[V.core]-1+N ];  }
	} // just a block of code

	// eliminate short segments, split long segments
	{ // just a block of code for hiding variables
	const double lower_lim = seg_len * seg_len / 2.1;
	const double upper_lim = seg_len * seg_len * 2.1;
	for ( std::list < Mesh > ::iterator it =  free_bdry_pieces .begin();
	                                    it != free_bdry_pieces .end();   it++ )
	{	Mesh piece = *it;
		split_elim_segs ( piece, lower_lim, upper_lim );  }
	for ( std::list < Mesh > ::iterator it =  Dirichlet_bdry_pieces .begin();
	                                    it != Dirichlet_bdry_pieces .end();   it++ )
	{	Mesh piece = *it;
		split_elim_segs ( piece, lower_lim, upper_lim );  }
	} // just a block of code
	
	// regularize vertices on free boundary
	// similar to baricenters but takes into account the curvature
	// tries to keep the previous curvature
	{ // just a block of code for hiding 'i'
	for ( std::list < Mesh > ::iterator i =  free_bdry_pieces .begin();
	                                    i != free_bdry_pieces .end();   i++ )
	{	Mesh piece = *i;  smoothen ( piece );  }
	} // just a block of code

	// regularize Dirichlet boundary by merely applying baricenters
	{ // just a block of code for hiding 'it'
	for ( std::list < Mesh > ::iterator i =  Dirichlet_bdry_pieces .begin();
	                                    i != Dirichlet_bdry_pieces .end(); i++ )
	{	Mesh piece = *i;
		Mesh::Iterator it = piece .iterator ( tag::over_vertices );
		for ( it .reset(); it .in_range(); it++ )
		{	Cell V = *it;
			if ( V .is_inner_to ( piece ) ) piece .baricenter ( V );  }  }
	} // just a block of code

	bdry = Mesh ( tag::join, bdry_pieces );
	free_bdry = Mesh ( tag::join, free_bdry_pieces );

	cantilever = Mesh ( tag::frontal, tag::boundary, bdry, tag::desired_length, seg_len );

	cantilever .draw_ps ("cantilever.eps");

	}  // end of  for iter
		
	// write to file coordinates of boundary vertices
	{ // just a block of code for hiding 'f'
	std::ofstream f ("bdry.dat");

	f << seg_len << std::endl;

	// begin at lower left corner
	{ // just a block of code for hiding 'it'
	f << AB .number_of ( tag::segments ) << std::endl;
	Mesh::Iterator it = AB .iterator ( tag::over_vertices, tag::require_order );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;
		f << x(V) << " " << y(V) << std::endl;  }
	} // just a block of code

	// three pieces on the right hand side
	{ // just a block of code for hiding 'it'
	f << BE .number_of ( tag::segments ) << std::endl;
	Mesh::Iterator it = BE .iterator ( tag::over_vertices, tag::require_order );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;
		f << x(V) << " " << y(V) << std::endl;  }
	} { // just a block of code for hiding 'it'
	f << EF .number_of ( tag::segments ) << std::endl;
	Mesh::Iterator it = EF .iterator ( tag::over_vertices, tag::require_order );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;
		f << x(V) << " " << y(V) << std::endl;  }
	} { // just a block of code for hiding 'it'
	f << FC .number_of ( tag::segments ) << std::endl;
	Mesh::Iterator it = FC .iterator ( tag::over_vertices, tag::require_order );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;
		f << x(V) << " " << y(V) << std::endl;  }
	} // just a block of code
	
	// upper side
	{ // just a block of code for hiding 'it'
	f << CD .number_of ( tag::segments ) << std::endl;
	Mesh::Iterator it = CD .iterator ( tag::over_vertices, tag::require_order );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;
		f << x(V) << " " << y(V) << std::endl;  }
	} // just a block of code

	// three pieces on the left hand side
	{ // just a block of code for hiding 'it'
	f << DG .number_of ( tag::segments ) << std::endl;
	Mesh::Iterator it = DG .iterator ( tag::over_vertices, tag::require_order );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;
		f << x(V) << " " << y(V) << std::endl;  }
	} { // just a block of code for hiding 'it'
	f << semi_circle .number_of ( tag::segments ) << std::endl;
	Mesh::Iterator it = semi_circle .iterator ( tag::over_vertices, tag::require_order );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;
		f << x(V) << " " << y(V) << std::endl;  }
	} { // just a block of code for hiding 'it'
	f << HA .number_of ( tag::segments ) << std::endl;
	Mesh::Iterator it = HA .iterator ( tag::over_vertices, tag::require_order );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell V = *it;
		f << x(V) << " " << y(V) << std::endl;  }
	} // just a block of code
	
	// several inner curves
	f << inner .size() << std::endl;
	for ( std::list < Mesh > ::iterator il = inner .begin();
				il != inner .end(); il++ )
	{	Mesh circle = *il;
		f << circle .number_of ( tag::segments ) << std::endl;
		Mesh::Iterator it = circle .iterator ( tag::over_vertices, tag::require_order );
		for ( it .reset(); it .in_range(); it++ )
		{	Cell V = *it;
			f << x(V) << " " << y(V) << std::endl;  }                                      }

	} // just a block of code

	{ // just a block of code for hiding 'f'
	std::ofstream f ("bdry.eps");

	// begin at lower left corner
	{ // just a block of code for hiding 'it'
	f << "0.01 setlinewidth 0 0 1 setrgbcolor" << std::endl;
	Mesh::Iterator it = AB .iterator ( tag::over_segments, tag::require_order );
	it .reset();  assert ( it .in_range() );
	Cell seg = *it, V = seg .base() .reverse();
	f << x(V) << " " << y(V) << " moveto ";
	for ( ; it .in_range(); it++ )
	{	seg = *it, V = seg .tip();
		f << x(V) << " " << y(V) << " lineto"<< std::endl;  }
	f << "stroke" << std::endl;
	} // just a block of code

	// three pieces on the right hand side
	{ // just a block of code for hiding 'it'
	// keep previous color blue
	Mesh::Iterator it = BE .iterator ( tag::over_segments, tag::require_order );
	it .reset();  assert ( it .in_range() );
	Cell seg = *it, V = seg .base() .reverse();
	f << x(V) << " " << y(V) << " moveto ";
	for ( ; it .in_range(); it++ )
	{	seg = *it, V = seg .tip();
		f << x(V) << " " << y(V) << " lineto"<< std::endl;  }
	f << "stroke" << std::endl;
	} { // just a block of code for hiding 'it'
	f << "0.01 setlinewidth 1 0 0 setrgbcolor" << std::endl;
	Mesh::Iterator it = EF .iterator ( tag::over_segments, tag::require_order );
	it .reset();  assert ( it .in_range() );
	Cell seg = *it, V = seg .base() .reverse();
	f << x(V) << " " << y(V) << " moveto ";
	for ( ; it .in_range(); it++ )
	{	seg = *it, V = seg .tip();
		f << x(V) << " " << y(V) << " lineto"<< std::endl;  }
	f << "stroke" << std::endl;
	} { // just a block of code for hiding 'it'
	f << "0.01 setlinewidth 0 0 1 setrgbcolor" << std::endl;
	Mesh::Iterator it = FC .iterator ( tag::over_segments, tag::require_order );
	it .reset();  assert ( it .in_range() );
	Cell seg = *it, V = seg .base() .reverse();
	f << x(V) << " " << y(V) << " moveto ";
	for ( ; it .in_range(); it++ )
	{	seg = *it, V = seg .tip();
		f << x(V) << " " << y(V) << " lineto"<< std::endl;  }
	f << "stroke" << std::endl;
	} // just a block of code
	
	// upper side
	{ // just a block of code for hiding 'it'
	// keep previous color blue
	Mesh::Iterator it = CD .iterator ( tag::over_segments, tag::require_order );
	it .reset();  assert ( it .in_range() );
	Cell seg = *it, V = seg .base() .reverse();
	f << x(V) << " " << y(V) << " moveto ";
	for ( ; it .in_range(); it++ )
	{	seg = *it, V = seg .tip();
		f << x(V) << " " << y(V) << " lineto"<< std::endl;  }
	f << "stroke" << std::endl;
	} // just a block of code

	// three pieces on the left hand side
	{ // just a block of code for hiding 'it'
	f << "0.01 setlinewidth 0 0 1 setrgbcolor" << std::endl;
	Mesh::Iterator it = semi_circle .iterator ( tag::over_segments, tag::require_order );
	it .reset();  assert ( it .in_range() );
	Cell seg = *it, V = seg .base() .reverse();
	f << x(V) << " " << y(V) << " moveto ";
	for ( ; it .in_range(); it++ )
	{	seg = *it, V = seg .tip();
		f << x(V) << " " << y(V) << " lineto"<< std::endl;  }
	f << "stroke" << std::endl;
	f << "0.02 setlinewidth 0.6 0.2 0 setrgbcolor" << std::endl;
	f << "0. 1. moveto 0. 0.673205 lineto stroke" << std::endl;
	f << "0.02 setlinewidth 0.6 0.2 0 setrgbcolor" << std::endl;
	f << "0. 0.326795 moveto 0. 0. lineto stroke" << std::endl;
	} // just a block of code
	
	// three inner curves
	f << "0.01 setlinewidth 0 0 1 setrgbcolor" << std::endl;
	for ( std::list < Mesh > ::iterator il = inner .begin();
				il != inner .end(); il++ )
	{	Mesh circle = *il;
		Mesh::Iterator it = circle .iterator ( tag::over_segments, tag::require_order );
		it .reset();  assert ( it .in_range() );
		Cell seg = *it, V = seg .base() .reverse();
		f << x(V) << " " << y(V) << " moveto ";
		for ( ; it .in_range(); it++ )
		{	seg = *it, V = seg .tip();
			f << x(V) << " " << y(V) << " lineto"<< std::endl;  }
		f << "stroke" << std::endl;                                                      }
	
	} // just a block of code

	std::cout << "produced file cantilever.eps, bdry.dat and bdry.eps" << std::endl;

	return 0;
	
}  // end of  main

//------------------------------------------------------------------------------------------------------//


void split_elim_segs_closed ( Mesh & chain, const double lower_lim, const double upper_lim )

{	assert ( upper_lim > lower_lim );
	assert ( chain .dim() == 1 );

	Function xy = Manifold::working .coordinates();
	Function x = xy [0], y = xy [1];

	Mesh::Connected::OneDim * chain_core = tag::Util::assert_cast
		< Mesh::Core *, Mesh::Connected::OneDim * > ( chain .core );

	size_t n = chain_core->nb_of_segs;
	
	// first split long segments :
	Cell A = chain_core->first_ver .reverse();
	assert ( A .is_positive() );

	for ( size_t i = 0; i < n; i++ )
	{	Cell seg = chain .cell_in_front_of ( A, tag::surely_exists );
		assert ( seg .base() == A .reverse() );
		Cell B = seg .tip();
		if ( sq_length ( seg ) < upper_lim )
		{	A = B;  continue;  }
		Cell new_ver ( tag::vertex );
		Manifold::working .interpolate ( new_ver, 0.5, A, 0.5, B );
		seg .remove_from_mesh ( chain, tag::do_not_bother );
		Cell seg1 ( tag::segment, A .reverse(), new_ver );
		Cell seg2 ( tag::segment, new_ver .reverse(), B );
		seg1 .add_to_mesh ( chain, tag::do_not_bother );
		seg2 .add_to_mesh ( chain, tag::do_not_bother );
		A = B;  i++;  n++;  // skip both seg1 and seg2
	}  // end of for

	// now eliminate short segments :
	A = chain_core->first_ver .reverse();
	assert ( A .is_positive() );
	Cell AB = chain .cell_in_front_of ( A, tag::surely_exists );
	Cell B = AB .tip();
	assert ( B != chain_core->last_ver );
	Cell BC = chain .cell_in_front_of ( B, tag::surely_exists );
	Cell C = BC .tip();
	assert ( C != chain_core->last_ver );
	for ( size_t i = 0; i < n; i++ )
	{	Cell CD = chain .cell_in_front_of ( C, tag::surely_exists );
		// eliminate inner segment
		Cell D = CD .tip();
		if ( sq_length ( BC ) > lower_lim )
		{	A = B;  B = C;  C = D;  AB = BC;  BC = CD;  continue;  }
		// Manifold::working .interpolate ( B, 0.5, B, 0.5, C );
		Cell BD ( tag::segment, B .reverse(), D );
		// std::cout << x(B) << " " << y(B) << " moveto "
		// << x(D) << " " << y(D) << " lineto stroke" << std::endl;
		BC .remove_from_mesh ( chain, tag::do_not_bother );
		CD .remove_from_mesh ( chain, tag::do_not_bother );
		BD .add_to_mesh      ( chain, tag::do_not_bother );
		n--;
		C = D;  BC = BD;  continue;                                  }
	// assert ( A == chain_core->last_ver );

	chain_core->nb_of_segs = n;

}  // end of  split_elim_segs_closed
	

void split_elim_segs_open ( Mesh & chain, const double lower_lim, const double upper_lim )
	
{	assert ( upper_lim > lower_lim );
	assert ( chain .dim() == 1 );

	Function xy = Manifold::working .coordinates();
	Function x = xy [0], y = xy [1];

	Mesh::Connected::OneDim * chain_core = tag::Util::assert_cast
		< Mesh::Core *, Mesh::Connected::OneDim * > ( chain .core );

	// first split long segments :
	Cell A = chain_core->first_ver .reverse();
	assert ( A .is_positive() );

	while ( A != chain_core->last_ver )
	{	Cell seg = chain .cell_in_front_of ( A, tag::surely_exists );
		assert ( seg .base() == A .reverse() );
		Cell B = seg .tip();
		if ( sq_length ( seg ) < upper_lim )
		{	A = B;  continue;  }
		Cell new_ver ( tag::vertex );
		Manifold::working .interpolate ( new_ver, 0.5, A, 0.5, B );
		seg .remove_from_mesh ( chain, tag::do_not_bother );
		Cell seg1 ( tag::segment, A .reverse(), new_ver );
		Cell seg2 ( tag::segment, new_ver .reverse(), B );
		seg1 .add_to_mesh ( chain, tag::do_not_bother );
		seg2 .add_to_mesh ( chain, tag::do_not_bother );
		chain_core->nb_of_segs ++;
		A = B;  // skip both seg1 and seg2
	}  // end of while

	// now eliminate short segments :
	A = chain_core->first_ver .reverse();
	assert ( A .is_positive() );
	Cell AB = chain .cell_in_front_of ( A, tag::surely_exists );
	Cell B = AB .tip();
	assert ( B != chain_core->last_ver );
	Cell BC = chain .cell_in_front_of ( B, tag::surely_exists );
	Cell C = BC .tip();
	if ( C == chain_core->last_ver ) return;
	if ( sq_length ( AB ) < lower_lim )
	// eliminate first segment
	{	Cell AC ( tag::segment, A .reverse(), C );
		AB .remove_from_mesh ( chain, tag::do_not_bother );
		BC .remove_from_mesh ( chain, tag::do_not_bother );
		AC .add_to_mesh      ( chain, tag::do_not_bother );
		chain_core->nb_of_segs --;
		// first_ver remains the same
		AB = AC;  B = C;
		BC = chain .cell_in_front_of ( B, tag::surely_exists );
		C = BC .tip();                                          }
	while ( C != chain_core->last_ver )
	{	Cell CD = chain .cell_in_front_of ( C, tag::surely_exists );
		Cell D = CD .tip();
		if ( sq_length ( BC ) > lower_lim )
		{	A = B;  B = C;  C = D;  AB = BC;  BC = CD;  continue;  }
		// eliminate inner segment
		Manifold::working .interpolate ( B, 0.5, B, 0.5, C );
		Cell BD ( tag::segment, B .reverse(), D );
		BC .remove_from_mesh ( chain, tag::do_not_bother );
		CD .remove_from_mesh ( chain, tag::do_not_bother );
		BD .add_to_mesh      ( chain, tag::do_not_bother );
		chain_core->nb_of_segs --;
		C = D;  BC = BD;  continue;                                  }
	assert ( C == chain_core->last_ver );
	if ( sq_length ( BC ) < lower_lim )
	// eliminate last segment
	{	Cell AC ( tag::segment, A .reverse(), C );
		AB .remove_from_mesh ( chain, tag::do_not_bother );
		BC .remove_from_mesh ( chain, tag::do_not_bother );
		AC .add_to_mesh      ( chain, tag::do_not_bother );
		chain_core->nb_of_segs --;                          }
		// last_ver remains the same

}  // end of  split_elim_segs_open
	

void split_elim_segs ( Mesh & chain, const double lower_lim, const double upper_lim )

// eliminates, from a 1d mesh, segments shorter than sqrt ( lower_lim )
// splits segments longer than sqrt ( upper_lim )

{	Mesh::Connected::OneDim * chain_core = tag::Util::assert_cast
		< Mesh::Core *, Mesh::Connected::OneDim * > ( chain .core );

	if ( chain_core->first_ver .reverse() == chain_core->last_ver )  // closed loop
		split_elim_segs_closed ( chain, lower_lim, upper_lim );

	else  // open chain
		split_elim_segs_open ( chain, lower_lim, upper_lim );          }
	
//------------------------------------------------------------------------------------------------------//


// regularize vertices on a give curve
// similar to baricenters but takes into account the curvature
// tries to keep the previous curvature

inline double improve ( const double & lambda, const double & norm )
{	return lambda / 3. + 0.00001 * norm * norm / ( norm + lambda )
	  - 0.08 * lambda * lambda / norm - 0.5 * lambda * lambda * lambda / norm / norm;  }


void smoothen_open ( Mesh & chain )

{	Function xy = Manifold::working .coordinates();
	Function x = xy [0], y = xy [1];

	Mesh::Iterator it = chain .iterator ( tag::over_vertices, tag::require_order );
	it .reset();  assert ( it .in_range() );
	Cell B = *it;  // first vertex in chain
	Cell seg = chain .cell_in_front_of ( B, tag::surely_exists );
	Cell C = seg .tip();
	seg = chain .cell_in_front_of ( C, tag::surely_exists );
	Cell D = seg .tip();
	seg = chain .cell_in_front_of ( D, tag::may_not_exist );
	if ( not seg .exists() ) return;
	Cell E = seg .tip();
	double tanx = x(D) - x(B), tany = y(D) - y(B);
	double norm = std::sqrt ( tanx*tanx + tany*tany );
	tanx /= norm;  tany /= norm;
	double nx = - tany, ny = tanx;
	double old_lam = ( x(C) - x(B) ) * nx + ( y(C) - y(B) ) * ny;
	double lambda = improve ( ( x(D) - x(E) ) * nx + ( y(D) - y(E) ) * ny, norm );
	x(C) = x(B) + tanx * norm / 2. + nx * ( lambda + old_lam ) / 2.;
	y(C) = y(B) + tany * norm / 2. + ny * ( lambda + old_lam ) / 2.;
	Cell A ( tag::non_existent );
	while ( true )
	{	A = B;  B = C;  C = D;  D = E;
		seg = chain .cell_in_front_of ( E, tag::may_not_exist );
		if ( not seg .exists() ) break;
		E = seg .tip();
		tanx = x(D) - x(B), tany = y(D) - y(B);
		norm = std::sqrt ( tanx*tanx + tany*tany );
		tanx /= norm;  tany /= norm;
		nx = - tany, ny = tanx;
		old_lam = ( x(C) - x(B) ) * nx + ( y(C) - y(B) ) * ny;
		double lambdaA = ( x(B) - x(A) ) * nx + ( y(B) - y(A) ) * ny;
		double lambdaE = ( x(D) - x(E) ) * nx + ( y(D) - y(E) ) * ny;
		lambda = improve ( ( lambdaA + lambdaE ) / 2., norm );
		x(C) = x(B) + tanx * norm / 2. + nx * ( lambda + old_lam ) / 2.;
		y(C) = y(B) + tany * norm / 2. + ny * ( lambda + old_lam ) / 2.;     }
	assert ( A .exists() );   // E does not exist
	tanx = x(D) - x(B), tany = y(D) - y(B);
	norm = std::sqrt ( tanx*tanx + tany*tany );
	tanx /= norm;  tany /= norm;
	nx = - tany, ny = tanx;
	old_lam = ( x(C) - x(B) ) * nx + ( y(C) - y(B) ) * ny;
	lambda = improve ( ( x(B) - x(A) ) * nx + ( y(B) - y(A) ) * ny, norm );
	x(C) = x(B) + tanx * norm / 2. + nx * ( lambda + old_lam ) / 2.;
	y(C) = y(B) + tany * norm / 2. + ny * ( lambda + old_lam ) / 2.;               }

	
void smoothen_loop ( Mesh & chain )

{	Function xy = Manifold::working .coordinates();
	Function x = xy [0], y = xy [1];

	Mesh::Iterator it = chain .iterator ( tag::over_vertices, tag::require_order );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell C = *it;
		Cell seg = chain .cell_behind ( C, tag::surely_exists );
		Cell B = seg .base() .reverse();
		seg = chain .cell_behind ( B, tag::surely_exists );
		Cell A = seg .base() .reverse();
		seg = chain .cell_in_front_of ( C, tag::surely_exists );
		Cell D = seg .tip();
		seg = chain .cell_in_front_of ( D, tag::surely_exists );
		Cell E = seg .tip();
		double tanx = x(D) - x(B), tany = y(D) - y(B);
		double norm = std::sqrt ( tanx*tanx + tany*tany );
		tanx /= norm;  tany /= norm;
		double nx = - tany, ny = tanx;
		double old_lam = ( x(C) - x(B) ) * nx + ( y(C) - y(B) ) * ny;
		double lambdaA = ( x(B) - x(A) ) * nx + ( y(B) - y(A) ) * ny;
		double lambdaE = ( x(D) - x(E) ) * nx + ( y(D) - y(E) ) * ny;
		double lambda = improve ( ( lambdaA + lambdaE ) / 2., norm );
		x(C) = x(B) + tanx * norm / 2. + nx * ( lambda + old_lam ) / 2.;
		y(C) = y(B) + tany * norm / 2. + ny * ( lambda + old_lam ) / 2.;  }                }


void smoothen ( Mesh & chain )

{	Mesh::Connected::OneDim * chain_core = tag::Util::assert_cast
		< Mesh::Core *, Mesh::Connected::OneDim * > ( chain .core );

	if ( chain_core->first_ver .reverse() == chain_core->last_ver )  // closed loop
		smoothen_loop ( chain );

	else  // open chain
		smoothen_open ( chain );                                       }
