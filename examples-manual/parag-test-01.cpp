
// solve a celullar problem, elasticity, given macroscopic strain
// instead of zero average we impose sum zero
// square periodicity, triangular elements, circular hole
// identification of opposite sides based order (when exporting msh)


#include "maniFEM.h"
#include <fstream>
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
		// elements.reserve(total_dim);
		elements.insert(elements.end(),total_dim,0.0);
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


int main ( )

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0], y = xy[1];

	size_t n = 10;
	double d = 2.6 / double(n);

	Cell A ( tag::vertex );  x(A) = -1.3;  y(A) = -1.3;
	Cell B ( tag::vertex );  x(B) =  1.3;  y(B) = -1.3;
	Cell C ( tag::vertex );  x(C) =  1.3;  y(C) =  1.3;
	Cell D ( tag::vertex );  x(D) = -1.3;  y(D) =  1.3;

	Mesh AB ( tag::segment, A.reverse(), B, tag::divided_in, n );
	Mesh BC ( tag::segment, B.reverse(), C, tag::divided_in, n );
	Mesh CD ( tag::segment, C.reverse(), D, tag::divided_in, n );
	Mesh DA ( tag::segment, D.reverse(), A, tag::divided_in, n );

	Manifold circle = RR2.implicit ( x*x + y*y == 0.7 );
	Mesh inner ( tag::progressive, tag::entire_manifold, circle, tag::desired_length, d );

	Mesh bdry ( tag::join, AB, BC, CD, DA, inner.reverse() );

	RR2.set_as_working_manifold();
	Mesh square ( tag::progressive, tag::boundary, bdry, tag::desired_length, d );

	Mesh torus = square.fold ( tag::identify, AB, tag::with, CD.reverse(),
	                           tag::identify, BC, tag::with, DA.reverse(),
	                           tag::use_existing_vertices                 );
							   
	// Hooke's Law 
	double lambda = 1., mu = 3.;

	myTensor <double> Hooke(2,2,2,2);
	Hooke(0,0,0,0) = 2*mu+lambda;
	Hooke(0,0,0,1) = lambda;
	Hooke(0,0,1,0) = lambda;
	Hooke(0,0,1,1) = lambda;
	Hooke(0,1,0,0) = 0.;
	Hooke(0,1,0,1) = mu;
	Hooke(0,1,1,0) = mu;
	Hooke(0,1,1,1) = 0.;
	Hooke(1,0,0,0) = 0.;
	Hooke(1,0,0,1) = mu;
	Hooke(1,0,1,0) = mu;
	Hooke(1,0,1,1) = 0.;
	Hooke(1,1,0,0) = lambda;
	Hooke(1,1,0,1) = lambda;
	Hooke(1,1,1,0) = lambda;
	Hooke(1,1,1,1) = 2*mu+lambda;

	// declare the type of finite element
	FiniteElement fe ( tag::with_master, tag::triangle, tag::Lagrange, tag::of_degree, 1 );
	fe.set_integrator ( tag::Gauss, tag::tri_4 );

	// we number all nodes in 'square', not only those belonging to 'torus'
	std::map < Cell, size_t > numbering;
	{ // just a block of code for hiding 'it' and 'counter'
	CellIterator it = square.iterator ( tag::over_vertices );
	size_t counter = 0;
	for ( it.reset() ; it.in_range(); it++ )
	{	Cell V = *it;  numbering[V] = counter;  ++counter;  }
	assert ( counter == numbering.size() );
	} // just a block of code

	size_t number_dofs = numbering.size();
	std::cout << "global matrix " << 2*number_dofs + 3 << "x" << 2*number_dofs << std::endl;
	Eigen::SparseMatrix < double > matrix_A ( 2*number_dofs + 3, 2*number_dofs );
	
	matrix_A.reserve ( Eigen::VectorXi::Constant ( 2*number_dofs, 17 ) );
	// since we will be working with a mesh of triangles,
	// there will be, in average, 17=2*(6+1)+3 non-zero elements per column
	// the diagonal entry plus six neighbour vertices plus the last three equations

	// we fill the main diagonal with ones
	// then we put zero for vertices belonging to 'torus'
	{ // just a block of code for hiding 'it'
	for ( size_t i = 0; i < 2*number_dofs; i++ ) matrix_A.coeffRef ( i, i ) = 1.;
	CellIterator it = torus.iterator ( tag::over_vertices );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell V = *it;
		matrix_A.coeffRef ( 2*numbering[V], 2*numbering[V] ) = 0.;
		matrix_A.coeffRef ( 2*numbering[V]+1, 2*numbering[V]+1 ) = 0.;		}
	} // just a block of code for hiding 'it'
	
	Eigen::VectorXd vector_b ( 2*number_dofs + 3 ), vector_sol ( 2*number_dofs );
	vector_b.setZero();

	xy = Manifold::working.coordinates();
	x = xy[0];  y = xy[1];

	// macroscopic temperature gradient
	myTensor < double > macro_strain(2,2);
	macro_strain(0,0)=1.;
	macro_strain(0,1)=0.;
	macro_strain(1,0)=0.;
	macro_strain(1,1)=1.;
	Function::Jump jump_of_u_1 = macro_strain(0,0) * x.jump() + macro_strain(0,1) * y.jump();
	Function::Jump jump_of_u_2 = macro_strain(1,0) * x.jump() + macro_strain(1,1) * y.jump();
	// run over all square cells composing 'torus'
	{ // just a block of code for hiding 'it'
	CellIterator it = torus.iterator ( tag::over_cells_of_max_dim );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell small_tri = *it;
		fe.dock_on ( small_tri, tag::winding );
		// run twice over the four vertices of 'small_tri'
		CellIterator it_V = small_tri.boundary().iterator ( tag::over_vertices );
		for ( it_V.reset(); it_V.in_range(); it_V++ )
		{	Cell V = *it_V;
			// perhaps implement an interator returning a vertex and a segment
			Cell seg = small_tri .boundary(). cell_in_front_of ( V );
			Cell W = V;
			Function psi_V = fe .basis_function ( V ),
			         d_psiV_dx = psi_V .deriv ( x ),
			         d_psiV_dy = psi_V .deriv ( y );
			double jump_V_W_1 = 0.;
			double jump_V_W_2 = 0.;
			while ( true )
			{	assert ( W == seg.base().reverse() );
				// V may be the same as W, no problem about that
				Function psi_W = fe .basis_function ( W ),
				         d_psiW_dx = psi_W .deriv ( x ),
				         d_psiW_dy = psi_W .deriv ( y );
				// 'fe' is already docked on 'small_tri' so this will be the domain of integration
				double int_d_psiW_dx_d_psiV_dx = fe.integrate ( d_psiW_dx * d_psiV_dx );
				double int_d_psiW_dx_d_psiV_dy = fe.integrate ( d_psiW_dx * d_psiV_dy );
				double int_d_psiW_dy_d_psiV_dx = fe.integrate ( d_psiW_dy * d_psiV_dx );
				double int_d_psiW_dy_d_psiV_dy = fe.integrate ( d_psiW_dy * d_psiV_dy );
				myTensor<double> energy(2,2);
				energy(0,0) = Hooke(0,0,0,0) * int_d_psiW_dx_d_psiV_dx +
				              Hooke(0,0,0,1) * int_d_psiW_dy_d_psiV_dx +
				              Hooke(0,1,0,0) * int_d_psiW_dx_d_psiV_dy +
				              Hooke(0,1,0,1) * int_d_psiW_dy_d_psiV_dy ;
				energy(0,1) = Hooke(0,0,1,0) * int_d_psiW_dx_d_psiV_dx +
				              Hooke(0,0,1,1) * int_d_psiW_dy_d_psiV_dx +
				              Hooke(0,1,1,0) * int_d_psiW_dx_d_psiV_dy +
				              Hooke(0,1,1,1) * int_d_psiW_dy_d_psiV_dy ;
				energy(1,0) = Hooke(1,0,0,0) * int_d_psiW_dx_d_psiV_dx +
				              Hooke(1,0,0,1) * int_d_psiW_dy_d_psiV_dx +
				              Hooke(1,1,0,0) * int_d_psiW_dx_d_psiV_dy +
				              Hooke(1,1,0,1) * int_d_psiW_dy_d_psiV_dy ;
				energy(1,1) = Hooke(1,0,1,0) * int_d_psiW_dx_d_psiV_dx +
				              Hooke(1,0,1,1) * int_d_psiW_dy_d_psiV_dx +
				              Hooke(1,1,1,0) * int_d_psiW_dx_d_psiV_dy +
				              Hooke(1,1,1,1) * int_d_psiW_dy_d_psiV_dy ;
				matrix_A.coeffRef ( 2*numbering[V], 2*numbering[W] ) += energy(0,0);
				matrix_A.coeffRef ( 2*numbering[V], 2*numbering[W]+1 ) += energy(0,1);
				matrix_A.coeffRef ( 2*numbering[V]+1, 2*numbering[W] ) += energy(1,0);
				matrix_A.coeffRef ( 2*numbering[V]+1, 2*numbering[W]+1 ) += energy(1,1);
				vector_b ( 2*numbering[V] ) -= jump_V_W_1 * energy(0,0) + jump_V_W_2 * energy(0,1);
				vector_b ( 2*numbering[V]+1 ) -= jump_V_W_1 * energy(1,0) + jump_V_W_2 * energy(1,1);
				jump_V_W_1 += jump_of_u_1 ( seg.winding() );// a*g actualiza
				jump_V_W_2 += jump_of_u_2 ( seg.winding() );
				W = seg.tip();
				if ( V == W ) break;
				seg = small_tri .boundary() .cell_in_front_of ( seg.tip() );                          }
			// here  jump_V_W_1 and jump_V_W_2  should be zero again
			// but we do not assert that, rounding errors may mess up things
		}  }
	} // just a block of code for hiding 'it'

	// we add, as the last three equations, the condition of zero translation and zero rotation
	// actually, a more rudimentary condition : two lines of  zeros and ones alternated 
	for ( size_t i = 0; i < number_dofs; i++ )
	{	matrix_A.coeffRef ( 2*number_dofs, 2*i ) = 1.;
		matrix_A.coeffRef ( 2*number_dofs+1, 2*i+1 ) = 1.;  }
	// on the last line we put y alternated with -x (rotation)
	{ // just a block of code for hiding 'it'
	CellIterator it = torus.iterator ( tag::over_vertices );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell V = *it;
		matrix_A.coeffRef ( 2*number_dofs+2, 2*numbering[V] ) = y(V, tag::winding,0);
		matrix_A.coeffRef ( 2*number_dofs+2, 2*numbering[V]+1 ) = -x(V, tag::winding,0);  }
	} // just a block of code for hiding 'it'

	matrix_A .makeCompressed();

	Eigen::SparseQR < Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;

	solver.compute ( matrix_A );
	if ( solver.info() != Eigen::Success )
	{	std::cout << "Eigen solver.compute failed" << std::endl;
		exit ( 0 );                                              }

	vector_sol = solver.solve ( vector_b );
	if ( solver.info() != Eigen::Success )
	{	std::cout << "Eigen solver.solve failed" << std::endl;
		exit ( 0 );                                             }
	
		
	myTensor<double> macro_stress(2,2);
	for(size_t i=0; i<2; i++ )
	for(size_t j=0; j<2; j++)
		macro_stress(i,j) = 0.;
	{ // just a block of code for hiding 'it'
	CellIterator it = torus.iterator ( tag::over_cells_of_max_dim );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell small_tri = *it;
		fe.dock_on ( small_tri, tag::winding );
		// run twice over the four vertices of 'small_tri'
		CellIterator it_V = small_tri.boundary().iterator ( tag::over_vertices );
		double jump_V_1 = 0., jump_V_2 = 0.; 
		for ( it_V .reset(); it_V .in_range(); it_V++ )
		{	Cell V = *it_V;
			Function psi_V = fe .basis_function ( V ),
			         d_psi_V_dx = psi_V .deriv ( x ),
			         d_psi_V_dy = psi_V .deriv ( y );
			for ( size_t i=0; i<2; i++ )
			for ( size_t j=0; j<2; j++ )
				macro_stress(i,j) +=
					( vector_sol(2*numbering[V]) + jump_V_1 ) *
					  ( Hooke(i,j,0,0) * fe.integrate(d_psi_V_dx) +
					    Hooke(i,j,0,1) * fe.integrate(d_psi_V_dy)  ) +
					( vector_sol(2*numbering[V]+1) + jump_V_2 ) *
					  ( Hooke(i,j,1,0) * fe.integrate(d_psi_V_dx) +
					    Hooke(i,j,1,1) * fe.integrate(d_psi_V_dy)  )  ;
		Cell seg=small_tri .boundary() .cell_in_front_of(V);
		jump_V_1 += jump_of_u_1 ( seg .winding() );
		jump_V_2 += jump_of_u_2 ( seg .winding() );                 }             }
	} // just a block of code for hiding 'it'
	
	cout << "macro strain " << macro_strain(0,0) << " "<< macro_strain(0,1) << " "
	                        << macro_strain(1,0) << " "<< macro_strain(1,1) << " " << endl;
	cout << "macro stress " << macro_stress(0,0) << " "<< macro_stress(0,1) << " "
	                        << macro_stress(1,0) << " "<< macro_stress(1,1) << " " << endl;

	RR2 .set_as_working_manifold();
	xy = Manifold::working.coordinates();
	x = xy[0];  y = xy[1];
	square.export_msh ("cell-elast-strain.msh", numbering );
	
	{ // just a block of code for hiding variables
	std::ofstream solution_file ("cell-elast-strain.msh", std::fstream::app );
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
	CellIterator it = torus.iterator ( tag::over_vertices );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell P = *it;
		size_t j = numbering [ P ];
		solution_file << j + 1 << " " << vector_sol[2*j] << " "
	                            << vector_sol[2*j+1] << " 0. "<< std::endl;  }
	size_t j = numbering [ B ];
	assert ( not A .belongs_to ( torus ) );
	solution_file << numbering[A] + 1 << " "
		<< vector_sol[2*j] + macro_strain(0,0) * ( x ( A ) - x ( B ) )
		                 + macro_strain(0,1) * ( y ( A ) - y ( B ) ) <<" " 
		<< vector_sol[2*j+1] + macro_strain(1,0) * ( x ( A ) - x ( B ) )
		                 + macro_strain(1,1) * ( y ( A ) - y ( B ) ) << " 0." << std::endl;
	assert ( not C .belongs_to ( torus ) );
	solution_file << numbering[C] + 1 << " "
		<< vector_sol[2*j] + macro_strain(0,0) * ( x ( C ) - x ( B ) )
		                 + macro_strain(0,1) * ( y ( C ) - y ( B ) ) <<" " 
		<< vector_sol[2*j+1] + macro_strain(1,0) * ( x ( C ) - x ( B ) )
		                 + macro_strain(1,1) * ( y ( C ) - y ( B ) ) << " 0." << std::endl;
	assert ( not D .belongs_to ( torus ) );
	solution_file << numbering[D] + 1 << " "
		<< vector_sol[2*j] + macro_strain(0,0) * ( x ( D ) - x ( B ) )
		                 + macro_strain(0,1) * ( y ( D ) - y ( B ) ) <<" " 
		<< vector_sol[2*j+1] + macro_strain(1,0) * ( x ( D ) - x ( B ) )
		                 + macro_strain(1,1) * ( y ( D ) - y ( B ) ) << " 0." << std::endl;
	CellIterator it_AB = AB .iterator ( tag::over_vertices, tag::require_order );
	CellIterator it_CD = CD .iterator ( tag::over_vertices, tag::backwards );
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
		solution_file << numbering[W] + 1 << " "
	    << vector_sol[2*j] + macro_strain(0,0) * ( x ( W ) - x ( V ) )
		                 + macro_strain(0,1) * ( y ( W ) - y ( V ) ) <<" " 
		  << vector_sol[2*j+1] + macro_strain(1,0) * ( x ( W ) - x ( V ) )
		                 + macro_strain(1,1) * ( y ( W ) - y ( V ) ) << " 0." << std::endl;  }
	CellIterator it_BC = BC .iterator ( tag::over_vertices, tag::require_order );
	CellIterator it_DA = DA .iterator ( tag::over_vertices, tag::backwards );
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
		solution_file << numbering[W] + 1 << " "
		  << vector_sol[2*j] + macro_strain(0,0) * ( x ( W ) - x ( V ) )
		                 + macro_strain(0,1) * ( y ( W ) - y ( V ) ) <<" " 
		  << vector_sol[2*j+1] + macro_strain(1,0) * ( x ( W ) - x ( V ) )
		                 + macro_strain(1,1) * ( y ( W ) - y ( V ) ) << " 0." << std::endl;  }
	} // just a block of code

	std::cout << "produced file cell-elast-strain.msh" << std::endl;

	return 0;

}
