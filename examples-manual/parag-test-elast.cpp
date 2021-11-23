
// example presented in section 7.2 of the manual
// rudimentary use of finite elements

#include "maniFEM.h"
#include "math.h"

#include <fstream>
#include <Eigen/Sparse>

using namespace maniFEM;
using namespace std;

#include <list>
#include <vector>
#include <iostream>
#include <cstring>
#include "assert.h"

template <typename T>
class myTensor
{ // data:
  public:
  vector<T> elements;
  list<unsigned int> dimensions, cumulative_dims;
  unsigned int total_dim;
  // constructors:
  myTensor () {};
  myTensor (const list<unsigned int> dims)
  { dimensions = dims;
  allocate_space (); };
  myTensor (const list<char> dims)
  { dimensions = str2list(dims);
  allocate_space (); };
  myTensor (const char dims[])
  { dimensions = str2list(dims);
  allocate_space (); };
  myTensor (unsigned int i, unsigned int j,
            unsigned int k, unsigned int l)
            { dimensions.push_back(i);
    dimensions.push_back(j);
    dimensions.push_back(k);
    dimensions.push_back(l);
    allocate_space (); };
  myTensor (unsigned int i, unsigned int j, unsigned int k)
  { dimensions.push_back(i);
    dimensions.push_back(j);
    dimensions.push_back(k);
    allocate_space (); };
  myTensor (unsigned int i, unsigned int j)
  { dimensions.push_back(i);
    dimensions.push_back(j);
    allocate_space (); };
  myTensor (int i)
  { dimensions.push_back(i);
  allocate_space (); };
  ~myTensor () { };
  // methods:
  list<unsigned int> str2list (const list<char> lc)
  { const unsigned int izero = int('0');
    list<char>::iterator i;
    for (i=lc.begin();i!=lc.end();i++)
    { assert (*i >= '0');
    assert (*i <= '9'); }
    list<unsigned int> li;
    for (i=lc.begin();i!=lc.end();i++)
      li.push_back(int(*i)-izero);
      return li; }
  list<unsigned int> str2list (const char lc[])
  { const unsigned int izero = int('0');
    for (int i=0;i<strlen(lc);i++)
    { assert (lc[i] >= '0');
    assert (lc[i] <= '9'); }
    list<unsigned int> li;
    for (int i=0;i<strlen(lc);i++)
      li.push_back(int(lc[i])-izero);
      return li; }
  int allocate_space ()
  { total_dim = 1;
    list<unsigned int>::iterator k;
    for (k=dimensions.begin();k!=dimensions.end();k++)
    { cumulative_dims.push_back(total_dim);
    total_dim *= *k; }
    //elements.reserve(total_dim);
    elements.insert(elements.end(),total_dim,0.0);
    assert (elements.size() == total_dim); };
  T& operator()(list<unsigned int> index)
  { assert (index.size() == dimensions.size());
    unsigned int pointer = 0;
    list<unsigned int>::iterator i,d,cd;
    for (i=index.begin(),d=dimensions.begin(),
         cd=cumulative_dims.begin();
         i!=index.end();i++,d++,cd++)
         { assert (*i >= 0);
      assert (*i < *d);
      pointer += (*i)*(*cd); }
      return (elements[pointer]); }
  T& operator()(const char index[])
  { return operator()(str2list(index)); }
  T& operator()(unsigned int i, unsigned int j,
                unsigned int k, unsigned int l)
                { assert (dimensions.size() == 4);
    list<unsigned int> index;
    index.push_back(i);
    index.push_back(j);
    index.push_back(k);
    index.push_back(l);
    return operator()(index); }
  T& operator()(unsigned int i, unsigned int j, unsigned int k)
  { assert (dimensions.size() == 3);
    list<unsigned int> index;
    index.push_back(i);
    index.push_back(j);
    index.push_back(k);
    return operator()(index); }
  T& operator()(unsigned int i, unsigned int j)
  { assert (dimensions.size() == 2);
    list<unsigned int> index;
    index.push_back(i);
    index.push_back(j);
    return operator()(index); }
  T& operator()(int i)
  { assert (dimensions.size() == 1);
    list<unsigned int> index(1,i);
    return operator()(index); }
    };

void impose_value_of_unknown
(	Eigen::SparseMatrix <double> & matrix_A, Eigen::VectorXd & vector_b,
	size_t i, double val                                                 )

// in a system of linear equations, destroy equation 'i' and impose u(i) = val
// change also column 'i' of the matrix, just to preserve symmetry of the matrix

// used for imposing Dirichlet boundary conditions

{	size_t size_matrix = matrix_A.innerSize();
	vector_b(i) = val;
	for ( size_t j = 0; j < size_matrix; j++ )
		matrix_A.coeffRef ( i, j ) = 0.;
	matrix_A.coeffRef ( i, i ) = 1.;
	for ( size_t j = 0; j < size_matrix; j++ )
	{	if ( i == j ) continue;
		vector_b(j) -= matrix_A.coeffRef ( j, i ) * val;
		matrix_A.coeffRef ( j, i ) = 0.;                          }  }

	
int main ()

{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0],  y = xy[1];

	// build a 20x10 quadrangular mesh
	Cell A ( tag::vertex );  x(A) = 0.;   y(A) = 0.;
	Cell B ( tag::vertex );  x(B) = 1.;   y(B) = 0.;
//	Cell F ( tag::vertex );  x(F) = 2.;  y(F) = 0.5;
	Cell C ( tag::vertex );  x(C) = 1.;   y(C) = 1.;
	Cell D ( tag::vertex );  x(D) = 0.;   y(D) = 1.;
	Mesh AB ( tag::segment, A.reverse(), B, tag::divided_in, 10 );
	Mesh BC ( tag::segment, B.reverse(), C, tag::divided_in, 10 );
	Mesh CD ( tag::segment, C.reverse(), D, tag::divided_in, 10 );
	Mesh DA ( tag::segment, D.reverse(), A, tag::divided_in, 10 );
	Mesh ABCD ( tag::rectangle, AB, BC, CD, DA );
	Mesh Front_ABCD ( tag::join, AB, BC, CD, DA );

	// declare the type of finite element
	FiniteElement fe ( tag::with_master, tag::quadrangle, tag::Lagrange, tag::of_degree, 1 );
	Integrator integ = fe.set_integrator ( tag::Gauss, tag::quad_4 );
	
	// Hook's Law 
	double lambda = 1., mu = 3.;

    myTensor <double> Hook(2,2,2,2);
		Hook(0,0,0,0) = 2*mu+lambda;
		Hook(0,0,0,1) = lambda;
		Hook(0,0,1,0) = lambda;
		Hook(0,0,1,1) = lambda;
		Hook(0,1,0,0) = 0.;
		Hook(0,1,0,1) = mu;
		Hook(0,1,1,0) = mu;
		Hook(0,1,1,1) = 0.;
		Hook(1,0,0,0) = 0.;
		Hook(1,0,0,1) = mu;
		Hook(1,0,1,0) = mu;
		Hook(1,0,1,1) = 0.;
		Hook(1,1,0,0) = lambda;
		Hook(1,1,0,1) = lambda;
		Hook(1,1,1,0) = lambda;
		Hook(1,1,1,1) = 2*mu+lambda;
	
	// there will be a more elegant and efficient way of producing the numbering
	std::map < Cell::Core *, size_t > numbering;
	{ // just a block of code for hiding 'it' and 'counter'
	CellIterator it = ABCD.iterator ( tag::over_vertices );
	size_t counter = 0;
	for ( it.reset() ; it.in_range(); it++ )
	{	Cell p = *it;  ++counter;  numbering [p.core] = counter;  }
	} // just a block of code 
	// in each node P with number i>=0 ,  u_x(P)=vector_sol(i) , u_y(P)=vector_sol(i+N) where N is the number of nodes
    //  i=numbering(P)-1	
	size_t N = ABCD.number_of ( tag::vertices );
	assert ( N == numbering.size() );
	size_t size_matrix = 2*N;
	std::cout << "global matrix " << size_matrix << "x" << size_matrix << std::endl;
	Eigen::SparseMatrix <double> matrix_A ( size_matrix, size_matrix );
	Eigen::VectorXd vector_b ( size_matrix ), vector_sol ( size_matrix );
	vector_b.setZero();

	// run over all square cells composing ABCD
	{ // just a block of code for hiding 'it'
	double vol_force_x = 5. ;
	double vol_force_y = 10. ;
	CellIterator it = ABCD.iterator ( tag::over_cells_of_dim, 2 );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell small_square = *it;
		fe.dock_on ( small_square );
		// run twice over the four vertices of 'small_square'
		CellIterator it1 = small_square.boundary().iterator ( tag::over_vertices );
		CellIterator it2 = small_square.boundary().iterator ( tag::over_vertices );
		for ( it1.reset(); it1.in_range(); it1++ ) {
			Cell V = *it1 ;
			Function psiV = fe.basis_function(V) ;
//			Function volumic_force = sin(xy[0]) * sin(xy[1]) ;
		for ( it2.reset(); it2.in_range(); it2++ )
		{	Cell W = *it2;  // V may be the same as W, no problem about that
			// std::cout << "vertices V=(" << x(V) << "," << y(V) << ") " << numbering[V.core] << ", W=("
			// 					<< x(W) << "," << y(W) << ") " << numbering[W.core] << std::endl;
			Function psiW = fe.basis_function(W),
			         d_psiV_dx = psiV.deriv(x),
		           d_psiV_dy = psiV.deriv(y),
		           d_psiW_dx = psiW.deriv(x),
		           d_psiW_dy = psiW.deriv(y);
				   double int_d_psiW_dx_d_psiV_dx = fe.integrate ( d_psiW_dx * d_psiV_dx );
				   double int_d_psiW_dx_d_psiV_dy = fe.integrate ( d_psiW_dx * d_psiV_dy );
				   double int_d_psiW_dy_d_psiV_dx = fe.integrate ( d_psiW_dy * d_psiV_dx );
				   double int_d_psiW_dy_d_psiV_dy = fe.integrate ( d_psiW_dy * d_psiV_dy );
			// 'fe' is already docked on 'small_square' so this will be the domain of integration
			matrix_A.coeffRef ( numbering[V.core]-1, numbering[W.core]-1 ) += 
			Hook(0,0,0,0) *	int_d_psiW_dx_d_psiV_dx +
			Hook(0,0,0,1) *	int_d_psiW_dy_d_psiV_dx +
			Hook(0,1,0,0) *	int_d_psiW_dx_d_psiV_dy +
			Hook(0,1,0,1) *	int_d_psiW_dy_d_psiV_dy ;
			matrix_A.coeffRef ( numbering[V.core]-1, numbering[W.core]-1+N ) += 
			Hook(0,0,1,0) *	int_d_psiW_dx_d_psiV_dx +
			Hook(0,0,1,1) *	int_d_psiW_dy_d_psiV_dx +
			Hook(0,1,1,0) *	int_d_psiW_dx_d_psiV_dy +
			Hook(0,1,1,1) *	int_d_psiW_dy_d_psiV_dy ;
			matrix_A.coeffRef ( numbering[V.core]-1+N, numbering[W.core]-1 ) += 
			Hook(1,0,0,0) *	int_d_psiW_dx_d_psiV_dx +
			Hook(1,0,0,1) *	int_d_psiW_dy_d_psiV_dx +
			Hook(1,1,0,0) *	int_d_psiW_dx_d_psiV_dy +
			Hook(1,1,0,1) *	int_d_psiW_dy_d_psiV_dy ;
			matrix_A.coeffRef ( numbering[V.core]-1+N, numbering[W.core]-1+N ) += 
			Hook(1,0,1,0) *	int_d_psiW_dx_d_psiV_dx +
			Hook(1,0,1,1) *	int_d_psiW_dy_d_psiV_dx +
			Hook(1,1,1,0) *	int_d_psiW_dx_d_psiV_dy +
			Hook(1,1,1,1) *	int_d_psiW_dy_d_psiV_dy ;
			// std::cout << fe.integrate ( psiV * psiW ) << " "
			//          << fe.integrate ( d_psiV_dx * d_psiW_dx ) << " "
			//					<< fe.integrate ( d_psiV_dx * d_psiW_dy ) << " "
			//					<< fe.integrate ( d_psiV_dy * d_psiW_dx ) << " "
			//					<< fe.integrate ( d_psiV_dy * d_psiW_dy ) << " " << std::endl;
		}  
		double int_psiV = fe.integrate ( psiV );
		vector_b.coeffRef ( numbering[V.core]-1 ) += vol_force_x * int_psiV ;
		vector_b.coeffRef ( numbering[V.core]-1+N ) += vol_force_y * int_psiV ;
		}
	}
	} // just a block of code 

	// impose Dirichlet boundary conditions  u = xy
	{ // just a block of code for hiding 'it'
//	CellIterator it = Front_ABCD.iterator ( tag::over_vertices );
    CellIterator it = AB.iterator ( tag::over_vertices );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell P = *it;
		size_t i = numbering[P.core]-1;
		impose_value_of_unknown ( matrix_A, vector_b, i, 0. ); 
        impose_value_of_unknown ( matrix_A, vector_b, i+N, 0. ); 		}
	}  // just a block of code for hiding 'it' 
	
	
	// concentrated load in F
///	size_t i = numbering[F.core]-1;
//	vector_b(i+N) += -1.;
	
//	for ( size_t i = 0; i < size_matrix; i++ )
//	{	std::cout << i << " ";
//		for ( size_t j = 0; j < size_matrix; j++ )
//			std::cout << matrix_A.coeffRef (i,j) << " ";
//		std::cout << std::endl;                           }

//	std::cout << std::endl;
//	for ( size_t j = 0; j < size_matrix; j++ )
//		std::cout << vector_b (j) << " ";
//	std::cout << std::endl;
//	for ( size_t i = 0; i < size_matrix; i++ )
//		std::cout << vector_sol (i) << " ";
//	std::cout << std::endl;

	// solve the system of linear equations
	Eigen::ConjugateGradient < Eigen::SparseMatrix<double>,
	                           Eigen::Lower|Eigen::Upper    > cg;
	cg.compute ( matrix_A );
	vector_sol = cg.solve ( vector_b );

	ABCD.export_msh ("square_elast.msh", numbering );

	{ // just a block of code for hiding variables
	ofstream solution_file ("square_elast.msh", fstream::app );
	solution_file << "$NodeData" << endl;
	solution_file << "1" << endl;   // one string follows
	solution_file << "\"elastic displacement\"" << endl;
	solution_file << "1" << endl;   //  one real follows
	solution_file << "0.0" << endl;  // time [??]
	solution_file << "3" << endl;   // three integers follow
	solution_file << "0" << endl;   // time step [??]
	solution_file << "3" << endl;  // components of the gradient
	solution_file << ABCD.number_of ( tag::vertices ) << endl;  // number of values listed below
	CellIterator it = ABCD.iterator ( tag::over_vertices );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell P = *it;
		size_t i = numbering[P.core];
		solution_file << i << " " << vector_sol[i-1]<< " " << vector_sol[i-1+N] << " " << 0. << std::endl;   }
	} // just a block of code

	std::cout << "produced file square_elast.msh" << std::endl;
	
	CellIterator it = ABCD.iterator ( tag::over_vertices );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell P = *it;
		size_t i = numbering[P.core];
	    x(P) += vector_sol[i-1];	
		y(P) += vector_sol[i-1+N];
	}

	ABCD.export_msh ("square_elast_deform.msh", numbering );

	return 0;
}

/*

		FiniteElement::WithMaster * fe_core =
			dynamic_cast < FiniteElement::WithMaster * > ( fe.core );
		assert ( fe_core );
		Function::Diffeomorphism * tran = Function::core_to_diffeom ( fe_core->transf.core );
		Function::name[tran->master_coords[0].core] = "xi";
		Function::name[tran->master_coords[1].core] = "eta";

			std::cout << "psiV : " << psiV.repr() << std::endl;
			std::cout << "psiW : " << psiW.repr() << std::endl;

	std::cout << "d psiV / dx : " << d_psiV_dx.repr() << std::endl;
	std::cout << "d psiV / dy : " << d_psiV_dy.repr() << std::endl;
	std::cout << "d psiW / dx : " << d_psiW_dx.repr() << std::endl;
	std::cout << "d psiW / dy : " << d_psiW_dy.repr() << std::endl;


*/	

