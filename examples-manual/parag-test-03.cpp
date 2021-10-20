
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/OrderingMethods>

int main ()

{
	Eigen::SparseMatrix < double > A ( 6, 5 );
	A.coeffRef ( 0, 0 ) =  2.;
	A.coeffRef ( 0, 1 ) = -1.;
	A.coeffRef ( 0, 4 ) = -1.;
	A.coeffRef ( 1, 0 ) = -1.;
	A.coeffRef ( 1, 1 ) =  2.;
	A.coeffRef ( 1, 2 ) = -1.;
	A.coeffRef ( 2, 1 ) = -1.;
	A.coeffRef ( 2, 2 ) =  2.;
	A.coeffRef ( 2, 3 ) = -1.;
	A.coeffRef ( 3, 2 ) = -1.;
	A.coeffRef ( 3, 3 ) =  2.;
	A.coeffRef ( 3, 4 ) = -1.;
	A.coeffRef ( 4, 0 ) = -1.;
	A.coeffRef ( 4, 3 ) = -1.;
	A.coeffRef ( 4, 4 ) =  2.;
	A.coeffRef ( 5, 0 ) =  1.;
	A.coeffRef ( 5, 1 ) =  1.;
	A.coeffRef ( 5, 2 ) =  1.;
	A.coeffRef ( 5, 3 ) =  1.;
	A.coeffRef ( 5, 4 ) =  1.;

	A.makeCompressed();

	Eigen::VectorXd b(6), x(5);
	b(0) =  1.;
	b(1) =  1.;
	b(2) =  1.;
	b(3) =  1.;
	b(4) = -4.;
	b(5) =  0.;

	Eigen::SparseQR < Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;

	solver.compute ( A );
	if ( solver.info() != Eigen::Success )
	{	std::cout << "compute failed" << std::endl;
		exit ( 0 );                                  }

	x = solver.solve ( b );
	if ( solver.info() != Eigen::Success )
	{	std::cout << "solve failed" << std::endl;
		exit ( 0 );                                  }

	std::cout << "x" << std::endl << x << std::endl;
	std::cout << " Ax" << std::endl << A*x << std::endl;
}
	
