

#include <iostream>
#include <cstring>
#include <vector>
#include <list>

#include "assert.h"

template <typename T>
class Tensor

// indices begin at 0
// constructors fill in values 0. automatically

{	// data:
	public:
	std::vector < T > elements;
	std::list < size_t > dimensions, cumulative_dims;
	size_t total_dim;
	
	// constructors:
	inline Tensor ( ) { };
	inline Tensor ( const std::list < size_t > dims )
	{	dimensions = dims;
		allocate_space (); };
	inline Tensor (const std::list < char > dims)
	{	dimensions = str2list ( dims );
		allocate_space ();             };
	inline Tensor ( const char dims[] )
	{	dimensions = str2list ( dims );
		allocate_space ();             };
	inline Tensor ( size_t i, size_t j,
	         size_t k, size_t l )
	{	dimensions .push_back (i);
		dimensions .push_back (j);
		dimensions .push_back (k);
		dimensions .push_back (l);
		allocate_space ();        };
	inline Tensor ( size_t i, size_t j, size_t k )
	{	dimensions .push_back (i);
		dimensions .push_back (j);
		dimensions .push_back (k);
		allocate_space ();        };
	inline Tensor ( size_t i, size_t j )
	{	dimensions .push_back (i);
		dimensions .push_back (j);
		allocate_space ();         };
	inline Tensor ( size_t i )
	{	dimensions .push_back (i);
		allocate_space ();        };
	inline ~Tensor () { };

	// methods:
	inline std::list < size_t > str2list ( const std::list < char > lc )
	{	const size_t izero = size_t ('0');
		std::list < char > ::iterator i;
		for ( i = lc .begin(); i != lc .end(); i++)
		{	assert ( *i >= '0' );
			assert ( *i <= '9' ); }
		std::list < size_t > li;
		for ( i = lc.begin(); i != lc .end(); i++)
			li .push_back ( size_t (*i) - izero );
		return li;                                   }
	inline std::list < size_t > str2list ( const char lc[] )
	{	const size_t izero = size_t ('0');
		for ( size_t i=0 ; i < std::strlen (lc); i++)
		{	assert ( lc[i] >= '0' );
			assert ( lc[i] <= '9' ); }
		std::list < size_t > li;
		for ( size_t i=0; i < std::strlen (lc); i++ )
			li .push_back ( size_t ( lc[i] ) - izero );
		return li;                                 }
	inline void allocate_space ()
	{	total_dim = 1;
		std::list < size_t > ::iterator k;
		for ( k = dimensions .begin(); k != dimensions .end(); k++)
		{	cumulative_dims .push_back ( total_dim );
			total_dim *= *k;                         }
		// elements .reserve ( total_dim );
		elements .insert ( elements .end(), total_dim, 0.0 );
		assert ( elements .size() == total_dim );                    };
	inline size_t pointer ( std::list < size_t > index ) const
	{	assert ( index .size() == dimensions .size() );
		size_t pointer = 0;
		std::list < size_t > ::const_iterator i, d, cd;
		for ( i = index .begin(), d = dimensions .begin(),
		      cd =cumulative_dims.begin();
		      i != index.end(); i++, d++, cd++)
		{	assert ( *i < *d );
			pointer += (*i) * (*cd); }
		return pointer;                                   }
	inline T& operator() ( std::list < size_t > index )
	{ return ( elements [ pointer ( index ) ] );  }
	inline const T& operator() ( std::list < size_t > index ) const
	{ return ( elements [ pointer ( index ) ] );  }
	inline T& operator[] ( const char index[] )
	{ return operator() ( str2list ( index ) );  }
	inline const T& operator() ( const char index[] ) const
	{ return operator() ( str2list ( index ) );  }
	inline T& operator() ( size_t i, size_t j,
	               size_t k, size_t l )
	{	assert ( dimensions .size() == 4 );
		std::list < size_t > index;
		index .push_back (i);
		index .push_back (j);
		index .push_back (k);
		index .push_back (l);
		return operator() ( index );        }
	inline const T& operator()( size_t i, size_t j,
	                     size_t k, size_t l ) const
	{	assert ( dimensions .size() == 4 );
		std::list < size_t > index;
		index .push_back (i);
		index .push_back (j);
		index .push_back (k);
		index .push_back (l);
		return operator() ( index );        }
	inline T& operator() ( size_t i, size_t j, size_t k )
	{	assert ( dimensions .size() == 3 );
		std::list < size_t > index;
		index .push_back (i);
		index .push_back (j);
		index .push_back (k);
		return operator() ( index );        }
	inline const T& operator() ( size_t i, size_t j, size_t k ) const
	{	assert ( dimensions .size() == 3 );
		std::list < size_t > index;
		index .push_back (i);
		index .push_back (j);
		index .push_back (k);
		return operator() ( index );        }
	inline T& operator() ( size_t i, size_t j )
	{	assert ( dimensions .size() == 2 );
		std::list < size_t > index;
		index .push_back (i);
		index .push_back (j);
		return operator() ( index );        }
	inline const T& operator() ( size_t i, size_t j ) const
	{	assert ( dimensions .size() == 2 );
		std::list < size_t > index;
		index .push_back (i);
		index .push_back (j);
		return operator() ( index );        }
	inline const T& operator() ( size_t i ) const
	{	assert ( dimensions .size() == 1 );
		std::list < size_t > index ( 1, i );
		return operator() ( index );       }
	inline T& operator[] ( size_t i )
	{	assert ( dimensions .size() == 1 );
	  std::list < size_t > index ( 1, i );
		return operator() ( index );       }

};  // end of  class Tensor



void minimize_residual
( const Tensor < double > & A, std::vector < double > & x, const std::vector < double > & b )

{	size_t n = x .size();
	assert ( A .dimensions == std::list < size_t > ( { n, n } ) );

	std::vector < double > bmAx = b, delta ( n, 0. ), Adelta ( n, 0. );
	assert ( bmAx .size() == n );   //  -Ax+b
	assert ( delta .size() == n );   //  -At (Ax-b)
	assert ( Adelta .size() == n );   //  -A At (Ax-b)

	for ( size_t i = 0; i < n; i++ )
	for ( size_t j = 0; j < n; j++ )
		bmAx [i] -= A (i,j) * x [j];

	for ( size_t i = 0; i < n; i++ )
	for ( size_t j = 0; j < n; j++ )
		delta [i] += A (j,i) * bmAx [j];

	std::cout << "delta ";
	for ( size_t i = 0; i < n; i++ )
		std::cout << " " << delta [i];
	std::cout << std::endl;
	
	for ( size_t i = 0; i < n; i++ )
	for ( size_t j = 0; j < n; j++ )
		Adelta [i] += A (i,j) * delta [j];

	double up = 0., down = 0.;
	for ( size_t i = 0; i < n; i++ )
	{	up   += Adelta [i] * bmAx [i];
		down += Adelta [i] * Adelta [i];  }

	if ( std::abs ( down ) < 1.e-10 ) return;

	const double s = up / down;

	for ( size_t i = 0; i < n; i++ )
		x [i] +=  s * delta [i];
}


void conjugate_gradient_1  // Fletcher Reeves
( std::vector < double > & x, const std::vector < double > constr,
	const Tensor < double > & grad_constr, std::vector < double > & grad,
	std::vector < double > & direc, bool first_step                      )

// performs one step for minimizing a sum of squares of m constraints, in n variables

// x is the current point, will be changed

// 'constr' has the values of each constraint at point x
// 'grad_constr' has the gradients of each individual constraint at point x

// grad  receives the gradient of the sum of squares at previous step
// and also returns the gradient of the sum of squares at x
	
// uses previous direction 'direc' and returns the new direction in the same vector
// if 'first step', do not use information about previous direction

// norm of 'grad' or of 'direc' can be used as stopping criterion
	
{	size_t n = x .size();
	assert ( direc .size() == n );
	assert ( grad .size() == n );
	size_t m = constr .size();
	assert ( grad_constr .dimensions == std::list < size_t > ( { m, n } ) );

	std::vector < double > new_grad ( n, 0. );
	assert ( new_grad .size() == n );
	Tensor < double > grad_grad ( n, n );
	for ( size_t i = 0; i <	n; i++ )
	for ( size_t k = 0; k <	m; k++ )
	{	new_grad [i] += constr [k] * grad_constr ( k, i );
		for ( size_t j = 0; j <	n; j++ )
			grad_grad ( i, j ) += grad_constr ( k, i ) * grad_constr ( k, j );  }
	
	if ( first_step ) direc = new_grad;
	else
	{	double up = 0., down = 0.;
		for ( size_t i = 0; i < n; i++ )
		{	up += new_grad [i] * new_grad [i];
			down += grad [i] * grad [i];                      }
		double alpha = up / down;
		for ( size_t i = 0; i < n; i++ ) direc [i] *= alpha;
		for ( size_t i = 0; i < n; i++ ) direc [i] += new_grad [i];  }

	double down = 0.;
	for ( size_t i = 0; i < n; i++ )
	for ( size_t j = 0; j < n; j++ )
		down += grad_grad ( i, j ) * direc [i] * direc [j];
	
	for ( size_t i = 0; i < n; i++ )
	for ( size_t j = 0; j < n; j++ )
		x [i] -= new_grad [j] * direc [j] * direc [i] / down;
	
	grad .swap ( new_grad );
}

//---------------------------------------------------------------------------------------------------//

void conjugate_gradient   // Polak Ribiere
( std::vector < double > & x, const std::vector < double > constr,
	const Tensor < double > & grad_constr, std::vector < double > & grad,
	std::vector < double > & direc, bool first_step                      )

// performs one step for minimizing a sum of squares of m constraints, in n variables

// x is the current point, will be changed

// 'constr' has the values of each constraint at point x
// 'grad_constr' has the gradients of each individual constraint at point x

// grad  receives the gradient of the sum of squares at previous step
// and also returns the gradient of the sum of squares at x
	
// uses previous direction 'direc' and returns the new direction in the same vector
// if 'first step', do not use information about previous direction

// norm of 'grad' or of 'direc' can be used as stopping criterion
	
{	size_t n = x .size();
	assert ( direc .size() == n );
	assert ( grad .size() == n );
	size_t m = constr .size();
	assert ( grad_constr .dimensions == std::list < size_t > ( { m, n } ) );

	std::vector < double > new_grad ( n, 0. );
	assert ( new_grad .size() == n );
	Tensor < double > grad_grad ( n, n );
	for ( size_t i = 0; i <	n; i++ )
	for ( size_t k = 0; k <	m; k++ )
	{	new_grad [i] += constr [k] * grad_constr ( k, i );
		for ( size_t j = 0; j <	n; j++ )
			grad_grad ( i, j ) += grad_constr ( k, i ) * grad_constr ( k, j );  }
	
	if ( first_step ) direc = new_grad;
	else
	{	double up = 0., down = 0.;
		for ( size_t i = 0; i < n; i++ )
		{	up += new_grad [i] * ( new_grad [i] - grad [i] );
			down += grad [i] * grad [i];                      }
		double alpha = up / down;
		for ( size_t i = 0; i < n; i++ ) direc [i] *= alpha;
		for ( size_t i = 0; i < n; i++ ) direc [i] += new_grad [i];  }

	double down = 0.;
	for ( size_t i = 0; i < n; i++ )
	for ( size_t j = 0; j < n; j++ )
		down += grad_grad ( i, j ) * direc [i] * direc [j];
	
	for ( size_t i = 0; i < n; i++ )
	for ( size_t j = 0; j < n; j++ )
		x [i] -= new_grad [j] * direc [j] * direc [i] / down;
	
	grad .swap ( new_grad );
}

//---------------------------------------------------------------------------------------------------//


void simple_gradient
( std::vector < double > & x, const std::vector < double > constr,
	const Tensor < double > & grad_constr, std::vector < double > & grad,
	std::vector < double > & direc, bool first_step                      )

// performs one step for minimizing a sum of squares of m constraints, in n variables

// x is the current point, will be changed

// 'constr' has the values of each constraint at point x
// 'grad_constr' has the gradients of each individual constraint at point x

// grad  receives the gradient of the sum of squares at previous step
// and also returns the gradient of the sum of squares at x
	
// returns the direction in 'direc'
// 'first step' is not used

// norm of 'grad' or of 'direc' can be used as stopping criterion
	
{	size_t n = x .size();
	assert ( direc .size() == n );
	assert ( grad .size() == n );
	size_t m = constr .size();
	assert ( grad_constr .dimensions == std::list < size_t > ( { m, n } ) );

	std::vector < double > new_grad ( n, 0. );
	assert ( new_grad .size() == n );
	Tensor < double > grad_grad ( n, n );
	for ( size_t i = 0; i <	n; i++ )
	for ( size_t k = 0; k <	m; k++ )
	{	new_grad [i] += constr [k] * grad_constr ( k, i );
		for ( size_t j = 0; j <	n; j++ )
			grad_grad ( i, j ) += grad_constr ( k, i ) * grad_constr ( k, j );  }
	
	direc = new_grad;

	double down = 0.;
	for ( size_t i = 0; i < n; i++ )
	for ( size_t j = 0; j < n; j++ )
		down += grad_grad ( i, j ) * direc [i] * direc [j];
	
	for ( size_t i = 0; i < n; i++ )
	for ( size_t j = 0; j < n; j++ )
		x [i] -= new_grad [j] * direc [j] * direc [i] / down;
	
	grad .swap ( new_grad );
}

//---------------------------------------------------------------------------------------------------//

int main ()

{	Tensor < double > points ( 4, 3 );
	points (0,0) =  1.;   points (0,1) = 0.;    points (0,2) = 0.;
	points (1,0) =  0.5;  points (1,1) = 0.86;  points (1,2) = 0.;
	points (2,0) =  0.5;  points (2,1) = 0.3;   points (2,2) = 0.6;
	points (3,0) = -0.6;  points (3,1) = 0.3;   points (3,2) = 0.5;

	std::list < size_t > dims = points .dimensions;
	assert ( dims .size() == 2 );
	std::list < size_t > ::iterator it = dims .begin();
	assert ( it != dims .end() );
	size_t n_points = *it;
	it ++;  assert ( it != dims .end() );
	size_t geom_dim = *it;

	std::vector < double > grad ( geom_dim ), direc ( geom_dim ), x ( geom_dim, 0.2 );
	
	std::vector < double > constr ( n_points );
	for ( size_t i = 0; i < n_points; i++ )
	{	constr [i] = -1.;
		for ( size_t j = 0; j < geom_dim; j++ )
			constr [i] += ( x [j] - points ( i, j ) ) * ( x [j] - points ( i, j ) );  }

	for ( size_t i = 0; i < geom_dim; i++ ) std::cout << " " << x [i];
	std::cout << " ---> ";
	for ( size_t i = 0; i < n_points; i++ ) std::cout << " " << constr [i] + 1.;
	std::cout << std::endl;

	Tensor < double > grad_constr ( n_points, geom_dim );
	for ( size_t i = 0; i < n_points; i++ )
	for ( size_t j = 0; j < geom_dim; j++ )
		grad_constr ( i, j ) = 2. * ( x [j] - points ( i, j ) );

	// keep old x for stopping criterion
	std::vector < double > old_x = x;
	conjugate_gradient ( x, constr, grad_constr, grad, direc, true );
	// last argument true means "first step"

	for ( size_t iter = 0; iter < 20; iter++ )
	{	// we need to compute again the values of the constraints
		for ( size_t i = 0; i < n_points; i++ )
		{	constr [i] = -1.;
			for ( size_t j = 0; j < geom_dim; j++ )
				constr [i] += ( x [j] - points ( i, j ) ) * ( x [j] - points ( i, j ) );  }
		for ( size_t i = 0; i < geom_dim; i++ ) std::cout << " " << x [i];
		std::cout << " ---> ";
		for ( size_t i = 0; i < n_points; i++ ) std::cout << " " << constr [i] + 1.;
		std::cout << std::endl;

		// test stopping criterion :
		double norm = 0.;
		for ( size_t i = 0; i < geom_dim; i++ )
		{	double tmp = std::abs ( x [i] - old_x [i] );
			if ( tmp > norm ) norm = tmp;                }
		if ( norm < 1.e-4 ) break;

		// we need to compute again the gradients
		for ( size_t i = 0; i < n_points; i++ )
		for ( size_t j = 0; j < geom_dim; j++ )
			grad_constr ( i, j ) = 2. * ( x [j] - points ( i, j ) );

		old_x = x;
		conjugate_gradient ( x, constr, grad_constr, grad, direc, false );                }
		// last argument false means "not first step"

}



void main1 ()

{	Tensor < double > A (3,3);
	A (0,0) = 1.;  A (1,0) =  0.5;  A (0,1) = 0.5;  A (1,1) = 1.;
	A (0,2) = 3.;  A (1,2) = -0.2;  A (2,0) = 0.;   A (2,1) = 0.;  A (2,2) = 2.;

	std::vector < double > b {0.,1.,2.}, x {0.,0.,0.};

	size_t n = x .size();
	size_t m = b .size();
	assert ( A .dimensions == std::list < size_t > ( { m, n } ) );

	std::vector < double > resid (m);
	for ( size_t i = 0; i < m; i++ )
	{	resid [i] = -b [i];
		for ( size_t j = 0; j < n; j++ )
			resid [i] += A ( i, j ) * x [j];  }

	std::vector < double > grad (n), direc (n);
	
	for ( size_t i = 0; i < n; i++ )
		std::cout << " " << x [i];
	std::cout << std::endl;
	conjugate_gradient ( x, resid, A, grad, direc, true );
	// last argument true means "first step"

	for ( size_t iter = 0; iter < 20; iter++ )
	{	for ( size_t i = 0; i < n; i++ )
			std::cout << " " << x [i];
		std::cout << std::endl;
		// we need to compute again the values of the constraints (resid)
		// no need to update gradients, they are constant
		double norm = 0.;
		for ( size_t i = 0; i < m; i++ )
		{	resid [i] = -b [i];
			for ( size_t j = 0; j < n; j++ )
				resid [i] += A ( i, j ) * x [j];
			double tmp = std::abs ( resid [i] );
			if ( tmp > norm ) norm = tmp;        }
		if ( norm < 1.e-7 ) break;
		conjugate_gradient ( x, resid, A, grad, direc, false );  }
		// last argument false means "not first step"

	for ( size_t i = 0; i < n; i++ )
		std::cout << " " << x [i];
	std::cout << std::endl;

	std::cout << "---------------------------" << std::endl;
	for ( size_t i = 0; i < m; i++ )
		std::cout << " " << b [i];
	std::cout << std::endl;

	for ( size_t i = 0; i < m; i++ )
	{	double tmp = 0;
		for ( size_t j = 0; j < n; j++ )
			tmp += A (i,j) * x [j];
		std::cout << " " << tmp;         }
	std::cout << std::endl;

}
