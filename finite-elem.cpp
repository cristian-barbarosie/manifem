
// finite-elem.cpp 2021.12.14

//   This file is part of maniFEM, a C++ library for meshes and finite elements on manifolds.

//   Copyright 2019, 2020, 2021 Cristian Barbarosie cristian.barbarosie@gmail.com
//   https://github.com/cristian-barbarosie/manifem

//   ManiFEM is free software: you can redistribute it and/or modify it
//   under the terms of the GNU Lesser General Public License as published
//   by the Free Software Foundation, either version 3 of the License
//   or (at your option) any later version.

//   ManiFEM is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//   See the GNU Lesser General Public License for more details.

//   You should have received a copy of the GNU Lesser General Public License
//   along with maniFEM.  If not, see <https://www.gnu.org/licenses/>.

#include "math.h"
#include <sstream>
#include <typeinfo>

#include "iterator.h"
#include "function.h"
#include "finite-elem.h"

using namespace maniFEM;


Integrator::Gauss::Gauss ( const tag::gauss_quadrature & q,
                           const tag::FromFiniteElementWithMaster &, FiniteElement & fe )

: Integrator::Core ()

// J.E. Flaherty, Finite Element Analysis, Lecture Notes : Spring 2000
// http://manifem.rd.ciencias.ulisboa.pt/flaherty-06.pdf
	
// E.B. Becker, G.F. Carey, J.T. Oden, Finite Elements, an introduction, vol 1

{ FiniteElement::WithMaster * fe_core = tag::Util::assert_cast
		< FiniteElement::Core*, FiniteElement::WithMaster* > ( fe.core );
	Manifold master_manifold = fe_core->master_manif;
	switch ( q )
	{	case tag::seg_2 :  // Gauss quadrature with two points on a segment
		{	assert ( dynamic_cast < FiniteElement::WithMaster::Segment* > ( fe.core ) );
			Function t = master_manifold .coordinates();
			assert ( t .nb_of_components() == 1 );
			this->points .reserve (2);  this->weights .reserve (2);
			const double one_over_sqrt_three = 1. / std::sqrt (3.);
			Cell Gauss_2_A ( tag::vertex );  t ( Gauss_2_A ) = - one_over_sqrt_three;
			Cell Gauss_2_B ( tag::vertex );  t ( Gauss_2_B ) =   one_over_sqrt_three;
			this->points .push_back ( Gauss_2_A );
			this->weights .push_back ( 1. );
			this->points .push_back ( Gauss_2_B );
			this->weights .push_back ( 1. );
			break;                                                                    }
		case tag::seg_3 :  // Gauss quadrature with three points on a segment
		{	assert ( dynamic_cast < FiniteElement::WithMaster::Segment* > ( fe.core ) );
			Function t = master_manifold .coordinates();
			assert ( t .nb_of_components() == 1 );
			this->points .reserve (3);  this->weights .reserve (3);
			const double sqrt3over5 = std::sqrt (0.6);
			Cell Gauss_3_A ( tag::vertex );  t ( Gauss_3_A ) = - sqrt3over5;
			Cell Gauss_3_B ( tag::vertex );  t ( Gauss_3_B ) =   0.;
			Cell Gauss_3_C ( tag::vertex );  t ( Gauss_3_C ) =   sqrt3over5;
			this->points .push_back ( Gauss_3_A );
			this->weights .push_back ( 5./9. );
			this->points .push_back ( Gauss_3_B );
			this->weights .push_back ( 8./9. );
			this->points .push_back ( Gauss_3_C );
			this->weights .push_back ( 5./9. );
			break;                                                                    }
		case tag::seg_4 :  // Gauss quadrature with four points on a segment
		{	assert ( dynamic_cast < FiniteElement::WithMaster::Segment* > ( fe.core ) );
			Function t = master_manifold .coordinates();
			assert ( t .nb_of_components() == 1 );
			this->points .reserve (4);  this->weights .reserve (4);
			Cell Gauss_4_A ( tag::vertex );  t ( Gauss_4_A ) = - 0.861136311594053;
			Cell Gauss_4_B ( tag::vertex );  t ( Gauss_4_B ) = - 0.339981043584856;
			Cell Gauss_4_C ( tag::vertex );  t ( Gauss_4_C ) =   0.339981043584856;
			Cell Gauss_4_D ( tag::vertex );  t ( Gauss_4_D ) =   0.861136311594053;
			this->points .push_back ( Gauss_4_A );
			this->weights .push_back ( 0.347854845137454 );
			this->points .push_back ( Gauss_4_B );
			this->weights .push_back ( 0.652145154862546 );
			this->points .push_back ( Gauss_4_C );
			this->weights .push_back ( 0.652145154862546 );
			this->points .push_back ( Gauss_4_D );
			this->weights .push_back ( 0.347854845137454 );
			break;                                                                    }
		case tag::seg_5 :  // Gauss quadrature with five points on a segment
		{	assert ( dynamic_cast < FiniteElement::WithMaster::Segment* > ( fe.core ) );
			Function t = master_manifold .coordinates();
			assert ( t .nb_of_components() == 1 );
			this->points .reserve (5);  this->weights .reserve (5);
			Cell Gauss_5_A ( tag::vertex );  t ( Gauss_5_A ) = - 0.906179845938664;
			Cell Gauss_5_B ( tag::vertex );  t ( Gauss_5_B ) = - 0.538469310105683;
			Cell Gauss_5_C ( tag::vertex );  t ( Gauss_5_C ) =   0.;
			Cell Gauss_5_D ( tag::vertex );  t ( Gauss_5_D ) =   0.538469310105683;
			Cell Gauss_5_E ( tag::vertex );  t ( Gauss_5_E ) =   0.906179845938664;
			this->points .push_back ( Gauss_5_A );
			this->weights .push_back ( 0.236926885056189 );
			this->points .push_back ( Gauss_5_B );
			this->weights .push_back ( 0.478628670499366 );
			this->points .push_back ( Gauss_5_C );
			this->weights .push_back ( 0.568888888888889 );
			this->points .push_back ( Gauss_5_D );
			this->weights .push_back ( 0.478628670499366 );
			this->points .push_back ( Gauss_5_E );
			this->weights .push_back ( 0.236926885056189 );
			break;                                                                    }
		case tag::seg_6 :  // Gauss quadrature with six points on a segment
		{	assert ( dynamic_cast < FiniteElement::WithMaster::Segment* > ( fe.core ) );
			Function t = master_manifold .coordinates();
			assert ( t .nb_of_components() == 1 );
			this->points .reserve (6);  this->weights .reserve (6);
			Cell Gauss_6_A ( tag::vertex );  t ( Gauss_6_A ) = - 0.932469514203152;
			Cell Gauss_6_B ( tag::vertex );  t ( Gauss_6_B ) = - 0.661209386466265;
			Cell Gauss_6_C ( tag::vertex );  t ( Gauss_6_C ) = - 0.238619186083197;
			Cell Gauss_6_D ( tag::vertex );  t ( Gauss_6_D ) =   0.238619186083197;
			Cell Gauss_6_E ( tag::vertex );  t ( Gauss_6_E ) =   0.661209386466265;
			Cell Gauss_6_F ( tag::vertex );  t ( Gauss_6_F ) =   0.932469514203152;
			this->points .push_back ( Gauss_6_A );
			this->weights .push_back ( 0.171324492379170 );
			this->points .push_back ( Gauss_6_B );
			this->weights .push_back ( 0.360761573048139 );
			this->points .push_back ( Gauss_6_C );
			this->weights .push_back ( 0.467913934572691 );
			this->points .push_back ( Gauss_6_D );
			this->weights .push_back ( 0.467913934572691 );
			this->points .push_back ( Gauss_6_E );
			this->weights .push_back ( 0.360761573048139 );
			this->points .push_back ( Gauss_6_F );
			this->weights .push_back ( 0.171324492379170 );
			break;                                                                    }
		case tag::tri_3 :  // Gauss quadrature with three points on a triangle
		{	assert ( dynamic_cast < FiniteElement::WithMaster::Triangle* > ( fe.core ) );
			Function xi_eta = master_manifold .coordinates();
			assert ( xi_eta .nb_of_components() == 2 );
			Function xi = xi_eta [0], eta = xi_eta [1];
			this->points .reserve (3);  this->weights .reserve (3);
			// points
			Cell Gauss_3_O ( tag::vertex );  xi ( Gauss_3_O ) = 1./6.;  eta ( Gauss_3_O ) = 1./6.;
			Cell Gauss_3_A ( tag::vertex );  xi ( Gauss_3_A ) = 2./3.;  eta ( Gauss_3_A ) = 1./6.;
			Cell Gauss_3_B ( tag::vertex );  xi ( Gauss_3_B ) = 1./6.;  eta ( Gauss_3_B ) = 2./3.;
			// weights
			double Gw_3_O = 1./6., Gw_3_A = 1./6., Gw_3_B = 1./6.;
			// in the book of Oden we have the same weights but a different distribution of points
			// instead of  2/3 1/6 1/6  Oden has  middle of segments 0 0.5 0.5
			// maybe the points are not unique ?  see below
			this->points .push_back ( Gauss_3_O );
			this->weights .push_back ( Gw_3_O );
			this->points .push_back ( Gauss_3_A );
			this->weights .push_back ( Gw_3_A );
			this->points .push_back ( Gauss_3_B );
			this->weights .push_back ( Gw_3_B );
			break;                                                                                 }
		case tag::tri_3_Oden :
		// Gauss quadrature with three points on a triangle, points presented in the book of Oden
		{	assert ( dynamic_cast < FiniteElement::WithMaster::Triangle* > ( fe.core ) );
			Function xi_eta = master_manifold .coordinates();
			assert ( xi_eta .nb_of_components() == 2 );
			Function xi = xi_eta [0], eta = xi_eta [1];
			this->points .reserve (3);  this->weights .reserve (3);
			// points
			Cell Gauss_3_O ( tag::vertex );  xi ( Gauss_3_O ) = 0.5;  eta ( Gauss_3_O ) = 0.5;
			Cell Gauss_3_A ( tag::vertex );  xi ( Gauss_3_A ) = 0.;   eta ( Gauss_3_A ) = 0.5;
			Cell Gauss_3_B ( tag::vertex );  xi ( Gauss_3_B ) = 0.5;  eta ( Gauss_3_B ) = 0.;
			// weights
			double Gw_3_O = 1./6., Gw_3_A = 1./6., Gw_3_B = 1./6.;
			this->points .push_back ( Gauss_3_O );
			this->weights .push_back ( Gw_3_O );
			this->points .push_back ( Gauss_3_A );
			this->weights .push_back ( Gw_3_A );
			this->points .push_back ( Gauss_3_B );
			this->weights .push_back ( Gw_3_B );
			break;                                                                             }
		case tag::tri_4 :  // Gauss quadrature with four points on a triangle
		{	assert ( dynamic_cast < FiniteElement::WithMaster::Triangle* > ( fe.core ) );
			Function xi_eta = master_manifold .coordinates();
			assert ( xi_eta .nb_of_components() == 2 );
			Function xi = xi_eta [0], eta = xi_eta [1];
			this->points .reserve (4);  this->weights .reserve (4);
			// points
			Cell Gauss_4_c ( tag::vertex );  xi ( Gauss_4_c ) = 1./3.;  eta ( Gauss_4_c ) = 1./3.;
			Cell Gauss_4_O ( tag::vertex );  xi ( Gauss_4_O ) = 0.2;    eta ( Gauss_4_O ) = 0.2;
			Cell Gauss_4_A ( tag::vertex );  xi ( Gauss_4_A ) = 0.6;    eta ( Gauss_4_A ) = 0.2;
			Cell Gauss_4_B ( tag::vertex );  xi ( Gauss_4_B ) = 0.2;    eta ( Gauss_4_B ) = 0.6;
			// weights
			double Gw_4_c = -27./96.;  // yes, negative weight
			double Gw_4_O = 25./96., Gw_4_A = 25./96., Gw_4_B = 25./96.;
			// in the book of Oden we have the same weights but a different distribution of points
			// instead of  0.2 0.2 0.6  we have  2/15 2/15 11/15
			// maybe the points are not unique ?  see below
			this->points .push_back ( Gauss_4_c );
			this->weights .push_back ( Gw_4_c );
			this->points .push_back ( Gauss_4_O );
			this->weights .push_back ( Gw_4_O );
			this->points .push_back ( Gauss_4_A );
			this->weights .push_back ( Gw_4_A );
			this->points .push_back ( Gauss_4_B );
			this->weights .push_back ( Gw_4_B );
			break;                                                                                 }
		case tag::tri_4_Oden :
		// Gauss quadrature with four points on a triangle, points presented in the book of Oden
		{	assert ( dynamic_cast < FiniteElement::WithMaster::Triangle* > ( fe.core ) );
			Function xi_eta = master_manifold .coordinates();
			assert ( xi_eta .nb_of_components() == 2 );
			Function xi = xi_eta [0], eta = xi_eta [1];
			this->points .reserve (4);  this->weights .reserve (4);
			// points
			Cell Gauss_4_c ( tag::vertex );  xi(Gauss_4_c) =  1./ 3.;  eta(Gauss_4_c) =  1./ 3.;
			Cell Gauss_4_O ( tag::vertex );  xi(Gauss_4_O) =  2./15.;  eta(Gauss_4_O) =  2./15.;
			Cell Gauss_4_A ( tag::vertex );  xi(Gauss_4_A) = 11./15.;  eta(Gauss_4_A) =  2./15.;
			Cell Gauss_4_B ( tag::vertex );  xi(Gauss_4_B) =  2./15.;  eta(Gauss_4_B) = 11./15.;
			// weights
			double Gw_4_c = -27./96.;  // yes, negative weight
			double Gw_4_O = 25./96., Gw_4_A = 25./96., Gw_4_B = 25./96.;
			this->points .push_back ( Gauss_4_c );
			this->weights .push_back ( Gw_4_c );
			this->points .push_back ( Gauss_4_O );
			this->weights .push_back ( Gw_4_O );
			this->points .push_back ( Gauss_4_A );
			this->weights .push_back ( Gw_4_A );
			this->points .push_back ( Gauss_4_B );
			this->weights .push_back ( Gw_4_B );
			break;                                                                                }
		case tag::tri_6 :  // Gauss quadrature with six points on a triangle
		{	assert ( dynamic_cast < FiniteElement::WithMaster::Triangle* > ( fe.core ) );
			Function xi_eta = master_manifold .coordinates();
			assert ( xi_eta .nb_of_components() == 2 );
			Function xi = xi_eta [0], eta = xi_eta [1];
			this->points .reserve (6);  this->weights .reserve (6);
			// points
			Cell Gauss_6_O ( tag::vertex );   xi  ( Gauss_6_O )  = 0.091576213509771;
			                                  eta ( Gauss_6_O )  = 0.091576213509771;
			Cell Gauss_6_A ( tag::vertex );   xi  ( Gauss_6_A )  = 0.816847572980459;
			                                  eta ( Gauss_6_A )  = 0.091576213509771;
			Cell Gauss_6_B ( tag::vertex );   xi  ( Gauss_6_B )  = 0.091576213509771;
			                                  eta ( Gauss_6_B )  = 0.816847572980459;
			Cell Gauss_6_OA ( tag::vertex );  xi  ( Gauss_6_OA ) = 0.445948490915965;
			                                  eta ( Gauss_6_OA ) = 0.108103018168070;
			Cell Gauss_6_OB ( tag::vertex );  xi  ( Gauss_6_OB ) = 0.108103018168070;
			                                  eta ( Gauss_6_OB ) = 0.445948490915965;
			Cell Gauss_6_AB ( tag::vertex );  xi  ( Gauss_6_AB ) = 0.445948490915965;
			                                  eta ( Gauss_6_AB ) = 0.445948490915965;
			// weights
			double Gw_6_O = 0.054975871827661, Gw_6_A = 0.054975871827661,
			      Gw_6_B = 0.054975871827661;
			double Gw_6_OA = 0.1116907948390055, Gw_6_OB = 0.1116907948390055,
			       Gw_6_AB = 0.1116907948390055;
			this->points .push_back ( Gauss_6_O );
			this->weights .push_back ( Gw_6_O );
			this->points .push_back ( Gauss_6_A );
			this->weights .push_back ( Gw_6_A );
			this->points .push_back ( Gauss_6_B );
			this->weights .push_back ( Gw_6_B );
			this->points .push_back ( Gauss_6_OA );
			this->weights .push_back ( Gw_6_OA );
			this->points .push_back ( Gauss_6_OB );
			this->weights .push_back ( Gw_6_OB );
			this->points .push_back ( Gauss_6_AB );
			this->weights .push_back ( Gw_6_AB );
			break;                                                                      }
		case tag::quad_4 :  // Gauss quadrature with four points on a quadrangle
		{	assert ( dynamic_cast < FiniteElement::WithMaster::Quadrangle* > ( fe.core ) );
			Function xi_eta = master_manifold .coordinates();
			assert ( xi_eta .nb_of_components() == 2 );
			Function xi = xi_eta [0], eta = xi_eta [1];
			this->points .reserve (4);  this->weights .reserve (4);
			const double sqrt_of_one_third = std::sqrt (1./3.);
			Cell Gauss_4_a ( tag::vertex );   xi  ( Gauss_4_a ) = -sqrt_of_one_third;
			                                  eta ( Gauss_4_a ) = -sqrt_of_one_third;
			Cell Gauss_4_b ( tag::vertex );   xi  ( Gauss_4_b ) = -sqrt_of_one_third;
			                                  eta ( Gauss_4_b ) =  sqrt_of_one_third;
			Cell Gauss_4_c ( tag::vertex );   xi  ( Gauss_4_c ) =  sqrt_of_one_third;
			                                  eta ( Gauss_4_c ) = -sqrt_of_one_third;
			Cell Gauss_4_d ( tag::vertex );   xi  ( Gauss_4_d ) =  sqrt_of_one_third;
			                                  eta ( Gauss_4_d ) =  sqrt_of_one_third;
			this->points .push_back ( Gauss_4_a );
			this->weights .push_back ( 1. );
			this->points .push_back ( Gauss_4_b );
			this->weights .push_back ( 1. );
			this->points .push_back ( Gauss_4_c );
			this->weights .push_back ( 1. );
			this->points .push_back ( Gauss_4_d );
			this->weights .push_back ( 1. );
			break;                                                                    }
		case tag::quad_9 :  // Gauss quadrature with nine points on a quadrangle
		{	assert ( dynamic_cast < FiniteElement::WithMaster::Quadrangle* > ( fe.core ) );
			Function xi_eta = master_manifold .coordinates();
			assert ( xi_eta .nb_of_components() == 2 );
			Function xi = xi_eta [0], eta = xi_eta [1];
			this->points .reserve (9);  this->weights .reserve( 9);
			// points :
			const double sqrt3over5 = std::sqrt (0.6);
			Cell Gauss_9_a ( tag::vertex );   xi  ( Gauss_9_a ) = -sqrt3over5;
			                                  eta ( Gauss_9_a ) = -sqrt3over5;
			Cell Gauss_9_b ( tag::vertex );   xi  ( Gauss_9_b ) =  0.;
			                                  eta ( Gauss_9_b ) = -sqrt3over5;
			Cell Gauss_9_c ( tag::vertex );   xi  ( Gauss_9_c ) =  sqrt3over5;
			                                  eta ( Gauss_9_c ) = -sqrt3over5;
			Cell Gauss_9_d ( tag::vertex );   xi  ( Gauss_9_d ) = -sqrt3over5;
			                                  eta ( Gauss_9_d ) =  0.;
			Cell Gauss_9_e ( tag::vertex );   xi  ( Gauss_9_e ) =  0.;
			                                  eta ( Gauss_9_e ) =  0.;
			Cell Gauss_9_f ( tag::vertex );   xi  ( Gauss_9_f ) =  sqrt3over5;
			                                  eta ( Gauss_9_f ) =  0.;
			Cell Gauss_9_g ( tag::vertex );   xi  ( Gauss_9_g ) = -sqrt3over5;
			                                  eta ( Gauss_9_g ) =  sqrt3over5;
			Cell Gauss_9_h ( tag::vertex );   xi  ( Gauss_9_h ) =  0.;
			                                  eta ( Gauss_9_h ) =  sqrt3over5;
			Cell Gauss_9_i ( tag::vertex );   xi  ( Gauss_9_i ) =  sqrt3over5;
			                                  eta ( Gauss_9_i ) =  sqrt3over5;
			//  a {-sqrt3over5,-sqrt3over5}, b {0, -sqrt3over5}, c {sqrt3over5, -sqrt3over5},
			//  d {-sqrt3over5,0}, e {0,0}, f {sqrt3over5, 0}, g {-sqrt3over5,sqrt3over5},
			//  h {0,sqrt3over5}, i {sqrt3over5, sqrt3over5}
			// weights :
			//  25./81., 40./81., 25./81., 40./81, 64./81., 40./81, 25./81., 40./81, 25./81.
			this->points .push_back ( Gauss_9_a );
			this->weights .push_back ( 25./81. );
			this->points .push_back ( Gauss_9_b );
			this->weights .push_back ( 40./81. );
			this->points .push_back ( Gauss_9_c );
			this->weights .push_back ( 25./81. );
			this->points .push_back ( Gauss_9_d );
			this->weights .push_back ( 40./81. );
			this->points .push_back ( Gauss_9_e );
			this->weights .push_back ( 64./81. );
			this->points .push_back ( Gauss_9_f );
			this->weights .push_back ( 40./81. );
			this->points .push_back ( Gauss_9_g );
			this->weights .push_back ( 25./81. );
			this->points .push_back ( Gauss_9_h );
			this->weights .push_back ( 40./81. );
			this->points .push_back ( Gauss_9_i );
			this->weights .push_back ( 25./81. );
			break;                                                                         }
		default :	assert ( false );

	}  // end of  switch ( q )

}  // end of  Integrator::Gauss::Gauss
	
//-----------------------------------------------------------------------------------------//


void Integrator::Gauss::pre_compute ( const std::vector < Function > & bf,
                                      const std::vector < Function > & v  )
// virtual from Integrator::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "pre_compute is not necessary for Gauss integrators" << std::endl;
	exit ( 1 );                                                                                 }


std::vector < double > Integrator::Gauss::retrieve_precomputed
( const Function & bf, const Function & psi )  // virtual from Integrator::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Gauss integrators do not pre-compute expressions" << std::endl;
	exit ( 1 );                                                                                 }


std::vector < double > Integrator::Gauss::retrieve_precomputed
( const Function & bf1, const Function & psi1,
  const Function & bf2, const Function & psi2 )  // virtual from Integrator::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Gauss integrators do not pre-compute expressions" << std::endl;
	exit ( 1 );                                                                                 }


void Integrator::UFL_FFC::pre_compute ( const std::vector < Function > & bf,
                                        const std::vector < Function > & v  )
// virtual from Integrator::Core

// bf is an arbitrary basis function (Function::MereSymbol) in the finite element
// we prepare computations for fast evaluation of integrals of expressions listed in 'v'

// we delegate this analysis to the finite element

{	this->fe .core->pre_compute ( bf, v );  }


std::vector < double > Integrator::UFL_FFC::retrieve_precomputed
( const Function & bf, const Function & psi )  // virtual from Integrator::Core

{	FiniteElement::StandAlone::TypeOne * fesa = tag::Util::assert_cast
		< FiniteElement::Core *, FiniteElement::StandAlone::TypeOne * > ( this->fe .core);
	assert ( fesa->dummy_bf .size() >= 1 );
	assert ( fesa->dummy_bf [0] .core == bf .core );
	assert ( fesa->result_of_integr .size() == 1 );
	assert ( fesa->basis_numbering [ psi .core ] < fesa->result_of_integr [0] .size() );
	std::vector < double > res ( fesa->selector .size() );
	std::vector < double > & r = fesa->result_of_integr
		[0] [ fesa->basis_numbering [ psi .core ] ];
	for ( size_t k = 0; k < fesa->selector .size(); k++ )
		res [k] = r [ fesa->selector [k] ];
	return res;                                                                         }
	// return Integrator::Result ( fesa->selector,
	//  	    fesa->result_of_integr [ 0 ] [ fesa->basis_numbering [ psi .core ] ]  )


std::vector < double > Integrator::UFL_FFC::retrieve_precomputed
( const Function & bf1, const Function & bf2,
  const Function & psi1, const Function & psi2 )  // virtual from Integrator::Core

{	FiniteElement::StandAlone::TypeOne * fesa = tag::Util::assert_cast
		< FiniteElement::Core *, FiniteElement::StandAlone::TypeOne * > ( this->fe .core);
	assert ( fesa->dummy_bf .size() >= 2 );
	assert ( fesa->dummy_bf [0] .core == bf1 .core );
	assert ( fesa->dummy_bf [1] .core == bf2 .core );
	assert ( fesa->basis_numbering [ psi1 .core ] < fesa->result_of_integr .size()  );
	assert ( fesa->basis_numbering [ psi2 .core ] <
					 fesa->result_of_integr [ fesa->basis_numbering [ psi1 .core ] ] .size() );
	// std::cout << "case " << fesa->cas << ", bf " << fesa->basis_numbering [ psi1 .core ] << " "
	//           << fesa->basis_numbering [ psi2 .core ] << ", sel " << fesa->selector [0]
	//           << std::endl;
	std::vector < double > res ( fesa->selector .size() );
	std::vector < double > & r = fesa->result_of_integr
		[ fesa->basis_numbering [ psi1 .core ] ] [ fesa->basis_numbering [ psi2 .core ] ];
	for ( size_t k = 0; k < fesa->selector .size(); k++ )
		res [k] = r [ fesa->selector [k] ];
	return res;                                                                         }
	// return Integrator::Result ( fesa->selector,
	//  	    fesa->result_of_integr [ fesa->basis_numbering [ bf1 .core ] ]
	//                              [ fesa->basis_numbering [ bf2 .core ] ]  )
	
//-----------------------------------------------------------------------------------------//


double Integrator::Gauss::action ( Function f, const FiniteElement & fe )
// virtual from Integrator::Core

// assumes the finite element is already docked on a cell
// thus, fe_core->transf is well defined

{	FiniteElement::WithMaster * fe_core = tag::Util::assert_cast
		< FiniteElement::Core*, FiniteElement::WithMaster* > ( fe .core );
	assert ( fe_core->transf .core );
	Function::Map * tran = dynamic_cast < Function::Map* > ( fe_core->transf .core );
	assert ( tran );
	// in the above we must use dynamic_cast
	// below, with -DNDEBUG, assert_cast calls a static_cast which does not compile
	//   Function::Map * tran = tag::Util::assert_cast
	// 	   < Function::Core*, Function::Map* > ( fe_core->transf .core );
	// classes Function::Core and Function::Map are not directly related
	// perhaps we should have used virtual inheritance ?

	size_t geom_dim = tran->geom_coords .nb_of_components();
	assert ( geom_dim == tran->back_geom_coords .nb_of_components() );
	for ( size_t i = 0; i < geom_dim; i++ )
		f = f .replace ( tran->geom_coords [i], tran->back_geom_coords [i] );

	double res = 0.;
	std::vector < double > ::iterator it_weight = this->weights .begin();
	for ( std::vector < Cell > ::iterator it_point = this->points .begin();
	      it_point != this->points .end(); it_point++                      )
	{	assert ( it_weight != this->weights .end() );
		Cell Gauss_point = *it_point;
		double w = *it_weight;
		res += w * f ( Gauss_point ) * tran->det ( Gauss_point );
		it_weight++;                                              }
	assert ( it_weight == this->weights .end() );
	return res;                                                                        }

//-----------------------------------------------------------------------------------------//

double Integrator::UFL_FFC::action ( Function f, const FiniteElement & )
// virtual from Integrator::Core

// assumes the finite element is already docked on a cell

{	assert ( false );  return 0.;   }

//-----------------------------------------------------------------------------------------//


void FiniteElement::WithMaster::Segment::dock_on ( const Cell & cll )
// virtual from FiniteElement::Core

{	assert ( cll .dim() == 1 );
	this->docked_on = cll;
	Cell P = cll .base() .reverse();
	Cell Q = cll .tip();

	Function t = this->master_manif  .coordinates();
	assert ( t .nb_of_components() == 1 );

	Function xyz = Manifold::working .coordinates();
	size_t geom_dim = xyz .nb_of_components();
	if ( geom_dim == 1 )
	{	double xP = xyz (P), xQ = xyz (Q);
		Function x_c = ( xP * (1.-t) + xQ * (1.+t) ) / 2.;
		this->transf = Function ( tag::diffeomorphism, tag::one_dim, xyz, t, x_c ); }
	else  // geometric dimension >= 2
	{	Function x = xyz [0];
		double xP = x(P), xQ = x(Q);
		Function xyz_c = ( xP * (1.-t) + xQ * (1.+t) ) / 2.;
		for ( size_t d = 1; d < geom_dim; d++ )
		{	x = xyz [d];
		  xP = x(P), xQ = x(Q);
			xyz_c = xyz_c && ( ( xP * (1.-t) + xQ * (1.+t) ) / 2. );  }
		assert ( xyz_c .nb_of_components() == geom_dim );
		this->transf = Function ( tag::immersion, xyz, t, xyz_c );     }

	this->base_fun_1 .clear();
	this->base_fun_1 .insert ( std::pair < Cell::Core*, Function >
		( P.core, Function ( (1.-t)/2., tag::composed_with, this->transf ) ) );
	this->base_fun_1 .insert ( std::pair < Cell::Core*, Function >
		( Q.core, Function ( (1.+t)/2., tag::composed_with, this->transf ) ) );
	
}  // end of  FiniteElement::WithMaster::Segment::dock_on

//-----------------------------------------------------------------------------------------//

void FiniteElement::WithMaster::Segment::dock_on ( const Cell & cll, const tag::Winding & )
// virtual from FiniteElement::Core

{	assert ( cll. dim() == 1 );
	this->docked_on = cll;
	Cell P = cll .base() .reverse();
	Cell Q = cll .tip();

	Function t = this->master_manif .coordinates();
	assert ( t .nb_of_components() == 1 );

	Function xyz = Manifold::working .coordinates();
	size_t geom_dim = xyz .nb_of_components();
	if ( geom_dim == 1 )
	{	double xP = xyz (P), xQ = xyz (Q);
		Function x_c = ( xP * (1.-t) + xQ * (1.+t) ) / 2.;
		this->transf = Function ( tag::diffeomorphism, tag::one_dim, xyz, t, x_c ); }
	else  // geometric dimension >= 2
	{	Function x = xyz [0];
		double xP = x (P), xQ = x (Q);
		Function xyz_c = ( xP * (1.-t) + xQ * (1.+t) ) / 2.;
		for ( size_t d = 1; d < geom_dim; d++ )
		{	x = xyz [d];
		  xP = x (P), xQ = x (Q);
			xyz_c = xyz_c && ( ( xP * (1.-t) + xQ * (1.+t) ) / 2. );  }
		assert ( xyz_c .nb_of_components() == geom_dim );
		this->transf = Function ( tag::immersion, xyz, t, xyz_c );     }

	this->base_fun_1 .clear();
	this->base_fun_1 .insert ( std::pair < Cell::Core*, Function >
		( P.core, Function ( (1.-t)/2., tag::composed_with, this->transf ) ) );
	this->base_fun_1 .insert ( std::pair < Cell::Core*, Function >
		( Q.core, Function ( (1.+t)/2., tag::composed_with, this->transf ) ) );
	
}  // end of  FiniteElement::WithMaster::Segment::dock_on  with tag::winding

//-----------------------------------------------------------------------------------------//


void FiniteElement::WithMaster::Triangle::dock_on ( const Cell & cll )
// virtual from FiniteElement::Core

{	assert ( cll .dim() == 2 );
	this->docked_on = cll;
	Mesh::Iterator it = cll .boundary() .iterator ( tag::over_vertices, tag::require_order );
	it .reset();  assert ( it .in_range() );  Cell P = *it;
	it++;  assert ( it .in_range() );  Cell Q = *it;
	it++;  assert ( it .in_range() );  Cell R = *it;
	it++;  assert ( not it .in_range() );
	
	Function xi_eta = this->master_manif .coordinates();
	assert ( xi_eta .nb_of_components() == 2 );
	Function xi = xi_eta [0], eta = xi_eta [1];
	Function one_m_xi_m_eta = 1. - xi - eta;

	// xi, eta and 1-xi-eta are a base of functions defined on the master element
	// we may need to differentiate them with respect to x and y (physical coordinates)
	// which is equivalent to differentiating x_c and y_c with respect to xi and eta
	// and then taking the inverse matrix

	Function xyz = Manifold::working .coordinates();
	size_t geom_dim = xyz .nb_of_components();
	assert ( geom_dim >= 2 );

	if ( geom_dim == 2 )

	{	Function x = xyz [0],  y = xyz [1];
		double xP = x (P), xQ = x (Q), xR = x (R);
		double yP = y (P), yQ = y (Q), yR = y (R);

		Function x_c = xP * one_m_xi_m_eta + xQ * xi + xR * eta;
		Function y_c = yP * one_m_xi_m_eta + yQ * xi + yR * eta;

		this->transf =
			Function ( tag::diffeomorphism, tag::high_dim, xyz, xi_eta, x_c && y_c );  }
	
	else  // geometric dimension >= 3
		
	{	Function x = xyz [0];
		double xP = x(P), xQ = x(Q), xR = x(R);
		Function xyz_c = xP * one_m_xi_m_eta + xQ * xi + xR * eta;

		for ( size_t d = 1; d < geom_dim; d++ )
		{	x = xyz [d];
		  xP = x (P), xQ = x (Q), xR = x (R);
			xyz_c = xyz_c && ( xP * one_m_xi_m_eta + xQ * xi + xR * eta );  }
		assert ( xyz_c .nb_of_components() == geom_dim );

		this->transf = Function ( tag::immersion, xyz, xi_eta, xyz_c );            }

	this->base_fun_1 .clear();
	this->base_fun_1 .insert ( std::pair < Cell::Core*, Function >
	       ( P.core, Function ( one_m_xi_m_eta, tag::composed_with, this->transf ) ) );
	this->base_fun_1 .insert ( std::pair < Cell::Core*, Function >
	       ( Q.core, Function ( xi, tag::composed_with, this->transf ) ) );
	this->base_fun_1 .insert ( std::pair < Cell::Core*, Function >
	       ( R.core, Function ( eta, tag::composed_with, this->transf ) ) );

	// the only use of composing is to allow the calling code to
	// differentiate base functions with respect to geometric coordinates
	// and anyway this only works for a diffeomorpsm, not for an immersion
	
}  // end of  FiniteElement::WithMaster::Triangle::dock_on

//-----------------------------------------------------------------------------------------//


void FiniteElement::WithMaster::Triangle::dock_on ( const Cell & cll, const tag::Winding & )
// virtual from FiniteElement::Core

{	assert ( cll .dim() == 2 );
	this->docked_on = cll;
	// perhaps implement a special iterator returning points and segments
	Mesh::Iterator it = cll .boundary() .iterator ( tag::over_vertices, tag::require_order );
	it .reset();  assert ( it .in_range() );  Cell P = *it;
	Cell PQ = cll .boundary() .cell_in_front_of ( P );
	it++;  assert ( it .in_range() );  Cell Q = *it;
	Cell QR = cll .boundary() .cell_in_front_of ( Q );
	it++;  assert ( it .in_range() );  Cell R = *it;
	#ifndef NDEBUG
	Cell RP = cll .boundary() .cell_in_front_of ( R );
	#endif
	it++;  assert ( not it .in_range() );

	Function::Action winding_Q = PQ .winding(), winding_R = winding_Q + QR .winding();
	assert ( winding_R + RP .winding() == 0 );
	
	Function xi_eta = this->master_manif .coordinates();
	assert ( xi_eta .nb_of_components() == 2 );
	Function xi = xi_eta [0], eta = xi_eta [1];
	Function one_m_xi_m_eta = 1.- xi - eta;

	// xi, eta and 1-xi-eta are a base of functions defined on the master element
	// we may need to differentiate them with respect to x and y (physical coordinates)
	// which is equivalent to differentiating x_c and y_c with respect to xi and eta
	// and then taking the inverse matrix

	Function xyz = Manifold::working .coordinates();
	size_t geom_dim = xyz .nb_of_components();
	assert ( geom_dim >= 2 );

	if ( geom_dim == 2 )

	{	std::vector < double > xyz_P = xyz ( P ),
		                       xyz_Q = xyz ( Q, tag::winding, winding_Q ),
		                       xyz_R = xyz ( R, tag::winding, winding_R );

		Function x_c = xyz_P [0] * one_m_xi_m_eta + xyz_Q [0] * xi + xyz_R [0] * eta;
		Function y_c = xyz_P [1] * one_m_xi_m_eta + xyz_Q [1] * xi + xyz_R [1] * eta;

		this->transf =
			Function ( tag::diffeomorphism, tag::high_dim, xyz, xi_eta, x_c && y_c );   }
	
	else  // geometric dimension >= 3
		
	{	Function x = xyz [0];
		double xP = x (P), xQ = x (Q), xR = x (R);
		Function xyz_c = xP * one_m_xi_m_eta + xQ * xi + xR * eta;

		for ( size_t d = 1; d < geom_dim; d++ )
		{	x = xyz [d];
		  xP = x (P), xQ = x (Q), xR = x (R);
			xyz_c = xyz_c && ( xP * one_m_xi_m_eta + xQ * xi + xR * eta );  }
		assert ( xyz_c .nb_of_components() == geom_dim );

		this->transf = Function ( tag::immersion, xyz, xi_eta, xyz_c );            }

	this->base_fun_1 .clear();
	this->base_fun_1 .insert ( std::pair < Cell::Core*, Function >
	       ( P.core, Function ( one_m_xi_m_eta, tag::composed_with, this->transf ) ) );
	this->base_fun_1 .insert ( std::pair < Cell::Core*, Function >
	       ( Q.core, Function ( xi, tag::composed_with, this->transf ) ) );
	this->base_fun_1 .insert ( std::pair < Cell::Core*, Function >
	       ( R.core, Function ( eta, tag::composed_with, this->transf ) ) );

	// the only use of composing is to allow the calling code to
	// differentiate base functions with respect to geometric coordinates
	// and anyway this only works for a diffeomorpsm, not for an immersion
	
}  // end of  FiniteElement::WithMaster::Triangle::dock_on  with tag::winding

//-----------------------------------------------------------------------------------------//


void FiniteElement::WithMaster::Quadrangle::dock_on ( const Cell & cll )
// virtual from FiniteElement::Core

{	assert ( cll .dim() == 2 );
	this->docked_on = cll;
	Mesh::Iterator it = cll .boundary() .iterator ( tag::over_vertices, tag::require_order );
	it .reset();  assert ( it .in_range() );  Cell P = *it;
	it++;  assert ( it .in_range() );  Cell Q = *it;
	it++;  assert ( it .in_range() );  Cell R = *it;
	it++;  assert ( it .in_range() );  Cell S = *it;
	it++;  assert ( not it .in_range() );
	
	Function xi_eta = this->master_manif .coordinates();
	assert ( xi_eta .nb_of_components() == 2 );
	Function xi = xi_eta [0], eta = xi_eta [1];

	Function psiP = (1.-xi) * (1.-eta), psiQ = (1.+xi) * (1.-eta),
		psiR = (1.+xi) * (1.+eta), psiS = (1.-xi) * (1.+eta);
	// psi* are a base of functions defined on the master element
	// each takes value 4 in the associated vertex, zero on the others
	// we may need to differentiate them with respect to x and y (physical coordinates)
	// which is equivalent to differentiating x_c and y_c with respect to xi and eta
	// and then taking the inverse matrix
	
	Function xyz = Manifold::working .coordinates();
	size_t geom_dim = xyz .nb_of_components();
	assert ( geom_dim >= 2 );

	if ( geom_dim == 2 )

	{	Function x = xyz [0],  y = xyz [1];
		double xP = x (P), xQ = x (Q), xR = x (R), xS = x (S);
		double yP = y (P), yQ = y (Q), yR = y (R), yS = y (S);

		Function x_c = ( xP * psiP + xQ * psiQ + xR * psiR + xS * psiS ) / 4.;
		Function y_c = ( yP * psiP + yQ * psiQ + yR * psiR + yS * psiS ) / 4.;

		this->transf =
			Function ( tag::diffeomorphism, tag::high_dim, xyz, xi_eta, x_c && y_c );  }

	else  // geometric dimension >= 3
		
	{	Function x = xyz [0];
		double xP = x(P), xQ = x(Q), xR = x(R), xS = x(S);
		Function xyz_c = ( xP * psiP + xQ * psiQ + xR * psiR + xS * psiS ) / 4.;

		for ( size_t d = 1; d < geom_dim; d++ )
		{	x = xyz [d];
		  xP = x(P), xQ = x(Q), xR = x(R);
			xyz_c = xyz_c && ( ( xP * psiP + xQ * psiQ + xR * psiR + xS * psiS ) / 4. );  }
		assert ( xyz_c .nb_of_components() == geom_dim );

		this->transf = Function ( tag::immersion, xyz, xi_eta, xyz_c );                    }

	this->base_fun_1 .clear();
	this->base_fun_1 .insert ( std::pair < Cell::Core*, Function >
	       ( P.core, Function ( psiP/4., tag::composed_with, this->transf ) ) );
	this->base_fun_1 .insert ( std::pair < Cell::Core*, Function >
	       ( Q.core, Function ( psiQ/4., tag::composed_with, this->transf ) ) );
	this->base_fun_1 .insert ( std::pair < Cell::Core*, Function >
	       ( R.core, Function ( psiR/4., tag::composed_with, this->transf ) ) );
	this->base_fun_1 .insert ( std::pair < Cell::Core*, Function >
	       ( S.core, Function ( psiS/4., tag::composed_with, this->transf ) ) );

	// the only use of composing is to allow the calling code to
	// differentiate base functions with respect to geometric coordinates
	// and anyway this only works for a diffeomorphism, not for an immersion
	
}  // end of  FiniteElement::WithMaster::Quadrangle::dock_on

//-----------------------------------------------------------------------------------------//


void FiniteElement::WithMaster::Quadrangle::dock_on ( const Cell & cll, const tag::Winding & )
// virtual from FiniteElement::Core

{	assert ( cll .dim() == 2 );
	this->docked_on = cll;
	Mesh::Iterator it = cll .boundary() .iterator ( tag::over_vertices, tag::require_order );
	it. reset();  assert ( it .in_range() );  Cell P = *it;
	Cell PQ = cll .boundary() .cell_in_front_of ( P );
	it++;  assert ( it .in_range() );  Cell Q = *it;
	Cell QR = cll .boundary() .cell_in_front_of ( Q );
	it++;  assert ( it .in_range() );  Cell R = *it;
	Cell RS = cll .boundary() .cell_in_front_of ( R );
	it++;  assert ( it .in_range() );  Cell S = *it;
	#ifndef NDEBUG
	Cell SP = cll .boundary() .cell_in_front_of ( S );
	#endif
	it++;  assert ( not it .in_range() );

	Function::Action winding_Q = PQ .winding(),
	                 winding_R = winding_Q + QR .winding(),
	                 winding_S = winding_R + RS .winding();
	assert ( winding_S + SP .winding() == 0 );
	
	Function xi_eta = this->master_manif .coordinates();
	assert ( xi_eta .nb_of_components() == 2 );
	Function xi = xi_eta [0], eta = xi_eta [1];

	Function psi_P = (1.-xi) * (1.-eta), psi_Q = (1.+xi) * (1.-eta),
	         psi_R = (1.+xi) * (1.+eta), psi_S = (1.-xi) * (1.+eta);
	// psi_* are a base of functions defined on the master element
	// each takes value 4 in the associated vertex, zero on the others
	// we may need to differentiate them with respect to x and y (physical coordinates)
	// which is equivalent to differentiating x_c and y_c with respect to xi and eta
	// and then taking the inverse matrix
	
	Function xyz = Manifold::working .coordinates();
	size_t geom_dim = xyz .nb_of_components();
	assert ( geom_dim >= 2 );

	if ( geom_dim == 2 )

	{	std::vector < double > xyz_P = xyz ( P ),
		                       xyz_Q = xyz ( Q, tag::winding, winding_Q ),
		                       xyz_R = xyz ( R, tag::winding, winding_R ),
		                       xyz_S = xyz ( S, tag::winding, winding_S );

		Function x_c = ( xyz_P[0] * psi_P + xyz_Q[0] * psi_Q +
	                   xyz_R[0] * psi_R + xyz_S[0] * psi_S ) / 4.;
		Function y_c = ( xyz_P[1] * psi_P + xyz_Q[1] * psi_Q +
	                   xyz_R[1] * psi_R + xyz_S[1] * psi_S ) / 4.;

	 this->transf =
		 Function ( tag::diffeomorphism, tag::high_dim, xyz, xi_eta, x_c && y_c );   }

	else  // geometric dimension >= 3
		
	{	Function x = xyz [0];
		double xP = x (P), xQ = x (Q), xR = x (R), xS = x (S);
		Function xyz_c = ( xP * psi_P + xQ * psi_Q + xR * psi_R + xS * psi_S ) / 4.;

		for ( size_t d = 1; d < geom_dim; d++ )
		{	x = xyz [d];
		  xP = x(P), xQ = x(Q), xR = x(R);
			xyz_c = xyz_c && ( ( xP * psi_P + xQ * psi_Q + xR * psi_R + xS * psi_S ) / 4. );  }
		assert ( xyz_c .nb_of_components() == geom_dim );

		this->transf = Function ( tag::immersion, xyz, xi_eta, xyz_c );                    }

	this->base_fun_1 .clear();
	this->base_fun_1 .insert ( std::pair < Cell::Core*, Function >
	       ( P.core, Function ( psi_P/4., tag::composed_with, this->transf ) ) );
	this->base_fun_1 .insert ( std::pair < Cell::Core*, Function >
	       ( Q.core, Function ( psi_Q/4., tag::composed_with, this->transf ) ) );
	this->base_fun_1 .insert ( std::pair < Cell::Core*, Function >
	       ( R.core, Function ( psi_R/4., tag::composed_with, this->transf ) ) );
	this->base_fun_1 .insert ( std::pair < Cell::Core*, Function >
	       ( S.core, Function ( psi_S/4., tag::composed_with, this->transf ) ) );

	// the only use of composing is to allow the calling code to
	// differentiate base functions with respect to geometric coordinates
	// and anyway this only works for a diffeomorphism, not for an immersion

}  // end of  FiniteElement::WithMaster::Quadrangle::dock_on  with tag::winding
	
//-----------------------------------------------------------------------------------------//


namespace { // anonymous namespace, mimics static linkage

inline void dock_on_common ( const double & xP, const double & yP,
                             const double & xQ, const double & yQ,
                             const double & xR, const double & yR,
														 const size_t & c,  // case
		  std::vector < std::vector < std::vector < double > > > & result )

// this function is speed-critical
	
{	double J_c0 = xQ - xP, J_c1 = xR - xP,
	       J_c2 = yQ - yP, J_c3 = yR - yP;

	// below we use a switch statement
	// we trust that the compiler implements it by means of a list of addresses
	// and not as a cascade of ifs (which would be rather slow)
	// we could use instead a vector of pointers to functions
			
	switch ( c )
		
	{	case  0 :
			std::cout << "UFL-FFC integrators require pre_compute" << std::endl;
			exit ( 1 );

		case  1 :  // { int psi }
			
		{	const double det = J_c0 * J_c3 - J_c1 * J_c2;
			assert ( det > 0. );  // det = 2. * area

			result [0][0][0] =                    // int psi^P
			result [0][1][0] =                    // int psi^Q
			result [0][2][0] = det / 6.;   }      // int psi^R
			break;  // end of case 1

		case  2 :  // { int psi1 * psi2 }
			
		{	const double det = J_c0 * J_c3 - J_c1 * J_c2;
			assert ( det > 0. );  // det = 2. * area

			result [0][0][0] =                 // int psi^P * psi^P
			result [1][1][0] =                 // int psi^Q * psi^Q
			result [2][2][0] =                 // int psi^R * psi^R
				0.08333333333333333 * det;
			result [0][1][0] =                 // int psi^P * psi^Q
			result [0][2][0] =                 // int psi^P * psi^R
			result [1][0][0] =                 // int psi^Q * psi^P
			result [1][2][0] =                 // int psi^Q * psi^R
			result [2][0][0] =                 // int psi^R * psi^P
			result [2][1][0] =                 // int psi^R * psi^Q
				0.04166666666666666 * det;                    }
			break;  // end of case 2

		case  3 :  // { int psi .deriv(x), int psi .deriv(y) }
		{
			#ifndef NDEBUG
			const double det = J_c0 * J_c3 - J_c1 * J_c2;
			assert ( det > 0. );  // det = 2. * area
			#endif

			J_c0 *= 0.5;  J_c1 *= 0.5;  J_c2 *= 0.5;  J_c3 *= 0.5;
	
			result [0][0][0] =   J_c2 - J_c3;          // int psi^P,x
			result [0][0][1] =   J_c1 - J_c0;          // int psi^P,y
			result [0][1][0] =   J_c3;                 // int psi^Q,x
			result [0][1][1] = - J_c1;                 // int psi^Q,y
			result [0][2][0] = - J_c2;                 // int psi^R,x
			result [0][2][1] =   J_c0;             }   // int psi^R,y
			break;  // end of case 3

		case  4 :  // { int psi1 .deriv(x) * psi2 .deriv(y) }
			
		{	double det = J_c0 * J_c3 - J_c1 * J_c2;
			assert ( det > 0. );  // det = 2. * area
			det += det;  // we double det to avoid later divisions by two
			const double p00 = J_c0 * J_c0 / det,
				           p01 = J_c0 * J_c1 / det,
				           p02 = J_c0 * J_c2 / det,
				           p03 = J_c0 * J_c3 / det,
				           p11 = J_c1 * J_c1 / det,
				           p12 = J_c1 * J_c2 / det,
				           p13 = J_c1 * J_c3 / det,
				           p22 = J_c2 * J_c2 / det,
				           p23 = J_c2 * J_c3 / det,
				           p33 = J_c3 * J_c3 / det;
			result [0][0][0] =   p33 - p23 - p23 + p22;    // int psi^P,x * psi^P,x
			result [0][1][0] =                             // int psi^P,x * psi^Q,x
			result [1][0][0] =   p23 - p33;                // int psi^Q,x * psi^P,x
			result [1][1][0] =   p33;                      // int psi^Q,x * psi^Q,x
			result [0][2][0] =                             // int psi^P,x * psi^R,x
			result [2][0][0] =   p23 - p22;                // int psi^R,x * psi^P,x
			result [1][2][0] =                             // int psi^Q,x * psi^R,x
			result [2][1][0] = - p23;                      // int psi^R,x * psi^Q,x
			result [2][2][0] =   p22;                      // int psi^R,x * psi^R,x
			result [0][0][1] =                             // int psi^P,x * psi^P,y
			result [0][0][2] =   p12 + p03 - p13 - p02;    // int psi^P,y * psi^P,x
			result [0][1][1] =                             // int psi^P,x * psi^Q,y
			result [1][0][2] =   p13 - p12;                // int psi^Q,y * psi^P,x
			result [0][2][1] =                             // int psi^P,x * psi^R,y
			result [2][0][2] =   p02 - p03;                // int psi^R,y * psi^P,x
			result [1][0][1] =                             // int psi^Q,x * psi^P,y
			result [0][1][2] =   p13 - p03;                // int psi^P,y * psi^Q,x
			result [1][1][1] =                             // int psi^Q,x * psi^Q,y
			result [1][1][2] = - p13;                      // int psi^Q,y * psi^Q,x
			result [1][2][1] =                             // int psi^Q,x * psi^R,y
			result [2][1][2] =   p03;                      // int psi^R,y * psi^Q,x
			result [2][0][1] =                             // int psi^R,x * psi^P,y
			result [0][2][2] =   p02 - p12;                // int psi^P,y * psi^R,x
			result [2][1][1] =                             // int psi^R,x * psi^Q,y
			result [1][2][2] =   p12;                      // int psi^Q,y * psi^R,x
			result [2][2][1] =                             // int psi^R,x * psi^R,y
			result [2][2][2] = - p02;                      // int psi^R,y * psi^R,x
			result [0][0][3] =   p11 - p01 - p01 + p00;    // int psi^P,y * psi^P,y
			result [0][1][3] =                             // int psi^P,y * psi^Q,y
			result [1][0][3] =   p01 - p11;                // int psi^Q,y * psi^P,y
			result [1][1][3] =   p11;                      // int psi^Q,y * psi^Q,y
			result [0][2][3] =                             // int psi^P,y * psi^R,y
			result [2][0][3] =   p01 - p00;                // int psi^R,y * psi^P,y
			result [1][2][3] =                             // int psi^Q,y * psi^R,y
			result [2][1][3] = - p01;                      // int psi^R,y * psi^Q,y
			result [2][2][3] =   p00;                   }  // int psi^R,y * psi^R,y
			break;  // end of case 4

		case  5 :  // { int grad psi grad psi }
		// { int psi1 .deriv(x) * psi2 .dervi(x) + psi1 .deriv(y) * psi2 .deriv(y) }

		{	double det = J_c0 * J_c3 - J_c1 * J_c2;
			assert ( det > 0. );  // det = 2. * area
			det += det;  // we double det to avoid later divisions by two
			const double sp17 =   ( J_c0 * J_c0 + J_c2 * J_c2 ) / det,
			             sp18 = - ( J_c0 * J_c1 + J_c3 * J_c2 ) / det,
			             sp19 =   ( J_c1 * J_c1 + J_c3 * J_c3 ) / det;
		
			result [0][0][0] =   sp19 + sp18 + sp18 + sp17;      // int grad psi^P grad psi^P
			result [1][0][0] =                                   // int grad psi^Q grad psi^P
			result [0][1][0] = - sp19 - sp18;                    // int grad psi^P grad psi^Q
			result [1][1][0] =   sp19;                           // int grad psi^Q grad psi^Q
			result [2][0][0] =                                   // int grad psi^R grad psi^P
			result [0][2][0] = - sp18 - sp17;                    // int grad psi^P grad psi^R
			result [2][1][0] =                                   // int grad psi^R grad psi^Q
			result [1][2][0] =   sp18;                           // int grad psi^Q grad psi^R
			result [2][2][0] =   sp17;                       }   // int grad psi^R grad psi^R
			break;  // end of case 5

		case  6 :  // { int psi2 * psi2 .deriv(x) }
		{
			#ifndef NDEBUG
			const double det = J_c0 * J_c3 - J_c1 * J_c2;
			assert ( det > 0. );  // det = 2. * area
			#endif

			J_c0 /= 6.;  J_c1 /= 6.;  J_c2 /= 6.;  J_c3 /= 6;
	
			result [0][0][0] =                      // int psi^P psi^P,x
			result [0][0][1] =                      // int psi^P psi^P,x
			result [1][0][0] =                      // int psi^Q psi^P,x
			result [0][1][1] =                      // int psi^Q psi^P,x
			result [2][0][0] =                      // int psi^R psi^P,x
			result [0][2][1] =   J_c2 - J_c3;       // int psi^R psi^P,x
			result [0][1][0] =                      // int psi^P psi^Q,x
			result [1][0][1] =                      // int psi^P psi^Q,x
			result [1][1][0] =                      // int psi^Q psi^Q,x
			result [1][1][1] =                      // int psi^Q psi^Q,x
			result [2][1][0] =                      // int psi^R psi^Q,x
			result [1][2][1] =   J_c3;              // int psi^R psi^Q,x
			result [0][2][0] =                      // int psi^P psi^R,x
			result [2][0][1] =                      // int psi^P psi^R,x
			result [1][2][0] =                      // int psi^Q psi^R,x
			result [2][1][1] =                      // int psi^Q psi^R,x
			result [2][2][0] =                      // int psi^R psi^R,x
			result [2][2][1] = - J_c2;              // int psi^R psi^R,x
			result [0][0][2] =                      // int psi^P psi^P,y
			result [0][0][3] =                      // int psi^P psi^P,y
			result [1][0][2] =                      // int psi^Q psi^P,y
			result [0][1][3] =                      // int psi^Q psi^P,y
			result [2][0][2] =                      // int psi^R psi^P,y
			result [0][2][3] =   J_c1 - J_c0;       // int psi^R psi^P,y
			result [0][1][2] =                      // int psi^P psi^Q,y
			result [1][0][3] =                      // int psi^P psi^Q,y
			result [1][1][2] =                      // int psi^Q psi^Q,y
			result [1][1][3] =                      // int psi^Q psi^Q,y
			result [2][1][2] =                      // int psi^R psi^Q,y
			result [1][2][3] = - J_c1;              // int psi^R psi^Q,y
			result [0][2][2] =                      // int psi^P psi^R,y
			result [2][0][3] =                      // int psi^P psi^R,y
			result [1][2][2] =                      // int psi^Q psi^R,y
			result [2][1][3] =                      // int psi^Q psi^R,y
			result [2][2][2] =                      // int psi^R psi^R,y
			result [2][2][3] =   J_c0;          }   // int psi^R psi^R,y
			break;  // end of case 6

		case  7 :  // { int psi, int psi1 * psi2 }

		{	double det = J_c0 * J_c3 - J_c1 * J_c2;
			assert ( det > 0. );  // det = 2. * area
			result [0][0][0] =                 // int psi^P
			result [1][0][0] =                 // int psi^P
			result [2][0][0] =                 // int psi^P
			result [0][0][1] =                 // int psi^P
			result [0][1][1] =                 // int psi^P
			result [0][2][1] =                 // int psi^P
			result [0][1][0] =                 // int psi^Q
			result [1][1][0] =                 // int psi^Q
			result [2][1][0] =                 // int psi^Q
			result [1][0][1] =                 // int psi^Q
			result [1][1][1] =                 // int psi^Q
			result [1][2][1] =                 // int psi^Q
			result [0][2][0] =                 // int psi^R
			result [1][2][0] =                 // int psi^R
			result [2][2][0] =                 // int psi^R
			result [2][0][1] =                 // int psi^R
			result [2][1][1] =                 // int psi^R
			result [2][2][1] = det / 6.;       // int psi^R
			result [0][0][2] =                 // int psi^P * psi^P
			result [1][1][2] =                 // int psi^Q * psi^Q
			result [2][2][2] =                 // int psi^R * psi^R
				0.08333333333333333 * det;
			result [0][1][2] =                 // int psi^P * psi^Q
			result [0][2][2] =                 // int psi^P * psi^R
			result [1][0][2] =                 // int psi^Q * psi^P
			result [1][2][2] =                 // int psi^Q * psi^R
			result [2][0][2] =                 // int psi^R * psi^P
			result [2][1][2] =                 // int psi^R * psi^Q
				0.04166666666666666 * det;                    }
			break;  // end of case 7

		case  8 :  // { int psi, int psi .deriv(x) }
			
		{	const double det = J_c0 * J_c3 - J_c1 * J_c2;
			assert ( det > 0. );  // det = 2. * area
			
			J_c0 *= 0.5;  J_c1 *= 0.5;  J_c2 *= 0.5;  J_c3 *= 0.5;

			result [0][0][0] =                    // int psi^P
			result [0][1][0] =                    // int psi^Q
			result [0][2][0] =   det / 6.;        // int psi^R	
			result [0][0][1] =   J_c2 - J_c3;          // int psi^P,x
			result [0][0][2] =   J_c1 - J_c0;          // int psi^P,y
			result [0][1][1] =   J_c3;                 // int psi^Q,x
			result [0][1][2] = - J_c1;                 // int psi^Q,y
			result [0][2][1] = - J_c2;                 // int psi^R,x
			result [0][2][2] =   J_c0;             }   // int psi^R,y
			break;  // end of case 8

		case  9 :  // { int psi1 .deriv(x), int psi1 .deriv(x) * psi2 .deriv(y) }
			
		{	double det = J_c0 * J_c3 - J_c1 * J_c2;
			assert ( det > 0. );  // det = 2. * area
			det += det;  // we double det to avoid later divisions by two

			const double p00 = J_c0 * J_c0 / det,
				           p01 = J_c0 * J_c1 / det,
				           p02 = J_c0 * J_c2 / det,
				           p03 = J_c0 * J_c3 / det,
				           p11 = J_c1 * J_c1 / det,
				           p12 = J_c1 * J_c2 / det,
				           p13 = J_c1 * J_c3 / det,
				           p22 = J_c2 * J_c2 / det,
				           p23 = J_c2 * J_c3 / det,
				           p33 = J_c3 * J_c3 / det;

			J_c0 *= 0.5;  J_c1 *= 0.5;  J_c2 *= 0.5;  J_c3 *= 0.5;

			result [0][0][0] =                         // int psi^P,x
			result [0][1][0] =                         // int psi^P,x
			result [0][2][0] =                         // int psi^P,x
			result [0][0][1] =                         // int psi^P,x
			result [1][0][1] =                         // int psi^P,x
			result [2][0][1] =   J_c2 - J_c3;          // int psi^P,x
			result [0][0][2] =                         // int psi^P,y
			result [0][1][2] =                         // int psi^P,y
			result [0][2][2] =                         // int psi^P,y
			result [0][0][3] =                         // int psi^P,y
			result [1][0][3] =                         // int psi^P,y
			result [2][0][3] =   J_c1 - J_c0;          // int psi^P,y
			result [1][0][0] =                         // int psi^Q,x
			result [1][1][0] =                         // int psi^Q,x
			result [1][2][0] =                         // int psi^Q,x
			result [0][1][1] =                         // int psi^Q,x
			result [1][1][1] =                         // int psi^Q,x
			result [2][1][1] =   J_c3;                 // int psi^Q,x
			result [1][0][2] =                         // int psi^Q,y
			result [1][1][2] =                         // int psi^Q,y
			result [1][2][2] =                         // int psi^Q,y
			result [0][1][3] =                         // int psi^Q,y
			result [1][1][3] =                         // int psi^Q,y
			result [2][1][3] = - J_c1;                 // int psi^Q,y
			result [2][0][0] =                         // int psi^R,x
			result [2][1][0] =                         // int psi^R,x
			result [2][2][0] =                         // int psi^R,x
			result [0][2][1] =                         // int psi^R,x
			result [1][2][1] =                         // int psi^R,x
			result [2][2][1] = - J_c2;                 // int psi^R,x
			result [2][0][2] =                         // int psi^R,y
			result [2][1][2] =                         // int psi^R,y
			result [2][2][2] =                         // int psi^R,y
			result [0][2][3] =                         // int psi^R,y
			result [1][2][3] =                         // int psi^R,y
			result [2][2][3] =   J_c0;                 // int psi^R,y
			result [0][0][4] =   p33 - p23 - p23 + p22;    // int psi^P,x * psi^P,x
			result [0][1][4] =                             // int psi^P,x * psi^Q,x
			result [1][0][4] =   p23 - p33;                // int psi^Q,x * psi^P,x
			result [1][1][4] =   p33;                      // int psi^Q,x * psi^Q,x
			result [0][2][4] =                             // int psi^P,x * psi^R,x
			result [2][0][4] =   p23 - p22;                // int psi^R,x * psi^P,x
			result [1][2][4] =                             // int psi^Q,x * psi^R,x
			result [2][1][4] = - p23;                      // int psi^R,x * psi^Q,x
			result [2][2][4] =   p22;                      // int psi^R,x * psi^R,x
			result [0][0][5] =                             // int psi^P,x * psi^P,y
			result [0][0][6] =   p12 + p03 - p13 - p02;    // int psi^P,y * psi^P,x
			result [0][1][5] =                             // int psi^P,x * psi^Q,y
			result [1][0][6] =   p13 - p12;                // int psi^Q,y * psi^P,x
			result [0][2][5] =                             // int psi^P,x * psi^R,y
			result [2][0][6] =   p02 - p03;                // int psi^R,y * psi^P,x
			result [1][0][5] =                             // int psi^Q,x * psi^P,y
			result [0][1][6] =   p13 - p03;                // int psi^P,y * psi^Q,x
			result [1][1][5] =                             // int psi^Q,x * psi^Q,y
			result [1][1][6] = - p13;                      // int psi^Q,y * psi^Q,x
			result [1][2][5] =                             // int psi^Q,x * psi^R,y
			result [2][1][6] =   p03;                      // int psi^R,y * psi^Q,x
			result [2][0][5] =                             // int psi^R,x * psi^P,y
			result [0][2][6] =   p02 - p12;                // int psi^P,y * psi^R,x
			result [2][1][5] =                             // int psi^R,x * psi^Q,y
			result [1][2][6] =   p12;                      // int psi^Q,y * psi^R,x
			result [2][2][5] =                             // int psi^R,x * psi^R,y
			result [2][2][6] = - p02;                      // int psi^R,y * psi^R,x
			result [0][0][7] =   p11 - p01 - p01 + p00;    // int psi^P,y * psi^P,y
			result [0][1][7] =                             // int psi^P,y * psi^Q,y
			result [1][0][7] =   p01 - p11;                // int psi^Q,y * psi^P,y
			result [1][1][7] =   p11;                      // int psi^Q,y * psi^Q,y
			result [0][2][7] =                             // int psi^P,y * psi^R,y
			result [2][0][7] =   p01 - p00;                // int psi^R,y * psi^P,y
			result [1][2][7] =                             // int psi^Q,y * psi^R,y
			result [2][1][7] = - p01;                      // int psi^R,y * psi^Q,y
			result [2][2][7] =   p00;                   }  // int psi^R,y * psi^R,y
			break;  // end of case 9

		case  10 :  // { int psi1 * psi2, int psi1 .deriv(x) * psi2 .deriv(y) }

		{	double det = J_c0 * J_c3 - J_c1 * J_c2;
			assert ( det > 0. );  // det = 2. * area
			
			result [0][0][0] =                 // int psi^P * psi^P
			result [1][1][0] =                 // int psi^Q * psi^Q
			result [2][2][0] =                 // int psi^R * psi^R
				0.08333333333333333 * det;
			result [0][1][0] =                 // int psi^P * psi^Q
			result [0][2][0] =                 // int psi^P * psi^R
			result [1][0][0] =                 // int psi^Q * psi^P
			result [1][2][0] =                 // int psi^Q * psi^R
			result [2][0][0] =                 // int psi^R * psi^P
			result [2][1][0] =                 // int psi^R * psi^Q
				0.04166666666666666 * det;
			
			det += det;  // we double det to avoid later divisions by two
			const double p00 = J_c0 * J_c0 / det,
				           p01 = J_c0 * J_c1 / det,
				           p02 = J_c0 * J_c2 / det,
				           p03 = J_c0 * J_c3 / det,
				           p11 = J_c1 * J_c1 / det,
				           p12 = J_c1 * J_c2 / det,
				           p13 = J_c1 * J_c3 / det,
				           p22 = J_c2 * J_c2 / det,
				           p23 = J_c2 * J_c3 / det,
				           p33 = J_c3 * J_c3 / det;
			result [0][0][1] =   p33 - p23 - p23 + p22;    // int psi^P,x * psi^P,x
			result [0][1][1] =                             // int psi^P,x * psi^Q,x
			result [1][0][1] =   p23 - p33;                // int psi^Q,x * psi^P,x
			result [1][1][1] =   p33;                      // int psi^Q,x * psi^Q,x
			result [0][2][1] =                             // int psi^P,x * psi^R,x
			result [2][0][1] =   p23 - p22;                // int psi^R,x * psi^P,x
			result [1][2][1] =                             // int psi^Q,x * psi^R,x
			result [2][1][1] = - p23;                      // int psi^R,x * psi^Q,x
			result [2][2][1] =   p22;                      // int psi^R,x * psi^R,x
			result [0][0][2] =                             // int psi^P,x * psi^P,y
			result [0][0][3] =   p12 + p03 - p13 - p02;    // int psi^P,y * psi^P,x
			result [0][1][2] =                             // int psi^P,x * psi^Q,y
			result [1][0][3] =   p13 - p12;                // int psi^Q,y * psi^P,x
			result [0][2][2] =                             // int psi^P,x * psi^R,y
			result [2][0][3] =   p02 - p03;                // int psi^R,y * psi^P,x
			result [1][0][2] =                             // int psi^Q,x * psi^P,y
			result [0][1][3] =   p13 - p03;                // int psi^P,y * psi^Q,x
			result [1][1][2] =                             // int psi^Q,x * psi^Q,y
			result [1][1][3] = - p13;                      // int psi^Q,y * psi^Q,x
			result [1][2][2] =                             // int psi^Q,x * psi^R,y
			result [2][1][3] =   p03;                      // int psi^R,y * psi^Q,x
			result [2][0][2] =                             // int psi^R,x * psi^P,y
			result [0][2][3] =   p02 - p12;                // int psi^P,y * psi^R,x
			result [2][1][2] =                             // int psi^R,x * psi^Q,y
			result [1][2][3] =   p12;                      // int psi^Q,y * psi^R,x
			result [2][2][2] =                             // int psi^R,x * psi^R,y
			result [2][2][3] = - p02;                      // int psi^R,y * psi^R,x
			result [0][0][4] =   p11 - p01 - p01 + p00;    // int psi^P,y * psi^P,y
			result [0][1][4] =                             // int psi^P,y * psi^Q,y
			result [1][0][4] =   p01 - p11;                // int psi^Q,y * psi^P,y
			result [1][1][4] =   p11;                      // int psi^Q,y * psi^Q,y
			result [0][2][4] =                             // int psi^P,y * psi^R,y
			result [2][0][4] =   p01 - p00;                // int psi^R,y * psi^P,y
			result [1][2][4] =                             // int psi^Q,y * psi^R,y
			result [2][1][4] = - p01;                      // int psi^R,y * psi^Q,y
			result [2][2][4] =   p00;                   }  // int psi^R,y * psi^R,y
			break;  // end of case 10

		case  11 :  // { int psi1 * psi2, int grad psi1 * grad psi2 }

		{	double det = J_c0 * J_c3 - J_c1 * J_c2;
			assert ( det > 0. );  // det = 2. * area
			
			result [0][0][0] =                 // int psi^P * psi^P
			result [1][1][0] =                 // int psi^Q * psi^Q
			result [2][2][0] =                 // int psi^R * psi^R
				0.08333333333333333 * det;
			result [0][1][0] =                 // int psi^P * psi^Q
			result [0][2][0] =                 // int psi^P * psi^R
			result [1][0][0] =                 // int psi^Q * psi^P
			result [1][2][0] =                 // int psi^Q * psi^R
			result [2][0][0] =                 // int psi^R * psi^P
			result [2][1][0] =                 // int psi^R * psi^Q
				0.04166666666666666 * det;
			
			det += det;  // we double det to avoid later divisions by two
			const double sp17 =   ( J_c0 * J_c0 + J_c2 * J_c2 ) / det,
			             sp18 = - ( J_c0 * J_c1 + J_c3 * J_c2 ) / det,
			             sp19 =   ( J_c1 * J_c1 + J_c3 * J_c3 ) / det;
		
			result [0][0][1] =   sp19 + sp18 + sp18 + sp17;      // int grad psi^P grad psi^P
			result [1][0][1] =                                   // int grad psi^Q grad psi^P
			result [0][1][1] = - sp19 - sp18;                    // int grad psi^P grad psi^Q
			result [1][1][1] =   sp19;                           // int grad psi^Q grad psi^Q
			result [2][0][1] =                                   // int grad psi^R grad psi^P
			result [0][2][1] = - sp18 - sp17;                    // int grad psi^P grad psi^R
			result [2][1][1] =                                   // int grad psi^R grad psi^Q
			result [1][2][1] =   sp18;                           // int grad psi^Q grad psi^R
			result [2][2][1] =   sp17;                       }   // int grad psi^R grad psi^R
			break;  // end of case 11

		case  12 :  // 1, 2, 3 and 4

		std::cout << "**** ";
			
		{	//  case  1 :  { int psi }

			double det = J_c0 * J_c3 - J_c1 * J_c2;
			assert ( det > 0. );  // det = 2. * area

			result [0][0][0] =                 // int psi^P
			result [1][0][0] =                 // int psi^P
			result [2][0][0] =                 // int psi^P
			result [0][0][1] =                 // int psi^P
			result [0][1][1] =                 // int psi^P
			result [0][2][1] =                 // int psi^P
			result [0][1][0] =                 // int psi^Q
			result [1][1][0] =                 // int psi^Q
			result [2][1][0] =                 // int psi^Q
			result [1][0][1] =                 // int psi^Q
			result [1][1][1] =                 // int psi^Q
			result [1][2][1] =                 // int psi^Q
			result [0][2][0] =                 // int psi^R
			result [1][2][0] =                 // int psi^R
			result [2][2][0] =                 // int psi^R
			result [2][0][1] =                 // int psi^R
			result [2][1][1] =                 // int psi^R
			result [2][2][1] = det / 6.;       // int psi^R

			// case  2 :  { int psi1 * psi2 }

			result [0][0][2] =                 // int psi^P * psi^P
			result [1][1][2] =                 // int psi^Q * psi^Q
			result [2][2][2] =                 // int psi^R * psi^R
				0.08333333333333333 * det;
			result [0][1][2] =                 // int psi^P * psi^Q
			result [0][2][2] =                 // int psi^P * psi^R
			result [1][0][2] =                 // int psi^Q * psi^P
			result [1][2][2] =                 // int psi^Q * psi^R
			result [2][0][2] =                 // int psi^R * psi^P
			result [2][1][2] =                 // int psi^R * psi^Q
				0.04166666666666666 * det;

			// case  4 :  { int psi1 .deriv(x) * psi2 .deriv(y) }
			
			det += det;  // we double det to avoid later divisions by two
			const double p00 = J_c0 * J_c0 / det,
				           p01 = J_c0 * J_c1 / det,
				           p02 = J_c0 * J_c2 / det,
				           p03 = J_c0 * J_c3 / det,
				           p11 = J_c1 * J_c1 / det,
				           p12 = J_c1 * J_c2 / det,
				           p13 = J_c1 * J_c3 / det,
				           p22 = J_c2 * J_c2 / det,
				           p23 = J_c2 * J_c3 / det,
				           p33 = J_c3 * J_c3 / det;
			result [0][0][3] =   p33 - p23 - p23 + p22;    // int psi^P,x * psi^P,x
			result [0][1][3] =                             // int psi^P,x * psi^Q,x
			result [1][0][3] =   p23 - p33;                // int psi^Q,x * psi^P,x
			result [1][1][3] =   p33;                      // int psi^Q,x * psi^Q,x
			result [0][2][3] =                             // int psi^P,x * psi^R,x
			result [2][0][3] =   p23 - p22;                // int psi^R,x * psi^P,x
			result [1][2][3] =                             // int psi^Q,x * psi^R,x
			result [2][1][3] = - p23;                      // int psi^R,x * psi^Q,x
			result [2][2][3] =   p22;                      // int psi^R,x * psi^R,x
			result [0][0][4] =                             // int psi^P,x * psi^P,y
			result [0][0][5] =   p12 + p03 - p13 - p02;    // int psi^P,y * psi^P,x
			result [0][1][4] =                             // int psi^P,x * psi^Q,y
			result [1][0][5] =   p13 - p12;                // int psi^Q,y * psi^P,x
			result [0][2][4] =                             // int psi^P,x * psi^R,y
			result [2][0][5] =   p02 - p03;                // int psi^R,y * psi^P,x
			result [1][0][4] =                             // int psi^Q,x * psi^P,y
			result [0][1][5] =   p13 - p03;                // int psi^P,y * psi^Q,x
			result [1][1][4] =                             // int psi^Q,x * psi^Q,y
			result [1][1][5] = - p13;                      // int psi^Q,y * psi^Q,x
			result [1][2][4] =                             // int psi^Q,x * psi^R,y
			result [2][1][5] =   p03;                      // int psi^R,y * psi^Q,x
			result [2][0][4] =                             // int psi^R,x * psi^P,y
			result [0][2][5] =   p02 - p12;                // int psi^P,y * psi^R,x
			result [2][1][4] =                             // int psi^R,x * psi^Q,y
			result [1][2][5] =   p12;                      // int psi^Q,y * psi^R,x
			result [2][2][4] =                             // int psi^R,x * psi^R,y
			result [2][2][5] = - p02;                      // int psi^R,y * psi^R,x
			result [0][0][6] =   p11 - p01 - p01 + p00;    // int psi^P,y * psi^P,y
			result [0][1][6] =                             // int psi^P,y * psi^Q,y
			result [1][0][6] =   p01 - p11;                // int psi^Q,y * psi^P,y
			result [1][1][6] =   p11;                      // int psi^Q,y * psi^Q,y
			result [0][2][6] =                             // int psi^P,y * psi^R,y
			result [2][0][6] =   p01 - p00;                // int psi^R,y * psi^P,y
			result [1][2][6] =                             // int psi^Q,y * psi^R,y
			result [2][1][6] = - p01;                      // int psi^R,y * psi^Q,y
			result [2][2][6] =   p00;                      // int psi^R,y * psi^R,y
	
			// case  3 :  { int psi .deriv(x), int psi .deriv(y) }

			J_c0 *= 0.5;  J_c1 *= 0.5;  J_c2 *= 0.5;  J_c3 *= 0.5;
	
			result [0][0][7] =                         // int psi^P,x
			result [0][1][7] =                         // int psi^P,x
			result [0][2][7] =                         // int psi^P,x
			result [0][0][8] =                         // int psi^P,x
			result [1][0][8] =                         // int psi^P,x
			result [2][0][8] =   J_c2 - J_c3;          // int psi^P,x
			result [0][0][9] =                         // int psi^P,y
			result [0][1][9] =                         // int psi^P,y
			result [0][2][9] =                         // int psi^P,y
			result [0][0][10] =                         // int psi^P,y
			result [1][0][10] =                         // int psi^P,y
			result [2][0][10] =   J_c1 - J_c0;          // int psi^P,y
			result [1][0][7] =                         // int psi^Q,x
			result [1][1][7] =                         // int psi^Q,x
			result [1][2][7] =                         // int psi^Q,x
			result [0][1][8] =                         // int psi^Q,x
			result [1][1][8] =                         // int psi^Q,x
			result [2][1][8] =   J_c3;                 // int psi^Q,x
			result [1][0][9] =                         // int psi^Q,y
			result [1][1][9] =                         // int psi^Q,y
			result [1][2][9] =                         // int psi^Q,y
			result [0][1][10] =                         // int psi^Q,y
			result [1][1][10] =                         // int psi^Q,y
			result [2][1][10] = - J_c1;                 // int psi^Q,y
			result [2][0][7] =                         // int psi^R,x
			result [2][1][7] =                         // int psi^R,x
			result [2][2][7] =                         // int psi^R,x
			result [0][2][8] =                         // int psi^R,x
			result [1][2][8] =                         // int psi^R,x
			result [2][2][8] = - J_c2;                 // int psi^R,x
			result [2][0][9] =                         // int psi^R,y
			result [2][1][9] =                         // int psi^R,y
			result [2][2][9] =                         // int psi^R,y
			result [0][2][10] =                         // int psi^R,y
			result [1][2][10] =                         // int psi^R,y
			result [2][2][10] =   J_c0;              }  // int psi^R,y

			break;  // end of case 12

		default : assert ( false );
	}  // end of  switch  statement
}  // end of dock_on_common
	
}  // anonymous namespace


void FiniteElement::StandAlone::TypeOne::Triangle::dock_on ( const Cell & cll )
// virtual from FiniteElement::Core

{	assert ( cll .dim() == 2 );
	this->docked_on = cll;
	Mesh::Iterator it = cll .boundary() .iterator ( tag::over_vertices, tag::require_order );
	it .reset();  assert ( it .in_range() );  Cell P = *it;
	it++;  assert ( it .in_range() );  Cell Q = *it;
	it++;  assert ( it .in_range() );  Cell R = *it;
	it++;  assert ( not it .in_range() );
	
	Function xyz = Manifold::working .coordinates();
	size_t geom_dim = xyz .nb_of_components();
	assert ( geom_dim >= 2 );

	assert ( geom_dim == 2 );
	Function x = xyz [0], y = xyz [1];

	// using Cell::Core::short_int_heap for local numbering of vertices
	// should be slightly faster

	dock_on_common ( x(P), y(P), x(Q), y(Q), x(R), y(R),
	                 this->cas, this->result_of_integr  );

	// const double xP = x (P), xQ = x (Q), xR = x (R);
	// const double yP = y (P), yQ = y (Q), yR = y (R);

	// code below can be viewed as a local numbering of vertices P, Q, R
	this->base_fun_1 .clear();
	this->base_fun_1 .insert
		( std::pair < Cell::Core*, Function > ( P .core, this->bf1 ) );
	this->base_fun_1 .insert
		( std::pair < Cell::Core*, Function > ( Q .core, this->bf2 ) );
	this->base_fun_1 .insert
		( std::pair < Cell::Core*, Function > ( R .core, this->bf3 ) );

}  // end of  FiniteElement::StandAlone::TypeOne::Triangle::dock_on

//-----------------------------------------------------------------------------------------//


void FiniteElement::StandAlone::TypeOne::Triangle::dock_on  
( const Cell & cll, const tag::Winding & )  // virtual from FiniteElement::Core


{	
}  // end of  FiniteElement::StandAlone::TypeOne::Triangle::dock_on  with tag::winding

//-----------------------------------------------------------------------------------------//


void FiniteElement::WithMaster::pre_compute  // virtual from FiniteElement::Core
( const std::vector < Function > & bf, const std::vector < Function > & result )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "pre-compute does not apply to elements with master" << std::endl;
	exit ( 1 );                                                                     }


void FiniteElement::StandAlone::TypeOne::Triangle::pre_compute
( const std::vector < Function > & bf, const std::vector < Function > & result )
// virtual from FiniteElement::Core

// 'bf' contains one or two dummy basis functions, of type Function::MereSymbol
// 'result' contains a list of expressions whose integral the user wants
// the numeric values will be only computed later, at 'dock_on'
// and retrieved by 'integrate' which calls 'retrieve_precomputed'
	
{	this->result_of_integr .clear();
	this->selector .clear();

	this->dummy_bf = bf;
	
	// we analyse syntactically expressions listed in 'result'
	// and identify cases previously studied and hard-coded
	// if none matches, stop with error message

	// possible situations :
	//  1. integral of psi
	//  2. integral of psi * psi
	//  3. integral of d_psi_dx
	//  4. integral of d_psi_dx * d_psi_dy
	//  5. integral of grad psi * grad psi
	//  6. integral of psi * d_psi_dx
	//  7. 1 and 2
	//  8. 1 and 3
	//  9. 3 and 4
	// 10. 2 and 4
	// 11. 2 and 5
	// 12. 1, 2, 3 and 4
	// 13. 1, 2, 3, 4 and 5
	// 14. everything

	bool int_psi = false,       int_psi_psi = false,   int_dpsi = false,
	     int_dpsi_dpsi = false, int_grad_grad = false, int_psi_dpsi = false;

	for ( std::vector < Function > ::const_iterator it = result .begin();
	      it != result .end(); it++                                      )
	{	Function expr = *it;
		Function::MereSymbol * ms = dynamic_cast < Function::MereSymbol* > ( expr .core );
		if ( ms ) { int_psi = true;  continue;  }
		Function::DelayedDerivative * dd =
			dynamic_cast < Function::DelayedDerivative* > ( expr .core );
		if ( dd )
		{ if ( dynamic_cast < Function::MereSymbol* > ( dd->base .core ) == nullptr )
				goto error;
			int_dpsi = true;  continue;                                                 }
		Function::Product * prod = dynamic_cast < Function::Product* > ( expr .core );
		if ( prod )
		{	std::forward_list < Function > ::const_iterator itt = prod->factors .begin();
			assert ( itt != prod->factors .end() );
			ms = dynamic_cast < Function::MereSymbol* > ( itt->core );
			if ( ms )
			{	itt++;  assert ( itt != prod->factors .end() );
				ms = dynamic_cast < Function::MereSymbol* > ( itt->core );
				if ( ms )
				{	itt++;  if ( itt != prod->factors .end() ) goto error;
					int_psi_psi = true;  continue;                         }
				dd = dynamic_cast < Function::DelayedDerivative* > ( itt->core );
				if ( dd == nullptr ) goto error;
				if ( dynamic_cast < Function::MereSymbol* > ( dd->base .core ) == nullptr )
					goto error;
				itt++;  if ( itt != prod->factors .end() ) goto error;
				int_psi_dpsi = true;  continue;                                              }
			dd = dynamic_cast < Function::DelayedDerivative* > ( itt->core );
			if ( dd == nullptr ) goto error;
			if ( dynamic_cast < Function::MereSymbol* > ( dd->base .core ) == nullptr )
				goto error;
			itt++;  assert ( itt != prod->factors .end() );
			ms = dynamic_cast < Function::MereSymbol* > ( itt->core );
			if ( ms )
			{	itt++;  if ( itt != prod->factors .end() ) goto error;
				int_psi_dpsi = true;  continue;                       }
			dd = dynamic_cast < Function::DelayedDerivative* > ( itt->core );
			if ( dd == nullptr ) goto error;
			if ( dynamic_cast < Function::MereSymbol* > ( dd->base .core ) == nullptr )
				goto error;
			itt++;  if ( itt != prod->factors .end() ) goto error;
			int_dpsi_dpsi = true;  continue;                                                   }
		Function::Sum * s = dynamic_cast < Function::Sum* > ( expr .core );
		if ( s == nullptr ) goto error;
		for ( std::forward_list < Function > ::const_iterator itt = s->terms .begin();
					itt != s->terms .end(); itt++                                           )
			assert ( dynamic_cast < Function::Product* > ( itt->core ) );
		int_grad_grad = true;  continue;
	}  // end of for over vector 'result'

	if ( int_psi and not int_psi_psi and not int_dpsi and
	     not int_dpsi_dpsi and not int_grad_grad and not int_psi_dpsi )
	{	this->cas = 1;  // case one : int_psi
		assert ( bf .size() == 1 );
		assert ( result .size() == 1 );
		this->result_of_integr .resize ( 1 );
		this->result_of_integr [0] .resize ( 3 );
		for ( size_t i = 0; i < 3; i++ ) this->result_of_integr [0] [i] .resize ( 1 );
		this->selector .push_back ( 0 ); // selector = { 0 }
		return;                                                                        }

	if ( not int_psi and int_psi_psi and not int_dpsi and
	     not int_dpsi_dpsi and not int_grad_grad and not int_psi_dpsi )
	{	this->cas = 2;  // case two : int_psi_psi
		assert ( bf .size() == 2 );
		assert ( result .size() == 1 );
		this->result_of_integr .resize ( 3 );
		for ( size_t i = 0; i < 3; i++ )
		{	this->result_of_integr [i] .resize ( 3 );
			for ( size_t j = 0; j < 3; j++ )
				this->result_of_integr [i] [j] .resize ( 1 );  }
		Function::Product * prod = tag::Util::assert_cast
			< Function::Core*, Function::Product* > ( result [0] .core );
		std::forward_list < Function > ::const_iterator
			itt = prod->factors .begin();
		assert ( itt != prod->factors .end() );
		assert ( itt->core == bf [0] .core );
		itt++;  assert ( itt != prod->factors .end() );
		assert ( itt->core == bf [1] .core );
		itt++;  assert ( itt == prod->factors .end() );
		this->selector .push_back ( 0 ); // selector = { 0 }
		return;                                                             }
	
	if ( not int_psi and not int_psi_psi and int_dpsi and
	     not int_dpsi_dpsi and not int_grad_grad and not int_psi_dpsi )
	{	this->cas = 3;  // case three : int_dpsi
		Function xy = Manifold::working .coordinates();
		assert ( xy .nb_of_components() == 2 );
		assert ( bf .size() == 1 );
		this->result_of_integr .resize ( 1 );
		this->result_of_integr [0] .resize ( 3 );
		for ( size_t i = 0; i < 3; i++ )
			this->result_of_integr [0] [i] .resize ( 2 );
		// we analyse again the requested 'result' to see which expressions are wanted
		assert ( result .size() > 0 );
		assert ( result .size() <= 2 );
		Function::DelayedDerivative * dd0 = tag::Util::assert_cast
			< Function::Core*, Function::DelayedDerivative* > ( result [0] .core );
		if ( dd0->variable .core == xy[0] .core )  selector .push_back ( 0 );
		else
		{	assert ( dd0->variable .core == xy[1] .core );
			this->selector .push_back ( 1 );                }
		if ( result .size() == 1 ) return;
		Function::DelayedDerivative * dd1 = tag::Util::assert_cast
			< Function::Core*, Function::DelayedDerivative* > ( result [1] .core );
		if ( dd0->variable .core == dd1->variable .core ) goto error;
		if ( dd1->variable .core == xy[0] .core )  selector .push_back ( 0 );
		else
		{	assert ( dd1->variable .core == xy[1] .core );
			this->selector .push_back ( 1 );                }
		return;                                                                    }
	
	if ( not int_psi and not int_psi_psi and not int_dpsi and
	     int_dpsi_dpsi and not int_grad_grad and not int_psi_dpsi )
	{	this->cas = 4;  // case four : int_dpsi_dpsi
		assert ( bf .size() == 2 );
		this->result_of_integr .resize ( 3 );
		for ( size_t i = 0; i < 3; i++ )
		{	this->result_of_integr [i] .resize ( 3 );
			for ( size_t j = 0; j < 3; j++ )
				this->result_of_integr [i] [j] .resize ( 4 );  }
		// we analyse again the requested 'result' to see which expressions are wanted
		assert ( result .size() >= 1 );
		Function xy = Manifold::working .coordinates();
		assert ( xy .nb_of_components() == 2 );
		for ( size_t i = 0; i < result .size(); i++ )
		{	Function f = result [i];
			Function::Product * prod = tag::Util::assert_cast
				< Function::Core*, Function::Product* > ( f .core );
			std::forward_list < Function > ::const_iterator
				itt = prod->factors .begin();
			assert ( itt != prod->factors .end() );
			Function::DelayedDerivative * d1 = tag::Util::assert_cast
				< Function::Core*, Function::DelayedDerivative* > ( itt->core );
			itt++;  assert ( itt != prod->factors .end() );
			Function::DelayedDerivative * d2 = tag::Util::assert_cast
				< Function::Core*, Function::DelayedDerivative* > ( itt->core );
			itt++;  assert ( itt == prod->factors .end() );
			Function::Core * v1 = d1->variable .core;
			Function::Core * v2 = d2->variable .core;
			if ( v1 == xy [0] .core )
				if ( v2 == xy [0] .core )
					this->selector .push_back ( 0 );
				else
				{	assert ( v2 == xy [1] .core );
					this->selector .push_back ( 1 );     }
			else
			{	assert ( v1 == xy [1] .core );
				if ( v2 == xy [0] .core )
					this->selector .push_back ( 2 );
				else
				{	assert ( v2 == xy [1] .core );
					this->selector .push_back ( 3 );     }  }                        }
		return;                                                                   }
	
	if ( not int_psi and not int_psi_psi and not int_dpsi and
	     not int_dpsi_dpsi and int_grad_grad and not int_psi_dpsi )
	{	this->cas = 5;  // case five : int_grad_grad
		assert ( bf .size() == 2 );
		this->result_of_integr .resize ( 3 );
		for ( size_t i = 0; i < 3; i++ )
		{	this->result_of_integr [i] .resize ( 3 );
			for ( size_t j = 0; j < 3; j++ )
				this->result_of_integr [i] [j] .resize ( 1 );  }
		Function xy = Manifold::working .coordinates();
		assert ( xy .nb_of_components() == 2 );
		assert ( result .size() == 1 );
		Function f = result [0];
		Function::Sum * s = dynamic_cast < Function::Sum* > ( f .core );
		if ( s == nullptr ) goto error;
		size_t i = 0;
		for ( std::forward_list < Function > ::const_iterator it_s = s->terms .begin();
					it_s != s->terms .end(); it_s++                                          )
		{	Function::Product * fp = tag::Util::assert_cast
				< Function::Core*, Function::Product* > ( it_s->core );
			std::forward_list < Function > ::const_iterator it_p = fp->factors .begin();
			assert ( it_p != fp->factors .end() );
			Function::DelayedDerivative * dd1 = tag::Util::assert_cast
				< Function::Core*, Function::DelayedDerivative* > ( it_p->core );
			assert ( dynamic_cast < Function::MereSymbol* > ( dd1->base .core ) );
			assert ( dd1->base .core == bf [0] .core );
			assert ( dd1->variable .core == xy [i] .core );
			it_p ++;  assert ( it_p != fp->factors .end() );
			Function::DelayedDerivative * dd2 = tag::Util::assert_cast
				< Function::Core*, Function::DelayedDerivative* > ( it_p->core );
			assert ( dynamic_cast < Function::MereSymbol* > ( dd2->base .core ) );
			assert ( dd2->base .core == bf [1] .core );
			assert ( dd2->variable .core == xy [i] .core );
			it_p ++;  assert ( it_p == fp->factors .end() );
			i++;                                                                           }
		assert ( i == 2 );
		this->selector .push_back ( 0 ); // selector = { 0 }
		return;                                                                            }
		
	if ( not int_psi and not int_psi_psi and not int_dpsi and
	     not int_dpsi_dpsi and not int_grad_grad and int_psi_dpsi )
	{	this->cas = 6;  // case six : int_psi_dpsi
		assert ( result .size() > 0 );
		assert ( result .size() <= 2 );
		Function xy = Manifold::working .coordinates();
		assert ( xy .nb_of_components() == 2 );
		this->result_of_integr .resize ( 3 );
		for ( size_t i = 0; i < 3; i++ )
		{	this->result_of_integr [i] .resize ( 3 );
			for ( size_t j = 0; j < 3; j++ )
				this->result_of_integr [i] [j] .resize ( 4 );  }
		// we analyse again the requested 'result' to see which expressions are wanted
		assert ( result .size() > 0 );
		assert ( result .size() <= 2 );
		for ( size_t i = 0; i < result .size(); i++ )
		{	Function::Product * prod = tag::Util::assert_cast
				< Function::Core*, Function::Product* > ( result [i] .core );
			std::forward_list < Function > ::const_iterator itt = prod->factors .begin();
			assert ( itt != prod->factors .end() );
			Function::DelayedDerivative * dd = dynamic_cast
				< Function::DelayedDerivative* > ( itt->core );
			if ( dd )
			{	assert ( dd->base .core == bf [0] .core );
				itt++;  assert ( itt != prod->factors .end() );
				Function::MereSymbol * ms = dynamic_cast
					< Function::MereSymbol* > ( itt->core );
				assert ( ms == bf [1] .core );
				Function::Core * var = dd->variable .core;
				if ( var == xy [0] .core )
					this->selector .push_back ( 1 );
				else
				{	assert ( var == xy [1] .core );
					this->selector .push_back ( 3 );  }
				continue;                                        }
			Function::MereSymbol * ms = tag::Util::assert_cast
				< Function::Core*, Function::MereSymbol* > ( itt->core );
			assert ( ms == bf [0] .core );
			itt++;  assert ( itt != prod->factors .end() );
			dd = tag::Util::assert_cast
				< Function::Core*, Function::DelayedDerivative* > ( itt->core );
			assert ( dd->base .core == bf [1] .core );
			Function::Core * var = dd->variable .core;
			if ( var == xy [0] .core )
				this->selector .push_back ( 0 );
			else
			{	assert ( var == xy [1] .core );
				this->selector .push_back ( 2 );  }
			itt++;  assert ( itt == prod->factors .end() );                                }
		return;                                                                             }
	
	if ( int_psi and int_psi_psi and not int_dpsi and
	     not int_dpsi_dpsi and not int_grad_grad and not int_psi_dpsi )
	{	this->cas = 7;  // case seven : int_psi and int_psi_psi
		this->result_of_integr .resize ( 3 );
		for ( size_t i = 0; i < 3; i++ )
		{	this->result_of_integr [i] .resize ( 3 );
			for ( size_t j = 0; j < 3; j++ )
				this->result_of_integr [i] [j] .resize ( 3 );  }
		// we analyse again the requested 'result' to see which expressions are wanted
		assert ( result .size() == 2 );
		Function f0 = result [0], f1 = result [1];
		Function::MereSymbol * ff0 = dynamic_cast < Function::MereSymbol* > ( f0 .core );
		if ( ff0 )
		{	assert ( dynamic_cast < Function::Product* > ( f1 .core ) );
			#ifndef NDEBUG
			Function::Product * prod = tag::Util::assert_cast
				< Function::Core*, Function::Product* > ( f1 .core );
			std::forward_list < Function > ::const_iterator itt = prod->factors .begin();
			assert ( itt != prod->factors .end() );
			Function::MereSymbol * ms = dynamic_cast < Function::MereSymbol* > ( itt->core );
			assert ( ms );
			assert ( ms == bf [0] .core );
			itt++;  assert ( itt != prod->factors .end() );
			ms = dynamic_cast < Function::MereSymbol* > ( itt->core );
			assert ( ms );
			assert ( ms == bf [1] .core );
			itt++;  assert ( itt == prod->factors .end() );
			#endif
			if ( ff0 == bf [0] .core ) this->selector = { 0, 2 };
			else
			{	assert ( ff0 == bf [1] .core );  this->selector = { 1, 2 };  }                 }
		else
		{	assert ( dynamic_cast < Function::Product* > ( f0 .core ) );
			#ifndef NDEBUG
			Function::Product * prod = tag::Util::assert_cast
				< Function::Core*, Function::Product* > ( f0 .core );
			std::forward_list < Function > ::const_iterator itt = prod->factors .begin();
			assert ( itt != prod->factors .end() );
			Function::MereSymbol * ms = dynamic_cast < Function::MereSymbol* > ( itt->core );
			assert ( ms );
			assert ( ms == bf [0] .core );
			itt++;  assert ( itt != prod->factors .end() );
			ms = dynamic_cast < Function::MereSymbol* > ( itt->core );
			assert ( ms );
			assert ( ms == bf [1] .core );
			itt++;  assert ( itt == prod->factors .end() );
			#endif
			Function::MereSymbol * ff1 = dynamic_cast < Function::MereSymbol* > ( f1 .core );
			if ( ff1 == bf [0] .core ) this->selector = { 2, 0 };
			else
			{	assert ( ff1 == bf [1] .core );  this->selector = { 2, 1 };  }                  }
		return;                                                                                }

	if ( int_psi and not int_psi_psi and int_dpsi and
	     not int_dpsi_dpsi and not int_grad_grad and not int_psi_dpsi )
	{	this->cas = 8;  // case eight : int_psi and int_dpsi
		assert ( bf .size() == 1 );
		this->result_of_integr .resize ( 1 );
		this->result_of_integr [0] .resize ( 3 );
		for ( size_t i = 0; i < 3; i++ )
				this->result_of_integr [0] [i] .resize ( 3 );
		Function xy = Manifold::working .coordinates();
		assert ( xy .nb_of_components() == 2 );
		// we analyse again the requested 'result' to see which expressions are wanted
		assert ( result .size() >= 2 );
		assert ( result .size() <= 3 );
		for ( size_t i = 0; i < result .size(); i++ )
		{	Function f = result [i];
			Function::MereSymbol * ms =
				dynamic_cast < Function::MereSymbol* > ( f .core );
			if ( ms )
			{	selector .push_back ( 0 );  continue;  }
			Function::DelayedDerivative * dd = tag::Util::assert_cast
				< Function::Core*, Function::DelayedDerivative* > ( f .core );
			Function::Core * v = dd->variable .core;
			if ( v == xy [0] .core )
				this->selector .push_back ( 1 );
			else
			{	assert ( v == xy [1] .core );
				this->selector .push_back ( 2 );     }                             }
		return;                                                                   }

	if ( not int_psi and not int_psi_psi and int_dpsi and
	     int_dpsi_dpsi and not int_grad_grad and not int_psi_dpsi )
	{	this->cas = 9;  // case nine : int_dpsi and int_dpsi_dpsi
		assert ( bf .size() == 2 );
		this->result_of_integr .resize ( 3 );
		for ( size_t i = 0; i < 3; i++ )
		{	this->result_of_integr [i] .resize ( 3 );
			for ( size_t j = 0; j < 3; j++ )
				this->result_of_integr [i] [j] .resize ( 8 );  }
		Function xy = Manifold::working .coordinates();
		assert ( xy .nb_of_components() == 2 );
		// we analyse again the requested 'result' to see which expressions are wanted
		assert ( result .size() >= 2 );
		for ( size_t i = 0; i < result .size(); i++ )
		{	Function f = result [i];
			Function::DelayedDerivative * dd = 
				dynamic_cast < Function::DelayedDerivative* > ( f .core );
			if ( dd )
			{	assert ( dynamic_cast < Function::MereSymbol* > ( dd->base .core ) );
				if ( dd->variable .core == xy[0] .core )
					if ( dd->base .core == bf [0] .core )
						selector .push_back ( 0 );
					else  selector .push_back ( 1 );
				else
				{	assert ( dd->variable .core == xy[1] .core );
					if ( dd->base .core == bf [0] .core )
						selector .push_back ( 2 );
					else  selector .push_back ( 3 );              }
				continue;                                                              }
			Function::Product * prod = tag::Util::assert_cast
				< Function::Core*, Function::Product* > ( f .core );
			std::forward_list < Function > ::const_iterator
				itt = prod->factors .begin();
			assert ( itt != prod->factors .end() );
			Function::DelayedDerivative * d1 = tag::Util::assert_cast
				< Function::Core*, Function::DelayedDerivative* > ( itt->core );
			itt++;  assert ( itt != prod->factors .end() );
			Function::DelayedDerivative * d2 = tag::Util::assert_cast
				< Function::Core*, Function::DelayedDerivative* > ( itt->core );
			itt++;  assert ( itt == prod->factors .end() );
			Function::Core * v1 = d1->variable .core;
			Function::Core * v2 = d2->variable .core;
			if ( v1 == xy [0] .core )
				if ( v2 == xy [0] .core )
					this->selector .push_back ( 4 );
				else
				{	assert ( v2 == xy [1] .core );
					this->selector .push_back ( 5 );     }
			else
			{	assert ( v1 == xy [1] .core );
				if ( v2 == xy [0] .core )
					this->selector .push_back ( 6 );
				else
				{	assert ( v2 == xy [1] .core );
					this->selector .push_back ( 7 );     }  }                        }
		return;                                                                   }
	
	if ( not int_psi and int_psi_psi and not int_dpsi and
	     int_dpsi_dpsi and not int_grad_grad and not int_psi_dpsi )
	{	this->cas = 10;  // case ten : int_psi_psi and int_dpsi_dpsi
		assert ( bf .size() == 2 );
		this->result_of_integr .resize ( 3 );
		for ( size_t i = 0; i < 3; i++ )
		{	this->result_of_integr [i] .resize ( 3 );
			for ( size_t j = 0; j < 3; j++ )
				this->result_of_integr [i] [j] .resize ( 5 );  }
		// we analyse again the requested 'result' to see which expressions are wanted
		assert ( result .size() >= 2 );
		Function xy = Manifold::working .coordinates();
		assert ( xy .nb_of_components() == 2 );
		for ( size_t i = 0; i < result .size(); i++ )
		{	Function f = result [i];
			Function::Product * prod = tag::Util::assert_cast
				< Function::Core*, Function::Product* > ( f .core );
			std::forward_list < Function > ::const_iterator
				itt = prod->factors .begin();
			assert ( itt != prod->factors .end() );
			Function::MereSymbol * ms = dynamic_cast
				< Function::MereSymbol* > ( itt->core );
			if ( ms )
			{	assert ( ms == bf [0] .core );
				itt++;  assert ( itt != prod->factors .end() );
				ms = tag::Util::assert_cast
					< Function::Core*, Function::MereSymbol* > ( itt->core );
				assert ( ms == bf [1] .core );
				this->selector .push_back ( 0 );
				continue;                                                   }
			Function::DelayedDerivative * d1 = tag::Util::assert_cast
				< Function::Core*, Function::DelayedDerivative* > ( itt->core );
			itt++;  assert ( itt != prod->factors .end() );
			Function::DelayedDerivative * d2 = tag::Util::assert_cast
				< Function::Core*, Function::DelayedDerivative* > ( itt->core );
			itt++;  assert ( itt == prod->factors .end() );
			assert ( d1->base .core == bf [0] .core );
			assert ( d2->base .core == bf [1] .core );
			Function::Core * v1 = d1->variable .core;
			Function::Core * v2 = d2->variable .core;
			if ( v1 == xy [0] .core )
				if ( v2 == xy [0] .core )
					this->selector .push_back ( 1 );
				else
				{	assert ( v2 == xy [1] .core );
					this->selector .push_back ( 2 );     }
			else
			{	assert ( v1 == xy [1] .core );
				if ( v2 == xy [0] .core )
					this->selector .push_back ( 3 );
				else
				{	assert ( v2 == xy [1] .core );
					this->selector .push_back ( 4 );     }  }                        }
		return;                                                                   }
	
	if ( not int_psi and int_psi_psi and not int_dpsi and
	     not int_dpsi_dpsi and int_grad_grad and not int_psi_dpsi )
	{	this->cas = 11;  // case eleven : int_psi_psi and int_grad_grad
		assert ( bf .size() == 2 );
		this->result_of_integr .resize ( 3 );
		for ( size_t i = 0; i < 3; i++ )
		{	this->result_of_integr [i] .resize ( 3 );
			for ( size_t j = 0; j < 3; j++ )
				this->result_of_integr [i] [j] .resize ( 2 );  }
		// we analyse again the requested 'result' to see which expressions are wanted
		assert ( result .size() == 2 );
		Function xy = Manifold::working .coordinates();
		assert ( xy .nb_of_components() == 2 );
		for ( size_t i = 0; i < result .size(); i++ )
		{	Function f = result [i];
			Function::Product * prod = 
				dynamic_cast < Function::Product* > ( f .core );
			if ( prod )
			{	std::forward_list < Function > ::const_iterator
					itt = prod->factors .begin();
				assert ( itt != prod->factors .end() );
				Function::MereSymbol * ms = tag::Util::assert_cast
					< Function::Core*, Function::MereSymbol* > ( itt->core );
				assert ( ms == bf [0] .core );
				itt++;  assert ( itt != prod->factors .end() );
				ms = tag::Util::assert_cast
					< Function::Core*, Function::MereSymbol* > ( itt->core );
				assert ( ms == bf [1] .core );
				this->selector .push_back ( 0 );
				continue;                                                   }
			Function::Sum * s = dynamic_cast < Function::Sum* > ( f .core );
			if ( s == nullptr ) goto error;
			size_t ii = 0;
			for ( std::forward_list < Function > ::const_iterator it_s = s->terms .begin();
						it_s != s->terms .end(); it_s++                                          )
			{	Function::Product * fp = tag::Util::assert_cast
					< Function::Core*, Function::Product* > ( it_s->core );
				std::forward_list < Function > ::const_iterator it_p = fp->factors .begin();
				assert ( it_p != fp->factors .end() );
				Function::DelayedDerivative * dd1 = tag::Util::assert_cast
					< Function::Core*, Function::DelayedDerivative* > ( it_p->core );
				assert ( dynamic_cast < Function::MereSymbol* > ( dd1->base .core ) );
				assert ( dd1->base .core == bf [0] .core );
				assert ( dd1->variable .core == xy [ii] .core );
				it_p ++;  assert ( it_p != fp->factors .end() );
				Function::DelayedDerivative * dd2 = tag::Util::assert_cast
					< Function::Core*, Function::DelayedDerivative* > ( it_p->core );
				assert ( dynamic_cast < Function::MereSymbol* > ( dd2->base .core ) );
				assert ( dd2->base .core == bf [1] .core );
				assert ( dd2->variable .core == xy [ii] .core );
				it_p ++;  assert ( it_p == fp->factors .end() );
				ii++;                                                                        }
			assert ( ii == 2 );
			this->selector .push_back ( 1 );                                                 }
		return;                                                                              }
	
	if ( int_psi and int_psi_psi and int_dpsi and
	     int_dpsi_dpsi and not int_grad_grad and not int_psi_dpsi )
	{	this->cas = 12;  // case twelve : int_psi, int_psi_psi, int_dpsi and int_dpsi_dpsi
		assert ( bf .size() == 2 );
		this->result_of_integr .resize ( 3 );
		for ( size_t i = 0; i < 3; i++ )
		{	this->result_of_integr [i] .resize ( 3 );
			for ( size_t j = 0; j < 3; j++ )
				this->result_of_integr [i] [j] .resize ( 11 );  }
		// we analyse again the requested 'result' to see which expressions are wanted
		assert ( result .size() >= 4 );
		Function xy = Manifold::working .coordinates();
		assert ( xy .nb_of_components() == 2 );
		for ( size_t i = 0; i < result .size(); i++ )
		{	Function f = result [i];
			Function::MereSymbol * ms =
				dynamic_cast < Function::MereSymbol* > ( f .core );
			if ( ms )
			{	if ( ms == bf [0] .core ) this->selector .push_back ( 0 );
				else
					{	assert ( ms == bf [1] .core );  this->selector .push_back ( 1 );  }
				continue;                                                               }
			Function::Product * prod = dynamic_cast < Function::Product* > ( f .core );
			if ( prod )
			{	std::forward_list < Function > ::const_iterator
					itt = prod->factors .begin();
				assert ( itt != prod->factors .end() );
				ms = dynamic_cast < Function::MereSymbol* > ( itt->core );
				if ( ms )
				{	assert ( ms == bf [0] .core );
					itt++;  assert ( itt != prod->factors .end() );
					ms = tag::Util::assert_cast
						< Function::Core*, Function::MereSymbol* > ( itt->core );
					assert ( ms == bf [1] .core );
					this->selector .push_back ( 2 );
					continue;                                                   }
				Function::DelayedDerivative * d1 = tag::Util::assert_cast
					< Function::Core*, Function::DelayedDerivative* > ( itt->core );
				itt++;  assert ( itt != prod->factors .end() );
				Function::DelayedDerivative * d2 = tag::Util::assert_cast
					< Function::Core*, Function::DelayedDerivative* > ( itt->core );
				itt++;  assert ( itt == prod->factors .end() );
				assert ( d1->base .core == bf [0] .core );
				assert ( d2->base .core == bf [1] .core );
				Function::Core * v1 = d1->variable .core;
				Function::Core * v2 = d2->variable .core;
				if ( v1 == xy [0] .core )
					if ( v2 == xy [0] .core )
						this->selector .push_back ( 3 );
					else
						{	assert ( v2 == xy [1] .core );
							this->selector .push_back ( 4 );     }
				else
				{	assert ( v1 == xy [1] .core );
					if ( v2 == xy [0] .core )
						this->selector .push_back ( 5 );
					else
						{	assert ( v2 == xy [1] .core );
							this->selector .push_back ( 6 );     }  }
				continue;                                                             }
			Function::DelayedDerivative * dd = tag::Util::assert_cast
				< Function::Core*, Function::DelayedDerivative* > ( f .core );
			if ( dd->base .core == bf [0] .core )
				if ( dd->variable .core == xy[0] .core )  selector .push_back ( 7 );
				else
				{	assert ( dd->variable .core == xy[1] .core );
					this->selector .push_back ( 9 );               }
			else
			{	assert ( dd->base .core == bf [1] .core );
				if ( dd->variable .core == xy[0] .core )  selector .push_back ( 8 );
				else
				{	assert ( dd->variable .core == xy[1] .core );
					this->selector .push_back ( 10 );              }                     }  }
		return;                                                                          }
	
	std::cout << "boolean variables at end of pre_compute : ";
	std::cout << int_psi << int_psi_psi << int_dpsi <<
		int_dpsi_dpsi << int_grad_grad << int_psi_dpsi << std::endl;
	
	error :
	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "pre_compute unsuccessful" << std::endl;
	exit ( 1 );

}  // end of  FiniteElement::StandAlone::TypeOne::Triangle::pre_compute


