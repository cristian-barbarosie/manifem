
// example presented in paragraph 3.22 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// a sharp cone

// the code shown in the manual does not work (yet)
// we fake the result by building by hand the triangles around the vertex

#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;
using namespace std;


void update_info_connected_one_dim ( const Mesh msh, const Cell start, const Cell stop )

// 'start' and 'stop' are positive vertices (may be one and the same)

{	assert ( start .dim() == 0 );
	assert ( stop .dim() == 0 );
	assert ( start .is_positive() );
	assert ( stop .is_positive() );

	Mesh::Connected::OneDim * msh_core = tag::Util::assert_cast
		< Mesh::Core*, Mesh::Connected::OneDim* > ( msh .core );
	msh_core->first_ver = start .reverse();
	msh_core->last_ver = stop;
	// now we can use an iterator

	Mesh::Iterator it = msh .iterator ( tag::over_segments, tag::require_order );
	size_t n = 0;
	for ( it .reset(); it .in_range(); it++ ) n++;
	msh_core->nb_of_segs = n;                                                   }


int main ()

{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz [0], y = xyz [1], z = xyz [2];
	double seg_size = 0.13;
	// std::cout << "segment size : ";
	// std::cin >> seg_size;

	Manifold cone_manif = RR3 .implicit ( x*x + y*y == z*z );
	
	Cell O ( tag::vertex );  x (O) = 0.;  y (O) = 0.;  z (O) = 0.;
	Cell A ( tag::vertex );  x (A) = 0.7*seg_size;  y (A) = 0.;  z (A) = 0.7*seg_size;
	cone_manif .project (A);
	Cell B ( tag::vertex );  x (B) = 0.2*seg_size;  y (B) = 0.6*seg_size;  z (B) = 0.7*seg_size;
	cone_manif .project (B);
	Cell C ( tag::vertex );  x (C) = -0.5*seg_size;  y (C) = 0.4*seg_size;  z (C) = 0.7*seg_size;
	cone_manif .project (C);
	Cell D ( tag::vertex );  x (D) = -0.5*seg_size;  y (D) = -0.4*seg_size;  z (D) = 0.7*seg_size;
	cone_manif .project (D);
	Cell E ( tag::vertex );  x (E) = 0.2*seg_size;  y (E) = -0.6*seg_size;  z (E) = 0.7*seg_size;
	cone_manif .project (E);
	Cell OA ( tag::segment, O .reverse(), A );
	Cell OB ( tag::segment, O .reverse(), B );
	Cell OC ( tag::segment, O .reverse(), C );
	Cell OD ( tag::segment, O .reverse(), D );
	Cell OE ( tag::segment, O .reverse(), E );
	Cell AB ( tag::segment, A .reverse(), B );
	Cell BC ( tag::segment, B .reverse(), C );
	Cell CD ( tag::segment, C .reverse(), D );
	Cell DE ( tag::segment, D .reverse(), E );
	Cell EA ( tag::segment, E .reverse(), A );
	Cell OAB ( tag::triangle, OA, AB, OB .reverse() );
	Cell OBC ( tag::triangle, OB, BC, OC .reverse() );
	Cell OCD ( tag::triangle, OC, CD, OD .reverse() );
	Cell ODE ( tag::triangle, OD, DE, OE .reverse() );
	Cell OEA ( tag::triangle, OE, EA, OA .reverse() );
	Mesh ponta ( tag::whose_core_is,
	       new Mesh::Fuzzy ( tag::of_dim, 3, tag::minus_one, tag::one_dummy_wrapper ),
	       tag::freshly_created, tag::is_positive                                      );
	OAB .add_to_mesh ( ponta );
	OBC .add_to_mesh ( ponta );
	OCD .add_to_mesh ( ponta );
	ODE .add_to_mesh ( ponta );
	OEA .add_to_mesh ( ponta );

	Mesh small_circle ( tag::whose_core_is,
         new Mesh::Connected::OneDim ( tag::with, 5, tag::segments, tag::one_dummy_wrapper ),
         tag::freshly_created, tag::is_positive                                               );
	AB .add_to_mesh ( small_circle, tag::do_not_bother );
	BC .add_to_mesh ( small_circle, tag::do_not_bother );
	CD .add_to_mesh ( small_circle, tag::do_not_bother );
	DE .add_to_mesh ( small_circle, tag::do_not_bother );
	EA .add_to_mesh ( small_circle, tag::do_not_bother );
	// the meaning of tag::do_not_bother is explained at the end of paragraph 11.6 in the manual
	update_info_connected_one_dim ( small_circle, A, A );

	Cell V ( tag::vertex );  x (V) = 1.;  y (V) = 0.;  z (V) = 1.;
	cone_manif .implicit ( z == 1. );
	std::vector < double > tau = { 0., 1., 0. };
	Mesh big_circle ( tag::frontal, tag::start_at, V, tag::towards, tau,
	                  tag::desired_length, seg_size                     );
	Mesh two_circles ( tag::join, small_circle .reverse(), big_circle );

	cone_manif .set_as_working_manifold();
	tau = { -1., 0., -1. };
	Mesh cone_up ( tag::frontal, tag::boundary, two_circles,
	               tag::start_at, V, tag::towards, tau,
	               tag::desired_length, seg_size            );

	Mesh cone ( tag::join, ponta, cone_up );
	cone .export_to_file ( tag::msh, "cone.msh");
	cout << "produced file cone.msh" << endl;

}  // end of main
