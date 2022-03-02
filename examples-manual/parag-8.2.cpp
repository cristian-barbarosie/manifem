
// example presented in paragraph 8.2 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// exchanging information with gmsh


//  https://stackoverflow.com/questions/216823/how-to-trim-a-stdstring

#include <algorithm> 
#include <cctype>
#include <locale>

// trim from start (in place)
static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

// trim from start (copying)
static inline std::string ltrim_copy(std::string s) {
    ltrim(s);
    return s;
}

// trim from end (copying)
static inline std::string rtrim_copy(std::string s) {
    rtrim(s);
    return s;
}

// trim from both ends (copying)
static inline std::string trim_copy(std::string s) {
    trim(s);
    return s;
}

//----------------------------------------------------------------------------------//


#include "maniFEM.h"
#include "math.h"
#include "fstream"

using namespace maniFEM;

int main ()

{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz [0], y = xyz [1], z = xyz [2];

	const double rs = 1.;    // radius of the sphere
	const double rc = 0.45;  // radius of the cylinder
	const double seg_size = 0.1;
	
	Manifold cylinder = RR3 .implicit ( y*y + (z-0.5)*(z-0.5) == rc*rc );
	Manifold intersection = cylinder .implicit ( x*x + y*y + z*z == rs*rs );
	Cell start1 ( tag::vertex );  x ( start1 ) = 1.;  y ( start1 ) = 0.;  z ( start1 ) = 0.5 - rc;
	intersection .project ( start1 );
	Mesh circle_1 ( tag::frontal, tag::start_at, start1, tag::towards, { 0., 1., 0. },
	                tag::desired_length, seg_size                                          );

	Cell start2 ( tag::vertex );  x ( start2 ) = -1.;  y ( start2 ) = 0.;  z ( start2 ) = 0.5 - rc;
	intersection .project ( start2 );
	Mesh circle_2 ( tag::frontal, tag::start_at, start2, tag::towards, { 0., 1., 0. },
                  tag::desired_length, seg_size                                          );

	Mesh two_circles ( tag::join, circle_1 .reverse(), circle_2 );
	cylinder .set_as_working_manifold();
	Mesh cyl ( tag::frontal, tag::boundary, two_circles, tag::start_at, start1,
						 tag::towards, { -1., 0., 0. }, tag::desired_length, seg_size         );

	Mesh two_circles_rev ( tag::join, circle_1, circle_2 .reverse() );
	RR3 .implicit ( x*x + y*y + z*z == rs*rs );
	Mesh sph ( tag::frontal, tag::boundary, two_circles_rev, tag::start_at, start1,
             tag::towards, { 0., 0., -1. }, tag::desired_length, seg_size             );

	Mesh all ( tag::join, cyl, sph );

	std::map < int, Mesh > composed_mesh
		{ { 1, all }, { 2, sph }, { 3, cyl }, { 4, circle_1 }, { 5, circle_2 } };
	
	// export_msh (
	all .export_to_file ( tag::msh, "sphere-tunnel.msh");

	std::cout << "produced file sphere-tunnel.msh" << std::endl;

}  // end of main

//-----------------------------------------------------------------------------------------------


