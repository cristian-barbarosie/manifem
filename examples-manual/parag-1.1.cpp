
// example presented in paragraph 1.1 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// builds a twisted rectangle in 3D

#include "maniFEM.h"

using namespace maniFEM;
using namespace std;

int main ()

{	// we choose our (geometric) space dimension :
	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	
	// xyz is a map defined on our future mesh with values in RR3 :
	Function xyz = RR3 .build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );

	// we can extract components of xyz using the [] operator :
	Function x = xyz [0],  y = xyz [1],  z = xyz [2];

	// Let's build a rectangular mesh. First, the four corners :
	Cell SW ( tag::vertex );  x(SW) = -1.;  y(SW) =  0.;  z(SW) =  0.;
	Cell SE ( tag::vertex );  x(SE) =  1.;  y(SE) =  0.;  z(SE) =  0.;
	Cell NE ( tag::vertex );  x(NE) =  1.;  y(NE) =  1.;  z(NE) =  0.;
	Cell NW ( tag::vertex );  x(NW) = -1.;  y(NW) =  1.;  z(NW) =  1.;

	// alternative syntax for declaring a vertex :
	// Cell SW ( tag::vertex, tag::of_coords, { -1., 0., 0. } )
	
	// we access the coordinates of a point using the () operator :
	cout << "coordinates of NW : " << x(NW) << " " << y(NW) << " " << z(NW) << endl;
	
	// now build the four sides of the rectangle :
	Mesh south ( tag::segment, SW .reverse(), SE, tag::divided_in, 10 );
	Mesh east  ( tag::segment, SE .reverse(), NE, tag::divided_in, 10 );
	Mesh north ( tag::segment, NE .reverse(), NW, tag::divided_in, 10 );
	Mesh west  ( tag::segment, NW .reverse(), SW, tag::divided_in, 10 );

	// and now the rectangle :
	Mesh rect_mesh ( tag::rectangle, south, east, north, west );
	 
	// we may want to visualize the resulting mesh
	// here is one way to export the mesh in the "msh" format :
	rect_mesh .export_to_file ( tag::msh, "rectangle.msh");

	// or, directly draw a postscript file :
	// rect_mesh.draw_ps ("rectangle.eps")
	
	// let's define two symbolic functions
	Function f = x*x + 1/(5+y), g = x*y;

	// and compute their integral on the rectangle, using Gauss quadrature with 9 points :
	Integrator integr ( tag::Gauss, tag::quad_9 );

	cout << "integral of f " << integr ( f, tag::on, rect_mesh ) << endl;
	cout << "integral of g " << integr ( g, tag::on, rect_mesh ) << endl;
		
	cout << "produced file rectangle.msh" << endl;

}  // end of main
