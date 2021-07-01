
// example presented in paragraph 1.1 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// builds a twisted rectangle in 3D

#include "maniFEM.h"

using namespace maniFEM;
using namespace std;

int main () {

	// we choose our (geometric) space dimension :
	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	
	// xyz is a map defined on our future mesh with values in RR3 :
	Function xyz = RR3.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );

	// we can extract components of xyz using the [] operator :
	Function x = xyz[0],  y = xyz[1],  z = xyz[2];

	// Let's build a rectangular mesh. First, the four corners :
	Cell SW ( tag::vertex );  x(SW) = -1.;  y(SW) =  0.;  z(SW) =  0.;
	Cell SE ( tag::vertex );  x(SE) =  1.;  y(SE) =  0.;  z(SE) =  0.;
	Cell NE ( tag::vertex );  x(NE) =  1.;  y(NE) =  1.;  z(NE) =  0.;
	Cell NW ( tag::vertex );  x(NW) = -1.;  y(NW) =  1.;  z(NW) =  1.;

	// we access the coordinates of a point using the () operator :
	cout << "coordinates of NW : " << x(NW) << " " << y(NW) << " " << z(NW) << endl;
	
	// now build the four sides of the rectangle :
	Mesh south ( tag::segment, SW.reverse(), SE, tag::divided_in, 10 );
	Mesh east  ( tag::segment, SE.reverse(), NE, tag::divided_in, 10 );
	Mesh north ( tag::segment, NE.reverse(), NW, tag::divided_in, 10 );
	Mesh west  ( tag::segment, NW.reverse(), SW, tag::divided_in, 10 );

	// and now the rectangle :
	Mesh rect_mesh ( tag::pretty, tag::rectangle, south, east, north, west );
	 
	// We may want to visualize the resulting mesh.
	// Here is one way to export the mesh in the "msh" format :
	rect_mesh.export_msh ("rectangle.msh");
	// rect_mesh.draw_ps ("rectangle.eps");
	
	// let's define two symbolic functions
	Function f = x*x + 1./(5.+y), g = x*y;

	// and compute their integral on the rectangle, using Gauss quadrature with 9 points :
	FiniteElement fe ( tag::with_master, tag::quadrangle, tag::Lagrange, tag::of_degree, 1 );
	fe.set_integrator ( tag::Gauss, tag::quad_9 );
	// cout << "integral of " << f.repr() << " = " << fe.integrate ( f, tag::over, rect_mesh ) << endl;
	// cout << "integral of " << g.repr() << " = " << fe.integrate ( g, tag::over, rect_mesh ) << endl;

	// the above code does not work yet, but we are very close to getting it right ;-)
	
	// abstract variational formulation :
	// auto & uu = FunctionOnMesh::unknown ( u, "Lagrange degree one");
	// The values of 'uu' are related to the field 'u'.
	// On segments, 'uu' varies linearly. On a triangle, it varies
	// also linearly. On a rectangle, it is a polynomial of degree 1,
	// so it has a linear part and a bilinear one.
	// auto & w = FunctionOnMesh::test ( uu );
	// 'w' will be used for defining variational formulations in tandem with 'uu'

	// auto & integ = Integrator::gauss ("Q9");
	// integ nao recebe as coordenadas espaciais x e y
	// fe (mais abaixo) extrai essa informacao da formulacao variacional
	// cout << ( xx*yy ) .integrate ( malha, integ ) << endl;
	
/*
	auto & var_pb =
		( uu.deriv(xx)*w.deriv(xx) + uu.deriv(yy)*w.deriv(yy) ) .integrate(malha)
		== w.integrate(malha) + w.integrate(south);

	var_pb.prescribe_on (north);  uu == 0.;          w == 0.;
	var_pb.prescribe_on (east);   uu == xx*(1.-yy);  w == 0.;
	var_pb.prescribe_on (west);   uu == 0.;          w == 0.;

	auto & fe = FiniteElement::lagrange ("Q1", malha );
	// fe.set_integrator ("on cells of dimension", 1, "gauss", 3, "nodes");
	fe.set_integrator ("on cells of dimension", 2, "gauss", 9, "nodes");
	var_pb.set_finite_element ( fe );
	// A finite element should perform the following tasks.
	// Build a base in the (discrete) Hilbert space and index it;
	// the index should be somehow related to an enumeration of cells.
	// Take a variational problem and express it as a system of linear equations,
	// by replacing the unknown and the test functions by functions in the base
	// and by computing the respective integrals.
	// Take the solution of the system and interpret it as a function on the mesh.

	var_pb.discretize ();
*/
//	do_test ();
	
	 cout << "produced file rectangle.msh" << endl;
}
