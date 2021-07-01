
#include "maniFEM.h"

using namespace maniFEM;

int main () {

	// we choose our (geometric) space dimension :
	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	
	// xyz is a map defined on our future mesh with values in RR3 :
	Function xyz = RR3.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );

	// we can extract components of xyz using the [] operator :
	Function x = xyz[0],  y = xyz[1],  z = xyz[2];

	Cell A ( tag::vertex );  x(A) = 0.;
	Cell B ( tag::vertex );  x(B) = 1.;

	Mesh AB ( tag::segment, A.reverse(), B, tag::divided_in, 10 );

	Mesh::Connected::OneDim * AB1d = tag::Util::assert_cast
		< Mesh::Core*, Mesh::Connected::OneDim* > ( AB.core );

	CellIterator::Over::VerticesOfConnectedOneDimMesh::ReverseOrder it ( AB1d );
	// it.reset ( );
	std::cout << "main " << B.core << " " << B.core->reverse_attr.core << std::endl;
  it.reset ( tag::start_at, B.core );
  assert ( it.in_range() );
	std::cout << x(it.deref()) << std::endl;
	it.advance();  assert ( it.in_range() );
	std::cout << x(it.deref()) << std::endl;
	it.advance();  assert ( it.in_range() );
	std::cout << x(it.deref()) << std::endl;

}
	
