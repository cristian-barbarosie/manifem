

#include "maniFEM.h"

using namespace maniFEM;

int main () {

	// we choose our (geometric) space dimension :
	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	
	// xyz is a map defined on our future mesh with values in RR3 :
	Function xyz = RR3.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );

	// we can extract components of xyz using the [] operator :
	Function x = xyz[0],  y = xyz[1],  z = xyz[2];

	// Let's build a rectangular mesh. First, the four corners :
	Cell SW ( tag::vertex );  x(SW) =  0.;  y(SW) =  0.;  z(SW) =  0.;
	Cell SE ( tag::vertex );  x(SE) =  1.;  y(SE) =  0.;  z(SE) =  0.;
	Cell NE ( tag::vertex );  x(NE) =  1.;  y(NE) =  1.;  z(NE) =  0.;
	Cell NW ( tag::vertex );  x(NW) =  0.;  y(NW) =  1.;  z(NW) =  1.;

	// now build the four sides of the rectangle :
	Mesh south ( tag::segment, SW.reverse(), SE, tag::divided_in, 10 );
	Mesh east  ( tag::segment, SE.reverse(), NE, tag::divided_in, 10 );
	Mesh north ( tag::segment, NE.reverse(), NW, tag::divided_in, 10 );
	Mesh west  ( tag::segment, NW.reverse(), SW, tag::divided_in, 10 );

	// and now the rectangle :
	Mesh rect_mesh ( tag::rectangle, south, east, north, west );

	Cell PQ = south.cell_in_front_of ( SW );
	Cell Q = PQ.tip();
	Cell square = rect_mesh.cell_behind ( PQ );
	Cell QR = square.boundary().cell_in_front_of ( Q );
	Cell R = QR.tip();
	Cell RS = square.boundary().cell_in_front_of ( R );
	Cell A ( tag::vertex );  x(A) = 0.;  y(A) = -0.1;
	Cell B ( tag::vertex );  x(B) = 0.1; y(B) = -0.1;
	Cell PA ( tag::segment, SW.reverse(), A );
	Cell AB ( tag::segment, A.reverse(), B );
	Cell BQ ( tag::segment, B.reverse(), Q );
	Cell PQBA ( tag::rectangle, PQ, BQ.reverse(), AB.reverse(), PA.reverse() );
	PQBA.reverse().add_to_mesh ( rect_mesh );

	std::cout << "main CellIterator::AroundCell::OfCodimTwo::OverVertices::NormalOrder::BuildReverseCells::AsTheyAre"	 << std::endl;

	{  // just a block of code for hiding 'it'

	std::cout << "--------------------------------------------------------" << std::endl;	

	Cell::Positive * SW_core = tag::Util::assert_cast < Cell::Core*, Cell::Positive* > ( SW.core );
  CellIterator::AroundCell::OfCodimTwo::OverVertices::NormalOrder::BuildReverseCells::AsTheyAre
		it ( rect_mesh.core, SW_core );

	std::cout << "around SW :" << std::endl;
	Cell SP = west.cell_in_front_of ( RS.tip() );
	for ( it.reset ( ); it.in_range(); it.advance() )
	{	Cell seg = it.deref();
		std::cout << "[("  << x(seg.base().reverse()) << "," << y(seg.base().reverse()) << ") -> (" << x(seg.tip()) << "," << y(seg.tip()) << ")]" << std::endl << std::flush;  }

	} {  // just a block of code for hiding 'it'

	std::cout << "--------------------------------------------------------" << std::endl;	

	Cell::Positive * Q_core = tag::Util::assert_cast < Cell::Core*, Cell::Positive* > ( Q.core );
  CellIterator::AroundCell::OfCodimTwo::OverVertices::NormalOrder::BuildReverseCells::AsTheyAre
		it ( rect_mesh.core, Q_core );

	std::cout << "around Q :" << std::endl;
	Cell QQ = south.cell_in_front_of ( Q );
	for ( it.reset ( ); it.in_range(); it.advance() )
	{	Cell seg = it.deref();
		std::cout << "[("  << x(seg.base().reverse()) << "," << y(seg.base().reverse()) << ") -> (" << x(seg.tip()) << "," << y(seg.tip()) << ")]" << std::endl << std::flush;  }

	} {  // just a block of code for hiding 'it'

	std::cout << "--------------------------------------------------------" << std::endl;	

	Cell::Positive * R_core = tag::Util::assert_cast < Cell::Core*, Cell::Positive* > ( R.core );
  CellIterator::AroundCell::OfCodimTwo::OverVertices::NormalOrder::BuildReverseCells::AsTheyAre
		it ( rect_mesh.core, R_core );

	std::cout << "around R :" << std::endl;
	for ( it.reset ( ); it.in_range(); it.advance() )
	{	Cell seg = it.deref();
		std::cout << "[("  << x(seg.base().reverse()) << "," << y(seg.base().reverse()) << ") -> (" << x(seg.tip()) << "," << y(seg.tip()) << ")]" << std::endl << std::flush;  }

	}  // just a block of code for hiding 'it'

	std::cout << "end of main" << std::endl << std::flush;
}
