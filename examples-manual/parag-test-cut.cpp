
#include "maniFEM.h"
#include <set>
#include <fstream>

using namespace maniFEM;


//-----------------------------------------------------------------------------------//


int main ()
	
{	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0],  y = xy[1];

	Cell X ( tag::vertex );  x(X) = 0.;  y(X) = 0.;
	Cell Y ( tag::vertex );  x(Y) = 1.;  y(Y) = 0.;
	Cell Z ( tag::vertex );  x(Z) = 1.;  y(Z) = 1.;
	Cell T ( tag::vertex );  x(T) = 0.;  y(T) = 1.;

	Mesh XY ( tag::segment, X.reverse(), Y, tag::divided_in, 2 );
	Mesh YZ ( tag::segment, Y.reverse(), Z, tag::divided_in, 2 );
	Mesh ZT ( tag::segment, Z.reverse(), T, tag::divided_in, 2 );
	Mesh TX ( tag::segment, T.reverse(), X, tag::divided_in, 2 );
	XY.core->name = "XY";
	YZ.core->name = "YZ";
	ZT.core->name = "ZT";
	TX.core->name = "TX";

	Mesh rect_mesh ( tag::rectangle, XY, YZ, ZT, TX );
	rect_mesh.core->name = "rect_mesh";

	Cell A = X;
	Cell AB = XY .cell_in_front_of ( A, tag::surely_exists );
	Cell ABCD = rect_mesh .cell_behind ( AB, tag::surely_exists );
	Cell B = AB .tip();
	Cell BC = ABCD .boundary() .cell_in_front_of ( B );
	Cell C = BC .tip();
	Cell CD = ABCD .boundary() .cell_in_front_of ( C );
	Cell D = CD .tip();
	Cell DA = ABCD .boundary() .cell_in_front_of ( D );
	assert ( DA.tip() == A );

	A.core->name = "A";
	B.core->name = "B";
	C.core->name = "C";
	D.core->name = "D";
	if ( AB .is_positive() ) AB .core->name = "AB";
	else AB .reverse() .core->name = "BA";
	if ( BC .is_positive() ) BC .core->name = "BC";
	else BC .reverse() .core->name = "CB";
	if ( CD .is_positive() ) CD .core->name = "CD";
	else CD .reverse() .core->name = "DC";
	if ( DA .is_positive() ) DA .core->name = "DA";
	else DA .reverse() .core->name = "AD";
	ABCD.core->name = "ABCD";

	{ // just a block of code
	std::cout << "meshes above A" << std::endl;
	Cell::Positive::Vertex * A_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	auto segs = A_p->segments;
	for ( auto it = segs.begin(); it != segs.end(); it++ )
	{	std::cout << "bdry_of_" << it->first->get_name() << " " << it->second << std::endl;  }
	auto meshes = A_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above B" << std::endl;
	Cell::Positive::Vertex * B_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	auto segs = B_p->segments;
	for ( auto it = segs.begin(); it != segs.end(); it++ )
	{	std::cout << "bdry_of_" << it->first->get_name() << " " << it->second << std::endl;  }
	auto meshes = B_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above C" << std::endl;
	Cell::Positive::Vertex * C_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( C.core );
	auto segs = C_p->segments;
	for ( auto it = segs.begin(); it != segs.end(); it++ )
	{	std::cout << "bdry_of_" << it->first->get_name() << " " << it->second << std::endl;  }
	auto meshes = C_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above AB" << std::endl;
	Cell::Positive::NotVertex * AB_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( AB.core );
	auto meshes_same_dim = AB_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = AB_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above BC" << std::endl;
	Cell::Positive::NotVertex * BC_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( BC.core );
	auto meshes_same_dim = BC_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = BC_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above CD" << std::endl;
	Cell::Positive::NotVertex * CD_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( CD.core );
	auto meshes_same_dim = CD_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = CD_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} // just a block of code

	AB .cut_from_bdry_of ( ABCD );
	BC .cut_from_bdry_of ( ABCD );
	
	std::cout << "*********************************************************" << std::endl;
	std::cout << "removed two sides of square" << std::endl;
	{ // just a block of code
	std::cout << "meshes above A" << std::endl;
	Cell::Positive::Vertex * A_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	auto segs = A_p->segments;
	for ( auto it = segs.begin(); it != segs.end(); it++ )
	{	std::cout << "bdry_of_" << it->first->get_name() << " " << it->second << std::endl;  }
	auto meshes = A_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above B" << std::endl;
	Cell::Positive::Vertex * B_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	auto segs = B_p->segments;
	for ( auto it = segs.begin(); it != segs.end(); it++ )
	{	std::cout << "bdry_of_" << it->first->get_name() << " " << it->second << std::endl;  }
	auto meshes = B_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above C" << std::endl;
	Cell::Positive::Vertex * C_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( C.core );
	auto segs = C_p->segments;
	for ( auto it = segs.begin(); it != segs.end(); it++ )
	{	std::cout << "bdry_of_" << it->first->get_name() << " " << it->second << std::endl;  }
	auto meshes = C_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above AB" << std::endl;
	Cell::Positive::NotVertex * AB_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( AB.core );
	auto meshes_same_dim = AB_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = AB_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above BC" << std::endl;
	Cell::Positive::NotVertex * BC_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( BC.core );
	auto meshes_same_dim = BC_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = BC_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above CD" << std::endl;
	Cell::Positive::NotVertex * CD_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( CD.core );
	auto meshes_same_dim = CD_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = CD_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} // just a block of code

	Cell AC ( tag::segment, A.reverse(), C );
	AC .core->name = "AC";
	AC .glue_on_bdry_of ( ABCD );
	ABCD .core->name = "ACD";
	
	std::cout << "*********************************************************" << std::endl;
	std::cout << "created diagonal AC" << std::endl;
	{ // just a block of code
	std::cout << "meshes above A" << std::endl;
	Cell::Positive::Vertex * A_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	auto segs = A_p->segments;
	for ( auto it = segs.begin(); it != segs.end(); it++ )
	{	std::cout << "bdry_of_" << it->first->get_name() << " " << it->second << std::endl;  }
	auto meshes = A_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above B" << std::endl;
	Cell::Positive::Vertex * B_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	auto segs = B_p->segments;
	for ( auto it = segs.begin(); it != segs.end(); it++ )
	{	std::cout << "bdry_of_" << it->first->get_name() << " " << it->second << std::endl;  }
	auto meshes = B_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above C" << std::endl;
	Cell::Positive::Vertex * C_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( C.core );
	auto segs = C_p->segments;
	for ( auto it = segs.begin(); it != segs.end(); it++ )
	{	std::cout << "bdry_of_" << it->first->get_name() << " " << it->second << std::endl;  }
	auto meshes = C_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above AB" << std::endl;
	Cell::Positive::NotVertex * AB_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( AB.core );
	auto meshes_same_dim = AB_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = AB_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above BC" << std::endl;
	Cell::Positive::NotVertex * BC_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( BC.core );
	auto meshes_same_dim = BC_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = BC_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above CD" << std::endl;
	Cell::Positive::NotVertex * CD_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( CD.core );
	auto meshes_same_dim = CD_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = CD_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} // just a block of code

	Cell ABC ( tag::triangle, AB, BC, AC.reverse() );
	ABC.core->name = "ABC";
	
	std::cout << "*********************************************************" << std::endl;
	std::cout << "created triangle ABC" << std::endl;
	{ // just a block of code
	std::cout << "meshes above A" << std::endl;
	Cell::Positive::Vertex * A_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	auto segs = A_p->segments;
	for ( auto it = segs.begin(); it != segs.end(); it++ )
	{	std::cout << "bdry_of_" << it->first->get_name() << " " << it->second << std::endl;  }
	auto meshes = A_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above B" << std::endl;
	Cell::Positive::Vertex * B_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	auto segs = B_p->segments;
	for ( auto it = segs.begin(); it != segs.end(); it++ )
	{	std::cout << "bdry_of_" << it->first->get_name() << " " << it->second << std::endl;  }
	auto meshes = B_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above C" << std::endl;
	Cell::Positive::Vertex * C_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( C.core );
	auto segs = C_p->segments;
	for ( auto it = segs.begin(); it != segs.end(); it++ )
	{	std::cout << "bdry_of_" << it->first->get_name() << " " << it->second << std::endl;  }
	auto meshes = C_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above AB" << std::endl;
	Cell::Positive::NotVertex * AB_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( AB.core );
	auto meshes_same_dim = AB_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = AB_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above BC" << std::endl;
	Cell::Positive::NotVertex * BC_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( BC.core );
	auto meshes_same_dim = BC_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = BC_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above CD" << std::endl;
	Cell::Positive::NotVertex * CD_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( CD.core );
	auto meshes_same_dim = CD_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = CD_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} // just a block of code

	ABC .add_to_mesh ( rect_mesh );

	std::cout << "*********************************************************" << std::endl;
	std::cout << "included triangle in big mesh" << std::endl;
	{ // just a block of code
	std::cout << "meshes above A" << std::endl;
	Cell::Positive::Vertex * A_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	auto segs = A_p->segments;
	for ( auto it = segs.begin(); it != segs.end(); it++ )
	{	std::cout << "bdry_of_" << it->first->get_name() << " " << it->second << std::endl;  }
	auto meshes = A_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above B" << std::endl;
	Cell::Positive::Vertex * B_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	auto segs = B_p->segments;
	for ( auto it = segs.begin(); it != segs.end(); it++ )
	{	std::cout << "bdry_of_" << it->first->get_name() << " " << it->second << std::endl;  }
	auto meshes = B_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above C" << std::endl;
	Cell::Positive::Vertex * C_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( C.core );
	auto segs = C_p->segments;
	for ( auto it = segs.begin(); it != segs.end(); it++ )
	{	std::cout << "bdry_of_" << it->first->get_name() << " " << it->second << std::endl;  }
	auto meshes = C_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above AB" << std::endl;
	Cell::Positive::NotVertex * AB_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( AB.core );
	auto meshes_same_dim = AB_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = AB_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above BC" << std::endl;
	Cell::Positive::NotVertex * BC_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( BC.core );
	auto meshes_same_dim = BC_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = BC_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above CD" << std::endl;
	Cell::Positive::NotVertex * CD_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( CD.core );
	auto meshes_same_dim = CD_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = CD_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} // just a block of code

	C .cut_from_bdry_of ( AC );

	std::cout << "*********************************************************" << std::endl;
	std::cout << "C .cut_from_bdry_of ( AC );" << std::endl;
	{ // just a block of code
	std::cout << "meshes above A" << std::endl;
	Cell::Positive::Vertex * A_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	auto segs = A_p->segments;
	for ( auto it = segs.begin(); it != segs.end(); it++ )
	{	std::cout << "bdry_of_" << it->first->get_name() << " " << it->second << std::endl;  }
	auto meshes = A_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above B" << std::endl;
	Cell::Positive::Vertex * B_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	auto segs = B_p->segments;
	for ( auto it = segs.begin(); it != segs.end(); it++ )
	{	std::cout << "bdry_of_" << it->first->get_name() << " " << it->second << std::endl;  }
	auto meshes = B_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above C" << std::endl;
	Cell::Positive::Vertex * C_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( C.core );
	auto segs = C_p->segments;
	for ( auto it = segs.begin(); it != segs.end(); it++ )
	{	std::cout << "bdry_of_" << it->first->get_name() << " " << it->second << std::endl;  }
	auto meshes = C_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above AB" << std::endl;
	Cell::Positive::NotVertex * AB_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( AB.core );
	auto meshes_same_dim = AB_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = AB_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above BC" << std::endl;
	Cell::Positive::NotVertex * BC_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( BC.core );
	auto meshes_same_dim = BC_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = BC_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above CD" << std::endl;
	Cell::Positive::NotVertex * CD_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( CD.core );
	auto meshes_same_dim = CD_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = CD_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " "
									<< itt->second.counter_neg << std::endl;  }  }
	} // just a block of code

	Cell E ( tag::vertex );
	E.core->name = "E";
	
}  // end of  main

