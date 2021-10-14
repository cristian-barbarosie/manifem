

#include "maniFEM.h"

using namespace maniFEM;


	
int main ( )

{	// begin with the usual two-dimensional space
	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0], y = xy[1];

	Cell A ( tag::vertex );  x(A) =  0.  ;  y(A) = -0.35;
	Cell B ( tag::vertex );  x(B) =  0.  ;  y(B) = -1.05;
	Cell C ( tag::vertex );  x(C) =  1.15;  y(C) = -2.15;
	Cell D ( tag::vertex );  x(C) =  1.15;  y(C) = -2.15;

	A.core->name = "A";
	B.core->name = "B";
	C.core->name = "C";
	D.core->name = "D";

	Cell BA ( tag::segment, B.reverse(), A );
	Cell BC ( tag::segment, B.reverse(), C );
	Cell CA ( tag::segment, C.reverse(), A );
	Cell AD ( tag::segment, A.reverse(), D );
	Cell BD ( tag::segment, B.reverse(), D );

	BA.core->name = "BA";
	BC.core->name = "BC";
	CA.core->name = "CA";
	AD.core->name = "AD";
	BD.core->name = "BD";

	Cell BAD ( tag::triangle, BA, AD, BD.reverse() );
	Cell ABC ( tag::triangle, BA.reverse(), BC, CA );
	BAD.core->name = "BAD";
	ABC.core->name = "ABC";

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
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
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
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
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
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above BA" << std::endl;
	Cell::Positive::NotVertex * BA_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( BA.core );
	auto meshes_same_dim = BA_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = BA_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above CA" << std::endl;
	Cell::Positive::NotVertex * CA_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( CA.core );
	auto meshes_same_dim = CA_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = CA_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
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
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
	} // just a block of code

	Mesh msh ( tag::whose_core_is,
     new Mesh::Fuzzy ( tag::of_dimension, 3, tag::minus_one,
											 tag::one_dummy_wrapper ),
						 tag::freshly_created, tag::is_positive           );
	msh.core->name = "msh";

	ABC .add_to_mesh ( msh );
	
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
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
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
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
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
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above BA" << std::endl;
	Cell::Positive::NotVertex * BA_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( BA.core );
	auto meshes_same_dim = BA_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = BA_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above CA" << std::endl;
	Cell::Positive::NotVertex * CA_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( CA.core );
	auto meshes_same_dim = CA_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = CA_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
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
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
	} // just a block of code

	BAD .add_to_mesh ( msh );
	
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
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
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
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
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
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above BA" << std::endl;
	Cell::Positive::NotVertex * BA_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( BA.core );
	auto meshes_same_dim = BA_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = BA_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above CA" << std::endl;
	Cell::Positive::NotVertex * CA_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( CA.core );
	auto meshes_same_dim = CA_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = CA_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
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
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
	} // just a block of code

  BA .reverse() .cut_from_bdry_of ( ABC );

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
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
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
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
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
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above BA" << std::endl;
	Cell::Positive::NotVertex * BA_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( BA.core );
	auto meshes_same_dim = BA_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = BA_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above CA" << std::endl;
	Cell::Positive::NotVertex * CA_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( CA.core );
	auto meshes_same_dim = CA_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = CA_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
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
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
	} // just a block of code

  BA .reverse() .glue_on_bdry_of ( ABC );

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
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
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
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
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
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above BA" << std::endl;
	Cell::Positive::NotVertex * BA_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( BA.core );
	auto meshes_same_dim = BA_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = BA_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
	} { // just a block of code
	std::cout << "meshes above CA" << std::endl;
	Cell::Positive::NotVertex * CA_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( CA.core );
	auto meshes_same_dim = CA_p->meshes_same_dim;
	for ( auto it = meshes_same_dim.begin(); it != meshes_same_dim.end(); it++ )
	{	std::cout << it->first->get_name() << " " << it->second.sign << std::endl;  }
	auto meshes = CA_p->meshes;
	for ( auto it = meshes.begin(); it != meshes.end(); it++ )
	{	auto m = *it;
		for ( auto itt = m.begin(); itt != m.end(); itt++ )
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
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
			{	std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << std::endl;  }  }
	} // just a block of code

	std::cout << "reached end of main" << std::endl;
}

