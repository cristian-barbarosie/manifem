

// see global.cpp fold_common_1

#include "maniFEM.h"
using namespace maniFEM;


int main ( )

{	Mesh m ( tag::fuzzy, tag::of_dim, 2 );

	Cell A ( tag::vertex);
	Cell B ( tag::vertex);
	Cell C ( tag::vertex);

	Cell AB ( tag::segment, A .reverse(), B );  AB .core->name = "AB";
	Cell BC ( tag::segment, B .reverse(), C );  BC .core->name = "BC";
	Cell CA ( tag::segment, C .reverse(), A );  CA .core->name = "CA";

	Cell::Positive::HighDim * new_cll_ptr = new Cell::Positive::HighDim
		( tag::whose_boundary_is,
			Mesh ( tag::whose_core_is,
		         new Mesh::Connected::OneDim ( tag::with, 3,
		           tag::segments, tag::one_dummy_wrapper       ),
         tag::freshly_created                                 ),
			tag::one_dummy_wrapper                                     );
	Cell ABC ( tag::whose_core_is, new_cll_ptr, tag::freshly_created );
	Mesh ABC_bdry = ABC .boundary();

	AB .core->add_to_mesh ( ABC_bdry .core, tag::do_not_bother );
	BC .core->add_to_mesh ( ABC_bdry .core, tag::do_not_bother );
	CA .core->add_to_mesh ( ABC_bdry .core, tag::do_not_bother );

	ABC_bdry .closed_loop ( A );

	ABC .add_to_mesh ( m );

	Cell D ( tag::vertex );

	Cell AD ( tag::segment, A .reverse(), D );
	Cell BD ( tag::segment, B .reverse(), D );

	new_cll_ptr = new Cell::Positive::HighDim
		( tag::whose_boundary_is,
			Mesh ( tag::whose_core_is,
		         new Mesh::Connected::OneDim ( tag::with, 3,
		           tag::segments, tag::one_dummy_wrapper       ),
         tag::freshly_created                                 ),
			tag::one_dummy_wrapper                                     );
	Cell BAD ( tag::whose_core_is, new_cll_ptr, tag::freshly_created );
	Mesh BAD_bdry = BAD .boundary();

	AB .reverse() .core->add_to_mesh ( BAD_bdry .core, tag::do_not_bother );
	AD .core->add_to_mesh ( BAD_bdry .core, tag::do_not_bother );
	BD .reverse() .core->add_to_mesh ( BAD_bdry .core, tag::do_not_bother );

	BAD_bdry .closed_loop ( A );

	BAD .add_to_mesh ( m );

	CellIterator it = ABC_bdry .iterator
		( tag::over_cells, tag::of_max_dim, tag::force_positive );
	for ( it .reset(); it .in_range(); it ++ )
		std::cout << ( * it ) .core->get_name() << std::endl;	
}
