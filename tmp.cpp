

void frontal_construct          // hidden in anonymous namespace
( Mesh & msh, const tag::StartWithInconsistentMesh &,
  const tag::StartAt &, const Cell & start,
  const tag::StopAt &, const Cell & stop,
  const tag::Orientation &, const tag::OrientationChoice & oc )

// 'start' and 'stop' are positive vertices (may be one and the same)

// here the working manifold is not a quotient manifold

{	assert ( start .dim() == 0 );
	assert ( stop  .dim() == 0 );
	assert ( start .is_positive() );
	assert ( stop  .is_positive() );
	assert ( msh .dim() == 1 );

	// here the working manifold is not a quotient manifold
	#ifndef NDEBUG
	{ // just a block of code for hiding 'm'
	Manifold::Quotient * m = dynamic_cast < Manifold::Quotient * > ( Manifold::working .core );
	assert ( m == nullptr );
	} // just a block of code
	#endif

	if ( oc == tag::not_provided )
	{	std::cout << "when starting and stopping points are provided," << std::endl;
		std::cout << "maniFEM needs to know how to choose the orientation of the curve;" << std::endl;
		std::cout << "please specify either tag::orientation or tag::shortest_path" << std::endl;
		exit (1);                                                                                      }
	
	std::vector < double > best_tangent = compute_tangent_vec ( tag::at_point, start );

	if ( oc == tag::geodesic )   // shortest path

	{	assert ( start != stop );
		
		// start walking along the manifold from 'start' in the direction of best_tangent
		// and, simultaneously, in the opposite direction, given by -best_tangent
		std::vector < double > tan1 = best_tangent, tan2 = best_tangent;
		for ( size_t i = 0; i < frontal_nb_of_coords; i++ ) tan2 [i] *= -1.;
		Cell ver1 ( tag::vertex ), ver2 ( tag::vertex );
		for ( size_t i = 0; i < frontal_nb_of_coords; i++ )
		{	Function x = Manifold::working .coordinates() [i];
			x ( ver1 ) = x ( start );  x ( ver2 ) = x ( start );  }
		int winner = 0;  //  will be 1 or -1
		while ( true )
		{	double augm_length = desired_length(ver1) * 1.5,
			       augm_len_sq = augm_length * augm_length;
			for ( size_t i = 0; i < frontal_nb_of_coords; i++ )
			{	Function x = Manifold::working .coordinates() [i];
				x ( temporary_vertex ) = x ( ver1 ) + tan1 [i];  }
			Manifold::working .project ( temporary_vertex );
			for ( size_t i = 0; i < frontal_nb_of_coords; i++ )
			{	Function x = Manifold::working .coordinates() [i];
				tan1[i] = x ( temporary_vertex ) - x ( ver1 );
				x ( ver1 ) = x ( temporary_vertex );            }
			double d = Manifold::working .dist_sq ( ver1, stop );
			if ( d < augm_len_sq )
			{	double prod = 0.;
				for ( size_t i = 0; i < frontal_nb_of_coords; i++ )
				{	Function x = Manifold::working .coordinates() [i];
					prod += tan1[i] * ( x (stop) - x (ver1) );         }
				if ( prod > 0. )  { winner = 1;  break;  }             }
			for ( size_t i = 0; i < frontal_nb_of_coords; i++ )
			{	Function x = Manifold::working .coordinates() [i];
				x ( temporary_vertex ) = x ( ver2 ) + tan2 [i];  }
			Manifold::working .project ( temporary_vertex );
			for ( size_t i = 0; i < frontal_nb_of_coords; i++ )
			{	Function x = Manifold::working .coordinates() [i];
				tan2[i] = x ( temporary_vertex ) - x ( ver2 );
				x ( ver2 ) = x ( temporary_vertex );            }
			d = Manifold::working .dist_sq ( ver2, stop );
			if ( d < augm_len_sq )
			{	double prod = 0.;
				for ( size_t i = 0; i < frontal_nb_of_coords; i++ )
				{	Function x = Manifold::working .coordinates() [i];
					prod += tan2 [i] * ( x (stop) - x (ver2) );          }
				if ( prod > 0. )  { winner = -1;  break;  }           }
		}  // end of  while true

		assert ( ( winner == 1 ) or ( winner == -1 ) );
		for ( size_t i = 0; i < frontal_nb_of_coords; i++ ) best_tangent [i] *= winner;
		frontal_construct ( msh, tag::start_with_inconsistent_mesh,
		                    tag::start_at, start, tag::towards, best_tangent,
		                    tag::stop_at, stop                                    );
		return;                                                                          }

	if ( ( oc == tag::inherent ) or ( oc == tag::random ) )

	{	if ( start == stop )
		{	frontal_construct ( msh, tag::start_with_inconsistent_mesh,
			                    tag::start_at, start, tag::orientation, oc );
			return;                                                           }

		frontal_construct ( msh, tag::start_with_inconsistent_mesh,
		                    tag::start_at, start, tag::towards, best_tangent,
 		                    tag::stop_at, stop                               );

		if ( oc == tag::random ) return;

		for ( size_t i = 0; i < frontal_nb_of_coords; i++ )  best_tangent[i] *= -1.;
		// the number of segments does not count, and we don't know it yet
		Mesh msh2 ( tag::whose_core_is,
		    new Mesh::Connected::OneDim ( tag::with, 1, tag::segments, tag::one_dummy_wrapper ),
	  	  tag::freshly_created, tag::is_positive                                              );
		// the number of segments does not count, and we don't know it yet
		// we compute it after the mesh is built, by counting segments
		// but we count segments using an iterator, and the iterator won't work
		// if this->msh->nb_of_segs == 0, so we set nb_of_segs to 1 (dirty trick)
		// see Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder::reset
		// in iterator.cpp
		frontal_construct ( msh2, tag::start_with_inconsistent_mesh,
		                    tag::start_at, start, tag::towards, best_tangent,
		                    tag::stop_at, stop                               );

		switch_orientation_direct ( msh2 );
		update_info_connected_one_dim ( msh2, stop, start );
		Mesh whole ( tag::join, msh, msh2 );

		if ( not correctly_oriented ( whole, tag::orientation, oc ) )
		{	switch_orientation_direct ( msh2 );  msh = msh2;  }

		return;                                                                                  }

	assert ( oc == tag::intrinsic );
	assert ( false );
	std::cout << "intrinsic orientation does not make sense for one-dimensional manifolds"
	          << std::endl;   // ???
	exit (1);

} // end of  frontal_construct

//------------------------------------------------------------------------------------------------------//


void Manifold::Core::frontal_method  // virtual, overridden by Manifold::Quotient
( Mesh & msh, const tag::StartWithInconsistentMesh &,
  const tag::StartAt &, const Cell & start,
  const tag::StopAt &,  const Cell & stop, const tag::ShortestPath & )

{	frontal_construct < ManifoldNoWinding >  // line ???
	( msh, tag::start_with_inconsistent_mesh,
	  tag::start_at, start, tag::stop_at, stop,
		tag::orientation, tag::geodesic          );  }


void Manifold::Quotient::frontal_method  // virtual from Manifold::Core, here overridden
( Mesh & msh, const tag::StartWithInconsistentMesh &,
  const tag::StartAt &, const Cell & start,
  const tag::StopAt &,  const Cell & stop, const tag::ShortestPath & )

{	assert ( false );  }
	
//-------------------------------------------------------------------------------------------------



