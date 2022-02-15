
// mesh.cpp 2022.02.15

//   This file is part of maniFEM, a C++ library for meshes and finite elements on manifolds.

//   Copyright 2019 - 2022 Cristian Barbarosie cristian.barbarosie@gmail.com

//   http://manifem.rd.ciencias.ulisboa.pt/
//   https://github.com/cristian-barbarosie/manifem

//   ManiFEM is free software: you can redistribute it and/or modify it
//   under the terms of the GNU Lesser General Public License as published
//   by the Free Software Foundation, either version 3 of the License
//   or (at your option) any later version.

//   ManiFEM is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//   See the GNU Lesser General Public License for more details.

//   You should have received a copy of the GNU Lesser General Public License
//   along with maniFEM.  If not, see <https://www.gnu.org/licenses/>.

#include <forward_list>

#include "mesh.h"
#include "iterator.h"
#include "manifold.h"

using namespace maniFEM;


size_t Mesh::maximum_dimension_plus_one { 4 };  // static data member

// we keep here the topological dimension of the largest mesh we intend to build
// '4' means three-dimensional meshes (cubes, tetrahedra, etc)
// '3' means two-dimensional meshes, including surfaces in R^3
// '2' would be for just polygonal lines
// '1' doesn't make much sense - just points ?
// '0' what ?? not even points ?! then perhaps you dont't need maniFEM at all

// see method Mesh::set_max_dim and paragraph 11.7 in the manual

// static data members :

std::vector < size_t > Cell::Positive::double_heap_size ( Mesh::maximum_dimension_plus_one, 0 );
std::vector < size_t > Cell::Negative::double_heap_size ( Mesh::maximum_dimension_plus_one, 0 );
std::vector < size_t > Cell::Positive::size_t_heap_size ( Mesh::maximum_dimension_plus_one, 0 );
std::vector < size_t > Cell::Negative::size_t_heap_size ( Mesh::maximum_dimension_plus_one, 0 );
std::vector < size_t > Cell::Positive::short_int_heap_size ( Mesh::maximum_dimension_plus_one, 0 );
std::vector < size_t > Cell::Negative::short_int_heap_size ( Mesh::maximum_dimension_plus_one, 0 );
std::vector < std::vector < void(*)(Cell::Core*,void*) > >
   Cell::init_pos_cell ( Mesh::maximum_dimension_plus_one );
std::vector < std::vector < void(*)(Cell::Core*,void*) > >
   Cell::init_neg_cell ( Mesh::maximum_dimension_plus_one );
std::vector < std::vector < void* > > Cell::data_for_init_pos ( Mesh::maximum_dimension_plus_one );
std::vector < std::vector < void* > > Cell::data_for_init_neg ( Mesh::maximum_dimension_plus_one );

const double tag::Util::one_third = 1. / 3.,
             tag::Util::minus_one_third = - tag::Util::minus_one_third,
             tag::Util::one_sixth = 1. / 6.,
             tag::Util::minus_one_sixth = - tag::Util::one_sixth,
             tag::Util::two_thirds = 2. / 3.,
             tag::Util::minus_two_thirds = - tag::Util::two_thirds,
             tag::Util::sqrt_2 = std::sqrt (2.),
             tag::Util::sqrt_half = 1./ tag::Util::sqrt_2,
             tag::Util::sqrt_3 = std::sqrt (3.),
             tag::Util::sqrt_third = 1./ tag::Util::sqrt_3,
             tag::Util::sqrt_two_thirds = std::sqrt ( tag::Util::two_thirds ),
             tag::Util::sqrt_three_quarters = std::sqrt (0.75);
const std::vector < std::vector < double > >
	tag::Util::ortho_basis_in_R2 { { 1., 0. }, { 0., 1. } },
	tag::Util::ortho_basis_in_R3 { { 1., 0., 0. }, { 0., 1., 0 }, { 0., 0., 1. } },
	tag::Util::four_directions_in_R2 { { 1., 0. }, { 0., 1. }, { -1., 0. }, { 0., -1. } },
	tag::Util::six_directions_in_R3 { { 1., 0., 0. },  { 0., 1., 0. }, { -1., 0., 0. },
	                                  { 0., -1., 0. }, { 0., 0., 1. }, { 0, 0., -1. }  },
	tag::Util::eight_directions_in_R2 { { 1., 0. },  {  tag::Util::sqrt_half,  tag::Util::sqrt_half },
	                                    { 0., 1. },  { -tag::Util::sqrt_half,  tag::Util::sqrt_half },
	                                    { -1., 0. }, { -tag::Util::sqrt_half, -tag::Util::sqrt_half },
	                                    { 0., -1. }, {  tag::Util::sqrt_half, -tag::Util::sqrt_half } },
	tag::Util::twentysix_directions_in_R3
				{ { 1., 0., 0. }, { -1., 0., 0. }, { 0., 1., 0. }, 
				  { 0., -1., 0. }, { 0., 0., 1. }, { 0., 0., -1. },
				  { tag::Util::sqrt_half, tag::Util::sqrt_half, 0. },
			  	{ tag::Util::sqrt_half, -tag::Util::sqrt_half, 0. },
				  { -tag::Util::sqrt_half, tag::Util::sqrt_half, 0. },
				  { -tag::Util::sqrt_half, -tag::Util::sqrt_half, 0. },
				  { tag::Util::sqrt_half, 0., tag::Util::sqrt_half },
				  { tag::Util::sqrt_half, 0., -tag::Util::sqrt_half },
			  	{ -tag::Util::sqrt_half, 0., tag::Util::sqrt_half },
				  { -tag::Util::sqrt_half, 0., -tag::Util::sqrt_half },
				  { 0., tag::Util::sqrt_half, tag::Util::sqrt_half },
				  { 0., tag::Util::sqrt_half, -tag::Util::sqrt_half },
				  { 0., -tag::Util::sqrt_half, tag::Util::sqrt_half },
			  	{ 0., -tag::Util::sqrt_half, -tag::Util::sqrt_half },
					{  tag::Util::sqrt_third,  tag::Util::sqrt_third,  tag::Util::sqrt_third },
					{  tag::Util::sqrt_third,  tag::Util::sqrt_third, -tag::Util::sqrt_third },
					{  tag::Util::sqrt_third, -tag::Util::sqrt_third,  tag::Util::sqrt_third },
					{  tag::Util::sqrt_third, -tag::Util::sqrt_third, -tag::Util::sqrt_third },
					{ -tag::Util::sqrt_third,  tag::Util::sqrt_third,  tag::Util::sqrt_third },
					{ -tag::Util::sqrt_third,  tag::Util::sqrt_third, -tag::Util::sqrt_third },
					{ -tag::Util::sqrt_third, -tag::Util::sqrt_third,  tag::Util::sqrt_third },
					{ -tag::Util::sqrt_third, -tag::Util::sqrt_third, -tag::Util::sqrt_third } };


// int Cell::counter { 0 };

//-----------------------------------------------------------------------------//


bool tag::Util::Core::default_dispose_query ( tag::Util::Core::DelegateDispose * that )
// static method
{	tag::Util::Core * that_one = static_cast < tag::Util::Core* > ( that );
	return that_one->dispose_query ();  // calls tag::Util::Core::dispose_query
}

bool tag::Util::Core::dispose_query_cell_with_reverse ( tag::Util::Core::DelegateDispose * that )
// static method
// a cell may still have one wrapper kept by its reverse which is a positive cell
#ifdef MANIFEM_COLLECT_CM	
{	tag::Util::Core * that_one = static_cast < tag::Util::Core* > ( that );
	bool res = that_one->dispose_query ();
	// calls tag::Util::Core::dispose_query which decrements nb_of_wrappers
	assert ( not res );  // 'that' cell is kept alive by its reverse
	if ( that->nb_of_wrappers == 1 )
	{	Cell::Core * that_cell = tag::Util::assert_cast
			< tag::Util::Core::DelegateDispose*, Cell::Core* > ( that );
		assert ( that_cell->reverse_attr.exists() );
		if ( that_cell->reverse_attr.core->nb_of_wrappers == 1 )
		// these two cells keep each other alive, both must be killed
		// it suffices to ask the calling code to kill 'that'
		// its reverse will be destroyed, too, in the process
		{	that_cell->reverse_attr.core->reverse_attr.core = nullptr;
			that_cell->reverse_attr.core->dispose_query_p =
				& tag::Util::Core::default_dispose_query;
			return true;                                                 }  }
	return false;                                                          }
#else  // no MANIFEM_COLLECT_CM
{	assert ( false );  // we return false just to avoid compilation errors
	return false;     }
#endif  // MANIFEM_COLLECT_CM	

//-----------------------------------------------------------------------------//


namespace { // anonymous namespace, mimics static linkage
	
// this should work well even if 'msh' is the boundary of some cell
// if the cell is about to be destroyed, then it does not belong to any
// higher-dimensional mesh
	
inline void shred_mesh ( Mesh::Core * msh )
	
// this is called before the destructor, so the object is still entire and functional
// however, we should only do this for terminal classes
// (no other class should inherit from here)
// but this is not so strict
// for instance, we call 'shred_mesh' on an ex-STSI mesh, partially destroyed,
// now acting like Fuzzy mesh
// it's OK since this kind of iterator is defined in Fuzzy, not in STSI
// mas como raio iremos percorrer a malha se ela tiver efectivamente auto-interseccoes ?
// conseguimos porque temos as listas do Fuzzy ainda ...
// mas vamos ter problemas com o cell_behind_within
// o melhor sera' certificar-nos que, no fim da sua vida,
// uma malha STSI ou esta vazia ou ja nao tem autointerseccoes
// but we should not call it from the destructor of a STSI mesh
// (if we do so, it will be executed twice on the same mesh)

// we build a list of cells
// it is not obvious whether a cell iterator would behave well
// on a mesh which is being torn down to nothing ...
	
{	std::forward_list < Cell > l;
	// these Cell wrappers are useful, they prevent premature destruction of cells
	{ // just a block of code for hiding 'it'
	Mesh::Iterator::Core * it = msh->iterator
		( tag::over_cells_of_max_dim, tag::as_they_are, tag::this_mesh_is_positive );
	// iterates over cells of 'msh'
	for ( it->reset(); it->in_range(); it->advance() ) l .push_front ( it->deref() );
	} { // just a block of code for hiding 'it'
	std::forward_list<Cell>::iterator it;
	for ( it = l .begin(); it != l .end(); it++ )
		(*it) .core->remove_from_mesh ( msh );	    }                                    }

}  // anonymous namespace


Mesh::Connected::OneDim::~OneDim ()  // virtual
{	shred_mesh ( this );  }


Mesh::Fuzzy::~Fuzzy ()  // virtual
{	shred_mesh ( this );  }


// Mesh::STSI::~STSI ()  // virtual
// {	// do not shred_mesh ! it will be called from ~Fuzzy
// 	// just check that there are no self-intersections
// 	// because if there are self-intersections,
// 	// the destructor of Mesh::Fuzzy will get in trouble with cell_behind_within
// }


//-----------------------------------------------------------------------------//


bool Cell::Positive::is_positive ( ) const  // virtual from Cell::Core
{	return true;  }

bool Cell::Negative::is_positive ( ) const  // virtual from Cell::Core
{	return false;  }

bool tag::Util::return_true ( ) // static
{	return true;  }

bool tag::Util::return_false ( ) // static
{	return false;  }

Cell Cell::Positive::get_positive ( )  // virtual from Cell::Core
{	return Cell ( tag::whose_core_is, this, tag::previously_existing, tag::surely_not_null );  }

Cell Cell::Negative::get_positive ( )  // virtual from Cell::Core
{	assert ( this->reverse_attr.exists() );
	assert ( this->reverse_attr.is_positive() );
	return this->reverse_attr;                   }


Cell::Negative * Cell::Positive::Vertex::build_reverse ( const tag::OneDummyWrapper & )
// virtual from Cell::Positive
{	assert ( not this->reverse_attr.exists() );
	#ifdef MANIFEM_COLLECT_CM	
	this->dispose_query_p = & tag::Util::Core::dispose_query_cell_with_reverse;
	#endif  // MANIFEM_COLLECT_CM	
	return new Cell::Negative::Vertex ( tag::reverse_of, this, tag::one_dummy_wrapper );  }

Cell::Negative * Cell::Positive::Segment::build_reverse ( const tag::OneDummyWrapper & )
// virtual from Cell::Positive
{	assert ( not this->reverse_attr.exists() );
	#ifdef MANIFEM_COLLECT_CM	
	this->dispose_query_p = & tag::Util::Core::dispose_query_cell_with_reverse;
	#endif  // MANIFEM_COLLECT_CM
	return new Cell::Negative::Segment ( tag::reverse_of, this, tag::one_dummy_wrapper );  }

Cell::Negative * Cell::Positive::HighDim::build_reverse ( const tag::OneDummyWrapper & )
// virtual from Cell::Positive
{	assert ( not this->reverse_attr.exists() );
	#ifdef MANIFEM_COLLECT_CM	
	this->dispose_query_p = & tag::Util::Core::dispose_query_cell_with_reverse;
	#endif  // MANIFEM_COLLECT_CM	
	return new Cell::Negative::HighDim ( tag::reverse_of, this, tag::one_dummy_wrapper );  }

//-----------------------------------------------------------------------------//


size_t Cell::Positive::Vertex::get_dim ( ) const  // virtual from Cell::Core
{	return 0;  }

size_t Cell::Negative::Vertex::get_dim ( ) const  // virtual from Cell::Core
{	return 0;  }

size_t Cell::Positive::Segment::get_dim ( ) const  // virtual from Cell::Core
{	return 1;  }

size_t Cell::Negative::Segment::get_dim ( ) const  // virtual from Cell::Core
{	return 1;  }

size_t Cell::Positive::HighDim::get_dim ( ) const  // virtual from Cell::Core
{	return this->boundary_attr .core->get_dim_plus_one();  }

size_t Cell::Negative::HighDim::get_dim ( ) const  // virtual from Cell::Core
{	assert ( this->reverse_attr .exists() );
	return this->reverse_attr .dim();        }


size_t Mesh::ZeroDim::get_dim_plus_one ( )  // virtual from Mesh::Core
{	return 1;  }

size_t Mesh::Connected::OneDim::get_dim_plus_one ( )  // virtual from Mesh::Core
{	return 2;  }

size_t Mesh::Connected::HighDim::get_dim_plus_one ( )  // virtual from Mesh::Core
{	return this->nb_of_cells .size();  }

size_t Mesh::MultiplyConnected::OneDim::get_dim_plus_one ( )  // virtual from Mesh::Core
{	return 2;  }

size_t Mesh::MultiplyConnected::HighDim::get_dim_plus_one ( )  // virtual from Mesh::Core
{	return this->nb_of_cells .size();  }

size_t Mesh::Fuzzy::get_dim_plus_one ( )  // virtual from Mesh::Core
{	return this->cells .size();  }

	
size_t Mesh::ZeroDim::number_of ( const tag::Vertices & )
// virtual from Mesh::Core
{	return 2;  }
// what if the segment is incomplete ?
	
size_t Mesh::ZeroDim::number_of ( const tag::Segments & )
// virtual from Mesh::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "zero-dimensional meshes have have no segments" << std::endl;
	exit ( 1 );                                                                                     }
	
size_t Mesh::ZeroDim::number_of ( const tag::CellsOfDim &, const size_t d )
// virtual from Mesh::Core
{	assert ( d == 0 );  return 2;  }
// what if the segment is incomplete ?
	
size_t Mesh::ZeroDim::number_of ( const tag::CellsOfMaxDim & )
// virtual from Mesh::Core
{	return 2;  }
// what if the segment is incomplete ?
	
size_t Mesh::Connected::OneDim::number_of ( const tag::Vertices & )
// virtual from Mesh::Core
{	if ( this->first_ver .reverse() == this->last_ver )  // closed loop
		return this->nb_of_segs;
	return this->nb_of_segs + 1;                          }
	
size_t Mesh::Connected::OneDim::number_of ( const tag::Segments & )
// virtual from Mesh::Core
{	return this->nb_of_segs;  }
	
size_t Mesh::Connected::OneDim::number_of ( const tag::CellsOfDim &, const size_t d )
// virtual from Mesh::Core
{	if ( d == 1 ) return this->nb_of_segs;
	assert ( d == 0 );
	if ( this->first_ver .reverse() == this->last_ver )  // closed loop
		return this->nb_of_segs;
	return this->nb_of_segs + 1;              }
	
size_t Mesh::Connected::OneDim::number_of ( const tag::CellsOfMaxDim & )
// virtual from Mesh::Core
{	return this->nb_of_segs;  }
	
size_t Mesh::Connected::HighDim::number_of ( const tag::Vertices & )
// virtual from Mesh::Core
{	assert ( this->nb_of_cells .size() > 0 );
	return this->nb_of_cells [0];             }
	
size_t Mesh::Connected::HighDim::number_of ( const tag::Segments & )
// virtual from Mesh::Core
{	assert ( this->nb_of_cells .size() > 1 );
	return this->nb_of_cells [1];             }
	
size_t Mesh::Connected::HighDim::number_of ( const tag::CellsOfDim &, const size_t d )
// virtual from Mesh::Core
{	assert ( d < this->get_dim_plus_one() );
	assert ( this->nb_of_cells .size() > d );
	return this->nb_of_cells [d];             }
	
size_t Mesh::Connected::HighDim::number_of ( const tag::CellsOfMaxDim & )
// virtual from Mesh::Core
{	return this->nb_of_cells .back();             }
	
// size_t Mesh::Connected::number_of ( const tag::CellsOfDim &, const size_t d )
// virtual from Mesh::Core
// {	assert ( d < this->get_dim_plus_one() );
// 	assert ( this->nb_of_cells .size() > d );
//	return this->nb_of_cells [d];             }
	
size_t Mesh::Fuzzy::number_of ( const tag::Vertices & )
// virtual from Mesh::Core
{	assert ( this->cells .size() > 0 );
	return this->cells [0] .size();       }
	
size_t Mesh::Fuzzy::number_of ( const tag::Segments & )
// virtual from Mesh::Core
{	assert ( 1 < this->get_dim_plus_one() );
	assert ( this->cells .size() > 1 );
	return this->cells [1] .size();            }
	
size_t Mesh::Fuzzy::number_of ( const tag::CellsOfDim &, const size_t d )
// virtual from Mesh::Core
{	assert ( d < this->get_dim_plus_one() );
	assert ( this->cells .size() > d );
	return this->cells [d] .size();            }
	
size_t Mesh::Fuzzy::number_of ( const tag::CellsOfMaxDim & )
// virtual from Mesh::Core
{	return this->cells .back() .size();  }
	

Cell Mesh::Core::first_vertex ( )  // virtual
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "only one-dimensional meshes have first vertex" << std::endl;
	exit ( 1 );                                                                                     }

Cell Mesh::Core::last_vertex ( )  // virtual
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "only one-dimensional meshes have last vertex" << std::endl;
	exit ( 1 );                                                                                     }

Cell Mesh::Core::first_segment ( )  // virtual
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "only one-dimensional meshes have first segment" << std::endl;
	exit ( 1 );                                                                                     }

Cell Mesh::Core::last_segment ( )  // virtual
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "only one-dimensional meshes have last segment" << std::endl;
	exit ( 1 );                                                                                     }

Cell Mesh::Connected::OneDim::first_vertex ( )
// virtual from Mesh::Core, here overriden
// returns a negative vertex
{	assert ( this->first_ver .exists() );
	assert ( not this->first_ver .is_positive() );
	assert ( this->last_ver .exists() );
	assert ( this->last_ver .is_positive() );
	return this->first_ver;                        }

Cell Mesh::Connected::OneDim::last_vertex ( )
// virtual from Mesh::Core, here overriden
{	assert ( this->first_ver .exists() );
	assert ( not this->first_ver .is_positive() );
	assert ( this->last_ver .exists() );
	assert ( this->last_ver .is_positive() );
	return this->last_ver;                         }

Cell Mesh::Connected::OneDim::first_segment ( )
// virtual from Mesh::Core, here overriden
{	assert ( this->first_ver.exists() );
	assert ( not this->first_ver .is_positive() );
	assert ( this->last_ver .exists() );
	assert ( this->last_ver .is_positive() );
	Mesh this_mesh ( tag::whose_core_is, this, tag::previously_existing, tag::is_positive );
	return this_mesh .cell_in_front_of ( this->first_ver );                                   }

Cell Mesh::Connected::OneDim::last_segment ( )
// virtual from Mesh::Core, here overriden
{	assert ( this->first_ver .exists() );
	assert ( not this->first_ver .is_positive() );
	assert ( this->last_ver .exists() );
	assert ( this->last_ver .is_positive() );
	Mesh this_mesh ( tag::whose_core_is, this, tag::previously_existing, tag::is_positive );
	return this_mesh .cell_behind ( this->last_ver );                                         }

//-----------------------------------------------------------------------------//


Mesh Cell::Positive::Vertex::boundary ( )  // virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "vertices have no boundary" << std::endl;
	exit ( 1 );                                             }

Mesh Cell::Negative::Vertex::boundary ( )  // virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "vertices have no boundary" << std::endl;
	exit ( 1 );                                             }

Mesh Cell::Positive::Segment::boundary ( )  // virtual from Cell::Core
{	return Mesh ( tag::whose_core_is,
	              new Mesh::ZeroDim ( tag::boundary_of, tag::positive, tag::segment,
                                    this, tag::one_dummy_wrapper                   ),
	              tag::freshly_created                                                  );  }

Mesh Cell::Negative::Segment::boundary ( )  // virtual from Cell::Core
{	assert ( this->reverse_attr.exists() );
	Cell::Positive::Segment * rev_seg = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Segment* > ( this->reverse_attr.core );
	// return reversed mesh, faces surely exist
	return Mesh ( tag::whose_core_is,
	              new Mesh::ZeroDim ( tag::boundary_of, tag::positive, tag::segment,
                                    rev_seg, tag::one_dummy_wrapper                ),
	              tag::freshly_created, tag::is_negative, tag::do_not_build_cells      );  }

Mesh Cell::Positive::HighDim::boundary ( )  // virtual from Cell::Core
{	return this->boundary_attr;  }

Mesh Cell::Negative::HighDim::boundary ( )  // virtual from Cell::Core
{	Cell pos_cll = this->reverse_attr;
	assert ( pos_cll .exists() );
	assert ( pos_cll .is_positive() );
	return ( pos_cll .boundary() .reverse() );  }

//-----------------------------------------------------------------------------//


// the two methods below are only relevant for STSI meshes
// so we forbid execution in Mesh::Core and then override them in Mesh::STSI

Cell::Core * Mesh::Core::cell_in_front_of  // virtual
(	const Cell::Core * face_p, const tag::SeenFrom &, const Cell neighbour,
	const tag::SurelyExists & se                                            ) const
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "cell_in_front as seen_from : use on STSI meshes only" << std::endl;
	exit ( 1 );                                                                                     }

Cell::Core * Mesh::Core::cell_in_front_of  // virtual
(	const Cell::Core * face_p, const tag::SeenFrom &, const Cell neighbour,
	const tag::MayNotExist &                                                ) const
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "cell_in_front as seen_from : use on STSI meshes only" << std::endl;
	exit ( 1 );                                                                                     }

Cell::Core * Mesh::Core::cell_behind  // virtual
(	const Cell::Core * face_p, const tag::SeenFrom &, const Cell neighbour,
	const tag::SurelyExists & se                                            ) const
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "cell_behind as seen_from : use on STSI meshes only" << std::endl;
	exit ( 1 );                                                                                     }

Cell::Core * Mesh::Core::cell_behind  // virtual
(	const Cell::Core * face_p, const tag::SeenFrom &, const Cell neighbour,
	const tag::MayNotExist &                                                ) const
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "cell_behind as seen_from : use on STSI meshes only" << std::endl;
	exit ( 1 );                                                                                     }

/*
// the four methods below are only relevant for STSI meshes

Cell::Core * Mesh::STSI::cell_in_front_of
// virtual from Mesh::Core, here overriden
( const Cell::Core * face_p, const tag::SeenFrom &, const Cell neighbour,
  const tag::MayNotExist                                                        ) const

{	assert ( this->get_dim_plus_one() == face_p->get_dim() + 2 );
	Cell::Core * face_rev_p = face_p->reverse_p;
	assert ( face_rev_p );
	std::map<Mesh::Core*,Cell::Core*>::const_iterator
		it = face_rev_p->cell_behind_within.find ( this->core );
	assert ( it != face_rev_p->cell_behind_within .end() );  //  to finish !
	Cell::Core * cll_rev_p = it->second;
	assert ( cll_rev_p );  assert ( cll_rev_p->reverse_p );
	return cll_rev_p->reverse_p;                                        }  }

}	

Cell::Core * Mesh::STSI::cell_in_front_of
// virtual from Mesh::Core, here overriden
( const Cell::Core * face_p, const tag::SeenFrom &, const Cell neighbour,
  const tag::SurelyExists & se                                                  ) const
// 'se' defaults to tag::surely_exists, so method may be called with only one argument

{	assert ( this->get_dim_plus_one() == face_p->get_dim() + 2 );
	Cell::Core * face_rev_p = face_p->reverse_p;
	assert ( face_rev_p );
	std::map < Mesh::Core*, Cell::Core* > ::const_iterator
		it = face_rev_p->cell_behind_within .find ( this->core );
	assert ( it != face_rev_p->cell_behind_within .end() );  //  to finish !
	Cell::Core * cll_rev_p = it->second;
	assert ( cll_rev_p );  assert ( cll_rev_p->reverse_p );
	return cll_rev_p->reverse_p;                                        }  }

}	

Cell::Core * Mesh::STSI::cell_behind
// virtual from Mesh::Core, here overriden
( const Cell face, const tag::SeenFrom &, const Cell neighbour,
  const tag::MayNotExist                                        ) const
	
{	assert ( this->get_dim_plus_one() == face.dim() + 2 );
	std::map < Mesh::Core*, Cell::Core* > ::const_iterator
		it = face .core->cell_behind_within .find ( this->core );
	if ( it == face .core->cell_behind_within .end() ) return nullptr;
		// nothing behind us, we are touching the boundary
	Cell::Core * cll_p = it->second;
	if ( cll_p ) return cll_p;
	// now cll_p is null which means face_p is a singular point,
	// i.e. a touching point or a self-intersection
	// so its neighbours are stored separately
	std::vector < std::pair < Cell, Cell > > ::iterator itt = this->singular .begin();
	for ( itt .reset(); itt .in_range(); itt++ );
	{	std::pair < Cell, Cell > & p = *itt;
		Cell c1 & = p .first;
		Cell c2 & = p .second;
		Cell face_rev = face .reverse ( tag::surely_exists );
		if ( face_rev .belongs_to ( c1 .boundary(), tag::oriented ) ) return c2 .core;
		if ( face_rev .belongs_to ( c2 .boundary(), tag::oriented ) ) return c1 .core;
		assert ( false );                                                                 }  }


Cell::Core * Mesh::STSI::cell_behind
// virtual from Mesh::Core, here overriden
( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core neighbour,
  const tag::SurelyExists & se                                                  ) const
// 'se' defaults to tag::surely_exists, so method may be called with only one argument
	
{	assert ( this->get_dim_plus_one() == face_p->get_dim() + 2 );
	std::map<Mesh::Core*,Cell::Core*>::const_iterator
		it = face_p->cell_behind_within .find ( this->core );
	assert ( it != face_p->cell_behind_within .end() );
	Cell::Core * cll_p = it->second;
	if ( cll_p ) return cll_p;
	// now cll_p == nullptr which means face_p is a singular point,
	// i.e. a touching point or a self-intersection
	// so its neighbours are stored separately

}

*/

//-----------------------------------------------------------------------------//

namespace { // anonymous namespace, mimics static linkage

inline bool belongs_to_common
(	const Cell::Positive * cll, Mesh::Core * msh, const tag::CellHasLowDim & )

{	size_t dim = tag::Util::assert_diff ( msh->get_dim_plus_one(), cll->get_dim() + 1 );
	assert ( dim > 0 );
	const std::map < Mesh::Core *, Cell::field_to_meshes > & mmap =	cll->meshes [ dim ];
	std::map < Mesh::Core*, Cell::field_to_meshes > ::const_iterator it = mmap .find ( msh );
	return ( it != mmap .end() );                                                             }

}  // anonymous namespace


bool Cell::Positive::belongs_to  // virtual from Cell::Core
(	Mesh::Core * msh, const tag::CellHasLowDim &, const tag::NotOriented & ) const
// third argument defaults to tag::not_oriented, so method can be called with two arguments only

{	return belongs_to_common ( this, msh, tag::cell_has_low_dim );   }


bool Cell::Negative::belongs_to  // virtual from Cell::Core
(	Mesh::Core * msh, const tag::CellHasLowDim &, const tag::NotOriented & ) const
// third argument defaults to tag::not_oriented, so method can be called with two arguments only

{	assert ( this->reverse_attr.exists() );
	Cell::Positive * cll = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive* > ( this->reverse_attr.core );
	return belongs_to_common ( cll, msh, tag::cell_has_low_dim );   }


bool Cell::Positive::Vertex::belongs_to  // virtual from Cell::Core
(	Mesh::Core * msh, const tag::SameDim &, const tag::Oriented & ) const
// third argument defaults to tag::not_oriented, so method can be called with two arguments only

{	Cell::Positive::Segment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::Positive::Segment* > ( msh->cell_enclosed );
	assert ( seg->get_dim() == 1 );
	std::map < Cell::Positive::Segment*, short int > ::const_iterator it =
		this->segments .find ( seg );
	assert ( ( seg->tip_attr .core == this ) ==
					 ( ( it != this->segments .end() ) && ( it->second == 1 ) ) );
	return ( ( it != this->segments .end() ) && ( it->second == 1 ) );       }
	

bool Cell::Negative::Vertex::belongs_to  // virtual from Cell::Core
(	Mesh::Core * msh, const tag::SameDim &, const tag::Oriented & ) const
// third argument defaults to tag::not_oriented, so method can be called with two arguments only

{	Cell::Positive::Segment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::Positive::Segment* > ( msh->cell_enclosed );
	assert ( seg->get_dim() == 1 );
	Cell::Positive::Vertex * rev_this = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( this->reverse_attr.core );
	std::map<Cell::Positive::Segment*,short int>::const_iterator
		it = rev_this->segments .find ( seg );
	assert ( ( seg->base_attr .core == this ) ==
					 ( ( it != rev_this->segments .end() ) && ( it->second == -1 ) ) );
	return ( ( it != rev_this->segments .end() ) && ( it->second == -1 ) );     }
	
		
bool Cell::Positive::NotVertex::belongs_to  // virtual from Cell::Core
(	Mesh::Core * msh, const tag::SameDim &, const tag::Oriented & ) const
// third argument defaults to tag::not_oriented, so method can be called with two arguments only

{	assert ( this->get_dim() + 1 == msh->get_dim_plus_one() );
	std::map < Mesh::Core*, Cell::field_to_meshes_same_dim > ::const_iterator
		it = this->meshes_same_dim .find ( msh );
	return ( ( it != this->meshes_same_dim .end() ) && ( it->second .sign == 1 ) );  }
	

bool Cell::Negative::NotVertex::belongs_to  // virtual from Cell::Core
(	Mesh::Core * msh, const tag::SameDim &, const tag::Oriented & ) const
// third argument defaults to tag::not_oriented, so method can be called with two arguments only

{	assert ( this->get_dim() + 1 == msh->get_dim_plus_one() );
	Cell::Positive::NotVertex * rev_this = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( this->reverse_attr.core );
	std::map < Mesh::Core*, Cell::field_to_meshes_same_dim > ::const_iterator
		it = rev_this->meshes_same_dim .find ( msh );
	return ( ( it != rev_this->meshes_same_dim .end() ) && ( it->second .sign == -1 ) );  }
	

bool Cell::Core::belongs_to
( Mesh::Core * msh, const tag::SameDim &, const tag::NotOriented & ) const

{	if ( this->belongs_to ( msh, tag::same_dim, tag::oriented ) ) return true;
	if ( this->reverse_attr .exists() )
		return this->reverse_attr .core->belongs_to ( msh, tag::same_dim, tag::oriented );
	return false;                                                                        }
	
//-----------------------------------------------------------------------------//


Cell tag::Util::CellCore::tip ()  // virtual 
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Only segments have extremities." << std::endl;
	exit ( 1 );                                                                                   }

Cell tag::Util::CellCore::base ()  // virtual
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Only segments have extremities." << std::endl;
	exit ( 1 );                                                                                   }

Cell Cell::Positive::Segment::tip () { return this->tip_attr;  }
// virtual from Cell::Core, overridden here

Cell Cell::Positive::Segment::base () { return this->base_attr;  }
// virtual from Cell::Core, overridden here

Cell Cell::Negative::Segment::tip ()  // virtual from Cell::Core, overridden here
{	Cell pos_seg = this->reverse_attr;
	assert ( pos_seg .exists() );
	assert ( pos_seg .is_positive() );
	return pos_seg .base() .reverse ( tag::surely_exists );  }

Cell Cell::Negative::Segment::base ()  // virtual from Cell::Core, overridden here
{	Cell pos_seg = this->reverse_attr;
	assert ( pos_seg.exists() );
	assert ( pos_seg.is_positive() );
	return pos_seg.tip().reverse ( tag::surely_exists );   }

//-----------------------------------------------------------------------------//


#ifndef NDEBUG

std::string Cell::Positive::get_name ()  // virtual from Cell::Core
{	return this->name;  }

std::string Cell::Negative::get_name ()  // virtual from Cell::Core
{	return "r" + this->reverse_attr .core->get_name();  }


void Cell::Positive::Vertex::print_everything ( )  // virtual from Cell::Core

{	std::cout << "this is PositiveVertex " << this->name << std::endl;
	if ( this->meshes .size() > 0 )
	{	if ( this->meshes [0] .size() > 0 )
				std::cout << "meshes of index 0, dim 0 (segments disguised as meshes)" << std::endl;
		std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
			it = this->meshes [0] .begin(), it_e = this->meshes [0] .end();
		for ( ; it != it_e; it++ )
		{	Cell::Positive::Segment * seg = (Cell::Positive::Segment*) it->first;  // !!
			std::cout << seg->name << " " << it->second.counter_pos
		                         << " " << it->second.counter_neg << "  ";  }
			if ( this->meshes [0] .size() > 0 )  std::cout << std::endl;
		for ( size_t d = 1; d < this->meshes .size(); d++ )
		{	if ( this->meshes [d] .size() > 0 )
				std::cout << "meshes of index " << d << ", that is, dimension " << d << std::endl;
			std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
							itt = this->meshes [d] .begin(), itt_e = this->meshes [d] .end();
			for ( ; itt != itt_e; itt++ )
				std::cout << itt->first->get_name() << " " << itt->second.counter_pos
			                                      << " " << itt->second.counter_neg << "  ";
			if ( this->meshes [d] .size() > 0 )  std::cout << std::endl;                     }     }
	if ( this->reverse_attr .exists() ) std::cout << "has reverse" << std::endl;                  }


void Cell::Negative::Vertex::print_everything ( )  // virtual from Cell::Core

{	std::cout << "this is Cell::Negative::Vertex, reverse of "
            << this->reverse_attr .core->name << std::endl;   }


void Cell::Positive::Segment::print_everything ( )  // virtual from Cell::Core

{	std::cout << "this is PositiveSegment " << this->name << std::endl;
	std::cout << "base : " << this->base_attr .core->get_name() << std::endl;
	std::cout << "tip :  " << this->tip_attr .core->get_name()  << std::endl;
	if ( this->meshes .size() > 0 )
	{	for ( size_t d = 0; d < this->meshes .size(); d++ )
		{	if ( this->meshes [d] .size() > 0 )
				std::cout << "meshes of index " << d << ", that is, dimension " << d+1 << std::endl;
			std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
							it = this->meshes [d] .begin(), it_e = this->meshes [d] .end();
			for ( ; it != it_e; it++ )
				std::cout << it->first->get_name() << " " << it->second .counter_pos
			                                     << " " << it->second .counter_neg << "  ";
			if ( this->meshes [d] .size() > 0 )  std::cout << std::endl;                     }  }
	if ( this->reverse_attr .exists() ) std::cout << "has reverse" << std::endl;                }


void Cell::Negative::Segment::print_everything ( )  // virtual from Cell::Core

{	std::cout << "this is NegativeSegment " << this->get_name() << std::endl;     }


void Cell::Positive::HighDim::print_everything ( )  // virtual from Cell::Core
	
{	size_t dim = this->get_dim();
	std::cout << "this is Cell::Positive::HighDim of dim "
	          << dim << " " << this->name << std::endl;
	if ( this->meshes .size() > 0 )
	{	std::cout << "meshes above me :" << std::endl;
		for ( size_t d = 0; d < this->meshes .size(); d++ )
		{	if ( this->meshes [d] .size() > 0 )
				std::cout << "meshes of index " << d
			            << ", that is, dimension " << d+dim << std::endl;
			std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
							it = this->meshes [d] .begin(), it_e = this->meshes [d] .end();
			for ( ; it != it_e; it++ )
				std::cout << it->first->get_name() << " " << it->second.counter_pos
		                                       << " " << it->second.counter_neg << "  ";
			if ( this->meshes [d] .size() > 0 )  std::cout << std::endl;                   }  }
	if ( this->reverse_attr .exists() ) std::cout << "has reverse" << std::endl;            }


void Cell::Negative::HighDim::print_everything ( )  // virtual from Cell::Core
	
{	std::cout << "this is Cell::Negative::HighDim "
						<< this->reverse_attr.core->get_name() << std::endl;  }


std::string Mesh::Core::get_name ( )

{	if ( this->cell_enclosed ) return "bdry_of_" + this->cell_enclosed->name;
	else return this->name;                                                   }

void Mesh::ZeroDim::print_everything ( )  // virtual from Mesh::Core
{	assert ( false );  }

void Mesh::Connected::OneDim::print_everything ( )  // virtual from Mesh::Core
{	assert ( false );  }

void Mesh::Connected::HighDim::print_everything ( )  // virtual from Mesh::Core
{	assert ( false );  }

void Mesh::Fuzzy::print_everything ( )  // virtual from Mesh::Core

{	std::cout << "this is Mesh::OneDim::Positive " << this->get_name() << std::endl;
	for ( size_t d = 0; d < this->cells.size(); d++ )
	{	std::cout << "cells of dim " << d << " :" << std::endl;
		int counter = 0;
		std::list < Cell > ::iterator it = this->cells [d] .begin(),
		                              it_e = this->cells [d] .end();
		for ( ; it != it_e; it++, counter++ ) std::cout << (*it) .core->get_name() << " ";
		if ( counter > 0 )  std::cout << std::endl;                                         }  }

#endif  // NDEBUG


///////////////////////////////////////////////////////////////////
////////////    add/remove a cell to/from a mesh    ///////////////
///////////////////////////////////////////////////////////////////
////////////     this is where the magic happens     //////////////
//////////// (and also where we get our hands dirty) //////////////
//////////////////////////////////////////////////////////////////


namespace { // anonymous namespace, mimics static linkage


inline void add_cell_behind_above ( Cell::Positive::NotVertex * cll, Cell::Core * face_p )

// hidden in anonymous namespace
// called from Cell::Positive::***::glue_on_my_bdry

{	typedef std::map < Mesh::Core*, Cell::field_to_meshes_same_dim > maptype_f;
	maptype_f & tm0 = cll->meshes_same_dim;
	maptype_f::iterator it;
	for ( it = tm0.begin(); it != tm0.end(); ++it )
	{	Mesh::Core * msh = it->first;
		// std::list<Cell>::iterator wh = it->second.where;
		// const Cell other_cell = *wh;  assert ( other_cell.exists() );
		// we used to inquire whether  other_cell.core == cll
		// but 'wh' is well defined only for Fuzzy (and STSI) meshes
		// so we use the sign instead
		short int s = it->second.sign;
		if ( s == 1 )  // orientations match
//////////////////////////////////////////////////////////////////////////////////
		// inspired in item 24 of the book : Scott Meyers, Effective STL            //
		{	typedef std::map < Mesh::Core *, Cell > maptype;                          //
			maptype & cmd = face_p->cell_behind_within;                               //
			maptype::iterator lb = cmd.lower_bound(msh);                              //
			assert ( ( lb == cmd.end() ) or ( cmd.key_comp()(msh,lb->first) ) );      //
			cmd.emplace_hint ( lb, std::piecewise_construct,                          //
			      std::forward_as_tuple(msh),                                         //
			      std::forward_as_tuple(Cell(tag::whose_core_is,cll,                 //
																			 tag::previously_existing,                //
																			 tag::surely_not_null      )) );       }  //
/////////  code below is conceptually equivalent to the above  ///////////////////
//		face_p->cell_behind_within[msh] =                                 //
//			Cell ( tag::whose_core_is, cll, tag::previously_existing ) ;   //
//////////////////////////////////////////////////////////////////////////
		else  // mismatched orientations
		{	assert ( s == -1 );
			Cell::Core * rev_face_p = face_p->reverse_attr.core;
			Cell & other_cell = cll->reverse_attr;
			assert ( rev_face_p );
/////////////////////////////////////////////////////////////////////////////////////
			// inspired in item 24 of the book : Scott Meyers, Effective STL             //
			typedef std::map < Mesh::Core *, Cell > maptype;                             //
			maptype & cmd = rev_face_p->cell_behind_within;                              //
			maptype::iterator lb = cmd.lower_bound(msh);                                 //
			assert ( ( lb == cmd.end() ) or ( cmd.key_comp()(msh,lb->first) ) );         //
			cmd.emplace_hint ( lb, std::piecewise_construct,                             //
			      std::forward_as_tuple(msh), std::forward_as_tuple(other_cell) );  }  }  }  
/////////  code below is conceptually equivalent to the above  //////////////////////
//		rev_face_p->cell_behind_within[msh] = other_cell;          //
///////////////////////////////////////////////////////////////////
	

inline void remove_cell_behind_above ( Cell::Positive::NotVertex * cll, Cell::Core * face_p )
	
// hidden in anonymous namespace
// called from Cell::Positive::***::cut_from_my_bdry

{	if ( cll->meshes.size() == 0 ) return;
	typedef std::map < Mesh::Core*, Cell::field_to_meshes_same_dim > maptype_f;
	maptype_f & tm0 = cll->meshes_same_dim;
	maptype_f::iterator it;
	for ( it = tm0.begin(); it != tm0.end(); ++it )
	{	Mesh::Core * msh = it->first;
		// std::list<Cell>::iterator wh = it->second.where;
		// const Cell other_cell = *wh;  assert ( other_cell.exists() );
		// we used to inquire whether  other_cell.core == cll
		// but 'wh' is well defined only for Fuzzy (and STSI) meshes
		// so we use the sign instead
		short int s = it->second.sign;
		if ( s == 1 )  // orientations match
		{	assert ( face_p->cell_behind_within.find(msh) !=
		           face_p->cell_behind_within.end()        );
			assert ( face_p->cell_behind_within[msh].core == cll );
			face_p->cell_behind_within .erase ( msh );               }
		else  // mismatched orientations
		{	assert ( s == -1 );
			Cell::Core * rev_face_p { face_p->reverse_attr.core };
			assert ( rev_face_p );
			assert ( rev_face_p->cell_behind_within .find ( msh ) !=
		           rev_face_p->cell_behind_within .end()          );
			assert ( rev_face_p->cell_behind_within [ msh ] == cll->reverse_attr );
			rev_face_p->cell_behind_within .erase ( msh );                           }  }  }


}  // anonymous namespace

//-----------------------------------------------------------------------------//


void Cell::Positive::Vertex::glue_on_my_bdry ( Cell::Core * )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Vertices have no boundary." << std::endl;
	exit ( 1 );                                               }


void Cell::Positive::Vertex::glue_on_my_bdry ( Cell::Core *, const tag::DoNotBother & )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Vertices have no boundary." << std::endl;
	exit ( 1 );                                               }


void Cell::Positive::Segment::glue_on_my_bdry ( Cell::Core * ver )
// virtual from Cell::Core

{	assert ( ver->get_dim() == 0 );
	ver->add_to_seg ( this );
	// 'add_to_seg' is virtual, so the computer will choose the right version
	// (Cell::Positive::Vertex::add_to_seg or Cell::Negative::Vertex::add_to_seg)
	add_cell_behind_above ( this, ver );        }
	

void Cell::Positive::Segment::glue_on_my_bdry ( Cell::Core * ver, const tag::DoNotBother & )
// virtual from Cell::Core

{	assert ( ver->get_dim() == 0 );
	ver->add_to_seg ( this, tag::do_not_bother );
	// 'add_to_seg' is virtual, so the computer will choose the right version
	// (Cell::Positive::Vertex::add_to_seg or Cell::Negative::Vertex::add_to_seg)
	add_cell_behind_above ( this, ver );                    }
	

void Cell::Positive::HighDim::glue_on_my_bdry ( Cell::Core * face )
// virtual from Cell::Core

{	assert ( this->get_dim() == face->get_dim() + 1 );
	assert ( this->boundary() .is_positive() );
	face->add_to_bdry ( this->boundary() .core );
	// 'add_to_bdry' is virtual, so the computer will choose the right version
	add_cell_behind_above ( this, face );                  }


void Cell::Positive::HighDim::glue_on_my_bdry ( Cell::Core * face, const tag::DoNotBother & )
// virtual from Cell::Core

{	assert ( this->get_dim() == face->get_dim() + 1 );
	assert ( this->boundary() .is_positive() );
	face->add_to_bdry ( this->boundary() .core, tag::do_not_bother );
	// 'add_to_bdry' is virtual, so the computer will choose the right version
	add_cell_behind_above ( this, face );                                      }


void Cell::Negative::glue_on_my_bdry ( Cell::Core * cll )
// virtual from Cell::Core

{	assert ( cll->reverse_attr .exists() );
	assert ( this->reverse_attr .exists() );
	this->reverse_attr .core->glue_on_my_bdry ( cll->reverse_attr .core );  }


void Cell::Negative::glue_on_my_bdry ( Cell::Core * cll, const tag::DoNotBother & )
// virtual from Cell::Core

{	assert ( cll->reverse_attr .exists() );
	assert ( this->reverse_attr .exists() );
	this->reverse_attr .core->glue_on_my_bdry ( cll->reverse_attr .core, tag::do_not_bother );  }


void Cell::Positive::Vertex::cut_from_my_bdry ( Cell::Core * )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Vertices have no boundary." << std::endl;
	exit ( 1 );                                               }


void Cell::Positive::Vertex::cut_from_my_bdry ( Cell::Core *, const tag::DoNotBother & )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Vertices have no boundary." << std::endl;
	exit ( 1 );                                               }


void Cell::Positive::Segment::cut_from_my_bdry ( Cell::Core * ver )
// virtual from Cell::Core

{	assert ( ver->get_dim() == 0 );
	remove_cell_behind_above ( this, ver );
	// 'remove_from_seg' is virtual, so the computer will choose the right version
	// (Cell::Positive::Vertex::remove_from_seg or Cell::Negative::Vertex::remove_from_seg)
	ver->remove_from_seg ( this );  }
	

void Cell::Positive::Segment::cut_from_my_bdry ( Cell::Core * ver, const tag::DoNotBother & )
// virtual from Cell::Core

{	assert ( ver->get_dim() == 0 );
	remove_cell_behind_above ( this, ver );
	// 'remove_from_seg' is virtual, so the computer will choose the right version
	// (Cell::Positive::Vertex::remove_from_seg or Cell::Negative::Vertex::remove_from_seg)
	ver->remove_from_seg ( this, tag::do_not_bother );  }
	

void Cell::Positive::HighDim::cut_from_my_bdry ( Cell::Core * face )
// virtual from Cell::Core

{	assert ( this->get_dim() == face->get_dim() + 1 );
	assert ( this->boundary_attr .is_positive() );
	remove_cell_behind_above ( this, face );
	// 'remove_from_bdry' is virtual, so the computer will choose the right version
	face->remove_from_bdry ( this->boundary_attr .core );   }


void Cell::Positive::HighDim::cut_from_my_bdry ( Cell::Core * face, const tag::DoNotBother & )
// virtual from Cell::Core

{	assert ( this->get_dim() == face->get_dim() + 1 );
	assert ( this->boundary_attr .is_positive() );
	remove_cell_behind_above ( this, face );
	// 'remove_from_bdry' is virtual, so the computer will choose the right version
	face->remove_from_bdry ( this->boundary_attr .core, tag::do_not_bother );  }


void Cell::Negative::cut_from_my_bdry ( Cell::Core * cll )
// virtual from Cell::Core

{	assert ( cll->reverse_attr .exists() );
	assert ( this->reverse_attr .exists() );
	this->reverse_attr .core->cut_from_my_bdry ( cll->reverse_attr .core );  }


void Cell::Negative::cut_from_my_bdry ( Cell::Core * cll, const tag::DoNotBother & )
// virtual from Cell::Core

{	assert ( cll->reverse_attr .exists() );
	assert ( this->reverse_attr .exists() );
	this->reverse_attr .core->cut_from_my_bdry ( cll->reverse_attr .core, tag::do_not_bother );  }

//-----------------------------------------------------------------------------//


void Cell::Negative::compute_sign
( short int & cp, short int & cn, Mesh::Core * const cell_bdry )
// we only use this method for positive cells
{	assert ( false );  }
	

void Cell::Positive::Vertex::compute_sign
( short int & cp, short int & cn, Mesh::Core * const cell_bdry )
	
// 'this' is a face of 'cll', so it is a cell in 'cell_bdry' of maximum dimension
// we just want to know if it is positively oriented or not

{	assert ( this->get_dim() + 1 == cell_bdry->get_dim_plus_one() );
	std::map < Cell::Positive::Segment*, short int > &
		cemd = this->segments;
	Cell::Positive::Segment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::Positive::Segment* > ( cell_bdry->cell_enclosed );
	std::map < Cell::Positive::Segment*, short int > ::iterator
		map_iter = cemd .find(seg);
	assert ( map_iter != cemd.end() );
	short int sign = map_iter->second;
	if ( sign == 1 )  {  cp = 1;  cn = 0;  }
	else  {  assert ( sign == -1 );  cp = 0;  cn = 1;  }                           }


void Cell::Positive::NotVertex::compute_sign
( short int & cp, short int & cn, Mesh::Core * const cell_bdry )
	
// 'this' is a face of 'cll', so it is a cell in 'cell_bdry' of maximum dimension
// we just want to know if it is positively oriented or not
	
{	assert ( this->get_dim() + 1 == cell_bdry->get_dim_plus_one() );
	std::map < Mesh::Core*, Cell::field_to_meshes_same_dim > &
		cemd = this->meshes_same_dim;
	std::map < Mesh::Core*, Cell::field_to_meshes_same_dim > ::iterator
		map_iter = cemd .find(cell_bdry);
	assert ( map_iter != cemd .end() );
	short int sign = map_iter->second.sign;
	if ( sign == 1 )  {  cp = 1;  cn = 0;  }
	else  {  assert ( sign == -1 );  cp = 0;  cn = 1;  }                }

//-----------------------------------------------------------------------------//
	

namespace {  // anonymous namespace, mimics static linkage

// here is where the low-level linking between cells and meshes happens

	
inline void add_link_zero_dim  // hidden in anonymous namespace
(	Cell::Positive::Vertex * const ver, Cell::Positive::Segment * const seg )

// This function "makes the link" between a cell and a mesh.
// This "link" asserts that the cell belongs to the mesh.

// this function updates ver->segments

{	assert ( ver );  assert ( seg );
/////////////////////////////////////////////////////////////////////////////
	// inspired in item 24 of the book : Scott Meyers, Effective STL         //
	typedef std::map < Cell::Positive::Segment *, short int > maptype;       //
	maptype & cmd = ver->segments;                                           //
	maptype::iterator lb = cmd.lower_bound(seg);                             //
	assert ( ( lb == cmd.end() ) or ( cmd.key_comp()(seg,lb->first) ) );     //
	cmd.emplace_hint ( lb, std::piecewise_construct,                         //
	      std::forward_as_tuple(seg), std::forward_as_tuple(1) );            //
/////////  code below is conceptually equivalent to the above  //////////////
//	assert ( ver->segments.find(seg) == ver->segments.end() );     //
//	ver->segments[seg] = 1;                                        //
/////////////////////////////////////////////////////////////////////

} // end of add_link_zero_dim


inline void add_link_zero_dim_rev  // hidden in anonymous namespace
(	Cell::Positive::Vertex * const ver, Cell::Positive::Segment * const seg )

// This function "makes the link" between a cell and a mesh.
// This "link" asserts that the cell belongs to the mesh.

// this function updates ver->segments

{	assert ( ver );  assert ( seg );
/////////////////////////////////////////////////////////////////////////////
	// inspired in item 24 of the book : Scott Meyers, Effective STL         //
	typedef std::map < Cell::Positive::Segment *, short int > maptype;       //
	maptype & cmd = ver->segments;                                           //
	maptype::iterator lb = cmd.lower_bound(seg);                             //
	assert ( ( lb == cmd.end() ) or ( cmd.key_comp()(seg,lb->first) ) );     //
	cmd.emplace_hint ( lb, std::piecewise_construct,                         //
	       std::forward_as_tuple(seg), std::forward_as_tuple(-1) );          //
/////////  code below is conceptually equivalent to the above  //////////////
//	assert ( ver->segments.find(seg) == ver->segments.end() );     //
//	ver->segments[seg] = -1;                                       //
/////////////////////////////////////////////////////////////////////

} // end of add_link_zero_dim_rev


inline void add_link_same_dim  // hidden in anonymous namespace
(	Cell::Positive::NotVertex * const cll, Mesh::Core * const msh )

// This function "makes the link" between a cell and a mesh.
// This "link" asserts that the cell belongs to the mesh.
// Here the dimensions are equal and the cell is positive.

// this function updates cll->meshes_same_dim
// also, calls the virtual method msh->add_to_my_cells
// which, for fuzzy and stsi meshes, updates msh->cells
// returning an iterator within msh->cells
// for other kinds of meshes does nothing, returns garbage

{	assert ( cll );  assert ( msh );
	size_t cll_dim = cll->get_dim();
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );
/////////////////////////////////////////////////////////////////////////////
	// inspired in item 24 of the book : Scott Meyers, Effective STL         //
	typedef std::map<Mesh::Core*,Cell::field_to_meshes_same_dim> maptype;    //
	maptype & cmd = cll->meshes_same_dim;                                    //
	maptype::iterator lb = cmd.lower_bound(msh);                             //
	assert ( ( lb == cmd.end() ) or ( cmd.key_comp()(msh,lb->first) ) );     //
	cmd.emplace_hint ( lb, std::piecewise_construct,                         //
	      std::forward_as_tuple(msh),                                        // 
	      std::forward_as_tuple(1,msh->add_to_my_cells(cll,cll_dim)) );      //
/////////  code below is conceptually equivalent to the above  //////////////////
//	assert ( cll->meshes_same_dim.find(msh) == cll->meshes_same_dim.end() );   //
//	msh->cells[cll->dim].push_front(cll);                                      //
//	Cell::field_to_meshes_same_dim field;                                      //
//	field.sign = 1;                                                            //
//	field.where = msh->cells[cll->dim]->begin();                               //
//	cll->meshes_same_dim[msh] = field;                                         //
/////////////////////////////////////////////////////////////////////////////////

} // end of add_link_same_dim


inline void add_link_same_dim  // hidden in anonymous namespace
(	Cell::Positive::NotVertex * const cll, Mesh::Core * const msh, const tag::DoNotBother & )

// This function "makes the link" between a cell and a mesh.
// This "link" asserts that the cell belongs to the mesh.
// Here the dimensions are equal and the cell is positive.

// this function updates cll->meshes_same_dim
// also, calls the virtual method msh->add_to_my_cells
// which, for fuzzy and stsi meshes, updates msh->cells
// returning an iterator within msh->cells
// for other kinds of meshes does nothing, returns garbage

{	assert ( cll );  assert ( msh );
	size_t cll_dim = cll->get_dim();
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );
////////////////////////////////////////////////////////////////////////////////////////
	// inspired in item 24 of the book : Scott Meyers, Effective STL                    //
	typedef std::map<Mesh::Core*,Cell::field_to_meshes_same_dim> maptype;               //
	maptype & cmd = cll->meshes_same_dim;                                               //
	maptype::iterator lb = cmd.lower_bound(msh);                                        //
	assert ( ( lb == cmd.end() ) or ( cmd.key_comp()(msh,lb->first) ) );                //
	cmd.emplace_hint ( lb, std::piecewise_construct,                                    //
	   std::forward_as_tuple(msh),                                                      // 
	   std::forward_as_tuple(1,msh->add_to_my_cells(cll,cll_dim,tag::do_not_bother)) ); //
/////////  code below is conceptually equivalent to the above  //////////////////>//////
//	assert ( cll->meshes_same_dim.find(msh) == cll->meshes_same_dim.end() );   //
//	msh->cells[cll->dim].push_front(cll);                                      //
//	Cell::field_to_meshes_same_dim field;                                      //
//	field.sign = 1;                                                            //
//	field.where = msh->cells[cll->dim]->begin();                               //
//	cll->meshes_same_dim[msh] = field;                                         //
/////////////////////////////////////////////////////////////////////////////////

} // end of add_link_same_dim with tag::do_not_bother


inline void add_link_same_dim_rev  // hidden in anonymous namespace
(	Cell::Core * const o_cll, Cell::Positive::NotVertex * const cll, Mesh::Core * const msh )

// This function "makes the link" between a cell and a mesh.
// This "link" asserts that the cell belongs to the mesh.
// Here the dimensions are equal and the cell is negative.

// this function updates cll->meshes_same_dim
// also, calls the virtual method msh->add_to_my_cells
// which, for fuzzy and stsi meshes, updates msh->cells
// returning an iterator within msh->cells
// for other kinds of meshes does nothing, returns garbage

{	assert ( cll );  assert ( o_cll );  assert ( msh );
	assert ( ! o_cll->is_positive() );
	size_t cll_dim = cll->get_dim();
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );
/////////////////////////////////////////////////////////////////////////////
	// inspired in item 24 of the book : Scott Meyers, Effective STL         //
	typedef std::map<Mesh::Core*,Cell::field_to_meshes_same_dim> maptype;    //
	maptype & cmd = cll->meshes_same_dim;                                    //
	maptype::iterator lb = cmd.lower_bound(msh);                             //
	assert ( ( lb == cmd.end() ) or ( cmd.key_comp()(msh,lb->first) ) );     //
	cmd.emplace_hint ( lb, std::piecewise_construct,                         //
	      std::forward_as_tuple(msh),                                        // 
	      std::forward_as_tuple(-1,msh->add_to_my_cells(o_cll,cll_dim)) );   //
/////////  code below is conceptually equivalent to the above  //////////////////
//	assert ( cll->meshes_same_dim.find(msh) == cll->meshes_same_dim.end() );   //
//	msh->cells[cll->dim].push_front(o_cll);                                    //
//	Cell::field_to_meshes_same_dim field;                                      //
//	field.sign = -1;                                                           //
//	field.where = msh->cells[cll->dim]->begin();                               //
//	cll->meshes_same_dim[msh] = field;                                         //
/////////////////////////////////////////////////////////////////////////////////

} // end of add_link_same_dim_rev


inline void add_link_same_dim_rev  // hidden in anonymous namespace
(	Cell::Core * const o_cll, Cell::Positive::NotVertex * const cll,
  Mesh::Core * const msh, const tag::DoNotBother &                )

// This function "makes the link" between a cell and a mesh.
// This "link" asserts that the cell belongs to the mesh.
// Here the dimensions are equal and the cell is negative.

// this function updates cll->meshes_same_dim
// also, calls the virtual method msh->add_to_my_cells
// which, for fuzzy and stsi meshes, updates msh->cells
// returning an iterator within msh->cells
// for other kinds of meshes does nothing, returns garbage

{	assert ( cll );  assert ( o_cll );  assert ( msh );
	assert ( ! o_cll->is_positive() );
	size_t cll_dim = cll->get_dim();
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );
////////////////////////////////////////////////////////////////////////////////////////////
	// inspired in item 24 of the book : Scott Meyers, Effective STL                        //
	typedef std::map<Mesh::Core*,Cell::field_to_meshes_same_dim> maptype;                   //
	maptype & cmd = cll->meshes_same_dim;                                                   //
	maptype::iterator lb = cmd.lower_bound(msh);                                            //
	assert ( ( lb == cmd.end() ) or ( cmd.key_comp()(msh,lb->first) ) );                    //
	cmd.emplace_hint ( lb, std::piecewise_construct,                                        //
	    std::forward_as_tuple(msh),                                                         // 
	    std::forward_as_tuple(-1,msh->add_to_my_cells(o_cll,cll_dim,tag::do_not_bother)) ); //
/////////  code below is conceptually equivalent to the above  /////////////////////////////
//	assert ( cll->meshes_same_dim.find(msh) == cll->meshes_same_dim.end() );   //
//	msh->cells[cll->dim].push_front(o_cll);                                    //
//	Cell::field_to_meshes_same_dim field;                                      //
//	field.sign = -1;                                                           //
//	field.where = msh->cells[cll->dim]->begin();                               //
//	cll->meshes_same_dim[msh] = field;                                         //
/////////////////////////////////////////////////////////////////////////////////

} // end of add_link_same_dim_rev with tag::do_not_bother


inline void add_link  // hidden in anonymous namespace
(	Cell::Positive * const cll, Mesh::Core * const msh, const short int cp, const short int cn )

// This function "makes the link" between a cell and a mesh.
// This "link" asserts that the cell belongs to the mesh.
// Here 'cll' has lower dimension, so it may appear several times in the mesh 'msh',
// with both orientations. The counter 'cp' asserts how many times 'cll' appears
// in the mesh with positive orientation,
// while 'cn' says the same thing for negative orientation.

// this function updates cll->meshes
// also, calls the virtual method msh->add_to_my_cells
// which, for fuzzy and stsi meshes, updates msh->cells
// returning an iterator within msh->cells
// for other kinds of meshes does nothing, returns garbage

{	assert ( cll );  assert ( msh );
	size_t cll_dim = cll->get_dim(),
	       msh_dim_p1 = msh->get_dim_plus_one();
	assert ( msh_dim_p1 > cll_dim + 1 );
	size_t dif_dim = msh_dim_p1 - cll_dim - 1;
	assert ( cll->meshes.size() > dif_dim );
///////////////////////////////////////////////////////////////////////////////
	// inspired in item 24 of the book : Scott Meyers, Effective STL           //
	typedef std::map<Mesh::Core*,Cell::field_to_meshes> maptype;               //
	maptype & cmd = cll->meshes[dif_dim];                                      //
	maptype::iterator lb = cmd.lower_bound(msh);                               //
	if ( ( lb == cmd.end() ) or ( cmd.key_comp()(msh,lb->first) ) )            //
		cmd.emplace_hint ( lb, std::piecewise_construct,                         //
		      std::forward_as_tuple(msh),                                        // 
		      std::forward_as_tuple(cp,cn,msh->add_to_my_cells(cll,cll_dim)) );  //
	else                                                                       //
	{ lb->second.counter_pos += cp;                                            //
	  lb->second.counter_neg += cn;   }                                        //
////////// code below is conceptually equivalent to the above /////////////////
//	if ( cll->meshes[dif_dim].find(msh)==cll->meshes[dif_dim].end() )   //
//	{ msh->cells[cll->dim].push_front(cll);                             //
//	  Cell::field_to_meshes field;                                      //
//	  field.counter_pos = cp;                                           //
//	  field.counter_neg = cn;                                           //
//	  field.where = msh->cells[cll->dim]->begin();                      //
//	  cll->meshes[dif_dim][msh] = field;             }                  //
//	else                                                                //
//	{ cll->meshes[dif_dim][msh].counter_pos += cp;                      //
//	  cll->meshes[dif_dim][msh].counter_neg += cn;   }                  //
//////////////////////////////////////////////////////////////////////////

} // end of add_link


inline void add_link  // hidden in anonymous namespace
(	Cell::Positive * const cll, Mesh::Core * const msh,
  const short int cp, const short int cn, const tag::DoNotBother & )

// This function "makes the link" between a cell and a mesh.
// This "link" asserts that the cell belongs to the mesh.
// Here 'cll' has lower dimension, so it may appear several times in the mesh 'msh',
// with both orientations. The counter 'cp' asserts how many times 'cll' appears
// in the mesh with positive orientation,
// while 'cn' says the same thing for negative orientation.

// this function updates cll->meshes
// also, calls the virtual method msh->add_to_my_cells
// which, for fuzzy and stsi meshes, updates msh->cells
// returning an iterator within msh->cells
// for other kinds of meshes does nothing, returns garbage

{	assert ( cll );  assert ( msh );
	size_t cll_dim = cll->get_dim(),
	       msh_dim_p1 = msh->get_dim_plus_one();
	assert ( msh_dim_p1 > cll_dim + 1 );
	size_t dif_dim = msh_dim_p1 - cll_dim - 1;
	assert ( cll->meshes.size() > dif_dim );
//////////////////////////////////////////////////////////////////////////////////////////////////
	// inspired in item 24 of the book : Scott Meyers, Effective STL                              //
	typedef std::map<Mesh::Core*,Cell::field_to_meshes> maptype;                                  //
	maptype & cmd = cll->meshes[dif_dim];                                                         //
	maptype::iterator lb = cmd.lower_bound(msh);                                                  //
	if ( ( lb == cmd.end() ) or ( cmd.key_comp()(msh,lb->first) ) )                               //
		cmd.emplace_hint ( lb, std::piecewise_construct,                                            //
		      std::forward_as_tuple(msh),                                                           // 
					std::forward_as_tuple(cp,cn,msh->add_to_my_cells(cll,cll_dim,tag::do_not_bother)) );  //
	else                                                                                          //
	{ lb->second.counter_pos += cp;                                                               //
	  lb->second.counter_neg += cn;   }                                                           //
////////// code below is conceptually equivalent to the above ////////////////////////////////////
//	if ( cll->meshes[dif_dim].find(msh)==cll->meshes[dif_dim].end() )   //
//	{ msh->cells[cll->dim].push_front(cll);                             //
//	  Cell::field_to_meshes field;                                      //
//	  field.counter_pos = cp;                                           //
//	  field.counter_neg = cn;                                           //
//	  field.where = msh->cells[cll->dim]->begin();                      //
//	  cll->meshes[dif_dim][msh] = field;             }                  //
//	else                                                                //
//	{ cll->meshes[dif_dim][msh].counter_pos += cp;                      //
//	  cll->meshes[dif_dim][msh].counter_neg += cn;   }                  //
//////////////////////////////////////////////////////////////////////////

} // end of add_link with tag::do_not_bother


inline void link_face_to_msh  // hidden in anonymous namespace
( Cell::Core * const face, Cell::Positive::NotVertex * cll,
  Mesh::Core * const msh, const short int cp, const short int cn )
// just a block of code for make_deep_connections

{	assert ( face );  assert ( cll );
	assert ( face->get_dim() < cll->get_dim() );
	assert ( cll->get_dim() + 1 == msh->get_dim_plus_one() );
	Cell::Positive * face_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive* > ( face );
	typedef std::map < Mesh::Core*, Cell::field_to_meshes_same_dim > maptype;
	maptype & cemd = cll->meshes_same_dim;
	maptype::iterator map_iter = cemd.find ( msh );
	assert ( map_iter != cemd.end() );
	short int mis = map_iter->second.sign;
	if ( mis == 1 ) add_link ( face_p, msh, cp, cn );
	else  // we switch the two counters
	{	assert ( mis == -1 ); add_link ( face_p, msh, cn, cp );  }              }
	

inline void link_face_to_msh  // hidden in anonymous namespace
( Cell::Core * const face, Cell::Positive::NotVertex * cll,
  Mesh::Core * const msh, const short int cp, const short int cn,
  const tag::DoNotBother &                                       )
// just a block of code for make_deep_connections

{	assert ( face );  assert ( cll );
	assert ( face->get_dim() < cll->get_dim() );
	assert ( cll->get_dim() + 1 == msh->get_dim_plus_one() );
	Cell::Positive * face_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive* > ( face );
	typedef std::map < Mesh::Core*, Cell::field_to_meshes_same_dim > maptype;
	maptype & cemd = cll->meshes_same_dim;
	maptype::iterator map_iter = cemd.find ( msh );
	assert ( map_iter != cemd.end() );
	short int mis = map_iter->second.sign;
	if ( mis == 1 ) add_link ( face_p, msh, cp, cn, tag::do_not_bother );
	else  // we switch the two counters
	{	assert ( mis == -1 ); add_link ( face_p, msh, cn, cp, tag::do_not_bother );  }  }
	

inline void link_face_to_higher  // hidden in anonymous namespace
(	Cell::Core * const face, Cell::Positive::NotVertex * const pmce,
  const short int cp, const short int cn                           )
// just a block of code for make_deep_connections

{	assert ( face );  assert ( pmce );
	assert ( face->get_dim() + 1 <= pmce->get_dim() );
	Cell::Positive * face_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive* > ( face );
	{ // just a block of code for hiding maptype, cemd and map_iter
	// pmce->meshes[0] is empty, we use pmce->meshes_same_dim instead
	typedef std::map < Mesh::Core*, Cell::field_to_meshes_same_dim > maptype;
	maptype & cemd = pmce->meshes_same_dim;
	maptype::iterator map_iter;
	// we now loop over all meshes of given dimension
	for ( map_iter = cemd.begin(); map_iter != cemd.end(); ++map_iter )
	{	assert ( map_iter->first->get_dim_plus_one() == pmce->get_dim() + 1 );
		Cell::field_to_meshes_same_dim & mis = map_iter->second;
		if ( mis.sign == 1 ) add_link ( face_p, map_iter->first, cp, cn );
		else  // we just switch cp and cn
		{	assert ( mis.sign == -1 );
			add_link ( face_p, map_iter->first, cn, cp );  }                      }
	} // just a block of code for hiding maptype, cemd and map_iter
	for ( size_t dif_dim = 1; dif_dim < pmce->meshes.size(); dif_dim++ )
	{	typedef std::map < Mesh::Core*, Cell::field_to_meshes > maptype;
		maptype & cemd = pmce->meshes[dif_dim];
		maptype::iterator map_iter;
		// we now loop over all meshes of given dimension
		for ( map_iter = cemd.begin(); map_iter != cemd.end(); ++map_iter )
		{	assert ( map_iter->first->get_dim_plus_one() == dif_dim + pmce->get_dim() + 1 );
			Cell::field_to_meshes & mis = map_iter->second;
			add_link ( face_p, map_iter->first,
			           cp*mis.counter_pos + cn*mis.counter_neg,
			           cn*mis.counter_pos + cp*mis.counter_neg );                              }  }  }


inline void link_face_to_higher  // hidden in anonymous namespace
(	Cell::Core * const face, Cell::Positive::NotVertex * const pmce,
  const short int cp, const short int cn, const tag::DoNotBother & )
// just a block of code for make_deep_connections

{	assert ( face );  assert ( pmce );
	assert ( face->get_dim() + 1 <= pmce->get_dim() );
	Cell::Positive * face_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive* > ( face );
	{ // just a block of code for hiding maptype, cemd and map_iter
	// pmce->meshes[0] is empty, we use pmce->meshes_same_dim instead
	typedef std::map < Mesh::Core*, Cell::field_to_meshes_same_dim > maptype;
	maptype & cemd = pmce->meshes_same_dim;
	maptype::iterator map_iter;
	// we now loop over all meshes of given dimension
	for ( map_iter = cemd.begin(); map_iter != cemd.end(); ++map_iter )
	{	assert ( map_iter->first->get_dim_plus_one() == pmce->get_dim() + 1 );
		Cell::field_to_meshes_same_dim & mis = map_iter->second;
		if ( mis.sign == 1 ) add_link ( face_p, map_iter->first, cp, cn, tag::do_not_bother );
		else  // we just switch cp and cn
		{	assert ( mis.sign == -1 );
			add_link ( face_p, map_iter->first, cn, cp, tag::do_not_bother );  }                }
	} // just a block of code for hiding maptype, cemd and map_iter
	for ( size_t dif_dim = 1; dif_dim < pmce->meshes.size(); dif_dim++ )
	{	typedef std::map < Mesh::Core*, Cell::field_to_meshes > maptype;
		maptype & cemd = pmce->meshes[dif_dim];
		maptype::iterator map_iter;
		// we now loop over all meshes of given dimension
		for ( map_iter = cemd.begin(); map_iter != cemd.end(); ++map_iter )
		{	assert ( map_iter->first->get_dim_plus_one() == dif_dim + pmce->get_dim() + 1 );
			Cell::field_to_meshes & mis = map_iter->second;
			add_link ( face_p, map_iter->first,
			           cp*mis.counter_pos + cn*mis.counter_neg,
			           cn*mis.counter_pos + cp*mis.counter_neg, tag::do_not_bother );         }  }  }


inline void compute_cp_cn  // hidden in anonymous namespace
(	short int & cp, short int & cn, Cell::Core * const face, Mesh::Core * const cell_bdry )
// just a block of code for make_deep_connections

// 'face' is a vertex or side of 'cll', so it is a cell in 'cell_bdry', of lower dimension
// we want to know how many times it appears as positive or as negative within 'cell_bdry'

{	assert ( face );
	size_t dif_dim = tag::Util::assert_diff
		( cell_bdry->get_dim_plus_one(), face->get_dim() + 1 );
	Cell::Positive * const face_p = tag::Util::assert_cast
		< Cell::Core*const, Cell::Positive*const > ( face );
	typedef std::map<Mesh::Core*,Cell::field_to_meshes> maptype;
	maptype & cemd = face_p->meshes[dif_dim];
	maptype::iterator map_iter = cemd.find(cell_bdry);
	assert ( map_iter != cemd.end() );
	Cell::field_to_meshes & mis = map_iter->second;
	cp = mis.counter_pos;  cn = mis.counter_neg;                                       }

	
inline void make_deep_connections_0d  // hidden in anonymous namespace
(	Cell::Positive::Vertex * const ver, Cell::Positive::Segment * const seg )

// make far connections when adding a positive vertex

{	assert ( ver );  assert ( seg );
	add_link_zero_dim ( ver, seg );
	link_face_to_higher ( ver, seg, 1, 0 );  }
	

inline void make_deep_connections_0d  // hidden in anonymous namespace
(	Cell::Positive::Vertex * const ver, Cell::Positive::Segment * const seg,
  const tag::DoNotBother &                                                )

// make far connections when adding a positive vertex

{	assert ( ver );  assert ( seg );
	add_link_zero_dim ( ver, seg );
	link_face_to_higher ( ver, seg, 1, 0, tag::do_not_bother );  }
	

inline void make_deep_connections_0d_rev  // hidden in anonymous namespace
(	Cell::Positive::Vertex * const ver, Cell::Positive::Segment * const seg )

// make far connections when adding a negative vertex

{	assert ( ver );  assert ( seg );
	add_link_zero_dim_rev ( ver, seg );
	link_face_to_higher ( ver, seg, 0, 1 );  }
	

inline void make_deep_connections_0d_rev  // hidden in anonymous namespace
(	Cell::Positive::Vertex * const ver, Cell::Positive::Segment * const seg,
  const tag::DoNotBother &                                                )

// make far connections when adding a negative vertex

{	assert ( ver );  assert ( seg );
	add_link_zero_dim_rev ( ver, seg );
	link_face_to_higher ( ver, seg, 0, 1, tag::do_not_bother );  }
	

inline void make_deep_connections_1d  // hidden in anonymous namespace
(	Cell::Positive::Segment * const seg, Mesh::Core * const msh, const tag::MeshIsNotBdry & )

// make far connections when adding a positive segment
// see paragraph 11.9 in the manual
	
{	assert ( msh->get_dim_plus_one() == 2 );
	assert ( seg );  assert ( seg->get_dim() == 1 );
	// we could assert  msh->cell_enclosed == nullptr  (i.e. mesh is not a boundary)
	// however, we sometimes call this method on the boundary of a
	// freshly created cell (pretending it is not a boundary)

	add_link_same_dim ( seg, msh );

	Cell A = seg->base_attr.reverse ( tag::surely_exists );
	Cell::Positive::Vertex * Ap = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	link_face_to_msh ( Ap, seg, msh, 0, 1 );
	Cell B = seg->tip_attr;
	Cell::Positive::Vertex * Bp = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	link_face_to_msh ( Bp, seg, msh, 1, 0 );

} // end of make_deep_connections_1d with tag::mesh_is_not_bdry


inline void make_deep_connections_1d  // hidden in anonymous namespace
(	Cell::Positive::Segment * const seg, Mesh::Core * const msh,
  const tag::MeshIsNotBdry &, const tag::DoNotBother &        )

// make far connections when adding a positive segment
// see paragraph 11.9 in the manual

// tag::do_not_bother is only relevant for Mesh::Connected::***Dim
	
{	assert ( msh->get_dim_plus_one() == 2 );
	assert ( seg );  assert ( seg->get_dim() == 1 );
	// we could assert  msh->cell_enclosed == nullptr  (i.e. mesh is not a boundary)
	// however, we sometimes call this method on the boundary of a
	// freshly created cell (pretending it is not a boundary)

	add_link_same_dim ( seg, msh, tag::do_not_bother );

	Cell A = seg->base_attr.reverse ( tag::surely_exists );
	Cell::Positive::Vertex * Ap = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	link_face_to_msh ( Ap, seg, msh, 0, 1, tag::do_not_bother );
	Cell B = seg->tip_attr;
	Cell::Positive::Vertex * Bp = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	link_face_to_msh ( Bp, seg, msh, 1, 0, tag::do_not_bother );

} // end of make_deep_connections_1d with tag::mesh_is_not_bdry, tag::do_not_bother


inline void make_deep_connections_1d  // hidden in anonymous namespace
(	Cell::Positive::Segment * const seg, Mesh::Core * const msh, const tag::MeshIsBdry & )

// make far connections when adding a positive cell
// see paragraph 11.9 in the manual
	
{	assert ( msh->get_dim_plus_one() == 2 );
	assert ( seg );
	assert ( seg->get_dim() == 1 );

	add_link_same_dim ( seg, msh );

	// for all meshes strictly above msh
	Cell::Core * mce = msh->cell_enclosed;
	assert ( mce );
	Cell::Positive::NotVertex * pmce = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( mce );
	// link 'cll' to all meshes above 'this->cell_enclosed' (of all dimensions)
	link_face_to_higher ( seg, pmce, 1, 0 );

	Cell A = seg->base_attr.reverse ( tag::surely_exists );
	Cell::Positive::Vertex * Ap = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	link_face_to_msh ( Ap, seg, msh, 0, 1 );
	link_face_to_higher ( Ap, pmce, 0, 1 );
	Cell B = seg->tip_attr;
	Cell::Positive::Vertex * Bp = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	link_face_to_msh ( Bp, seg, msh, 1, 0 );
	link_face_to_higher ( Bp, pmce, 1, 0 );

} // end of make_deep_connections_1d with tag::mesh_is_bdry


inline void make_deep_connections_1d  // hidden in anonymous namespace
(	Cell::Positive::Segment * const seg, Mesh::Core * const msh,
  const tag::MeshIsBdry &, const tag::DoNotBother &           )

// make far connections when adding a positive cell
// see paragraph 11.9 in the manual
	
// tag::do_not_bother is only relevant for Mesh::Connected::***Dim

{	assert ( msh->get_dim_plus_one() == 2 );
	assert ( seg );
	assert ( seg->get_dim() == 1 );

	add_link_same_dim ( seg, msh, tag::do_not_bother );

	// for all meshes strictly above msh
	Cell::Core * mce = msh->cell_enclosed;
	assert ( mce );
	Cell::Positive::NotVertex * pmce = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( mce );
	// link 'cll' to all meshes above 'this->cell_enclosed' (of all dimensions)
	link_face_to_higher ( seg, pmce, 1, 0, tag::do_not_bother );

	Cell A = seg->base_attr.reverse ( tag::surely_exists );
	Cell::Positive::Vertex * Ap = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	link_face_to_msh ( Ap, seg, msh, 0, 1, tag::do_not_bother );
	link_face_to_higher ( Ap, pmce, 0, 1, tag::do_not_bother );
	Cell B = seg->tip_attr;
	Cell::Positive::Vertex * Bp = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	link_face_to_msh ( Bp, seg, msh, 1, 0, tag::do_not_bother );
	link_face_to_higher ( Bp, pmce, 1, 0, tag::do_not_bother );

} // end of make_deep_connections_1d with tag::mesh_is_bdry, tag::do_not_bother


inline void make_deep_connections_1d_rev  // hidden in anonymous namespace
(	Cell::Core * const o_seg, Cell::Positive::Segment * const seg,
  Mesh::Core * const msh, const tag::MeshIsNotBdry &             )

// make far connections when adding a negative cell
// see paragraph 11.9 in the manual
	
{	assert ( seg != o_seg );
	assert ( seg );  assert ( o_seg );
	assert ( msh->get_dim_plus_one() == 2 );
	assert ( seg->get_dim() == 1 );
	assert ( o_seg->get_dim() == 1 );
	// we could assert  msh->cell_enclosed == nullptr  (i.e. mesh is not a boundary)
	// however, we sometimes call this method on the boundary of a
	// freshly created cell (pretending it is not a boundary)

	add_link_same_dim_rev ( o_seg, seg, msh );
	
	Cell A = seg->base_attr.reverse ( tag::surely_exists );
	Cell::Positive::Vertex * Ap = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	link_face_to_msh ( Ap, seg, msh, 0, 1 );
	Cell B = seg->tip_attr;
	Cell::Positive::Vertex * Bp = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	link_face_to_msh ( Bp, seg, msh, 1, 0 );

} // end of make_deep_connections_1d_rev with tag::mesh_is_not_bdry


inline void make_deep_connections_1d_rev  // hidden in anonymous namespace
(	Cell::Core * const o_seg, Cell::Positive::Segment * const seg, Mesh::Core * const msh,
  const tag::MeshIsNotBdry &, const tag::DoNotBother &                                  )

// make far connections when adding a negative cell
// see paragraph 11.9 in the manual
	
// tag::do_not_bother is only relevant for Mesh::Connected::***Dim

{	assert ( seg != o_seg );
	assert ( seg );  assert ( o_seg );
	assert ( msh->get_dim_plus_one() == 2 );
	assert ( seg->get_dim() == 1 );
	assert ( o_seg->get_dim() == 1 );
	// we could assert  msh->cell_enclosed == nullptr  (i.e. mesh is not a boundary)
	// however, we sometimes call this method on the boundary of a
	// freshly created cell (pretending it is not a boundary)

	add_link_same_dim_rev ( o_seg, seg, msh, tag::do_not_bother );
	
	Cell A = seg->base_attr.reverse ( tag::surely_exists );
	Cell::Positive::Vertex * Ap = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	link_face_to_msh ( Ap, seg, msh, 0, 1, tag::do_not_bother );
	Cell B = seg->tip_attr;
	Cell::Positive::Vertex * Bp = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	link_face_to_msh ( Bp, seg, msh, 1, 0, tag::do_not_bother );

} // end of make_deep_connections_1d_rev with tag::mesh_is_not_bdry, tag::do_not_bother


inline void make_deep_connections_1d_rev  // hidden in anonymous namespace
(	Cell::Core * o_seg, Cell::Positive::Segment * seg,
  Mesh::Core * msh, const tag::MeshIsBdry &          )

// make far connections when adding a negative cell
// see paragraph 11.9 in the manual
	
{	assert ( seg != o_seg );
	assert ( seg );  assert ( o_seg );
	assert ( msh->get_dim_plus_one() == 2 );
	assert ( seg->get_dim() == 1 );
	assert ( o_seg->get_dim() == 1 );

	add_link_same_dim_rev ( o_seg, seg, msh );
	
	// for all meshes strictly above msh
	Cell::Core * mce = msh->cell_enclosed;
	assert ( mce );
	Cell::Positive::NotVertex * pmce = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( mce );
	// link 'cll' to all meshes above 'this->cell_enclosed' (of all dimensions)
	link_face_to_higher ( seg, pmce, 0, 1 );  // we switch the two counters

	Cell A = seg->base_attr.reverse ( tag::surely_exists );
	Cell::Positive::Vertex * Ap = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	link_face_to_msh ( Ap, seg, msh, 0, 1 );
	link_face_to_higher ( Ap, pmce, 1, 0 );  // we switch the two counters
	Cell B = seg->tip_attr;
	Cell::Positive::Vertex * Bp = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	link_face_to_msh ( Bp, seg, msh, 1, 0 );
	link_face_to_higher ( Bp, pmce, 0, 1 );  // we switch the two counters

} // end of make_deep_connections_1d_rev with tag::mesh_is_bdry


inline void make_deep_connections_1d_rev  // hidden in anonymous namespace
(	Cell::Core * o_seg, Cell::Positive::Segment * seg, Mesh::Core * msh,
  const tag::MeshIsBdry &, const tag::DoNotBother &                   )

// make far connections when adding a negative cell
// see paragraph 11.9 in the manual
	
// tag::do_not_bother is only relevant for Mesh::Connected::***Dim

{	assert ( seg != o_seg );
	assert ( seg );  assert ( o_seg );
	assert ( msh->get_dim_plus_one() == 2 );
	assert ( seg->get_dim() == 1 );
	assert ( o_seg->get_dim() == 1 );

	add_link_same_dim_rev ( o_seg, seg, msh, tag::do_not_bother );
	
	// for all meshes strictly above msh
	Cell::Core * mce = msh->cell_enclosed;
	assert ( mce );
	Cell::Positive::NotVertex * pmce = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( mce );
	// link 'cll' to all meshes above 'this->cell_enclosed' (of all dimensions)
	link_face_to_higher ( seg, pmce, 0, 1, tag::do_not_bother );  // we switch the two counters

	Cell A = seg->base_attr.reverse ( tag::surely_exists );
	Cell::Positive::Vertex * Ap = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	link_face_to_msh ( Ap, seg, msh, 0, 1, tag::do_not_bother );
	link_face_to_higher ( Ap, pmce, 1, 0, tag::do_not_bother );  // we switch the two counters
	Cell B = seg->tip_attr;
	Cell::Positive::Vertex * Bp = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	link_face_to_msh ( Bp, seg, msh, 1, 0, tag::do_not_bother );
	link_face_to_higher ( Bp, pmce, 0, 1, tag::do_not_bother );  // we switch the two counters

} // end of make_deep_connections_1d_rev with tag::mesh_is_bdry, tag::do_not_bother


inline void make_deep_connections_hd  // hidden in anonymous namespace
(	Cell::Positive::HighDim * const cll, Mesh::Core * const msh, const tag::MeshIsNotBdry & )

// make far connections when adding a positive cell
// see paragraph 11.9 in the manual
	
{	assert ( cll );

	assert ( msh->get_dim_plus_one() > 1 );
	// make_deep_connections_0d deals with the case of a vertex cll
	// make_deep_connections_1d deals with the case of a segment cll
	size_t cll_dim = cll->get_dim();
	size_t cll_dim_m1 = tag::Util::assert_diff ( cll_dim, 1 );
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );
	assert ( cll );
	// we could assert  msh->cell_enclosed == nullptr  (i.e. mesh is not a boundary)
	// however, we sometimes call this method on the boundary of a
	// freshly created cell (pretending it is not a boundary)

	add_link_same_dim ( cll, msh );

	Mesh cll_bdry = cll->boundary();
	Mesh::Iterator it = cll_bdry .iterator
		( tag::over_cells, tag::of_max_dim, tag::force_positive );
	// talvez implementar um iterador especial que devolva cp e cn
	for ( it .reset(); it .in_range(); it ++ )
	{	Cell::Core * face = ( * it ) .core;  // add link from face to 'msh'
		short int cp, cn;
		face->compute_sign ( cp, cn, cll_bdry .core );
		assert ( ( ( cp == 1 ) and ( cn == 0 ) ) or ( ( cp == 0 ) and ( cn == 1 ) ) );
	  link_face_to_msh ( face, cll, msh, cp, cn );                                   }

	for ( size_t d = 0; d < cll_dim_m1 ; d++ )
	{	Mesh::Iterator itt = cll_bdry .iterator  // as they are : positive
			( tag::over_cells, tag::of_dim, d, tag::as_they_are );
		// talvez implementar um iterador especial que devolva cp e cn
		for ( itt .reset(); itt .in_range(); itt ++ )
		{	Cell::Core * fface = ( * itt ) .core;  // add link from face to 'msh'
			short int cp, cn;
			compute_cp_cn ( cp, cn, fface, cll_bdry .core );
			link_face_to_msh ( fface, cll, msh, cp, cn );   }     }

} // end of make_deep_connections_hd with tag::mesh_is_not_bdry


inline void make_deep_connections_hd  // hidden in anonymous namespace
(	Cell::Positive::HighDim * const cll, Mesh::Core * const msh,
  const tag::MeshIsNotBdry &, const tag::DoNotBother &        )

// make far connections when adding a positive cell
// see paragraph 11.9 in the manual
	
// tag::do_not_bother is only relevant for Mesh::Connected::***Dim

{	assert ( cll );
	assert ( msh->get_dim_plus_one() > 1 );
	// make_deep_connections_0d deals with the case of a vertex cll
	// make_deep_connections_1d deals with the case of a segment cll
	size_t cll_dim = cll->get_dim();
	size_t cll_dim_m1 = tag::Util::assert_diff ( cll_dim, 1 );
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );
	assert ( cll );
	// we could assert  msh->cell_enclosed == nullptr  (i.e. mesh is not a boundary)
	// however, we sometimes call this method on the boundary of a
	// freshly created cell (pretending it is not a boundary)

	add_link_same_dim ( cll, msh, tag::do_not_bother );

	Mesh cll_bdry = cll->boundary();
	Mesh::Iterator it = cll_bdry.iterator
		( tag::over_cells, tag::of_max_dim, tag::force_positive );
	// talvez implementar um iterador especial que devolva cp e cn
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face = (*it).core;  // add link from face to 'msh'
		short int cp, cn;
		face->compute_sign ( cp, cn, cll_bdry.core );
		assert ( ( ( cp == 1 ) and ( cn == 0 ) ) or ( ( cp == 0 ) and ( cn == 1 ) ) );
	  link_face_to_msh ( face, cll, msh, cp, cn, tag::do_not_bother );               }

	for ( size_t d = 0; d < cll_dim_m1 ; d++ )
	{	Mesh::Iterator itt = cll_bdry.iterator  // as they are : positive
			( tag::over_cells, tag::of_dim, d, tag::as_they_are );
		// talvez implementar um iterador especial que devolva cp e cn
		for ( itt.reset(); itt.in_range(); itt++ )
		{	Cell::Core * fface = (*itt).core;  // add link from face to 'msh'
			short int cp, cn;
			compute_cp_cn ( cp, cn, fface, cll_bdry.core );
			link_face_to_msh ( fface, cll, msh, cp, cn, tag::do_not_bother );  }  }

} // end of make_deep_connections_hd with tag::mesh_is_not_bdry, tag::do_not_bother


inline void make_deep_connections_hd  // hidden in anonymous namespace
(	Cell::Positive::HighDim * const cll, Mesh::Core * const msh, const tag::MeshIsBdry & )

// make far connections when adding a positive cell
// see paragraph 11.9 in the manual
	
{	assert ( cll );
	assert ( msh->get_dim_plus_one() > 2 );
	// make_deep_connections_0d deals with the case of a vertex cll
	// make_deep_connections_1d deals with the case of a segment cll
	size_t cll_dim = cll->get_dim();
	size_t cll_dim_m1 = tag::Util::assert_diff ( cll_dim, 1 );
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );

	add_link_same_dim ( cll, msh );
	
	// for all meshes strictly above msh
	Cell::Core * mce = msh->cell_enclosed;
	assert ( mce );
	Cell::Positive::NotVertex * pmce = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( mce );
	// link 'cll' to all meshes above 'this->cell_enclosed' (of all dimensions)
	link_face_to_higher ( cll, pmce, 1, 0 );

	Mesh cll_bdry = cll->boundary();
	Mesh::Iterator it = cll_bdry.iterator
		( tag::over_cells, tag::of_max_dim, tag::force_positive );
	// talvez implementar um iterador especial que devolva cp e cn
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face = (*it).core;  // add link from face to 'msh' and
		// to all meshes above 'msh->cell_enclosed' (of all dimensions)
		short int cp, cn;
		face->compute_sign ( cp, cn, cll_bdry.core );
		assert ( ( ( cp == 1 ) and ( cn == 0 ) ) or ( ( cp == 0 ) and ( cn == 1 ) ) );
	  link_face_to_msh ( face, cll, msh, cp, cn );
	  link_face_to_higher ( face, pmce, cp, cn );                                    }

	for ( size_t d = 0; d < cll_dim_m1 ; d++ )
	{	Mesh::Iterator itt = cll_bdry.iterator  // as they are : positive
			( tag::over_cells, tag::of_dim, d, tag::as_they_are );
		// talvez implementar um iterador especial que devolva cp e cn
		for ( itt.reset(); itt.in_range(); itt++ )
		{	Cell::Core * fface = (*itt).core;  // add link from face to 'msh' and
			// to all meshes above 'msh->cell_enclosed' (of all dimensions)
			short int cp, cn;
			compute_cp_cn ( cp, cn, fface, cll_bdry.core );
			link_face_to_msh ( fface, cll, msh, cp, cn );
			link_face_to_higher ( fface, pmce, cp, cn );            }  }

} // end of make_deep_connections_hd with tag::mesh_is_bdry


inline void make_deep_connections_hd  // hidden in anonymous namespace
(	Cell::Positive::HighDim * const cll, Mesh::Core * const msh,
  const tag::MeshIsBdry &, const tag::DoNotBother &           )

// make far connections when adding a positive cell
// see paragraph 11.9 in the manual
	
// tag::do_not_bother is only relevant for Mesh::Connected::***Dim

{	assert ( cll );
	assert ( msh->get_dim_plus_one() > 2 );
	// make_deep_connections_0d deals with the case of a vertex cll
	// make_deep_connections_1d deals with the case of a segment cll
	size_t cll_dim = cll->get_dim();
	size_t cll_dim_m1 = tag::Util::assert_diff ( cll_dim, 1 );
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );

	add_link_same_dim ( cll, msh, tag::do_not_bother );
	
	// for all meshes strictly above msh
	Cell::Core * mce = msh->cell_enclosed;
	assert ( mce );
	Cell::Positive::NotVertex * pmce = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( mce );
	// link 'cll' to all meshes above 'this->cell_enclosed' (of all dimensions)
	link_face_to_higher ( cll, pmce, 1, 0, tag::do_not_bother );

	Mesh cll_bdry = cll->boundary();
	Mesh::Iterator it = cll_bdry.iterator
		( tag::over_cells, tag::of_max_dim, tag::force_positive );
	// talvez implementar um iterador especial que devolva cp e cn
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face = (*it).core;  // add link from face to 'msh' and
		// to all meshes above 'msh->cell_enclosed' (of all dimensions)
		short int cp, cn;
		face->compute_sign ( cp, cn, cll_bdry.core );
		assert ( ( ( cp == 1 ) and ( cn == 0 ) ) or ( ( cp == 0 ) and ( cn == 1 ) ) );
	  link_face_to_msh ( face, cll, msh, cp, cn, tag::do_not_bother );
	  link_face_to_higher ( face, pmce, cp, cn, tag::do_not_bother );                }

	for ( size_t d = 0; d < cll_dim_m1 ; d++ )
	{	Mesh::Iterator itt = cll_bdry.iterator  // as they are : positive
			( tag::over_cells, tag::of_dim, d, tag::as_they_are );
		// talvez implementar um iterador especial que devolva cp e cn
		for ( itt.reset(); itt.in_range(); itt++ )
		{	Cell::Core * fface = (*itt).core;  // add link from face to 'msh' and
			// to all meshes above 'msh->cell_enclosed' (of all dimensions)
			short int cp, cn;
			compute_cp_cn ( cp, cn, fface, cll_bdry.core );
			link_face_to_msh ( fface, cll, msh, cp, cn, tag::do_not_bother );
			link_face_to_higher ( fface, pmce, cp, cn, tag::do_not_bother );  }  }

} // end of make_deep_connections_hd with tag::mesh_is_bdry, tag::do_not_bother


inline void make_deep_connections_hd_rev  // hidden in anonymous namespace
(	Cell::Core * const o_cll, Cell::Positive::HighDim * const cll,
  Mesh::Core * const msh, const tag::MeshIsNotBdry &             )

// make far connections when adding a negative cell
// see paragraph 11.9 in the manual
	
{	assert ( cll != o_cll );
	assert ( cll );  assert ( o_cll );  assert ( msh );
	assert ( msh->get_dim_plus_one() > 2 );
	// make_deep_connections_0d deals with the case of a vertex cll
	// make_deep_connections_1d deals with the case of a segment cll
	size_t cll_dim = cll->get_dim();
	size_t cll_dim_m1 = tag::Util::assert_diff ( cll_dim, 1 );
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );
	// we could assert  msh->cell_enclosed == nullptr  (i.e. mesh is not a boundary)
	// however, we sometimes call this method on the boundary of a
	// freshly created cell (pretending it is not a boundary)

	add_link_same_dim_rev ( o_cll, cll, msh );
	
	Mesh cll_bdry = cll->boundary();
	Mesh::Iterator it = cll_bdry.iterator
		( tag::over_cells, tag::of_max_dim, tag::force_positive );
	// talvez implementar um iterador especial que devolva cp e cn
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face = (*it).core;  // add link from face to 'msh'
		short int cp, cn;
		face->compute_sign ( cp, cn, cll_bdry.core );
		assert ( ( ( cp == 1 ) and ( cn == 0 ) ) or ( ( cp == 0 ) and ( cn == 1 ) ) );
	  link_face_to_msh ( face, cll, msh, cp, cn );                                   }

	for ( size_t d = 0; d < cll_dim_m1 ; d++ )
	{	Mesh::Iterator itt = cll_bdry.iterator  // as they are : positive
			( tag::over_cells, tag::of_dim, d, tag::as_they_are );
		// talvez implementar um iterador especial que devolva cp e cn
		for ( itt.reset(); itt.in_range(); itt++ )
		{	Cell::Core * fface = (*itt).core;  // add link from face to 'msh'
			short int cp, cn;
			compute_cp_cn ( cp, cn, fface, cll_bdry.core );
			link_face_to_msh ( fface, cll, msh, cp, cn );   }       }

} // end of make_deep_connections_hd_rev with tag::mesh_is_not_bdry


inline void make_deep_connections_hd_rev  // hidden in anonymous namespace
(	Cell::Core * const o_cll, Cell::Positive::HighDim * const cll, Mesh::Core * const msh,
  const tag::MeshIsNotBdry &, const tag::DoNotBother &                                  )

// make far connections when adding a negative cell
// see paragraph 11.9 in the manual
	
{	assert ( cll != o_cll );
	assert ( cll );  assert ( o_cll );  assert ( msh );
	assert ( msh->get_dim_plus_one() > 2 );
	// make_deep_connections_0d deals with the case of a vertex cll
	// make_deep_connections_1d deals with the case of a segment cll
	size_t cll_dim = cll->get_dim();
	size_t cll_dim_m1 = tag::Util::assert_diff ( cll_dim, 1 );
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );
	// we could assert  msh->cell_enclosed == nullptr  (i.e. mesh is not a boundary)
	// however, we sometimes call this method on the boundary of a
	// freshly created cell (pretending it is not a boundary)

	add_link_same_dim_rev ( o_cll, cll, msh, tag::do_not_bother );
	
	Mesh cll_bdry = cll->boundary();
	Mesh::Iterator it = cll_bdry.iterator
		( tag::over_cells, tag::of_max_dim, tag::force_positive );
	// talvez implementar um iterador especial que devolva cp e cn
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face = (*it).core;  // add link from face to 'msh'
		short int cp, cn;
		face->compute_sign ( cp, cn, cll_bdry.core );
		assert ( ( ( cp == 1 ) and ( cn == 0 ) ) or ( ( cp == 0 ) and ( cn == 1 ) ) );
	  link_face_to_msh ( face, cll, msh, cp, cn, tag::do_not_bother );              }

	for ( size_t d = 0; d < cll_dim_m1 ; d++ )
	{	Mesh::Iterator itt = cll_bdry.iterator  // as they are : positive
			( tag::over_cells, tag::of_dim, d, tag::as_they_are );
		// talvez implementar um iterador especial que devolva cp e cn
		for ( itt.reset(); itt.in_range(); itt++ )
		{	Cell::Core * fface = (*itt).core;  // add link from face to 'msh'
			short int cp, cn;
			compute_cp_cn ( cp, cn, fface, cll_bdry.core );
			link_face_to_msh ( fface, cll, msh, cp, cn, tag::do_not_bother );  }  }

} // end of make_deep_connections_hd_rev with tag::mesh_is_not_bdry, tag::do_not_bother


inline void make_deep_connections_hd_rev  // hidden in anonymous namespace
(	Cell::Core * o_cll, Cell::Positive::HighDim * cll,
  Mesh::Core * msh, const tag::MeshIsBdry &          )

// make far connections when adding a negative cell
// see paragraph 11.9 in the manual
	
{	assert ( cll );  assert ( o_cll );  assert ( msh );
	assert ( cll != o_cll );
	assert ( msh->get_dim_plus_one() > 1 );
	// make_deep_connections_0d deals with the case of a vertex cll
	// make_deep_connections_1d deals with the case of a segment cll
	size_t cll_dim = cll->get_dim();
	size_t cll_dim_m1 = tag::Util::assert_diff ( cll_dim, 1 );
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );

	add_link_same_dim_rev ( o_cll, cll, msh );
	
	// for all meshes strictly above msh
	Cell::Core * mce = msh->cell_enclosed;
	assert ( mce );
	Cell::Positive::NotVertex * pmce = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( mce );
	// link 'cll' to all meshes above 'this->cell_enclosed' (of all dimensions)
	link_face_to_higher ( cll, pmce, 0, 1 );  // we switch the two counters

	Mesh cll_bdry = cll->boundary();
	Mesh::Iterator it = cll_bdry.iterator
		( tag::over_cells, tag::of_max_dim, tag::force_positive );
	// talvez implementar um iterador especial que devolva cp e cn
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face = (*it).core;  // add link from face to 'msh' and
		// to all meshes above 'msh->cell_enclosed' (of all dimensions)
		short int cp, cn;
		face->compute_sign ( cp, cn, cll_bdry.core );
		assert ( ( ( cp == 1 ) and ( cn == 0 ) ) or ( ( cp == 0 ) and ( cn == 1 ) ) );
	  link_face_to_msh ( face, cll, msh, cp, cn );
		// we switch the two counters
	  link_face_to_higher ( face, pmce, cn, cp );                                    }

	for ( size_t d = 0; d < cll_dim_m1 ; d++ )
	{	Mesh::Iterator itt = cll_bdry.iterator  // as they are : positive
			( tag::over_cells, tag::of_dim, d, tag::as_they_are );
		// talvez implementar um iterador especial que devolva cp e cn
		for ( itt.reset(); itt.in_range(); itt++ )
		{	Cell::Core * fface = (*itt).core;  // add link from face to 'msh' and
			// to all meshes above 'msh->cell_enclosed' (of all dimensions)
			short int cp, cn;
			compute_cp_cn ( cp, cn, fface, cll_bdry.core );
			link_face_to_msh ( fface, cll, msh, cp, cn );
			// we switch the two counters
			link_face_to_higher ( fface, pmce, cn, cp );    }  }

} // end of make_deep_connections_hd_rev with tag::mesh_is_bdry


inline void make_deep_connections_hd_rev  // hidden in anonymous namespace
(	Cell::Core * o_cll, Cell::Positive::HighDim * cll, Mesh::Core * msh,
  const tag::MeshIsBdry &, const tag::DoNotBother                      )

// make far connections when adding a negative cell
// see paragraph 11.9 in the manual
	
// tag::do_not_bother is only relevant for Mesh::Connected::***Dim

{	assert ( cll );  assert ( o_cll );  assert ( msh );
	assert ( cll != o_cll );
	assert ( msh->get_dim_plus_one() > 1 );
	// make_deep_connections_0d deals with the case of a vertex cll
	// make_deep_connections_1d deals with the case of a segment cll
	size_t cll_dim = cll->get_dim();
	size_t cll_dim_m1 = tag::Util::assert_diff ( cll_dim, 1 );
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );

	add_link_same_dim_rev ( o_cll, cll, msh, tag::do_not_bother );
	
	// for all meshes strictly above msh
	Cell::Core * mce = msh->cell_enclosed;
	assert ( mce );
	Cell::Positive::NotVertex * pmce = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( mce );
	// link 'cll' to all meshes above 'this->cell_enclosed' (of all dimensions)
	link_face_to_higher ( cll, pmce, 0, 1, tag::do_not_bother );  // we switch the two counters

	Mesh cll_bdry = cll->boundary();
	Mesh::Iterator it = cll_bdry.iterator
		( tag::over_cells, tag::of_max_dim, tag::force_positive );
	// talvez implementar um iterador especial que devolva cp e cn
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face = (*it).core;  // add link from face to 'msh' and
		// to all meshes above 'msh->cell_enclosed' (of all dimensions)
		short int cp, cn;
		face->compute_sign ( cp, cn, cll_bdry.core );
		assert ( ( ( cp == 1 ) and ( cn == 0 ) ) or ( ( cp == 0 ) and ( cn == 1 ) ) );
	  link_face_to_msh ( face, cll, msh, cp, cn, tag::do_not_bother );
		// we switch the two counters
	  link_face_to_higher ( face, pmce, cn, cp, tag::do_not_bother );               }

	for ( size_t d = 0; d < cll_dim_m1 ; d++ )
	{	Mesh::Iterator itt = cll_bdry.iterator  // as they are : positive
			( tag::over_cells, tag::of_dim, d, tag::as_they_are );
		// talvez implementar um iterador especial que devolva cp e cn
		for ( itt.reset(); itt.in_range(); itt++ )
		{	Cell::Core * fface = (*itt).core;  // add link from face to 'msh' and
			// to all meshes above 'msh->cell_enclosed' (of all dimensions)
			short int cp, cn;
			compute_cp_cn ( cp, cn, fface, cll_bdry.core );
			link_face_to_msh ( fface, cll, msh, cp, cn, tag::do_not_bother );
			// we switch the two counters
			link_face_to_higher ( fface, pmce, cn, cp, tag::do_not_bother );  }  }

} // end of make_deep_connections_hd_rev with tag::mesh_is_bdry, tag::do_not_bother


inline void remove_link_zero_dim  // hidden in anonymous namespace
(	Cell::Positive::Vertex * const ver, Cell::Positive::Segment * const seg )

// this function "removes the link" between a vertex and a segment
// this "link" asserts that the vertex belongs to the (boundary of the) segment

// this function updates ver->segments

{	assert ( ver );  assert ( seg );
	typedef std::map < Cell::Positive::Segment *, short int > maptype;
	maptype & vs = ver->segments;
	maptype::iterator vsfs = vs.find(seg);
	assert ( vsfs != vs.end() );
	assert ( vsfs->second == 1 );
	vs.erase ( vsfs );                                                   }


inline void remove_link_zero_dim_rev  // hidden in anonymous namespace
(	Cell::Positive::Vertex * const ver, Cell::Positive::Segment * const seg )

// this function "removes the link" between a vertex and a segment
// this "link" asserts that the vertex belongs to the (boundary of the) segment

// this function updates ver->segments

{	assert ( ver );  assert ( seg );
	typedef std::map < Cell::Positive::Segment *, short int > maptype;
	maptype & vs = ver->segments;
	maptype::iterator vsfs = vs.find(seg);
	assert ( vsfs != vs.end() );
	assert ( vsfs->second == -1 );
	vs.erase ( vsfs );                                                   }


inline void remove_link_same_dim  // hidden in anonymous namespace
(	Cell::Positive::NotVertex * cll, Mesh::Core * msh )

// this function "removes the link" between a cell and a mesh
// this "link" asserts that the cell belongs to the mesh
// here the dimensions are equal and the cell is positive

// this function updates cll->meshes_same_dim
// also, calls the virtual method msh->remove_from_my_cells
// which, for fuzzy and stsi meshes, updates msh->cells
// removing cll with the aid of the iterator provided
// for other kinds of meshes, does nothing

{	assert ( cll );  assert ( msh );
	assert ( cll->is_positive() ); // assert ( msh->is_positive() );
	size_t cll_dim = cll->get_dim();
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );
	typedef std::map <Mesh::Core*, Cell::field_to_meshes_same_dim> maptype;
	maptype & cmd = cll->meshes_same_dim;
	maptype::iterator cmdfm = cmd.find(msh);
	assert ( cmdfm != cmd.end() );
	assert ( cmdfm->second.sign == 1 );
	msh->remove_from_my_cells ( cll, cll_dim, cmdfm->second.where );
	cmd.erase ( cmdfm );

} // end of remove_link_same_dim


inline void remove_link_same_dim  // hidden in anonymous namespace
(	Cell::Positive::NotVertex * cll, Mesh::Core * msh, const tag::DoNotBother & )

// this function "removes the link" between a cell and a mesh
// this "link" asserts that the cell belongs to the mesh
// here the dimensions are equal and the cell is positive

// this function updates cll->meshes_same_dim
// also, calls the virtual method msh->remove_from_my_cells
// which, for fuzzy and stsi meshes, updates msh->cells
// removing cll with the aid of the iterator provided
// for other kinds of meshes, does nothing

{	assert ( cll );  assert ( msh );
	assert ( cll->is_positive() ); // assert ( msh->is_positive() );
	size_t cll_dim = cll->get_dim();
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );
	typedef std::map <Mesh::Core*, Cell::field_to_meshes_same_dim> maptype;
	maptype & cmd = cll->meshes_same_dim;
	maptype::iterator cmdfm = cmd.find(msh);
	assert ( cmdfm != cmd.end() );
	assert ( cmdfm->second.sign == 1 );
	msh->remove_from_my_cells ( cll, cll_dim, cmdfm->second.where, tag::do_not_bother );
	cmd.erase ( cmdfm );

} // end of remove_link_same_dim with tag::do_not_bother


inline void remove_link_same_dim_rev  // hidden in anonymous namespace
(	Cell::Core * o_cll, Cell::Positive::NotVertex * cll, Mesh::Core * msh )

// this function "removes the link" between a cell and a mesh
// this "link" asserts that the cell belongs to the mesh
// here the dimensions are equal and the cell is negative

// this function updates cll->meshes_same_dim
// also, calls the virtual method msh->remove_from_my_cells
// which, for fuzzy and stsi meshes, updates msh->cells
// removing cll with the aid of the iterator provided
// for other kinds of meshes, does nothing

{	assert ( cll );  assert ( msh );
	assert ( cll->is_positive() ); // assert ( msh->is_positive() );
	size_t cll_dim = cll->get_dim();
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );
	typedef std::map <Mesh::Core*, Cell::field_to_meshes_same_dim> maptype;
	maptype & cmd = cll->meshes_same_dim;
	maptype::iterator cmdfm = cmd.find(msh);
	assert ( cmdfm != cmd.end() );
	assert ( cmdfm->second.sign == -1 );
	msh->remove_from_my_cells ( o_cll, cll_dim, cmdfm->second.where );
	cmd.erase ( cmdfm );

} // end of remove_link_same_dim


inline void remove_link_same_dim_rev  // hidden in anonymous namespace
(	Cell::Core * o_cll, Cell::Positive::NotVertex * cll, Mesh::Core * msh,
  const tag::DoNotBother &                                              )

// this function "removes the link" between a cell and a mesh
// this "link" asserts that the cell belongs to the mesh
// here the dimensions are equal and the cell is negative

// this function updates cll->meshes_same_dim
// also, calls the virtual method msh->remove_from_my_cells
// which, for fuzzy and stsi meshes, updates msh->cells
// removing cll with the aid of the iterator provided
// for other kinds of meshes, does nothing

{	assert ( cll );  assert ( msh );
	assert ( cll->is_positive() ); // assert ( msh->is_positive() );
	size_t cll_dim = cll->get_dim();
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );
	typedef std::map <Mesh::Core*, Cell::field_to_meshes_same_dim> maptype;
	maptype & cmd = cll->meshes_same_dim;
	maptype::iterator cmdfm = cmd.find(msh);
	assert ( cmdfm != cmd.end() );
	assert ( cmdfm->second.sign == -1 );
	msh->remove_from_my_cells ( o_cll, cll_dim, cmdfm->second.where, tag::do_not_bother );
	cmd.erase ( cmdfm );

} // end of remove_link_same_dim with tag::do_not_bother


inline void remove_link  // hidden in anonymous namespace
(	Cell::Positive * const cll, Mesh::Core * const msh, const short int cp, const short int cn )

// this function "removes the link" between a cell and a mesh
// this "link" asserts that the cell belongs to the mesh
// here 'cll' has lower dimension, so it may appear several times in the mesh 'msh',
// with both orientations
// the counter 'cp' asserts how many times 'cll' appears
// in the mesh with positive orientation,
// while 'cn' says the same thing for negative orientation

// this function updates cll->meshes
// also, calls the virtual method msh->remove_from_my_cells
// which, for fuzzy and stsi meshes, updates msh->cells
// removing cll with the aid of the iterator provided
// for other kinds of meshes, does nothing

{	assert ( cll );  assert ( msh );
	assert ( cll->is_positive() ); // assert ( msh->is_positive() );
	size_t cll_dim = cll->get_dim(),
	       msh_dim_p1 = msh->get_dim_plus_one();
	assert ( msh_dim_p1 > cll_dim + 1 );
	size_t dif_dim = msh_dim_p1 - cll_dim - 1;
	assert ( cll->meshes.size() > dif_dim );
	typedef std::map <Mesh::Core*, Cell::field_to_meshes> maptype;
	maptype & cmd = cll->meshes[dif_dim];
	maptype::iterator cmdfm = cmd.find(msh);
	assert ( cmdfm != cmd.end() );
	short int c_p = cmdfm->second.counter_pos -= cp;
	short int c_n = cmdfm->second.counter_neg -= cn;

	assert ( ( c_p >= 0 ) and ( c_n >= 0 ) );
	if ( ( c_p == 0 ) and ( c_n == 0 ) )
	{	msh->remove_from_my_cells ( cll, cll_dim, cmdfm->second.where );
		cmd.erase ( cmdfm );                                             }

} // end of remove_link


inline void remove_link  // hidden in anonymous namespace
( Cell::Positive * const cll, Mesh::Core * const msh,
  const short int cp, const short int cn, const tag::DoNotBother & )

// this function "removes the link" between a cell and a mesh
// this "link" asserts that the cell belongs to the mesh
// here 'cll' has lower dimension, so it may appear several times in the mesh 'msh',
// with both orientations
// the counter 'cp' asserts how many times 'cll' appears
// in the mesh with positive orientation,
// while 'cn' says the same thing for negative orientation

// this function updates cll->meshes
// also, calls the virtual method msh->remove_from_my_cells
// which, for fuzzy and stsi meshes, updates msh->cells
// removing cll with the aid of the iterator provided
// for other kinds of meshes, does nothing

{	assert ( cll );  assert ( msh );
	assert ( cll->is_positive() ); // assert ( msh->is_positive() );
	size_t cll_dim = cll->get_dim(),
	       msh_dim_p1 = msh->get_dim_plus_one();
	assert ( msh_dim_p1 > cll_dim + 1 );
	size_t dif_dim = msh_dim_p1 - cll_dim - 1;
	assert ( cll->meshes.size() > dif_dim );
	typedef std::map <Mesh::Core*, Cell::field_to_meshes> maptype;
	maptype & cmd = cll->meshes[dif_dim];
	maptype::iterator cmdfm = cmd.find(msh);
	assert ( cmdfm != cmd.end() );
	short int c_p = cmdfm->second.counter_pos -= cp;
	short int c_n = cmdfm->second.counter_neg -= cn;
	assert ( ( c_p >= 0 ) and ( c_n >= 0 ) );
	if ( ( c_p == 0 ) and ( c_n == 0 ) )
		{	msh->remove_from_my_cells ( cll, cll_dim, cmdfm->second.where, tag::do_not_bother );
		cmd.erase ( cmdfm );                                                                   }

} // end of remove_link with tag::do_not_bother


inline void unlink_face_from_msh  // hidden in anonymous namespace
(	Cell::Core * const face, Cell::Positive::NotVertex * cll,
  Mesh::Core * const msh, const short int cp, const short int cn )
// just a block of code for break_deep_connections

{	assert ( face );  assert ( cll );
	assert ( face->get_dim() < cll->get_dim() );
	assert ( cll->get_dim() + 1 == msh->get_dim_plus_one() );
	Cell::Positive * face_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive* > ( face );
	// size_t dif_dim = msh->get_dim_plus_one() - face->get_dim() - 1;
	assert ( msh->get_dim_plus_one() < face_p->meshes.size() + face->get_dim() + 1 );
	typedef std::map<Mesh::Core*,Cell::field_to_meshes_same_dim> maptype;
	maptype & cemd = cll->meshes_same_dim;
	maptype::iterator map_iter = cemd.find ( msh );
	assert ( map_iter != cemd.end() );
	short int mis = map_iter->second.sign;
	if ( mis == 1 ) remove_link ( face_p, msh, cp, cn );
	else  // we switch the two counters
	{	assert ( mis == -1 ); remove_link ( face_p, msh, cn, cp );  }       }


inline void unlink_face_from_msh  // hidden in anonymous namespace
(	Cell::Core * const face, Cell::Positive::NotVertex * cll,
  Mesh::Core * const msh, const short int cp, const short int cn,
  const tag::DoNotBother &                                        )
// just a block of code for break_deep_connections

{	assert ( face );  assert ( cll );
	assert ( face->get_dim() < cll->get_dim() );
	assert ( cll->get_dim() + 1 == msh->get_dim_plus_one() );
	Cell::Positive * face_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive* > ( face );
	// size_t dif_dim = msh->get_dim_plus_one() - face->get_dim() - 1;
	assert ( msh->get_dim_plus_one() < face_p->meshes.size() + face->get_dim() + 1 );
	typedef std::map<Mesh::Core*,Cell::field_to_meshes_same_dim> maptype;
	maptype & cemd = cll->meshes_same_dim;
	maptype::iterator map_iter = cemd.find ( msh );
	assert ( map_iter != cemd.end() );
	short int mis = map_iter->second.sign;
	if ( mis == 1 ) remove_link ( face_p, msh, cp, cn, tag::do_not_bother );
	else  // we switch the two counters
	{	assert ( mis == -1 ); remove_link ( face_p, msh, cn, cp, tag::do_not_bother );  }  }


inline void unlink_face_from_higher  // hidden in anonymous namespace
(	Cell::Core * const face, Cell::Positive::NotVertex * const pmce,
  const short int cp, const short int cn                          )
// just a block of code for break_deep_connections

{	assert ( face );  assert ( pmce );
	assert ( face->get_dim() + 1 <= pmce->get_dim() );
	Cell::Positive * face_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive* > ( face );
	{ // just a block of code for hiding maptype, cemd and map_iter
	// pmce->meshes[0] is empty, we use pmce->meshes_same_dim instead
	typedef std::map<Mesh::Core*,Cell::field_to_meshes_same_dim> maptype;
	maptype & cemd = pmce->meshes_same_dim;
	maptype::iterator map_iter;
	// we now loop over all meshes of given dimension
	for ( map_iter = cemd.begin(); map_iter != cemd.end(); ++map_iter )
	{	assert ( map_iter->first->get_dim_plus_one() == pmce->get_dim() + 1 );
		Cell::field_to_meshes_same_dim & mis = map_iter->second;
		if ( mis.sign == 1 ) remove_link ( face_p, map_iter->first, cp, cn );
		else  // we just switch cp and cn
		{	assert ( mis.sign == -1 );
			remove_link ( face_p, map_iter->first, cn, cp );  }                   }
	} // just a block of code for hiding maptype, cemd and map_iter
	for ( size_t dif_dim = 1; dif_dim < pmce->meshes.size(); dif_dim++ )
	{	typedef std::map<Mesh::Core*,Cell::field_to_meshes> maptype;
		maptype & cemd = pmce->meshes[dif_dim];
		maptype::iterator map_iter;
		// we now loop over all meshes of given dimension
		for ( map_iter = cemd.begin(); map_iter != cemd.end(); ++map_iter )
		{	assert ( map_iter->first->get_dim_plus_one() == dif_dim + pmce->get_dim() + 1 );
			Cell::field_to_meshes & mis = map_iter->second;
			remove_link ( face_p, map_iter->first,
			              cp*mis.counter_pos + cn*mis.counter_neg,
			              cn*mis.counter_pos + cp*mis.counter_neg );                          }  }  }


inline void unlink_face_from_higher  // hidden in anonymous namespace
(	Cell::Core * const face, Cell::Positive::NotVertex * const pmce,
  const short int cp, const short int cn, const tag::DoNotBother & )
// just a block of code for break_deep_connections

{	assert ( face );  assert ( pmce );
	assert ( face->get_dim() + 1 <= pmce->get_dim() );
	Cell::Positive * face_p = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive* > ( face );
	{ // just a block of code for hiding maptype, cemd and map_iter
	// pmce->meshes[0] is empty, we use pmce->meshes_same_dim instead
	typedef std::map<Mesh::Core*,Cell::field_to_meshes_same_dim> maptype;
	maptype & cemd = pmce->meshes_same_dim;
	maptype::iterator map_iter;
	// we now loop over all meshes of given dimension
	for ( map_iter = cemd.begin(); map_iter != cemd.end(); ++map_iter )
	{	assert ( map_iter->first->get_dim_plus_one() == pmce->get_dim() + 1 );
		Cell::field_to_meshes_same_dim & mis = map_iter->second;
		if ( mis.sign == 1 ) remove_link ( face_p, map_iter->first, cp, cn, tag::do_not_bother );
		else  // we just switch cp and cn
		{	assert ( mis.sign == -1 );
			remove_link ( face_p, map_iter->first, cn, cp, tag::do_not_bother );  }                 }
	} // just a block of code for hiding maptype, cemd and map_iter
	for ( size_t dif_dim = 1; dif_dim < pmce->meshes.size(); dif_dim++ )
	{	typedef std::map<Mesh::Core*,Cell::field_to_meshes> maptype;
		maptype & cemd = pmce->meshes[dif_dim];
		maptype::iterator map_iter;
		// we now loop over all meshes of given dimension
		for ( map_iter = cemd.begin(); map_iter != cemd.end(); ++map_iter )
		{	assert ( map_iter->first->get_dim_plus_one() == dif_dim + pmce->get_dim() + 1 );
			Cell::field_to_meshes & mis = map_iter->second;
			remove_link ( face_p, map_iter->first,
			              cp*mis.counter_pos + cn*mis.counter_neg,
			              cn*mis.counter_pos + cp*mis.counter_neg, tag::do_not_bother );      }  }  }


inline void break_deep_connections_0d  // hidden in anonymous namespace
(	Cell::Positive::Vertex * const ver, Cell::Positive::Segment * const seg )

// break far connections when adding a positive vertex
// see paragraph 11.9 in the manual

{	assert ( ver );  assert ( seg );
	unlink_face_from_higher ( ver, seg, 1, 0 );
	remove_link_zero_dim ( ver, seg );          }
	

inline void break_deep_connections_0d  // hidden in anonymous namespace
(	Cell::Positive::Vertex * const ver, Cell::Positive::Segment * const seg,
  const tag::DoNotBother &                                                )

// break far connections when adding a positive vertex
// see paragraph 11.9 in the manual

{	assert ( ver );  assert ( seg );
	unlink_face_from_higher ( ver, seg, 1, 0, tag::do_not_bother );
	remove_link_zero_dim ( ver, seg );                              }
	

inline void break_deep_connections_0d_rev  // hidden in anonymous namespace
(	Cell::Positive::Vertex * const ver, Cell::Positive::Segment * const seg )

// break far connections when adding a negative vertex
// see paragraph 11.9 in the manual

{	assert ( ver );  assert ( seg );
	unlink_face_from_higher ( ver, seg, 0, 1 );    // we switch the two counters
	remove_link_zero_dim_rev ( ver, seg );      }
	

inline void break_deep_connections_0d_rev  // hidden in anonymous namespace
(	Cell::Positive::Vertex * const ver, Cell::Positive::Segment * const seg,
  const tag::DoNotBother &                                                )

// break far connections when adding a negative vertex
// see paragraph 11.9 in the manual

{	assert ( ver );  assert ( seg );  // we switch the two counters
	unlink_face_from_higher ( ver, seg, 0, 1, tag::do_not_bother );
	remove_link_zero_dim_rev ( ver, seg );                         }
	

inline void break_deep_connections_1d  // hidden in anonymous namespace
(	Cell::Positive::Segment * const seg, Mesh::Core * const msh, const tag::MeshIsNotBdry & )

// break far connections when removing a positive segment
// see paragraph 11.9 in the manual
	
{	assert ( msh->get_dim_plus_one() == 2 );
	assert ( seg );  assert ( seg->get_dim() == 1 );
	// we could assert  msh->cell_enclosed == nullptr  (i.e. mesh is not a boundary)
	// however, we sometimes call this method on the boundary of a
	// cell being destroyed (pretending it is not a boundary)

	Cell A = seg->base_attr.reverse ( tag::surely_exists );
	Cell::Positive::Vertex * Ap = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	unlink_face_from_msh ( Ap, seg, msh, 0, 1 );
	Cell B = seg->tip_attr;
	Cell::Positive::Vertex * Bp = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	unlink_face_from_msh ( Bp, seg, msh, 1, 0 );

	remove_link_same_dim ( seg, msh );

} // end of break_deep_connections_1d with tag::mesh_is_not_bdry


inline void break_deep_connections_1d  // hidden in anonymous namespace
(	Cell::Positive::Segment * const seg, Mesh::Core * const msh,
  const tag::MeshIsNotBdry &, const tag::DoNotBother &        )

// break far connections when removing a positive segment
// see paragraph 11.9 in the manual
	
// tag::do_not_bother is only relevant for Mesh::Connected::***Dim

{	assert ( msh->get_dim_plus_one() == 2 );
	assert ( seg );  assert ( seg->get_dim() == 1 );
	// we could assert  msh->cell_enclosed == nullptr  (i.e. mesh is not a boundary)
	// however, we sometimes call this method on the boundary of a
	// cell being destroyed (pretending it is not a boundary)

	Cell A = seg->base_attr.reverse ( tag::surely_exists );
	Cell::Positive::Vertex * Ap = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	unlink_face_from_msh ( Ap, seg, msh, 0, 1, tag::do_not_bother );
	Cell B = seg->tip_attr;
	Cell::Positive::Vertex * Bp = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	unlink_face_from_msh ( Bp, seg, msh, 1, 0, tag::do_not_bother );

	remove_link_same_dim ( seg, msh, tag::do_not_bother );

} // end of break_deep_connections_1d with tag::mesh_is_not_bdry, tag::do_not_bother


inline void break_deep_connections_1d  // hidden in anonymous namespace
(	Cell::Positive::Segment * const seg, Mesh::Core * const msh, const tag::MeshIsBdry & )

// break far connections when removing a positive segment
// see paragraph 11.9 in the manual
	
{	assert ( msh->get_dim_plus_one() == 2 );
	assert ( seg );
	assert ( seg->get_dim() == 1 );

	// unlink 'seg' from all meshes above 'msh->cell_enclosed' (of all dimensions)
	Cell::Core * mce = msh->cell_enclosed;
	assert ( mce );
	Cell::Positive::NotVertex * pmce = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( mce );
	Cell A = seg->base_attr.reverse ( tag::surely_exists );
	Cell::Positive::Vertex * Ap = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	unlink_face_from_msh ( Ap, seg, msh, 0, 1 );
	unlink_face_from_higher ( Ap, pmce, 0, 1 );
	Cell B = seg->tip_attr;
	Cell::Positive::Vertex * Bp = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	unlink_face_from_msh ( Bp, seg, msh, 1, 0 );
	unlink_face_from_higher ( Bp, pmce, 1, 0 );

	unlink_face_from_higher ( seg, pmce, 1, 0 );

	remove_link_same_dim ( seg, msh );

} // end of break_deep_connections_1d with tag::mesh_is_bdry


inline void break_deep_connections_1d  // hidden in anonymous namespace
(	Cell::Positive::Segment * const seg, Mesh::Core * const msh,
  const tag::MeshIsBdry &, const tag::DoNotBother &           )

// break far connections when removing a positive segment
// see paragraph 11.9 in the manual
	
// tag::do_not_bother is only relevant for Mesh::Connected::***Dim

{	assert ( msh->get_dim_plus_one() == 2 );
	assert ( seg );
	assert ( seg->get_dim() == 1 );

	// unlink 'seg' from all meshes above 'msh->cell_enclosed' (of all dimensions)
	Cell::Core * mce = msh->cell_enclosed;
	assert ( mce );
	Cell::Positive::NotVertex * pmce = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( mce );
	Cell A = seg->base_attr.reverse ( tag::surely_exists );
	Cell::Positive::Vertex * Ap = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	unlink_face_from_msh ( Ap, seg, msh, 0, 1, tag::do_not_bother );
	unlink_face_from_higher ( Ap, pmce, 0, 1, tag::do_not_bother );
	Cell B = seg->tip_attr;
	Cell::Positive::Vertex * Bp = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	unlink_face_from_msh ( Bp, seg, msh, 1, 0, tag::do_not_bother );
	unlink_face_from_higher ( Bp, pmce, 1, 0, tag::do_not_bother );

	unlink_face_from_higher ( seg, pmce, 1, 0, tag::do_not_bother );

	remove_link_same_dim ( seg, msh, tag::do_not_bother );

} // end of break_deep_connections_1d with tag::mesh_is_bdry, tag::do_not_bother


inline void break_deep_connections_1d_rev  // hidden in anonymous namespace
(	Cell::Core * const o_seg, Cell::Positive::Segment * const seg,
  Mesh::Core * const msh, const tag::MeshIsNotBdry &               )

// break far connections when removing a negative cell
// see paragraph 11.9 in the manual
	
{	assert ( seg != o_seg );
	assert ( seg );  assert ( o_seg );
	assert ( msh->get_dim_plus_one() == 2 );
	assert ( seg->get_dim() == 1 );
	assert ( o_seg->get_dim() == 1 );
	// we could assert  msh->cell_enclosed == nullptr  (i.e. mesh is not a boundary)
	// however, we sometimes call this method on the boundary of a
	// cell being destroyed (pretending it is not a boundary)

	Cell A = seg->base_attr.reverse ( tag::surely_exists );
	Cell::Positive::Vertex * Ap = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	unlink_face_from_msh ( Ap, seg, msh, 0, 1 );
	Cell B = seg->tip_attr;
	Cell::Positive::Vertex * Bp = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	unlink_face_from_msh ( Bp, seg, msh, 1, 0 );

	remove_link_same_dim_rev ( o_seg, seg, msh );
	
} // end of break_deep_connections_1d_rev with tag::mesh_is_not_bdry


inline void break_deep_connections_1d_rev  // hidden in anonymous namespace
(	Cell::Core * const o_seg, Cell::Positive::Segment * const seg, Mesh::Core * const msh,
  const tag::MeshIsNotBdry &, const tag::DoNotBother &                                  )

// break far connections when removing a negative cell
// see paragraph 11.9 in the manual
	
// tag::do_not_bother is only relevant for Mesh::Connected::***Dim

{	assert ( seg != o_seg );
	assert ( seg );  assert ( o_seg );
	assert ( msh->get_dim_plus_one() == 2 );
	assert ( seg->get_dim() == 1 );
	assert ( o_seg->get_dim() == 1 );
	// we could assert  msh->cell_enclosed == nullptr  (i.e. mesh is not a boundary)
	// however, we sometimes call this method on the boundary of a
	// cell being destroyed (pretending it is not a boundary)

	Cell A = seg->base_attr.reverse ( tag::surely_exists );
	Cell::Positive::Vertex * Ap = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	unlink_face_from_msh ( Ap, seg, msh, 0, 1, tag::do_not_bother );
	Cell B = seg->tip_attr;
	Cell::Positive::Vertex * Bp = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	unlink_face_from_msh ( Bp, seg, msh, 1, 0, tag::do_not_bother );

	remove_link_same_dim_rev ( o_seg, seg, msh, tag::do_not_bother );
	
} // end of break_deep_connections_1d_rev with tag::mesh_is_not_bdry, tag::do_not_bother


inline void break_deep_connections_1d_rev  // hidden in anonymous namespace
(	Cell::Core * o_seg, Cell::Positive::Segment * seg,
  Mesh::Core * msh, const tag::MeshIsBdry &          )

// break far connections when removing a negative cell
// see paragraph 11.9 in the manual
	
{	assert ( seg != o_seg );
	assert ( seg );  assert ( o_seg );
	assert ( msh->get_dim_plus_one() == 2 );
	assert ( seg->get_dim() == 1 );
	assert ( o_seg->get_dim() == 1 );

	// unlink 'seg' from meshes above 'msh->cell_enclosed' (of all dimensions)
	Cell::Core * mce = msh->cell_enclosed;
	assert ( mce );
	Cell::Positive::NotVertex * pmce = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( mce );
	Cell A = seg->base_attr.reverse ( tag::surely_exists );
	Cell::Positive::Vertex * Ap = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	unlink_face_from_msh ( Ap, seg, msh, 0, 1 );
	unlink_face_from_higher ( Ap, pmce, 1, 0 );  // we switch the two counters
	Cell B = seg->tip_attr;
	Cell::Positive::Vertex * Bp = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	unlink_face_from_msh ( Bp, seg, msh, 1, 0 );
	unlink_face_from_higher ( Bp, pmce, 0, 1 );  // we switch the two counters

	unlink_face_from_higher ( seg, pmce, 0, 1 );  // we switch the two counters

	remove_link_same_dim_rev ( o_seg, seg, msh );
	
} // end of break_deep_connections_1d_rev with tag::mesh_is_not_bdry


inline void break_deep_connections_1d_rev  // hidden in anonymous namespace
(	Cell::Core * o_seg, Cell::Positive::Segment * seg, Mesh::Core * msh,
  const tag::MeshIsBdry &, const tag::DoNotBother &                   )

// break far connections when removing a negative cell
// see paragraph 11.9 in the manual
	
// tag::do_not_bother is only relevant for Mesh::Connected::***Dim

{	assert ( seg != o_seg );
	assert ( seg );  assert ( o_seg );
	assert ( msh->get_dim_plus_one() == 2 );
	assert ( seg->get_dim() == 1 );
	assert ( o_seg->get_dim() == 1 );

	// unlink 'seg' from meshes above 'msh->cell_enclosed' (of all dimensions)
	Cell::Core * mce = msh->cell_enclosed;
	assert ( mce );
	Cell::Positive::NotVertex * pmce = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( mce );
	Cell A = seg->base_attr.reverse ( tag::surely_exists );
	Cell::Positive::Vertex * Ap = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( A.core );
	unlink_face_from_msh ( Ap, seg, msh, 0, 1, tag::do_not_bother );
	unlink_face_from_higher ( Ap, pmce, 1, 0, tag::do_not_bother );  // we switch the two counters
	Cell B = seg->tip_attr;
	Cell::Positive::Vertex * Bp = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B.core );
	unlink_face_from_msh ( Bp, seg, msh, 1, 0, tag::do_not_bother );
	unlink_face_from_higher ( Bp, pmce, 0, 1, tag::do_not_bother );  // we switch the two counters

	unlink_face_from_higher ( seg, pmce, 0, 1, tag::do_not_bother );  // we switch the two counters

	remove_link_same_dim_rev ( o_seg, seg, msh, tag::do_not_bother );
	
} // end of break_deep_connections_1d_rev with tag::mesh_is_not_bdry, tag::do_not_bother


inline void break_deep_connections_hd  // hidden in anonymous namespace
(	Cell::Positive::HighDim * const cll, Mesh::Core * const msh, const tag::MeshIsNotBdry & )

// break far connections when removing a positive cell
// see paragraph 11.9 in the manual
	
{	assert ( cll );
	assert ( msh->get_dim_plus_one() > 1 );
	// break_deep_connections_0d deals with the case of a vertex cll
	// break_deep_connections_1d deals with the case of a segment cll
	size_t cll_dim = cll->get_dim();
	size_t cll_dim_m1 = tag::Util::assert_diff ( cll_dim, 1 );
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );
	assert ( cll );
	// we could assert  msh->cell_enclosed == nullptr  (i.e. mesh is not a boundary)
	// however, we sometimes call this method on the boundary of a
	// cell being destroyed (pretending it is not a boundary)

	Mesh cll_bdry = cll->boundary();
	for ( size_t d = 0; d < cll_dim_m1 ; d++ )
	{	Mesh::Iterator itt = cll_bdry.iterator  // as they are : positive
			( tag::over_cells, tag::of_dim, d, tag::as_they_are );
		// talvez implementar um iterador especial que devolva cp e cn
		for ( itt.reset(); itt.in_range(); itt++ )
		{	Cell::Core * fface = (*itt).core;  // add link from face to 'msh'
			short int cp, cn;
			compute_cp_cn ( cp, cn, fface, cll_bdry.core );
			unlink_face_from_msh ( fface, cll, msh, cp, cn );   }     }

	Mesh::Iterator it = cll_bdry.iterator
		( tag::over_cells, tag::of_max_dim, tag::force_positive );
	// talvez implementar um iterador especial que devolva cp e cn
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face = (*it).core;  // add link from face to 'msh'
		short int cp, cn;
		face->compute_sign ( cp, cn, cll_bdry.core );
		assert ( ( ( cp == 1 ) and ( cn == 0 ) ) or ( ( cp == 0 ) and ( cn == 1 ) ) );
	  unlink_face_from_msh ( face, cll, msh, cp, cn );                                }

	remove_link_same_dim ( cll, msh );

} // end of break_deep_connections_hd with tag::mesh_is_not_bdry


inline void break_deep_connections_hd  // hidden in anonymous namespace
(	Cell::Positive::HighDim * const cll, Mesh::Core * const msh,
  const tag::MeshIsNotBdry &, const tag::DoNotBother &        )

// break far connections when removing a positive cell
// see paragraph 11.9 in the manual
	
// tag::do_not_bother is only relevant for Mesh::Connected::***Dim

{	assert ( cll );
	assert ( msh->get_dim_plus_one() > 1 );
	// break_deep_connections_0d deals with the case of a vertex cll
	// break_deep_connections_1d deals with the case of a segment cll
	size_t cll_dim = cll->get_dim();
	size_t cll_dim_m1 = tag::Util::assert_diff ( cll_dim, 1 );
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );
	assert ( cll );
	// we could assert  msh->cell_enclosed == nullptr  (i.e. mesh is not a boundary)
	// however, we sometimes call this method on the boundary of a
	// cell being destroyed (pretending it is not a boundary)

	Mesh cll_bdry = cll->boundary();
	for ( size_t d = 0; d < cll_dim_m1 ; d++ )
	{	Mesh::Iterator itt = cll_bdry.iterator  // as they are : positive
			( tag::over_cells, tag::of_dim, d, tag::as_they_are );
		// talvez implementar um iterador especial que devolva cp e cn
		for ( itt.reset(); itt.in_range(); itt++ )
		{	Cell::Core * fface = (*itt).core;  // add link from face to 'msh'
			short int cp, cn;
			compute_cp_cn ( cp, cn, fface, cll_bdry.core );
			unlink_face_from_msh ( fface, cll, msh, cp, cn, tag::do_not_bother );  }  }

	Mesh::Iterator it = cll_bdry.iterator
		( tag::over_cells, tag::of_max_dim, tag::force_positive );
	// talvez implementar um iterador especial que devolva cp e cn
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face = (*it).core;  // add link from face to 'msh'
		short int cp, cn;
		face->compute_sign ( cp, cn, cll_bdry.core );
		assert ( ( ( cp == 1 ) and ( cn == 0 ) ) or ( ( cp == 0 ) and ( cn == 1 ) ) );
	  unlink_face_from_msh ( face, cll, msh, cp, cn, tag::do_not_bother );           }

	remove_link_same_dim ( cll, msh, tag::do_not_bother );

} // end of break_deep_connections_hd with tag::mesh_is_not_bdry, tag::do_not_bother


inline void break_deep_connections_hd  // hidden in anonymous namespace
(	Cell::Positive::HighDim * const cll, Mesh::Core * const msh, const tag::MeshIsBdry & )

// break far connections when removing a positive cell
// see paragraph 11.9 in the manual
	
{	assert ( cll );
	assert ( msh->get_dim_plus_one() > 2 );
	// break_deep_connections_0d deals with the case of a vertex cll
	// break_deep_connections_1d deals with the case of a segment cll
	size_t cll_dim = cll->get_dim();
	size_t cll_dim_m1 = tag::Util::assert_diff ( cll_dim, 1 );
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );

	// unlink 'cll' from all meshes above 'msh->cell_enclosed' (of all dimensions)
	Cell::Core * mce = msh->cell_enclosed;
	assert ( mce );
	Cell::Positive::NotVertex * pmce = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( mce );
	Mesh cll_bdry = cll->boundary();
	for ( size_t d = 0; d < cll_dim_m1 ; d++ )
	{	Mesh::Iterator itt = cll_bdry.iterator  // as they are : positive
			( tag::over_cells, tag::of_dim, d, tag::as_they_are );
		// talvez implementar um iterador especial que devolva cp e cn
		for ( itt.reset(); itt.in_range(); itt++ )
		{	Cell::Core * fface = (*itt).core;  // add link from face to 'msh' and
			// to all meshes above 'msh->cell_enclosed' (of all dimensions)
			short int cp, cn;
			compute_cp_cn ( cp, cn, fface, cll_bdry.core );
			unlink_face_from_msh ( fface, cll, msh, cp, cn );
			unlink_face_from_higher ( fface, pmce, cp, cn );            }  }

	Mesh::Iterator it = cll_bdry.iterator
		( tag::over_cells, tag::of_max_dim, tag::force_positive );
	// talvez implementar um iterador especial que devolva cp e cn
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face = (*it).core;  // add link from face to 'msh' and
		// to all meshes above 'msh->cell_enclosed' (of all dimensions)
		short int cp, cn;
		face->compute_sign ( cp, cn, cll_bdry.core );
		assert ( ( ( cp == 1 ) and ( cn == 0 ) ) or ( ( cp == 0 ) and ( cn == 1 ) ) );
	  unlink_face_from_msh ( face, cll, msh, cp, cn );
	  unlink_face_from_higher ( face, pmce, cp, cn );                                    }

	unlink_face_from_higher ( cll, pmce, 1, 0 );

	remove_link_same_dim ( cll, msh );
	
} // end of break_deep_connections_hd with tag::mesh_is_bdry


inline void break_deep_connections_hd  // hidden in anonymous namespace
(	Cell::Positive::HighDim * const cll, Mesh::Core * const msh,
  const tag::MeshIsBdry &, const tag::DoNotBother &           )

// break far connections when removing a positive cell
// see paragraph 11.9 in the manual
	
// tag::do_not_bother is only relevant for Mesh::Connected::***Dim

{	assert ( cll );
	assert ( msh->get_dim_plus_one() > 2 );
	// break_deep_connections_0d deals with the case of a vertex cll
	// break_deep_connections_1d deals with the case of a segment cll
	size_t cll_dim = cll->get_dim();
	size_t cll_dim_m1 = tag::Util::assert_diff ( cll_dim, 1 );
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );

	// unlink 'cll' from all meshes above 'msh->cell_enclosed' (of all dimensions)
	Cell::Core * mce = msh->cell_enclosed;
	assert ( mce );
	Cell::Positive::NotVertex * pmce = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( mce );
	Mesh cll_bdry = cll->boundary();
	for ( size_t d = 0; d < cll_dim_m1 ; d++ )
	{	Mesh::Iterator itt = cll_bdry.iterator  // as they are : positive
			( tag::over_cells, tag::of_dim, d, tag::as_they_are );
		// talvez implementar um iterador especial que devolva cp e cn
		for ( itt.reset(); itt.in_range(); itt++ )
		{	Cell::Core * fface = (*itt).core;  // add link from face to 'msh' and
			// to all meshes above 'msh->cell_enclosed' (of all dimensions)
			short int cp, cn;
			compute_cp_cn ( cp, cn, fface, cll_bdry.core );
			unlink_face_from_msh ( fface, cll, msh, cp, cn, tag::do_not_bother );
			unlink_face_from_higher ( fface, pmce, cp, cn, tag::do_not_bother );  }  }

	Mesh::Iterator it = cll_bdry.iterator
		( tag::over_cells, tag::of_max_dim, tag::force_positive );
	// talvez implementar um iterador especial que devolva cp e cn
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face = (*it).core;  // add link from face to 'msh' and
		// to all meshes above 'msh->cell_enclosed' (of all dimensions)
		short int cp, cn;
		face->compute_sign ( cp, cn, cll_bdry.core );
		assert ( ( ( cp == 1 ) and ( cn == 0 ) ) or ( ( cp == 0 ) and ( cn == 1 ) ) );
	  unlink_face_from_msh ( face, cll, msh, cp, cn, tag::do_not_bother );
	  unlink_face_from_higher ( face, pmce, cp, cn, tag::do_not_bother );            }

	unlink_face_from_higher ( cll, pmce, 1, 0, tag::do_not_bother );

	remove_link_same_dim ( cll, msh, tag::do_not_bother );
	
} // end of break_deep_connections_hd with tag::mesh_is_bdry, tag::do_not_bother


inline void break_deep_connections_hd_rev  // hidden in anonymous namespace
(	Cell::Core * const o_cll, Cell::Positive::HighDim * const cll,
  Mesh::Core * const msh, const tag::MeshIsNotBdry &               )

// break far connections when removing a negative cell
// see paragraph 11.9 in the manual
	
{	assert ( cll != o_cll );
	assert ( cll );  assert ( o_cll );  assert ( msh );
	assert ( msh->get_dim_plus_one() > 2 );
	// break_deep_connections_0d deals with the case of a vertex cll
	// break_deep_connections_1d deals with the case of a segment cll
	size_t cll_dim = cll->get_dim();
	size_t cll_dim_m1 = tag::Util::assert_diff ( cll_dim, 1 );
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );
	// we could assert  msh->cell_enclosed == nullptr  (i.e. mesh is not a boundary)
	// however, we sometimes call this method on the boundary of a
	// cell being destroyed (pretending it is not a boundary)

	Mesh cll_bdry = cll->boundary();
	for ( size_t d = 0; d < cll_dim_m1 ; d++ )
	{	Mesh::Iterator itt = cll_bdry.iterator  // as they are : positive
			( tag::over_cells, tag::of_dim, d, tag::as_they_are );
		// talvez implementar um iterador especial que devolva cp e cn
		for ( itt.reset(); itt.in_range(); itt++ )
		{	Cell::Core * fface = (*itt).core;  // add link from face to 'msh'
			short int cp, cn;
			compute_cp_cn ( cp, cn, fface, cll_bdry.core );
			unlink_face_from_msh ( fface, cll, msh, cp, cn );   }       }

	Mesh::Iterator it = cll_bdry.iterator
		( tag::over_cells, tag::of_max_dim, tag::force_positive );
	// talvez implementar um iterador especial que devolva cp e cn
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face = (*it).core;  // add link from face to 'msh'
		short int cp, cn;
		face->compute_sign ( cp, cn, cll_bdry.core );
		assert ( ( ( cp == 1 ) and ( cn == 0 ) ) or ( ( cp == 0 ) and ( cn == 1 ) ) );
	  unlink_face_from_msh ( face, cll, msh, cp, cn );                                   }

	remove_link_same_dim_rev ( o_cll, cll, msh );
	
} // end of break_deep_connections_hd_rev with tag::mesh_is_not_bdry


inline void break_deep_connections_hd_rev  // hidden in anonymous namespace
(	Cell::Core * const o_cll, Cell::Positive::HighDim * const cll, Mesh::Core * const msh,
  const tag::MeshIsNotBdry &, const tag::DoNotBother &                                  )

// break far connections when removing a negative cell
// see paragraph 11.9 in the manual
	
// tag::do_not_bother is only relevant for Mesh::Connected::***Dim

{	assert ( cll != o_cll );
	assert ( cll );  assert ( o_cll );  assert ( msh );
	assert ( msh->get_dim_plus_one() > 2 );
	// break_deep_connections_0d deals with the case of a vertex cll
	// break_deep_connections_1d deals with the case of a segment cll
	size_t cll_dim = cll->get_dim();
	size_t cll_dim_m1 = tag::Util::assert_diff ( cll_dim, 1 );
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );
	// we could assert  msh->cell_enclosed == nullptr  (i.e. mesh is not a boundary)
	// however, we sometimes call this method on the boundary of a
	// cell being destroyed (pretending it is not a boundary)

	Mesh cll_bdry = cll->boundary();
	for ( size_t d = 0; d < cll_dim_m1 ; d++ )
	{	Mesh::Iterator itt = cll_bdry.iterator  // as they are : positive
			( tag::over_cells, tag::of_dim, d, tag::as_they_are );
		// talvez implementar um iterador especial que devolva cp e cn
		for ( itt.reset(); itt.in_range(); itt++ )
		{	Cell::Core * fface = (*itt).core;  // add link from face to 'msh'
			short int cp, cn;
			compute_cp_cn ( cp, cn, fface, cll_bdry.core );
			unlink_face_from_msh ( fface, cll, msh, cp, cn, tag::do_not_bother );  }  }

	Mesh::Iterator it = cll_bdry.iterator
		( tag::over_cells, tag::of_max_dim, tag::force_positive );
	// talvez implementar um iterador especial que devolva cp e cn
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face = (*it).core;  // add link from face to 'msh'
		short int cp, cn;
		face->compute_sign ( cp, cn, cll_bdry.core );
		assert ( ( ( cp == 1 ) and ( cn == 0 ) ) or ( ( cp == 0 ) and ( cn == 1 ) ) );
	  unlink_face_from_msh ( face, cll, msh, cp, cn, tag::do_not_bother );          }

	remove_link_same_dim_rev ( o_cll, cll, msh, tag::do_not_bother );
	
} // end of break_deep_connections_hd_rev with tag::mesh_is_not_bdry, tag::do_not_bother


inline void break_deep_connections_hd_rev  // hidden in anonymous namespace
(	Cell::Core * o_cll, Cell::Positive::HighDim * cll,
  Mesh::Core * msh, const tag::MeshIsBdry &             )

// break far connections when removing a negative cell
// see paragraph 11.9 in the manual
	
{	assert ( cll );  assert ( o_cll );  assert ( msh );
	assert ( cll != o_cll );
	assert ( msh->get_dim_plus_one() > 1 );
	// break_deep_connections_0d deals with the case of a vertex cll
	// break_deep_connections_1d deals with the case of a segment cll
	size_t cll_dim = cll->get_dim();
	size_t cll_dim_m1 = tag::Util::assert_diff ( cll_dim, 1 );
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );

	// unlink 'cll' from all meshes above 'msh->cell_enclosed' (of all dimensions)
	Cell::Core * mce = msh->cell_enclosed;
	assert ( mce );
	Cell::Positive::NotVertex * pmce = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( mce );
	Mesh cll_bdry = cll->boundary();
	for ( size_t d = 0; d < cll_dim_m1 ; d++ )
	{	Mesh::Iterator itt = cll_bdry.iterator  // as they are : positive
			( tag::over_cells, tag::of_dim, d, tag::as_they_are );
		// talvez implementar um iterador especial que devolva cp e cn
		for ( itt.reset(); itt.in_range(); itt++ )
		{	Cell::Core * fface = (*itt).core;  // add link from face to 'msh' and
			// to all meshes above 'msh->cell_enclosed' (of all dimensions)
			short int cp, cn;
			compute_cp_cn ( cp, cn, fface, cll_bdry.core );
			unlink_face_from_msh ( fface, cll, msh, cp, cn );
			// we switch the two counters
			unlink_face_from_higher ( fface, pmce, cn, cp );    }  }

	Mesh::Iterator it = cll_bdry.iterator
		( tag::over_cells, tag::of_max_dim, tag::force_positive );
	// talvez implementar um iterador especial que devolva cp e cn
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face = (*it).core;  // add link from face to 'msh' and
		// to all meshes above 'msh->cell_enclosed' (of all dimensions)
		short int cp, cn;
		face->compute_sign ( cp, cn, cll_bdry.core );
		assert ( ( ( cp == 1 ) and ( cn == 0 ) ) or ( ( cp == 0 ) and ( cn == 1 ) ) );
	  unlink_face_from_msh ( face, cll, msh, cp, cn );
		// we switch the two counters
	  unlink_face_from_higher ( face, pmce, cn, cp );                                }

	unlink_face_from_higher ( cll, pmce, 0, 1 );  // we switch the two counters

	remove_link_same_dim_rev ( o_cll, cll, msh );
	
} // end of break_deep_connections_hd_rev with tag::mesh_is_bdry

	
inline void break_deep_connections_hd_rev  // hidden in anonymous namespace
(	Cell::Core * o_cll, Cell::Positive::HighDim * cll, Mesh::Core * msh,
  const tag::MeshIsBdry &, const tag::DoNotBother &                   )

// break far connections when removing a negative cell
// see paragraph 11.9 in the manual
	
// tag::do_not_bother is only relevant for Mesh::Connected::***Dim

{	assert ( cll );  assert ( o_cll );  assert ( msh );
	assert ( cll != o_cll );
	assert ( msh->get_dim_plus_one() > 1 );
	// break_deep_connections_0d deals with the case of a vertex cll
	// break_deep_connections_1d deals with the case of a segment cll
	size_t cll_dim = cll->get_dim();
	size_t cll_dim_m1 = tag::Util::assert_diff ( cll_dim, 1 );
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );

	// unlink 'cll' from all meshes above 'msh->cell_enclosed' (of all dimensions)
	Cell::Core * mce = msh->cell_enclosed;
	assert ( mce );
	Cell::Positive::NotVertex * pmce = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::NotVertex* > ( mce );
	Mesh cll_bdry = cll->boundary();
	for ( size_t d = 0; d < cll_dim_m1 ; d++ )
	{	Mesh::Iterator itt = cll_bdry.iterator  // as they are : positive
			( tag::over_cells, tag::of_dim, d, tag::as_they_are );
		// talvez implementar um iterador especial que devolva cp e cn
		for ( itt.reset(); itt.in_range(); itt++ )
		{	Cell::Core * fface = (*itt).core;  // add link from face to 'msh' and
			// to all meshes above 'msh->cell_enclosed' (of all dimensions)
			short int cp, cn;
			compute_cp_cn ( cp, cn, fface, cll_bdry.core );
			unlink_face_from_msh ( fface, cll, msh, cp, cn, tag::do_not_bother );
			// we switch the two counters
			unlink_face_from_higher ( fface, pmce, cn, cp, tag::do_not_bother );    }  }

	Mesh::Iterator it = cll_bdry.iterator
		( tag::over_cells, tag::of_max_dim, tag::force_positive );
	// talvez implementar um iterador especial que devolva cp e cn
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face = (*it).core;  // add link from face to 'msh' and
		// to all meshes above 'msh->cell_enclosed' (of all dimensions)
		short int cp, cn;
		face->compute_sign ( cp, cn, cll_bdry.core );
		assert ( ( ( cp == 1 ) and ( cn == 0 ) ) or ( ( cp == 0 ) and ( cn == 1 ) ) );
	  unlink_face_from_msh ( face, cll, msh, cp, cn, tag::do_not_bother );
		// we switch the two counters
	  unlink_face_from_higher ( face, pmce, cn, cp, tag::do_not_bother );            }

	unlink_face_from_higher ( cll, pmce, 0, 1, tag::do_not_bother );  // we switch the two counters

	remove_link_same_dim_rev ( o_cll, cll, msh, tag::do_not_bother );
	
} // end of break_deep_connections_hd_rev with tag::mesh_is_bdry, tag::do_not_bother

	
inline void add_cell_behind_below_pos_seg  // hidden in anonymous namespace
( Cell::Positive::Segment * const seg, Mesh::NotZeroDim * const that )

// just a block of code called from four versions of Mesh::NotZeroDim::add_pos_seg
	
///////////////////////////////////////////////////////////////////////////////////////////
	// inspired in item 24 of the book : Scott Meyers, Effective STL                       //
{	typedef std::map < Mesh::Core *, Cell > maptype;                                       //
	{  // just a block of code for hiding cmd and lb                                       //
	maptype & cmd = seg->base_attr.core->cell_behind_within;                               //
	maptype::iterator lb = cmd.lower_bound(that);                                          //
	assert ( ( lb == cmd.end() ) or ( cmd.key_comp()(that,lb->first) ) );                  //
	cmd.emplace_hint ( lb, std::piecewise_construct,                                       //
	      std::forward_as_tuple(that),                                                     //
	      std::forward_as_tuple(Cell(tag::whose_core_is,seg,                               //
	                                 tag::previously_existing,tag::surely_not_null)) );    //
	} {  // just a block of code for hiding cmd and lb                                     //
	maptype & cmd = seg->tip_attr.core->cell_behind_within;                                //
	maptype::iterator lb = cmd.lower_bound(that);                                          //
	assert ( ( lb == cmd.end() ) or ( cmd.key_comp()(that,lb->first) ) );                  //
	cmd.emplace_hint ( lb, std::piecewise_construct,                                       //
	      std::forward_as_tuple(that),                                                     //
	      std::forward_as_tuple(Cell(tag::whose_core_is,seg,                               //
	                                 tag::previously_existing,tag::surely_not_null)) );  }  }
/////////  code below is conceptually equivalent to the above  ////////////////////////////
//seg->base_attr.core->cell_behind_within[that] =                   //
//	Cell ( tag::whose_core_is, seg, tag::previously_existing );     //
//seg->tip_attr.core->cell_behind_within[that] =                    //
//	Cell ( tag::whose_core_is, seg, tag::previously_existing );     //
//////////////////////////////////////////////////////////////////////


inline void add_cell_behind_below_neg_seg  // hidden in anonymous namespace
( Cell::Negative::Segment * const seg, Cell::Positive::Segment * const pos_seg,
  Mesh::NotZeroDim * const that                                                )

// just a block of code called from four versions of Mesh::NotZeroDim::add_pos_seg
	
///////////////////////////////////////////////////////////////////////////////////////////
	// inspired in item 24 of the book : Scott Meyers, Effective STL                       //
{	typedef std::map < Mesh::Core *, Cell > maptype;                                       //
	{  // just a block of code for hiding cmd and lb                                       //
	maptype & cmd = pos_seg->base_attr.core->reverse_attr.core->cell_behind_within;        //
	maptype::iterator lb = cmd.lower_bound(that);                                          //
	assert ( ( lb == cmd.end() ) or ( cmd.key_comp()(that,lb->first) ) );                  //
	cmd.emplace_hint ( lb, std::piecewise_construct,                                       //
	      std::forward_as_tuple(that),                                                     //
	      std::forward_as_tuple(Cell(tag::whose_core_is,seg,                               //
	                                 tag::previously_existing,tag::surely_not_null)) );    //
	} {  // just a block of code for hiding cmd and lb                                     //
	maptype & cmd = pos_seg->tip_attr.core->reverse_attr.core->cell_behind_within;         //
	maptype::iterator lb = cmd.lower_bound(that);                                          //
	assert ( ( lb == cmd.end() ) or ( cmd.key_comp()(that,lb->first) ) );                  //
	cmd.emplace_hint ( lb, std::piecewise_construct,                                       //
	      std::forward_as_tuple(that),                                                     //
	      std::forward_as_tuple(Cell(tag::whose_core_is,seg,                               //
	                                 tag::previously_existing,tag::surely_not_null)) );  }  }
/////////  code below is conceptually equivalent to the above  ////////////////////////////
//pos_seg->base_attr.core->reverse_attr.core->cell_behind_within[that] =   //
//	Cell ( tag::whose_core_is, seg, tag::previously_existing );            //
//pos_seg->tip_attr.core->reverse_attr.core->cell_behind_within[that] =    //
//	Cell ( tag::whose_core_is, seg, tag::previously_existing );            //
/////////////////////////////////////////////////////////////////////////////


inline void add_cell_behind_below_pos_hd  // hidden in anonymous namespace
( Cell::Positive::HighDim * const cll, Mesh::NotZeroDim * const that )

// just a block of code called from four versions of Mesh::NotZeroDim::add_pos_hd_cell
	
{	Mesh bdry = cll->boundary_attr;
	assert ( bdry.core->get_dim_plus_one() + 1 == that->get_dim_plus_one() );
	Mesh::Iterator it = bdry.iterator ( tag::over_cells_of_max_dim, tag::as_they_are );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face_p = (*it).core;
/////////////////////////////////////////////////////////////////////////////////////////////
		// inspired in item 24 of the book : Scott Meyers, Effective STL                       //
		typedef std::map < Mesh::Core *, Cell > maptype;                                       //
		maptype & cmd = face_p->cell_behind_within;                                            //
		maptype::iterator lb = cmd.lower_bound(that);                                          //
		assert ( ( lb == cmd.end() ) or ( cmd.key_comp()(that,lb->first) ) );                  //
		cmd.emplace_hint ( lb, std::piecewise_construct,                                       //
		      std::forward_as_tuple(that),                                                     //
		      std::forward_as_tuple(Cell(tag::whose_core_is,cll,                               //
	                                   tag::previously_existing,tag::surely_not_null)) );  }  }
/////////  code below is conceptually equivalent to the above  //////////////////////////////
//	face_p->cell_behind_within[that] =                              //
//		Cell ( tag::whose_core_is, cll, tag::previously_existing );   //
//////////////////////////////////////////////////////////////////////


inline void add_cell_behind_below_neg_hd  // hidden in anonymous namespace
( Cell::Negative::HighDim * const cll, Cell::Positive::HighDim * const pos_cll,
  Mesh::NotZeroDim * const that                                                )

// just a block of code called from four versions of Mesh::NotZeroDim::add_neg_hd_cell
	
{	Mesh bdry = pos_cll->boundary_attr;
	assert ( bdry.core->get_dim_plus_one() + 1 == that->get_dim_plus_one() );
	Mesh::Iterator it = bdry.iterator ( tag::over_cells_of_max_dim, tag::as_they_are );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face_p = (*it).core;
		Cell::Core * rev_face = face_p->reverse_attr.core;
		assert ( rev_face );
/////////////////////////////////////////////////////////////////////////////////////////////
		// inspired in item 24 of the book : Scott Meyers, Effective STL                       //
		typedef std::map < Mesh::Core *, Cell > maptype;                                       //
		maptype & cmd = rev_face->cell_behind_within;                                          //
		maptype::iterator lb = cmd.lower_bound(that);                                          //
		assert ( ( lb == cmd.end() ) or ( cmd.key_comp()(that,lb->first) ) );                  //
		cmd.emplace_hint ( lb, std::piecewise_construct,                                       //
		      std::forward_as_tuple(that),                                                     //
		      std::forward_as_tuple(Cell(tag::whose_core_is,cll,                               //
	                                   tag::previously_existing,tag::surely_not_null)) );  }  }
/////////  code below is conceptually equivalent to the above  //////////////////////////////
//	rev_face->cell_behind_within[that] =                            //
//		Cell ( tag::whose_core_is, cll, tag::previously_existing );   //
//////////////////////////////////////////////////////////////////////

}  // anonymous namespace


void Cell::Positive::Vertex::add_to_seg ( Cell::Positive::Segment * seg )
// virtual from Cell::Core	
{	assert ( seg );
	assert ( not seg->tip_attr .exists() );
	seg->tip_attr = Cell ( tag::whose_core_is, this,
	                       tag::previously_existing, tag::surely_not_null );
	make_deep_connections_0d ( this, seg );                                   }


void Cell::Positive::Vertex::add_to_seg  // virtual from Cell::Core
( Cell::Positive::Segment * seg, const tag::DoNotBother & )
{	assert ( seg );
	assert ( not seg->tip_attr .exists() );
	seg->tip_attr = Cell ( tag::whose_core_is, this,
	                       tag::previously_existing, tag::surely_not_null );
	make_deep_connections_0d ( this, seg, tag::do_not_bother );               }


void Cell::Positive::Vertex::remove_from_seg ( Cell::Positive::Segment * seg )
// virtual from Cell::Core
{	assert ( seg );
	break_deep_connections_0d ( this, seg );
	assert ( seg->tip_attr .core == this );
	seg->tip_attr = Cell ( tag::non_existent );  }


void Cell::Positive::Vertex::remove_from_seg  // virtual from Cell::Core
( Cell::Positive::Segment * seg, const tag::DoNotBother & )
{	assert ( seg );
	break_deep_connections_0d ( this, seg, tag::do_not_bother );
	assert ( seg->tip_attr .core == this );
	seg->tip_attr = Cell ( tag::non_existent );                  }


void Cell::Negative::Vertex::add_to_seg ( Cell::Positive::Segment * seg )
// virtual from Cell::Core
{ assert ( seg );
	assert ( this->reverse_attr .exists() );
	Cell::Positive::Vertex * pos_ver = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( this->reverse_attr.core );
	// assert that 'this' vertex does not belong yet to the boundary of 'seg'
	assert ( not seg->base_attr .exists() );
	seg->base_attr = Cell ( tag::whose_core_is, this,
                          tag::previously_existing, tag::surely_not_null );
	make_deep_connections_0d_rev ( pos_ver, seg );                             }


void Cell::Negative::Vertex::add_to_seg  // virtual from Cell::Core
( Cell::Positive::Segment * seg, const tag::DoNotBother & )
{ assert ( seg );
	assert ( this->reverse_attr .exists() );
	Cell::Positive::Vertex * pos_ver = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( this->reverse_attr.core );
	// assert that 'this' vertex does not belong yet to the boundary of 'seg'
	assert ( not seg->base_attr .exists() );
	seg->base_attr = Cell ( tag::whose_core_is, this,
                          tag::previously_existing, tag::surely_not_null );
	make_deep_connections_0d_rev ( pos_ver, seg, tag::do_not_bother );        }


void Cell::Negative::Vertex::remove_from_seg ( Cell::Positive::Segment * seg )
// virtual from Cell::Core
{ assert ( seg );
	assert ( this->reverse_attr .exists() );
	Cell::Positive::Vertex * pos_ver = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( this->reverse_attr.core );
	break_deep_connections_0d_rev ( pos_ver, seg );
	assert ( seg->base_attr .core == this );
	seg->base_attr = Cell ( tag::non_existent );                             }


void Cell::Negative::Vertex::remove_from_seg  // virtual from Cell::Core
( Cell::Positive::Segment * seg, const tag::DoNotBother & )
{ assert ( seg );
	assert ( this->reverse_attr .exists() );
	Cell::Positive::Vertex * pos_ver = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( this->reverse_attr.core );
	break_deep_connections_0d_rev ( pos_ver, seg, tag::do_not_bother );
	assert ( seg->base_attr .core == this );
	seg->base_attr = Cell ( tag::non_existent );                             }


void Cell::Positive::NotVertex::add_to_seg ( Cell::Positive::Segment * seg )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
            << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot add a cell to a segment, use add_to_mesh." << std::endl;
	exit ( 1 );                                                                   }


void Cell::Positive::NotVertex::add_to_seg   // virtual from Cell::Core
( Cell::Positive::Segment * seg, const tag::DoNotBother & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
            << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot add a cell to a segment, use add_to_mesh." << std::endl;
	exit ( 1 );                                                                   }


void Cell::Negative::NotVertex::add_to_seg ( Cell::Positive::Segment * seg )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
            << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot add a cell to a segment, use add_to_mesh." << std::endl;
	exit ( 1 );                                                                   }


void Cell::Negative::NotVertex::add_to_seg   // virtual from Cell::Core
( Cell::Positive::Segment * seg, const tag::DoNotBother & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
            << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot add a cell to a segment, use add_to_mesh." << std::endl;
	exit ( 1 );                                                                   }


void Cell::Positive::NotVertex::remove_from_seg ( Cell::Positive::Segment * seg )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot remove a cell from a segment, use remove_from_mesh." << std::endl;
	exit ( 1 );                                                                             }


void Cell::Positive::NotVertex::remove_from_seg   // virtual from Cell::Core
( Cell::Positive::Segment * seg, const tag::DoNotBother & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot remove a cell from a segment, use remove_from_mesh." << std::endl;
	exit ( 1 );                                                                             }


void Cell::Negative::NotVertex::remove_from_seg ( Cell::Positive::Segment * seg )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot remove a cell from a segment, use remove_from_mesh." << std::endl;
	exit ( 1 );                                                                             }


void Cell::Negative::NotVertex::remove_from_seg   // virtual from Cell::Core
( Cell::Positive::Segment * seg, const tag::DoNotBother & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot remove a cell from a segment, use remove_from_mesh." << std::endl;
	exit ( 1 );                                                                             }

//-----------------------------------------------------------------------------//


void Cell::Positive::Vertex::add_to_mesh ( Mesh::Core * )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot add a vertex to a mesh, use add_to_seg." << std::endl;
	exit ( 1 );                                                                  }


void Cell::Positive::Vertex::add_to_mesh ( Mesh::Core *, const tag::DoNotBother & )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot add a vertex to a mesh, use add_to_seg." << std::endl;
	exit ( 1 );                                                                  }


void Cell::Positive::Vertex::remove_from_mesh ( Mesh::Core * )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
            << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot remove a vertex from a mesh, use remove_from_seg." << std::endl;
	exit ( 1 );                                                                            }


void Cell::Positive::Vertex::remove_from_mesh ( Mesh::Core *, const tag::DoNotBother & )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
            << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot remove a vertex from a mesh, use remove_from_seg." << std::endl;
	exit ( 1 );                                                                            }


void Cell::Negative::Vertex::add_to_mesh ( Mesh::Core * )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot add a vertex to a mesh, use add_to_seg." << std::endl;
	exit ( 1 );                                                                 }


void Cell::Negative::Vertex::add_to_mesh ( Mesh::Core *, const tag::DoNotBother & )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot add a vertex to a mesh, use add_to_seg." << std::endl;
	exit ( 1 );                                                                 }


void Cell::Negative::Vertex::remove_from_mesh ( Mesh::Core * )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot remove a vertex from a mesh, use remove_from_seg." << std::endl;
	exit ( 1 );                                                                            }


void Cell::Negative::Vertex::remove_from_mesh ( Mesh::Core *, const tag::DoNotBother & )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot remove a vertex from a mesh, use remove_from_seg." << std::endl;
	exit ( 1 );                                                                            }


void Cell::Positive::Segment::add_to_mesh ( Mesh::Core * msh )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh );  msh->add_pos_seg ( this, tag::mesh_is_not_bdry );  }


void Cell::Positive::Segment::add_to_mesh ( Mesh::Core * msh, const tag::DoNotBother & )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh );  msh->add_pos_seg ( this, tag::mesh_is_not_bdry, tag::do_not_bother );  }


void Cell::Positive::Segment::remove_from_mesh ( Mesh::Core * msh )
// virtual from Cell::Core

// remove 'this' cell from the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'cut_from_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh );  msh->remove_pos_seg ( this, tag::mesh_is_not_bdry );  }


void Cell::Positive::Segment::remove_from_mesh ( Mesh::Core * msh, const tag::DoNotBother & )
// virtual from Cell::Core

// remove 'this' cell from the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'cut_from_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh );
	msh->remove_pos_seg ( this, tag::mesh_is_not_bdry, tag::do_not_bother );  }


void Cell::Negative::Segment::add_to_mesh ( Mesh::Core * msh )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh );  msh->add_neg_seg ( this, tag::mesh_is_not_bdry );  }


void Cell::Negative::Segment::add_to_mesh ( Mesh::Core * msh, const tag::DoNotBother & )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh );  msh->add_neg_seg ( this, tag::mesh_is_not_bdry, tag::do_not_bother );  }


void Cell::Negative::Segment::remove_from_mesh ( Mesh::Core * msh )
// virtual from Cell::Core

// remove 'this' cell from the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'cut_from_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh );  msh->remove_neg_seg ( this, tag::mesh_is_not_bdry );  }


void Cell::Negative::Segment::remove_from_mesh ( Mesh::Core * msh, const tag::DoNotBother & )
// virtual from Cell::Core

// remove 'this' cell from the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'cut_from_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh );  msh->remove_neg_seg ( this, tag::mesh_is_not_bdry, tag::do_not_bother );  }


void Cell::Positive::HighDim::add_to_mesh ( Mesh::Core * msh )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh );  msh->add_pos_hd_cell ( this, tag::mesh_is_not_bdry );  }


void Cell::Positive::HighDim::add_to_mesh ( Mesh::Core * msh, const tag::DoNotBother & )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh );  msh->add_pos_hd_cell ( this, tag::mesh_is_not_bdry, tag::do_not_bother );  }


void Cell::Positive::HighDim::remove_from_mesh ( Mesh::Core * msh )
// virtual from Cell::Core

// remove 'this' cell from the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'cut_from_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh );  msh->remove_pos_hd_cell ( this, tag::mesh_is_not_bdry );  }


void Cell::Positive::HighDim::remove_from_mesh ( Mesh::Core * msh, const tag::DoNotBother & )
// virtual from Cell::Core

// remove 'this' cell from the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'cut_from_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh );  msh->remove_pos_hd_cell ( this, tag::mesh_is_not_bdry, tag::do_not_bother );  }


void Cell::Negative::HighDim::add_to_mesh ( Mesh::Core * msh )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh );  msh->add_neg_hd_cell ( this, tag::mesh_is_not_bdry );  }


void Cell::Negative::HighDim::add_to_mesh ( Mesh::Core * msh, const tag::DoNotBother & )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh );  msh->add_neg_hd_cell ( this, tag::mesh_is_not_bdry, tag::do_not_bother );  }


void Cell::Negative::HighDim::remove_from_mesh ( Mesh::Core * msh )
// virtual from Cell::Core

// remove 'this' cell from the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'cut_from_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh );  msh->remove_neg_hd_cell ( this, tag::mesh_is_not_bdry );  }


void Cell::Negative::HighDim::remove_from_mesh ( Mesh::Core * msh,const tag::DoNotBother & )
// virtual from Cell::Core

// remove 'this' cell from the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'cut_from_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh );  msh->remove_neg_hd_cell ( this, tag::mesh_is_not_bdry, tag::do_not_bother );  }


void Cell::Positive::Vertex::add_to_bdry ( Mesh::Core * )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Please use add_to_seg." << std::endl;
	exit ( 1 );                                               }


void Cell::Positive::Vertex::add_to_bdry ( Mesh::Core *, const tag::DoNotBother & )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Please use add_to_seg." << std::endl;
	exit ( 1 );                                               }


void Cell::Positive::Vertex::remove_from_bdry ( Mesh::Core * )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
            << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Please use remove_from_seg." << std::endl;
	exit ( 1 );                                               }


void Cell::Positive::Vertex::remove_from_bdry ( Mesh::Core *, const tag::DoNotBother & )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
            << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Please use remove_from_seg." << std::endl;
	exit ( 1 );                                               }


void Cell::Negative::Vertex::add_to_bdry ( Mesh::Core * )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Please use add_to_seg." << std::endl;
	exit ( 1 );                                               }


void Cell::Negative::Vertex::add_to_bdry ( Mesh::Core *, const tag::DoNotBother & )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Please use add_to_seg." << std::endl;
	exit ( 1 );                                               }


void Cell::Negative::Vertex::remove_from_bdry ( Mesh::Core * )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Please use remove_from_seg." << std::endl;
	exit ( 1 );                                               }


void Cell::Negative::Vertex::remove_from_bdry ( Mesh::Core *, const tag::DoNotBother & )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Please use remove_from_seg." << std::endl;
	exit ( 1 );                                               }


void Cell::Positive::Segment::add_to_bdry ( Mesh::Core * msh )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh' which is the boundary of some cell
// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->add_pos_seg ( this, tag::mesh_is_bdry );  }


void Cell::Positive::Segment::add_to_bdry ( Mesh::Core * msh, const tag::DoNotBother & )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh' which is the boundary of some cell
// there are many types of meshes, so we call a virtual method
// tag::do_not_bother is only relevant for Mesh:Connected::***Dim

{	assert ( msh ); msh->add_pos_seg ( this, tag::mesh_is_bdry, tag::do_not_bother );  }


void Cell::Positive::Segment::remove_from_bdry ( Mesh::Core * msh )
// virtual from Cell::Core

// remove 'this' cell from the mesh 'msh' which is the boundary of some cell
// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->remove_pos_seg ( this, tag::mesh_is_bdry );  }


void Cell::Positive::Segment::remove_from_bdry ( Mesh::Core * msh, const tag::DoNotBother & )
// virtual from Cell::Core

// remove 'this' cell from the mesh 'msh' which is the boundary of some cell
// there are many types of meshes, so we call a virtual method
// tag::do_not_bother is only relevant for Mesh:Connected::***Dim

{	assert ( msh ); msh->remove_pos_seg ( this, tag::mesh_is_bdry, tag::do_not_bother );  }


void Cell::Negative::Segment::add_to_bdry ( Mesh::Core * msh )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh' which is the boundary of some cell
// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->add_neg_seg ( this, tag::mesh_is_bdry );  }


void Cell::Negative::Segment::add_to_bdry ( Mesh::Core * msh, const tag::DoNotBother & )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh' which is the boundary of some cell
// there are many types of meshes, so we call a virtual method
// tag::do_not_bother is only relevant for Mesh:Connected::***Dim

{	assert ( msh ); msh->add_neg_seg ( this, tag::mesh_is_bdry, tag::do_not_bother );  }


void Cell::Negative::Segment::remove_from_bdry ( Mesh::Core * msh )
// virtual from Cell::Core

// remove 'this' cell from the mesh 'msh' which is the boundary of some cell
// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->remove_neg_seg ( this, tag::mesh_is_bdry );  }


void Cell::Negative::Segment::remove_from_bdry ( Mesh::Core * msh, const tag::DoNotBother & )
// virtual from Cell::Core

// remove 'this' cell from the mesh 'msh' which is the boundary of some cell
// there are many types of meshes, so we call a virtual method
// tag::do_not_bother is only relevant for Mesh:Connected::***Dim

{	assert ( msh ); msh->remove_neg_seg ( this, tag::mesh_is_bdry, tag::do_not_bother );  }


void Cell::Positive::HighDim::add_to_bdry ( Mesh::Core * msh )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh' which is the boundary of some cell
// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->add_pos_hd_cell ( this, tag::mesh_is_bdry );  }


void Cell::Positive::HighDim::add_to_bdry ( Mesh::Core * msh, const tag::DoNotBother & )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh' which is the boundary of some cell
// there are many types of meshes, so we call a virtual method
// tag::do_not_bother is only relevant for Mesh:Connected::***Dim

{	assert ( msh ); msh->add_pos_hd_cell ( this, tag::mesh_is_bdry, tag::do_not_bother );  }


void Cell::Positive::HighDim::remove_from_bdry ( Mesh::Core * msh )
// virtual from Cell::Core

// remove 'this' cell from the mesh 'msh' which is the boundary of some cell
// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->remove_pos_hd_cell ( this, tag::mesh_is_bdry );  }


void Cell::Positive::HighDim::remove_from_bdry ( Mesh::Core * msh, const tag::DoNotBother & )
// virtual from Cell::Core

// remove 'this' cell from the mesh 'msh' which is the boundary of some cell
// there are many types of meshes, so we call a virtual method
// tag::do_not_bother is only relevant for Mesh:Connected::***Dim

{	assert ( msh ); msh->remove_pos_hd_cell ( this, tag::mesh_is_bdry, tag::do_not_bother );  }


void Cell::Negative::HighDim::add_to_bdry ( Mesh::Core * msh )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh' which is the boundary of some cell
// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->add_neg_hd_cell ( this, tag::mesh_is_bdry );  }


void Cell::Negative::HighDim::add_to_bdry ( Mesh::Core * msh, const tag::DoNotBother & )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh' which is the boundary of some cell
// there are many types of meshes, so we call a virtual method
// tag::do_not_bother is only relevant for Mesh:Connected::***Dim

{	assert ( msh ); msh->add_neg_hd_cell ( this, tag::mesh_is_bdry, tag::do_not_bother );  }


void Cell::Negative::HighDim::remove_from_bdry ( Mesh::Core * msh )
// virtual from Cell::Core

// remove 'this' cell from the mesh 'msh' which is the boundary of some cell
// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->remove_neg_hd_cell ( this, tag::mesh_is_bdry );  }


void Cell::Negative::HighDim::remove_from_bdry ( Mesh::Core * msh, const tag::DoNotBother & )
// virtual from Cell::Core

// remove 'this' cell from the mesh 'msh' which is the boundary of some cell
// there are many types of meshes, so we call a virtual method
// tag::do_not_bother is only relevant for Mesh:Connected::***Dim

{	assert ( msh ); msh->remove_neg_hd_cell ( this, tag::mesh_is_bdry, tag::do_not_bother );  }

//-----------------------------------------------------------------------------//


std::list<Cell>::iterator Mesh::ZeroDim::add_to_my_cells
(	Cell::Core * const cll, const size_t d )
// virtual from Mesh::Core, here execution forbidden

// called from add_link_same_dim and add_link (both hidden in anonymous namespace above)
	
{	assert ( d == 0 );  assert ( false );
	return static_cast < std::list<Cell>::iterator > ( nullptr );     }


std::list<Cell>::iterator Mesh::ZeroDim::add_to_my_cells
(	Cell::Core * const cll, const size_t d, const tag::DoNotBother & )
// virtual from Mesh::Core, here execution forbidden

// called from add_link_same_dim and add_link (both hidden in anonymous namespace above)
	
{	assert ( d == 0 );  assert ( false );
	return static_cast < std::list<Cell>::iterator > ( nullptr );     }


std::list<Cell>::iterator Mesh::Connected::OneDim::add_to_my_cells
(	Cell::Core * const cll, const size_t d )
// virtual from Mesh::Core, here returns garbage
// updates nb_of_segs, first_ver, last_ver

// called from add_link_same_dim and add_link (both hidden in anonymous namespace above)
	
{	assert ( d == cll->get_dim() );
	assert ( d <= 1 );
	assert ( this->last_ver != this->first_ver .reverse() );
	if ( d == 1 )
	{	// we cannot add a segment to a closed loop :
		assert ( this->last_ver != this->first_ver .reverse() );
		assert ( ( cll->base() .reverse() == this->last_ver ) or
	           ( cll->tip() == this->first_ver .reverse() )    );
		this->nb_of_segs++;
		if ( cll->base() .reverse() == this->last_ver )
			this->last_ver = cll->tip();
		else  // cll->tip() == this->first_ver.reverse()
			this->first_ver = cll->base();                           }
	return static_cast < std::list<Cell>::iterator > ( nullptr );     }


std::list<Cell>::iterator Mesh::Connected::OneDim::add_to_my_cells
(	Cell::Core * const cll, const size_t d, const tag::DoNotBother & )
// virtual from Mesh::Core, here returns garbage
// do not bother with nb_of_segs, first_ver, last_ver

// called from add_link_same_dim and add_link (both hidden in anonymous namespace above)
	
{	assert ( d == cll->get_dim() );
	assert ( d <= 1 );
	return static_cast < std::list<Cell>::iterator > ( nullptr );  }


std::list<Cell>::iterator Mesh::Fuzzy::add_to_my_cells
(	Cell::Core * cll, const size_t d )
// virtual from Mesh::Core, later overriden by Mesh::STSI

// called from add_link_same_dim and add_link (both hidden in anonymous namespace above)

// add a cell to 'this->cells [d]' list, return iterator into that list
		
{	assert ( d == cll->get_dim() );
	assert ( d < this->get_dim_plus_one() );
	std::list < Cell > & mcd = this->cells [d];
	mcd .push_front ( Cell ( tag::whose_core_is, cll,
                           tag::previously_existing, tag::surely_not_null ) );
	return mcd .begin();                                                          }


std::list<Cell>::iterator Mesh::Fuzzy::add_to_my_cells
(	Cell::Core * cll, const size_t d, const tag::DoNotBother & )
// virtual from Mesh::Core
// tag::do_not_bother makes no difference here, the body is the same as above

// called from add_link_same_dim and add_link (both hidden in anonymous namespace above)

// add a cell to 'this->cells [d]' list, return iterator into that list
		
{	assert ( d == cll->get_dim() );
	assert ( d < this->get_dim_plus_one() );
	std::list <Cell> & mcd = this->cells [d];
	mcd .push_front ( Cell ( tag::whose_core_is, cll,
                           tag::previously_existing, tag::surely_not_null ) );
	return mcd .begin();                                                          }


std::list<Cell>::iterator Mesh::STSI::add_to_my_cells
(	Cell::Core * cll, const size_t d )
// virtual from Mesh::Core, defined by Mesh::Fuzzy, here overriden

// called from add_link_same_dim and add_link (both hidden in anonymous namespace above)

// add a cell to 'this->cells[d]' list, return iterator into that list
	
{	return Mesh::Fuzzy::add_to_my_cells ( cll, d );  }  // will change


void Mesh::ZeroDim::remove_from_my_cells
(	Cell::Core * cll, const size_t d, std::list<Cell>::iterator )
// virtual from Mesh::Core, here execution forbidden

// called from remove_link_same_dim and remove_link (both hidden in anonymous namespace above)

{	assert ( d == 0 );  assert ( false );  }


void Mesh::ZeroDim::remove_from_my_cells
(	Cell::Core * cll, const size_t d, std::list<Cell>::iterator, const tag::DoNotBother & )
// virtual from Mesh::Core, here execution forbidden

// called from remove_link_same_dim and remove_link (both hidden in anonymous namespace above)

{	assert ( d == 0 );  assert ( false );  }


void Mesh::Connected::OneDim::remove_from_my_cells
(	Cell::Core * cll, const size_t d, std::list<Cell>::iterator )
// virtual from Mesh::Core, here updates nb_of_segs, first_ver, last_ver

// called from remove_link_same_dim and remove_link (both hidden in anonymous namespace above)

{	assert ( d == cll->get_dim() );
	assert ( d <= 1 );
	if ( d == 1 )
	{	assert ( this->nb_of_segs >= 1 );
		// more often than not, this->nb_of_segs >= 2
		// thus, after eliminating the 'cll' segment, 'this' mesh wont become empty
		// however, destructors of meshes also call this method
		// in which case the mesh may become empty
		this->nb_of_segs--;
		// if we remove a segment from a closed loop, it will become open chain :
		if ( this->last_ver == this->first_ver .reverse() )
		{	this->last_ver = cll->base() .reverse();
			this->first_ver = cll->tip() .reverse();  }
		else  // open chain
		{	assert ( ( cll->base() == this->first_ver ) or
			         ( cll->tip() == this->last_ver )      );
			if ( cll->base() == this->first_ver )
				this->first_ver = cll->tip() .reverse();
			else  // cll->tip() == this->last_ver
				this->last_ver = cll->base() .reverse();            }  }  }


void Mesh::Connected::OneDim::remove_from_my_cells
(	Cell::Core * cll, const size_t d, std::list<Cell>::iterator, const tag::DoNotBother & )
// virtual from Mesh::Core, here does nothing
// do not bother with nb_of_segs, first_ver, last_ver

// called from remove_link_same_dim and remove_link (both hidden in anonymous namespace above)

{	assert ( d == cll->get_dim() );
	assert ( d <= 1 );               }


void Mesh::Fuzzy::remove_from_my_cells
(	Cell::Core * cll, const size_t d, std::list<Cell>::iterator it )
// virtual from Cell::Core, later overriden by Mesh::STSI

// called from remove_link_same_dim and remove_link (both hidden in anonymous namespace above)

// remove a cell from 'this->cells [d]' list using the provided iterator

{	assert ( d == cll->get_dim() );
	assert ( d < this->get_dim_plus_one() );
	assert ( it != this->cells [d] .end() );
	this->cells [d] .erase(it);              }


void Mesh::Fuzzy::remove_from_my_cells
(	Cell::Core * cll, const size_t d, std::list<Cell>::iterator it, const tag::DoNotBother & )
// virtual from Cell::Core
// tag::do_not_bother makes no difference here, the body is the same as above

// called from remove_link_same_dim and remove_link (both hidden in anonymous namespace above)

// remove a cell from 'this->cells[d]' list using the provided iterator

{	assert ( d == cll->get_dim() );
	assert ( d < this->get_dim_plus_one() );
	assert ( it != this->cells [d] .end() );
	this->cells [d] .erase(it);              }


void Mesh::STSI::remove_from_my_cells
(	Cell::Core * cll, const size_t d, std::list<Cell>::iterator it )
// virtual from Cell::Core, overriden here a second time

// remove a cell from 'this->cells [d]' list using the provided iterator

{	Mesh::Fuzzy::remove_from_my_cells ( cll, d, it );  }

//-----------------------------------------------------------------------------//


// see paragraph 11.9 in the manual

void Mesh::ZeroDim::add_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsNotBdry & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::add_pos_seg
( Cell::Positive::Segment *, const tag::MeshIsNotBdry &, const tag::DoNotBother & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::add_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsBdry & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::add_pos_seg
( Cell::Positive::Segment *, const tag::MeshIsBdry &, const tag::DoNotBother & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::remove_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsNotBdry & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::remove_pos_seg
( Cell::Positive::Segment *, const tag::MeshIsNotBdry &, const tag::DoNotBother & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::remove_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsBdry & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::remove_pos_seg
( Cell::Positive::Segment *, const tag::MeshIsBdry &, const tag::DoNotBother & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::add_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsNotBdry & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::add_neg_seg
( Cell::Negative::Segment *, const tag::MeshIsNotBdry &, const tag::DoNotBother & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::add_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsBdry & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::add_neg_seg
( Cell::Negative::Segment *, const tag::MeshIsBdry &, const tag::DoNotBother & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::remove_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsNotBdry & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::remove_neg_seg
( Cell::Negative::Segment *, const tag::MeshIsNotBdry &, const tag::DoNotBother & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::remove_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsBdry & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::remove_neg_seg
( Cell::Negative::Segment *, const tag::MeshIsBdry &, const tag::DoNotBother & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::add_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsNotBdry & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::add_pos_hd_cell
( Cell::Positive::HighDim *, const tag::MeshIsNotBdry &, const tag::DoNotBother & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::add_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsBdry & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::add_pos_hd_cell
( Cell::Positive::HighDim *, const tag::MeshIsBdry &, const tag::DoNotBother & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::remove_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsNotBdry & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::remove_pos_hd_cell
( Cell::Positive::HighDim *, const tag::MeshIsNotBdry &, const tag::DoNotBother & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::remove_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsBdry & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::remove_pos_hd_cell
( Cell::Positive::HighDim *, const tag::MeshIsBdry &, const tag::DoNotBother & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::add_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsNotBdry & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::add_neg_hd_cell
( Cell::Negative::HighDim *, const tag::MeshIsNotBdry &, const tag::DoNotBother & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::add_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsBdry & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::add_neg_hd_cell
( Cell::Negative::HighDim *, const tag::MeshIsBdry &, const tag::DoNotBother & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden
	
void Mesh::ZeroDim::remove_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsNotBdry & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden

void Mesh::ZeroDim::remove_neg_hd_cell
( Cell::Negative::HighDim *, const tag::MeshIsNotBdry &, const tag::DoNotBother & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden

void Mesh::ZeroDim::remove_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsBdry & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden

void Mesh::ZeroDim::remove_neg_hd_cell
( Cell::Negative::HighDim *, const tag::MeshIsBdry &, const tag::DoNotBother & )
{	assert ( false ); }  // virtual from Mesh::Core, here execution forbidden


void Mesh::NotZeroDim::add_pos_seg ( Cell::Positive::Segment * seg, const tag::MeshIsNotBdry & )
// virtual from Mesh::Core

{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->base_attr .exists() );
	assert ( seg->tip_attr .exists() );
	// assert that 'seg' does not belong yet to 'this' mesh
	assert ( seg->meshes_same_dim .find (this) == seg->meshes_same_dim .end() );

	make_deep_connections_1d ( seg, this, tag::mesh_is_not_bdry );

	add_cell_behind_below_pos_seg ( seg, this );                                 }


void Mesh::NotZeroDim::add_pos_seg
( Cell::Positive::Segment * seg, const tag::MeshIsNotBdry &, const tag::DoNotBother & )
// virtual from Mesh::Core

{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->base_attr .exists() );
	assert ( seg->tip_attr .exists() );
	// assert that 'seg' does not belong yet to 'this' mesh
	assert ( seg->meshes_same_dim .find (this) == seg->meshes_same_dim .end() );

	make_deep_connections_1d ( seg, this, tag::mesh_is_not_bdry, tag::do_not_bother );

	add_cell_behind_below_pos_seg ( seg, this );                                        }


void Mesh::NotZeroDim::add_pos_seg ( Cell::Positive::Segment * seg, const tag::MeshIsBdry & )
// virtual from Mesh::Core

{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->base_attr .exists() );
	assert ( seg->tip_attr .exists() );
	// assert that 'seg' does not belong yet to 'this' mesh
	assert ( seg->meshes_same_dim .find (this) == seg->meshes_same_dim .end() );

	make_deep_connections_1d ( seg, this, tag::mesh_is_bdry );

	add_cell_behind_below_pos_seg ( seg, this );                                  }


void Mesh::NotZeroDim::add_pos_seg
( Cell::Positive::Segment * seg, const tag::MeshIsBdry &, const tag::DoNotBother & )
// virtual from Mesh::Core

{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->base_attr .exists() );
	assert ( seg->tip_attr .exists() );
	// assert that 'seg' does not belong yet to 'this' mesh
	assert ( seg->meshes_same_dim .find (this) == seg->meshes_same_dim .end() );

	make_deep_connections_1d ( seg, this, tag::mesh_is_bdry, tag::do_not_bother );

	add_cell_behind_below_pos_seg ( seg, this );                                    }


void Mesh::NotZeroDim::remove_pos_seg  // virtual from Mesh::Core
( Cell::Positive::Segment * seg, const tag::MeshIsNotBdry & )

{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->base_attr .exists() );
	assert ( seg->tip_attr .exists() );
	// assert that 'seg' belongs to 'this' mesh
	assert ( seg->meshes_same_dim .find (this) != seg->meshes_same_dim .end() );

	assert ( seg->base_attr .core->cell_behind_within .find (this) !=
	         seg->base_attr .core->cell_behind_within .end()          );
	seg->base_attr.core->cell_behind_within .erase (this);

	assert ( seg->tip_attr .core->cell_behind_within .find (this) !=
	         seg->tip_attr .core->cell_behind_within .end()          );
	seg->tip_attr.core->cell_behind_within .erase (this);

  break_deep_connections_1d ( seg, this, tag::mesh_is_not_bdry );

}  // end of Mesh::NotZeroDim::remove_pos_seg with tag::mesh_is_not_bdry


void Mesh::NotZeroDim::remove_pos_seg  // virtual from Mesh::Core
( Cell::Positive::Segment * seg, const tag::MeshIsNotBdry &, const tag::DoNotBother & )

{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->base_attr .exists() );
	assert ( seg->tip_attr .exists() );
	// assert that 'seg' belongs to 'this' mesh
	assert ( seg->meshes_same_dim .find (this) != seg->meshes_same_dim .end() );

	assert ( seg->base_attr .core->cell_behind_within .find (this) !=
	         seg->base_attr .core->cell_behind_within .end()         );
	seg->base_attr .core->cell_behind_within .erase (this);

	assert ( seg->tip_attr .core->cell_behind_within .find (this) !=
	         seg->tip_attr .core->cell_behind_within .end()         );
	seg->tip_attr .core->cell_behind_within .erase (this);

  break_deep_connections_1d ( seg, this, tag::mesh_is_not_bdry, tag::do_not_bother );

}  // end of Mesh::NotZeroDim::remove_pos_seg with tag::mesh_is_not_bdry, tag::do_not_bother


void Mesh::NotZeroDim::remove_pos_seg ( Cell::Positive::Segment * seg, const tag::MeshIsBdry & )
// virtual from Mesh::Core

{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->base_attr .exists() );
	assert ( seg->tip_attr .exists() );
	// assert that 'seg' belongs to 'this' mesh
	assert ( seg->meshes_same_dim .find (this) != seg->meshes_same_dim .end() );

	assert ( seg->base_attr .core->cell_behind_within .find (this) !=
	         seg->base_attr .core->cell_behind_within .end()         );
	seg->base_attr .core->cell_behind_within .erase (this);
	assert ( seg->tip_attr .core->cell_behind_within .find (this) !=
	         seg->tip_attr .core->cell_behind_within .end()         );
	seg->tip_attr .core->cell_behind_within .erase (this);

	break_deep_connections_1d ( seg, this, tag::mesh_is_bdry );

}  // end of Mesh::NotZeroDim::remove_pos_seg with tag::mesh_is_bdry


void Mesh::NotZeroDim::remove_pos_seg  // virtual from Mesh::Core
( Cell::Positive::Segment * seg, const tag::MeshIsBdry &, const tag::DoNotBother & )

{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->base_attr .exists() );
	assert ( seg->tip_attr .exists() );
	// assert that 'seg' belongs to 'this' mesh
	assert ( seg->meshes_same_dim .find (this) != seg->meshes_same_dim .end() );

	assert ( seg->base_attr .core->cell_behind_within .find (this) !=
	         seg->base_attr .core->cell_behind_within .end()         );
	seg->base_attr .core->cell_behind_within.erase(this);
	assert ( seg->tip_attr .core->cell_behind_within .find (this) !=
	         seg->tip_attr .core->cell_behind_within .end()         );
	seg->tip_attr .core->cell_behind_within .erase (this);

	break_deep_connections_1d ( seg, this, tag::mesh_is_bdry, tag::do_not_bother );

}  // end of Mesh::NotZeroDim::remove_pos_seg with tag::mesh_is_bdry, tag::do_not_bother


void Mesh::NotZeroDim::add_neg_seg ( Cell::Negative::Segment * seg, const tag::MeshIsNotBdry & )
// virtual from Mesh::Core

{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->reverse_attr .exists() );
	Cell::Positive::Segment * pos_seg = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Segment* > ( seg->reverse_attr .core );
	// assert that 'pos_seg' does not belong yet to 'this' mesh
	assert ( pos_seg->meshes_same_dim .find (this) == pos_seg->meshes_same_dim .end() );

	assert ( pos_seg->base_attr .exists() );
	assert ( pos_seg->tip_attr .exists() );
	assert ( pos_seg->base_attr .core->reverse_attr .exists() );
	assert ( pos_seg->tip_attr .core->reverse_attr .exists() );

	make_deep_connections_1d_rev ( seg, pos_seg, this, tag::mesh_is_not_bdry );

	add_cell_behind_below_neg_seg ( seg, pos_seg, this );                                  }

	
void Mesh::NotZeroDim::add_neg_seg  // virtual from Mesh::Core
( Cell::Negative::Segment * seg, const tag::MeshIsNotBdry &, const tag::DoNotBother & )

{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->reverse_attr .exists() );
	Cell::Positive::Segment * pos_seg = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Segment* > ( seg->reverse_attr .core );
	// assert that 'pos_seg' does not belong yet to 'this' mesh
	assert ( pos_seg->meshes_same_dim .find(this) == pos_seg->meshes_same_dim .end() );

	assert ( pos_seg->base_attr .exists() );
	assert ( pos_seg->tip_attr .exists() );
	assert ( pos_seg->base_attr .core->reverse_attr .exists() );
	assert ( pos_seg->tip_attr .core->reverse_attr .exists() );

	make_deep_connections_1d_rev
		( seg, pos_seg, this, tag::mesh_is_not_bdry, tag::do_not_bother );

	add_cell_behind_below_neg_seg ( seg, pos_seg, this );                                 }

	
void Mesh::NotZeroDim::add_neg_seg ( Cell::Negative::Segment * seg, const tag::MeshIsBdry & )
// virtual from Mesh::Core

{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->reverse_attr .exists() );
	Cell::Positive::Segment * pos_seg = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Segment* > ( seg->reverse_attr .core );
	// assert that 'pos_seg' does not belong yet to 'this' mesh
	assert ( pos_seg->meshes_same_dim .find (this) == pos_seg->meshes_same_dim .end() );

	assert ( pos_seg->base_attr .exists() );
	assert ( pos_seg->tip_attr .exists() );
	assert ( pos_seg->base_attr .core->reverse_attr .exists() );
	assert ( pos_seg->tip_attr .core->reverse_attr .exists() );

	make_deep_connections_1d_rev ( seg, pos_seg, this, tag::mesh_is_bdry );

	add_cell_behind_below_neg_seg ( seg, pos_seg, this );                                 }

	
void Mesh::NotZeroDim::add_neg_seg  // virtual from Mesh::Core
( Cell::Negative::Segment * seg, const tag::MeshIsBdry &, const tag::DoNotBother & )

{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->reverse_attr .exists() );
	Cell::Positive::Segment * pos_seg = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Segment* > ( seg->reverse_attr .core );
	// assert that 'pos_seg' does not belong yet to 'this' mesh
	assert ( pos_seg->meshes_same_dim .find (this) == pos_seg->meshes_same_dim .end() );

	assert ( pos_seg->base_attr .exists() );
	assert ( pos_seg->tip_attr  .exists() );
	assert ( pos_seg->base_attr .core->reverse_attr .exists() );
	assert ( pos_seg->tip_attr  .core->reverse_attr .exists() );

	make_deep_connections_1d_rev ( seg, pos_seg, this, tag::mesh_is_bdry, tag::do_not_bother );

	add_cell_behind_below_neg_seg ( seg, pos_seg, this );                                       }

	
void Mesh::NotZeroDim::remove_neg_seg  // virtual from Mesh::Core
( Cell::Negative::Segment * seg, const tag::MeshIsNotBdry & )
	
{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->reverse_attr .exists() );
	Cell::Positive::Segment * pos_seg = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Segment* > ( seg->reverse_attr.core );
	// assert that 'pos_seg' belongs to 'this' mesh
	assert ( pos_seg->meshes_same_dim .find (this) != pos_seg->meshes_same_dim .end() );

	assert ( pos_seg->base_attr .exists() );
	assert ( pos_seg->tip_attr  .exists() );
	assert ( pos_seg->base_attr .core->reverse_attr .exists() );
	assert ( pos_seg->tip_attr  .core->reverse_attr .exists() );

	assert ( pos_seg->base_attr .core->reverse_attr .core->cell_behind_within .find (this) !=
	         pos_seg->base_attr .core->reverse_attr .core->cell_behind_within .end()         );
	assert
		( pos_seg->base_attr .core->reverse_attr .core->cell_behind_within [this] .core == seg );
	pos_seg->base_attr .core->reverse_attr .core->cell_behind_within .erase (this);
	assert ( pos_seg->tip_attr .core->reverse_attr .core->cell_behind_within .find (this) !=
	         pos_seg->tip_attr .core->reverse_attr .core->cell_behind_within .end()         );
	assert
		( pos_seg->tip_attr .core->reverse_attr .core->cell_behind_within [this] .core == seg );
	pos_seg->tip_attr .core->reverse_attr .core->cell_behind_within .erase (this);
	
	break_deep_connections_1d_rev ( seg, pos_seg, this, tag::mesh_is_not_bdry );

}  // end of Mesh::NotZeroDim::remove_neg_seg with tag::mesh_is_not_bdry


void Mesh::NotZeroDim::remove_neg_seg  // virtual from Mesh::Core
( Cell::Negative::Segment * seg, const tag::MeshIsNotBdry &, const tag::DoNotBother & )
	
{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->reverse_attr .exists() );
	Cell::Positive::Segment * pos_seg = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Segment* > ( seg->reverse_attr.core );
	// assert that 'pos_seg' belongs to 'this' mesh
	assert ( pos_seg->meshes_same_dim .find (this) != pos_seg->meshes_same_dim .end() );

	assert ( pos_seg->base_attr .exists() );
	assert ( pos_seg->tip_attr  .exists() );
	assert ( pos_seg->base_attr .core->reverse_attr .exists() );
	assert ( pos_seg->tip_attr  .core->reverse_attr .exists() );

	assert ( pos_seg->base_attr .core->reverse_attr .core->cell_behind_within .find (this) !=
	         pos_seg->base_attr .core->reverse_attr .core->cell_behind_within .end()         );
	assert
		( pos_seg->base_attr .core->reverse_attr .core->cell_behind_within [this] .core == seg );
	pos_seg->base_attr .core->reverse_attr .core->cell_behind_within .erase (this);
	assert ( pos_seg->tip_attr .core->reverse_attr .core->cell_behind_within .find (this) !=
	         pos_seg->tip_attr .core->reverse_attr .core->cell_behind_within .end()         );
	assert
		( pos_seg->tip_attr .core->reverse_attr .core->cell_behind_within [this] .core == seg );
	pos_seg->tip_attr .core->reverse_attr .core->cell_behind_within .erase (this);
	
	break_deep_connections_1d_rev ( seg, pos_seg, this, tag::mesh_is_not_bdry, tag::do_not_bother );

}  // end of Mesh::NotZeroDim::remove_neg_seg with tag::mesh_is_not_bdry, tag::do_not_bother


void Mesh::NotZeroDim::remove_neg_seg ( Cell::Negative::Segment * seg, const tag::MeshIsBdry & )
// virtual from Mesh::Core
	
{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->reverse_attr .exists() );
	Cell::Positive::Segment * pos_seg = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Segment* > ( seg->reverse_attr.core );
	// assert that 'pos_seg' belongs to 'this' mesh
	assert ( pos_seg->meshes_same_dim .find (this) != pos_seg->meshes_same_dim .end() );
	
	assert ( pos_seg->base_attr .exists() );
	assert ( pos_seg->tip_attr  .exists() );
	assert ( pos_seg->base_attr .core->reverse_attr .exists() );
	assert ( pos_seg->tip_attr  .core->reverse_attr .exists() );

	assert ( pos_seg->base_attr .core->reverse_attr .core->cell_behind_within .find (this) !=
	         pos_seg->base_attr .core->reverse_attr .core->cell_behind_within .end()         );
	assert
		( pos_seg->base_attr .core->reverse_attr .core->cell_behind_within [this] .core == seg );
	pos_seg->base_attr .core->reverse_attr .core->cell_behind_within .erase (this);
	assert ( pos_seg->tip_attr .core->reverse_attr .core->cell_behind_within .find (this) !=
	         pos_seg->tip_attr .core->reverse_attr .core->cell_behind_within .end()         );
	assert
		( pos_seg->tip_attr .core->reverse_attr .core->cell_behind_within [this] .core == seg );
	pos_seg->tip_attr .core->reverse_attr .core->cell_behind_within .erase (this);

	break_deep_connections_1d_rev ( seg, pos_seg, this, tag::mesh_is_bdry );

}  // end of Mesh::NotZeroDim::remove_neg_seg with tag::mesh_is_bdry
	
	
void Mesh::NotZeroDim::remove_neg_seg  // virtual from Mesh::Core
( Cell::Negative::Segment * seg, const tag::MeshIsBdry &, const tag::DoNotBother & )
	
{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->reverse_attr .exists() );
	Cell::Positive::Segment * pos_seg = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Segment* > ( seg->reverse_attr .core );
	// assert that 'pos_seg' belongs to 'this' mesh
	assert ( pos_seg->meshes_same_dim .find (this) != pos_seg->meshes_same_dim .end() );

	assert ( pos_seg->base_attr .exists() );
	assert ( pos_seg->tip_attr  .exists() );
	assert ( pos_seg->base_attr .core->reverse_attr .exists() );
	assert ( pos_seg->tip_attr  .core->reverse_attr .exists() );

	assert ( pos_seg->base_attr .core->reverse_attr .core->cell_behind_within .find (this) !=
	         pos_seg->base_attr .core->reverse_attr .core->cell_behind_within .end()         );
	assert
		( pos_seg->base_attr .core->reverse_attr .core->cell_behind_within [this] .core == seg );
	pos_seg->base_attr .core->reverse_attr .core->cell_behind_within .erase (this);
	assert ( pos_seg->tip_attr .core->reverse_attr .core->cell_behind_within .find (this) !=
	         pos_seg->tip_attr .core->reverse_attr .core->cell_behind_within .end()         );
	assert
		( pos_seg->tip_attr .core->reverse_attr .core->cell_behind_within [this] .core == seg );
	pos_seg->tip_attr .core->reverse_attr .core->cell_behind_within .erase (this);

	break_deep_connections_1d_rev ( seg, pos_seg, this, tag::mesh_is_bdry, tag::do_not_bother );

}  // end of Mesh::NotZeroDim::remove_neg_seg with tag::mesh_is_bdry, tag::do_not_bother
	
	
void Mesh::NotZeroDim::add_pos_hd_cell  // virtual from Mesh::Core
( Cell::Positive::HighDim * cll, const tag::MeshIsNotBdry & )
	
{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	// assert that 'cll' does not belong yet to 'this' mesh
	assert ( cll->meshes_same_dim .find (this) == cll->meshes_same_dim .end() );

	make_deep_connections_hd ( cll, this, tag::mesh_is_not_bdry );
	
	add_cell_behind_below_pos_hd ( cll, this );                                }


void Mesh::NotZeroDim::add_pos_hd_cell  // virtual from Mesh::Core
( Cell::Positive::HighDim * cll, const tag::MeshIsNotBdry &, const tag::DoNotBother & )
	
{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	// assert that 'cll' does not belong yet to 'this' mesh
	assert ( cll->meshes_same_dim .find (this) == cll->meshes_same_dim .end() );
	
	make_deep_connections_hd ( cll, this, tag::mesh_is_not_bdry, tag::do_not_bother );
	
	add_cell_behind_below_pos_hd ( cll, this );                                             }


void Mesh::NotZeroDim::add_pos_hd_cell  // virtual from Mesh::Core
( Cell::Positive::HighDim * cll, const tag::MeshIsBdry & )
	
{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	// assert that 'cll' does not belong yet to 'this' mesh
	assert ( cll->meshes_same_dim .find (this) == cll->meshes_same_dim .end() );
	
  make_deep_connections_hd ( cll, this, tag::mesh_is_bdry );
	
	add_cell_behind_below_pos_hd ( cll, this );                                    }


void Mesh::NotZeroDim::add_pos_hd_cell  // virtual from Mesh::Core
( Cell::Positive::HighDim * cll, const tag::MeshIsBdry &, const tag::DoNotBother & )
	
{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	// assert that 'cll' does not belong yet to 'this' mesh
	assert ( cll->meshes_same_dim .find (this) == cll->meshes_same_dim .end() );
	
  make_deep_connections_hd ( cll, this, tag::mesh_is_bdry, tag::do_not_bother );
	
	add_cell_behind_below_pos_hd ( cll, this );                                         }


void Mesh::NotZeroDim::remove_pos_hd_cell  // virtual from Mesh::Core
( Cell::Positive::HighDim * cll, const tag::MeshIsNotBdry & )

{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	// assert that 'cll' belongs to 'this' mesh
	assert ( cll->meshes_same_dim .find (this) != cll->meshes_same_dim .end() );

	Mesh bdry = cll->boundary_attr;
	assert ( bdry .core->get_dim_plus_one() + 1 == this->get_dim_plus_one() );
	Mesh::Iterator it = bdry .iterator ( tag::over_cells_of_max_dim, tag::as_they_are );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell::Core * face_p = ( *it ) .core;
		assert ( face_p->cell_behind_within .find (this) !=
		         face_p->cell_behind_within .end()         );
		assert ( face_p->cell_behind_within [this] .core == cll );
		// optimize map access !!
		face_p->cell_behind_within .erase(this);                   }

	break_deep_connections_hd ( cll, this, tag::mesh_is_not_bdry );
	
}  // end of Mesh::NotZeroDim::remove_pos_hd_cell with tag::mesh_is_not_bdry
	
	
void Mesh::NotZeroDim::remove_pos_hd_cell  // virtual from Mesh::Core
( Cell::Positive::HighDim * cll, const tag::MeshIsNotBdry &, const tag::DoNotBother & )

{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	// assert that 'cll' belongs to 'this' mesh
	assert ( cll->meshes_same_dim .find (this) != cll->meshes_same_dim .end() );

	Mesh bdry = cll->boundary_attr;
	assert ( bdry .core->get_dim_plus_one() + 1 == this->get_dim_plus_one() );
	Mesh::Iterator it = bdry .iterator ( tag::over_cells_of_max_dim, tag::as_they_are );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell::Core * face_p = ( *it ) .core;
		assert ( face_p->cell_behind_within .find (this) !=
		         face_p->cell_behind_within .end()         );
		assert ( face_p->cell_behind_within [this] .core == cll );
		// optimize map access !!
		face_p->cell_behind_within .erase (this);                  }

	break_deep_connections_hd ( cll, this, tag::mesh_is_not_bdry, tag::do_not_bother );
	
}  // end of Mesh::NotZeroDim::remove_pos_hd_cell with tag::mesh_is_not_bdry, tag::do_not_bother
	
	
void Mesh::NotZeroDim::remove_pos_hd_cell  // virtual from Mesh::Core
( Cell::Positive::HighDim * cll, const tag::MeshIsBdry & )

{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	// assert that 'cll' belongs to 'this' mesh
	assert ( cll->meshes_same_dim .find (this) != cll->meshes_same_dim .end() );

	Mesh bdry = cll->boundary_attr;
	assert ( bdry .core->get_dim_plus_one() + 1 == this->get_dim_plus_one() );
	Mesh::Iterator it = bdry .iterator ( tag::over_cells_of_max_dim, tag::as_they_are );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell::Core * face_p = ( *it ) .core;
		assert ( face_p->cell_behind_within .find (this) !=
		         face_p->cell_behind_within .end()         );
		assert ( face_p->cell_behind_within [this] .core == cll );
		// optimize map access !!
		face_p->cell_behind_within .erase (this);                  }

	break_deep_connections_hd ( cll, this, tag::mesh_is_bdry );
	
}  // end of Mesh::NotZeroDim::remove_pos_hd_cell with tag::mesh_is_bdry

	
void Mesh::NotZeroDim::remove_pos_hd_cell  // virtual from Mesh::Core
( Cell::Positive::HighDim * cll, const tag::MeshIsBdry &, const tag::DoNotBother & )

{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	// assert that 'cll' belongs to 'this' mesh
	assert ( cll->meshes_same_dim.find(this) != cll->meshes_same_dim.end() );

	Mesh bdry = cll->boundary_attr;
	assert ( bdry .core->get_dim_plus_one() + 1 == this->get_dim_plus_one() );
	Mesh::Iterator it = bdry .iterator ( tag::over_cells_of_max_dim, tag::as_they_are );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell::Core * face_p = ( *it ) .core;
		assert ( face_p->cell_behind_within .find (this) !=
		         face_p->cell_behind_within .end()         );
		assert ( face_p->cell_behind_within [this] .core == cll );
		// optimize map access !!
		face_p->cell_behind_within .erase (this);                  }

	break_deep_connections_hd ( cll, this, tag::mesh_is_bdry, tag::do_not_bother );
	
}  // end of Mesh::NotZeroDim::remove_pos_hd_cell with tag::mesh_is_bdry, tag::do_not_bother

	
void Mesh::NotZeroDim::add_neg_hd_cell  // virtual from Mesh::Core
( Cell::Negative::HighDim * cll, const tag::MeshIsNotBdry & )

{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	assert ( cll->reverse_attr .exists() );
	Cell::Positive::HighDim * pos_cll = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::HighDim* > ( cll->reverse_attr.core );
	// assert that 'pos_cll' does not belong yet to 'this' mesh
	assert ( pos_cll->meshes_same_dim .find (this) == pos_cll->meshes_same_dim .end() );

	make_deep_connections_hd_rev ( cll, pos_cll, this, tag::mesh_is_not_bdry );
	
	add_cell_behind_below_neg_hd ( cll, pos_cll, this );                                    }


void Mesh::NotZeroDim::add_neg_hd_cell  // virtual from Mesh::Core
( Cell::Negative::HighDim * cll, const tag::MeshIsNotBdry &, const tag::DoNotBother & )

{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	assert ( cll->reverse_attr .exists() );
	Cell::Positive::HighDim * pos_cll = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::HighDim* > ( cll->reverse_attr.core );
	// assert that 'pos_cll' does not belong yet to 'this' mesh
	assert ( pos_cll->meshes_same_dim .find (this) == pos_cll->meshes_same_dim .end() );

	make_deep_connections_hd_rev ( cll, pos_cll, this, tag::mesh_is_not_bdry, tag::do_not_bother );
	
	add_cell_behind_below_neg_hd ( cll, pos_cll, this );                                            }


void Mesh::NotZeroDim::add_neg_hd_cell ( Cell::Negative::HighDim * cll, const tag::MeshIsBdry & )
// virtual from Mesh::Core

{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	assert ( cll->reverse_attr .exists() );
	Cell::Positive::HighDim * pos_cll = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::HighDim* > ( cll->reverse_attr.core );
	// assert that 'pos_cll' does not belong yet to 'this' mesh
	assert ( pos_cll->meshes_same_dim .find (this) == pos_cll->meshes_same_dim .end() );

	make_deep_connections_hd_rev ( cll, pos_cll, this, tag::mesh_is_bdry );
	
	add_cell_behind_below_neg_hd ( cll, pos_cll, this );                                    }

	
void Mesh::NotZeroDim::add_neg_hd_cell  // virtual from Mesh::Core
( Cell::Negative::HighDim * cll, const tag::MeshIsBdry &, const tag::DoNotBother & )

{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	assert ( cll->reverse_attr .exists() );
	Cell::Positive::HighDim * pos_cll = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::HighDim* > ( cll->reverse_attr.core );
	// assert that 'pos_cll' does not belong yet to 'this' mesh
	assert ( pos_cll->meshes_same_dim .find (this) == pos_cll->meshes_same_dim .end() );

	make_deep_connections_hd_rev ( cll, pos_cll, this, tag::mesh_is_bdry, tag::do_not_bother );
	
	add_cell_behind_below_neg_hd ( cll, pos_cll, this );                                             }

	
void Mesh::NotZeroDim::remove_neg_hd_cell  // virtual from Mesh::Core
( Cell::Negative::HighDim * cll, const tag::MeshIsNotBdry & )
	
{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	assert ( cll->reverse_attr .exists() );
	Cell::Positive::HighDim * pos_cll = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::HighDim* > ( cll->reverse_attr .core );
	// assert that 'cll' belongs to 'this' mesh
	assert ( pos_cll->meshes_same_dim .find (this) != pos_cll->meshes_same_dim .end() );

	Mesh bdry = pos_cll->boundary_attr;
	assert ( bdry .core->get_dim_plus_one() + 1 == this->get_dim_plus_one() );
	Mesh::Iterator it = bdry .iterator ( tag::over_cells_of_max_dim, tag::as_they_are );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell::Core * face_p = ( *it ) .core;
		Cell::Core * rev_face = face_p->reverse_attr .core;
		assert ( rev_face );
		assert ( rev_face->cell_behind_within [this] .core == cll );
		rev_face->cell_behind_within .erase (this);                  }

	break_deep_connections_hd_rev ( cll, pos_cll, this, tag::mesh_is_not_bdry );

}  // end of Mesh::NotZeroDim::remove_neg_hd_cell with tag::mesh_is_not_bdry

	
void Mesh::NotZeroDim::remove_neg_hd_cell  // virtual from Mesh::Core
( Cell::Negative::HighDim * cll, const tag::MeshIsNotBdry &, const tag::DoNotBother & )
	
{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	assert ( cll->reverse_attr .exists() );
	Cell::Positive::HighDim * pos_cll = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::HighDim* > ( cll->reverse_attr .core );
	// assert that 'cll' belongs to 'this' mesh
	assert ( pos_cll->meshes_same_dim .find (this) != pos_cll->meshes_same_dim .end() );

	Mesh bdry = pos_cll->boundary_attr;
	assert ( bdry .core->get_dim_plus_one() + 1 == this->get_dim_plus_one() );
	Mesh::Iterator it = bdry .iterator ( tag::over_cells_of_max_dim, tag::as_they_are );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell::Core * face_p = ( *it ) .core;
		Cell::Core * rev_face = face_p->reverse_attr .core;
		assert ( rev_face );
		assert ( rev_face->cell_behind_within [this] .core == cll );
		rev_face->cell_behind_within .erase (this);                  }

	break_deep_connections_hd_rev ( cll, pos_cll, this, tag::mesh_is_not_bdry, tag::do_not_bother );

}  // end of Mesh::NotZeroDim::remove_neg_hd_cell with tag::mesh_is_not_bdry, tag::do_not_bother

	
void Mesh::NotZeroDim::remove_neg_hd_cell  // virtual from Mesh::Core
( Cell::Negative::HighDim * cll, const tag::MeshIsBdry & )
	
{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	assert ( cll->reverse_attr .exists() );
	Cell::Positive::HighDim * pos_cll = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::HighDim* > ( cll->reverse_attr .core );
	// assert that 'cll' belongs to 'this' mesh
	assert ( pos_cll->meshes_same_dim .find (this) != pos_cll->meshes_same_dim .end() );

	Mesh bdry = pos_cll->boundary_attr;
	assert ( bdry .core->get_dim_plus_one() + 1 == this->get_dim_plus_one() );
	Mesh::Iterator it = bdry .iterator ( tag::over_cells_of_max_dim, tag::as_they_are );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell::Core * face_p = ( *it ) .core;
		Cell::Core * rev_face = face_p->reverse_attr .core;
		assert ( rev_face );
		assert ( rev_face->cell_behind_within [this] .core == cll );
		rev_face->cell_behind_within .erase (this);                 }

	break_deep_connections_hd_rev ( cll, pos_cll, this, tag::mesh_is_bdry );

}  // end of Mesh::NotZeroDim::remove_neg_hd_cell with tag::mesh_is_bdry


void Mesh::NotZeroDim::remove_neg_hd_cell  // virtual from Mesh::Core
( Cell::Negative::HighDim * cll, const tag::MeshIsBdry &, const tag::DoNotBother & )
	
{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	assert ( cll->reverse_attr .exists() );
	Cell::Positive::HighDim * pos_cll = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::HighDim* > ( cll->reverse_attr .core );
	// assert that 'cll' belongs to 'this' mesh
	assert ( pos_cll->meshes_same_dim .find (this) != pos_cll->meshes_same_dim .end() );

	Mesh bdry = pos_cll->boundary_attr;
	assert ( bdry .core->get_dim_plus_one() + 1 == this->get_dim_plus_one() );
	Mesh::Iterator it = bdry .iterator ( tag::over_cells_of_max_dim, tag::as_they_are );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell::Core * face_p = ( *it ) .core;
		Cell::Core * rev_face = face_p->reverse_attr .core;
		assert ( rev_face );
		assert ( rev_face->cell_behind_within [this] .core == cll );
		rev_face->cell_behind_within .erase (this);              }

	break_deep_connections_hd_rev ( cll, pos_cll, this, tag::mesh_is_bdry, tag::do_not_bother );

}  // end of Mesh::NotZeroDim::remove_neg_hd_cell with tag::mesh_is_bdry, tag::do_not_bother

//-----------------------------------------------------------------------------------------//


void Mesh::ZeroDim::closed_loop ( const Cell & ver )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dim meshes cannot be closed loops" << std::endl;
	exit ( 1 );                                                     }


void Mesh::ZeroDim::closed_loop ( const Cell & ver, size_t n )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dim meshes cannot be closed loops" << std::endl;
	exit ( 1 );                                                     }


void Mesh::Connected::OneDim::closed_loop ( const Cell & ver )

// check it is indeed a closed loop !!

{	assert ( ver .is_positive() );
	this->first_ver = ver .reverse();
	this->last_ver = ver;             }
	

void Mesh::Connected::OneDim::closed_loop ( const Cell & ver, size_t n )
	
{	this->closed_loop ( ver );
	this->nb_of_segs = n;      }


void Mesh::Fuzzy::closed_loop ( const Cell & ver )
// do nothing !!
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Fuzzy meshes cannot be closed loops" << std::endl;
	exit ( 1 );                                                     }


void Mesh::Fuzzy::closed_loop ( const Cell & ver, size_t n )
// do nothing !!
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Fuzzy meshes cannot be closed loops" << std::endl;
	exit ( 1 );                                                     }

//-----------------------------------------------------------------------------------------//


Mesh::Core * Mesh::ZeroDim::build_deep_copy ( )
// virtual from Mesh::Core, here execution forbidden
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "No deep copy for zero-dim meshes" << std::endl;
	exit ( 1 );                                                     }
	
Mesh::Core * Mesh::NotZeroDim::build_deep_copy ( )
// virtual from Mesh::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "build_deep_copy not implemented yet" << std::endl;
	exit ( 1 );                                                     }

//-----------------------------------------------------------------------------------------//


Mesh Mesh::convert_to
( const tag::Connected &, const tag::OneDim &, const tag::SurelyExists & ) const

{	assert ( this->dim() == 1 );
	assert ( this->number_of ( tag::segments ) > 0 );

	Mesh::Iterator it = this->iterator ( tag::over_vertices );
	it .reset();  assert ( it .in_range() );
	Cell A = *it, start = A;
	Cell stop ( tag::non_existent );
	size_t nb_of_segs = 0;
	bool closed_loop = false;
	while ( true )
	{	Cell seg = this->cell_in_front_of ( A, tag::may_not_exist );
		if ( not seg .exists() )  // found the "stop"
		{	stop = A;  break;  }
		nb_of_segs ++;
		Cell B = seg .tip();
		if ( B == start )  {  closed_loop = true;  break;  }
		A = B;                                                       }

	if ( not closed_loop )
	{	assert ( stop .exists() );
		while ( true )
		{	Cell seg = this->cell_behind ( start, tag::may_not_exist );
			if ( not seg .exists() )  break;  // this is the true "start"
			nb_of_segs ++;
			start = seg .base() .reverse();                              }  }

	assert ( nb_of_segs == this->number_of ( tag::segments ) );
	
	if ( closed_loop )
	{	assert ( not stop .exists() );
		assert ( nb_of_segs == this->number_of ( tag::segments ) );
		stop = start;                                               }
		
	Mesh::Connected::OneDim * res_core =
	    new Mesh::Connected::OneDim ( tag::with, nb_of_segs, tag::segments,
	                                  tag::one_dummy_wrapper               );
	Mesh result ( tag::whose_core_is, res_core, tag::freshly_created, tag::is_positive );

	Mesh::Iterator itt = this->iterator ( tag::over_segments );
	for ( itt .reset(); itt .in_range(); itt++ )
		(*itt) .add_to_mesh ( result, tag::do_not_bother );

	res_core->first_ver = start .reverse();
	res_core->last_ver = stop;
	return result;                                                                         }

//-----------------------------------------------------------------------------------------//


Mesh Mesh::convert_to
( const tag::Connected &, const tag::OneDim &, const tag::MayNotExist & ) const

{	assert ( this->dim() == 1 );
	assert ( this->number_of ( tag::segments ) > 0 );

	Mesh::Iterator it = this->iterator ( tag::over_vertices );
	it .reset();  assert ( it .in_range() );
	Cell A = *it, start = A;
	Cell stop ( tag::non_existent );
	size_t nb_of_segs = 0;
	bool closed_loop = false;
	while ( true )
	{	Cell seg = this->cell_in_front_of ( A, tag::may_not_exist );
		if ( not seg .exists() )  // found the "stop"
		{	stop = A;  break;  }
		nb_of_segs ++;
		Cell B = seg .tip();
		if ( A == B )  // closed loop
		{	closed_loop = true;  break;  }
		A = B;                                                       }

	if ( closed_loop )
	{	assert ( not stop .exists() );
		if ( nb_of_segs != this->number_of ( tag::segments ) )
			return Mesh ( tag::non_existent );
		Mesh::Connected::OneDim * res_core =
		    new Mesh::Connected::OneDim ( tag::with, nb_of_segs, tag::segments,
		                                  tag::one_dummy_wrapper               );
		Mesh result ( tag::whose_core_is, res_core,
		    tag::freshly_created, tag::is_positive );
		res_core->first_ver = A;
		res_core->last_ver = A;
		return result;                                                             }
	assert ( stop .exists() );

	while ( true )
	{	Cell seg = this->cell_behind ( start, tag::may_not_exist );
		if ( not seg .exists() )  break;  // this is the true "start"
		nb_of_segs ++;
		start = seg .base() .reverse();                             }

	if ( nb_of_segs != this->number_of ( tag::segments ) )
		return Mesh ( tag::non_existent );

	Mesh::Connected::OneDim * res_core =
	    new Mesh::Connected::OneDim ( tag::with, nb_of_segs, tag::segments,
	                                  tag::one_dummy_wrapper               );
	Mesh result ( tag::whose_core_is, res_core, tag::freshly_created, tag::is_positive );

	Mesh::Iterator itt = this->iterator ( tag::over_segments );
	for ( itt .reset(); itt .in_range(); itt++ )
		(*itt) .add_to_mesh ( result, tag::do_not_bother );
	
	res_core->first_ver = start .reverse();
	res_core->last_ver = stop;
	return result;                                                                        }


//-----------------------------------------------------------------------------------------//


// for one-dimensional meshes, if you call 'baricenter' on an extremity
// nothing will happen

// in contrast, for two or more dimensions, 'baricenter' will act even
// on a vertex on the boundary of 'this' mesh
// depending on the current working manifold, the resulting coordinates
// may be projected on the boundary or may lean towards the interior of the mesh

void Mesh::baricenter ( const Cell & ver )

// 'ver' is a vertex in 'this' mesh

{	assert ( ver .dim() == 0 );
	std::vector < Cell > neighbours;  // vertices
	size_t n = 0;
	if ( this->dim() == 1 )
	{	Cell front = this->cell_in_front_of ( ver, tag::may_not_exist );
		if ( not front .exists() ) return;
		assert ( front .base() == ver .reverse() );
		neighbours .push_back ( front .tip() );
		Cell back = this->cell_behind ( ver, tag::may_not_exist );
		if ( not back .exists() ) return;
		assert ( back .tip() == ver );
		neighbours .push_back ( back .base() .reverse() );
		n = 2;                                                           }
	else
	{	assert ( this->dim() >= 2 );
		Mesh::Iterator it = this->iterator ( tag::over_vertices, tag::around, ver );
		for ( it .reset(); it .in_range(); it++ )
		{	n++;  neighbours .push_back ( *it );  }                                   }
	assert ( n == neighbours .size() );
	assert ( n >= 2 );
	std::vector < double > coefs ( n, 1. / double(n) );
	Manifold::working .interpolate ( ver, coefs, neighbours );                      }

//-----------------------------------------------------------------------------------------//


void Mesh::baricenter ( const Cell & ver, const tag::Winding & )

// 'ver' is a vertex in 'this' mesh
// tag::winding is a mere indication that we are on a quotient manifold
// and, as such, neighbour segments may have winding

{	Manifold space = Manifold::working;
	assert ( space .exists() );  // we use the current (quotient) manifold
	Manifold::Quotient * mani_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space .core );
	Function coords_q = space .coordinates();
	Manifold mani_Eu = mani_q->base_space;  // underlying Euclidian manifold
	Function coords_Eu = mani_Eu .coordinates();

	assert ( ver .dim() == 0 );
	std::vector < Cell > neighbours;  // vertices
	size_t n = 0;
	if ( this->dim() == 1 )
	{	Cell front = this->cell_in_front_of ( ver, tag::may_not_exist );
		if ( not front .exists() ) return;
		assert ( front .base() == ver.reverse() );
		Cell back = this->cell_behind ( ver, tag::may_not_exist );
		if ( not back .exists() ) return;
		assert ( back .tip() == ver );
		Cell shadow_front ( tag::vertex );
		Cell shadow_back ( tag::vertex );
		if ( coords_Eu.nb_of_components() == 1 )
		{	double new_co = coords_q ( front.tip(), tag::winding, front.winding() );
			coords_Eu ( shadow_front ) = new_co;
			new_co = coords_q
				( back .base() .reverse(), tag::winding, back .reverse() .winding() );
			coords_Eu ( shadow_back ) = new_co;                                       }
		else
		{	assert ( coords_Eu .nb_of_components() > 1 );
			std::vector < double > new_co =
				coords_q ( front .tip(), tag::winding, front.winding() );
			coords_Eu ( shadow_front ) = new_co;
			new_co = coords_q
				( back .base() .reverse(), tag::winding, back .reverse() .winding() );
			coords_Eu ( shadow_back ) = new_co;                                       }
		neighbours .push_back ( shadow_front );
		neighbours .push_back ( shadow_back );
		n = 2;                                                                                 }
	else
	{	assert ( this->dim() >= 2 );
		Mesh::Iterator it = this->iterator ( tag::over_segments, tag::around, ver.reverse() );
		for ( it .reset(); it .in_range(); it++ )
		{	n++;
			Cell seg = *it;
			Cell shadow ( tag::vertex );
			if ( coords_Eu .nb_of_components() == 1 )
			{	double new_co = coords_q ( seg .tip(), tag::winding, seg .winding() );
				coords_Eu ( shadow ) = new_co;                                         }
			else
			{	assert ( coords_Eu .nb_of_components() > 1 );
				std::vector < double > new_co =
					coords_q ( seg .tip(), tag::winding, seg .winding() );
				coords_Eu ( shadow ) = new_co;                           }
			neighbours .push_back ( shadow );                                                  }  }
	std::vector < double > co = coords_Eu ( ver );
	assert ( n == neighbours .size() );
	assert ( n >= 2 );
	std::vector < double > coefs ( n, 1. / double(n) );
	mani_Eu.interpolate ( ver, coefs, neighbours );
}

//-----------------------------------------------------------------------------------------//


void Mesh::baricenter ( const Cell & ver, const tag::Winding &,
                        const tag::ShadowVertices &, const std::vector < Cell > & vec_cll )

// 'ver' is a vertex in 'this' mesh
// tag::winding is a mere indication that we are on a quotient manifold
// and, as such, neighbour segments may be winding

// the above version (without shadow vertices as argument)
// builds new vertices each time it is invoked,
// then destroys them (ifdef MANIFEM_COLLECT_CM)
// this version is more efficient

{	Manifold space = Manifold::working;
	assert ( space .exists() );  // we use the current (quotient) manifold
	Manifold::Quotient * mani_q = tag::Util::assert_cast
		< Manifold::Core*, Manifold::Quotient* > ( space .core );
	Function coords_q = space .coordinates();
	Manifold mani_Eu = mani_q->base_space;  // underlying Euclidian manifold
	Function coords_Eu = mani_Eu .coordinates();

	assert ( ver .dim() == 0 );
	std::vector < Cell > neighbours;  // vertices
	size_t n = 0;
	if ( this->dim() == 1 )
	{	Cell front = this->cell_in_front_of ( ver, tag::may_not_exist );
		if ( not front .exists() ) return;
		assert ( front .base() == ver .reverse() );
		Cell back = this->cell_behind ( ver, tag::may_not_exist );
		if ( not back .exists() ) return;
		assert ( back .tip() == ver );
		assert ( vec_cll .size() >= 2 );
		Cell shadow_front = vec_cll [0];
		Cell shadow_back  = vec_cll [1];
		if ( coords_Eu .nb_of_components() == 1 )
		{	double new_co = coords_q ( front .tip(), tag::winding, front .winding() );
			coords_Eu ( shadow_front ) = new_co;
			new_co = coords_q
				( back .base() .reverse(), tag::winding, back .reverse() .winding() );
			coords_Eu ( shadow_back ) = new_co;                                        }
		else
		{	assert ( coords_Eu .nb_of_components() > 1 );
			std::vector < double > new_co =
				coords_q ( front .tip(), tag::winding, front .winding() );
			coords_Eu ( shadow_front ) = new_co;
			new_co = coords_q
				( back .base() .reverse(), tag::winding, back .reverse() .winding() );
			coords_Eu ( shadow_back ) = new_co;                                       }
		neighbours .push_back ( shadow_front );
		neighbours .push_back ( shadow_back );
		n = 2;                                                                          }
	else
	{	assert ( this->dim() >= 2 );
		Mesh::Iterator it = this->iterator ( tag::over_segments, tag::around, ver .reverse() );
		for ( it .reset(); it .in_range(); it++ )
		{	Cell seg = *it;
			assert ( vec_cll .size() > n );
			Cell shadow = vec_cll [n];
			n++;
			if ( coords_Eu .nb_of_components() == 1 )
			{	double new_co = coords_q ( seg .tip(), tag::winding, seg .winding() );
				coords_Eu ( shadow ) = new_co;                                         }
			else
			{	assert ( coords_Eu .nb_of_components() > 1 );
				std::vector < double > new_co =
					coords_q ( seg .tip(), tag::winding, seg .winding() );
				coords_Eu ( shadow ) = new_co;                            }
			neighbours .push_back ( shadow );                                            }        }
	assert ( n == neighbours .size() );
	assert ( n >= 2 );
	std::vector < double > coefs ( n, 1. / double(n) );
	mani_Eu .interpolate ( ver, coefs, neighbours );                                             }


//-----------------------------------------------------------------------------//


