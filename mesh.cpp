
// maniFEM mesh.cpp 2020.01.03

//    This file is part of maniFEM, a C++ library for meshes on manifolds and finite elements.

//    ManiFEM is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Lesser General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.

//    ManiFEM is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.

//    You should have received a copy of the GNU Lesser General Public License
//    along with maniFEM.  If not, see <https://www.gnu.org/licenses/>.

//    Copyright 2019, 2020 Cristian Barbarosie cristian.barbarosie@gmail.com
//    https://github.com/cristian-barbarosie/manifem

#include <forward_list>

#include "mesh.h"

using namespace maniFEM;


size_t Mesh::maximum_dimension_plus_one { 4 };  // static data member

// we keep here the topological dimension of the largest mesh we intend to build
// '4' means three-dimensional meshes (cubes, tetrahedra, etc)
// '3' means two-dimensional meshes, including surfaces in R^3
// '2' would be for just polygonal lines
// '1' doesn't make much sense - just points ?

// static data members :
std::vector < size_t > Cell::double_heap_size_pos ( Mesh::maximum_dimension_plus_one, 0. );
std::vector < size_t > Cell::double_heap_size_neg ( Mesh::maximum_dimension_plus_one, 0. );
std::vector < size_t > Cell::size_t_heap_size_pos ( Mesh::maximum_dimension_plus_one, 0 );
std::vector < size_t > Cell::size_t_heap_size_neg ( Mesh::maximum_dimension_plus_one, 0 );
std::vector < size_t > Cell::short_int_heap_size_pos ( Mesh::maximum_dimension_plus_one, 0 );
std::vector < size_t > Cell::short_int_heap_size_neg ( Mesh::maximum_dimension_plus_one, 0 );

	
Cell::Core * const Cell::ghost { new Cell::Negative::Vertex ( tag::ghost ) };
// static data member, see paragraph 8.14 in the manual

//-----------------------------------------------------------------------------//


const Mesh::Methods Mesh::Positive::methods_pos  // static data member
{	& Mesh::Positive::is_positive, & Mesh::Positive::reverse  };

const Mesh::Methods & Mesh::Positive::get_meth_pos ( ) // virtual from Mesh::Core
{	return Mesh::Positive::methods_pos;  }
	
const Mesh::Methods Mesh::OneDim::Positive::methods_pos  // static data member
{	& Mesh::Positive::is_positive, & Mesh::Positive::reverse  };

const Mesh::Methods & Mesh::OneDim::Positive::get_meth_pos ( ) // virtual from Mesh::Core
{	return Mesh::OneDim::Positive::methods_pos;  }
	
const Mesh::Methods Mesh::Positive::methods_neg  // static data member
{ & Mesh::Negative::is_positive, & Mesh::Negative::reverse  };

const Mesh::Methods & Mesh::Positive::get_meth_neg ( ) // virtual from Mesh::Core
{	return Mesh::Positive::methods_neg;  }
	
const Mesh::Methods Mesh::OneDim::Positive::methods_neg  // static data member
{ & Mesh::Negative::is_positive, & Mesh::Negative::reverse  };

const Mesh::Methods & Mesh::OneDim::Positive::get_meth_neg ( ) // virtual from Mesh::Core
{	return Mesh::OneDim::Positive::methods_neg;  }
	
//-----------------------------------------------------------------------------//


bool Cell::Core::Positive::is_positive ( ) const  // virtual from Cell::Core
{	return true;  }

bool Cell::Core::Negative::is_positive ( ) const  // virtual from Cell::Core
{	return false;  }

bool Mesh::Positive::is_positive ( ) // static
{	return true;  }

bool Mesh::Negative::is_positive ( ) // NegativeMesh is a namespace
{	return false;  }

Cell::Core::Positive * Cell::Core::Positive::get_positive ( )  // virtual from Cell::Core
{	return this;  }

Cell::Core::Positive * Cell::Core::Negative::get_positive ( )  // virtual from Cell::Core
{	assert ( this->reverse_p );
	assert ( this->reverse_p->is_positive() );
	return (Cell::Core::Positive*) this->reverse_p;     }

//-----------------------------------------------------------------------------//


size_t Cell::Positive::Vertex::get_dim ( ) const  // virtual from Cell::Core
{	return 0;  }

size_t Cell::Negative::Vertex::get_dim ( ) const  // virtual from Cell::Core
{	return 0;  }

size_t Cell::Positive::Segment::get_dim ( ) const  // virtual from Cell::Core
{	return 1;  }

size_t Cell::Negative::Segment::get_dim ( ) const  // virtual from Cell::Core
{	return 1;  }

size_t Cell::Positive::get_dim ( ) const  // virtual from Cell::Core
{	assert ( this->boundary_p );
	return this->boundary_p->get_dim_plus_one();  }

size_t Cell::Negative::get_dim ( ) const  // virtual from Cell::Core
{	assert ( this->reverse_p );
	return this->reverse_p->get_dim();  }


size_t Mesh::Positive::get_dim_plus_one ( )  // virtual from Mesh::Core
{	return this->cells.size();  }

size_t Mesh::OneDim::Positive::get_dim_plus_one ( )  // virtual from Mesh::Core
{	return 2;  }

	
size_t Mesh::Positive::number_of ( const tag::CellsOfDim &, size_t d )
// virtual from Mesh::Core
{	assert ( d < this->get_dim_plus_one() );
	return this->cells[d].size();             }
	
size_t Mesh::OneDim::Positive::number_of ( const tag::CellsOfDim &, size_t d )
// virtual from Mesh::Core
{	assert ( d <= 1 );
	return this->cells[d].size();  }  // will change


Cell::Core * Mesh::Positive::first_vertex ( )  // virtual from Mesh::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "only one-dimensional meshes have first vertex" << std::endl;
	exit ( 1 );                                                                                     }

Cell::Core * Mesh::Positive::last_vertex ( )  // virtual from Mesh::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "only one-dimensional meshes have last vertex" << std::endl;
	exit ( 1 );                                                                                     }

Cell::Core * Mesh::Positive::first_segment ( )  // virtual from Mesh::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "only one-dimensional meshes have first segment" << std::endl;
	exit ( 1 );                                                                                     }

Cell::Core * Mesh::Positive::last_segment ( )  // virtual from Mesh::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "only one-dimensional meshes have last segment" << std::endl;
	exit ( 1 );                                                                                     }

Cell::Core * Mesh::OneDim::Positive::first_vertex ( )  // virtual from Mesh::Core
// returns a positive vertex
{	this->order();
	assert ( this->first_ver );
	assert ( this->first_ver != Cell::ghost );
	assert ( not this->first_ver->is_positive() );
	assert ( this->last_ver );
	assert ( this->last_ver->is_positive() );
	return this->first_ver->reverse_p;              }

Cell::Core * Mesh::OneDim::Positive::last_vertex ( )  // virtual from Mesh::Core
{	this->order();
	assert ( this->first_ver );
	assert ( this->first_ver != Cell::ghost );
	assert ( not this->first_ver->is_positive() );
	assert ( this->last_ver );
	assert ( this->last_ver->is_positive() );
	return this->last_ver;                          }

Cell::Core * Mesh::OneDim::Positive::first_segment ( )  // virtual from Mesh::Core
{	this->order();
	assert ( this->first_ver );
	assert ( this->first_ver != Cell::ghost );
	assert ( not this->first_ver->is_positive() );
	assert ( this->last_ver );
	assert ( this->last_ver->is_positive() );
	Cell::Core * neg_ver = this->first_ver;
	assert ( neg_ver );  assert ( not neg_ver->is_positive() );
	std::map<Mesh::Core*,Cell::Core*>::const_iterator
		it = neg_ver->cell_behind_within.find ( this );
	assert ( it != neg_ver->cell_behind_within.end() );
	assert ( it->second );  // check what happens for an empty mesh
	return it->second;                                               }

Cell::Core * Mesh::OneDim::Positive::last_segment ( )  // virtual from Mesh::Core
{	Cell::Core * ver = this->last_vertex();
	assert ( ver );  assert ( ver->is_positive() );
	std::map<Mesh::Core*,Cell::Core*>::const_iterator
		it = ver->cell_behind_within.find ( this );
	assert ( it != ver->cell_behind_within.end() );
	assert ( it->second );  // check what happens for an empty mesh
	return it->second;                                               }


	
//-----------------------------------------------------------------------------//


bool Cell::Core::Positive::belongs_to ( Mesh::Core * msh, const tag::Oriented & ) const
// virtual from Cell::Core

{	assert ( msh->get_dim_plus_one() == this->get_dim() + 1 );
	const std::map < Mesh::Core *, Cell::field_to_meshes > & mmap = this->meshes[0];
	std::map<Mesh::Core*,Cell::field_to_meshes>::const_iterator it = mmap.find(msh);
	if ( it == mmap.end() ) return false;
	const Cell::field_to_meshes & field = it->second;
	if ( field.counter_pos != 1 ) return false;
	assert ( field.counter_neg == 0 );
	return true;                                                                    }


bool Cell::Core::Positive::belongs_to ( Mesh::Core * msh, const tag::NotOriented & ) const
// virtual from Cell::Core

{	const std::map < Mesh::Core *, Cell::field_to_meshes > & mmap =
		this->meshes [ Mesh::diff ( msh->get_dim_plus_one(), this->get_dim() + 1 ) ];
	std::map<Mesh::Core*,Cell::field_to_meshes>::const_iterator it = mmap.find(msh);
	return ( it != mmap.end() );                                                     }


bool Cell::Core::Negative::belongs_to ( Mesh::Core * msh, const tag::Oriented & ) const
// virtual from Cell::Core

{	assert ( msh->get_dim_plus_one() == this->get_dim() + 1 );
	assert ( this->reverse_p );
	assert ( this->reverse_p->is_positive() );
	Cell::Core::Positive * rev = static_cast < Cell::Core::Positive * > ( this->reverse_p );
	const std::map < Mesh::Core *, Cell::field_to_meshes > & mmap = rev->meshes[0];
	std::map<Mesh::Core*,Cell::field_to_meshes>::const_iterator it = mmap.find(msh);
	if ( it == mmap.end() ) return false;
	const Cell::field_to_meshes & field = it->second;
	if ( field.counter_neg != 1 ) return false;
	assert ( field.counter_pos == 0 );
	return true;                                                                               }


bool Cell::Core::Negative::belongs_to ( Mesh::Core * msh, const tag::NotOriented & ) const
// virtual from Cell::Core

{	assert ( this->reverse_p );
	assert ( this->reverse_p->is_positive() );
	Cell::Core::Positive * rev = static_cast < Cell::Core::Positive * > ( this->reverse_p );
	std::map < Mesh::Core *, Cell::field_to_meshes > & mmap =
		rev->meshes [ Mesh::diff ( msh->get_dim_plus_one(), this->get_dim() + 1 ) ];
	std::map<Mesh::Core*,Cell::field_to_meshes>::iterator it = mmap.find(msh);
	return ( it != mmap.end() );                                                                }

//-----------------------------------------------------------------------------//


Cell::Core * Cell::Positive::Vertex::reverse ( const tag::BuildIfNotExists & )
// virtual from Cell::Core

{	if ( this->reverse_p == nullptr )
		this->reverse_p = new Cell::Negative::Vertex ( tag::reverse_of, this );
	return this->reverse_p;                                                    }


Cell::Core * Cell::Positive::Segment::reverse ( const tag::BuildIfNotExists & )
// virtual from Cell::Core

{	if ( this->reverse_p == nullptr )
		this->reverse_p = new Cell::Negative::Segment ( tag::reverse_of, this );
	return this->reverse_p;                                                     }


Cell::Core * Cell::Positive::reverse ( const tag::BuildIfNotExists & )
// virtual from Cell::Core

{	if ( this->reverse_p == nullptr )
		this->reverse_p = new Cell::Negative ( tag::reverse_of, this );
	return this->reverse_p;                                                   }


Cell::Core * Cell::Core::Negative::reverse ( const tag::BuildIfNotExists & )
// virtual from Cell::Core

{	assert ( this->reverse_p );
	assert ( this->reverse_p->is_positive() );
	return this->reverse_p;                      }


Mesh Mesh::Positive::reverse ( Mesh::Core * core )  // static

// negative meshes have no core, the core points to the (positive) reverse mesh

{	return Mesh ( tag::whose_core_is, core, tag::is_negative, tag::build_cells_if_necessary );  }


Mesh Mesh::Negative::reverse ( Mesh::Core * core )
// NegativeMesh is a namespace

// negative meshes have no core, the core points to the (positive) reverse mesh

{	return Mesh ( tag::whose_core_is, core, tag::is_positive );  }

//-----------------------------------------------------------------------------//


Cell::Core * Cell::Core::tip ()  // virtual
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Only segments have extremities." << std::endl;
	exit ( 1 );                                                                                   }

Cell::Core * Cell::Core::base ()  // virtual
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Only segments have extremities." << std::endl;
	exit ( 1 );                                                                                   }

Cell::Core * Cell::Positive::Segment::tip () { return this->tip_p;  }
// virtual from Cell::Core

Cell::Core * Cell::Positive::Segment::base () { return this->base_p;  }
// virtual from Cell::Core

Cell::Core * Cell::Negative::Segment::tip ()  // virtual, overrides definition by Cell::Core
{	assert ( this->reverse_p );
	Cell::Positive::Segment * pos_seg = (Cell::Positive::Segment*) this->reverse_p;
	assert ( pos_seg->base_p );
	// assert ( pos_seg->base_p->reverse_p );
	return pos_seg->base_p->reverse_p;                                               }

Cell::Core * Cell::Negative::Segment::base ()  // virtual, overrides definition by Cell::Core
{	assert ( this->reverse_p );
	Cell::Positive::Segment * pos_seg = (Cell::Positive::Segment*) this->reverse_p;
	assert ( pos_seg->tip_p );
	// assert ( pos_seg->tip_p->reverse_p );
	return pos_seg->tip_p->reverse_p;                                                }

//-----------------------------------------------------------------------------//

#ifndef NDEBUG


std::string Cell::Core::Positive::get_name ()  // virtual from Cell::Core
{	return this->name;  }

std::string Cell::Core::Negative::get_name ()  // virtual from Cell::Core
{	return "r" + this->reverse_p->name;  }


void Cell::Positive::Vertex::print_everything ( )  // virtual from Cell::Core

{	std::cout << "this is Positive::Vertex " << this->name << std::endl;
	if ( this->meshes.size() > 0 )
	{	if ( this->meshes[0].size() > 0 )
				std::cout << "meshes of index 0, dim 0 (segments disguised as meshes)" << std::endl;
		std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
			it = this->meshes[0].begin(), it_e = this->meshes[0].end();
		for ( ; it != it_e; it++ )
		{	Cell::Positive::Segment * seg = (Cell::Positive::Segment*) it->first;
			std::cout << seg->name << " " << it->second.counter_pos << " " << it->second.counter_neg << "  ";  }
			if ( this->meshes[0].size() > 0 )  std::cout << std::endl;
		for ( size_t d = 1; d < this->meshes.size(); d++ )
		{	if ( this->meshes[d].size() > 0 )
				std::cout << "meshes of index " << d << ", that is, dimension " << d << std::endl;
			std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
							itt = this->meshes[d].begin(), itt_e = this->meshes[d].end();
			for ( ; itt != itt_e; itt++ )
				std::cout << itt->first->get_name() << " " << itt->second.counter_pos << " " << itt->second.counter_neg << "  ";
			if ( this->meshes[d].size() > 0 )  std::cout << std::endl;                 }  }
	if ( this->reverse_p ) std::cout << "has reverse" << std::endl;                       }


void Cell::Negative::Vertex::print_everything ( )  // virtual from Cell::Core

{	std::cout << "this is Negative::Vertex, reverse of " << this->reverse_p->name << std::endl;    }


void Cell::Positive::Segment::print_everything ( )  // virtual from Cell::Core

{	std::cout << "this is Positive::Segment " << this->name << std::endl;
	std::cout << "base : " << this->base_p->get_name() << std::endl;
	std::cout << "tip :  " << this->tip_p->get_name()  << std::endl;
	if ( this->meshes.size() > 0 )
	{	for ( size_t d = 0; d < this->meshes.size(); d++ )
		{	if ( this->meshes[d].size() > 0 )
				std::cout << "meshes of index " << d << ", that is, dimension " << d+1 << std::endl;
			std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
							it = this->meshes[d].begin(), it_e = this->meshes[d].end();
			for ( ; it != it_e; it++ ) std::cout << it->first->get_name() << " " << it->second.counter_pos << " " << it->second.counter_neg << "  ";
			if ( this->meshes[d].size() > 0 )  std::cout << std::endl;                 }  }
	if ( this->reverse_p ) std::cout << "has reverse" << std::endl;                       }


void Cell::Negative::Segment::print_everything ( )  // virtual from Cell::Core

{	std::cout << "this is Negative::Segment " << this->get_name() << std::endl;     }


void Cell::Positive::print_everything ( )  // virtual from Cell::Core
	
{	size_t dim = this->get_dim();
	std::cout << "this is PositiveCell of dim " << dim << " " << this->name << std::endl;
	if ( this->meshes.size() > 0 )
	{	std::cout << "meshes above me :" << std::endl;
		for ( size_t d = 0; d < this->meshes.size(); d++ )
		{	if ( this->meshes[d].size() > 0 )
				std::cout << "meshes of index " << d << ", that is, dimension " << d+dim << std::endl;
			std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
							it = this->meshes[d].begin(), it_e = this->meshes[d].end();
			for ( ; it != it_e; it++ ) std::cout << it->first->get_name() << " " << it->second.counter_pos
		                                       << " " << it->second.counter_neg << "  ";
			if ( this->meshes[d].size() > 0 )  std::cout << std::endl;                 }  }
	if ( this->reverse_p ) std::cout << "has reverse" << std::endl;                       }


void Cell::Negative::print_everything ( )  // virtual from Cell::Core
	
{	std::cout << "this is NegativeCell " << this->reverse_p->get_name() << std::endl;   }


std::string Mesh::Positive::get_name ( )  // virtual from Mesh::Core

{	if ( this->cell_enclosed ) return "bdry_of_" + this->cell_enclosed->name;
	else return this->name;                                                   }

std::string Mesh::OneDim::Positive::get_name ( )  // virtual from Mesh::Core

{	if ( this->cell_enclosed ) return "bdry_of_" + this->cell_enclosed->name;
	else return this->name;                                                   }

void Mesh::Positive::print_everything ( )  // virtual from Mesh::Core

{	std::cout << "this is Mesh::Positive " << this->get_name() << std::endl;
	for ( size_t d = 0; d < this->cells.size(); d++ )
	{	std::cout << "cells of dim " << d << " :" << std::endl;
		int counter = 0;
		std::list<Cell::Core*>::iterator it = this->cells[d].begin(),
		                                 it_e = this->cells[d].end();
		for ( ; it != it_e; it++, counter++ ) std::cout << (*it)->get_name() << " ";
		if ( counter > 0 )  std::cout << std::endl;                                   }  }

void Mesh::OneDim::Positive::print_everything ( )  // virtual from Mesh::Core

{	std::cout << "this is Mesh::OneDim::Positive " << this->get_name() << std::endl;
	for ( size_t d = 0; d < this->cells.size(); d++ )
	{	std::cout << "cells of dim " << d << " :" << std::endl;
		int counter = 0;
		std::list<Cell::Core*>::iterator it = this->cells[d].begin(),
		                                 it_e = this->cells[d].end();
		for ( ; it != it_e; it++, counter++ ) std::cout << (*it)->get_name() << " ";
		if ( counter > 0 )  std::cout << std::endl;                                   }  }

#endif

 //////////////////////////////////////////////////////////////////
////////////    add/remove a cell to/from a mesh    ///////////////
///////////////////////////////////////////////////////////////////
////////////     this is where the magic happens     //////////////
//////////// (and also where we get our hands dirty) //////////////
//////////////////////////////////////////////////////////////////


void Cell::Positive::Vertex::glue_on_my_bdry ( Cell::Core * )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Vertices have no boundary." << std::endl;
	exit ( 1 );                                                                                   }


void Cell::Positive::Segment::glue_on_my_bdry ( Cell::Core * ver )
// virtual from Cell::Core

{	assert ( ver->get_dim() == 0 );
	ver->add_to ( (Mesh::Core*) this );
	// ver.meshes[0] contains segments disguised as zero-dimensional meshes
	// 'add_to' is virtual, so the computer will choose the right version
	// (Positive::Vertex::add_to or Negative::Vertex::add_to)
	this->glue_common ( ver );                        }
	

void Cell::Positive::glue_on_my_bdry ( Cell::Core * face )
// virtual from Cell::Core

{	assert ( this->get_dim() == face->get_dim() + 1 );
	face->add_to ( this->boundary_p );
	// 'add_to' is virtual, so the computer will choose the right version
	// ( Positive::Segment::add_to or Negative::Segment::add_to
	//   or PositiveCell::add_to or NegativeCell::add_to    )
	this->glue_common ( face );                           }


void Cell::Core::Negative::glue_on_my_bdry ( Cell::Core * cll )
// virtual from Cell::Core

{	assert ( cll->reverse_p );
	assert ( this->reverse_p );
	this->reverse_p->glue_on_my_bdry ( cll->reverse_p );  }


void Cell::Positive::Vertex::cut_from_my_bdry ( Cell::Core * )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Vertices have no boundary." << std::endl;
	exit ( 1 );                                                                                   }


void Cell::Positive::Segment::cut_from_my_bdry ( Cell::Core * ver )
// virtual from Cell::Core

{	assert ( ver->get_dim() == 0 );
	ver->remove_from ( (Mesh::Core*) this );
	// ver.meshes[0] contains segments disguised as zero-dimensional meshes
	// 'remove_from' is virtual, so the computer will choose the right version
	// (Positive::Vertex::remove_from or Negative::Vertex::remove_from)
	this->cut_common ( ver );                           }
	

void Cell::Positive::cut_from_my_bdry ( Cell::Core * face )
// virtual from Cell::Core

{	assert ( this->get_dim() == face->get_dim() + 1 );
	face->remove_from ( this->boundary_p );
	// 'remove_from' is virtual, so the computer will choose the right version
	// ( Positive::Segment::remove_from or Negative::Segment::remove_from
	//   or PositiveCell::remove_from or NegativeCell::remove_from    )
	this->cut_common ( face );                            }


void Cell::Core::Negative::cut_from_my_bdry ( Cell::Core * cll )
// virtual from Cell::Core

{	assert ( cll->reverse_p );
	assert ( this->reverse_p );
	this->reverse_p->cut_from_my_bdry ( cll->reverse_p );  }

//-----------------------------------------------------------------------------//


void Cell::Positive::Vertex::add_to ( Mesh::Core * msh ) // virtual from Cell::Core
	
{	assert ( msh );
	Cell::Positive::Segment * seg = ( Cell::Positive::Segment * ) msh;
	// ver.meshes[0] contains segments disguised as zero-dimensional meshes
	assert ( this->meshes.size() > 0 );
	// assert that 'this' vertex does not belong yet to the mesh 'msh'
	assert ( this->meshes[0].find(msh) == this->meshes[0].end() );
	// '[0]' means meshes of the same dimension as the cell
	assert ( seg->tip_p == nullptr );
	seg->tip_p = this;
	// this->meshes[0][msh] = Cell::field_to_meshes { 1, 0 };
	// the third component 'where' is irrelevant here
	this->meshes[0].emplace ( std::piecewise_construct,
		std::forward_as_tuple(msh), std::forward_as_tuple(1,0) );
	seg->deep_connections ( this, Mesh::action_add );                 }


void Cell::Positive::Vertex::remove_from ( Mesh::Core * msh ) // virtual from Cell::Core
	
{	assert ( msh );
	Cell::Positive::Segment * seg = ( Cell::Positive::Segment * ) msh;
	// ver.meshes[0] contains segments disguised as zero-dimensional meshes
	assert ( this->meshes.size() > 0 );
	// assert that 'this' segment belongs to the mesh 'msh'
	assert ( this->meshes[0].find(msh) != this->meshes[0].end() );
	// '[0]' means meshes of the same dimension as the cell
	assert ( seg->tip_p == this );
	seg->tip_p = nullptr;
	std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
		it = this->meshes[0].find ( msh );
	Cell::field_to_meshes & field = it->second;
	assert ( field.counter_pos == 1 );
	assert ( field.counter_neg == 0 );
	// field.where is meaningless
	this->meshes[0].erase ( it );
	seg->deep_connections ( this, Mesh::action_remove );                  }


void Cell::Negative::Vertex::add_to ( Mesh::Core * msh ) // virtual from Cell::Core

{ assert ( msh );
	Cell::Positive::Segment * seg = ( Cell::Positive::Segment * ) msh;
	// ver.meshes[0] contains segments disguised as zero-dimensional meshes
	assert ( this->reverse_p );
	Cell::Positive::Vertex * pos_ver = (Cell::Positive::Vertex*) this->reverse_p;
	assert ( pos_ver->meshes.size() > 0 );
	// assert that 'this' vertex does not belong yet to the mesh 'msh'
	assert ( pos_ver->meshes[0].find(msh) == pos_ver->meshes[0].end() );
	// '[0]' means meshes of the same dimension as the cell
	assert ( seg->base_p == nullptr );
	seg->base_p = this;
	// pos_ver->meshes[0][msh] = Cell::field_to_meshes { 0, 1 };
	// the third component 'where' is irrelevant here
	pos_ver->meshes[0].emplace ( std::piecewise_construct,
		std::forward_as_tuple(msh), std::forward_as_tuple(0,1) );
	seg->deep_connections ( pos_ver, Mesh::action_add_rev );                       }


void Cell::Negative::Vertex::remove_from ( Mesh::Core * msh ) // virtual from Cell::Core

{ assert ( msh );
	Cell::Positive::Segment * seg = ( Cell::Positive::Segment * ) msh;
	// ver.meshes[0] contains segments disguised as zero-dimensional meshes
	assert ( this->reverse_p );
	Cell::Positive::Vertex * pos_ver = (Cell::Positive::Vertex*) this->reverse_p;
	assert ( pos_ver->meshes.size() > 0 );
	// assert that 'this' vertex belongs to the mesh 'msh'
	assert ( pos_ver->meshes[0].find(msh) != pos_ver->meshes[0].end() );
	// '[0]' means meshes of the same dimension as the cell
	assert ( seg->base_p );
	seg->base_p = nullptr;
	std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
		it = pos_ver->meshes[0].find ( msh );
	Cell::field_to_meshes & field = it->second;
	assert ( field.counter_pos == 0 );
	assert ( field.counter_neg == 1 );
	// field.where is meaningless
	pos_ver->meshes[0].erase ( it );
	seg->deep_connections ( pos_ver, Mesh::action_remove_rev );                     }


void Cell::Positive::Segment::add_to ( Mesh::Core * mmsh ) // virtual from Cell::Core

{	assert ( mmsh );
	Mesh::OneDim::Positive * msh = (Mesh::OneDim::Positive*) mmsh;
	assert ( this->meshes.size() > 0 );
	// assert that 'this' segment does not belong yet to the mesh 'msh'
	assert ( this->meshes[0].find(msh) == this->meshes[0].end() );
	// '[0]' means meshes of the same dimension as the cell

	assert ( this->base_p );
	assert ( this->tip_p );
	// assert that the segment fits into the chain of segments
	assert ( msh->cells.size() > 1 );

	// msh->first_ver == nullptr  means no check, no ordering
	// msh->first_ver == Cell::ghost  means closed chain
	// see paragraph 8.14 in the manual
	msh->first_ver = nullptr;
	
	msh->deep_connections ( this, this, Mesh::action_add );
	
	assert ( this->base_p->cell_behind_within.find(msh) ==
	         this->base_p->cell_behind_within.end()         );
	this->base_p->cell_behind_within[msh] = this;
	assert ( this->tip_p->cell_behind_within.find(msh) ==
	         this->tip_p->cell_behind_within.end()         );
	this->tip_p->cell_behind_within[msh] = this;

} // end of Cell::Positive::Segment::add_to


void Cell::Positive::Segment::remove_from ( Mesh::Core * mmsh ) // virtual from Cell::Core

{	assert ( mmsh );
	Mesh::OneDim::Positive * msh = (Mesh::OneDim::Positive*) mmsh;
	assert ( this->meshes.size() > 0 );
	// assert that 'this' segment belongs to the mesh 'msh'
	assert ( this->meshes[0].find(msh) != this->meshes[0].end() );
	// '[0]' means meshes of the same dimension as the cell

	assert ( this->base_p );
	assert ( this->tip_p );
//	assert ( this->tip_p->reverse_p );
//	assert ( this->base_p->reverse_p );
	assert ( msh->cells.size() > 1 );

	// msh->first_ver == nullptr  means no check, no ordering
	// msh->first_ver == Cell::ghost  means closed chain
	// see paragraph 8.14 in the manual
	msh->first_ver = nullptr;

	msh->deep_connections ( this, this, Mesh::action_remove );
	
	assert ( this->base_p->cell_behind_within.find(msh) !=
	         this->base_p->cell_behind_within.end()         );
	assert ( this->base_p->cell_behind_within[msh] == this );
	this->base_p->cell_behind_within.erase(msh);
	assert ( this->tip_p->cell_behind_within.find(msh) !=
	         this->tip_p->cell_behind_within.end()         );
	assert ( this->tip_p->cell_behind_within[msh] == this );
	this->tip_p->cell_behind_within.erase(msh);

} // end of Cell::Positive::Segment::remove_from


void Cell::Negative::Segment::add_to ( Mesh::Core * mmsh ) // virtual from Cell::Core

{	assert ( mmsh );
	Mesh::OneDim::Positive * msh = (Mesh::OneDim::Positive*) mmsh;
	Cell::Positive::Segment * pos_seg = ( Cell::Positive::Segment * ) this->reverse_p;
	assert ( pos_seg );
	assert ( pos_seg->meshes.size() > 0 );
	// assert that 'this' segment does not belong yet to the mesh 'msh'
	assert ( pos_seg->meshes[0].find(msh) == pos_seg->meshes[0].end() );
	// '[0]' means meshes of the same dimension as the cell

	assert ( pos_seg->base_p );
	assert ( pos_seg->tip_p );
	assert ( pos_seg->base_p->reverse_p );
	assert ( pos_seg->tip_p->reverse_p );

	// assert that the segment fits into the chain of segments
	assert ( msh->cells.size() > 1 );

	// msh->first_ver == nullptr  means no check, no ordering
	// msh->first_ver == Cell::ghost  means closed chain
	// see paragraph 8.14 in the manual
	msh->first_ver = nullptr;
	
	msh->deep_connections ( pos_seg, this, Mesh::action_add_rev );

	assert ( pos_seg->base_p->reverse_p->cell_behind_within.find(msh) ==
	         pos_seg->base_p->reverse_p->cell_behind_within.end()         );
	pos_seg->base_p->reverse_p->cell_behind_within[msh] = this;
	assert ( pos_seg->tip_p->reverse_p->cell_behind_within.find(msh) ==
	         pos_seg->tip_p->reverse_p->cell_behind_within.end()         );
	pos_seg->tip_p->reverse_p->cell_behind_within[msh] = this;

} // end of Cell::Negative::Segment::add_to


void Cell::Negative::Segment::remove_from ( Mesh::Core * mmsh ) // virtual from Cell::Core

{	assert ( mmsh );
	Mesh::OneDim::Positive * msh = (Mesh::OneDim::Positive*) mmsh;
	assert ( this->reverse_p );
	Cell::Positive::Segment * pos_seg = ( Cell::Positive::Segment * ) this->reverse_p;
	assert ( pos_seg->meshes.size() > 0 );
	// assert that 'this' segment belongs to the mesh 'msh'
	assert ( pos_seg->meshes[0].find(msh) != pos_seg->meshes[0].end() );
	// '[0]' means meshes of the same dimension as the cell

	assert ( pos_seg->base_p );
	assert ( pos_seg->tip_p );
	assert ( pos_seg->base_p->reverse_p );
	assert ( pos_seg->tip_p->reverse_p );
	// assert that the segment fits into the chain of segments
	assert ( msh->cells.size() > 1 );

	// msh->first_ver == nullptr  means no check, no ordering
	// msh->first_ver == Cell::ghost  means closed chain
	// see paragraph 8.14 in the manual
	msh->first_ver = nullptr;
	
	msh->deep_connections ( pos_seg, this, Mesh::action_remove_rev );
	
	assert ( pos_seg->base_p->reverse_p->cell_behind_within.find(msh) !=
	         pos_seg->base_p->reverse_p->cell_behind_within.end()         );
	assert ( pos_seg->base_p->reverse_p->cell_behind_within[msh] == this );
	pos_seg->base_p->reverse_p->cell_behind_within.erase(msh);
	assert ( pos_seg->tip_p->reverse_p->cell_behind_within.find(msh) !=
	         pos_seg->tip_p->reverse_p->cell_behind_within.end()         );
	assert ( pos_seg->tip_p->reverse_p->cell_behind_within[msh] == this );
	pos_seg->tip_p->reverse_p->cell_behind_within.erase(msh);

} // end of Cell::Negative::Segment::remove_from


void Cell::Positive::add_to ( Mesh::Core * msh ) // virtual from Cell::Core

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

{	assert ( msh );
	assert ( msh->get_dim_plus_one() == this->get_dim() + 1 );
	assert ( this->meshes.size() > 0 );
	// assert that 'this' cell does not belong yet to the mesh 'msh'
	assert ( this->meshes[0].find(msh) == this->meshes[0].end() );
	// '[0]' means meshes of the same dimension as the cell
	msh->deep_connections ( this, this, Mesh::action_add );
	
	Mesh::Core * bdry = this->boundary_p;  assert ( bdry );
	assert ( bdry->get_dim_plus_one() == this->get_dim() );
	std::list<Cell::Core*> & bcd = bdry->cells.back(); // last elem, cells[d-1]
	std::list<Cell::Core*>::iterator it = bcd.begin(), it_e = bcd.end();
	for ( ; it != it_e; it++ )
	{	Cell::Core * face = *it;
		assert ( face->cell_behind_within.find(msh) == face->cell_behind_within.end() );
		face->cell_behind_within[msh] = this;  	                                          }

} // end of Cell::Positive::add_to


void Cell::Positive::remove_from ( Mesh::Core * msh ) // virtual from Cell::Core

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'cut_from_bdry_of'

{	assert ( msh );
	assert ( msh->get_dim_plus_one() == this->get_dim() + 1 );
	assert ( this->meshes.size() > 0 );
	// assert that 'this' cell belongs to the mesh 'msh'
	assert ( this->meshes[0].find(msh) != this->meshes[0].end() );
	// '[0]' means meshes of the same dimension as the cell
	msh->deep_connections ( this, this, Mesh::action_remove );
	
	Mesh::Core * bdry = this->boundary_p;  assert ( bdry );
	assert ( bdry->get_dim_plus_one() == this->get_dim() );
	std::list<Cell::Core*> & bcd = bdry->cells.back(); // last elem, cells[d-1]
	std::list<Cell::Core*>::iterator it = bcd.begin(), it_e = bcd.end();
	for ( ; it != it_e; it++ )
	{	Cell::Core * face = *it;
		assert ( face->cell_behind_within.find(msh) != face->cell_behind_within.end() );
		assert ( face->cell_behind_within[msh] == this );
		face->cell_behind_within.erase(msh);                                              }

} // end of Cell::Positive::remove_from


void Cell::Negative::add_to ( Mesh::Core * msh ) // virtual from Cell::Core

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

{	assert ( msh );
	assert ( msh->get_dim_plus_one() == this->get_dim() + 1 );
	assert ( this->reverse_p );
	Cell::Positive * pos_cll = ( Cell::Positive* ) this->reverse_p;
	assert ( pos_cll );
	assert ( pos_cll->meshes.size() > 0 );
	// assert that 'this' cell does not belong yet to the mesh 'msh'
	assert ( pos_cll->meshes[0].find(msh) == pos_cll->meshes[0].end() );
	// '[0]' means meshes of the same dimension as the cell
	msh->deep_connections ( pos_cll, this, Mesh::action_add_rev );
	
	Mesh::Core * bdry = pos_cll->boundary_p;  assert ( bdry );
	assert ( bdry->get_dim_plus_one() == this->get_dim() );
	std::list<Cell::Core*> & bcd = bdry->cells.back(); // last elem, cells[d-1]
	std::list<Cell::Core*>::iterator it = bcd.begin(), it_e = bcd.end();
	for ( ; it != it_e; it++ )
	{	Cell::Core * face = *it;
		Cell::Core * rev_face = face->reverse_p;
		assert ( rev_face );
		assert ( rev_face->cell_behind_within.find(msh) ==
		         rev_face->cell_behind_within.end() );
		rev_face->cell_behind_within[msh] = this;             }

} // end of Cell::Negative::add_to


void Cell::Negative::remove_from ( Mesh::Core * msh ) // virtual from Cell::Core

// remove 'this' cell from the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'cut_from_bdry_of'

{	assert ( msh );
	assert ( msh->get_dim_plus_one() == this->get_dim() + 1 );
	assert ( this->reverse_p );
	Cell::Positive * pos_cll = ( Cell::Positive* ) this->reverse_p;
	assert ( pos_cll );
	assert ( pos_cll->meshes.size() > 0 );
	// assert that 'this' cell belongs to the mesh 'msh'
	assert ( pos_cll->meshes[0].find(msh) != pos_cll->meshes[0].end() );
	// '[0]' means meshes of the same dimension as the cell
	msh->deep_connections ( pos_cll, this, Mesh::action_remove_rev );
	
	Mesh::Core * bdry = pos_cll->boundary_p;  assert ( bdry );
	assert ( bdry->get_dim_plus_one() == this->get_dim() );
	std::list<Cell::Core*> & bcd = bdry->cells.back(); // last elem, cells[d-1]
	std::list<Cell::Core*>::iterator it = bcd.begin(), it_e = bcd.end();
	for ( ; it != it_e; it++ )
	{	Cell::Core * face = *it;
		Cell::Core * rev_face = face->reverse_p;
		assert ( rev_face );
		assert ( rev_face->cell_behind_within.find(msh) !=
		         rev_face->cell_behind_within.end() );
		assert ( rev_face->cell_behind_within[msh] == this );
		rev_face->cell_behind_within.erase(msh);              }

} // end of Cell::Negative::remove_from

//-----------------------------------------------------------------------------//


void Mesh::action_add ( Cell::Core * cll, Cell::Core * o_cll,
                        Mesh::Core * msh, short int cp, short int cn )
	
// This function "makes the link" between a cell and a mesh.
// This "link" asserts that the cell belongs to the mesh.
// If the dimensions are equal, then there will be only one instance of the cell
// in the mesh, so either cp==1 and cn==0 or cp==0 and cn==1, depending on the
// orientations.
// If 'cll' has lower dimension, it may appear several times in the mesh 'msh',
// with both orientations. The counter 'cp' asserts how many times 'cll' appears
// in the mesh with direct orientation (the same as the mesh), while 'cn' says
// the same thing for the reverse orientation (opposite to the mesh).
// This function is used only in add_to, passed to Mesh::Core::deep_connections.

{	assert ( cll );  assert ( o_cll );  assert ( msh );
	assert ( cll->is_positive() ); // assert ( msh->is_positive() );
	Cell::Core::Positive * cll_pos = (Cell::Core::Positive*) cll;
	size_t cll_dim = cll_pos->get_dim(),
	       msh_dim_p1 = msh->get_dim_plus_one();
	assert ( msh_dim_p1 > cll_dim );
	size_t dif_dim = msh_dim_p1 - cll_dim - 1;
	assert ( cll_pos->meshes.size() > dif_dim );
/////////////////////////////////////////////////////////////////////////////////////
	// inspired in item 24 of the book : Scott Meyers, Effective STL           //
	typedef std::map<Mesh::Core*,Cell::field_to_meshes> maptype;               //
	maptype & cmd = cll_pos->meshes[dif_dim];                                  //
	maptype::iterator lb = cmd.lower_bound(msh);                               //
	if ( ( lb == cmd.end() ) or ( cmd.key_comp()(msh,lb->first) ) )            //
	{ std::list <Cell::Core*> & mcd = msh->cells[cll_dim];                     //
	  mcd.push_front ( o_cll );                                                //
	  cmd.emplace_hint ( lb, std::piecewise_construct,                         //
	                     std::forward_as_tuple(msh),                           //
	                     std::forward_as_tuple(cp,cn,mcd.begin()) );    }      //
	else                                                                       //
	{ lb->second.counter_pos += cp;                                            //
	  lb->second.counter_neg += cn;   }                                        //
////////// code below is conceptually equivalent to the above //////////////////////
//	if (cll->meshes[dif_dim].find(msh)==cll->meshes[dif_dim].end())     //
//	{ msh->cells[cll->dim].push_front(o_cll);                           //
//	  Cell::field_to_meshes field;                                      //
//	  field.counter_pos = cp;                                           //
//	  field.counter_neg = cn;                                           //
//	  field.where = msh->cells[cll->dim]->begin();                      //
//	  cll->meshes[dif_dim][msh] = field;             }                  //
//	else                                                                //
//	{ cll->meshes[dif_dim][msh].counter_pos += cp;                      //
//	  cll->meshes[dif_dim][msh].counter_neg += cn;   }                  //
/////////////////////////////////////////////////////////////////////////////

} // end of Mesh::action_add


void Mesh::action_add_rev
( Cell::Core * cll, Cell::Core * o_cll,
  Mesh::Core * msh, short int cp, short int cn )
// we just switch the two counters
{	Mesh::action_add ( cll, o_cll, msh, cn, cp );  }


void Mesh::action_remove ( Cell::Core * cll, Cell::Core * o_cll,
                           Mesh::Core * msh, short int cp, short int cn )

// used only in remove_from, passed to Mesh::Core::deep_connections

{	assert ( cll );  assert ( o_cll );  assert ( msh );
	assert ( cll->is_positive() ); // assert ( msh->is_positive() );
	Cell::Core::Positive * cll_pos = (Cell::Core::Positive*) cll;
	size_t cll_dim = cll_pos->get_dim(),
	       msh_dim_p1 = msh->get_dim_plus_one();
	assert ( msh_dim_p1 > cll_dim );
	size_t dif_dim = msh_dim_p1 - cll_dim - 1;
	assert ( cll_pos->meshes.size() > dif_dim );
	typedef std::map <Mesh::Core*, Cell::field_to_meshes> maptype;
	maptype & cmd = cll_pos->meshes[dif_dim];
	maptype::iterator cmdm = cmd.find(msh);
	assert ( cmdm != cmd.end() );
	short int c_p = cmdm->second.counter_pos -= cp;
	short int c_n = cmdm->second.counter_neg -= cn;
	assert ( ( c_p >= 0 ) and ( c_n >= 0 ) );
	if ( ( c_p == 0 ) and ( c_n == 0 ) )
	{	std::list<Cell::Core*>::iterator w = cmdm->second.where;
		msh->cells[cll_dim].erase(w);
		cmd.erase(cmdm);                                                 }
} // end of Mesh::action_remove


void Mesh::action_remove_rev
( Cell::Core * cll, Cell::Core * o_cll,
  Mesh::Core * msh, short int cp, short int cn )
// we just switch the two counters
{	Mesh::action_remove ( cll, o_cll, msh, cn, cp );    }

//-----------------------------------------------------------------------------//


void Cell::Positive::Segment::deep_connections
( Cell::Core::Positive * cll,
	void (*action) ( Cell::Core*, Cell::Core*, Mesh::Core*, short int, short int ) )

// make or destroy connections when adding or removing a cell,
// according to the 'action' argument
// this version (belonging to the class Positive::Segment)
// is used to add a vertex to a segment

{	assert ( cll );
	assert ( cll->get_dim() == 0 );

	// this->base and this->tip have been set in 'add_to' before calling 'deep_connections'
	// also, cll->meshes[0] have been updated
	
	// we build a list of meshes "above" 'this' segment
	// we then use this list to create or destroy all connections

	struct triplet_mesh
	{	Mesh::Core * obj;
		short int counter_pos;
		short int counter_neg;
		triplet_mesh ( Mesh::Core * o, short int cp, short int cn )
			: obj (o), counter_pos (cp), counter_neg (cn)  { }                   };
	std::forward_list < triplet_mesh > all_meshes;
	// we now loop over all meshes above 'this' segment (of all dimensions)
	for ( size_t dif_dim = 0; dif_dim < this->meshes.size(); dif_dim++ )
	{	std::map<Mesh::Core*,Cell::field_to_meshes> & cemd = this->meshes[dif_dim];
		std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
			map_iter = cemd.begin(), cemd_end = cemd.end();
		// we now loop over all meshes of given dimension 'd'
		for ( ; map_iter != cemd_end; ++map_iter)
		{	Cell::field_to_meshes & mis = map_iter->second;
			assert ( map_iter->first->get_dim_plus_one() == dif_dim + 2 );
			all_meshes.emplace_front ( map_iter->first, mis.counter_pos, mis.counter_neg );
		}   }  // end of for, end of for

	std::forward_list <triplet_mesh>::iterator
		meshes_iter = all_meshes.begin(),  meshes_end = all_meshes.end();
	// we loop over "superior" dimensions	
	for (; meshes_iter != meshes_end; ++meshes_iter)
	{	triplet_mesh & higher_mesh_tr = *meshes_iter;
		short int mesh_counter_pos = higher_mesh_tr.counter_pos;
		short int mesh_counter_neg = higher_mesh_tr.counter_neg;
		action ( cll, cll, higher_mesh_tr.obj, mesh_counter_pos, mesh_counter_neg  );  }

} // end of Cell::Positive::Segment::deep_connections


void Mesh::Core::deep_connections
( Cell::Core::Positive * cll,  Cell::Core * o_cll,
	void (*action) ( Cell::Core *, Cell::Core *,
	                 Mesh::Core *, short int, short int ) )

// make or destroy connections when adding or removing a cell,
// according to the 'action' argument

{	assert ( this->get_dim_plus_one() > 1 );
	// the case cll.get_dim() == 0 is dealt with separately,
	// see Cell::Positive::Segment::deep_connections
	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );

	// We build two lists, a list of cells "below" 'cll'
	// (that is, belonging to the boundary of 'cll')
	// and a list of meshes "above" 'this' mesh
	// (that is, meshes to which 'this->cell_enclosed' belongs).
	// We then use these two lists to create or destroy all connections.

	struct triplet_cell
	{	Cell::Core::Positive * obj;
		size_t dim;
		short int counter_pos;
		short int counter_neg;
		triplet_cell ( Cell::Core::Positive * o, size_t d, short int cp, short int cn )
			: obj {o}, dim {d}, counter_pos {cp}, counter_neg {cn}  { }                        };
	std::forward_list < triplet_cell > all_cells;
	size_t cll_dim = cll->get_dim();
	all_cells.emplace_front ( cll, cll_dim, 1, 0 );
	if ( cll_dim == 1 ) // 'cll' is a segment
	{	Cell::Positive::Segment * seg = (Cell::Positive::Segment*) cll;
		assert ( seg->base_p );  assert ( seg->tip_p );
		assert ( not ( seg->base_p->is_positive() ) );  assert ( seg->tip_p->is_positive() );
		assert ( seg->base_p->reverse_p );
		Cell::Positive::Vertex * pos_base = (Cell::Positive::Vertex*) seg->base_p->reverse_p;
		all_cells.emplace_front ( pos_base, 0, 0, 1 );
		all_cells.emplace_front ( seg->tip_p, 0, 1, 0 );                                        }
	else
	{	assert ( cll_dim > 1 );
		Cell::Positive * cll_loc = (Cell::Positive*) cll;
		Mesh::Core * cell_bdry = cll_loc->boundary_p;
		assert ( cell_bdry );
		size_t cell_bdry_dim = cell_bdry->get_dim_plus_one() - 1; // debug mode
		assert ( cell_bdry_dim == cll_dim - 1 );
		// we loop over all cells of cell_bdry (of all dimensions)	
		for ( size_t d = 0; d < cell_bdry_dim; d++ )
		{ // we now loop over all cells of given dimension 'd' (not maximum)
			std::list<Cell::Core*> & cbcd = cell_bdry->cells[d];
			std::list<Cell::Core*>::iterator
				list_iter = cbcd.begin(), cbcd_end = cbcd.end();
			for ( ; list_iter != cbcd_end; ++list_iter )
			{	Cell::Core::Positive* lower_cell = (Cell::Core::Positive*) *list_iter;
				assert ( lower_cell->get_dim() == d );
				assert ( cell_bdry_dim >= d );
				size_t dif_dim = cell_bdry_dim - d;
				assert ( lower_cell->meshes.size() > dif_dim );
				// code below should be equivalent to & fcb = lower_cell->meshes[dif_dim][cell_bdry]
				// but apparantly it's not, probably because operator[]
				// is designed to create a new element if it doesn't find the key
				std::map<Mesh::Core*,Cell::field_to_meshes >::iterator
					it = lower_cell->meshes[dif_dim].find(cell_bdry);
				assert ( it != lower_cell->meshes[dif_dim].end() );
				Cell::field_to_meshes & fcb = it->second ;
				all_cells.emplace_front ( lower_cell, d, fcb.counter_pos, fcb.counter_neg );  }  }
		// we now loop over all cells of maximum dimension
		std::list<Cell::Core*> & cbcd = cell_bdry->cells[cell_bdry_dim];
		std::list<Cell::Core*>::iterator
			list_iter = cbcd.begin(), cbcd_end = cbcd.end();
		for ( ; list_iter != cbcd_end; ++list_iter )
		{	Cell::Core* lower_cell = *list_iter;
			assert ( lower_cell );
			assert ( lower_cell->get_dim() == cell_bdry_dim );
			Cell::Core::Positive* pos_lower_cll = lower_cell->get_positive();
			assert ( pos_lower_cll );
			assert ( pos_lower_cll->meshes.size() > 0 );
			// code below should be equivalent to & fcb = pos_lower_cll->meshes[0][cell_bdry]
			// but apparantly it's not, probably because operator[] is designed to
			// create a new element if the key does not exist
			std::map<Mesh::Core*,Cell::field_to_meshes >::iterator
				it = pos_lower_cll->meshes[0].find(cell_bdry);
			// 'meshes[0]' means meshes of dimension equal to the dimension of the cell
			assert ( it != pos_lower_cll->meshes[0].end() );
			Cell::field_to_meshes & fcb = it->second ;
			all_cells.emplace_front ( pos_lower_cll, cell_bdry_dim,
		                            fcb.counter_pos, fcb.counter_neg );                 }
	} // end of if-else cll_dim
	
	struct triplet_mesh
	{	Mesh::Core * obj;
		size_t dim;
		short int counter_pos;
		short int counter_neg;
		triplet_mesh ( Mesh::Core * o, size_t d, short int cp, short int cn )
			: obj (o), dim (d), counter_pos (cp), counter_neg (cn)  { }                   };
	std::forward_list < triplet_mesh > all_meshes;
	all_meshes.emplace_front ( this, cll_dim, 1, 0 );
	if ( this->cell_enclosed )
	{	// we now loop over all meshes above 'this->cell_enclosed' (of all dimensions)
		for ( size_t dif_dim = 0; dif_dim < this->cell_enclosed->meshes.size(); dif_dim++ )
		{	std::map<Mesh::Core*,Cell::field_to_meshes> &
				cemd = this->cell_enclosed->meshes[dif_dim];
			std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
				map_iter = cemd.begin(), cemd_end = cemd.end();
			// we now loop over all meshes of given dimension
			for (; map_iter != cemd_end; ++map_iter)
			{	assert ( map_iter->first->get_dim_plus_one() ==
			           dif_dim + this->cell_enclosed->get_dim() + 1                );
				size_t current_mesh_dim = map_iter->first->get_dim_plus_one() - 1;
				Cell::field_to_meshes & mis = map_iter->second;
				all_meshes.emplace_front ( map_iter->first, current_mesh_dim,
			                             mis.counter_pos, mis.counter_neg );	            }  }         }

	// we loop over cells "below" 'cell'
	std::forward_list<triplet_cell>::iterator cells_iter, cells_end = all_cells.end();
	for ( cells_iter = all_cells.begin(); cells_iter != cells_end; ++cells_iter)
	{	std::forward_list <triplet_mesh>::iterator
			meshes_iter = all_meshes.begin(),  meshes_end = all_meshes.end();
		triplet_cell & lower_cell_tri = *cells_iter;
		short int cell_counter_pos = lower_cell_tri.counter_pos;
		short int cell_counter_neg = lower_cell_tri.counter_neg;
		// we loop over "superior" dimensions	
		for (; meshes_iter != meshes_end; ++meshes_iter)
		{	triplet_mesh & higher_mesh_tri = *meshes_iter;
			short int mesh_counter_pos = higher_mesh_tri.counter_pos;
			short int mesh_counter_neg = higher_mesh_tri.counter_neg;
			if ( lower_cell_tri.dim == higher_mesh_tri.dim )
			{	assert ( lower_cell_tri.obj == cll );
				action ( cll, o_cll, higher_mesh_tri.obj,
				         cell_counter_pos*mesh_counter_pos + cell_counter_neg*mesh_counter_neg,
				         cell_counter_pos*mesh_counter_neg + cell_counter_neg*mesh_counter_pos  );  }
			else
				action ( lower_cell_tri.obj, lower_cell_tri.obj, higher_mesh_tri.obj,
				         cell_counter_pos*mesh_counter_pos + cell_counter_neg*mesh_counter_neg,
				         cell_counter_pos*mesh_counter_neg + cell_counter_neg*mesh_counter_pos  );
		}  }

} // end of Mesh::Core::deep_connections

//-----------------------------------------------------------------------------//

