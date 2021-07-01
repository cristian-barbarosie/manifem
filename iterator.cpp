
// maniFEM iterator.cpp 2019.10.30

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

#include "iterator.h"

using namespace maniFEM;


void CellIterator::Over::CellsOfMesh::reset ( Cell::Core * cll )
// virtual from CoreCellIterator
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "This should never happen. Something is very wrong." << std::endl;
	exit ( 1 );                                                                                   }

void CellIterator::Over::CellsOfMesh::reset ( ) // virtual from CoreCellIterator
{	this->iter = this->list_p->begin();  }

void CellIterator::Over::CellsOfMesh::advance ( ) // virtual from CoreCellIterator
{	this->iter++;  }
	
bool CellIterator::Over::CellsOfMesh::in_range ( )	 // virtual from CoreCellIterator
{	return this->iter != this->list_p->end();  }


Cell::Core * CellIterator::Over::CellsOfPosMesh::deref ( )
// virtual from CoreCellIterator
{	return * ( this->iter );  }

Cell::Core * CellIterator::Over::CellsOfPosMesh::Positive::deref ( )
// virtual from CoreCellIterator
{	return ( * this->iter )->get_positive();  }

Cell::Core * CellIterator::Over::CellsOfNegMesh::deref ( )
// virtual from CoreCellIterator
{	return * ( this->iter );  }

Cell::Core * CellIterator::Over::CellsOfNegMesh::MaxDim::deref ( )
// virtual from CoreCellIterator
{	Cell::Core * cll_p = * this->iter;
	cll_p = cll_p->reverse_p;
	assert ( cll_p );
	return cll_p;                              }
	
Cell::Core * CellIterator::Over::CellsOfNegMesh::Positive::deref ( )
// virtual from CoreCellIterator
{	return ( * this->iter )->get_positive();  }

	
void CellIterator::Over::CellsOfMesh::order_mesh ( )  // virtual from CoreCellIterator
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "This should never happen. Something is very wrong." << std::endl;
	exit ( 1 );                                                                                   }
	
//-------------------------------------------------------------------------------//

	
void CellIterator::Over::SegsOfPosChain::reset ( Cell::Core * cll )
// virtual from CoreCellIterator,
// through IterOver::OrderedCells, IterOver::OrderedCells::Open, CellIterator::Over::SegsOfChain
{	assert ( cll->get_dim() == 1 );
	// assert cll is a segment in this->mesh !
	this->order_mesh();
	assert ( this->mesh->first_ver );
	assert ( this->mesh->first_ver != Cell::ghost );
	this->current_cell = cll;        }

void CellIterator::Over::SegsOfPosChain::Reverse::reset ( Cell::Core * cll )
// virtual from CoreCellIterator,
// through IterOver::OrderedCells, IterOver::OrderedCells::Open, CellIterator::Over::SegsOfChain
{	assert ( cll->get_dim() == 1 );
	// assert cll is a segment in this->mesh !
	this->order_mesh();
	assert ( this->mesh->first_ver );
	assert ( this->mesh->first_ver != Cell::ghost );
	this->current_cell = cll;        }

void CellIterator::Over::SegsOfNegChain::reset ( Cell::Core * cll )
// virtual, overrides definition by IterOver::SegsOfPosChain::Reverse
{	assert ( cll->get_dim() == 1 );
	// assert cll->reverse_p is a segment in this->mesh !
	this->order_mesh();
	assert ( this->mesh->first_ver );
	assert ( this->mesh->first_ver != Cell::ghost );
	assert ( cll->reverse_p );
	this->current_cell = cll->reverse_p;   }

void CellIterator::Over::SegsOfNegChain::Reverse::reset ( Cell::Core * cll )
// virtual, overrides definition by IterOver::SegsOfPosChain
{	assert ( cll->get_dim() == 1 );
	// assert cll->reverse_p is a segment in this->mesh !
	this->order_mesh();
	assert ( this->mesh->first_ver );
	assert ( this->mesh->first_ver != Cell::ghost );
	assert ( cll->reverse_p );
	this->current_cell = cll->reverse_p;   }

void CellIterator::Over::SegsOfLoop::reset ( Cell::Core * cll )
{	assert ( cll->get_dim() == 1 );
	// assert cll is a segment in this->mesh !
	this->order_mesh();
	assert ( this->mesh->first_ver == Cell::ghost );
	this->current_cell = this->first_cell = cll;
	this->fresh = true;                           }

void CellIterator::Over::SegsOfNegLoop::reset ( Cell::Core * cll )
// virtual, overrides definition by IterOver::SegsOfLoop
{	assert ( cll->get_dim() == 1 );
	// assert cll->reverse_p is a segment in this->mesh !
	this->order_mesh();
	assert ( this->mesh->first_ver == Cell::ghost );
	this->current_cell = this->first_cell = cll->reverse_p;
	this->fresh = true;                                      }

void CellIterator::Over::SegsOfNegLoop::Reverse::reset ( Cell::Core * cll )
// virtual, overrides definition by IterOver::SegsOfLoop
{	assert ( cll->get_dim() == 1 );
	// assert cll->reverse_p is a segment in this->mesh !
	this->order_mesh();
	assert ( this->mesh->first_ver == Cell::ghost );
	this->current_cell = this->first_cell = cll->reverse_p;
	this->fresh = true;                                      }

void CellIterator::Over::VerOfChain::reset ( Cell::Core * cll )
{	assert ( cll->get_dim() == 0 );
	assert ( cll->is_positive() );
	// assert cll is a vertex in this->mesh !
	this->order_mesh();
	assert ( this->mesh->first_ver );
	assert ( this->mesh->first_ver != Cell::ghost );
	this->current_cell = cll;         }

void CellIterator::Over::VerOfLoop::reset ( Cell::Core * cll )
{	assert ( cll->get_dim() == 0 );
	assert ( cll->is_positive() );
	// assert cll is a vertex in this->mesh !
	this->order_mesh();
	assert ( this->mesh->first_ver == Cell::ghost );
	this->current_cell = this->first_cell = cll;
	this->fresh = true;                           }

	
void CellIterator::Over::SegsOfPosChain::reset ( )
// virtual from CoreCellIterator, through IterOver::CellsOfOneDimMesh, 
//    IterOver::OrderedCells, IterOver::OrderedCells::Open, CellIterator::Over::SegsOfChain
{	this->current_cell = this->mesh->first_segment();  }

void CellIterator::Over::SegsOfPosChain::Reverse::reset ( )
// virtual from CoreCellIterator, through IterOver::CellsOfOneDimMesh, 
//    IterOver::OrderedCells, IterOver::OrderedCells::Open, CellIterator::Over::SegsOfChain
{	this->current_cell = this->mesh->last_segment();  }

void CellIterator::Over::SegsOfNegChain::reset ( )
// virtual, overrides definition by IterOver::SegsOfPosChain::Reverse
// the two functions are identical, but the compiler forces me to override both 'reset's
{	this->current_cell = this->mesh->last_segment();  }

void CellIterator::Over::SegsOfNegChain::Reverse::reset ( )
// virtual, overrides definition by IterOver::SegsOfPosChain
// the two functions are identical, but the compiler forces me to override both 'reset's
{	this->current_cell = this->mesh->first_segment();  }

void CellIterator::Over::VerOfPosChain::reset ( )
{	this->current_cell = this->mesh->first_vertex();  }

void CellIterator::Over::VerOfPosChain::Reverse::reset ( )
{	this->current_cell = this->mesh->last_vertex();  }
	
void CellIterator::Over::SegsOfLoop::reset ( )
{	this->order_mesh();
	assert ( this->mesh->first_ver == Cell::ghost );
	std::list<Cell::Core*>::iterator it0 = this->mesh->cells[1].begin();  // will change !
	if ( it0 == this->mesh->cells[1].end() )  // empty mesh
	{	this->current_cell = this->first_cell = nullptr;
		this->fresh = false;  return;                    }
	Cell::Core * seg = *it0;  assert ( seg );
	this->current_cell = this->first_cell = seg;
	this->fresh = true;                                                    }

void CellIterator::Over::SegsOfNegLoop::reset ( )
// virtual, overrides definition by IterOver::SegsOfLoop
// the two functions are identical, but the compiler forces me to override both 'reset's
{	this->order_mesh();
	assert ( this->mesh->first_ver == Cell::ghost );
	std::list<Cell::Core*>::iterator it0 = this->mesh->cells[1].begin();  // will change !
	if ( it0 == this->mesh->cells[1].end() )  // empty mesh
	{	this->current_cell = this->first_cell = nullptr;
		this->fresh = false;  return;                    }
	Cell::Core * seg = *it0;  assert ( seg );
	this->current_cell = this->first_cell = seg;
	this->fresh = true;                                                    }

void CellIterator::Over::SegsOfNegLoop::Reverse::reset ( )
// virtual, overrides definition by IterOver::SegsOfLoop
// the two functions are identical, but the compiler forces me to override both 'reset's
{	this->order_mesh();
	assert ( this->mesh->first_ver == Cell::ghost );
	std::list<Cell::Core*>::iterator it0 = this->mesh->cells[1].begin();  // will change !
	if ( it0 == this->mesh->cells[1].end() )  // empty mesh
	{	this->current_cell = this->first_cell = nullptr;
		this->fresh = false;  return;                    }
	Cell::Core * seg = *it0;  assert ( seg );
	this->current_cell = this->first_cell = seg;
	this->fresh = true;                                                    }

void CellIterator::Over::VerOfLoop::reset ( )
{	this->order_mesh();
	assert ( this->mesh->first_ver == Cell::ghost );
	std::list<Cell::Core*>::iterator it0 = this->mesh->cells[1].begin();  // will change !
	if ( it0 == this->mesh->cells[1].end() )  // empty mesh
	{	this->current_cell = this->first_cell = nullptr;
		this->fresh = false;
		return;                                            }
	Cell::Core * seg = *it0;
	this->current_cell = this->first_cell = seg->base()->reverse_p;
	assert ( this->current_cell );
	this->fresh = true;                                                   }
	

Cell::Core * CellIterator::Over::OrderedCells::deref ( )
// virtual from CoreCellIterator
{	assert ( this->current_cell );
	return this->current_cell;     }

Cell::Core * CellIterator::Over::SegsOfNegChain::deref ( )
// virtual, overrides definition by IterOver::OrderedCells
{	assert ( this->current_cell );
	assert ( this->current_cell->reverse_p );
	return this->current_cell->reverse_p;     }
	
Cell::Core * CellIterator::Over::SegsOfNegChain::Reverse::deref ( )
// virtual, overrides definition by IterOver::OrderedCells
{	assert ( this->current_cell );
	assert ( this->current_cell->reverse_p );
	return this->current_cell->reverse_p;     }
	
Cell::Core * CellIterator::Over::SegsOfNegLoop::deref ( )
// virtual, overrides definition by IterOver::OrderedCells
{	assert ( this->current_cell );
	assert ( this->current_cell->reverse_p );
	return this->current_cell->reverse_p;     }
	
Cell::Core * CellIterator::Over::SegsOfNegLoop::Reverse::deref ( )
// virtual, overrides definition by IterOver::OrderedCells
{	assert ( this->current_cell );
	assert ( this->current_cell->reverse_p );
	return this->current_cell->reverse_p;     }
	
Cell::Core * CellIterator::Over::SegsOfPosChain::Positive::deref ( )
// virtual, overrides definition by IterOver::OrderedCells
{	assert ( this->current_cell );
	assert ( this->current_cell->get_positive() );
	return this->current_cell->get_positive();     }
	
Cell::Core * CellIterator::Over::SegsOfPosChain::Reverse::Positive::deref ( )
// virtual, overrides definition by IterOver::OrderedCells
{	assert ( this->current_cell );
	assert ( this->current_cell->get_positive() );
	return this->current_cell->get_positive();     }
	
Cell::Core * CellIterator::Over::SegsOfNegChain::Positive::deref ( )
// virtual, overrides definition by IterOver::OrderedCells and by IterOver::SegsOfNegChain
{	assert ( this->current_cell );
	assert ( this->current_cell->get_positive() );
	return this->current_cell->get_positive();     }
	
Cell::Core * CellIterator::Over::SegsOfNegChain::Reverse::Positive::deref ( )
// virtual, overrides definition by IterOver::OrderedCells and by SegsOfNegChain::Reverse
{	assert ( this->current_cell );
	assert ( this->current_cell->get_positive() );
	return this->current_cell->get_positive();     }
	
Cell::Core * CellIterator::Over::SegsOfPosLoop::Positive::deref ( )
// virtual, overrides definition by IterOver::OrderedCells and by IterOver::SegsOfPosLoop
{	assert ( this->current_cell );
	assert ( this->current_cell->get_positive() );
	return this->current_cell->get_positive();     }
	
Cell::Core * CellIterator::Over::SegsOfPosLoop::Reverse::Positive::deref ( )
// virtual, overrides definition by IterOver::OrderedCells and by SegsOfPosLoop::Reverse
{	assert ( this->current_cell );
	assert ( this->current_cell->get_positive() );
	return this->current_cell->get_positive();     }
	
Cell::Core * CellIterator::Over::SegsOfNegLoop::Positive::deref ( )
// virtual, overrides definition by IterOver::OrderedCells and by IterOver::SegsOfNegLoop
{	assert ( this->current_cell );
	assert ( this->current_cell->get_positive() );
	return this->current_cell->get_positive();     }
	
Cell::Core * CellIterator::Over::SegsOfNegLoop::Reverse::Positive::deref ( )
// virtual, overrides definition by IterOver::OrderedCells and by SegsOfNegLoop::Reverse
{	assert ( this->current_cell );
	assert ( this->current_cell->get_positive() );
	return this->current_cell->get_positive();     }
	

bool CellIterator::Over::OrderedCells::Open::in_range ( )
// virtual from CoreCellIterator, through IterOver::OrderedCells
{	return this->current_cell != nullptr;  }
	
bool CellIterator::Over::OrderedCells::Loop::in_range ( )
// virtual from CoreCellIterator, through IterOver::OrderedCells
{	if ( this->fresh ) return this->current_cell != nullptr;
	// this allows for empty meshes
	return this->current_cell != this->first_cell;           }

	
void CellIterator::Over::SegsOfPosChain::advance ( )
// virtual from CoreCellIterator, through IterOver::CellsOfOneDimMesh, 
//    IterOver::OrderedCells, IterOver::OrderedCells::Open, CellIterator::Over::SegsOfChain

{ Cell::Core * ver = this->current_cell->tip();
	assert ( ver );  assert ( ver->is_positive() );
	assert ( ver );
	// we have no Mesh, only a core this->mesh, so we cannot use	
	// this->current_cell = this->mesh->cell_in_front_of ( ver, tag::may_not_exist )
	Cell::Core * neg_ver = ver->reverse_p;
	if ( neg_ver == nullptr )
	{	assert ( ver == this->mesh->last_ver );
		this->current_cell = nullptr;  return;   }
	assert ( not neg_ver->is_positive() );
	std::map<Mesh::Core*,Cell::Core*>::const_iterator
		it = neg_ver->cell_behind_within.find ( this->mesh );
	if ( it == neg_ver->cell_behind_within.end() )
	{	this->current_cell = nullptr;  return;  }
	assert ( it->second );
	this->current_cell = it->second;                                 }
	

void CellIterator::Over::SegsOfPosChain::Reverse::advance ( )
// virtual from CoreCellIterator, through IterOver::CellsOfOneDimMesh, 
//    IterOver::OrderedCells, IterOver::OrderedCells::Open, CellIterator::Over::SegsOfChain

{ Cell::Core * neg_ver = this->current_cell->base();
	assert ( neg_ver );  assert ( not neg_ver->is_positive() );
	// we have no Mesh, only a core this->mesh, so we cannot use	
	// this->current_cell = this->mesh->cell_in_front_of ( ver, tag::may_not_exist )
	Cell::Core * ver = neg_ver->reverse_p;
	assert ( ver );  assert ( ver->is_positive() );
	std::map<Mesh::Core*,Cell::Core*>::const_iterator
		it = ver->cell_behind_within.find ( this->mesh );
	if ( it == ver->cell_behind_within.end() )
	{	this->current_cell = nullptr;  return;  }
	assert ( it->second );
	this->current_cell = it->second;                                 }
	

void CellIterator::Over::SegsOfPosLoop::advance ( )  // virtual from CoreCellIterator,
  // through IterOver::CellsOfOneDimMesh, IterOver::OrderedCells, IterOver::OrderedCells::Open

{ Cell::Core * ver = this->current_cell->tip();
	assert ( ver );  assert ( ver->is_positive() );
	// we have no Mesh, only a core this->mesh, so we cannot use	
	// this->current_cell = this->mesh->cell_in_front_of ( ver, tag::may_not_exist )
	Cell::Core * neg_ver = ver->reverse_p;
	assert ( neg_ver );  assert ( not neg_ver->is_positive() );
	std::map<Mesh::Core*,Cell::Core*>::const_iterator
		it = neg_ver->cell_behind_within.find ( this->mesh );
	assert ( it != neg_ver->cell_behind_within.end() );
	assert ( it->second );
	this->current_cell = it->second;
	this->fresh = false;                                            }

	
void CellIterator::Over::SegsOfPosLoop::Reverse::advance ( )  // virtual from CoreCellIterator,
  // through IterOver::CellsOfOneDimMesh, IterOver::OrderedCells, IterOver::OrderedCells::Open

{ Cell::Core * neg_ver = this->current_cell->base();
	assert ( neg_ver );  assert ( not neg_ver->is_positive() );
	Cell::Core * ver = neg_ver->reverse_p;
	assert ( ver );  assert ( ver->is_positive() );
	// we have no Mesh, only a core this->mesh, so we cannot use	
	// this->current_cell = this->mesh->cell_behind ( ver, tag::may_not_exist )
	std::map<Mesh::Core*,Cell::Core*>::const_iterator
		it = ver->cell_behind_within.find ( this->mesh );
	assert ( it != ver->cell_behind_within.end() );
	assert ( it->second );
	this->current_cell = it->second;
	this->fresh = false;                                            }

	
void CellIterator::Over::VerOfPosChain::advance ( )  // virtual from CoreCellIterator,
//   through IterOver::CellsOfOneDimMesh, IterOver::OrderedCells,
//   IterOver::OrderedCells::Open, IterOver::VerOfChain

{	// we have no Mesh, only a core this->mesh, so we cannot use	
	// seg = this->mesh->cell_in_front_of ( this->current_cell, tag::may_not_exist )
	Cell::Core * neg_ver = this->current_cell->reverse_p;
	if ( neg_ver == nullptr )
	{	this->current_cell = nullptr;  return;  }
	std::map<Mesh::Core*,Cell::Core*>::const_iterator
		it = neg_ver->cell_behind_within.find ( this->mesh );
	if ( it == neg_ver->cell_behind_within.end() )
	{	this->current_cell = nullptr;  return;  }
	Cell::Core * seg = it->second;
	if ( seg )
	{	this->current_cell = seg->tip();
		assert ( this->current_cell );   }
  else this->current_cell = nullptr;                                            }

	
void CellIterator::Over::VerOfPosChain::Reverse::advance ( )  // virtual from CoreCellIterator,
//   through IterOver::CellsOfOneDimMesh, IterOver::OrderedCells,
//   IterOver::OrderedCells::Open, IterOver::VerOfChain

{	// we have no Mesh, only a core this->mesh, so we cannot use	
	// seg = this->mesh->cell_behind ( this->current_cell, tag::may_not_exist )
	Cell::Core * ver = this->current_cell;  assert ( ver );
	std::map<Mesh::Core*,Cell::Core*>::const_iterator
		it = ver->cell_behind_within.find ( this->mesh );
	if ( it == ver->cell_behind_within.end() )
	{	this->current_cell = nullptr;  return;  }
	Cell::Core * seg = it->second;
	assert ( seg );
	this->current_cell = seg->base()->reverse_p;
	assert ( this->current_cell );                                         }

	
void CellIterator::Over::VerOfPosLoop::advance ( )  // virtual from CoreCellIterator,
// through IterOver::CellsOfOneDimMesh, IterOver::OrderedCells, IterOver::OrderedCells::Open

{	// we have no Mesh, only a core this->mesh, so we cannot use	
	// seg = this->mesh->cell_in_front_of ( this->current_cell, tag::may_not_exist )
	Cell::Core * neg_ver = this->current_cell->reverse_p;
	assert ( neg_ver );
	std::map<Mesh::Core*,Cell::Core*>::const_iterator
		it = neg_ver->cell_behind_within.find ( this->mesh );
	assert ( it != neg_ver->cell_behind_within.end() );
	Cell::Core * seg = it->second;  assert ( seg );
	assert ( seg->base()->reverse_p == this->current_cell );
	this->current_cell = seg->tip();
	assert ( this->current_cell );                           
	this->fresh = false;                                                          }


void CellIterator::Over::VerOfPosLoop::Reverse::advance ( )  // virtual from CoreCellIterator,
// through IterOver::CellsOfOneDimMesh, OrderedCells, OrderedCells::Open, VerOfPosLoop

{	// we have no Mesh, only a core this->mesh, so we cannot use	
	// seg = this->mesh->cell_behind ( this->current_cell, tag::may_not_exist )
	std::map<Mesh::Core*,Cell::Core*>::const_iterator
		it = this->current_cell->cell_behind_within.find ( this->mesh );
	assert ( it != this->current_cell->cell_behind_within.end() );
	Cell::Core * seg = it->second;  assert ( seg );
	assert ( seg->tip() == this->current_cell );
	this->current_cell = seg->base()->reverse_p;
	assert ( this->current_cell );                           
	this->fresh = false;                                                          }
	

void CellIterator::Over::CellsOfOneDimMesh::order_mesh ( )  // virtual from CoreCellIterator
{	this->mesh->order();
	assert ( this->mesh->first_ver );  }

