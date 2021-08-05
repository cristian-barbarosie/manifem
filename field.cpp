
// field.cpp 2020.08.05

//   This file is part of maniFEM, a C++ library for meshes and finite elements on manifolds.

//   Copyright 2019, 2020 Cristian Barbarosie cristian.barbarosie@gmail.com
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

#include "field.h"

using namespace maniFEM;

size_t Field::ShortInt::Scalar::nb_of_components ( )  // virtual from Field::Core
{	return 1;  }

size_t Field::Double::Scalar::nb_of_components ( )  // virtual from Field::Core
{	return 1;  }

size_t Field::SizeT::nb_of_components ( )  // virtual from Field::Core
{	return 1;  }


size_t Field::Double::Block::nb_of_components ( )  // virtual from Field::Core
{	return tag::Util::assert_diff ( this->max_index_p1, this->min_index );  }
// tag::Util::assert_diff  provides a safe way to substract two size_t numbers


Field::ShortInt::Scalar * Field::ShortInt::Scalar::component ( size_t i )
// virtual from Field::Core
{	assert ( i == 0 );
	return this;        }

Field::Double::Scalar * Field::Double::Scalar::component ( size_t i )  // virtual from Field::Core
{	assert ( i == 0 );
	return this;        }


Field::ShortInt::Scalar * Field::ShortInt::Block::component ( size_t i )
// virtual from Field::Core
// we build a new Field::Scalar
// it would be nice to check whether a field for index i exists already ...
{	assert ( i < this->nb_of_components() );
	return new Field::ShortInt::Scalar ( tag::lives_on_positive_cells, tag::of_dim,
	      this->lives_on_cells_of_dim, tag::has_index_in_heap, this->min_index + i );  }

Field::Double::Scalar * Field::Double::Block::component ( size_t i )  // virtual from Field::Core
// we build a new Field::Scalar
// it would be nice to check whether a field for index i exists already ...
{	assert ( i < this->nb_of_components() );
	return new Field::Double::Scalar ( tag::lives_on_positive_cells, tag::of_dim,
	    this->lives_on_cells_of_dim, tag::has_index_in_heap, this->min_index + i );  }

//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------


size_t & Cell::Numbering::Field::operator() ( const Cell p )  // virtual from Cell::Numbering
{	return this->field.on_cell ( p.core );  }

void Cell::Numbering::Field::set_and_increment ( Cell::Core * cll, void * num )  // static
{	Cell::Numbering::Field * num_f = static_cast < Cell::Numbering::Field* > ( num );
	cll->size_t_heap[num_f->field.index_in_heap] = num_f->counter;
	num_f->counter++;                                                                 }
		
