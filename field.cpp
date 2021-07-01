
// maniFEM field.cpp 2019.10.30

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

#include "field.h"

using namespace maniFEM;

size_t Field::Scalar::nb_of_components ( )  // virtual from Field::Core
{	return 1;  }


size_t Field::Block::nb_of_components ( )  // virtual from Field::Corex
{	return Mesh::diff ( this->max_index_p1, this->min_index );  }
// Mesh::diff  provides a safe way to substract two size_t numbers

Field::Scalar * Field::Scalar::component ( size_t i )  // virtual from Field::Core
{	assert ( i == 0 );
	return this;        }

Field::Scalar * Field::Block::component ( size_t i )  // virtual from Field::Core
// we build a new Field::Scalar
// it would be nice to check whether a field for index i exists already ...
{	assert ( i < this->nb_of_components() );
	return new Field::Scalar ( tag::lives_on_positive_cells, tag::of_dim,
	  this->lives_on_cells_of_dim, tag::has_index_in_heap, this->min_index + i );  }





