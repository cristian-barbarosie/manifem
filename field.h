
// field.h 2020.01.03

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

#ifndef MANIFEM_FIELD_H
#define MANIFEM_FIELD_H

#include "mesh.h"

namespace maniFEM {

class Field;

namespace tag {
	struct LivesOnPositiveCells { };  static const LivesOnPositiveCells lives_on_positive_cells;
	struct LivesOnPoints { };  static const LivesOnPoints lives_on_points;
	struct HasIndexInHeap { };  static const HasIndexInHeap has_index_in_heap;
}

	
class Field

// wrapper for fields

{	public :

	class Core;
	
	Field::Core * core;

	inline Field ( const tag::WhoseCoreIs &, Field::Core * c )
	:	core { c }
	{	assert ( c );  }

	inline Field ( const tag::LivesOnPoints & );
	inline Field ( const tag::LivesOnPoints &, const tag::HasSize &, size_t s );

	inline Field operator[] ( size_t );
	class TakenOnCell;
	inline Field::TakenOnCell operator() ( Cell );
	
	class Block;  class Scalar;
};


class Field::Core

// base class for several different fields

{	public :

	size_t lives_on_cells_of_dim;

	inline Core ( const tag::LivesOnPositiveCells &, const tag::OfDimension &, size_t d )
	:	lives_on_cells_of_dim { d }
	{	}

	virtual ~Core ( ) { };

	virtual size_t nb_of_components ( ) = 0;
	virtual Field::Scalar * component ( size_t ) = 0;

	inline Field::TakenOnCell on_cell ( Cell::Core * cll );

};


class Field::Scalar : public Field::Core
	
{	public :

	size_t index_in_heap;

	inline Scalar ( const tag::LivesOnPositiveCells &, const tag::OfDimension &, size_t d )
	:	Field::Core ( tag::lives_on_positive_cells, tag::of_dim, d )
	{	index_in_heap = Cell::double_heap_size_pos[d];
		Cell::double_heap_size_pos[d] ++;               }
	
	inline Scalar ( const tag::LivesOnPositiveCells &, const tag::OfDimension &,
		size_t d, const tag::HasIndexInHeap, size_t i )
	:	Field::Core ( tag::lives_on_positive_cells, tag::of_dim, d )
	{	index_in_heap = i;  }
	
	size_t nb_of_components ( );  // virtual from Field::Core
	Field::Scalar * component ( size_t );  // virtual from Field::Core

};


// we should have a Field::Vector as parent for Field::Block and Field::Aggregate

class Field::Block : public Field::Core
	
{	public :

	size_t min_index, max_index_p1;
	
	inline Block ( const tag::LivesOnPositiveCells &, const tag::OfDimension &, size_t d,
	                    const tag::HasSize &, size_t s                                         )
	:	Field::Core ( tag::lives_on_positive_cells, tag::of_dim, d )
	{	min_index = Cell::double_heap_size_pos[d];
		Cell::double_heap_size_pos[d] += s;
		max_index_p1 = Cell::double_heap_size_pos[d];  }
	
	size_t nb_of_components ( );  // virtual from Field::Core
	Field::Scalar * component ( size_t );  // virtual from Field::Core

};

	
inline Field::Field ( const tag::LivesOnPoints & )
:	core { new Field::Scalar ( tag::lives_on_positive_cells, tag::of_dim, 0 ) }
{	assert ( core );  }

	
inline Field::Field ( const tag::LivesOnPoints &, const tag::HasSize &, size_t s )
:	core { new Field::Block ( tag::lives_on_positive_cells, tag::of_dim, 0,
                                  tag::has_size, s                               ) }
{	assert ( core );  }

	
// class Field::TakenOnCell holds a temporary object,
// the result of Field::operator(), or of Field::Core::on_cell
// it is versatile in the sense that the Field may have one components or several
// if the Field f has one component, we may use the Field::TakenOnCell object like in
// double x = f(cll)  or  f(cll) = 2.0
// if the Field f has several components, we may use the Field::TakenOnCell object like in
// vector<double> vec = f(cll)  or  f(cll) = vec

// Trebuie sa ne gandim si la chestii de tipul += *=
// si chiar  x = y = F(P)  sau  F(P) = F(Q) = 2  si chiar  F(P) = F(P)

class Field::TakenOnCell	

{	public :

	Field::Core * f;
	Cell::Core * cll;

	TakenOnCell ( const Field::TakenOnCell & ) = delete;
	TakenOnCell ( const Field::TakenOnCell && ) = delete;
	Field::TakenOnCell operator= ( const Field::TakenOnCell & ) = delete;
	Field::TakenOnCell operator= ( const Field::TakenOnCell && ) = delete;
	
	inline operator double ()
	// can be used like in  double x = f(cll)  or  cout << f(cll)
	{	return this->reference();  }
		
	inline double & operator= ( const double & x )
	// can be used like in  f(cll) = 2.0
	{	double & y = this->reference();
		y = x;  return y;                                   	}

	inline double & reference ( )
	// can be used like in  f(cll) = 2.0
	{	Field::Scalar * f_scalar = dynamic_cast<Field::Scalar*> ( f );
		assert ( f_scalar );
		return cll->double_heap[f_scalar->index_in_heap];               }

	inline operator std::vector<double> ()
	// can be used like in  vector<double> vec = f(cll)
	{	Field::Block * f_block = dynamic_cast<Field::Block*> ( f );
		assert ( f_block );
		return std::vector<double> { & cll->double_heap[f_block->min_index],
		                             & cll->double_heap[f_block->max_index_p1] };  }
	
	inline std::vector<double> operator= ( const std::vector<double> & x )
	// can be used like in  f(cll) = vec
	{	Field::Block * f_block = dynamic_cast<Field::Block*> ( f );
		assert ( f_block );
		size_t i_min    = f_block->min_index,
		       i_max_p1 = f_block->max_index_p1;
		for ( size_t i = i_min; i < i_max_p1; i++ )
			cll->double_heap[i] = x[i];  // there should be a faster solution
		return x;                                                       }

};


inline Field::TakenOnCell Field::Core::on_cell ( Cell::Core * cll )
{	assert ( this->lives_on_cells_of_dim == cll->get_dim() );
	return Field::TakenOnCell { this, cll };                 }

inline Field Field::operator[] ( size_t i )
{	return Field ( tag::whose_core_is, this->core->component(i) );  }
	
inline Field::TakenOnCell Field::operator() ( Cell cll )
{	return this->core->on_cell(cll.core);  }

	
} // namespace maniFEM

#endif  // ifndef MANIFEM_FIELD_H
