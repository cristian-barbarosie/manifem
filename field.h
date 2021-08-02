
//   field.h  2021.08.02

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
#include "iterator.h"

namespace maniFEM {

class Field;

namespace tag {
	struct LivesOnPositiveCells { };  static const LivesOnPositiveCells lives_on_positive_cells;
	struct LivesOnPoints { };  static const LivesOnPoints lives_on_points;
	struct HasIndexInHeap { };  static const HasIndexInHeap has_index_in_heap;
}

//---------------------------------------------------------------------------------------

	
// we have three types of Fields
// a Field::ShortInt can be used for labeling cells in different regions
//   or for holding the jump of a segment
// a Field::SizeT can be used for enumerating vertices
// a Field::Double can be used for storing coordinates of vertices or values of solutions

struct Field  {  class Core;  class ShortInt;  class SizeT;  class Double;  };


class Field::Core

// core field holding different types of values
// specialized in classes Field::ShortInt, Field::SizeT, Field::Double

{	public :

	size_t lives_on_cells_of_dim;

	inline Core ( const tag::LivesOnPositiveCells &, const tag::OfDimension &, size_t d )
	:	lives_on_cells_of_dim { d }
	{	}

	virtual ~Core ( ) { };

	virtual size_t nb_of_components ( ) = 0;  

};


//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------

// there will be just a few Field objects in a program,
// and they will usually be destroyed only at the end of the program,
// so we do not use here the mechanism of inheriting from tag::Util::Wrapper and tag::Util::Core
	
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------


// a Field::ShortInt can be used for labeling cells in different regions
//   or for holding the jump of a segment
	
class Field::ShortInt

// wrapper for fields holding short int values

{	public :

	class Core;
	
	Field::ShortInt::Core * core;

	inline ShortInt ( const tag::WhoseCoreIs &, Field::ShortInt::Core * c )
	:	core { c }
	{	assert ( c );  }

	inline ShortInt ( const tag::LivesOnPoints & );
	inline ShortInt ( const tag::LivesOnPoints &, const tag::HasSize &, size_t s );

	inline Field::ShortInt operator[] ( size_t );
	class TakenOnCell;
	inline Field::ShortInt::TakenOnCell operator() ( Cell );
	
	class Core;  class Block;  class Scalar;
};


class Field::ShortInt::Core : public Field::Core

// base class for several different fields

{	public :

	// attribute  lives_on_cells_of_dim  inherited from Field::Core

	inline Core ( const tag::LivesOnPositiveCells &, const tag::OfDimension &, size_t d )
	:	Field::Core ( tag::lives_on_positive_cells, tag::of_dim, d )
	{	}

	virtual ~Core ( ) { };

	// size_t nb_of_components  stays pure virtual from Field::Core
	
	virtual Field::ShortInt::Scalar * component ( size_t ) = 0;

	inline Field::ShortInt::TakenOnCell on_cell ( Cell::Core * cll );

};

	
class Field::ShortInt::Scalar : public Field::ShortInt::Core
	
{	public :

	// attribute  lives_on_cells_of_dim  inherited from Field::Core

	size_t index_in_heap;

	inline Scalar ( const tag::LivesOnPositiveCells &, const tag::OfDimension &, size_t d )
	:	Field::ShortInt::Core ( tag::lives_on_positive_cells, tag::of_dim, d )
	{	index_in_heap = Cell::double_heap_size_pos[d];
		Cell::double_heap_size_pos[d] ++;               }
	
	inline Scalar ( const tag::LivesOnPositiveCells &, const tag::OfDimension &,
		size_t d, const tag::HasIndexInHeap, size_t i )
	:	Field::ShortInt::Core ( tag::lives_on_positive_cells, tag::of_dim, d )
	{	index_in_heap = i;  }
	
	size_t nb_of_components ( );  // virtual from Field::Core

	Field::ShortInt::Scalar * component ( size_t );  // virtual from Field::ShortInt::Core

	// Field::ShortInt::TakenOnCell on_cell ( Cell::Core * cll )
	//   defined inline by Field::ShortInt::Core

};


class Field::ShortInt::Block : public Field::ShortInt::Core
	
{	public :

	// attribute  lives_on_cells_of_dim  inherited from Field::Core

	size_t min_index, max_index_p1;
	
	inline Block ( const tag::LivesOnPositiveCells &, const tag::OfDimension &, size_t d,
	                    const tag::HasSize &, size_t s                                         )
	:	Field::ShortInt::Core ( tag::lives_on_positive_cells, tag::of_dim, d )
	{	min_index = Cell::short_int_heap_size_pos[d];
		Cell::short_int_heap_size_pos[d] += s;
		max_index_p1 = Cell::short_int_heap_size_pos[d];  }
	
	size_t nb_of_components ( );  // virtual from Field::Core

	Field::ShortInt::Scalar * component ( size_t );  // virtual from Field::Core

	// Field::ShortInt::TakenOnCell on_cell ( Cell::Core * cll )
	//   defined inline by Field::ShortInt::Core

};

//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------

	
// a Field::SizeT can be used for enumerating vertices

class Field::SizeT

// wrapper for fields holding size_t values

{	public :

	class Core;
	
	Field::SizeT::Core * core;

	inline SizeT ( const tag::WhoseCoreIs &, Field::SizeT::Core * c )
	:	core { c }
	{	assert ( c );  }

	inline SizeT ( const tag::LivesOnPoints & );
	inline SizeT ( const tag::LivesOnPoints &, const tag::HasSize &, size_t s );

	inline Field::SizeT operator[] ( size_t );
	class TakenOnCell;
	inline Field::SizeT::TakenOnCell operator() ( Cell );
	
	class Core;  class Block;  class Scalar;
};


class Field::SizeT::Core : public Field::Core

// base class for several different fields

{	public :

	// attribute  lives_on_cells_of_dim  inherited from Field::Core

	inline Core ( const tag::LivesOnPositiveCells &, const tag::OfDimension &, size_t d )
	:	Field::Core ( tag::lives_on_positive_cells, tag::of_dim, d )
	{	}

	virtual ~Core ( ) { };

	// size_t nb_of_components  stays pure virtual from Field::Core

	virtual Field::SizeT::Scalar * component ( size_t ) = 0;

	inline Field::SizeT::TakenOnCell on_cell ( Cell::Core * cll );

	// Field::SizeT::TakenOnCell on_cell ( Cell::Core * cll )
	//   defined inline by Field::SizeT::Core

};


class Field::SizeT::Scalar : public Field::SizeT::Core
	
{	public :

	// attribute  lives_on_cells_of_dim  inherited from Field::Core

	size_t index_in_heap;

	inline Scalar ( const tag::LivesOnPositiveCells &, const tag::OfDimension &, size_t d )
	:	Field::SizeT::Core ( tag::lives_on_positive_cells, tag::of_dim, d )
	{	index_in_heap = Cell::double_heap_size_pos[d];
		Cell::double_heap_size_pos[d] ++;               }
	
	inline Scalar ( const tag::LivesOnPositiveCells &, const tag::OfDimension &,
		size_t d, const tag::HasIndexInHeap, size_t i )
	:	Field::SizeT::Core ( tag::lives_on_positive_cells, tag::of_dim, d )
	{	index_in_heap = i;  }
	
	size_t nb_of_components ( );  // virtual from Field::Core

	Field::SizeT::Scalar * component ( size_t );  // virtual from Field::SizeT::Core

};


class Field::SizeT::Block : public Field::SizeT::Core
	
{	public :

	// attribute  lives_on_cells_of_dim  inherited from Field::Core

	size_t min_index, max_index_p1;
	
	inline Block ( const tag::LivesOnPositiveCells &, const tag::OfDimension &, size_t d,
	                    const tag::HasSize &, size_t s                                         )
	:	Field::SizeT::Core ( tag::lives_on_positive_cells, tag::of_dim, d )
	{	min_index = Cell::size_t_heap_size_pos[d];
		Cell::size_t_heap_size_pos[d] += s;
		max_index_p1 = Cell::size_t_heap_size_pos[d];  }
	
	size_t nb_of_components ( );  // virtual from Field::Core

	Field::SizeT::Scalar * component ( size_t );  // virtual from Field::SizeT::Core

	// Field::SizeT::TakenOnCell on_cell ( Cell::Core * cll )
	//   defined inline by Field::SizeT::Core

};

//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------


// a Field::Double can be used for storing coordinates of vertices or values of solutions
	
class Field::Double

// wrapper for fields

{	public :

	class Core;
	
	Field::Double::Core * core;

	inline Double ( const tag::WhoseCoreIs &, Field::Double::Core * c )
	:	core { c }
	{	assert ( c );  }

	inline Double ( const tag::LivesOnPoints & );
	inline Double ( const tag::LivesOnPoints &, const tag::HasSize &, size_t s );

	inline Field::Double operator[] ( size_t );
	class TakenOnCell;
	inline Field::Double::TakenOnCell operator() ( Cell );
	
	class Core;  class Block;  class Scalar;
};


class Field::Double::Core : public Field::Core

// base class for several different fields

{	public :

	// attribute  lives_on_cells_of_dim  inherited from Field::Core

	inline Core ( const tag::LivesOnPositiveCells &, const tag::OfDimension &, size_t d )
	:	Field::Core ( tag::lives_on_positive_cells, tag::of_dim, d )
	{	}

	virtual ~Core ( ) { };

	// size_t nb_of_components  stays pure virtual from Field::Core

	virtual Field::Double::Scalar * component ( size_t ) = 0;

	inline Field::Double::TakenOnCell on_cell ( Cell::Core * cll );

};


class Field::Double::Scalar : public Field::Double::Core
	
{	public :

	// attribute  lives_on_cells_of_dim  inherited from Field::Core

	size_t index_in_heap;

	inline Scalar ( const tag::LivesOnPositiveCells &, const tag::OfDimension &, size_t d )
	:	Field::Double::Core ( tag::lives_on_positive_cells, tag::of_dim, d )
	{	index_in_heap = Cell::double_heap_size_pos[d];
		Cell::double_heap_size_pos[d] ++;               }
	
	inline Scalar ( const tag::LivesOnPositiveCells &, const tag::OfDimension &,
		size_t d, const tag::HasIndexInHeap, size_t i )
	:	Field::Double::Core ( tag::lives_on_positive_cells, tag::of_dim, d )
	{	index_in_heap = i;  }
	
	size_t nb_of_components ( );  // virtual from Field::Double::Core

	Field::Double::Scalar * component ( size_t );  // virtual from Field::Double::Core

	// Field::Double::TakenOnCell on_cell ( Cell::Core * cll )
	//   defined inline by Field::Double::Core

};


class Field::Double::Block : public Field::Double::Core
	
{	public :

	// attribute  lives_on_cells_of_dim  inherited from Field::Core

	size_t min_index, max_index_p1;
	
	inline Block ( const tag::LivesOnPositiveCells &, const tag::OfDimension &, size_t d,
	                    const tag::HasSize &, size_t s                                   )
	:	Field::Double::Core ( tag::lives_on_positive_cells, tag::of_dim, d )
	{	min_index = Cell::double_heap_size_pos[d];
		Cell::double_heap_size_pos[d] += s;
		max_index_p1 = Cell::double_heap_size_pos[d];  }
	
	size_t nb_of_components ( );  // virtual from Field::Double::Core
	
	Field::Double::Scalar * component ( size_t );  // virtual from Field::Double::Core

	// Field::Double::TakenOnCell on_cell ( Cell::Core * cll )
	//   defined inline by Field::Double::Core

};

	
inline Field::Double::Double ( const tag::LivesOnPoints & )
:	core { new Field::Double::Scalar ( tag::lives_on_positive_cells, tag::of_dim, 0 ) }
{	assert ( core );  }

	
inline Field::Double::Double ( const tag::LivesOnPoints &, const tag::HasSize &, size_t s )
:	core { new Field::Double::Block ( tag::lives_on_positive_cells, tag::of_dim, 0,
                                  tag::has_size, s                               ) }
{	assert ( core );  }

	
// class Field::Double::TakenOnCell holds a temporary object,
// the result of Field::Double::operator(), or of Field::Double::Core::on_cell
// it is versatile in the sense that the Field may have one components or several
// if the Field f has one component, we may use the Field::Double::TakenOnCell object like in
// double x = f(cll)  or  f(cll) = 2.0
// if the Field f has several components, we may use the Field::Double::TakenOnCell object like in
// vector<double> vec = f(cll)  or  f(cll) = vec

// Trebuie sa ne gandim si la chestii de tipul += *=
// si chiar  x = y = F(P)  sau  F(P) = F(Q) = 2  si chiar  F(P) = F(P)

class Field::Double::TakenOnCell	

{	public :

	Field::Double::Core * f;
	Cell::Core * cll;

	TakenOnCell ( const Field::Double::TakenOnCell & ) = delete;
	TakenOnCell ( const Field::Double::TakenOnCell && ) = delete;
	Field::Double::TakenOnCell operator= ( const Field::Double::TakenOnCell & ) = delete;
	Field::Double::TakenOnCell operator= ( const Field::Double::TakenOnCell && ) = delete;
	
	inline operator double ()
	// can be used like in  double x = f(cll)  or  cout << f(cll)
	{	return this->reference();  }
		
	inline double & operator= ( const double & x )
	// can be used like in  f(cll) = 2.0
	{	double & y = this->reference();
		y = x;  return y;               }

	inline double & reference ( )
	// can be used like in  f(cll) = 2.0
	{	Field::Double::Scalar * f_scalar = tag::Util::assert_cast
			< Field::Double::Core*, Field::Double::Scalar* > ( f );
		return cll->double_heap[f_scalar->index_in_heap];        }

	inline operator std::vector<double> ()
	// can be used like in  vector<double> vec = f(cll)
	{	Field::Double::Block * f_block = tag::Util::assert_cast
			< Field::Double::Core*, Field::Double::Block* > ( f );
		return std::vector<double> { & cll->double_heap[f_block->min_index],
		                             & cll->double_heap[f_block->max_index_p1] };  }
	
	inline std::vector<double> operator= ( const std::vector<double> & x )
	// can be used like in  f(cll) = vec
	{	Field::Double::Block * f_block = tag::Util::assert_cast
			< Field::Double::Core*, Field::Double::Block* > ( f );
		size_t i_min    = f_block->min_index,
		       i_max_p1 = f_block->max_index_p1;
		for ( size_t i = i_min; i < i_max_p1; i++ )
			cll->double_heap[i] = x[i];  // there should be a faster solution
		return x;                                               }

};


inline Field::Double::TakenOnCell Field::Double::Core::on_cell ( Cell::Core * cll )
{	assert ( this->lives_on_cells_of_dim == cll->get_dim() );
	return Field::Double::TakenOnCell { this, cll };          }

inline Field::Double Field::Double::operator[] ( size_t i )
{	return Field::Double ( tag::whose_core_is, this->core->component(i) );  }
	
inline Field::Double::TakenOnCell Field::Double::operator() ( Cell cll )
{	return this->core->on_cell(cll.core);  }

	
} // namespace maniFEM

#endif  // ifndef MANIFEM_FIELD_H
