
// iterator.h 2021.11.26

//   This file is part of maniFEM, a C++ library for meshes and finite elements on manifolds.

//   Copyright 2019, 2020, 2021 Cristian Barbarosie cristian.barbarosie@gmail.com
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

#ifndef MANIFEM_ITERATOR_H
#define MANIFEM_ITERATOR_H

#include "mesh.h"

namespace maniFEM {

namespace tag {
	struct OverTwoVerticesOfSeg { };
	static const OverTwoVerticesOfSeg over_two_vertices_of_seg;
	struct OverVerticesOf { };  static const OverVerticesOf over_vertices_of;
	struct OverSegmentsOf { };  static const OverSegmentsOf over_segments_of;
	struct OverCellsOf { };  static const OverCellsOf over_cells_of;
	struct FuzzyPosMesh { };  static const FuzzyPosMesh fuzzy_pos_mesh;
	struct Center { };  static const Center center;
	struct Inner { };  static const Inner inner;
}

//-----------------------------------------------------------------------------------------

// class Mesh::Iterator defined in mesh.h

//-----------------------------------------------------------------------------------------

class Mesh::Iterator::Core

// iterates over all cells of a given mesh (cells of given dimension)
// for maximum dimension (equal to the dimension of the mesh) returns oriented cells
// which may be positive or negative
// for lower dimension, returns positive cells

// or, iterates over all cells above a given cell (containing a given cell)

{	public :

	inline Core ( ) { };
	
	virtual ~Core ( ) { };

	virtual void reset ( ) = 0;
	virtual void reset ( const tag::StartAt &, Cell::Core * cll ) = 0;

	virtual Cell deref ( ) = 0;
	virtual void advance ( ) = 0;
	virtual bool in_range ( ) = 0;
	
};  // end of class Mesh::Iterator::Core


inline Mesh::Iterator::CellIterator ( const Mesh::Iterator & it)
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "ManiFEM doesn't know (yet) how to copy construct an iterator. Sorry." << std::endl;
	exit ( 1 );                                                                                }
// define virtual self-replicator for core iterators !!
	

// inline Mesh::Iterator::CellIterator ( Mesh::Iterator && it)
// {	this->core = std::move ( it.core );  }
	

inline Mesh::Iterator & Mesh::Iterator::operator= ( const Mesh::Iterator & it )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "ManiFEM doesn't know (yet) how to copy an iterator. Sorry." << std::endl;
	exit ( 1 );                                                                                   }
// define virtual self-replicator for core iterators !!


// inline Mesh::Iterator & Mesh::Iterator::operator= ( Mesh::Iterator && it )
// {	this->core = std::move ( it.core );  return *this;  }
	

inline void Mesh::Iterator::reset ( )  {  this->core->reset();  }


inline void Mesh::Iterator::reset ( const tag::StartAt &, const Cell & cll )
{  this->core->reset( tag::start_at, cll.core );  }
	

inline Cell Mesh::Iterator::operator* ( )
{	return this->core->deref();  }


inline Mesh::Iterator & Mesh::Iterator::operator++ ( )  {  return this->advance();  }


inline Mesh::Iterator & Mesh::Iterator::operator++ ( int )  {  return this->advance();  }


inline Mesh::Iterator & Mesh::Iterator::advance ( )
{ this->core->advance();  return *this;  }


inline bool Mesh::Iterator::in_range ( )  {  return this->core->in_range();  }
		
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::TwoVerticesOfSeg : public Mesh::Iterator::Core

{	public :

	Cell::PositiveSegment * seg_p;  // should an iterator keep the segment alive ?
	short unsigned int passage;

	inline TwoVerticesOfSeg ( Cell::PositiveSegment * seg )
	:	Mesh::Iterator::Core (), seg_p { seg } { };
	
	void reset ( );  // virtual from Mesh::Iterator::Core

	void reset ( const tag::StartAt &, Cell::Core * cll );
	// virtual from Mesh::Iterator::Core, here execution forbidden

	// Cell deref ( )  stays pure virtual from Mesh::Iterator::Core
	void advance ( );  // virtual from Mesh::Iterator::Core
	bool in_range ( );  // virtual from Mesh::Iterator::Core

	struct NormalOrder  {  class AsTheyAre;  class ForcePositive;
		struct ReverseEachCell  {  class  AssumeCellsExist;  class BuildCells;  };    };
	struct ReverseOrder  {  class AsTheyAre;  class ForcePositive;
		struct ReverseEachCell  {  class  AssumeCellsExist;  class BuildCells;  };    };
	
};  // end of class Mesh::Iterator::Over::TwoVerticesOfSeg

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::TwoVerticesOfSeg::NormalOrder::AsTheyAre
: public Mesh::Iterator::Over::TwoVerticesOfSeg
// iterate over the two vertices, first base (negative) then tip (positive)

{	public :

	// inherited from Mesh::Iterator::Over::TwoVerticesOfSeg :
	// Cell::PositiveSegment * seg_p;
	// short unsigned int passage;

	inline AsTheyAre ( Cell::PositiveSegment * seg )
	:	Mesh::Iterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( )  virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg, execution forbidden

	Cell deref ( );  // virtual from Mesh::Iterator::Core
	// void advance ( )  virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg

};  // end of class Mesh::Iterator::Over::TwoVerticesOfSeg::NormalOrder::AsTheyAre

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::TwoVerticesOfSeg::NormalOrder::ForcePositive
: public Mesh::Iterator::Over::TwoVerticesOfSeg
// iterate over the two vertices, first base then tip (both positive)

{	public :

	// inherited from Mesh::Iterator::Over::TwoVerticesOfSeg :
	// Cell::PositiveSegment * seg_p;
	// short unsigned int passage;

	inline ForcePositive ( Cell::PositiveSegment * seg )
	:	Mesh::Iterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg, execution forbidden

	Cell deref ( );  // virtual from Mesh::Iterator::Core
	// void advance ( )  virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg

};  // end of class Mesh::Iterator::Over::TwoVerticesOfSeg::NormalOrder::ForcePositive

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::TwoVerticesOfSeg::ReverseOrder::AsTheyAre
: public Mesh::Iterator::Over::TwoVerticesOfSeg
// iterate over the two vertices, first tip (positive) then base (negative)

{	public :

	// inherited from Mesh::Iterator::Over::TwoVerticesOfSeg :
	// Cell::PositiveSegment * seg_p;
	// short unsigned int passage;
	
	inline AsTheyAre ( Cell::PositiveSegment * seg )
	:	Mesh::Iterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( )  virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg, execution forbidden

	Cell deref ( );  // virtual from Mesh::Iterator::Core
	// void advance ( )  virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg

};  // end of class Mesh::Iterator::Over::TwoVerticesOfSeg::ReverseOrder::AsTheyAre

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::TwoVerticesOfSeg::ReverseOrder::ForcePositive
: public Mesh::Iterator::Over::TwoVerticesOfSeg
// iterate over the two vertices, first tip then base (both positive)

{	public :

	// inherited from Mesh::Iterator::Over::TwoVerticesOfSeg :
	// Cell::PositiveSegment * seg_p;
	// short unsigned int passage;

	inline ForcePositive ( Cell::PositiveSegment * seg )
	:	Mesh::Iterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( )  virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg, execution forbidden

	Cell deref ( );  // virtual from Mesh::Iterator::Core
	// void advance ( )  virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg

};  // end of class Mesh::Iterator::Over::TwoVerticesOfSeg::ReverseOrder::ForcePositive

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::TwoVerticesOfSeg::NormalOrder::ReverseEachCell::AssumeCellsExist
: public Mesh::Iterator::Over::TwoVerticesOfSeg
// iterate over the two vertices, first base (positive) then tip (negative)

{	public :

	// inherited from Mesh::Iterator::Over::TwoVerticesOfSeg :
	// Cell::PositiveSegment * seg_p;
	// short unsigned int passage;

	inline AssumeCellsExist ( Cell::PositiveSegment * seg )
	:	Mesh::Iterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( )  virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg, execution forbidden

	Cell deref ( );  // virtual from Mesh::Iterator::Core
	// void advance ( )  virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg

};  // end of class Mesh::Iterator::Over::TwoVerticesOfSeg
    //         ::NormalOrder::ReverseEachCell::AssumeCellsExist

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::TwoVerticesOfSeg::NormalOrder::ReverseEachCell::BuildCells
: public Mesh::Iterator::Over::TwoVerticesOfSeg
// iterate over the two vertices, first base (positive) then tip (negative)

{	public :

	// inherited from Mesh::Iterator::Over::TwoVerticesOfSeg :
	// Cell::PositiveSegment * seg_p;
	// short unsigned int passage;

	inline BuildCells ( Cell::PositiveSegment * seg )
	:	Mesh::Iterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( )  virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg, execution forbidden

	Cell deref ( );  // virtual from Mesh::Iterator::Core
	// void advance ( )  virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg

};  // end of class Mesh::Iterator::Over::TwoVerticesOfSeg
    //         ::NormalOrder::ReverseEachCell::BuildCells

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::TwoVerticesOfSeg::ReverseOrder::ReverseEachCell::AssumeCellsExist
: public Mesh::Iterator::Over::TwoVerticesOfSeg
// iterate over the two vertices, first tip (negative) then base (positive)

{	public :

	// inherited from Mesh::Iterator::Over::TwoVerticesOfSeg :
	// Cell::PositiveSegment * seg_p;
	// short unsigned int passage;

	inline AssumeCellsExist ( Cell::PositiveSegment * seg )
	:	Mesh::Iterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( )  virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg, execution forbidden

	Cell deref ( );  // virtual from Mesh::Iterator::Core
	// void advance ( )  virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg

};  // end of class Mesh::Iterator::Over::TwoVerticesOfSeg
    //         ::ReverseOrder::ReverseEachCell::AssumeCellsExist

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::TwoVerticesOfSeg::ReverseOrder::ReverseEachCell::BuildCells
: public Mesh::Iterator::Over::TwoVerticesOfSeg
// iterate over the two vertices, first tip (negative) then base (positive)

{	public :

	// inherited from Mesh::Iterator::Over::TwoVerticesOfSeg :
	// Cell::PositiveSegment * seg_p;
	// short unsigned int passage;

	inline BuildCells ( Cell::PositiveSegment * seg )
	:	Mesh::Iterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( )  virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg, execution forbidden

	Cell deref ( );  // virtual from Mesh::Iterator::Core
	// void advance ( )  virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  virtual, defined by Mesh::Iterator::Over::TwoVerticesOfSeg

};  // end of class Mesh::Iterator::Over::TwoVerticesOfSeg
    //         ::ReverseOrder::ReverseEachCell::BuildCells

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::CellsOfConnectedOneDimMesh : public Mesh::Iterator::Core

{	public :

	Mesh::Connected::OneDim * msh;  // should an iterator keep the mesh alive ?
	Cell::Core * last_vertex;
	// positive vertex for iterator over vertices
	// positive vertex for iterator over segments, forward
	// negative vertex for iterator over segments, back

	inline CellsOfConnectedOneDimMesh ( Mesh::Connected::OneDim * m )
	:	Mesh::Iterator::Core (), msh { m } { }
	
	// void reset ( )  stays pure virtual from Mesh::Iterator::Core

	// void reset ( tag::StartAt , Cell::Core * cll )  stays pure virtual from Mesh::Iterator::Core

	// Cell deref ( )  stays pure virtual from Mesh::Iterator::Core
	// void advance ( )  stays pure virtual from Mesh::Iterator::Core
	// bool in_range ( )  stays pure virtual from Mesh::Iterator::Core

};  // end of class Mesh::Iterator::Over::CellsOfConnectedOneDimMesh

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh :
public Mesh::Iterator::Over::CellsOfConnectedOneDimMesh

{	public :

	// attributes inherited from Mesh::Iterator::Over::CellsOfConnectedOneDimMesh :
	// Mesh::Connected::OneDim * msh
	// Cell::Core * last_vertex (here positive vertex)

	Cell::Positive::Vertex * current_vertex;

	inline VerticesOfConnectedOneDimMesh ( Mesh::Connected::OneDim * m )
	:	Mesh::Iterator::Over::CellsOfConnectedOneDimMesh ( m ) { }
	
	// void reset ( )  stays pure virtual from Mesh::Iterator::Core

	// void reset ( tag::StartAt, Cell::Core * cll )  stays pure virtual from Mesh::Iterator::Core

	Cell deref ( );  // virtual from Mesh::Iterator::Core
	// void advance ( )  stays pure virtual from Mesh::Iterator::Core
	bool in_range ( );  // virtual from Mesh::Iterator::Core

	class NormalOrder;  class ReverseOrder;
	
};  // end of class Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder
: public Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (here positive vertex)
	// inherited from Mesh::Iterator::Over::CellsOfConnectedOneDimMesh

	// attribute Cell::Positive::Vertex * current_vertex
	// inherited from Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh
	
	inline NormalOrder ( Mesh::Connected::OneDim * m )
	:	Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh ( m ) { };
	
	void reset ( );  // virtual from Mesh::Iterator::Core

	void reset ( const tag::StartAt &, Cell::Core * cll );  // virtual from Mesh::Iterator::Core

	// Cell deref ( )  defined by Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh
	void advance ( );  // virtual from Mesh::Iterator::Core
	// bool in_range ( )  virtual, defined by Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh

};  // end of class Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh::ReverseOrder
: public Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (here positive vertex)
	// inherited from Mesh::Iterator::Over::CellsOfConnectedOneDimMesh

	// attribute Cell::Positive::Vertex * current_vertex
	// inherited from Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh
	
	inline ReverseOrder ( Mesh::Connected::OneDim * m )
	:	Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh ( m ) { };
	
	void reset ( );  // virtual from Mesh::Iterator::Core

	void reset ( const tag::StartAt &, Cell::Core * cll );  // virtual from Mesh::Iterator::Core

	// Cell deref ( )  defined by Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh
	void advance ( );  // virtual from Mesh::Iterator::Core
	// bool in_range ( )  virtual, defined by Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh

};  // end of class Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh::ReverseOrder

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh :
public Mesh::Iterator::Over::CellsOfConnectedOneDimMesh

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (may be positive or negative)
	// inherited from Mesh::Iterator::Over::CellsOfConnectedOneDimMesh

	Cell::Core * current_segment;
	
	inline SegmentsOfConnectedOneDimMesh ( Mesh::Connected::OneDim * m )
	:	Mesh::Iterator::Over::CellsOfConnectedOneDimMesh ( m ) { }
	
	// void reset ( )  stays pure virtual from Mesh::Iterator::Core

	// void reset ( tag::StartAt, Cell::Core * cll )  stays pure virtual from Mesh::Iterator::Core

	// Cell deref ( )  stays pure virtual from Mesh::Iterator::Core
	// void advance ( )  stays pure virtual from Mesh::Iterator::Core
	bool in_range ( );  // virtual from Mesh::Iterator::Core

	class NormalOrder;  class ReverseOrder;
	
};  // end of class Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder
: public Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (here positive vertex)
	// inherited from Mesh::Iterator::Over::CellsOfConnectedOneDimMesh

	// attribute Cell::Core * current_segment
	// inherited from Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh

	inline NormalOrder ( Mesh::Connected::OneDim * m )
	:	Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh ( m ) { };
	
	void reset ( );  // virtual from Mesh::Iterator::Core

	void reset ( const tag::StartAt &, Cell::Core * cll );  // virtual from Mesh::Iterator::Core

	// Cell deref ( )  stays pure virtual from Mesh::Iterator::Core
	void advance ( );  // virtual from Mesh::Iterator::Core
	// bool in_range ( )  virtual, defined by Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh

	class AsTheyAre;  class ForcePositive;
	struct ReverseEachCell  {  class  AssumeCellsExist;  class BuildCells;  };

};  // end of class Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder::AsTheyAre
: public Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (here positive vertex)
	// inherited from Mesh::Iterator::Over::CellsOfConnectedOneDimMesh

	// attribute Cell::Core * current_segment
	// inherited from Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh

	inline AsTheyAre ( Mesh::Connected::OneDim * m )
	:	Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder ( m ) { };
	
	// void reset ( )  defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder

	Cell deref ( );  // virtual from Mesh::Iterator::Core
	// void advance ( )  defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder
	// bool in_range ( )  defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh

};  // end of class Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder::AsTheyAre

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder::ForcePositive
: public Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (here positive vertex)
	// inherited from Mesh::Iterator::Over::CellsOfConnectedOneDimMesh

	// attribute Cell::Core * current_segment
	// inherited from Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh

	inline ForcePositive ( Mesh::Connected::OneDim * m )
	:	Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder ( m ) { };
	
	// void reset ( )  defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder

	Cell deref ( );  // virtual from Mesh::Iterator::Core
	// void advance ( )  defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder
	// bool in_range ( )  defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh

};  // end of class Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder::ForcePositive

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh
          ::NormalOrder::ReverseEachCell::AssumeCellsExist
: public Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (here positive vertex)
	// inherited from Mesh::Iterator::Over::CellsOfConnectedOneDimMesh

	// attribute Cell::Core * current_segment
	// inherited from Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh

	inline AssumeCellsExist ( Mesh::Connected::OneDim * m )
	:	Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder ( m ) { };
	
	// void reset ( )  defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder

	Cell deref ( );  // virtual from Mesh::Iterator::Core
	// void advance ( )  defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder
	// bool in_range ( )  defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh

};  // end of class Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh
    //         ::NormalOrder::ReverseEachCell::AssumeCellsExist

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder::ReverseEachCell::BuildCells
: public Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (here positive vertex)
	// inherited from Mesh::Iterator::Over::CellsOfConnectedOneDimMesh

	// attribute Cell::Core * current_segment
	// inherited from Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh

	inline BuildCells ( Mesh::Connected::OneDim * m )
	:	Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder ( m ) { };
	
	// void reset ( )  defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder

	Cell deref ( );  // virtual from Mesh::Iterator::Core
	// void advance ( )  defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder
	// bool in_range ( )  defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh

};  // end of class Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh
    //         ::NormalOrder::BuildCells::AssumeCellsExist

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder
: public Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (here negative vertex)
	// inherited from Mesh::Iterator::Over::CellsOfConnectedOneDimMesh

	// attribute Cell::Core * current_segment
	// inherited from Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh

	inline ReverseOrder ( Mesh::Connected::OneDim * m )
	:	Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh ( m ) { };
	
	void reset ( );  // virtual from Mesh::Iterator::Core

	void reset ( const tag::StartAt &, Cell::Core * cll );  // virtual from Mesh::Iterator::Core

	// Cell deref ( )  stays pure virtual from Mesh::Iterator::Core
	void advance ( );  // virtual from Mesh::Iterator::Core
	// bool in_range ( )  virtual, defined by Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh

	class AsTheyAre;  class ForcePositive;
	struct ReverseEachCell  {  class  AssumeCellsExist;  class BuildCells;  };

};  // end of class Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder::AsTheyAre
: public Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (here negative vertex)
	// inherited from Mesh::Iterator::Over::CellsOfConnectedOneDimMesh

	// attribute Cell::Core * current_segment
	// inherited from Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh
	
	inline AsTheyAre ( Mesh::Connected::OneDim * m )
	:	Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder ( m ) { };
	
	// void reset ( )  defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder

	Cell deref ( );  // virtual from Mesh::Iterator::Core
	// void advance ( )  defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder
	// bool in_range ( )  virtual, defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh

	class ForcePositive;
	struct ReverseEachCell  {  class AssumeCellsExist;  };

};  // end of class Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder::AsTheyAre

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder::ForcePositive
: public Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (here negative vertex)
	// inherited from Mesh::Iterator::Over::CellsOfConnectedOneDimMesh

	// attribute Cell::Core * current_segment
	// inherited from Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh

	inline ForcePositive ( Mesh::Connected::OneDim * m )
	:	Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder ( m ) { };
	
	// void reset ( )  defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder

	Cell deref ( );  // virtual from Mesh::Iterator::Core
	// void advance ( )  defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder
	// bool in_range ( )  virtual, defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh

};  // end of class Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder::ForcePositive

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh
          ::ReverseOrder::ReverseEachCell::AssumeCellsExist
: public Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (here negative vertex)
	// inherited from Mesh::Iterator::Over::CellsOfConnectedOneDimMesh

	// attribute Cell::Core * current_segment
	// inherited from Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh

	inline AssumeCellsExist ( Mesh::Connected::OneDim * m )
	:	Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder ( m ) { };
	
	// void reset ( )  defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder

	Cell deref ( );  // virtual from Mesh::Iterator::Core
	// void advance ( )  defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder
	// bool in_range ( )  virtual, defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh

};  // end of class Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh
    //         ::ReverseOrder::ReverseEachCell::AssumeCellsExist

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh
          ::ReverseOrder::ReverseEachCell::BuildCells
: public Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (here negative vertex)
	// inherited from Mesh::Iterator::Over::CellsOfConnectedOneDimMesh

	// attribute Cell::Core * current_segment
	// inherited from Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh

	inline BuildCells ( Mesh::Connected::OneDim * m )
	:	Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder ( m ) { };
	
	// void reset ( )  defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder

	Cell deref ( );  // virtual from Mesh::Iterator::Core
	// void advance ( )  defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder
	// bool in_range ( )  virtual, defined by Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh

};  // end of class Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh
    //         ::ReverseOrder::ReverseEachCell::BuildCells


//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::CellsOfFuzzyMesh : public Mesh::Iterator::Core

{	public :

	std::list<Cell> & list;
  std::list<Cell>::iterator iter;

	inline CellsOfFuzzyMesh ( Mesh::Fuzzy * msh, const tag::CellsOfDim &, const size_t d )
	:	Mesh::Iterator::Core (), list { msh->cells[d] }
	{	}  // no need to initialize 'iter', 'reset' will do that
	
	void reset ( );  // virtual from Mesh::Iterator::Core

	void reset ( const tag::StartAt &, Cell::Core * cll );
	// virtual from Mesh::Iterator::Core, here execution forbidden

	// Cell deref ( )  stays pure virtual from Mesh::Iterator::Core
	void advance ( );  // virtual from Mesh::Iterator::Core
	bool in_range ( );  // virtual from Mesh::Iterator::Core

	class AsTheyAre;  class ForcePositive;
	struct ReverseEachCell  {  class AssumeCellsExist;  };
	
};  // end of class Mesh::Iterator::Over::CellsOfFuzzyMesh

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::CellsOfFuzzyMesh::AsTheyAre
: public Mesh::Iterator::Over::CellsOfFuzzyMesh

{	public :

	// inherited from Mesh::Iterator::Over::CellsOfFuzzyMesh :
	// std::list<Cell> & list;
  // std::list<Cell>::iterator iter;

	inline AsTheyAre ( Mesh::Fuzzy * msh, const tag::CellsOfDim &, const size_t d )
		:	Mesh::Iterator::Over::CellsOfFuzzyMesh ( msh, tag::cells_of_dim, d )
	{	}  // no need to initialize 'iter', 'reset' will do that
	
	// void reset ( )  virtual, defined by Mesh::Iterator::Over::CellsOfFuzzyMesh

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by Mesh::Iterator::Over::CellsOfFuzzyMesh, execution forbidden

	Cell deref ( );  // virtual from Mesh::Iterator::Core
	// void advance ( )  virtual, defined by Mesh::Iterator::Over::CellsOfFuzzyMesh
	// bool in_range ( )  virtual, defined by Mesh::Iterator::Over::CellsOfFuzzyMesh

};  // end of class Mesh::Iterator::Over::CellsOfFuzzyMesh::AsTheyAre

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::CellsOfFuzzyMesh::ForcePositive
: public Mesh::Iterator::Over::CellsOfFuzzyMesh

{	public :

	// inherited from Mesh::Iterator::Over::CellsOfFuzzyMesh :
	// std::list<Cell> & list;
  // std::list<Cell>::iterator iter;

	inline ForcePositive ( Mesh::Fuzzy * msh, const tag::CellsOfDim &, const size_t d )
		:	Mesh::Iterator::Over::CellsOfFuzzyMesh ( msh, tag::cells_of_dim, d )
	{	}  // no need to initialize 'iter', 'reset' will do that
	
	// void reset ( )  virtual, defined by Mesh::Iterator::Over::CellsOfFuzzyMesh

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by Mesh::Iterator::Over::CellsOfFuzzyMesh, execution forbidden

	Cell deref ( );  // virtual from Mesh::Iterator::Core
	// void advance ( )  virtual, defined by Mesh::Iterator::Over::CellsOfFuzzyMesh
	// bool in_range ( )  virtual, defined by Mesh::Iterator::Over::CellsOfFuzzyMesh

};  // end of class Mesh::Iterator::Over::CellsOfFuzzyMesh::ForcePositive

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Over::CellsOfFuzzyMesh::ReverseEachCell::AssumeCellsExist
: public Mesh::Iterator::Over::CellsOfFuzzyMesh

{	public :

	// inherited from Mesh::Iterator::Over::CellsOfFuzzyMesh :
	// std::list<Cell> & list;
  // std::list<Cell>::iterator iter;

	inline AssumeCellsExist ( Mesh::Fuzzy * msh, const tag::CellsOfDim &, const size_t d )
		:	Mesh::Iterator::Over::CellsOfFuzzyMesh ( msh, tag::cells_of_dim, d )
	{	}  // no need to initialize 'iter', 'reset' will do that
	
	// void reset ( )  virtual, defined by Mesh::Iterator::Over::CellsOfFuzzyMesh

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by Mesh::Iterator::Over::CellsOfFuzzyMesh, execution forbidden

	Cell deref ( );  // virtual from Mesh::Iterator::Core
	// we trust each cell has already a reverse
	// void advance ( )  virtual, defined by Mesh::Iterator::Over::CellsOfFuzzyMesh
	// bool in_range ( )  virtual, defined by Mesh::Iterator::Over::CellsOfFuzzyMesh

};  // end of class Mesh::Iterator::Over::CellsOfFuzzyMesh
    //                ::ReverseEachCell::AssumeCellsExist


//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

//    ITERATORS OVER CELLS AROUND A GIVEN CELL (CENTERED AT CELL)

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


// see paragraph 9.12 in the manual

class Mesh::Iterator::Around::OneCell : public Mesh::Iterator::Core

// we may want to run over all segments starting at a given vertex
// or over all squares sharing a given edge (in a mesh of cubes)

// the dimension of the cells to be returned varies
// e.g. we can return segments or squares (or cubes) around a vertex
// or squares (or cubes) around a segment

{	public :

	Mesh::Core * msh_p;  // should an iterator keep the mesh alive ?
	Cell::Positive * center;  // center->dim <= msh->dim - 2

	inline OneCell ( Mesh::Core * msh, Cell::Positive * const c )
	:	Mesh::Iterator::Core (), msh_p { msh }, center { c } { }

	// void reset ( )  stays pure virtual from Mesh::Iterator::Core
	// void reset ( tag::StartAt, Cell::Core * cll )  stays pure virtual from Mesh::Iterator::Core
	
	// Cell deref ( )  stays pure virtual from Mesh::Iterator::Core
	// void advance ( )  stays pure virtual from Mesh::Iterator::Core
	// bool in_range ( )  stays pure virtual from Mesh::Iterator::Core

	struct OfCodimTwo
	{ class Ordered; class OverVertices; class OverSegments;
		struct WorkAround2D { struct NormalOrder { class BuildReverseCells; };
		                      struct ReverseOrder { class BuildReverseCells; }; };	};
	class OfAnyCodim;
	
};  // end of class Mesh::Iterator::Around::OneCell

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered
: public Mesh::Iterator::Around::OneCell

// when cll->dim == msh->dim - 2, there is a linear order
// e.g. we rotate around a vertex in a 2D mesh
// or we rotate around a segment in a 3D mesh

// the dimension of the cells to be returned is a different matter
// e.g. we can return segments or squares around a vertex in a 2D mesh
// or squares or cubes around a segment in 3D mesh

// below, words "vertex" and "segment" are quite misleading
// they come from an analogy to iterators over chains of segments
// "vertex" means a cell of dimension center->dim + 1 == msh->dim - 1
// "segment" means a cell of dimension center->dim + 2 == msh->dim
// when we rotate around a vertex, "vertex" means segment, "segment" means square
// when we rotate around a segment, "vertex" means square, "segment" means cube

{	public :

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	Cell::Core * current_segment { nullptr };  // word "segment" is misleading here

	inline Ordered ( Mesh::Core * msh, Cell::Positive * const c )
	:	Mesh::Iterator::Around::OneCell ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  stays pure virtual from Mesh::Iterator::Core
	// void reset ( tag::StartAt, Cell::Core * cll )  stays pure virtual from Mesh::Iterator::Core
	
	// void reset ( tag::StartAt, Cell::Core * cll )  stays pure virtual from Mesh::Iterator::Core
	// Cell deref ( )  stays pure virtual from Mesh::Iterator::Core
	// void advance ( )  stays pure virtual from Mesh::Iterator::Core
	// bool in_range ( )  stays virtual from Mesh::Iterator::Core
	
};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices
: public Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered

// when cll->dim == msh->dim - 2, there is a linear order
// e.g. we rotate around a vertex in a 2D mesh
// or we rotate around a segment in a 3D mesh

// the dimension of the cells to be returned is a different matter
// here we return segments around a vertex or squares around a segment

// "Vertices" means we want to get cells of dimension center->dim + 1 == msh->dim - 1
// the procedure is analogous to Mesh::Iterator::Over::VerticesOfConnectedOneDimMesh

// here, words "vertex" and "segment" are quite misleading
// they come from an analogy to iterators over chains of segments
// "vertex" means a cell of dimension center->dim + 1 == msh->dim - 1
// "segment" means a cell of dimension center->dim + 2 == msh->dim
// when we rotate around a vertex, "vertex" means segment, "segment" means square
// when we rotate around a segment, "vertex" means square, "segment" means cube

{	public :

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	Cell::Core * current_vertex { nullptr };  // word "vertex" is misleading here
	Cell::Core * first_vertex { nullptr };  // word "vertex" is misleading here
	Cell::Core * last_segment { nullptr };  // word "segment" is misleading here

	inline OverVertices ( Mesh::Core * msh, Cell::Positive * const c )
	:	Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  both stay pure virtual from Mesh::Iterator::Core
	// void reset ( const tag::StartAt &, Cell::Core * cll ) 
	
	// Cell deref ( )  stays pure virtual from Mesh::Iterator::Core
	// void advance ( )  stays pure virtual from Mesh::Iterator::Core

	bool in_range ( );  // virtual from Mesh::Iterator::Core
	
	class NormalOrder;  class ReverseOrder;

};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices


//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder
: public Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices

// when cll->dim == msh->dim - 2, there is a linear order
// e.g. we rotate around a vertex in a 2D mesh
// or we rotate around a segment in a 3D mesh
// here we follow this order

// the dimension of the cells to be returned is a different matter
// here we return segments around a vertex or squares around a segment

{	public :

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	// we recognize a closed loop if  last_segment == msh.cell_behind ( first_vertex )
	// another possible criterion :  last_segment == nullptr  ==>  open chain

	inline NormalOrder ( Mesh::Core * msh, Cell::Positive * const c )
	:	Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices ( msh, c )
	// msh_p { msh }, center { c }
	{ }
	
	void reset ( );  // virtual from Mesh::Iterator::Core
	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   stays pure virtual from Mesh::Iterator::Core

	// Cell deref ( )  stays pure virtual from Mesh::Iterator::Core

	void advance ( );  // virtual from Mesh::Iterator::Core
	
	// bool in_range ( )  virtual, defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices
	
	class AssumeCellsExist;  class BuildReverseCells;

};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder
: public Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices

// when cll->dim == msh->dim - 2, there is a linear order
// e.g. we rotate around a vertex in a 2D mesh
// or we rotate around a segment in a 3D mesh
// here we follow the reversed order

// the dimension of the cells to be returned is a different matter
// here we return segments around a vertex or squares around a segment

{	public :

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	// we recognize a closed loop if  last_segment == msh.cell_in_front_of ( first_vertex )
	// another possible criterion :  last_segment == nullptr  ==>  open chain

	inline ReverseOrder ( Mesh::Core * msh, Cell::Positive * const c );
	
	void reset ( );  // virtual from Mesh::Iterator::Core
	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   stays pure virtual from Mesh::Iterator::Core
	
	// Cell deref ( )  stays pure virtual from Mesh::Iterator::Core

	// void advance ( )  stays pure virtual from Mesh::Iterator::Core
	
	// bool in_range ( )  virtual, defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices

	class AssumeCellsExist;  class BuildReverseCells;

};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder::AssumeCellsExist
: public Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder

{	public :

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline AssumeCellsExist ( Mesh::Core * msh, Cell::Positive * const c );
	
	// void reset ( )  virtual,
	//   defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder
	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   stays pure virtual from Mesh::Iterator::Core
	
	// Cell deref ( )  stays pure virtual from Mesh::Iterator::Core

	// void advance ( )  virtual from Mesh::Iterator::Core, defined by
	//   Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder

	// bool in_range ( )  virtual, defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices

	class AsTheyAre;  class ForcePositive;  class ReverseEachCell;

};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo
    //                   ::OverVertices::NormalOrder::AssumeCellsExist

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder::BuildReverseCells
: public Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder

{	public :

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline BuildReverseCells ( Mesh::Core * msh, Cell::Positive * const c );
	
	// void reset ( )  virtual,
	//   defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder
	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   stays pure virtual from Mesh::Iterator::Core
	
	// Cell deref ( )  stays pure virtual from Mesh::Iterator::Core

	// void advance ( )  virtual from Mesh::Iterator::Core, defined by
	//   Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder

	// bool in_range ( )  virtual, defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices

	class AsTheyAre;  class ForcePositive;  class ReverseEachCell;

};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo
    //                   ::OverVertices::NormalOrder::BuildReverseCells

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

namespace tag::local_functions { namespace iterator_h { }; }

namespace tag::local_functions::iterator_h {

inline Cell::Core * find_first_seg ( Cell::Positive * const pos_cen, Mesh::Core * const msh )

{	typedef std::map < Mesh::Core*, Cell::field_to_meshes > maptype;
	typedef std::map < Mesh::Core*, Cell::field_to_meshes_same_dim > maptype_sd;
	maptype & cm2 = pos_cen->meshes[1];
	typename maptype::iterator it = cm2.begin();
	for ( ; it != cm2.end(); it++ )
	{	Mesh::Core * m = it->first;
		// is 'm' the boundary of some cell ?
		Cell::Positive * mce = m->cell_enclosed;
		if ( mce == nullptr ) continue;
		#ifndef NDEBUG
		// a boundary has no boundary, so it will appear exactly twice, with opposite orientations
		// do we accept only 11 here, or are 20 and 02 acceptable too ?
		// Cell::field_to_meshes f = it->second;
		// assert ( f.counter_pos == 1 );
		// assert ( f.counter_neg == 1 );
		#endif
		// does mce belong to msh ? note that mce->dim == msh->dim and mce->is_positive
		Cell::Positive::NotVertex * mce_nv = tag::Util::assert_cast
			< Cell::Core*, Cell::Positive::NotVertex* > ( mce );
		maptype_sd & mce_msd = mce_nv->meshes_same_dim;
		typename maptype_sd::iterator itt = mce_msd.find(msh);
		if ( itt == mce_msd.end() ) continue;
		assert ( itt->first == msh );
		// yes ! mce may be current_seg, if correctly oriented
		Cell::field_to_meshes_same_dim ff = itt->second;
		if ( ff.sign == 1 ) return mce;
		else assert ( ff.sign == -1 );                               }
	return nullptr;                                                             }


inline void rotate_backwards_for_ver
( Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices * that )

// that->first_vertex should be reversed by calling function
	
{	Cell this_center ( tag::whose_core_is, that->center,
	                   tag::previously_existing, tag::surely_not_null );
	Mesh m ( tag::whose_core_is, that->msh_p, tag::previously_existing, tag::is_positive );
	assert ( that->current_segment );
	Cell::Core * vertex_stop = that->current_segment->boundary().cell_behind
		( this_center, tag::surely_exists ).core;
	while ( true )
	{	Cell::Core * vertex_p = that->current_segment->boundary().cell_in_front_of
			( this_center, tag::surely_exists ).core;
		Cell vertex ( tag::whose_core_is, vertex_p,
		              tag::previously_existing, tag::surely_not_null );
		Cell new_seg = m.cell_in_front_of ( vertex, tag::may_not_exist );
		if ( not new_seg.exists() )  // we have hit the boundary
		{	that->first_vertex = vertex_p;
			that->last_segment = nullptr;
			return;                          }
		if ( vertex_p->reverse_attr.core == vertex_stop )  // completed loop
		{	that->last_segment = m.cell_in_front_of ( vertex ).core;
			that->first_vertex = vertex_p;
			return;                                                  }
		that->current_segment = new_seg.core;                                                  }  }

// the above could get simpler and faster by using a method like
// inline Cell::Core * Mesh::Core::cell_behind_ptr
// ( const Cell::Core * face, const tag::MayNotExist & ) const
// which is not difficult to implement


inline void rotate_forward_for_ver
( Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices * that )
	
{	Cell this_center ( tag::whose_core_is, that->center,
	                   tag::previously_existing, tag::surely_not_null );
	Mesh m ( tag::whose_core_is, that->msh_p, tag::previously_existing, tag::is_positive );
	assert ( that->current_segment );
	Cell::Core * vertex_stop = that->current_segment->boundary().cell_in_front_of
		( this_center, tag::surely_exists ).core;
	while ( true )
	{	Cell::Core * vertex_p = that->current_segment->boundary().cell_behind
			( this_center, tag::surely_exists ).core;
		Cell vertex ( tag::whose_core_is, vertex_p,
		              tag::previously_existing, tag::surely_not_null );
		Cell new_seg = m.cell_in_front_of ( vertex, tag::may_not_exist );
		if ( not new_seg.exists() )  // we have hit the boundary
		{	that->first_vertex = vertex_p;
			that->last_segment = nullptr;
			return;                          }
		if ( vertex_p->reverse_attr.core == vertex_stop )  // completed loop
		{	that->last_segment = m.cell_in_front_of ( vertex ).core;
			that->first_vertex = vertex_p;
			return;                                                  }
		that->current_segment = new_seg.core;                                                 }  }

// the above could get simpler and faster by using a method like
// inline Cell::Core * Mesh::Core::cell_behind_ptr
// ( const Cell::Core * face, const tag::MayNotExist & ) const
// which is not difficult to implement
	

inline bool inner_cell ( Mesh::Core * const msh, Cell::Positive * const cen )

// we want to know whether we are in the interior or on the boundary of the mesh
// so we search for some first_segment, then rotate
	
{	assert ( msh );
	assert ( cen );
	assert ( cen->get_dim() + 3 == msh->get_dim_plus_one() );
	Cell::Core * first_segment =
		tag::local_functions::iterator_h::find_first_seg ( cen, msh );
	assert ( first_segment );

	// now we rotate forward until we meet the boundary or close the loop
	Cell this_center ( tag::whose_core_is, cen,
	                   tag::previously_existing, tag::surely_not_null );
	Mesh m ( tag::whose_core_is, msh, tag::previously_existing, tag::is_positive );
	Cell::Core * vertex_stop = first_segment->boundary().cell_in_front_of
		( this_center, tag::surely_exists ).core;

	while ( true )
	{	Cell::Core * vertex_p = first_segment->boundary().cell_behind
			( this_center, tag::surely_exists ).core;
		Cell vertex ( tag::whose_core_is, vertex_p,
		              tag::previously_existing, tag::surely_not_null );
		Cell new_seg = m.cell_in_front_of ( vertex, tag::may_not_exist );
		if ( not new_seg.exists() ) return false;  // we have hit the boundary
		if ( vertex_p->reverse_attr.core == vertex_stop )  // completed loop
			return true;
		first_segment = new_seg.core;                                                 }  }

// the above could get simpler and faster by using a method like
// inline Cell::Core * Mesh::Core::cell_behind_ptr
// ( const Cell::Core * face, const tag::MayNotExist & ) const
// which is not difficult to implement
	
}  // namespace tag::local_functions::iterator_h

//-----------------------------------------------------------------------------------------


inline bool Cell::is_inner_to ( const Mesh & msh ) const
// builds reverse of 'cen' !
{	assert ( msh.exists() );
	assert ( this->exists() );
	assert ( this->dim() < msh.dim() );
	if ( this->dim() + 2 == msh.core->get_dim_plus_one() )  // co-dimension one
		return ( msh.cell_in_front_of ( *this, tag::may_not_exist ) .exists() and
		         msh.cell_behind ( *this, tag::may_not_exist ) .exists()          );
	assert ( this->dim() + 3 == msh.core->get_dim_plus_one() );  // co-dimension two
	Cell::Positive * cll_pos = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive* > ( this->core->get_positive().core );
	return tag::local_functions::iterator_h::inner_cell ( msh.core, cll_pos );  }

// add code for co-dimension tree or higher ! how do we do that ?

//-----------------------------------------------------------------------------------------


inline Mesh::Iterator::Around::OneCell::OfCodimTwo
           ::OverVertices::NormalOrder::AssumeCellsExist::AssumeCellsExist
( Mesh::Core * msh, Cell::Positive * cen )
:	Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder ( msh, cen )
// msh_p { msh }, center { cen }

// we don't know whether we are in the interior or on the boundary of the mesh
// so we must search for this->first_vertex, rotating backwards
	
{	assert ( this->msh_p );
	assert ( cen );
	assert ( cen->get_dim() + 3 == this->msh_p->get_dim_plus_one() );
	this->center = cen;
	this->current_segment = tag::local_functions::iterator_h::find_first_seg ( cen, this->msh_p );
	assert ( this->current_segment );

	// now we rotate backwards until we meet the boundary or close the loop
	tag::local_functions::iterator_h::rotate_backwards_for_ver ( this );

	assert ( this->first_vertex->reverse_attr.core );
	this->first_vertex = this->first_vertex->reverse_attr.core;                           }

// the above could get simpler and faster by using a method like
// inline Cell::Core * Mesh::Core::cell_behind_ptr
// ( const Cell::Core * face, const tag::MayNotExist & ) const
// which is not difficult to implement
	

inline Mesh::Iterator::Around::OneCell::OfCodimTwo
           ::OverVertices::NormalOrder::BuildReverseCells::BuildReverseCells
( Mesh::Core * msh, Cell::Positive * cen )
:	Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder ( msh, cen )
// msh_p { msh }, center { cen }

// we don't know whether we are in the interior or on the boundary of the mesh
// so we must search for this->first_vertex, rotating backwards
	
{	assert ( this->msh_p );
	assert ( cen );
	assert ( cen->get_dim() + 3 == this->msh_p->get_dim_plus_one() );
	this->center = cen;
	this->current_segment = tag::local_functions::iterator_h::find_first_seg ( cen, this->msh_p );
	assert ( this->current_segment );

	// now we rotate backwards until we meet the boundary or close the loop
	tag::local_functions::iterator_h::rotate_backwards_for_ver ( this );
	// we do not reverse this->first_vertex
	// instead, this->reset() will reverse it while assigning to this->current_vertex
	// however, we build the reverse if it does not exist
	
	if ( not this->first_vertex->reverse_attr.exists() )
	//	this->first_vertex->reverse_attr =
	//		Cell ( tag::whose_core_is, this->first_vertex->build_reverse
	//	         ( tag::one_dummy_wrapper ), tag::freshly_created         );
	{	Cell::Positive * tfv = tag::Util::assert_cast
			< Cell::Core*, Cell::Positive* > ( this->first_vertex );
		tfv->reverse_attr.core = tfv->build_reverse ( tag::one_dummy_wrapper );  }
	
	assert ( this->first_vertex->reverse_attr.core );
	this->first_vertex = this->first_vertex->reverse_attr.core;                   }

// the above could get simpler and faster by using a method like
// inline Cell::Core * Mesh::Core::cell_behind_ptr
// ( const Cell::Core * face, const tag::MayNotExist & ) const
// which is not difficult to implement
	
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfCodimTwo
         ::OverVertices::NormalOrder::AssumeCellsExist::AsTheyAre
: public Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder::AssumeCellsExist

{	public :

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline AsTheyAre ( Mesh::Core * msh, Cell::Positive * const c )
	:	Mesh::Iterator::Around::OneCell::OfCodimTwo
	      ::OverVertices::NormalOrder::AssumeCellsExist ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual,
	//   defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder
	void reset ( const tag::StartAt &, Cell::Core * cll );  // virtual from Mesh::Iterator::Core
	
	Cell deref ( );  // virtual from Mesh::Iterator::Core

	// void advance ( )  virtual from Mesh::Iterator::Core, defined by
	//   Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder

	// bool in_range ( )  virtual, defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices

};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo
    //                   ::OverVertices::NormalOrder::AssumeCellsExist::AsTheyAre

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfCodimTwo
         ::OverVertices::NormalOrder::AssumeCellsExist::ForcePositive
: public Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder::AssumeCellsExist

{	public :

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline ForcePositive ( Mesh::Core * msh, Cell::Positive * const c )
	:	Mesh::Iterator::Around::OneCell::OfCodimTwo
	      ::OverVertices::NormalOrder::AssumeCellsExist ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual,
	//   defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder
	void reset ( const tag::StartAt &, Cell::Core * cll );  // virtual from Mesh::Iterator::Core
	
	Cell deref ( );  // virtual from Mesh::Iterator::Core

	// void advance ( )  virtual from Mesh::Iterator::Core, defined by
	//   Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder

	// bool in_range ( )  virtual, defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices
	
};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo
    //                  ::OverVertices::NormalOrder::AssumeCellsExist::ForcePositive

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfCodimTwo
         ::OverVertices::NormalOrder::AssumeCellsExist::ReverseEachCell
: public Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder::AssumeCellsExist

{	public :

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline ReverseEachCell ( Mesh::Core * msh, Cell::Positive * const c )
	:	Mesh::Iterator::Around::OneCell::OfCodimTwo
	      ::OverVertices::NormalOrder::AssumeCellsExist ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual,
	//   defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder
	void reset ( const tag::StartAt &, Cell::Core * cll );  // virtual from Mesh::Iterator::Core
	
	Cell deref ( );  // virtual from Mesh::Iterator::Core

	// void advance ( )  virtual from Mesh::Iterator::Core, defined by
	//   Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder

	// bool in_range ( )  virtual, defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices
	
};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo
    //                  ::OverVertices::NormalOrder::AssumeCellsExist::ReverseEachCell

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfCodimTwo
         ::OverVertices::NormalOrder::BuildReverseCells::AsTheyAre
: public Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder::BuildReverseCells

{	public :

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline AsTheyAre ( Mesh::Core * msh, Cell::Positive * const c )
	:	Mesh::Iterator::Around::OneCell::OfCodimTwo
	      ::OverVertices::NormalOrder::BuildReverseCells ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual,
	//   defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder
	void reset ( const tag::StartAt &, Cell::Core * cll );  // virtual from Mesh::Iterator::Core
	
	Cell deref ( );  // virtual from Mesh::Iterator::Core

	// void advance ( )  virtual from Mesh::Iterator::Core, defined by
	//   Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder

	// bool in_range ( )  virtual, defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices

};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo
    //                   ::OverVertices::NormalOrder::BuildReverseCells::AsTheyAre

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfCodimTwo
         ::OverVertices::NormalOrder::BuildReverseCells::ForcePositive
: public Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder::BuildReverseCells

{	public :

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline ForcePositive ( Mesh::Core * msh, Cell::Positive * const c )
	:	Mesh::Iterator::Around::OneCell::OfCodimTwo
	      ::OverVertices::NormalOrder::BuildReverseCells ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual,
	//   defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder
	void reset ( const tag::StartAt &, Cell::Core * cll );  // virtual from Mesh::Iterator::Core
	
	Cell deref ( );  // virtual from Mesh::Iterator::Core

	// void advance ( )  virtual from Mesh::Iterator::Core, defined by
	//   Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder

	// bool in_range ( )  virtual, defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices
	
};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo
    //                  ::OverVertices::NormalOrder::BuildReverseCells::ForcePositive

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfCodimTwo
         ::OverVertices::NormalOrder::BuildReverseCells::ReverseEachCell
: public Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder::BuildReverseCells

{	public :

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline ReverseEachCell ( Mesh::Core * msh, Cell::Positive * const c )
	:	Mesh::Iterator::Around::OneCell::OfCodimTwo
	      ::OverVertices::NormalOrder::BuildReverseCells ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual,
	//   defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder
	void reset ( const tag::StartAt &, Cell::Core * cll );  // virtual from Mesh::Iterator::Core
	
	Cell deref ( );  // virtual from Mesh::Iterator::Core

	// void advance ( )  virtual from Mesh::Iterator::Core, defined by
	//   Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder

	// bool in_range ( )  virtual, defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices
	
};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo
    //                  ::OverVertices::NormalOrder::BuildReverseCells::ReverseEachCell

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::AssumeCellsExist
: public Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder

{	public :

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline AssumeCellsExist ( Mesh::Core * msh, Cell::Positive * const c )
	:	Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder ( msh, c )
	// msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual,
	//   defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder
	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   stays pure virtual from Mesh::Iterator::Core
	
	// Cell deref ( )  stays pure virtual from Mesh::Iterator::Core

	void advance ( );  // virtual from Mesh::Iterator::Core

	// bool in_range ( )  virtual, defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices

	class AsTheyAre;  class ForcePositive;  class ReverseEachCell;

};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo
    //                   ::OverVertices::ReverseOrder::AssumeCellsExist

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::BuildReverseCells
: public Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder

{	public :

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline BuildReverseCells ( Mesh::Core * msh, Cell::Positive * const c )
	:	Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder ( msh, c )
	// msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual,
	//   defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverselOrder
	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   stays pure virtual from Mesh::Iterator::Core
	
	// Cell deref ( )  stays pure virtual from Mesh::Iterator::Core

	void advance ( );  // virtual from Mesh::Iterator::Core

	// bool in_range ( )  virtual, defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices

	class AsTheyAre;  class ForcePositive;  class ReverseEachCell;

};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo
    //                   ::OverVertices::ReverseOrder::BuildReverseCells

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


inline Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::ReverseOrder
( Mesh::Core * msh, Cell::Positive * cen )
:	Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices ( msh, cen )
// msh_p { msh }, center { cen }

// we don't know whether we are in the interior or on the boundary of the mesh
// so we must search for this->first_vertex, rotating backwards
	
{	assert ( this->msh_p );
	assert ( cen );
	assert ( cen->get_dim() + 3 == this->msh_p->get_dim_plus_one() );
	this->center = cen;
	this->current_segment = tag::local_functions::iterator_h::find_first_seg ( cen, this->msh_p );
	assert ( this->current_segment );

	// now we rotate forward until we meet the boundary or close the loop
	tag::local_functions::iterator_h::rotate_forward_for_ver ( this );                              }

// the above could get simpler and faster by using a method like
// inline Cell::Core * Mesh::Core::cell_behind_ptr
// ( const Cell::Core * face, const tag::MayNotExist & ) const
// which is not difficult to implement
	
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfCodimTwo
         ::OverVertices::ReverseOrder::AssumeCellsExist::AsTheyAre
: public Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::AssumeCellsExist

{	public :

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline AsTheyAre ( Mesh::Core * msh, Cell::Positive * const c )
	:	Mesh::Iterator::Around::OneCell::OfCodimTwo
		::OverVertices::ReverseOrder::AssumeCellsExist ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual,
	//   defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder
	void reset ( const tag::StartAt &, Cell::Core * cll );  // virtual from Mesh::Iterator::Core
	
	Cell deref ( );  // virtual from Mesh::Iterator::Core

	// void advance ( )  virtual from Mesh::Iterator::Core, defined by
	//   Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::AssumeCellsExist

	// bool in_range ( )  virtual, defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices

};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo
    //                   ::OverVertices::ReverseOrder::AssumeCellsExist::AsTheyAre

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfCodimTwo
         ::OverVertices::ReverseOrder::AssumeCellsExist::ForcePositive
: public Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::AssumeCellsExist

{	public :

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline ForcePositive ( Mesh::Core * msh, Cell::Positive * const c )
	:	Mesh::Iterator::Around::OneCell::OfCodimTwo
		::OverVertices::ReverseOrder::AssumeCellsExist ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual,
	//   defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder
	void reset ( const tag::StartAt &, Cell::Core * cll );  // virtual from Mesh::Iterator::Core
	
	Cell deref ( );  // virtual from Mesh::Iterator::Core

	// void advance ( )  virtual from Mesh::Iterator::Core, defined by
	//   Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::AssumeCellsExist

	// bool in_range ( )  virtual, defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices
	
};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo
    //                  ::OverVertices::ReverseOrder::AssumeCellsExist::ForcePositive

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfCodimTwo
         ::OverVertices::ReverseOrder::AssumeCellsExist::ReverseEachCell
: public Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::AssumeCellsExist

{	public :

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline ReverseEachCell ( Mesh::Core * msh, Cell::Positive * const c )
	:	Mesh::Iterator::Around::OneCell::OfCodimTwo
		::OverVertices::ReverseOrder::AssumeCellsExist ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual,
	//   defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder
	void reset ( const tag::StartAt &, Cell::Core * cll );  // virtual from Mesh::Iterator::Core
	
	Cell deref ( );  // virtual from Mesh::Iterator::Core

	// void advance ( )  virtual from Mesh::Iterator::Core, defined by
	//   Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::AssumeCellsExist

	// bool in_range ( )  virtual, defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices
	
};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo
    //                  ::OverVertices::ReverseOrder::AssumeCellsExist::ReverseEachCell

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfCodimTwo
         ::OverVertices::ReverseOrder::BuildReverseCells::AsTheyAre
: public Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::BuildReverseCells

{	public :

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline AsTheyAre ( Mesh::Core * msh, Cell::Positive * const c )
	:	Mesh::Iterator::Around::OneCell::OfCodimTwo
		::OverVertices::ReverseOrder::BuildReverseCells ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual,
	//   defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder
	void reset ( const tag::StartAt &, Cell::Core * cll );  // virtual from Mesh::Iterator::Core
	
	Cell deref ( );  // virtual from Mesh::Iterator::Core

	// void advance ( )  virtual from Mesh::Iterator::Core, defined by
	//   Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::BuildReverseCells

	// bool in_range ( )  virtual, defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices

};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo
    //                   ::OverVertices::ReverseOrder::BuildReverseCells::AsTheyAre

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfCodimTwo
         ::OverVertices::ReverseOrder::BuildReverseCells::ForcePositive
: public Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::BuildReverseCells

{	public :

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline ForcePositive ( Mesh::Core * msh, Cell::Positive * const c )
	:	Mesh::Iterator::Around::OneCell::OfCodimTwo
		::OverVertices::ReverseOrder::BuildReverseCells ( msh,c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual,
	//   defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder
	void reset ( const tag::StartAt &, Cell::Core * cll );  // virtual from Mesh::Iterator::Core
	
	Cell deref ( );  // virtual from Mesh::Iterator::Core

	// void advance ( )  virtual from Mesh::Iterator::Core, defined by
	//   Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::BuildReverseCells

	// bool in_range ( )  virtual, defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices
	
};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo
    //                  ::OverVertices::ReverseOrder::BuildReverseCells::ForcePositive

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfCodimTwo
         ::OverVertices::ReverseOrder::BuildReverseCells::ReverseEachCell
: public Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::BuildReverseCells

{	public :

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline ReverseEachCell ( Mesh::Core * msh, Cell::Positive * const c )
	:	Mesh::Iterator::Around::OneCell::OfCodimTwo
		::OverVertices::ReverseOrder::BuildReverseCells ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual,
	//   defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder
	void reset ( const tag::StartAt &, Cell::Core * cll );  // virtual from Mesh::Iterator::Core
	
	Cell deref ( );  // virtual from Mesh::Iterator::Core

	// void advance ( )  virtual from Mesh::Iterator::Core, defined by
	//   Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::BuildReverseCells

	// bool in_range ( )  virtual from Mesh::Iterator::Core,
  //   defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices
	
};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo
    //                  ::OverVertices::ReverseOrder::BuildReverseCells::ReverseEachCell

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments
: public Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered

// when cll->dim == msh->dim - 2, there is a linear order
// e.g. we rotate around a vertex in a 2D mesh
// or we rotate around a segment in a 3D mesh

// the dimension of the cells to be returned is a different matter
// here we return squares around a vertex or cubes around a segment

// "Segments" means we want to get cells of dimension center->dim + 2 == msh->dim
// the procedure is analogous to Mesh::Iterator::Over::SegmentsOfConnectedOneDimMesh

// here, words "vertex" and "segment" are quite misleading
// they come from an analogy to iterators over chains of segments
// "vertex" means a cell of dimension center->dim + 1 == msh->dim - 1
// "segment" means a cell of dimension center->dim + 2 == msh->dim
// when we rotate around a vertex, "vertex" means segment, "segment" means square
// when we rotate around a segment, "vertex" means square, "segment" means cube

{	public :

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here

	Cell::Core * first_segment { nullptr };  // word "segment" is misleading here
	Cell::Core * last_vertex { nullptr };  // word "vertex" is misleading here

	inline OverSegments ( Mesh::Core * msh, Cell::Positive * const c )
	:	Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered ( msh, c )  // msh_p { msh }, center { c }
	{ }

	void reset ( );  // virtual from Mesh::Iterator::Core
	
	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//      stays pure virtual from Mesh::Iterator::Core
	
	// Cell deref ( )  stays pure virtual from Mesh::Iterator::Core
	// void advance ( )  stays pure virtual from Mesh::Iterator::Core

	bool in_range ( );  // virtual from Mesh::Iterator::Core
	
	class NormalOrder;  class ReverseOrder;

};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments::NormalOrder
: public Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments

// when cll->dim == msh->dim - 2, there is a linear order
// e.g. we rotate around a vertex in a 2D mesh
// or we rotate around a segment in a 3D mesh
// here we follow this order

// the dimension of the cells to be returned is a different matter
// here we return squares around a vertex or cubes around a segment

{	public :

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here

	// attributes inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments :
	// Cell::Core * first_segment { nullptr }  -- word "segment" is misleading here
	// Cell::Core * last_vertex { nullptr }  -- word "vertex" is misleading here

	inline NormalOrder ( Mesh::Core * msh, Cell::Positive * const c );
	
	// void reset ( )  -- virtual,
	//    defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments

	void reset ( const tag::StartAt &, Cell::Core * cll );  // virtual from Mesh::Iterator::Core

	// Cell deref ( )  stays pure virtual from Mesh::Iterator::Core

	void advance ( );  // virtual from Mesh::Iterator::Core
	
	// bool in_range ( )  virtual, defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments
	
	class AsTheyAre;  class ForcePositive;  class ReverseEachCell;

};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments::NormalOrder

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments::ReverseOrder
: public Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments

// when cll->dim == msh->dim - 2, there is a linear order
// e.g. we rotate around a vertex in a 2D mesh
// or we rotate around a segment in a 3D mesh
// here we follow the reversed order

// the dimension of the cells to be returned is a different matter
// here we return squares around a vertex or cubes around a segment

{	public :

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here

	// attributes inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments :
	// Cell::Core * first_segment { nullptr }  -- word "segment" is misleading here
	// Cell::Core * last_vertex { nullptr }  -- word "vertex" is misleading here

	inline ReverseOrder ( Mesh::Core * msh, Cell::Positive * const c );
	
	// void reset ( )  -- virtual,
	//    defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments

	void reset ( const tag::StartAt &, Cell::Core * cll );  // virtual from Mesh::Iterator::Core
	
	// Cell deref ( )  stays pure virtual from Mesh::Iterator::Core

	void advance ( );  // virtual from Mesh::Iterator::Core
	
	// bool in_range ( )  virtual, defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments

	class AsTheyAre;  class ForcePositive;  class ReverseEachCell;

};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments::ReverseOrder

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


namespace tag::local_functions::iterator_h {
		
inline void rotate_backwards_for_seg ( Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments * that )

{	Cell this_center ( tag::whose_core_is, that->center,
	                   tag::previously_existing, tag::surely_not_null );
	Mesh m ( tag::whose_core_is, that->msh_p, tag::previously_existing, tag::is_positive );
	assert ( that->first_segment );
	Cell::Core * vertex_stop = that->first_segment->boundary().cell_behind
		( this_center, tag::surely_exists ).core;
	while ( true )
	{	Cell::Core * vertex_p = that->first_segment->boundary().cell_in_front_of
			( this_center, tag::surely_exists ).core;
		Cell vertex ( tag::whose_core_is, vertex_p,
	                tag::previously_existing, tag::surely_not_null );
		Cell new_seg = m.cell_in_front_of ( vertex, tag::may_not_exist );
		if ( not new_seg.exists() ) return;  // we have hit the boundary
		        // that->last_vertex remains nullptr
		if ( vertex_p->reverse_attr.core == vertex_stop )  // completed loop
		{	that->last_vertex = vertex_stop;  return; }
		that->first_segment = new_seg.core;                                         }           }

// the above could get simpler and faster by using a method like
// inline Cell::Core * Mesh::Core::cell_behind_ptr
// ( const Cell::Core * face, const tag::MayNotExist & ) const
// which is not difficult to implement


inline void rotate_forward_for_seg
( Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments * that )

{	Cell this_center ( tag::whose_core_is, that->center,
	                   tag::previously_existing, tag::surely_not_null );
	Mesh m ( tag::whose_core_is, that->msh_p, tag::previously_existing, tag::is_positive );
	assert ( that->first_segment );
	Cell::Core * vertex_stop = that->first_segment->boundary().cell_in_front_of
		( this_center, tag::surely_exists ).core;
	while ( true )
	{	Cell::Core * vertex_p = that->first_segment->boundary().cell_behind
			( this_center, tag::surely_exists ).core;
		Cell vertex ( tag::whose_core_is, vertex_p,
	                tag::previously_existing, tag::surely_not_null );
		Cell new_seg = m.cell_in_front_of ( vertex, tag::may_not_exist );
		if ( not new_seg.exists() ) return;  // we have hit the boundary
		        // that->last_vertex remains nullptr
		if ( vertex_p->reverse_attr.core == vertex_stop )  // completed loop
		{	that->last_vertex = vertex_stop;  return; }
		that->first_segment = new_seg.core;                                         }           }

// the above could get simpler and faster by using a method like
// inline Cell::Core * Mesh::Core::cell_behind_ptr
// ( const Cell::Core * face, const tag::MayNotExist & ) const
// which is not difficult to implement

}  // namespace tag::local_functions::iterator_h
	
//-----------------------------------------------------------------------------------------


inline Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments::NormalOrder::NormalOrder
( Mesh::Core * msh, Cell::Positive * cen )
:	Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments ( msh, cen )
// msh_p { msh }, center { cen }

// we don't know whether we are in the interior or on the boundary of the mesh
// so we must search for this->first_vertex, rotating backwards
	
{	assert ( this->msh_p );
	assert ( cen );
	assert ( cen->get_dim() + 3 == this->msh_p->get_dim_plus_one() );
	this->center = cen;
	this->first_segment = tag::local_functions::iterator_h::find_first_seg ( cen, this->msh_p );
	assert ( this->first_segment );

	// now we rotate backwards until we meet the boundary or close the loop
	tag::local_functions::iterator_h::rotate_backwards_for_seg ( this );
	assert ( this->first_segment );                                                              }
	

inline Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments::ReverseOrder::ReverseOrder
( Mesh::Core * msh, Cell::Positive * cen )
:	Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments ( msh, cen )
// msh_p { msh }, center { cen }

// we don't know whether we are in the interior or on the boundary of the mesh
// so we must search for this->first_vertex, rotating forward
	
{	assert ( this->msh_p );
	assert ( cen );
	assert ( cen->get_dim() + 3 == this->msh_p->get_dim_plus_one() );
	this->center = cen;
	this->first_segment = tag::local_functions::iterator_h::find_first_seg ( cen, this->msh_p );
	assert ( this->first_segment );

	// now we rotate backwards until we meet the boundary or close the loop
	tag::local_functions::iterator_h::rotate_forward_for_seg ( this );
	assert ( this->first_segment );                                                              }
 	
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments::NormalOrder::AsTheyAre
: public Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments::NormalOrder

{	public :

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here

	// attributes inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments :
	// Cell::Core * first_segment { nullptr }  -- word "segment" is misleading here
	// Cell::Core * last_vertex { nullptr }  -- word "vertex" is misleading here

	inline AsTheyAre ( Mesh::Core * msh, Cell::Positive * const c )
	:	Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments::NormalOrder ( msh, c )
	// msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  -- virtual,
	//    defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments

	// void reset ( const tag::StartAt &, Cell::Core * cll )  -- virtual,
	//    defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments::NormalOrder

	Cell deref ( );  // virtual from Mesh::Iterator::Core

	// void advance ( )  -- virtual,
	//    defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments::NormalOrder
	
	// bool in_range ( )  virtual, defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments
	
};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments::NormalOrder::AsTheyAre

//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments::ReverseOrder::AsTheyAre
: public Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments::ReverseOrder

{	public :

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here

	// attributes inherited from Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments :
	// Cell::Core * first_segment { nullptr }  -- word "segment" is misleading here
	// Cell::Core * last_vertex { nullptr }  -- word "vertex" is misleading here

	inline AsTheyAre ( Mesh::Core * msh, Cell::Positive * const c )
	:	Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments::ReverseOrder ( msh, c )
	// msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  -- virtual,
	//    defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments

	// void reset ( const tag::StartAt &, Cell::Core * cll )  -- virtual,
	//    defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments::ReverseOrder

	Cell deref ( );  // virtual from Mesh::Iterator::Core

	// void advance ( )  -- virtual,
	//    defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments::ReverseOrder
	
	// bool in_range ( )  virtual, defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments
	
};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo::OverSegments::ReverseOrder::AsTheyAre

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfCodimTwo::WorkAround2D::NormalOrder::BuildReverseCells
: public Mesh::Iterator::Around::OneCell::OfCodimTwo
              ::OverVertices::NormalOrder::BuildReverseCells::AsTheyAre

// this is a quick and ugly workaround for a specific need
// we are in a 2D mesh and we want all vertices around a given vertex
// so we inherit from ::OverVertices:: - which produces segments !
// and just override 'deref' to return the reversed base of the segment

{	public :

	inline BuildReverseCells ( Mesh::Core * msh, Cell::Positive * const c )
	:	Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices
	    ::NormalOrder::BuildReverseCells::AsTheyAre ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	Cell deref ( ) override;  // virtual from Mesh::Iterator::Core,
	// defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder::
	// ::BuildReverseCells::AsTheyAre, here overriden to return the reversed base of the segment

};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo::WorkAround2D::NormalOrder::BuildReverseCells


class Mesh::Iterator::Around::OneCell::OfCodimTwo::WorkAround2D::ReverseOrder::BuildReverseCells
: public Mesh::Iterator::Around::OneCell::OfCodimTwo
              ::OverVertices::ReverseOrder::BuildReverseCells::AsTheyAre

// this is a quick and ugly workaround for a specific need
// we are in a 2D mesh and we want all vertices around a given vertex
// so we inherit from ::OverVertices:: - which produces segments !
// and just override 'deref' to return the reversed base of the segment

{	public :

	Cell deref ( ) override;  // virtual from Mesh::Iterator::Core,
	// defined by Mesh::Iterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder::
	// ::BuildReverseCells::AsTheyAre, here overriden to return the reversed base of the segment

};  // end of class Mesh::Iterator::Around::OneCell::OfCodimTwo::WorkAround2D::ReverseOrder::BuildReverseCells

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

class Mesh::Iterator::Around::OneCell::OfAnyCodim
: public Mesh::Iterator::Around::OneCell

// cll->dim <= msh->dim - 2, no order assumed
// more often than not, cll->dim < msh->dim - 2

{	public :

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	size_t dim;  // dimension of cells to be returned
	// center->dim < dim <= msh->dim

	inline OfAnyCodim ( Mesh::Core * msh, size_t d, Cell::Positive * const c )
	:	Mesh::Iterator::Around::OneCell ( msh, c ), dim { d }  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  --  virtual from Mesh::Iterator::Core
	
	// void reset ( const tag::StartAt &, Cell::Core * cll )  --  virtual from Mesh::Iterator::Core
	
	// Cell deref ( )  stays pure virtual from Mesh::Iterator::Core

	void advance ( );  // virtual from Mesh::Iterator::Core

	bool in_range ( );  // virtual from Mesh::Iterator::Core

	class OverVertices;  class OverHDCells;

};  // end of class Mesh::Iterator::Around::OneCell::OfAnyCodim

//-------------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneVertex::OfAnyCodim::OverSegments
: public Mesh::Iterator::Around::OneCell::OfAnyCodim

// center->dim  == 0, no order assumed
// more often than not, msh->dim >= 3 (otherwise, an ordered iterator could be used)

// here we no longer use the misleading terms insstroduced in
// class Mesh::Iterator::Around::OneCell::OfCodimTwo::Ordered
// so, "Vertex" is really a vertex, "Segments" are really segments

{	public :

	typedef std::map < Cell::Positive::Segment*, short int > maptype_segs;

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attributes inherited from Mesh::Iterator::Around::OneCell::OfAnyCodim :
	// size_t dim  -- dimension of cells to be returned
	// center->dim + 1 == dim <= msh->dim

	maptype_segs::iterator it;;
	
	inline OverSegments ( Mesh::Core * msh, Cell::Positive * const c )
	:	Mesh::Iterator::Around::OneCell::OfAnyCodim ( msh, c->get_dim() + 1, c )
	// msh_p { msh }, center { c }, dim { c->get_dim() + 1 }
	{ }
	
	// void reset ( )  --  virtual from Mesh::Iterator::Core
	
	// void reset ( const tag::StartAt &, Cell::Core * cll )  --  virtual from Mesh::Iterator::Core
	
	// Cell deref ( )  stays pure virtual from Mesh::Iterator::Core

	void advance ( );  // virtual from Mesh::Iterator::Core

	bool in_range ( );  // virtual from Mesh::Iterator::Core

	struct ReverseEachCell  {  class  AssumeCellsExist;  };

};  // end of class Mesh::Iterator::Around::OneVertex::OfAnyCodim::OverSegments

//-------------------------------------------------------------------------------------------


class Mesh::Iterator::Around::OneCell::OfAnyCodim::OverHDCells
: public Mesh::Iterator::Around::OneCell::OfAnyCodim

// center->dim <= msh->dim - 2, no order assumed
// more often than not, center->dim < msh->dim - 2

// we return cells of dimension >= center->dim + 2
// but not of maximum dimension (maximum dimension means the dimension of the mesh)

// example : squares arond a vertex in mesh of cubes

{	public :

	typedef std::map < Mesh::Core*, Cell::field_to_meshes > maptype;

	// attributes inherited from Mesh::Iterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attributes inherited from Mesh::Iterator::Around::OneCell::OfAnyCodim :
	// size_t dim  -- dimension of cells to be returned
	// center->dim + 2 <= dim < msh->dim

	maptype::iterator it;;
	
	inline OverHDCells ( Mesh::Core * msh, Cell::Positive * const c )
	:	Mesh::Iterator::Around::OneCell::OfAnyCodim ( msh, c->get_dim() + 1, c )
	// msh_p { msh }, center { c }, dim { c->get_dim() + 1 }
	{ }
	
	// void reset ( )  --  virtual from Mesh::Iterator::Core
	
	// void reset ( const tag::StartAt &, Cell::Core * cll )  --  virtual from Mesh::Iterator::Core
	
	// Cell deref ( )  stays pure virtual from Mesh::Iterator::Core

	void advance ( );  // virtual from Mesh::Iterator::Core

	bool in_range ( );  // virtual from Mesh::Iterator::Core

	struct ReverseEachCell  {  class  AssumeCellsExist;  };

};  // end of class Mesh::Iterator::Around::OneCell::OfAnyCodim::OverHDCells

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------


class Mesh::Iterator::Adaptor::ForcePositive : public Mesh::Iterator::Core

// modifies the behaviour of another Mesh::Iterator
// forcing the resulting cells to be positive

{	public :

	Mesh::Iterator base;

	inline ForcePositive ( Mesh::Iterator && b )
	:	base { b } { };

	void reset ( );  // virtual from Mesh::Iterator::Core
	void reset ( const tag::StartAt &, Cell::Core * cll );
	// virtual from Mesh::Iterator::Core, here execution forbidden
	// Cell deref ( )  stays pure virtual from Mesh::Iterator::Core
	void advance ( );  // virtual from Mesh::Iterator::Core
	bool in_range ( );  // virtual from Mesh::Iterator::Core

};  // end of class Mesh::Iterator::Adaptor::ForcePositive


//-------------------------------------------------------------------------------------------

//      ITERATORS OVER VERTICES
//--------------------------------------------------------------------------------
		

inline Mesh::Iterator Mesh::iterator ( const tag::OverVertices & ) const
{	return this->iterator ( tag::over_vertices, tag::as_they_are );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverVertices &, const tag::AsTheyAre & ) const
{	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::as_they_are, tag::this_mesh_is_positive ) );
	// else
	if ( this->dim() == 0 )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::reverse_each_cell, tag::do_not_build_cells,
			  tag::reverse_order, tag::this_mesh_is_positive                      ) );
	// else : dim >= 1, all vertices are positive, no need to reverse them
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_vertices, tag::as_they_are,
		  tag::reverse_order_if_any, tag::this_mesh_is_positive )     );           }


inline Mesh::Iterator Mesh::iterator
( const tag::OverVertices &, const tag::RequireOrder & ) const
{	return this->iterator ( tag::over_vertices, tag::as_they_are, tag::require_order );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverVertices &, const tag::AsTheyAre &, const tag::RequireOrder & ) const
{	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::as_they_are,
			  tag::require_order, tag::this_mesh_is_positive )           );
	// else : negative mesh
	if ( this->dim() == 0 )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::reverse_each_cell, tag::do_not_build_cells,
			  tag::reverse_order, tag::this_mesh_is_positive                      ) );
	// else : dim >= 1, all vertices are positive, no need to reverse them
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_vertices, tag::as_they_are,
		  tag::reverse_order, tag::this_mesh_is_positive )           );           }


inline Mesh::Iterator Mesh::iterator
( const tag::OverVertices &, const tag::RequireOrder &, const tag::AsTheyAre & ) const
{	return this->iterator ( tag::over_vertices, tag::as_they_are, tag::require_order );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverVertices &, const tag::Backwards & ) const
{	return this->iterator ( tag::over_vertices, tag::as_they_are, tag::backwards );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverVertices &, const tag::AsTheyAre &, const tag::Backwards & ) const
{	assert ( this->dim() <= 1 );  // no order for dim >= 2
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::as_they_are,
			  tag::reverse_order, tag::this_mesh_is_positive )            );
	// else : negative cell, back to normal order
	// but reverse cells if mesh is zero-dimensional
	if ( this->core->get_dim_plus_one() == 1 )  // dimension zero
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::reverse_each_cell,
			  tag::do_not_build_cells, tag::this_mesh_is_positive )      );
	// else : dimension one
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
	  ( tag::over_vertices, tag::as_they_are,
		  tag::require_order, tag::this_mesh_is_positive )           );    }


inline Mesh::Iterator Mesh::iterator
( const tag::OverVertices &, const tag::Backwards &, const tag::AsTheyAre & ) const
{	return this->iterator ( tag::over_vertices, tag::as_they_are, tag::backwards );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverVertices &, const tag::ForcePositive & ) const
{	if ( this->dim() == 0 )
	{	if ( this->is_positive() )
			return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
				( tag::over_vertices, tag::force_positive, tag::this_mesh_is_positive ) );
		// else : negative mesh
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::force_positive,
			  tag::reverse_order, tag::this_mesh_is_positive ) );                        }
	// else : dim >= 1, all vertices are positive
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::as_they_are, tag::this_mesh_is_positive ) );
	// else : negative mesh
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_vertices, tag::as_they_are,
		  tag::reverse_order_if_any, tag::this_mesh_is_positive ) );                     }


inline Mesh::Iterator Mesh::iterator
( const tag::OverVertices &, const tag::ForcePositive &, const tag::RequireOrder & ) const
{	assert ( this->dim() <= 1 );  // no order for dim >= 2
	if ( this->core->get_dim_plus_one() == 1 )  // dimension zero
	{	if ( this->is_positive() )
			return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
				( tag::over_vertices, tag::force_positive, tag::this_mesh_is_positive ) );
		// else : negative mesh
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::force_positive,
			  tag::reverse_order, tag::this_mesh_is_positive ) );                        }
	// else : dimension one, all vertices are positive
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::as_they_are,
			  tag::require_order, tag::this_mesh_is_positive )            );
	// else : negative mesh
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_vertices, tag::as_they_are,
		  tag::reverse_order, tag::this_mesh_is_positive )            );                  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverVertices &, const tag::RequireOrder &, const tag::ForcePositive & ) const
{	return this->iterator ( tag::over_vertices, tag::force_positive, tag::require_order );  }


inline Mesh::Iterator Mesh:: iterator
( const tag::OverVertices &, const tag::ForcePositive &, const tag::Backwards & ) const
{ assert ( this->dim() <= 1 );  // no order for dim >= 2
	if ( this->dim() == 0 )
	{	if ( this->is_positive() )
			return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
				( tag::over_vertices, tag::force_positive,
				  tag::reverse_order, tag::this_mesh_is_positive )           );
		// else : negative mesh
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::force_positive, tag::this_mesh_is_positive ) ); }
	// else : dim == 1, all vertices are positive
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::as_they_are,
			  tag::reverse_order, tag::this_mesh_is_positive )           );
	else
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::as_they_are,
			  tag::require_order, tag::this_mesh_is_positive )            );                 }


inline Mesh::Iterator Mesh::iterator
( const tag::OverVertices &, const tag::Backwards &, const tag::ForcePositive & ) const
{	return this->iterator ( tag::over_vertices, tag::force_positive, tag::backwards );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == 0 );  // because reverse_each_cell
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::reverse_each_cell,
			  tag::do_not_build_cells, tag::this_mesh_is_positive )      );
	// else : negative mesh, reverse order
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
	  ( tag::over_vertices, tag::as_they_are, tag::reverse_order, tag::this_mesh_is_positive ) );  }
	

inline Mesh::Iterator Mesh::iterator
( const tag::OverVertices &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::RequireOrder &      ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == 0 );  // because reverse_each_cell
	return this->iterator ( tag::over_vertices, tag::reverse_each_cell, tag::do_not_build_cells );  }
	

inline Mesh::Iterator Mesh::iterator
( const tag::OverVertices &, const tag::RequireOrder &,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == 0 );  // because reverse_each_cell
	return this->iterator ( tag::over_vertices, tag::reverse_each_cell, tag::do_not_build_cells );  }
	

inline Mesh::Iterator Mesh::iterator
( const tag::OverVertices &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Backwards &        ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == 0 );  // because reverse_each_cell
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::reverse_each_cell, tag::do_not_build_cells,
			  tag::reverse_order, tag::this_mesh_is_positive                      ) );
	// else : negative mesh, back to normal order
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
	  ( tag::over_vertices, tag::as_they_are, tag::this_mesh_is_positive ) );   }


inline Mesh::Iterator Mesh::iterator
( const tag::OverVertices &, const tag::Backwards &,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_vertices, tag::reverse_each_cell,
                          tag::do_not_build_cells, tag::backwards    );  }

//-------------------------------------------------------------------------------------

//      ITERATORS OVER SEGMENTS
//--------------------------------------------------------------------------------
		

inline Mesh::Iterator Mesh::iterator ( const tag::OverSegments & ) const
{	return this->iterator ( tag::over_segments, tag::as_they_are );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverSegments &, const tag::AsTheyAre & ) const
{	assert ( this->dim() >= 1 );
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_segments, tag::as_they_are, tag::this_mesh_is_positive ) );
	// else
	if ( this->core->get_dim_plus_one() == 2 )  // dimension one
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_segments, tag::reverse_each_cell, tag::do_not_build_cells,
			  tag::reverse_order_if_any, tag::this_mesh_is_positive               ) );
	// else : dim >= 2, all segments are positive, no need to reverse them
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_segments, tag::as_they_are, tag::this_mesh_is_positive ) );   }


inline Mesh::Iterator Mesh::iterator
( const tag::OverSegments &, const tag::RequireOrder & ) const
{	return this->iterator ( tag::over_segments, tag::as_they_are, tag::require_order );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::RequireOrder & ) const
{	assert ( this->dim() == 1 );  // no order for dim >= 2
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		  ( tag::over_segments, tag::as_they_are, tag::require_order,
			  tag::this_mesh_is_positive                                ) );
	// else
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_segments, tag::reverse_each_cell, tag::do_not_build_cells,
		  tag::reverse_order, tag::this_mesh_is_positive                      ) ); }


inline Mesh::Iterator Mesh::iterator
( const tag::OverSegments &, const tag::RequireOrder &, const tag::AsTheyAre & ) const
{	return this->iterator ( tag::over_segments, tag::as_they_are, tag::require_order );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverSegments &, const tag::Backwards & ) const
{	return this->iterator ( tag::over_segments, tag::as_they_are, tag::backwards );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::Backwards & ) const
{	assert ( this->dim() == 1 );  // no order for dim >= 2
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_segments, tag::as_they_are,
			  tag::reverse_order, tag::this_mesh_is_positive )            );
	// else
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_segments, tag::reverse_each_cell, tag::do_not_build_cells,
		  tag::require_order, tag::this_mesh_is_positive                      ) ); }


inline Mesh::Iterator Mesh::iterator
( const tag::OverSegments &, const tag::Backwards &, const tag::AsTheyAre & ) const
{	return this->iterator ( tag::over_segments, tag::as_they_are, tag::backwards );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverSegments &, const tag::ForcePositive & ) const
{	assert ( this->dim() >= 1 );
	if ( this->core->get_dim_plus_one() == 2 )  // dimension one
	{	if ( this->is_positive() )
			return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
				( tag::over_segments, tag::force_positive, tag::this_mesh_is_positive ) );
		// else : negative mesh
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_segments, tag::force_positive,
				tag::reverse_order_if_any, tag::this_mesh_is_positive ) );                 }
	// else : dim >= 2, all segments are positive
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_segments, tag::as_they_are, tag::this_mesh_is_positive ) );
	// else : negative mesh
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_segments, tag::as_they_are,
		  tag::reverse_order_if_any, tag::this_mesh_is_positive )     );                  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverSegments &, const tag::ForcePositive &, const tag::RequireOrder & ) const
{	assert ( this->dim() == 1 );  // no order for dim >= 2
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_segments, tag::force_positive,
			  tag::require_order, tag::this_mesh_is_positive ) );
	else
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_segments, tag::force_positive,
			  tag::reverse_order, tag::this_mesh_is_positive )            );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverSegments &, const tag::RequireOrder &, const tag::ForcePositive & ) const
{	return this->iterator ( tag::over_segments, tag::force_positive, tag::require_order );  }


inline Mesh::Iterator Mesh:: iterator
( const tag::OverSegments &, const tag::ForcePositive &, const tag::Backwards & ) const
{ assert ( this->dim() == 1 );  // no order for dim >= 2
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_segments, tag::force_positive,
			  tag::reverse_order, tag::this_mesh_is_positive ) );
	// else : negative mesh
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_segments, tag::force_positive,
		  tag::require_order, tag::this_mesh_is_positive ) );           }


inline Mesh::Iterator Mesh::iterator
( const tag::OverSegments &, const tag::Backwards &, const tag::ForcePositive & ) const
{	return this->iterator ( tag::over_segments, tag::force_positive, tag::backwards );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == 1 );  // because reverse_each_cell
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_segments, tag::reverse_each_cell,
			  tag::do_not_build_cells, tag::this_mesh_is_positive )      );
	// else : negative mesh, reverse order
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_segments, tag::as_they_are,
		  tag::reverse_order_if_any, tag::this_mesh_is_positive )    );   }
	

inline Mesh::Iterator Mesh::iterator
( const tag::OverSegments &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::RequireOrder &      ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == 1 );  // because reverse_each_cell
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_segments, tag::reverse_each_cell, tag::do_not_build_cells,
			  tag::require_order, tag::this_mesh_is_positive                      ) );
	// else : negative mesh, reverse order
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_segments, tag::as_they_are,
		  tag::reverse_order, tag::this_mesh_is_positive )           );              }
	

inline Mesh::Iterator Mesh::iterator
( const tag::OverSegments &, const tag::RequireOrder &,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_segments, tag::reverse_each_cell,
                          tag::do_not_build_cells, tag::require_order );  }
	

inline Mesh::Iterator Mesh::iterator
( const tag::OverSegments &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Backwards &      ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == 1 );  // because reverse_each_cell
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_segments, tag::reverse_each_cell, tag::do_not_build_cells,
			  tag::reverse_order, tag::this_mesh_is_positive                      ) );
	// else : negative mesh, back to normal order
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_segments, tag::as_they_are,
		  tag::require_order, tag::this_mesh_is_positive )           );              }


inline Mesh::Iterator Mesh::iterator
( const tag::OverSegments &, const tag::Backwards &,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_segments, tag::reverse_each_cell,
                          tag::do_not_build_cells, tag::backwards    );  }


//-------------------------------------------------------------------------------------

//      ITERATORS OVER HIGH DIMENSIONAL CELLS
//--------------------------------------------------------------------------------
		

inline Mesh::Iterator Mesh::iterator ( const tag::OverCellsOfDim &, const size_t d ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre & ) const
{	assert ( this->dim() >= d );
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_dim, d, tag::as_they_are, tag::this_mesh_is_positive ) );
	// else : negative mesh
	if ( this->dim() == d )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::reverse_each_cell, tag::do_not_build_cells,
			  tag::reverse_order_if_any, tag::this_mesh_is_positive                       ) );
	// else : dim > d, all cells are positive, no need to reverse them
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_dim, d, tag::as_they_are,
		  tag::reverse_order_if_any, tag::this_mesh_is_positive )      );                    }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::RequireOrder & ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are, tag::require_order );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d,
  const tag::AsTheyAre &, const tag::RequireOrder & ) const
{	assert ( this->dim() <= 1 );  // no order for dim >= 2
	assert ( this->dim() >= d );
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		  ( tag::over_cells_of_dim, d, tag::as_they_are,
			  tag::require_order, tag::this_mesh_is_positive )           );
	// else : negative mesh
	if ( this->dim() == d )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::reverse_each_cell, tag::do_not_build_cells,
			  tag::reverse_order, tag::this_mesh_is_positive                              ) );
	// else : dim == 1 > d == 0, all cells are positive, no need to reverse them
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_vertices, tag::as_they_are,
		  tag::reverse_order, tag::this_mesh_is_positive )           );                       }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d,
  const tag::RequireOrder &, const tag::AsTheyAre & ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are, tag::require_order );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::Backwards & ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are, tag::backwards );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d,
  const tag::AsTheyAre &, const tag::Backwards & ) const
{	assert ( this->dim() <= 1 );  // no order for dim >= 2
	assert ( this->dim() >= d );
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_dim, d, tag::as_they_are,
			  tag::reverse_order, tag::this_mesh_is_positive )            );
	// else : negative mesh
	if ( this->dim() == d )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::reverse_each_cell, tag::do_not_build_cells,
			  tag::require_order, tag::this_mesh_is_positive                              ) );
	// else : dim == 1 > d == 0, all cells are positive, no need to reverse them
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_vertices, tag::as_they_are,
		  tag::require_order, tag::this_mesh_is_positive )           );                       }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d,
  const tag::Backwards &, const tag::AsTheyAre & ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are, tag::backwards );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive & ) const
{	assert ( this->dim() >= d );
	if ( this->dim() == d )
	{	if ( this->is_positive() )
			return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
				( tag::over_cells_of_max_dim, tag::force_positive, tag::this_mesh_is_positive ) );
		// else : negative mesh
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::force_positive,
			  tag::reverse_order_if_any, tag::this_mesh_is_positive )    );                      }
	// else : dim > d, all cells are positive
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_dim, d, tag::as_they_are, tag::this_mesh_is_positive ) );
	// else : negative mesh
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_dim, d, tag::as_they_are,
		  tag::reverse_order_if_any, tag::this_mesh_is_positive )     );                          }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d,
  const tag::ForcePositive &, const tag::RequireOrder & ) const
{	assert ( this->dim() <= 1 );  // no order for dim >= 2
	assert ( this->dim() >= d );
	if ( this->dim() == d )
	{	if ( this->is_positive() )
			return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
				( tag::over_cells_of_max_dim, tag::force_positive,
				  tag::require_order, tag::this_mesh_is_positive   )         );
		// else : negative mesh
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::force_positive,
			  tag::reverse_order, tag::this_mesh_is_positive   )          );  }
	// else : dim == 1 > d == 0, all cells are positive
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::as_they_are,
			  tag::require_order, tag::this_mesh_is_positive )           );
	// else : negative mesh
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_vertices, tag::as_they_are,
		  tag::reverse_order, tag::this_mesh_is_positive )           );       }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d,
  const tag::RequireOrder &, const tag::ForcePositive & ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::force_positive, tag::require_order );  }


inline Mesh::Iterator Mesh:: iterator
( const tag::OverCellsOfDim &, const size_t d,
  const tag::ForcePositive &, const tag::Backwards & ) const
{ assert ( this->dim() <= 1 );  // no order for dim >= 2
	assert ( this->dim() >= d );
	if ( this->dim() == d )
	{	if ( this->is_positive() )
			return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
				( tag::over_cells_of_max_dim, tag::force_positive,
				  tag::reverse_order, tag::this_mesh_is_positive   )         );
		// else : negative mesh
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::force_positive,
			  tag::require_order, tag::this_mesh_is_positive   )          );  }
	// else : dim == 1 > d == 0, all cells are positive
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::as_they_are,
			  tag::reverse_order, tag::this_mesh_is_positive )           );
	// else : negative mesh
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_vertices, tag::as_they_are,
		  tag::require_order, tag::this_mesh_is_positive )           );       }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d,
  const tag::Backwards &, const tag::ForcePositive & ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::force_positive, tag::backwards );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == d );  // because reverse_each_cell
	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
                          tag::do_not_build_cells                             ); }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::RequireOrder &                   ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == d );  // because reverse_each_cell
	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
	                        tag::do_not_build_cells, tag::require_order        );  }
	

inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::RequireOrder &,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells &             ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == d );  // because reverse_each_cell
	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
	                        tag::do_not_build_cells, tag::require_order          );  }
	

inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Backwards &                       ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == d );  // because reverse_each_cell
	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
	                        tag::do_not_build_cells, tag::backwards            );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::Backwards &,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells &          ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == d );  // because reverse_each_cell
	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
	                        tag::do_not_build_cells, tag::backwards            );  }

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------


//      ITERATORS OVER CELLS OF MAXIMUM DIMENSION
//------------------------------------------------------------------------------------
		

inline Mesh::Iterator Mesh::iterator ( const tag::OverCellsOfMaxDim & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::as_they_are );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre & ) const
{	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::as_they_are, tag::this_mesh_is_positive ) );
	// else : negative mesh
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_max_dim, tag::reverse_each_cell, tag::do_not_build_cells,
		  tag::reverse_order_if_any, tag::this_mesh_is_positive                       ) );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::RequireOrder & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::as_they_are, tag::require_order );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::RequireOrder & ) const
{	assert ( this->dim() <= 1 );  // no order for dim >= 2
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		  ( tag::over_cells_of_max_dim, tag::as_they_are, tag::require_order,
			  tag::this_mesh_is_positive                                        ) );
	// else : negative mesh
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_max_dim, tag::reverse_each_cell, tag::do_not_build_cells,
		  tag::reverse_order, tag::this_mesh_is_positive                              ) );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::RequireOrder &, const tag::AsTheyAre & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::as_they_are, tag::require_order );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::Backwards & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::as_they_are, tag::backwards );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::Backwards & ) const
{	assert ( this->dim() <= 1 );  // no order for dim >= 2
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::as_they_are,
			  tag::reverse_order, tag::this_mesh_is_positive )            );
	// else : negative mesh
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_max_dim, tag::reverse_each_cell, tag::do_not_build_cells,
		  tag::require_order, tag::this_mesh_is_positive                              ) );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::Backwards &, const tag::AsTheyAre & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::as_they_are, tag::backwards );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive & ) const
{	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::force_positive, tag::this_mesh_is_positive ) );
	// else : negative mesh
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_max_dim, tag::force_positive,
		  tag::reverse_order_if_any, tag::this_mesh_is_positive )    );                     }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &, const tag::RequireOrder & ) const
{	assert ( this->dim() <= 1 );  // no order for dim >= 2
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::force_positive,
			  tag::require_order, tag::this_mesh_is_positive  )          );
	// else : negative mesh, reverse order
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_max_dim, tag::force_positive,
		  tag::reverse_order, tag::this_mesh_is_positive  )          );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::RequireOrder &, const tag::ForcePositive & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::force_positive,
                          tag::require_order                               );  }


inline Mesh::Iterator Mesh:: iterator
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &, const tag::Backwards & ) const
{ assert ( this->dim() <= 1 );  // no order for dim >= 2
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::force_positive,
			  tag::reverse_order, tag::this_mesh_is_positive   )         );
	// else : negative mesh
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_max_dim, tag::force_positive,
		  tag::require_order, tag::this_mesh_is_positive   )          );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::Backwards &, const tag::ForcePositive & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::force_positive, tag::backwards );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
// do not bother whether reverse cells exist or not	
{	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::reverse_each_cell,
			  tag::do_not_build_cells, tag::this_mesh_is_positive )       );
	// else : negative mesh, reverse order
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_max_dim, tag::as_they_are,
		  tag::reverse_order_if_any, tag::this_mesh_is_positive )    );   }
	

inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::RequireOrder &      ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() <= 1 );  // no order for dim >= 2
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::reverse_each_cell, tag::do_not_build_cells,
			  tag::require_order, tag::this_mesh_is_positive                              ) );
	// else : negative mesh, reverse order
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_max_dim, tag::as_they_are,
		  tag::reverse_order, tag::this_mesh_is_positive )           );                     }
	

inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::RequireOrder &,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
                          tag::do_not_build_cells, tag::require_order        );  }
	

inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Backwards &         ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == 1 );  // no order for dim >= 2
	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::reverse_each_cell, tag::do_not_build_cells,
			  tag::reverse_order, tag::this_mesh_is_positive                              ) );
	// else : negative mesh, back to normal order
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_max_dim, tag::as_they_are,
		  tag::require_order, tag::this_mesh_is_positive )           );                     }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::Backwards &,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
                          tag::do_not_build_cells, tag::backwards             );  }
	
//-----------------------------------------------------------------------------//


inline Mesh::Iterator Mesh::iterator ( const tag::OverCells &, const tag::OfMaxDim & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::as_they_are );  }

inline Mesh::Iterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &, const tag::AsTheyAre & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::as_they_are );  }

inline Mesh::Iterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &, const tag::RequireOrder & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::as_they_are, tag::require_order );  }

inline Mesh::Iterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &,
  const tag::AsTheyAre &, const tag::RequireOrder & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::as_they_are, tag::require_order );  }

inline Mesh::Iterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &,
  const tag::RequireOrder &, const tag::AsTheyAre & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::as_they_are, tag::require_order );  }

inline Mesh::Iterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &, const tag::Backwards & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::as_they_are, tag::backwards );  }

inline Mesh::Iterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &,
  const tag::AsTheyAre &, const tag::Backwards & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::as_they_are, tag::backwards );  }

inline Mesh::Iterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &,
  const tag::Backwards &, const tag::AsTheyAre & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::as_they_are, tag::backwards );  }

inline Mesh::Iterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &, const tag::ForcePositive & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::force_positive );  }

inline Mesh::Iterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &,
  const tag::ForcePositive &, const tag::RequireOrder & ) const
{	return this->iterator ( tag::over_cells_of_max_dim,
                          tag::force_positive, tag::require_order );  }

inline Mesh::Iterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &,
  const tag::RequireOrder &, const tag::ForcePositive & ) const
{	return this->iterator ( tag::over_cells_of_max_dim,
                          tag::force_positive, tag::require_order );  }

inline Mesh::Iterator Mesh:: iterator
( const tag::OverCells &, const tag::OfMaxDim &,
  const tag::ForcePositive &, const tag::Backwards & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::force_positive, tag::backwards );  }

inline Mesh::Iterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &,
  const tag::Backwards &, const tag::ForcePositive & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::force_positive, tag::backwards );  }

inline Mesh::Iterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_max_dim,
                          tag::reverse_each_cell, tag::do_not_build_cells );  }
	
inline Mesh::Iterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::RequireOrder &                     ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
                          tag::do_not_build_cells, tag::require_order        );  }
	
inline Mesh::Iterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &, const tag::RequireOrder &,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells &               ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
                          tag::do_not_build_cells, tag::require_order        );  }
	
inline Mesh::Iterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Backwards &                         ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
                          tag::do_not_build_cells, tag::backwards             );  }

inline Mesh::Iterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &, const tag::Backwards &,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells &            ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
                          tag::do_not_build_cells, tag::backwards             );  }
	
//-----------------------------------------------------------------------------//

//      ITERATORS OVER SEGMENTS AROUND A GIVEN CELL
//------------------------------------------------------------------------------------


inline Mesh::Iterator Mesh::iterator
( const tag::OverVertices &, const tag::Around &, const Cell & c ) const
{	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
      ( tag::over_vertices, tag::as_they_are, tag::build_cells_if_necessary,
        tag::around, c.core, tag::this_mesh_is_positive                     ) );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverSegments &, const tag::Around &, const Cell & c ) const
{	return this->iterator ( tag::over_segments, tag::as_they_are,
                          tag::build_cells_if_necessary, tag::around, c );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::Around &, const Cell & c ) const
{	return this->iterator ( tag::over_segments, tag::as_they_are,
                          tag::build_cells_if_necessary, tag::around, c );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
  const tag::Around &, const Cell & cll                                           ) const
{	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		  ( tag::over_segments, tag::as_they_are, tag::build_cells_if_necessary,
	      tag::around, cll.core, tag::this_mesh_is_positive                   ) );
	// else : 'this' mesh is negative
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_segments, tag::reverse_each_cell, tag::build_cells_if_necessary,
		  tag::around, cll.core, tag::reverse_order_if_any, tag::this_mesh_is_positive ) );  }
// or, equivalently :
//	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
//		( tag::over_segments, tag::as_they_are, tag::build_cells_if_necessary,
//		  tag::reverse_order_if_any, tag::around, cll.core->reverse_attr.core,
//		  tag::this_mesh_is_positive                                           ) );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverSegments &, const tag::DoNotBuildCells &,
  const tag::Around &, const Cell & c                     ) const
{	return this->iterator ( tag::over_segments, tag::as_they_are,
                          tag::do_not_build_cells, tag::around, c );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
  const tag::Around &, const Cell & cll                                           ) const
{	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		  ( tag::over_segments, tag::as_they_are, tag::build_cells_if_necessary,
	      tag::around, cll.core, tag::this_mesh_is_positive                   ) );
	// else : 'this' mesh is negative
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_segments, tag::reverse_each_cell, tag::do_not_build_cells,
		  tag::around, cll.core, tag::reverse_order_if_any, tag::this_mesh_is_positive ) );  }
// or, equivalently :
//	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
//		( tag::over_segments, tag::as_they_are, tag::do_not_build_cells,
//		  tag::reverse_order_if_any, tag::around, cll.core->reverse_attr.core,
//		  tag::this_mesh_is_positive                                           ) );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::Around &, const Cell & c                              ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are,
                          tag::build_cells_if_necessary, tag::around, c );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::AsTheyAre &, const tag::Around &, const Cell & c      ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are,
                          tag::build_cells_if_necessary, tag::around, c );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
	const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
  const tag::Around &, const Cell & c                              ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are,
                          tag::build_cells_if_necessary, tag::around, c );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::DoNotBuildCells &, const tag::Around &, const Cell & c ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are,
                          tag::do_not_build_cells, tag::around, c     );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
	const tag::AsTheyAre &, const tag::DoNotBuildCells &,
  const tag::Around &, const Cell & c                              ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are,
                          tag::do_not_build_cells, tag::around, c     );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::Around &, const Cell & c ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are,
                          tag::build_cells_if_necessary, tag::around, c );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::Around &, const Cell & c                                 ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are,
                          tag::build_cells_if_necessary, tag::around, c );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d,
	const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
  const tag::Around &, const Cell & cll                ) const
{	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_dim, d, tag::as_they_are, tag::build_cells_if_necessary,
			  tag::around, cll.core, tag::this_mesh_is_positive                          ) );
	// else : 'this' mesh is negative
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_dim, d, tag::reverse_each_cell, tag::build_cells_if_necessary,
		  tag::around, cll.core, tag::this_mesh_is_positive                                ) );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::DoNotBuildCells &,
  const tag::Around &, const Cell & c                                       ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are,
                          tag::do_not_build_cells, tag::around, c     );  }


inline Mesh::Iterator Mesh::iterator
( const tag::OverCellsOfDim, const size_t d,
  const tag::AsTheyAre &, const tag::DoNotBuildCells &,
  const tag::Around &, const Cell & cll                ) const
{	if ( this->is_positive() )
		return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_dim, d, tag::as_they_are, tag::do_not_build_cells,
			  tag::around, cll.core, tag::this_mesh_is_positive                    ) );
	// else : 'this' mesh is negative
	return Mesh::Iterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_dim, d, tag::reverse_each_cell, tag::do_not_build_cells,
		  tag::around, cll.core, tag::this_mesh_is_positive                          ) );  }



}  // namespace maniFEM

#endif  // ifndef MANIFEM_ITERATOR_H
