
// iterator.h 2021.07.06

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

// class CellIterator defined in mesh.h

//-----------------------------------------------------------------------------------------

class CellIterator::Core

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
	
};  // end of class CellIterator::Core


inline CellIterator::CellIterator ( const CellIterator & it)
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "ManiFEM doesn't know (yet) how to copy construct an iterator. Sorry." << std::endl;
	exit ( 1 );                                                                                }
// define virtual self-replicator for core iterators !!
	

// inline CellIterator::CellIterator ( CellIterator && it)
// {	this->core = std::move ( it.core );  }
	

inline CellIterator & CellIterator::operator= ( const CellIterator & it )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "ManiFEM doesn't know (yet) how to copy an iterator. Sorry." << std::endl;
	exit ( 1 );                                                                                   }
// define virtual self-replicator for core iterators !!


// inline CellIterator & CellIterator::operator= ( CellIterator && it )
// {	this->core = std::move ( it.core );  return *this;  }
	

inline void CellIterator::reset ( )  {  this->core->reset();  }


inline void CellIterator::reset ( const tag::StartAt &, Cell & cll )
{  this->core->reset( tag::start_at, cll.core );  }
	

inline Cell CellIterator::operator* ( )
{	return this->core->deref();  }


inline CellIterator & CellIterator::operator++ ( )  {  return this->advance();  }


inline CellIterator & CellIterator::operator++ ( int )  {  return this->advance();  }


inline CellIterator & CellIterator::advance ( )
{ this->core->advance();  return *this;  }


inline bool CellIterator::in_range ( )  {  return this->core->in_range();  }
		
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class CellIterator::Over::TwoVerticesOfSeg : public CellIterator::Core

{	public :

	Cell::PositiveSegment * seg_p;  // should an iterator keep the segment alive ?
	short unsigned int passage;

	inline TwoVerticesOfSeg ( Cell::PositiveSegment * seg )
	:	CellIterator::Core (), seg_p { seg } { };
	
	void reset ( );  // virtual from CellIterator::Core

	void reset ( const tag::StartAt &, Cell::Core * cll );
	// virtual from CellIterator::Core, here execution forbidden

	// Cell deref ( )  stays pure virtual from CellIterator::Core
	void advance ( );  // virtual from CellIterator::Core
	bool in_range ( );  // virtual from CellIterator::Core

	struct NormalOrder  {  class AsTheyAre;  class ForcePositive;
		struct ReverseEachCell  {  class  AssumeCellsExist;  class BuildCells;  };    };
	struct ReverseOrder  {  class AsTheyAre;  class ForcePositive;
		struct ReverseEachCell  {  class  AssumeCellsExist;  class BuildCells;  };    };
	
};  // end of class CellIterator::Over::TwoVerticesOfSeg

//-----------------------------------------------------------------------------------------


class CellIterator::Over::TwoVerticesOfSeg::NormalOrder::AsTheyAre
: public CellIterator::Over::TwoVerticesOfSeg
// iterate over the two vertices, first base (negative) then tip (positive)

{	public :

	// inherited from CellIterator::Over::TwoVerticesOfSeg :
	// Cell::PositiveSegment * seg_p;
	// short unsigned int passage;

	inline AsTheyAre ( Cell::PositiveSegment * seg )
	:	CellIterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by CellIterator::Over::TwoVerticesOfSeg, execution forbidden

	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg

};  // end of class CellIterator::Over::TwoVerticesOfSeg::NormalOrder::AsTheyAre

//-----------------------------------------------------------------------------------------


class CellIterator::Over::TwoVerticesOfSeg::NormalOrder::ForcePositive
: public CellIterator::Over::TwoVerticesOfSeg
// iterate over the two vertices, first base then tip (both positive)

{	public :

	// inherited from CellIterator::Over::TwoVerticesOfSeg :
	// Cell::PositiveSegment * seg_p;
	// short unsigned int passage;

	inline ForcePositive ( Cell::PositiveSegment * seg )
	:	CellIterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by CellIterator::Over::TwoVerticesOfSeg, execution forbidden

	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg

};  // end of class CellIterator::Over::TwoVerticesOfSeg::NormalOrder::ForcePositive

//-----------------------------------------------------------------------------------------


class CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::AsTheyAre
: public CellIterator::Over::TwoVerticesOfSeg
// iterate over the two vertices, first tip (positive) then base (negative)

{	public :

	// inherited from CellIterator::Over::TwoVerticesOfSeg :
	// Cell::PositiveSegment * seg_p;
	// short unsigned int passage;
	
	inline AsTheyAre ( Cell::PositiveSegment * seg )
	:	CellIterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by CellIterator::Over::TwoVerticesOfSeg, execution forbidden

	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg

};  // end of class CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::AsTheyAre

//-----------------------------------------------------------------------------------------


class CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::ForcePositive
: public CellIterator::Over::TwoVerticesOfSeg
// iterate over the two vertices, first tip then base (both positive)

{	public :

	// inherited from CellIterator::Over::TwoVerticesOfSeg :
	// Cell::PositiveSegment * seg_p;
	// short unsigned int passage;

	inline ForcePositive ( Cell::PositiveSegment * seg )
	:	CellIterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by CellIterator::Over::TwoVerticesOfSeg, execution forbidden

	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg

};  // end of class CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::ForcePositive

//-----------------------------------------------------------------------------------------


class CellIterator::Over::TwoVerticesOfSeg::NormalOrder::ReverseEachCell::AssumeCellsExist
: public CellIterator::Over::TwoVerticesOfSeg
// iterate over the two vertices, first base (positive) then tip (negative)

{	public :

	// inherited from CellIterator::Over::TwoVerticesOfSeg :
	// Cell::PositiveSegment * seg_p;
	// short unsigned int passage;

	inline AssumeCellsExist ( Cell::PositiveSegment * seg )
	:	CellIterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by CellIterator::Over::TwoVerticesOfSeg, execution forbidden

	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg

};  // end of class CellIterator::Over::TwoVerticesOfSeg
    //         ::NormalOrder::ReverseEachCell::AssumeCellsExist

//-----------------------------------------------------------------------------------------


class CellIterator::Over::TwoVerticesOfSeg::NormalOrder::ReverseEachCell::BuildCells
: public CellIterator::Over::TwoVerticesOfSeg
// iterate over the two vertices, first base (positive) then tip (negative)

{	public :

	// inherited from CellIterator::Over::TwoVerticesOfSeg :
	// Cell::PositiveSegment * seg_p;
	// short unsigned int passage;

	inline BuildCells ( Cell::PositiveSegment * seg )
	:	CellIterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by CellIterator::Over::TwoVerticesOfSeg, execution forbidden

	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg

};  // end of class CellIterator::Over::TwoVerticesOfSeg
    //         ::NormalOrder::ReverseEachCell::BuildCells

//-----------------------------------------------------------------------------------------


class CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::ReverseEachCell::AssumeCellsExist
: public CellIterator::Over::TwoVerticesOfSeg
// iterate over the two vertices, first tip (negative) then base (positive)

{	public :

	// inherited from CellIterator::Over::TwoVerticesOfSeg :
	// Cell::PositiveSegment * seg_p;
	// short unsigned int passage;

	inline AssumeCellsExist ( Cell::PositiveSegment * seg )
	:	CellIterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by CellIterator::Over::TwoVerticesOfSeg, execution forbidden

	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg

};  // end of class CellIterator::Over::TwoVerticesOfSeg
    //         ::ReverseOrder::ReverseEachCell::AssumeCellsExist

//-----------------------------------------------------------------------------------------


class CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::ReverseEachCell::BuildCells
: public CellIterator::Over::TwoVerticesOfSeg
// iterate over the two vertices, first tip (negative) then base (positive)

{	public :

	// inherited from CellIterator::Over::TwoVerticesOfSeg :
	// Cell::PositiveSegment * seg_p;
	// short unsigned int passage;

	inline BuildCells ( Cell::PositiveSegment * seg )
	:	CellIterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by CellIterator::Over::TwoVerticesOfSeg, execution forbidden

	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg

};  // end of class CellIterator::Over::TwoVerticesOfSeg
    //         ::ReverseOrder::ReverseEachCell::BuildCells

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class CellIterator::Over::CellsOfConnectedOneDimMesh : public CellIterator::Core

{	public :

	Mesh::Connected::OneDim * msh;  // should an iterator keep the mesh alive ?
	Cell::Core * last_vertex;
	// positive vertex for iterator over vertices
	// positive vertex for iterator over segments, forward
	// negative vertex for iterator over segments, back

	inline CellsOfConnectedOneDimMesh ( Mesh::Connected::OneDim * m )
	:	CellIterator::Core (), msh { m } { }
	
	// void reset ( )  stays pure virtual from CellIterator::Core

	// void reset ( tag::StartAt , Cell::Core * cll )  stays pure virtual from CellIterator::Core

	// Cell deref ( )  stays pure virtual from CellIterator::Core
	// void advance ( )  stays pure virtual from CellIterator::Core
	// bool in_range ( )  stays pure virtual from CellIterator::Core

};  // end of class CellIterator::Over::CellsOfConnectedOneDimMesh

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class CellIterator::Over::VerticesOfConnectedOneDimMesh :
public CellIterator::Over::CellsOfConnectedOneDimMesh

{	public :

	// attributes inherited from CellIterator::Over::CellsOfConnectedOneDimMesh :
	// Mesh::Connected::OneDim * msh
	// Cell::Core * last_vertex (here positive vertex)

	Cell::Positive::Vertex * current_vertex;

	inline VerticesOfConnectedOneDimMesh ( Mesh::Connected::OneDim * m )
	:	CellIterator::Over::CellsOfConnectedOneDimMesh ( m ) { }
	
	// void reset ( )  stays pure virtual from CellIterator::Core

	// void reset ( tag::StartAt, Cell::Core * cll )  stays pure virtual from CellIterator::Core

	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  stays pure virtual from CellIterator::Core
	bool in_range ( );  // virtual from CellIterator::Core

	class NormalOrder;  class ReverseOrder;
	
};  // end of class CellIterator::Over::VerticesOfConnectedOneDimMesh

//-----------------------------------------------------------------------------------------


class CellIterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder
: public CellIterator::Over::VerticesOfConnectedOneDimMesh

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (here positive vertex)
	// inherited from CellIterator::Over::CellsOfConnectedOneDimMesh

	// attribute Cell::Positive::Vertex * current_vertex
	// inherited from CellIterator::Over::VerticesOfConnectedOneDimMesh
	
	inline NormalOrder ( Mesh::Connected::OneDim * m )
	:	CellIterator::Over::VerticesOfConnectedOneDimMesh ( m ) { };
	
	void reset ( );  // virtual from CellIterator::Core

	void reset ( const tag::StartAt &, Cell::Core * cll );  // virtual from CellIterator::Core

	// Cell deref ( )  defined by CellIterator::Over::VerticesOfConnectedOneDimMesh
	void advance ( );  // virtual from CellIterator::Core
	// bool in_range ( )  virtual, defined by CellIterator::Over::VerticesOfConnectedOneDimMesh

};  // end of class CellIterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder

//-----------------------------------------------------------------------------------------


class CellIterator::Over::VerticesOfConnectedOneDimMesh::ReverseOrder
: public CellIterator::Over::VerticesOfConnectedOneDimMesh

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (here positive vertex)
	// inherited from CellIterator::Over::CellsOfConnectedOneDimMesh

	// attribute Cell::Positive::Vertex * current_vertex
	// inherited from CellIterator::Over::VerticesOfConnectedOneDimMesh
	
	inline ReverseOrder ( Mesh::Connected::OneDim * m )
	:	CellIterator::Over::VerticesOfConnectedOneDimMesh ( m ) { };
	
	void reset ( );  // virtual from CellIterator::Core

	void reset ( const tag::StartAt &, Cell::Core * cll );  // virtual from CellIterator::Core

	// Cell deref ( )  defined by CellIterator::Over::VerticesOfConnectedOneDimMesh
	void advance ( );  // virtual from CellIterator::Core
	// bool in_range ( )  virtual, defined by CellIterator::Over::VerticesOfConnectedOneDimMesh

};  // end of class CellIterator::Over::VerticesOfConnectedOneDimMesh::ReverseOrder

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class CellIterator::Over::SegmentsOfConnectedOneDimMesh :
public CellIterator::Over::CellsOfConnectedOneDimMesh

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (may be positive or negative)
	// inherited from CellIterator::Over::CellsOfConnectedOneDimMesh

	Cell::Core * current_segment;
	
	inline SegmentsOfConnectedOneDimMesh ( Mesh::Connected::OneDim * m )
	:	CellIterator::Over::CellsOfConnectedOneDimMesh ( m ) { }
	
	// void reset ( )  stays pure virtual from CellIterator::Core

	// void reset ( tag::StartAt, Cell::Core * cll )  stays pure virtual from CellIterator::Core

	// Cell deref ( )  stays pure virtual from CellIterator::Core
	// void advance ( )  stays pure virtual from CellIterator::Core
	bool in_range ( );  // virtual from CellIterator::Core

	class NormalOrder;  class ReverseOrder;
	
};  // end of class CellIterator::Over::SegmentsOfConnectedOneDimMesh

//-----------------------------------------------------------------------------------------


class CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder
: public CellIterator::Over::SegmentsOfConnectedOneDimMesh

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (here positive vertex)
	// inherited from CellIterator::Over::CellsOfConnectedOneDimMesh

	// attribute Cell::Core * current_segment
	// inherited from CellIterator::Over::SegmentsOfConnectedOneDimMesh

	inline NormalOrder ( Mesh::Connected::OneDim * m )
	:	CellIterator::Over::SegmentsOfConnectedOneDimMesh ( m ) { };
	
	void reset ( );  // virtual from CellIterator::Core

	void reset ( const tag::StartAt &, Cell::Core * cll );  // virtual from CellIterator::Core

	// Cell deref ( )  stays pure virtual from CellIterator::Core
	void advance ( );  // virtual from CellIterator::Core
	// bool in_range ( )  virtual, defined by CellIterator::Over::VerticesOfConnectedOneDimMesh

	class AsTheyAre;  class ForcePositive;
	struct ReverseEachCell  {  class  AssumeCellsExist;  class BuildCells;  };

};  // end of class CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder

//-----------------------------------------------------------------------------------------


class CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder::AsTheyAre
: public CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (here positive vertex)
	// inherited from CellIterator::Over::CellsOfConnectedOneDimMesh

	// attribute Cell::Core * current_segment
	// inherited from CellIterator::Over::SegmentsOfConnectedOneDimMesh

	inline AsTheyAre ( Mesh::Connected::OneDim * m )
	:	CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder ( m ) { };
	
	// void reset ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder

	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder
	// bool in_range ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh

};  // end of class CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder::AsTheyAre

//-----------------------------------------------------------------------------------------


class CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder::ForcePositive
: public CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (here positive vertex)
	// inherited from CellIterator::Over::CellsOfConnectedOneDimMesh

	// attribute Cell::Core * current_segment
	// inherited from CellIterator::Over::SegmentsOfConnectedOneDimMesh

	inline ForcePositive ( Mesh::Connected::OneDim * m )
	:	CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder ( m ) { };
	
	// void reset ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder

	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder
	// bool in_range ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh

};  // end of class CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder::ForcePositive

//-----------------------------------------------------------------------------------------


class CellIterator::Over::SegmentsOfConnectedOneDimMesh
          ::NormalOrder::ReverseEachCell::AssumeCellsExist
: public CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (here positive vertex)
	// inherited from CellIterator::Over::CellsOfConnectedOneDimMesh

	// attribute Cell::Core * current_segment
	// inherited from CellIterator::Over::SegmentsOfConnectedOneDimMesh

	inline AssumeCellsExist ( Mesh::Connected::OneDim * m )
	:	CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder ( m ) { };
	
	// void reset ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder

	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder
	// bool in_range ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh

};  // end of class CellIterator::Over::SegmentsOfConnectedOneDimMesh
    //         ::NormalOrder::ReverseEachCell::AssumeCellsExist

//-----------------------------------------------------------------------------------------


class CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder::ReverseEachCell::BuildCells
: public CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (here positive vertex)
	// inherited from CellIterator::Over::CellsOfConnectedOneDimMesh

	// attribute Cell::Core * current_segment
	// inherited from CellIterator::Over::SegmentsOfConnectedOneDimMesh

	inline BuildCells ( Mesh::Connected::OneDim * m )
	:	CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder ( m ) { };
	
	// void reset ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder

	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder
	// bool in_range ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh

};  // end of class CellIterator::Over::SegmentsOfConnectedOneDimMesh
    //         ::NormalOrder::BuildCells::AssumeCellsExist

//-----------------------------------------------------------------------------------------


class CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder
: public CellIterator::Over::SegmentsOfConnectedOneDimMesh

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (here negative vertex)
	// inherited from CellIterator::Over::CellsOfConnectedOneDimMesh

	// attribute Cell::Core * current_segment
	// inherited from CellIterator::Over::SegmentsOfConnectedOneDimMesh

	inline ReverseOrder ( Mesh::Connected::OneDim * m )
	:	CellIterator::Over::SegmentsOfConnectedOneDimMesh ( m ) { };
	
	void reset ( );  // virtual from CellIterator::Core

	void reset ( const tag::StartAt &, Cell::Core * cll );  // virtual from CellIterator::Core

	// Cell deref ( )  stays pure virtual from CellIterator::Core
	void advance ( );  // virtual from CellIterator::Core
	// bool in_range ( )  virtual, defined by CellIterator::Over::VerticesOfConnectedOneDimMesh

	class AsTheyAre;  class ForcePositive;
	struct ReverseEachCell  {  class  AssumeCellsExist;  class BuildCells;  };

};  // end of class CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder

//-----------------------------------------------------------------------------------------


class CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder::AsTheyAre
: public CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (here negative vertex)
	// inherited from CellIterator::Over::CellsOfConnectedOneDimMesh

	// attribute Cell::Core * current_segment
	// inherited from CellIterator::Over::SegmentsOfConnectedOneDimMesh
	
	inline AsTheyAre ( Mesh::Connected::OneDim * m )
	:	CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder ( m ) { };
	
	// void reset ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder

	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder
	// bool in_range ( )  virtual, defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh

	class ForcePositive;
	struct ReverseEachCell  {  class AssumeCellsExist;  };

};  // end of class CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder::AsTheyAre

//-----------------------------------------------------------------------------------------


class CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder::ForcePositive
: public CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (here negative vertex)
	// inherited from CellIterator::Over::CellsOfConnectedOneDimMesh

	// attribute Cell::Core * current_segment
	// inherited from CellIterator::Over::SegmentsOfConnectedOneDimMesh

	inline ForcePositive ( Mesh::Connected::OneDim * m )
	:	CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder ( m ) { };
	
	// void reset ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder

	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder
	// bool in_range ( )  virtual, defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh

};  // end of class CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder::ForcePositive

//-----------------------------------------------------------------------------------------


class CellIterator::Over::SegmentsOfConnectedOneDimMesh
          ::ReverseOrder::ReverseEachCell::AssumeCellsExist
: public CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (here negative vertex)
	// inherited from CellIterator::Over::CellsOfConnectedOneDimMesh

	// attribute Cell::Core * current_segment
	// inherited from CellIterator::Over::SegmentsOfConnectedOneDimMesh

	inline AssumeCellsExist ( Mesh::Connected::OneDim * m )
	:	CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder ( m ) { };
	
	// void reset ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder

	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder
	// bool in_range ( )  virtual, defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh

};  // end of class CellIterator::Over::SegmentsOfConnectedOneDimMesh
    //         ::ReverseOrder::ReverseEachCell::AssumeCellsExist

//-----------------------------------------------------------------------------------------


class CellIterator::Over::SegmentsOfConnectedOneDimMesh
          ::ReverseOrder::ReverseEachCell::BuildCells
: public CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (here negative vertex)
	// inherited from CellIterator::Over::CellsOfConnectedOneDimMesh

	// attribute Cell::Core * current_segment
	// inherited from CellIterator::Over::SegmentsOfConnectedOneDimMesh

	inline BuildCells ( Mesh::Connected::OneDim * m )
	:	CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder ( m ) { };
	
	// void reset ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder

	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder
	// bool in_range ( )  virtual, defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh

};  // end of class CellIterator::Over::SegmentsOfConnectedOneDimMesh
    //         ::ReverseOrder::ReverseEachCell::BuildCells


//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class CellIterator::Over::CellsOfFuzzyMesh : public CellIterator::Core

{	public :

	std::list<Cell> & list;
  std::list<Cell>::iterator iter;

	inline CellsOfFuzzyMesh ( Mesh::Fuzzy * msh, const tag::CellsOfDim &, const size_t d )
	:	CellIterator::Core (), list { msh->cells[d] }
	{	}  // no need to initialize 'iter', 'reset' will do that
	
	void reset ( );  // virtual from CellIterator::Core

	void reset ( const tag::StartAt &, Cell::Core * cll );
	// virtual from CellIterator::Core, here execution forbidden

	// Cell deref ( )  stays pure virtual from CellIterator::Core
	void advance ( );  // virtual from CellIterator::Core
	bool in_range ( );  // virtual from CellIterator::Core

	class AsTheyAre;  class ForcePositive;
	struct ReverseEachCell  {  class AssumeCellsExist;  };
	
};  // end of class CellIterator::Over::CellsOfFuzzyMesh

//-----------------------------------------------------------------------------------------


class CellIterator::Over::CellsOfFuzzyMesh::AsTheyAre
: public CellIterator::Over::CellsOfFuzzyMesh

{	public :

	// inherited from CellIterator::Over::CellsOfFuzzyMesh :
	// std::list<Cell> & list;
  // std::list<Cell>::iterator iter;

	inline AsTheyAre ( Mesh::Fuzzy * msh, const tag::CellsOfDim &, const size_t d )
		:	CellIterator::Over::CellsOfFuzzyMesh ( msh, tag::cells_of_dim, d )
	{	}  // no need to initialize 'iter', 'reset' will do that
	
	// void reset ( )  virtual, defined by CellIterator::Over::CellsOfFuzzyMesh

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by CellIterator::Over::CellsOfFuzzyMesh, execution forbidden

	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  virtual, defined by CellIterator::Over::CellsOfFuzzyMesh
	// bool in_range ( )  virtual, defined by CellIterator::Over::CellsOfFuzzyMesh

};  // end of class CellIterator::Over::CellsOfFuzzyMesh::AsTheyAre

//-----------------------------------------------------------------------------------------


class CellIterator::Over::CellsOfFuzzyMesh::ForcePositive
: public CellIterator::Over::CellsOfFuzzyMesh

{	public :

	// inherited from CellIterator::Over::CellsOfFuzzyMesh :
	// std::list<Cell> & list;
  // std::list<Cell>::iterator iter;

	inline ForcePositive ( Mesh::Fuzzy * msh, const tag::CellsOfDim &, const size_t d )
		:	CellIterator::Over::CellsOfFuzzyMesh ( msh, tag::cells_of_dim, d )
	{	}  // no need to initialize 'iter', 'reset' will do that
	
	// void reset ( )  virtual, defined by CellIterator::Over::CellsOfFuzzyMesh

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by CellIterator::Over::CellsOfFuzzyMesh, execution forbidden

	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  virtual, defined by CellIterator::Over::CellsOfFuzzyMesh
	// bool in_range ( )  virtual, defined by CellIterator::Over::CellsOfFuzzyMesh

};  // end of class CellIterator::Over::CellsOfFuzzyMesh::ForcePositive

//-----------------------------------------------------------------------------------------


class CellIterator::Over::CellsOfFuzzyMesh::ReverseEachCell::AssumeCellsExist
: public CellIterator::Over::CellsOfFuzzyMesh

{	public :

	// inherited from CellIterator::Over::CellsOfFuzzyMesh :
	// std::list<Cell> & list;
  // std::list<Cell>::iterator iter;

	inline AssumeCellsExist ( Mesh::Fuzzy * msh, const tag::CellsOfDim &, const size_t d )
		:	CellIterator::Over::CellsOfFuzzyMesh ( msh, tag::cells_of_dim, d )
	{	}  // no need to initialize 'iter', 'reset' will do that
	
	// void reset ( )  virtual, defined by CellIterator::Over::CellsOfFuzzyMesh

	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//   virtual, defined by CellIterator::Over::CellsOfFuzzyMesh, execution forbidden

	Cell deref ( );  // virtual from CellIterator::Core
	// we trust each cell has already a reverse
	// void advance ( )  virtual, defined by CellIterator::Over::CellsOfFuzzyMesh
	// bool in_range ( )  virtual, defined by CellIterator::Over::CellsOfFuzzyMesh

};  // end of class CellIterator::Over::CellsOfFuzzyMesh
    //                ::ReverseEachCell::AssumeCellsExist


//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

//    ITERATORS OVER CELLS AROUND A GIVEN CELL (CENTERED AT CELL)

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


// see paragraph 9.12 in the manual

class CellIterator::Around::OneCell : public CellIterator::Core

// we may want to run over all segments starting at a given vertex
// or over all squares sharing a given edge (in a mesh of cubes)

// the dimension of the cells to be returned varies
// e.g. we can return segments or squares (or cubes) around a vertex
// or squares (or cubes) around a segment

{	public :

	Mesh::Core * msh_p;  // should an iterator keep the mesh alive ?
	Cell::Positive * center;  // center->dim <= msh->dim - 2

	inline OneCell ( Mesh::Core * msh, Cell::Positive * const c )
	:	CellIterator::Core (), msh_p { msh }, center { c } { }

	// void reset ( )  stays pure virtual from CellIterator::Core
	// void reset ( tag::StartAt, Cell::Core * cll )  stays pure virtual from CellIterator::Core
	
	// Cell deref ( )  stays pure virtual from CellIterator::Core
	// void advance ( )  stays pure virtual from CellIterator::Core
	// bool in_range ( )  stays pure virtual from CellIterator::Core

	struct OfCodimTwo
	{ class Ordered; class OverVertices; class OverSegments; };
	class OfAnyCodim;
	
};  // end of class CellIterator::Around::OneCell

//-----------------------------------------------------------------------------------------


class CellIterator::Around::OneCell::OfCodimTwo::Ordered
: public CellIterator::Around::OneCell

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

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	Cell::Core * current_segment { nullptr };  // word "segment" is misleading here

	inline Ordered ( Mesh::Core * msh, Cell::Positive * const c )
	:	CellIterator::Around::OneCell ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  stays pure virtual from CellIterator::Core
	// void reset ( tag::StartAt, Cell::Core * cll )  stays pure virtual from CellIterator::Core
	
	// void reset ( tag::StartAt, Cell::Core * cll )  stays pure virtual from CellIterator::Core
	// Cell deref ( )  stays pure virtual from CellIterator::Core
	// void advance ( )  stays pure virtual from CellIterator::Core
	// bool in_range ( )  stays virtual from CellIterator::Core
	
};  // end of class CellIterator::Around::OneCell::OfCodimTwo::Ordered

//-----------------------------------------------------------------------------------------


class CellIterator::Around::OneCell::OfCodimTwo::OverVertices
: public CellIterator::Around::OneCell::OfCodimTwo::Ordered

// when cll->dim == msh->dim - 2, there is a linear order
// e.g. we rotate around a vertex in a 2D mesh
// or we rotate around a segment in a 3D mesh

// the dimension of the cells to be returned is a different matter
// here we return segments around a vertex or squares around a segment

// "Vertices" means we want to get cells of dimension center->dim + 1 == msh->dim - 1
// the procedure is analogous to CellIterator::Over::VerticesOfConnectedOneDimMesh

// here, words "vertex" and "segment" are quite misleading
// they come from an analogy to iterators over chains of segments
// "vertex" means a cell of dimension center->dim + 1 == msh->dim - 1
// "segment" means a cell of dimension center->dim + 2 == msh->dim
// when we rotate around a vertex, "vertex" means segment, "segment" means square
// when we rotate around a segment, "vertex" means square, "segment" means cube

{	public :

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	Cell::Core * current_vertex { nullptr };  // word "vertex" is misleading here
	Cell::Core * first_vertex { nullptr };  // word "vertex" is misleading here
	Cell::Core * last_segment { nullptr };  // word "segment" is misleading here

	inline OverVertices ( Mesh::Core * msh, Cell::Positive * const c )
	:	CellIterator::Around::OneCell::OfCodimTwo::Ordered ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  both stay pure virtual from CellIterator::Core
	// void reset ( const tag::StartAt &, Cell::Core * cll ) 
	
	// Cell deref ( )  stays pure virtual from CellIterator::Core
	// void advance ( )  stays pure virtual from CellIterator::Core

	bool in_range ( );  // virtual from CellIterator::Core
	
	class NormalOrder;  class ReverseOrder;

};  // end of class CellIterator::Around::OneCell::OfCodimTwo::OverVertices


//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder
: public CellIterator::Around::OneCell::OfCodimTwo::OverVertices

// when cll->dim == msh->dim - 2, there is a linear order
// e.g. we rotate around a vertex in a 2D mesh
// or we rotate around a segment in a 3D mesh
// here we follow this order

// the dimension of the cells to be returned is a different matter
// here we return segments around a vertex or squares around a segment

{	public :

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from CellIterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from CellIterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	// we recognize a closed loop if  last_segment == msh.cell_behind ( first_vertex )
	// another possible criterion :  last_segment == nullptr  ==>  open chain

	inline NormalOrder ( Mesh::Core * msh, Cell::Positive * const c )
	:	CellIterator::Around::OneCell::OfCodimTwo::OverVertices ( msh, c )
	// msh_p { msh }, center { c }
	{ }
	
	void reset ( );  // virtual from CellIterator::Core
	void reset ( const tag::StartAt &, Cell::Core * cll );  // virtual from CellIterator::Core

	// Cell deref ( )  stays pure virtual from CellIterator::Core

	void advance ( );  // virtual from CellIterator::Core
	
	// bool in_range ( )  virtual, defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices
	
	class AssumeCellsExist;  class BuildReverseCells;

};  // end of class CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder

//-----------------------------------------------------------------------------------------


class CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder
: public CellIterator::Around::OneCell::OfCodimTwo::OverVertices

// when cll->dim == msh->dim - 2, there is a linear order
// e.g. we rotate around a vertex in a 2D mesh
// or we rotate around a segment in a 3D mesh
// here we follow the reversed order

// the dimension of the cells to be returned is a different matter
// here we return segments around a vertex or squares around a segment

{	public :

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from CellIterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from CellIterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	// we recognize a closed loop if  last_segment == msh.cell_in_front_of ( first_vertex )
	// another possible criterion :  last_segment == nullptr  ==>  open chain

	inline ReverseOrder ( Mesh::Core * msh, Cell::Positive * const c );
	
	void reset ( );  // virtual from CellIterator::Core
	void reset ( const tag::StartAt &, Cell::Core * cll );  // virtual from CellIterator::Core
	
	// Cell deref ( )  stays pure virtual from CellIterator::Core

	// void advance ( )  stays pure virtual from CellIterator::Core
	
	// bool in_range ( )  virtual, defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices

	class AssumeCellsExist;  class BuildReverseCells;

};  // end of class CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder::AssumeCellsExist
: public CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder

{	public :

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from CellIterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from CellIterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline AssumeCellsExist ( Mesh::Core * msh, Cell::Positive * const c );
	
	// void reset ( )  virtual
	// void reset ( const tag::StartAt &, Cell::Core * cll )  virtual
	//   both defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder
	
	// Cell deref ( )  stays pure virtual from CellIterator::Core

	// void advance ( )  virtual from CellIterator::Core, defined by
	//   CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder

	// bool in_range ( )  virtual, defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices

	class AsTheyAre;  class ForcePositive;  class ReverseEachCell;

};  // end of class CellIterator::Around::OneCell::OfCodimTwo
    //                   ::OverVertices::NormalOrder::AssumeCellsExist

//-----------------------------------------------------------------------------------------


class CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder::BuildReverseCells
: public CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder

{	public :

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from CellIterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from CellIterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline BuildReverseCells ( Mesh::Core * msh, Cell::Positive * const c );
	
	// void reset ( )  virtual
	// void reset ( const tag::StartAt &, Cell::Core * cll )  virtual
	//   both defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder
	
	// Cell deref ( )  stays pure virtual from CellIterator::Core

	// void advance ( )  virtual from CellIterator::Core, defined by
	//   CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder

	// bool in_range ( )  virtual, defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices

	class AsTheyAre;  class ForcePositive;  class ReverseEachCell;

};  // end of class CellIterator::Around::OneCell::OfCodimTwo
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
		// do we accept only 11 here, or is 20 02 acceptable too ?
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
( CellIterator::Around::OneCell::OfCodimTwo::OverVertices * that )

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
( CellIterator::Around::OneCell::OfCodimTwo::OverVertices * that )
	
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

}  // namespace tag::local_functions::iterator_h
	
//-----------------------------------------------------------------------------------------


inline CellIterator::Around::OneCell::OfCodimTwo
           ::OverVertices::NormalOrder::AssumeCellsExist::AssumeCellsExist
( Mesh::Core * msh, Cell::Positive * cen )
:	CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder ( msh, cen )
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
	

inline CellIterator::Around::OneCell::OfCodimTwo
           ::OverVertices::NormalOrder::BuildReverseCells::BuildReverseCells
( Mesh::Core * msh, Cell::Positive * cen )
:	CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder ( msh, cen )
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


class CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::NormalOrder::AssumeCellsExist::AsTheyAre
: public CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder::AssumeCellsExist

{	public :

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from CellIterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from CellIterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline AsTheyAre ( Mesh::Core * msh, Cell::Positive * const c )
	:	CellIterator::Around::OneCell::OfCodimTwo
	      ::OverVertices::NormalOrder::AssumeCellsExist ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual
	// void reset ( const tag::StartAt &, Cell::Core * cll )  virtual
	//   both defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder
	
	Cell deref ( );  // virtual from CellIterator::Core

	// void advance ( )  virtual from CellIterator::Core, defined by
	//   CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder

	// bool in_range ( )  virtual, defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices

};  // end of class CellIterator::Around::OneCell::OfCodimTwo
    //                   ::OverVertices::NormalOrder::AssumeCellsExist::AsTheyAre

//-----------------------------------------------------------------------------------------


class CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::NormalOrder::AssumeCellsExist::ForcePositive
: public CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder::AssumeCellsExist

{	public :

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from CellIterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from CellIterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline ForcePositive ( Mesh::Core * msh, Cell::Positive * const c )
	:	CellIterator::Around::OneCell::OfCodimTwo
	      ::OverVertices::NormalOrder::AssumeCellsExist ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual
	// void reset ( const tag::StartAt &, Cell::Core * cll )  virtual
	//   both defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder
	
	Cell deref ( );  // virtual from CellIterator::Core

	// void advance ( )  virtual from CellIterator::Core, defined by
	//   CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder

	// bool in_range ( )  virtual, defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices
	
};  // end of class CellIterator::Around::OneCell::OfCodimTwo
    //                  ::OverVertices::NormalOrder::AssumeCellsExist::ForcePositive

//-----------------------------------------------------------------------------------------


class CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::NormalOrder::AssumeCellsExist::ReverseEachCell
: public CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder::AssumeCellsExist

{	public :

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from CellIterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from CellIterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline ReverseEachCell ( Mesh::Core * msh, Cell::Positive * const c )
	:	CellIterator::Around::OneCell::OfCodimTwo
	      ::OverVertices::NormalOrder::AssumeCellsExist ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual
	// void reset ( const tag::StartAt &, Cell::Core * cll )  virtual
	//   both defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder
	
	Cell deref ( );  // virtual from CellIterator::Core

	// void advance ( )  virtual from CellIterator::Core, defined by
	//   CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder

	// bool in_range ( )  virtual, defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices
	
};  // end of class CellIterator::Around::OneCell::OfCodimTwo
    //                  ::OverVertices::NormalOrder::AssumeCellsExist::ReverseEachCell

//-----------------------------------------------------------------------------------------


class CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::NormalOrder::BuildReverseCells::AsTheyAre
: public CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder::BuildReverseCells

{	public :

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from CellIterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from CellIterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline AsTheyAre ( Mesh::Core * msh, Cell::Positive * const c )
	:	CellIterator::Around::OneCell::OfCodimTwo
	      ::OverVertices::NormalOrder::BuildReverseCells ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual
	// void reset ( const tag::StartAt &, Cell::Core * cll )  virtual
	//   both defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder
	
	Cell deref ( );  // virtual from CellIterator::Core

	// void advance ( )  virtual from CellIterator::Core, defined by
	//   CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder

	// bool in_range ( )  virtual, defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices

};  // end of class CellIterator::Around::OneCell::OfCodimTwo
    //                   ::OverVertices::NormalOrder::BuildReverseCells::AsTheyAre

//-----------------------------------------------------------------------------------------


class CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::NormalOrder::BuildReverseCells::ForcePositive
: public CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder::BuildReverseCells

{	public :

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from CellIterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from CellIterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline ForcePositive ( Mesh::Core * msh, Cell::Positive * const c )
	:	CellIterator::Around::OneCell::OfCodimTwo
	      ::OverVertices::NormalOrder::BuildReverseCells ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual
	// void reset ( const tag::StartAt &, Cell::Core * cll )  virtual
	//   both defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder
	
	Cell deref ( );  // virtual from CellIterator::Core

	// void advance ( )  virtual from CellIterator::Core, defined by
	//   CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder

	// bool in_range ( )  virtual, defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices
	
};  // end of class CellIterator::Around::OneCell::OfCodimTwo
    //                  ::OverVertices::NormalOrder::BuildReverseCells::ForcePositive

//-----------------------------------------------------------------------------------------


class CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::NormalOrder::BuildReverseCells::ReverseEachCell
: public CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder::BuildReverseCells

{	public :

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from CellIterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from CellIterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline ReverseEachCell ( Mesh::Core * msh, Cell::Positive * const c )
	:	CellIterator::Around::OneCell::OfCodimTwo
	      ::OverVertices::NormalOrder::BuildReverseCells ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual
	// void reset ( const tag::StartAt &, Cell::Core * cll )  virtual
	//   both defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder
	
	Cell deref ( );  // virtual from CellIterator::Core

	// void advance ( )  virtual from CellIterator::Core, defined by
	//   CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder

	// bool in_range ( )  virtual, defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices
	
};  // end of class CellIterator::Around::OneCell::OfCodimTwo
    //                  ::OverVertices::NormalOrder::BuildReverseCells::ReverseEachCell

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::AssumeCellsExist
: public CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder

{	public :

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from CellIterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from CellIterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline AssumeCellsExist ( Mesh::Core * msh, Cell::Positive * const c )
	:	CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder ( msh, c )
	// msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual
	// void reset ( const tag::StartAt &, Cell::Core * cll )  virtual
	//   both defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder
	
	// Cell deref ( )  stays pure virtual from CellIterator::Core

	void advance ( );  // virtual from CellIterator::Core

	// bool in_range ( )  virtual, defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices

	class AsTheyAre;  class ForcePositive;  class ReverseEachCell;

};  // end of class CellIterator::Around::OneCell::OfCodimTwo
    //                   ::OverVertices::ReverseOrder::AssumeCellsExist

//-----------------------------------------------------------------------------------------


class CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::BuildReverseCells
: public CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder

{	public :

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from CellIterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from CellIterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline BuildReverseCells ( Mesh::Core * msh, Cell::Positive * const c )
	:	CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder ( msh, c )
	// msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual
	// void reset ( const tag::StartAt &, Cell::Core * cll )  virtual
	//   both defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder
	
	// Cell deref ( )  stays pure virtual from CellIterator::Core

	void advance ( );  // virtual from CellIterator::Core

	// bool in_range ( )  virtual, defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices

	class AsTheyAre;  class ForcePositive;  class ReverseEachCell;

};  // end of class CellIterator::Around::OneCell::OfCodimTwo
    //                   ::OverVertices::ReverseOrder::BuildReverseCells

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


inline CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::ReverseOrder
( Mesh::Core * msh, Cell::Positive * cen )
:	CellIterator::Around::OneCell::OfCodimTwo::OverVertices ( msh, cen )
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


class CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::ReverseOrder::AssumeCellsExist::AsTheyAre
: public CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::AssumeCellsExist

{	public :

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from CellIterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from CellIterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline AsTheyAre ( Mesh::Core * msh, Cell::Positive * const c )
	:	CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::ReverseOrder::AssumeCellsExist ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual
	// void reset ( const tag::StartAt &, Cell::Core * cll )  virtual
	//   both defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder
	
	Cell deref ( );  // virtual from CellIterator::Core

	// void advance ( )  virtual from CellIterator::Core, defined by
	//   CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::AssumeCellsExist

	// bool in_range ( )  virtual, defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices

};  // end of class CellIterator::Around::OneCell::OfCodimTwo
    //                   ::OverVertices::ReverseOrder::AssumeCellsExist::AsTheyAre

//-----------------------------------------------------------------------------------------


class CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::ReverseOrder::AssumeCellsExist::ForcePositive
: public CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::AssumeCellsExist

{	public :

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from CellIterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from CellIterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline ForcePositive ( Mesh::Core * msh, Cell::Positive * const c )
	:	CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::ReverseOrder::AssumeCellsExist ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual
	// void reset ( const tag::StartAt &, Cell::Core * cll )  virtual
	//   both defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder
	
	Cell deref ( );  // virtual from CellIterator::Core

	// void advance ( )  virtual from CellIterator::Core, defined by
	//   CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::AssumeCellsExist

	// bool in_range ( )  virtual, defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices
	
};  // end of class CellIterator::Around::OneCell::OfCodimTwo
    //                  ::OverVertices::ReverseOrder::AssumeCellsExist::ForcePositive

//-----------------------------------------------------------------------------------------


class CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::ReverseOrder::AssumeCellsExist::ReverseEachCell
: public CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::AssumeCellsExist

{	public :

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from CellIterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from CellIterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline ReverseEachCell ( Mesh::Core * msh, Cell::Positive * const c )
	:	CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::ReverseOrder::AssumeCellsExist ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual
	// void reset ( const tag::StartAt &, Cell::Core * cll )  virtual
	//   both defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder
	
	Cell deref ( );  // virtual from CellIterator::Core

	// void advance ( )  virtual from CellIterator::Core, defined by
	//   CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::AssumeCellsExist

	// bool in_range ( )  virtual, defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices
	
};  // end of class CellIterator::Around::OneCell::OfCodimTwo
    //                  ::OverVertices::ReverseOrder::AssumeCellsExist::ReverseEachCell

//-----------------------------------------------------------------------------------------


class CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::ReverseOrder::BuildReverseCells::AsTheyAre
: public CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::BuildReverseCells

{	public :

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from CellIterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from CellIterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline AsTheyAre ( Mesh::Core * msh, Cell::Positive * const c )
	:	CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::ReverseOrder::BuildReverseCells ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual
	// void reset ( const tag::StartAt &, Cell::Core * cll )  virtual
	//   both defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder
	
	Cell deref ( );  // virtual from CellIterator::Core

	// void advance ( )  virtual from CellIterator::Core, defined by
	//   CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::BuildReverseCells

	// bool in_range ( )  virtual, defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices

};  // end of class CellIterator::Around::OneCell::OfCodimTwo
    //                   ::OverVertices::ReverseOrder::BuildReverseCells::AsTheyAre

//-----------------------------------------------------------------------------------------


class CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::ReverseOrder::BuildReverseCells::ForcePositive
: public CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::BuildReverseCells

{	public :

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from CellIterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from CellIterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline ForcePositive ( Mesh::Core * msh, Cell::Positive * const c )
	:	CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::ReverseOrder::BuildReverseCells ( msh,c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual
	// void reset ( const tag::StartAt &, Cell::Core * cll )  virtual
	//   both defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder
	
	Cell deref ( );  // virtual from CellIterator::Core

	// void advance ( )  virtual from CellIterator::Core, defined by
	//   CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::BuildReverseCells

	// bool in_range ( )  virtual, defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices
	
};  // end of class CellIterator::Around::OneCell::OfCodimTwo
    //                  ::OverVertices::ReverseOrder::BuildReverseCells::ForcePositive

//-----------------------------------------------------------------------------------------


class CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::ReverseOrder::BuildReverseCells::ReverseEachCell
: public CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::BuildReverseCells

{	public :

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from CellIterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here
	
	// attributes inherited from CellIterator::Around::OneCell::OfCodimTwo::OverVertices :
	// Cell::Core * current_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * first_vertex { nullptr }  -- word "vertex" is misleading here
	// Cell::Core * last_segment { nullptr }  -- word "segment" is misleading here

	inline ReverseEachCell ( Mesh::Core * msh, Cell::Positive * const c )
	:	CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::ReverseOrder::BuildReverseCells ( msh, c )  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  virtual
	// void reset ( const tag::StartAt &, Cell::Core * cll )  virtual
	//   both defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder
	
	Cell deref ( );  // virtual from CellIterator::Core

	// void advance ( )  virtual from CellIterator::Core, defined by
	//   CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::BuildReverseCells

	// bool in_range ( )  virtual, defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices
	
};  // end of class CellIterator::Around::OneCell::OfCodimTwo
    //                  ::OverVertices::ReverseOrder::BuildReverseCells::ReverseEachCell

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class CellIterator::Around::OneCell::OfCodimTwo::OverSegments
: public CellIterator::Around::OneCell::OfCodimTwo::Ordered

// when cll->dim == msh->dim - 2, there is a linear order
// e.g. we rotate around a vertex in a 2D mesh
// or we rotate around a segment in a 3D mesh

// the dimension of the cells to be returned is a different matter
// here we return squares around a vertex or cubes around a segment

// "Segments" means we want to get cells of dimension center->dim + 2 == msh->dim
// the procedure is analogous to CellIterator::Over::SegmentsOfConnectedOneDimMesh

// here, words "vertex" and "segment" are quite misleading
// they come from an analogy to iterators over chains of segments
// "vertex" means a cell of dimension center->dim + 1 == msh->dim - 1
// "segment" means a cell of dimension center->dim + 2 == msh->dim
// when we rotate around a vertex, "vertex" means segment, "segment" means square
// when we rotate around a segment, "vertex" means square, "segment" means cube

{	public :

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from CellIterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here

	Cell::Core * first_segment { nullptr };  // word "segment" is misleading here
	Cell::Core * last_vertex { nullptr };  // word "vertex" is misleading here

	inline OverSegments ( Mesh::Core * msh, Cell::Positive * const c )
	:	CellIterator::Around::OneCell::OfCodimTwo::Ordered ( msh, c )  // msh_p { msh }, center { c }
	{ }

	void reset ( );  // virtual from CellIterator::Core
	
	// void reset ( const tag::StartAt &, Cell::Core * cll )
	//      stays pure virtual from CellIterator::Core
	
	// Cell deref ( )  stays pure virtual from CellIterator::Core
	// void advance ( )  stays pure virtual from CellIterator::Core

	bool in_range ( );  // virtual from CellIterator::Core
	
	class NormalOrder;  class ReverseOrder;

};  // end of class CellIterator::Around::OneCell::OfCodimTwo::OverSegments

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class CellIterator::Around::OneCell::OfCodimTwo::OverSegments::NormalOrder
: public CellIterator::Around::OneCell::OfCodimTwo::OverSegments

// when cll->dim == msh->dim - 2, there is a linear order
// e.g. we rotate around a vertex in a 2D mesh
// or we rotate around a segment in a 3D mesh
// here we follow this order

// the dimension of the cells to be returned is a different matter
// here we return squares around a vertex or cubes around a segment

{	public :

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from CellIterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here

	// attributes inherited from CellIterator::Around::OneCell::OfCodimTwo::OverSegments :
	// Cell::Core * first_segment { nullptr }  -- word "segment" is misleading here
	// Cell::Core * last_vertex { nullptr }  -- word "vertex" is misleading here

	inline NormalOrder ( Mesh::Core * msh, Cell::Positive * const c );
	
	// void reset ( )  -- virtual,
	//    defined by CellIterator::Around::OneCell::OfCodimTwo::OverSegments

	void reset ( const tag::StartAt &, Cell::Core * cll );  // virtual from CellIterator::Core

	// Cell deref ( )  stays pure virtual from CellIterator::Core

	void advance ( );  // virtual from CellIterator::Core
	
	// bool in_range ( )  virtual, defined by CellIterator::Around::OneCell::OfCodimTwo::OverSegments
	
	class AsTheyAre;  class ForcePositive;  class ReverseEachCell;

};  // end of class CellIterator::Around::OneCell::OfCodimTwo::OverSegments::NormalOrder

//-----------------------------------------------------------------------------------------


class CellIterator::Around::OneCell::OfCodimTwo::OverSegments::ReverseOrder
: public CellIterator::Around::OneCell::OfCodimTwo::OverSegments

// when cll->dim == msh->dim - 2, there is a linear order
// e.g. we rotate around a vertex in a 2D mesh
// or we rotate around a segment in a 3D mesh
// here we follow the reversed order

// the dimension of the cells to be returned is a different matter
// here we return squares around a vertex or cubes around a segment

{	public :

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from CellIterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here

	// attributes inherited from CellIterator::Around::OneCell::OfCodimTwo::OverSegments :
	// Cell::Core * first_segment { nullptr }  -- word "segment" is misleading here
	// Cell::Core * last_vertex { nullptr }  -- word "vertex" is misleading here

	inline ReverseOrder ( Mesh::Core * msh, Cell::Positive * const c );
	
	// void reset ( )  -- virtual,
	//    defined by CellIterator::Around::OneCell::OfCodimTwo::OverSegments

	void reset ( const tag::StartAt &, Cell::Core * cll );  // virtual from CellIterator::Core
	
	// Cell deref ( )  stays pure virtual from CellIterator::Core

	void advance ( );  // virtual from CellIterator::Core
	
	// bool in_range ( )  virtual, defined by CellIterator::Around::OneCell::OfCodimTwo::OverSegments

	class AsTheyAre;  class ForcePositive;  class ReverseEachCell;

};  // end of class CellIterator::Around::OneCell::OfCodimTwo::OverSegments::ReverseOrder

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


namespace tag::local_functions::iterator_h {
		
inline void rotate_backwards_for_seg ( CellIterator::Around::OneCell::OfCodimTwo::OverSegments * that )

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


inline void rotate_forward_for_seg ( CellIterator::Around::OneCell::OfCodimTwo::OverSegments * that )

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


inline CellIterator::Around::OneCell::OfCodimTwo::OverSegments::NormalOrder::NormalOrder
( Mesh::Core * msh, Cell::Positive * cen )
:	CellIterator::Around::OneCell::OfCodimTwo::OverSegments ( msh, cen )
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
	

inline CellIterator::Around::OneCell::OfCodimTwo::OverSegments::ReverseOrder::ReverseOrder
( Mesh::Core * msh, Cell::Positive * cen )
:	CellIterator::Around::OneCell::OfCodimTwo::OverSegments ( msh, cen )
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


class CellIterator::Around::OneCell::OfCodimTwo::OverSegments::NormalOrder::AsTheyAre
: public CellIterator::Around::OneCell::OfCodimTwo::OverSegments::NormalOrder

{	public :

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from CellIterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here

	// attributes inherited from CellIterator::Around::OneCell::OfCodimTwo::OverSegments :
	// Cell::Core * first_segment { nullptr }  -- word "segment" is misleading here
	// Cell::Core * last_vertex { nullptr }  -- word "vertex" is misleading here

	inline AsTheyAre ( Mesh::Core * msh, Cell::Positive * const c )
	:	CellIterator::Around::OneCell::OfCodimTwo::OverSegments::NormalOrder ( msh, c )
	// msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  -- virtual,
	//    defined by CellIterator::Around::OneCell::OfCodimTwo::OverSegments

	// void reset ( const tag::StartAt &, Cell::Core * cll )  -- virtual,
	//    defined by CellIterator::Around::OneCell::OfCodimTwo::OverSegments::NormalOrder

	Cell deref ( );  // virtual from CellIterator::Core

	// void advance ( )  -- virtual,
	//    defined by CellIterator::Around::OneCell::OfCodimTwo::OverSegments::NormalOrder
	
	// bool in_range ( )  virtual, defined by CellIterator::Around::OneCell::OfCodimTwo::OverSegments
	
};  // end of class CellIterator::Around::OneCell::OfCodimTwo::OverSegments::NormalOrder::AsTheyAre

//-----------------------------------------------------------------------------------------


class CellIterator::Around::OneCell::OfCodimTwo::OverSegments::ReverseOrder::AsTheyAre
: public CellIterator::Around::OneCell::OfCodimTwo::OverSegments::ReverseOrder

{	public :

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attribute inherited from CellIterator::Around::OneCell::OfCodimTwo::Ordered :
	// Cell::Core * current_segment { nullptr }  -- word "segment" is misleading here

	// attributes inherited from CellIterator::Around::OneCell::OfCodimTwo::OverSegments :
	// Cell::Core * first_segment { nullptr }  -- word "segment" is misleading here
	// Cell::Core * last_vertex { nullptr }  -- word "vertex" is misleading here

	inline AsTheyAre ( Mesh::Core * msh, Cell::Positive * const c )
	:	CellIterator::Around::OneCell::OfCodimTwo::OverSegments::ReverseOrder ( msh, c )
	// msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  -- virtual,
	//    defined by CellIterator::Around::OneCell::OfCodimTwo::OverSegments

	// void reset ( const tag::StartAt &, Cell::Core * cll )  -- virtual,
	//    defined by CellIterator::Around::OneCell::OfCodimTwo::OverSegments::ReverseOrder

	Cell deref ( );  // virtual from CellIterator::Core

	// void advance ( )  -- virtual,
	//    defined by CellIterator::Around::OneCell::OfCodimTwo::OverSegments::ReverseOrder
	
	// bool in_range ( )  virtual, defined by CellIterator::Around::OneCell::OfCodimTwo::OverSegments
	
};  // end of class CellIterator::Around::OneCell::OfCodimTwo::OverSegments::ReverseOrder::AsTheyAre

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class CellIterator::Around::OneCell::OfAnyCodim
: public CellIterator::Around::OneCell

// cll->dim <= msh->dim - 2, no order assumed
// more often than not, cll->dim < msh->dim - 2

{	public :

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	size_t dim;  // dimension of cells to be returned
	// center->dim < dim <= msh->dim

	inline OfAnyCodim ( Mesh::Core * msh, size_t d, Cell::Positive * const c )
	:	CellIterator::Around::OneCell ( msh, c ), dim { d }  // msh_p { msh }, center { c }
	{ }
	
	// void reset ( )  --  virtual from CellIterator::Core
	
	// void reset ( const tag::StartAt &, Cell::Core * cll )  --  virtual from CellIterator::Core
	
	// Cell deref ( )  stays pure virtual from CellIterator::Core

	void advance ( );  // virtual from CellIterator::Core

	bool in_range ( );  // virtual from CellIterator::Core

	class OverVertices;  class OverHDCells;

};  // end of class CellIterator::Around::OneCell::OfAnyCodim

//-------------------------------------------------------------------------------------------


class CellIterator::Around::OneVertex::OfAnyCodim::OverSegments
: public CellIterator::Around::OneCell::OfAnyCodim

// center->dim  == 0, no order assumed
// more often than not, msh->dim >= 3 (otherwise, an ordered iterator could be used)

// here we no longer use the misleading terms insstroduced in
// class CellIterator::Around::OneCell::OfCodimTwo::Ordered
// so, "Vertex" is really a vertex, "Segments" are really segments

{	public :

	typedef std::map < Cell::Positive::Segment*, short int > maptype_segs;

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attributes inherited from CellIterator::Around::OneCell::OfAnyCodim :
	// size_t dim  -- dimension of cells to be returned
	// center->dim + 1 == dim <= msh->dim

	maptype_segs::iterator it;;
	
	inline OverSegments ( Mesh::Core * msh, Cell::Positive * const c )
	:	CellIterator::Around::OneCell::OfAnyCodim ( msh, c->get_dim() + 1, c )
	// msh_p { msh }, center { c }, dim { c->get_dim() + 1 }
	{ }
	
	// void reset ( )  --  virtual from CellIterator::Core
	
	// void reset ( const tag::StartAt &, Cell::Core * cll )  --  virtual from CellIterator::Core
	
	// Cell deref ( )  stays pure virtual from CellIterator::Core

	void advance ( );  // virtual from CellIterator::Core

	bool in_range ( );  // virtual from CellIterator::Core

	struct ReverseEachCell  {  class  AssumeCellsExist;  };

};  // end of class CellIterator::Around::OneVertex::OfAnyCodim::OverSegments

//-------------------------------------------------------------------------------------------


class CellIterator::Around::OneCell::OfAnyCodim::OverHDCells
: public CellIterator::Around::OneCell::OfAnyCodim

// center->dim <= msh->dim - 2, no order assumed
// more often than not, center->dim < msh->dim - 2

// we return cells of dimension >= center->dim + 2
// but not of maximum dimension (maximum dimension means the dimension of the mesh)

// example : squares arond a vertex in mesh of cubes

{	public :

	typedef std::map < Mesh::Core*, Cell::field_to_meshes > maptype;

	// attributes inherited from CellIterator::Around::OneCell :
	// Mesh::Connected::OneDim * msh_p
	// Cell::Positive * center
	
	// attributes inherited from CellIterator::Around::OneCell::OfAnyCodim :
	// size_t dim  -- dimension of cells to be returned
	// center->dim + 2 <= dim < msh->dim

	maptype::iterator it;;
	
	inline OverHDCells ( Mesh::Core * msh, Cell::Positive * const c )
	:	CellIterator::Around::OneCell::OfAnyCodim ( msh, c->get_dim() + 1, c )
	// msh_p { msh }, center { c }, dim { c->get_dim() + 1 }
	{ }
	
	// void reset ( )  --  virtual from CellIterator::Core
	
	// void reset ( const tag::StartAt &, Cell::Core * cll )  --  virtual from CellIterator::Core
	
	// Cell deref ( )  stays pure virtual from CellIterator::Core

	void advance ( );  // virtual from CellIterator::Core

	bool in_range ( );  // virtual from CellIterator::Core

	struct ReverseEachCell  {  class  AssumeCellsExist;  };

};  // end of class CellIterator::Around::OneCell::OfAnyCodim::OverHDCells

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------


class CellIterator::Adaptor::ForcePositive : public CellIterator::Core

// modifies the behaviour of another CellIterator
// forcing the resulting cells to be positive

{	public :

	CellIterator base;

	inline ForcePositive ( CellIterator && b )
	:	base { b } { };

	void reset ( );  // virtual from CellIterator::Core
	void reset ( const tag::StartAt &, Cell::Core * cll );
	// virtual from CellIterator::Core, here execution forbidden
	// Cell deref ( )  stays pure virtual from CellIterator::Core
	void advance ( );  // virtual from CellIterator::Core
	bool in_range ( );  // virtual from CellIterator::Core

};  // end of class CellIterator::Adaptor::ForcePositive


//-------------------------------------------------------------------------------------------

//      ITERATORS OVER VERTICES
//--------------------------------------------------------------------------------
		

inline CellIterator Mesh::iterator ( const tag::OverVertices & ) const
{	return this->iterator ( tag::over_vertices, tag::as_they_are );  }


inline CellIterator Mesh::iterator
( const tag::OverVertices &, const tag::AsTheyAre & ) const
{	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::as_they_are, tag::this_mesh_is_positive ) );
	// else
	if ( this->dim() == 0 )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::reverse_each_cell, tag::do_not_build_cells,
			  tag::reverse_order, tag::this_mesh_is_positive                      ) );
	// else : dim >= 1, all vertices are positive, no need to reverse them
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_vertices, tag::as_they_are,
		  tag::reverse_order_if_any, tag::this_mesh_is_positive )     );           }


inline CellIterator Mesh::iterator
( const tag::OverVertices &, const tag::RequireOrder & ) const
{	return this->iterator ( tag::over_vertices, tag::as_they_are, tag::require_order );  }


inline CellIterator Mesh::iterator
( const tag::OverVertices &, const tag::AsTheyAre &, const tag::RequireOrder & ) const
{	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::as_they_are,
			  tag::require_order, tag::this_mesh_is_positive )           );
	// else : negative mesh
	if ( this->dim() == 0 )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::reverse_each_cell, tag::do_not_build_cells,
			  tag::reverse_order, tag::this_mesh_is_positive                      ) );
	// else : dim >= 1, all vertices are positive, no need to reverse them
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_vertices, tag::as_they_are,
		  tag::reverse_order, tag::this_mesh_is_positive )           );           }


inline CellIterator Mesh::iterator
( const tag::OverVertices &, const tag::RequireOrder &, const tag::AsTheyAre & ) const
{	return this->iterator ( tag::over_vertices, tag::as_they_are, tag::require_order );  }


inline CellIterator Mesh::iterator
( const tag::OverVertices &, const tag::Backwards & ) const
{	return this->iterator ( tag::over_vertices, tag::as_they_are, tag::backwards );  }


inline CellIterator Mesh::iterator
( const tag::OverVertices &, const tag::AsTheyAre &, const tag::Backwards & ) const
{	assert ( this->dim() <= 1 );  // no order for dim >= 2
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::as_they_are,
			  tag::reverse_order, tag::this_mesh_is_positive )            );
	// else : negative cell, back to normal order
	// but reverse cells if mesh is zero-dimensional
	if ( this->core->get_dim_plus_one() == 1 )  // dimension zero
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::reverse_each_cell,
			  tag::do_not_build_cells, tag::this_mesh_is_positive )      );
	// else : dimension one
	return CellIterator ( tag::whose_core_is, this->core->iterator
	  ( tag::over_vertices, tag::as_they_are,
		  tag::require_order, tag::this_mesh_is_positive )           );    }


inline CellIterator Mesh::iterator
( const tag::OverVertices &, const tag::Backwards &, const tag::AsTheyAre & ) const
{	return this->iterator ( tag::over_vertices, tag::as_they_are, tag::backwards );  }


inline CellIterator Mesh::iterator
( const tag::OverVertices &, const tag::ForcePositive & ) const
{	if ( this->dim() == 0 )
	{	if ( this->is_positive() )
			return CellIterator ( tag::whose_core_is, this->core->iterator
				( tag::over_vertices, tag::force_positive, tag::this_mesh_is_positive ) );
		// else : negative mesh
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::force_positive,
			  tag::reverse_order, tag::this_mesh_is_positive ) );                        }
	// else : dim >= 1, all vertices are positive
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::as_they_are, tag::this_mesh_is_positive ) );
	// else : negative mesh
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_vertices, tag::as_they_are,
		  tag::reverse_order_if_any, tag::this_mesh_is_positive ) );                     }


inline CellIterator Mesh::iterator
( const tag::OverVertices &, const tag::ForcePositive &, const tag::RequireOrder & ) const
{	assert ( this->dim() <= 1 );  // no order for dim >= 2
	if ( this->core->get_dim_plus_one() == 1 )  // dimension zero
	{	if ( this->is_positive() )
			return CellIterator ( tag::whose_core_is, this->core->iterator
				( tag::over_vertices, tag::force_positive, tag::this_mesh_is_positive ) );
		// else : negative mesh
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::force_positive,
			  tag::reverse_order, tag::this_mesh_is_positive ) );                        }
	// else : dimension one, all vertices are positive
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::as_they_are,
			  tag::require_order, tag::this_mesh_is_positive )            );
	// else : negative mesh
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_vertices, tag::as_they_are,
		  tag::reverse_order, tag::this_mesh_is_positive )            );                  }


inline CellIterator Mesh::iterator
( const tag::OverVertices &, const tag::RequireOrder &, const tag::ForcePositive & ) const
{	return this->iterator ( tag::over_vertices, tag::force_positive, tag::require_order );  }


inline CellIterator Mesh:: iterator
( const tag::OverVertices &, const tag::ForcePositive &, const tag::Backwards & ) const
{ assert ( this->dim() <= 1 );  // no order for dim >= 2
	if ( this->dim() == 0 )
	{	if ( this->is_positive() )
			return CellIterator ( tag::whose_core_is, this->core->iterator
				( tag::over_vertices, tag::force_positive,
				  tag::reverse_order, tag::this_mesh_is_positive )           );
		// else : negative mesh
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::force_positive, tag::this_mesh_is_positive ) ); }
	// else : dim == 1, all vertices are positive
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::as_they_are,
			  tag::reverse_order, tag::this_mesh_is_positive )           );
	else
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::as_they_are,
			  tag::require_order, tag::this_mesh_is_positive )            );                 }


inline CellIterator Mesh::iterator
( const tag::OverVertices &, const tag::Backwards &, const tag::ForcePositive & ) const
{	return this->iterator ( tag::over_vertices, tag::force_positive, tag::backwards );  }


inline CellIterator Mesh::iterator
( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == 0 );  // because reverse_each_cell
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::reverse_each_cell,
			  tag::do_not_build_cells, tag::this_mesh_is_positive )      );
	// else : negative mesh, reverse order
	return CellIterator ( tag::whose_core_is, this->core->iterator
	  ( tag::over_vertices, tag::as_they_are, tag::reverse_order, tag::this_mesh_is_positive ) );  }
	

inline CellIterator Mesh::iterator
( const tag::OverVertices &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::RequireOrder &      ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == 0 );  // because reverse_each_cell
	return this->iterator ( tag::over_vertices, tag::reverse_each_cell, tag::do_not_build_cells );  }
	

inline CellIterator Mesh::iterator
( const tag::OverVertices &, const tag::RequireOrder &,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == 0 );  // because reverse_each_cell
	return this->iterator ( tag::over_vertices, tag::reverse_each_cell, tag::do_not_build_cells );  }
	

inline CellIterator Mesh::iterator
( const tag::OverVertices &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Backwards &        ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == 0 );  // because reverse_each_cell
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::reverse_each_cell, tag::do_not_build_cells,
			  tag::reverse_order, tag::this_mesh_is_positive                      ) );
	// else : negative mesh, back to normal order
	return CellIterator ( tag::whose_core_is, this->core->iterator
	  ( tag::over_vertices, tag::as_they_are, tag::this_mesh_is_positive ) );   }


inline CellIterator Mesh::iterator
( const tag::OverVertices &, const tag::Backwards &,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_vertices, tag::reverse_each_cell,
                          tag::do_not_build_cells, tag::backwards    );  }

//-------------------------------------------------------------------------------------

//      ITERATORS OVER SEGMENTS
//--------------------------------------------------------------------------------
		

inline CellIterator Mesh::iterator ( const tag::OverSegments & ) const
{	return this->iterator ( tag::over_segments, tag::as_they_are );  }


inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::AsTheyAre & ) const
{	assert ( this->dim() >= 1 );
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_segments, tag::as_they_are, tag::this_mesh_is_positive ) );
	// else
	if ( this->core->get_dim_plus_one() == 2 )  // dimension one
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_segments, tag::reverse_each_cell, tag::do_not_build_cells,
			  tag::reverse_order_if_any, tag::this_mesh_is_positive               ) );
	// else : dim >= 2, all segments are positive, no need to reverse them
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_segments, tag::as_they_are, tag::this_mesh_is_positive ) );   }


inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::RequireOrder & ) const
{	return this->iterator ( tag::over_segments, tag::as_they_are, tag::require_order );  }


inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::RequireOrder & ) const
{	assert ( this->dim() == 1 );  // no order for dim >= 2
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
		  ( tag::over_segments, tag::as_they_are, tag::require_order,
			  tag::this_mesh_is_positive                                ) );
	// else
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_segments, tag::reverse_each_cell, tag::do_not_build_cells,
		  tag::reverse_order, tag::this_mesh_is_positive                      ) ); }


inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::RequireOrder &, const tag::AsTheyAre & ) const
{	return this->iterator ( tag::over_segments, tag::as_they_are, tag::require_order );  }


inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::Backwards & ) const
{	return this->iterator ( tag::over_segments, tag::as_they_are, tag::backwards );  }


inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::Backwards & ) const
{	assert ( this->dim() == 1 );  // no order for dim >= 2
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_segments, tag::as_they_are,
			  tag::reverse_order, tag::this_mesh_is_positive )            );
	// else
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_segments, tag::reverse_each_cell, tag::do_not_build_cells,
		  tag::require_order, tag::this_mesh_is_positive                      ) ); }


inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::Backwards &, const tag::AsTheyAre & ) const
{	return this->iterator ( tag::over_segments, tag::as_they_are, tag::backwards );  }


inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::ForcePositive & ) const
{	assert ( this->dim() >= 1 );
	if ( this->core->get_dim_plus_one() == 2 )  // dimension one
	{	if ( this->is_positive() )
			return CellIterator ( tag::whose_core_is, this->core->iterator
				( tag::over_segments, tag::force_positive, tag::this_mesh_is_positive ) );
		// else : negative mesh
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_segments, tag::force_positive,
				tag::reverse_order_if_any, tag::this_mesh_is_positive ) );                 }
	// else : dim >= 2, all segments are positive
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_segments, tag::as_they_are, tag::this_mesh_is_positive ) );
	// else : negative mesh
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_segments, tag::as_they_are,
		  tag::reverse_order_if_any, tag::this_mesh_is_positive )     );                  }


inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::ForcePositive &, const tag::RequireOrder & ) const
{	assert ( this->dim() == 1 );  // no order for dim >= 2
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_segments, tag::force_positive,
			  tag::require_order, tag::this_mesh_is_positive ) );
	else
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_segments, tag::force_positive,
			  tag::reverse_order, tag::this_mesh_is_positive )            );  }


inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::RequireOrder &, const tag::ForcePositive & ) const
{	return this->iterator ( tag::over_segments, tag::force_positive, tag::require_order );  }


inline CellIterator Mesh:: iterator
( const tag::OverSegments &, const tag::ForcePositive &, const tag::Backwards & ) const
{ assert ( this->dim() == 1 );  // no order for dim >= 2
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_segments, tag::force_positive,
			  tag::reverse_order, tag::this_mesh_is_positive ) );
	// else : negative mesh
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_segments, tag::force_positive,
		  tag::require_order, tag::this_mesh_is_positive ) );           }


inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::Backwards &, const tag::ForcePositive & ) const
{	return this->iterator ( tag::over_segments, tag::force_positive, tag::backwards );  }


inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == 1 );  // because reverse_each_cell
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_segments, tag::reverse_each_cell,
			  tag::do_not_build_cells, tag::this_mesh_is_positive )      );
	// else : negative mesh, reverse order
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_segments, tag::as_they_are,
		  tag::reverse_order_if_any, tag::this_mesh_is_positive )    );   }
	

inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::RequireOrder &      ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == 1 );  // because reverse_each_cell
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_segments, tag::reverse_each_cell, tag::do_not_build_cells,
			  tag::require_order, tag::this_mesh_is_positive                      ) );
	// else : negative mesh, reverse order
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_segments, tag::as_they_are,
		  tag::reverse_order, tag::this_mesh_is_positive )           );              }
	

inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::RequireOrder &,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_segments, tag::reverse_each_cell,
                          tag::do_not_build_cells, tag::require_order );  }
	

inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Backwards &      ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == 1 );  // because reverse_each_cell
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_segments, tag::reverse_each_cell, tag::do_not_build_cells,
			  tag::reverse_order, tag::this_mesh_is_positive                      ) );
	// else : negative mesh, back to normal order
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_segments, tag::as_they_are,
		  tag::require_order, tag::this_mesh_is_positive )           );              }


inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::Backwards &,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_segments, tag::reverse_each_cell,
                          tag::do_not_build_cells, tag::backwards    );  }


//-------------------------------------------------------------------------------------

//      ITERATORS OVER HIGH DIMENSIONAL CELLS
//--------------------------------------------------------------------------------
		

inline CellIterator Mesh::iterator ( const tag::OverCellsOfDim &, const size_t d ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are );  }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre & ) const
{	assert ( this->dim() >= d );
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_dim, d, tag::as_they_are, tag::this_mesh_is_positive ) );
	// else : negative mesh
	if ( this->dim() == d )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::reverse_each_cell, tag::do_not_build_cells,
			  tag::reverse_order_if_any, tag::this_mesh_is_positive                       ) );
	// else : dim > d, all cells are positive, no need to reverse them
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_dim, d, tag::as_they_are,
		  tag::reverse_order_if_any, tag::this_mesh_is_positive )      );                    }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::RequireOrder & ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are, tag::require_order );  }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d,
  const tag::AsTheyAre &, const tag::RequireOrder & ) const
{	assert ( this->dim() <= 1 );  // no order for dim >= 2
	assert ( this->dim() >= d );
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
		  ( tag::over_cells_of_dim, d, tag::as_they_are,
			  tag::require_order, tag::this_mesh_is_positive )           );
	// else : negative mesh
	if ( this->dim() == d )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::reverse_each_cell, tag::do_not_build_cells,
			  tag::reverse_order, tag::this_mesh_is_positive                              ) );
	// else : dim == 1 > d == 0, all cells are positive, no need to reverse them
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_vertices, tag::as_they_are,
		  tag::reverse_order, tag::this_mesh_is_positive )           );                       }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d,
  const tag::RequireOrder &, const tag::AsTheyAre & ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are, tag::require_order );  }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::Backwards & ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are, tag::backwards );  }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d,
  const tag::AsTheyAre &, const tag::Backwards & ) const
{	assert ( this->dim() <= 1 );  // no order for dim >= 2
	assert ( this->dim() >= d );
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_dim, d, tag::as_they_are,
			  tag::reverse_order, tag::this_mesh_is_positive )            );
	// else : negative mesh
	if ( this->dim() == d )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::reverse_each_cell, tag::do_not_build_cells,
			  tag::require_order, tag::this_mesh_is_positive                              ) );
	// else : dim == 1 > d == 0, all cells are positive, no need to reverse them
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_vertices, tag::as_they_are,
		  tag::require_order, tag::this_mesh_is_positive )           );                       }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d,
  const tag::Backwards &, const tag::AsTheyAre & ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are, tag::backwards );  }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive & ) const
{	assert ( this->dim() >= d );
	if ( this->dim() == d )
	{	if ( this->is_positive() )
			return CellIterator ( tag::whose_core_is, this->core->iterator
				( tag::over_cells_of_max_dim, tag::force_positive, tag::this_mesh_is_positive ) );
		// else : negative mesh
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::force_positive,
			  tag::reverse_order_if_any, tag::this_mesh_is_positive )    );                      }
	// else : dim > d, all cells are positive
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_dim, d, tag::as_they_are, tag::this_mesh_is_positive ) );
	// else : negative mesh
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_dim, d, tag::as_they_are,
		  tag::reverse_order_if_any, tag::this_mesh_is_positive )     );                          }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d,
  const tag::ForcePositive &, const tag::RequireOrder & ) const
{	assert ( this->dim() <= 1 );  // no order for dim >= 2
	assert ( this->dim() >= d );
	if ( this->dim() == d )
	{	if ( this->is_positive() )
			return CellIterator ( tag::whose_core_is, this->core->iterator
				( tag::over_cells_of_max_dim, tag::force_positive,
				  tag::require_order, tag::this_mesh_is_positive   )         );
		// else : negative mesh
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::force_positive,
			  tag::reverse_order, tag::this_mesh_is_positive   )          );  }
	// else : dim == 1 > d == 0, all cells are positive
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::as_they_are,
			  tag::require_order, tag::this_mesh_is_positive )           );
	// else : negative mesh
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_vertices, tag::as_they_are,
		  tag::reverse_order, tag::this_mesh_is_positive )           );       }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d,
  const tag::RequireOrder &, const tag::ForcePositive & ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::force_positive, tag::require_order );  }


inline CellIterator Mesh:: iterator
( const tag::OverCellsOfDim &, const size_t d,
  const tag::ForcePositive &, const tag::Backwards & ) const
{ assert ( this->dim() <= 1 );  // no order for dim >= 2
	assert ( this->dim() >= d );
	if ( this->dim() == d )
	{	if ( this->is_positive() )
			return CellIterator ( tag::whose_core_is, this->core->iterator
				( tag::over_cells_of_max_dim, tag::force_positive,
				  tag::reverse_order, tag::this_mesh_is_positive   )         );
		// else : negative mesh
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::force_positive,
			  tag::require_order, tag::this_mesh_is_positive   )          );  }
	// else : dim == 1 > d == 0, all cells are positive
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::as_they_are,
			  tag::reverse_order, tag::this_mesh_is_positive )           );
	// else : negative mesh
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_vertices, tag::as_they_are,
		  tag::require_order, tag::this_mesh_is_positive )           );       }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d,
  const tag::Backwards &, const tag::ForcePositive & ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::force_positive, tag::backwards );  }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == d );  // because reverse_each_cell
	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
                          tag::do_not_build_cells                             ); }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::RequireOrder &                   ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == d );  // because reverse_each_cell
	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
	                        tag::do_not_build_cells, tag::require_order        );  }
	

inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::RequireOrder &,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells &             ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == d );  // because reverse_each_cell
	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
	                        tag::do_not_build_cells, tag::require_order          );  }
	

inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Backwards &                       ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == d );  // because reverse_each_cell
	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
	                        tag::do_not_build_cells, tag::backwards            );  }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::Backwards &,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells &          ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == d );  // because reverse_each_cell
	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
	                        tag::do_not_build_cells, tag::backwards            );  }

//------------------------------------------------------------------------------------


inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::AsTheyAre &                                           ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::RequireOrder &                                        ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are, tag::require_order );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::AsTheyAre &, const tag::RequireOrder &                ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are, tag::require_order );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::RequireOrder &, const tag::AsTheyAre &                ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are, tag::require_order );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::Backwards &                                        ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are, tag::backwards );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::AsTheyAre &, const tag::Backwards &                ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are, tag::backwards );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::Backwards &, const tag::AsTheyAre &                ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are, tag::backwards );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::ForcePositive &                                       ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::force_positive );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::ForcePositive &, const tag::RequireOrder & ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::force_positive, tag::require_order );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::RequireOrder &, const tag::ForcePositive & ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::force_positive, tag::require_order );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::ForcePositive &, const tag::Backwards & ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::force_positive, tag::backwards );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::Backwards &, const tag::ForcePositive & ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::force_positive, tag::backwards );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_dim, d,
                          tag::reverse_each_cell, tag::do_not_build_cells );  }
	
inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells &, const tag::RequireOrder & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_dim, d, tag::reverse_each_cell,
                          tag::do_not_build_cells, tag::require_order       );  }
	
inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::RequireOrder &, const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_dim, d, tag::reverse_each_cell,
                          tag::do_not_build_cells, tag::require_order       );  }
	
inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells &, const tag::Backwards & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_dim, d, tag::reverse_each_cell,
                          tag::do_not_build_cells, tag::backwards           );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::Backwards &, const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_dim, d, tag::reverse_each_cell,
                          tag::do_not_build_cells, tag::backwards           );  }


//------------------------------------------------------------------------------------

//      ITERATORS OVER CELLS OF MAXIMUM DIMENSION
//------------------------------------------------------------------------------------
		

inline CellIterator Mesh::iterator ( const tag::OverCellsOfMaxDim & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::as_they_are );  }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre & ) const
{	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::as_they_are, tag::this_mesh_is_positive ) );
	// else : negative mesh
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_max_dim, tag::reverse_each_cell, tag::do_not_build_cells,
		  tag::reverse_order_if_any, tag::this_mesh_is_positive                       ) );  }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::RequireOrder & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::as_they_are, tag::require_order );  }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::RequireOrder & ) const
{	assert ( this->dim() <= 1 );  // no order for dim >= 2
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
		  ( tag::over_cells_of_max_dim, tag::as_they_are, tag::require_order,
			  tag::this_mesh_is_positive                                        ) );
	// else : negative mesh
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_max_dim, tag::reverse_each_cell, tag::do_not_build_cells,
		  tag::reverse_order, tag::this_mesh_is_positive                              ) );  }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::RequireOrder &, const tag::AsTheyAre & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::as_they_are, tag::require_order );  }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::Backwards & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::as_they_are, tag::backwards );  }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::Backwards & ) const
{	assert ( this->dim() <= 1 );  // no order for dim >= 2
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::as_they_are,
			  tag::reverse_order, tag::this_mesh_is_positive )            );
	// else : negative mesh
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_max_dim, tag::reverse_each_cell, tag::do_not_build_cells,
		  tag::require_order, tag::this_mesh_is_positive                              ) );  }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::Backwards &, const tag::AsTheyAre & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::as_they_are, tag::backwards );  }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive & ) const
{	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::force_positive, tag::this_mesh_is_positive ) );
	// else : negative mesh
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_max_dim, tag::force_positive,
		  tag::reverse_order_if_any, tag::this_mesh_is_positive )    );                     }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &, const tag::RequireOrder & ) const
{	assert ( this->dim() <= 1 );  // no order for dim >= 2
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::force_positive,
			  tag::require_order, tag::this_mesh_is_positive  )          );
	// else : negative mesh, reverse order
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_max_dim, tag::force_positive,
		  tag::reverse_order, tag::this_mesh_is_positive  )          );  }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::RequireOrder &, const tag::ForcePositive & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::force_positive,
                          tag::require_order                               );  }


inline CellIterator Mesh:: iterator
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &, const tag::Backwards & ) const
{ assert ( this->dim() <= 1 );  // no order for dim >= 2
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::force_positive,
			  tag::reverse_order, tag::this_mesh_is_positive   )         );
	// else : negative mesh
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_max_dim, tag::force_positive,
		  tag::require_order, tag::this_mesh_is_positive   )          );  }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::Backwards &, const tag::ForcePositive & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::force_positive, tag::backwards );  }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
// do not bother whether reverse cells exist or not	
{	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::reverse_each_cell,
			  tag::do_not_build_cells, tag::this_mesh_is_positive )       );
	// else : negative mesh, reverse order
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_max_dim, tag::as_they_are,
		  tag::reverse_order_if_any, tag::this_mesh_is_positive )    );   }
	

inline CellIterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::RequireOrder &      ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() <= 1 );  // no order for dim >= 2
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::reverse_each_cell, tag::do_not_build_cells,
			  tag::require_order, tag::this_mesh_is_positive                              ) );
	// else : negative mesh, reverse order
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_max_dim, tag::as_they_are,
		  tag::reverse_order, tag::this_mesh_is_positive )           );                     }
	

inline CellIterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::RequireOrder &,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
                          tag::do_not_build_cells, tag::require_order        );  }
	

inline CellIterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Backwards &         ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == 1 );  // no order for dim >= 2
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::reverse_each_cell, tag::do_not_build_cells,
			  tag::reverse_order, tag::this_mesh_is_positive                              ) );
	// else : negative mesh, back to normal order
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_max_dim, tag::as_they_are,
		  tag::require_order, tag::this_mesh_is_positive )           );                     }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::Backwards &,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
                          tag::do_not_build_cells, tag::backwards             );  }
	
//-----------------------------------------------------------------------------//


inline CellIterator Mesh::iterator ( const tag::OverCells &, const tag::OfMaxDim & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::as_they_are );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &, const tag::AsTheyAre & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::as_they_are );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &, const tag::RequireOrder & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::as_they_are, tag::require_order );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &,
  const tag::AsTheyAre &, const tag::RequireOrder & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::as_they_are, tag::require_order );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &,
  const tag::RequireOrder &, const tag::AsTheyAre & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::as_they_are, tag::require_order );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &, const tag::Backwards & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::as_they_are, tag::backwards );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &,
  const tag::AsTheyAre &, const tag::Backwards & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::as_they_are, tag::backwards );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &,
  const tag::Backwards &, const tag::AsTheyAre & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::as_they_are, tag::backwards );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &, const tag::ForcePositive & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::force_positive );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &,
  const tag::ForcePositive &, const tag::RequireOrder & ) const
{	return this->iterator ( tag::over_cells_of_max_dim,
                          tag::force_positive, tag::require_order );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &,
  const tag::RequireOrder &, const tag::ForcePositive & ) const
{	return this->iterator ( tag::over_cells_of_max_dim,
                          tag::force_positive, tag::require_order );  }

inline CellIterator Mesh:: iterator
( const tag::OverCells &, const tag::OfMaxDim &,
  const tag::ForcePositive &, const tag::Backwards & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::force_positive, tag::backwards );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &,
  const tag::Backwards &, const tag::ForcePositive & ) const
{	return this->iterator ( tag::over_cells_of_max_dim, tag::force_positive, tag::backwards );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_max_dim,
                          tag::reverse_each_cell, tag::do_not_build_cells );  }
	
inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::RequireOrder &                     ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
                          tag::do_not_build_cells, tag::require_order        );  }
	
inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &, const tag::RequireOrder &,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells &               ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
                          tag::do_not_build_cells, tag::require_order        );  }
	
inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Backwards &                         ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
                          tag::do_not_build_cells, tag::backwards             );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &, const tag::Backwards &,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells &            ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
                          tag::do_not_build_cells, tag::backwards             );  }
	
//-----------------------------------------------------------------------------//

//      ITERATORS OVER SEGMENTS AROUND A GIVEN CELL
//------------------------------------------------------------------------------------


inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::Around &, const Cell & c ) const
{	return this->iterator ( tag::over_segments, tag::as_they_are,
                          tag::build_cells_if_necessary, tag::around, c );  }


inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::Around &, const Cell & c ) const
{	return this->iterator ( tag::over_segments, tag::as_they_are,
                          tag::build_cells_if_necessary, tag::around, c );  }


inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
  const tag::Around &, const Cell & cll                                           ) const
{	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
		  ( tag::over_segments, tag::as_they_are, tag::build_cells_if_necessary,
	      tag::around, cll.core, tag::this_mesh_is_positive                   ) );
	// else
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_segments, tag::reverse_each_cell, tag::build_cells_if_necessary,
		  tag::around, cll.core, tag::reverse_order_if_any, tag::this_mesh_is_positive ) );  }
// or, equivalently :
//	return CellIterator ( tag::whose_core_is, this->core->iterator
//		( tag::over_segments, tag::as_they_are, tag::build_cells_if_necessary,
//		  tag::reverse_order_if_any, tag::around, cll.core->reverse_attr.core,
//		  tag::this_mesh_is_positive                                           ) );  }


inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::DoNotBuildCells &,
  const tag::Around &, const Cell & c                     ) const
{	return this->iterator ( tag::over_segments, tag::as_they_are,
                          tag::do_not_build_cells, tag::around, c );  }


inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
  const tag::Around &, const Cell & cll                                           ) const
{	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
		  ( tag::over_segments, tag::as_they_are, tag::build_cells_if_necessary,
	      tag::around, cll.core, tag::this_mesh_is_positive                   ) );
	// else
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_segments, tag::reverse_each_cell, tag::do_not_build_cells,
		  tag::around, cll.core, tag::reverse_order_if_any, tag::this_mesh_is_positive ) );  }
// or, equivalently :
//	return CellIterator ( tag::whose_core_is, this->core->iterator
//		( tag::over_segments, tag::as_they_are, tag::do_not_build_cells,
//		  tag::reverse_order_if_any, tag::around, cll.core->reverse_attr.core,
//		  tag::this_mesh_is_positive                                           ) );  }


inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::Around &, const Cell & c                              ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are,
                          tag::build_cells_if_necessary, tag::around, c );  }


inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::AsTheyAre &, const tag::Around &, const Cell & c      ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are,
                          tag::build_cells_if_necessary, tag::around, c );  }


inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
	const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
  const tag::Around &, const Cell & c                              ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are,
                          tag::build_cells_if_necessary, tag::around, c );  }


inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::DoNotBuildCells &, const tag::Around &, const Cell & c ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are,
                          tag::do_not_build_cells, tag::around, c     );  }


inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
	const tag::AsTheyAre &, const tag::DoNotBuildCells &,
  const tag::Around &, const Cell & c                              ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are,
                          tag::do_not_build_cells, tag::around, c     );  }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::Around &, const Cell & c ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are,
                          tag::build_cells_if_necessary, tag::around, c );  }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::Around &, const Cell & c                                 ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are,
                          tag::build_cells_if_necessary, tag::around, c );  }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d,
	const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
  const tag::Around &, const Cell & cll                ) const
{	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_dim, d, tag::as_they_are, tag::build_cells_if_necessary,
			  tag::around, cll.core, tag::this_mesh_is_positive                          ) );
	// else
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_dim, d, tag::reverse_each_cell, tag::build_cells_if_necessary,
		  tag::around, cll.core, tag::this_mesh_is_positive                                ) );  }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::DoNotBuildCells &,
  const tag::Around &, const Cell & c                                       ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are,
                          tag::do_not_build_cells, tag::around, c     );  }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim, const size_t d,
  const tag::AsTheyAre &, const tag::DoNotBuildCells &,
  const tag::Around &, const Cell & cll                ) const
{	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_dim, d, tag::as_they_are, tag::do_not_build_cells,
			  tag::around, cll.core, tag::this_mesh_is_positive                    ) );
	// else
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_dim, d, tag::reverse_each_cell, tag::do_not_build_cells,
		  tag::around, cll.core, tag::this_mesh_is_positive                          ) );  }



}  // namespace maniFEM

#endif  // ifndef MANIFEM_ITERATOR_H
