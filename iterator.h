
// iterator.h 2021.06.12

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
}

// class CellIterator defined in mesh.h

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
	// Cell deref ( )  remains pure virtual from CellIterator::Core
	void advance ( );  // virtual from CellIterator::Core
	bool in_range ( );  // virtual from CellIterator::Core

	struct NormalOrder  {  class AsTheyAre;  class ForcePositive;
		struct ReverseEachCell  {  class  AssumeCellsExist;  };    };
	struct ReverseOrder  {  class AsTheyAre;  class ForcePositive;
		struct ReverseEachCell  {  class  AssumeCellsExist;  };    };
	
};  // end of class CellIterator::Over::TwoVerticesOfSeg


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
	// virtual, defined by CellIterator::Over::TwoVerticesOfSeg, execution forbidden
	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg

};  // end of class CellIterator::Over::TwoVerticesOfSeg::NormalOrder::AsTheyAre


class CellIterator::Over::TwoVerticesOfSeg::NormalOrder::ForcePositive
: public CellIterator::Over::TwoVerticesOfSeg
// iterate over the two vertices, first base then tip (both positive)

{	public :

	// inherited from CellIterator::Over::TwoVerticesOfSeg :
	// Cell::PositiveSegment * seg_p;
	// short unsigned int passage;

	inline ForcePositive ( Cell::PositiveSegment * seg )
	:	CellIterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg
	// void reset ( const tag::StartAt &, Cell::Core * cll )
	// virtual, defined by CellIterator::Over::TwoVerticesOfSeg, execution forbidden
	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg

};  // end of class CellIterator::Over::TwoVerticesOfSeg::NormalOrder::ForcePositive


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
	// virtual, defined by CellIterator::Over::TwoVerticesOfSeg, execution forbidden
	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg

};  // end of class CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::AsTheyAre


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
	// virtual, defined by CellIterator::Over::TwoVerticesOfSeg, execution forbidden
	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg

};  // end of class CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::ForcePositive


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
	// virtual, defined by CellIterator::Over::TwoVerticesOfSeg, execution forbidden
	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg

};  // end of class CellIterator::Over::TwoVerticesOfSeg
    //         ::NormalOrder::ReverseEachCell::AssumeCellsExist


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
	// virtual, defined by CellIterator::Over::TwoVerticesOfSeg, execution forbidden
	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  virtual, defined by CellIterator::Over::TwoVerticesOfSeg

};  // end of class CellIterator::Over::TwoVerticesOfSeg
    //         ::ReverseOrder::ReverseEachCell::AssumeCellsExist

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

	struct NormalOrder  {  class AsTheyAre;  class ForcePositive;
		struct ReverseEachCell  {  class  AssumeCellsExist;  };    };
	struct ReverseOrder  {  class AsTheyAre;  class ForcePositive;
		struct ReverseEachCell  {  class  AssumeCellsExist;  };    };
	
};  // end of class CellIterator::Over::CellsOfConnectedOneDimMesh

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


class CellIterator::Over::VerticesOfConnectedOneDimMesh :
public CellIterator::Over::CellsOfConnectedOneDimMesh

{	public :

	// attribute Mesh::Connected::OneDim * msh
	// attribute Cell::Core * last_vertex (here positive vertex)
	// inherited from CellIterator::Over::CellsOfConnectedOneDimMesh

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

};  // end of class CellIterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder::AsTheyAre


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

};  // end of class CellIterator::Over::VerticesOfConnectedOneDimMesh::ReverseOrder::AsTheyAre

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
	struct ReverseEachCell  {  class  AssumeCellsExist;  };

};  // end of class CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder


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
	// virtual, defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder
	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder
	// bool in_range ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh

};  // end of class CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder::AsTheyAre


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
	// virtual, defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder
	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder
	// bool in_range ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh

};  // end of class CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder::ForcePositive


class CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder::ReverseEachCell::AssumeCellsExist
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
	// virtual, defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder
	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder
	// bool in_range ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh

};  // end of class CellIterator::Over::SegmentsOfConnectedOneDimMesh
    //         ::NormalOrder::ReverseEachCell::AssumeCellsExist


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
	struct ReverseEachCell  {  class  AssumeCellsExist;  };

};  // end of class CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder


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
	// virtual, defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder
	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder
	// bool in_range ( )  virtual, defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh

	class ForcePositive;
	struct ReverseEachCell  {  class AssumeCellsExist;  };

};  // end of class CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder::AsTheyAre


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
	// virtual, defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder
	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder
	// bool in_range ( )  virtual, defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh

};  // end of class CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder::ForcePositive


class CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder::ReverseEachCell::AssumeCellsExist
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
	// virtual, defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder
	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder
	// bool in_range ( )  virtual, defined by CellIterator::Over::SegmentsOfConnectedOneDimMesh

};  // end of class CellIterator::Over::SegmentsOfConnectedOneDimMesh
    //         ::ReverseOrder::ReverseEachCell::AssumeCellsExist

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
	// Cell deref ( )  remains pure virtual from CellIterator::Core
	void advance ( );  // virtual from CellIterator::Core
	bool in_range ( );  // virtual from CellIterator::Core

	class AsTheyAre;  class ForcePositive;
	struct ReverseEachCell  {  class AssumeCellsExist;  };
	
};  // end of class CellIterator::Over::CellsOfFuzzyMesh


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
	// virtual, defined by CellIterator::Over::CellsOfFuzzyMesh, execution forbidden
	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  virtual, defined by CellIterator::Over::CellsOfFuzzyMesh
	// bool in_range ( )  virtual, defined by CellIterator::Over::CellsOfFuzzyMesh

};  // end of class CellIterator::Over::CellsOfFuzzyMesh::AsTheyAre


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
	// virtual, defined by CellIterator::Over::CellsOfFuzzyMesh, execution forbidden
	Cell deref ( );  // virtual from CellIterator::Core
	// void advance ( )  virtual, defined by CellIterator::Over::CellsOfFuzzyMesh
	// bool in_range ( )  virtual, defined by CellIterator::Over::CellsOfFuzzyMesh

};  // end of class CellIterator::Over::CellsOfFuzzyMesh::ForcePositive


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
	// virtual, defined by CellIterator::Over::CellsOfFuzzyMesh, execution forbidden
	Cell deref ( );  // virtual from CellIterator::Core
	// we trust each cell has already a reverse
	// void advance ( )  virtual, defined by CellIterator::Over::CellsOfFuzzyMesh
	// bool in_range ( )  virtual, defined by CellIterator::Over::CellsOfFuzzyMesh

};  // end of class CellIterator::Over::CellsOfFuzzyMesh
    //                ::ReverseEachCell::AssumeCellsExist

//--------------------------------------------------------------------------------

	
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
	// Cell deref ( )  remains pure virtual from CellIterator::Core
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
			( tag::over_vertices, tag::reverse_each_cell, tag::do_not_bother,
			  tag::reverse_order, tag::this_mesh_is_positive                  ) );
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
			( tag::over_vertices, tag::reverse_each_cell, tag::do_not_bother,
			  tag::reverse_order, tag::this_mesh_is_positive                  ) );
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
			  tag::do_not_bother, tag::this_mesh_is_positive )            );
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
( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBother & ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == 0 );  // because reverse_each_cell
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::reverse_each_cell,
			  tag::do_not_bother, tag::this_mesh_is_positive )            );
	// else : negative mesh, reverse order
	return CellIterator ( tag::whose_core_is, this->core->iterator
	  ( tag::over_vertices, tag::as_they_are, tag::reverse_order, tag::this_mesh_is_positive ) );  }
	

inline CellIterator Mesh::iterator
( const tag::OverVertices &, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::RequireOrder &      ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == 0 );  // because reverse_each_cell
	return this->iterator ( tag::over_vertices, tag::reverse_each_cell, tag::do_not_bother );  }
	

inline CellIterator Mesh::iterator
( const tag::OverVertices &, const tag::RequireOrder &,
  const tag::ReverseEachCell &, const tag::DoNotBother & ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == 0 );  // because reverse_each_cell
	return this->iterator ( tag::over_vertices, tag::reverse_each_cell, tag::do_not_bother );  }
	

inline CellIterator Mesh::iterator
( const tag::OverVertices &, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::Backwards &        ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == 0 );  // because reverse_each_cell
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_vertices, tag::reverse_each_cell, tag::do_not_bother,
			  tag::reverse_order, tag::this_mesh_is_positive                  ) );
	// else : negative mesh, back to normal order
	return CellIterator ( tag::whose_core_is, this->core->iterator
	  ( tag::over_vertices, tag::as_they_are, tag::this_mesh_is_positive ) );   }


inline CellIterator Mesh::iterator
( const tag::OverVertices &, const tag::Backwards &,
  const tag::ReverseEachCell &, const tag::DoNotBother & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_vertices, tag::reverse_each_cell,
                          tag::do_not_bother, tag::backwards      );  }

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
			( tag::over_segments, tag::reverse_each_cell, tag::do_not_bother,
			  tag::reverse_order_if_any, tag::this_mesh_is_positive            ) );
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
		( tag::over_segments, tag::reverse_each_cell, tag::do_not_bother,
		  tag::reverse_order, tag::this_mesh_is_positive                  ) );      }


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
		( tag::over_segments, tag::reverse_each_cell, tag::do_not_bother,
		  tag::require_order, tag::this_mesh_is_positive                  ) ); }


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
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBother & ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == 1 );  // because reverse_each_cell
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_segments, tag::reverse_each_cell,
			  tag::do_not_bother, tag::this_mesh_is_positive )            );
	// else : negative mesh, reverse order
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_segments, tag::as_they_are,
		  tag::reverse_order_if_any, tag::this_mesh_is_positive )    );   }
	

inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::RequireOrder &      ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == 1 );  // because reverse_each_cell
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_segments, tag::reverse_each_cell, tag::do_not_bother,
			  tag::require_order, tag::this_mesh_is_positive                  ) );
	// else : negative mesh, reverse order
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_segments, tag::as_they_are,
		  tag::reverse_order, tag::this_mesh_is_positive )           );           }
	

inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::RequireOrder &,
  const tag::ReverseEachCell &, const tag::DoNotBother & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_segments, tag::reverse_each_cell,
                          tag::do_not_bother, tag::require_order       );  }
	

inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::Backwards &      ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == 1 );  // because reverse_each_cell
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_segments, tag::reverse_each_cell, tag::do_not_bother,
			  tag::reverse_order, tag::this_mesh_is_positive                  ) );
	// else : negative mesh, back to normal order
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_segments, tag::as_they_are,
		  tag::require_order, tag::this_mesh_is_positive )           );           }


inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::Backwards &,
  const tag::ReverseEachCell &, const tag::DoNotBother & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_segments, tag::reverse_each_cell,
                          tag::do_not_bother, tag::backwards      );  }


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
			( tag::over_cells_of_max_dim, tag::reverse_each_cell, tag::do_not_bother,
			  tag::reverse_order_if_any, tag::this_mesh_is_positive                  ) );
	// else : dim > d, all cells are positive, no need to reverse them
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_dim, d, tag::as_they_are,
		  tag::reverse_order_if_any, tag::this_mesh_is_positive )      );                  }


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
			( tag::over_cells_of_max_dim, tag::reverse_each_cell, tag::do_not_bother,
			  tag::reverse_order, tag::this_mesh_is_positive                         ) );
	// else : dim == 1 > d == 0, all cells are positive, no need to reverse them
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_vertices, tag::as_they_are,
		  tag::reverse_order, tag::this_mesh_is_positive )           );                 }


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
			( tag::over_cells_of_max_dim, tag::reverse_each_cell, tag::do_not_bother,
			  tag::require_order, tag::this_mesh_is_positive                         ) );
	// else : dim == 1 > d == 0, all cells are positive, no need to reverse them
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_vertices, tag::as_they_are,
		  tag::require_order, tag::this_mesh_is_positive )           );                }


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
  const tag::ReverseEachCell &, const tag::DoNotBother & ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == d );  // because reverse_each_cell
	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell, tag::do_not_bother );  }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::RequireOrder &                       ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == d );  // because reverse_each_cell
	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
	                        tag::do_not_bother, tag::require_order              );  }
	

inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::RequireOrder &,
  const tag::ReverseEachCell &, const tag::DoNotBother &                 ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == d );  // because reverse_each_cell
	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
	                        tag::do_not_bother, tag::require_order              );  }
	

inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::Backwards &                       ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == d );  // because reverse_each_cell
	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
	                        tag::do_not_bother, tag::backwards                  );  }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::Backwards &,
  const tag::ReverseEachCell &, const tag::DoNotBother &                 ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == d );  // because reverse_each_cell
	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
	                        tag::do_not_bother, tag::backwards                  );  }

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

inline CellIterator Mesh:: iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::ForcePositive &, const tag::Backwards & ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::force_positive, tag::backwards );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::Backwards &, const tag::ForcePositive & ) const
{	return this->iterator ( tag::over_cells_of_dim, d, tag::force_positive, tag::backwards );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::ReverseEachCell &, const tag::DoNotBother & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_dim, d,
                          tag::reverse_each_cell, tag::do_not_bother );  }
	
inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::ReverseEachCell &, const tag::DoNotBother &, const tag::RequireOrder & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_dim, d, tag::reverse_each_cell,
                          tag::do_not_bother, tag::require_order            );  }
	
inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::RequireOrder &, const tag::ReverseEachCell &, const tag::DoNotBother & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_dim, d, tag::reverse_each_cell,
                          tag::do_not_bother, tag::require_order            );  }
	
inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::ReverseEachCell &, const tag::DoNotBother &, const tag::Backwards & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_dim, d, tag::reverse_each_cell,
                          tag::do_not_bother, tag::backwards                );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfDimension &, const size_t d,
  const tag::Backwards &, const tag::ReverseEachCell &, const tag::DoNotBother & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_dim, d, tag::reverse_each_cell,
                          tag::do_not_bother, tag::backwards                );  }


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
		( tag::over_cells_of_max_dim, tag::reverse_each_cell, tag::do_not_bother,
		  tag::reverse_order_if_any, tag::this_mesh_is_positive                   ) );     }


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
		( tag::over_cells_of_max_dim, tag::reverse_each_cell, tag::do_not_bother,
		  tag::reverse_order, tag::this_mesh_is_positive                          ) );  }


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
		( tag::over_cells_of_max_dim, tag::reverse_each_cell, tag::do_not_bother,
		  tag::require_order, tag::this_mesh_is_positive                          ) );  }


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
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &, const tag::DoNotBother & ) const
// do not bother whether reverse cells exist or not	
{	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::reverse_each_cell,
			  tag::do_not_bother, tag::this_mesh_is_positive      )       );
	// else : negative mesh, reverse order
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_max_dim, tag::as_they_are,
		  tag::reverse_order_if_any, tag::this_mesh_is_positive )    );   }
	

inline CellIterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::RequireOrder &      ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() <= 1 );  // no order for dim >= 2
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::reverse_each_cell, tag::do_not_bother,
			  tag::require_order, tag::this_mesh_is_positive                          ) );
	// else : negative mesh, reverse order
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_max_dim, tag::as_they_are,
		  tag::reverse_order, tag::this_mesh_is_positive )           );                 }
	

inline CellIterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::RequireOrder &,
  const tag::ReverseEachCell &, const tag::DoNotBother &     ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
                          tag::do_not_bother, tag::require_order              );  }
	

inline CellIterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::Backwards &           ) const
// do not bother whether reverse cells exist or not	
{	assert ( this->dim() == 1 );  // no order for dim >= 2
	if ( this->is_positive() )
		return CellIterator ( tag::whose_core_is, this->core->iterator
			( tag::over_cells_of_max_dim, tag::reverse_each_cell, tag::do_not_bother,
			  tag::reverse_order, tag::this_mesh_is_positive                          ) );
	// else : negative mesh, back to normal order
	return CellIterator ( tag::whose_core_is, this->core->iterator
		( tag::over_cells_of_max_dim, tag::as_they_are,
		  tag::require_order, tag::this_mesh_is_positive )           );                  }


inline CellIterator Mesh::iterator
( const tag::OverCellsOfMaxDim &, const tag::Backwards &,
  const tag::ReverseEachCell &, const tag::DoNotBother & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
                          tag::do_not_bother, tag::backwards             );  }
	
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
  const tag::ReverseEachCell &, const tag::DoNotBother & ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_max_dim,
                          tag::reverse_each_cell, tag::do_not_bother );  }
	
inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::RequireOrder &                         ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
                          tag::do_not_bother, tag::require_order              );  }
	
inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &, const tag::RequireOrder &,
  const tag::ReverseEachCell &, const tag::DoNotBother &                   ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
                          tag::do_not_bother, tag::require_order              );  }
	
inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::Backwards &                         ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
                          tag::do_not_bother, tag::backwards             );  }

inline CellIterator Mesh::iterator
( const tag::OverCells &, const tag::OfMaxDim &, const tag::Backwards &,
  const tag::ReverseEachCell &, const tag::DoNotBother &                   ) const
// do not bother whether reverse cells exist or not	
{	return this->iterator ( tag::over_cells_of_max_dim, tag::reverse_each_cell,
                          tag::do_not_bother, tag::backwards             );  }
	
//-----------------------------------------------------------------------------//


}  // namespace maniFEM

#endif  // ifndef MANIFEM_ITERATOR_H
