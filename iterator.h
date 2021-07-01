
// maniFEM iterator.h 2019.12.30

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

#ifndef MANIFEM_ITER_H
#define MANIFEM_ITER_H

#include "mesh.h"

namespace maniFEM {

namespace tag {
	struct DoNotOrderMesh { };  static const DoNotOrderMesh do_not_order_mesh;
}

	
class CellIterator

// a thin wrapper around a CellIterator::Core, with most methods delegated to 'core'

{	public :

	class Core;
	
	CellIterator::Core * core;

	inline CellIterator ( const CellIterator & it );
	inline CellIterator ( CellIterator && it );
	
	inline CellIterator
	( const tag::OverCellsOf &, const Mesh & msh, const tag::OfDimension &, size_t );
	
	inline CellIterator ( const tag::OverCellsOf &, const Mesh & msh,
                        const tag::OfDimension &, size_t, const tag::Reverse & );

	inline CellIterator ( const tag::OverCellsOf &, const Mesh & msh,
                        const tag::OfDimension &, size_t, const tag::ForcePositive & );

	inline CellIterator ( const tag::OverCellsOf &, const Mesh & msh,
	    const tag::OfDimension &, size_t, const tag::ForcePositive &, const tag::Reverse & );

	inline CellIterator ( const tag::OverCellsOf &, const Mesh & msh,
	    const tag::OfDimension &, size_t, const tag::Reverse &, const tag::ForcePositive & );

	inline ~CellIterator ( );

	inline CellIterator & operator= ( const CellIterator & it );
	inline CellIterator & operator= ( const CellIterator && it );
	
	inline void reset ( Cell & );	
	inline void reset ( );	
	inline Cell operator* ( );
	inline CellIterator & operator++ ( );
	inline CellIterator & operator++ ( int );
	inline CellIterator & advance ( );
	inline bool in_range ( );

	struct Over
	{	class CellsOfMesh;  class CellsOfPosMesh;  class CellsOfNegMesh;
		class CellsOfOneDimMesh;  class OrderedCells;
		class SegsOfChain; class SegsOfPosChain;  class SegsOfNegChain;
		class VerOfChain;  class VerOfPosChain;
		class SegsOfLoop; class SegsOfPosLoop;  class SegsOfNegLoop;
		class VerOfLoop;  class VerOfPosLoop;                              };
		
};  // end of class CellIterator


class CellIterator::Core

// iterates over all cells of a given mesh (cells of given dimension)
// for maximum dimension (equal to the dimension of the mesh) returns oriented cells
// which may be positive or negative
// for lower dimension, returns positive cells

// or, iterates over all cells above a given cell

{	public :

	inline Core ( ) { };
	
	virtual ~Core ( ) { };

	virtual void reset ( Cell::Core * cll ) = 0;
	virtual void reset ( ) = 0;
	virtual Cell::Core * deref ( ) = 0;
	virtual void advance ( ) = 0;
	virtual bool in_range ( ) = 0;
	
  virtual void order_mesh ( ) = 0;
			
};  // end of class CellIterator::Core

	
inline CellIterator::CellIterator ( const CellIterator & it)
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "ManiFEM doesn't know (yet) how to copy construct an iterator. Sorry." << std::endl;
	exit ( 1 );                                                                                        }
// define virtual self-replicator for core iterators ?
	
inline CellIterator::CellIterator ( CellIterator && it)
{	this->core = it.core;  } // it.~CellIterator();  }
// should we destroy 'it' ?  if we do,
// the core gets destroyed because that's how we define CellIterator::~CellIterator below
// then how does the new iterator work ? without a core ?
// if we don't, what happens to the memory slot occupied by 'it', more precisely by it->core ?
// is it automatically released ?  hope so ...
	
inline CellIterator & CellIterator::operator= ( const CellIterator & it )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "ManiFEM doesn't know (yet) how to copy an iterator. Sorry." << std::endl;
	exit ( 1 );                                                                                   }
// define virtual self-replicator for core iterators ?

inline CellIterator & CellIterator::operator= ( const CellIterator && it )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "ManiFEM doesn't know (yet) how to move an iterator. Sorry." << std::endl;
	exit ( 1 );                                                                                   }
// ownership transfer ? I don't know how to implement it ...

inline CellIterator::~CellIterator ( )  { delete core; }
	
inline void CellIterator::reset ( )  {  this->core->reset();  }

inline void CellIterator::reset ( Cell & cll )  {  this->core->reset(cll.core);  }
	
inline Cell CellIterator::operator* ( )
{	return Cell ( tag::whose_core_is, this->core->deref() );  }

inline CellIterator & CellIterator::operator++ ( )  {  return this->advance();  }

inline CellIterator & CellIterator::operator++ ( int )  {  return this->advance();  }

inline CellIterator & CellIterator::advance ( )
{	this->core->advance();  return *this;  }

inline bool CellIterator::in_range ( )  {  return this->core->in_range();  }
		

class CellIterator::Over::CellsOfMesh : public virtual CellIterator::Core

// iterates over all cells of a given mesh (cells of given dimension)

{	public :

	std::list<Cell::Core*> * list_p;
  std::list<Cell::Core*>::iterator iter;

	// in the future, give separate treatment to one-dimensional meshes
	
	inline CellsOfMesh ( Mesh::Core * msh, size_t dim )
	:	list_p { & ( msh->cells[dim] ) }, iter ( list_p->begin() )
	{	}

	void reset ( ); // virtual from CellIterator::Core
	void reset ( Cell::Core * cll ); // virtual from CellIterator::Core
	void advance ( ); // virtual from CellIterator::Core
	// Cell::Core * deref ( ) remains pure virtual from CellIterator::Core
	bool in_range ( ); // virtual from CellIterator::Core
	
  virtual void order_mesh ( );
			
};  // end of class CellIterator::Over::CellsOfMesh


class CellIterator::Over::CellsOfPosMesh : public CellIterator::Over::CellsOfMesh

// iterates over all cells of a positive mesh (cells of given dimension)
// for maximum dimension (equal to the dimension of the mesh) returns oriented cells
// which may be positive or negative
// for lower dimension, returns positive cells

{	public :
	
	inline CellsOfPosMesh ( Mesh::Core * msh, size_t dim )
	:	CellIterator::Over::CellsOfMesh ( msh, dim )
	{	}

	Cell::Core * deref ( ); //  virtual from CellIterator::Core

	class Positive;
};

	
class CellIterator::Over::CellsOfPosMesh::Positive : public virtual CellIterator::Over::CellsOfMesh

// iterates over all cells of a positive mesh (cells of given dimension)
// the dereference will produce always positive cells
// should be used only for cells of maximum dimension

{	public :

	inline Positive ( Mesh::Core * msh, size_t dim )
	:	CellIterator::Over::CellsOfMesh ( msh, dim )
	{	}
	
	Cell::Core * deref ( ); //  virtual from CellIterator::Core
};
	
	
class CellIterator::Over::CellsOfNegMesh : public virtual CellIterator::Over::CellsOfMesh

// iterates over all cells of a negative mesh (cells of given dimension)
// works only for cells of lower dimension, returns positive cells

{	public :

	inline CellsOfNegMesh ( Mesh::Core * msh, size_t dim )
	:	CellIterator::Over::CellsOfMesh ( msh, dim )
	{	}

	Cell::Core * deref ( ); //  virtual from CellIterator::Core

	class MaxDim;  class Positive;
};

	
class CellIterator::Over::CellsOfNegMesh::MaxDim : public virtual CellIterator::Over::CellsOfMesh

// iterates over all cells of a negative mesh (cells of maximum dimension)
// the dereference will produce reversed cells, according to the negative character of the mesh
// works only for cells of maximum dimension

{	public :
	
	inline MaxDim ( Mesh::Core * msh, size_t dim )
	:	CellIterator::Over::CellsOfMesh ( msh, dim )
	{	}

	Cell::Core * deref ( ); //  virtual from CellIterator::Core
};
	
	
class CellIterator::Over::CellsOfNegMesh::Positive : public virtual CellIterator::Over::CellsOfMesh

// iterates over all cells of a negative mesh (cells of given dimension)
// the dereference will produce always positive cells
// should be used only for cells of maximum dimension

{	public :

	inline Positive ( Mesh::Core * msh, size_t dim )
	:	CellIterator::Over::CellsOfMesh ( msh, dim )
	{	}

	Cell::Core * deref ( ); //  virtual from CellIterator::Core
};


class CellIterator::Over::OrderedCells : public virtual CellIterator::Core	
	
// iterates over a linearly ordered collection of cells
// maybe cells of a one-dimensional mesh
// or segments around a vertex or tetrahedra around a segment etc

{	public :

	Cell::Core * current_cell;

	Cell::Core * deref ( );  // virtual from CellIterator::Core

	// void reset ( Cell::Core * cll )  remains pure virtual from CellIterator::Core
	// void reset ( )  remains pure virtual from CellIterator::Core
	// bool in_range ( )  remains pure virtual from CellIterator::Core
	// void advance ( )  remains pure virtual from CellIterator::Core

	class Open;  class Loop;
};

	
class CellIterator::Over::OrderedCells::Open : public virtual CellIterator::Over::OrderedCells
	
// iterates over a linearly ordered collection of cells
// maybe cells of an (open) chain of segments
// or segments around a vertex on the boundary of a mesh
// or tetrahedra around a segment on the boundary of a mesh etc

{	public :

	bool in_range ( );  // virtual from CellIterator::Core, through IterOver::OrderedCells
	
	// void reset ( Cell::Core * cll )  remains pure virtual from CellIterator::Core,
	//    through IterOver::OrderedCells
	// void reset ( )  remains pure virtual from CellIterator::Core, through IterOver::OrderedCells
	// Cell::Core * deref ( )  defined by IterOver::OrderedCells
	// void advance ( )  remains pure virtual from CellIterator::Core, through IterOver::OrderedCells

};

	
class CellIterator::Over::OrderedCells::Loop : public virtual CellIterator::Over::OrderedCells
	
// iterates over a linearly ordered collection of cells
// maybe cells of a closed chain of segments (e.g. boundary of a triangle)
// or segments around a vertex inside a mesh
// or tetrahedra around a segment inside a mesh etc

{	public :

	Cell::Core * first_cell;
	bool fresh;
	// used to check whether we are just starting or we have already made a complete loop

	bool in_range ( );  // virtual from CellIterator::Core, through IterOver::OrderedCells

	// void reset ( Cell::Core * cll )  remains pure virtual from CellIterator::Core,
	//    through IterOver::OrderedCells
	// void reset ( )  remains pure virtual from CellIterator::Core, through IterOver::OrderedCells
	// Cell::Core * deref ( )  defined by IterOver::OrderedCells
	// void advance ( )  remains pure virtual from CellIterator::Core, through IterOver::OrderedCells

};

	
class CellIterator::Over::CellsOfOneDimMesh : public virtual CellIterator::Core

// iterates over all cells of a given open one-dimensional mesh (segments or vertices)

{	public :

	Mesh::OneDim::Positive * mesh;

	inline CellsOfOneDimMesh ( Mesh::OneDim::Positive * msh )
	: mesh { msh }
	{	msh->order();  }  // we need to order here and in each reset

	inline CellsOfOneDimMesh ( Mesh::OneDim::Positive * msh, const tag::DoNotOrderMesh & )
	: mesh { msh }
	{	}

	inline CellsOfOneDimMesh ( Mesh::OneDim::Positive * msh, const tag::Bizarre & )
	{	  }

	void order_mesh ( );  // virtual from CellIterator::Core
			
	// Cell::Core * deref ( )  remains pure virtual from CellIterator::Core
	// bool in_range ( )  remains pure virtual from CellIterator::Core
	// void advance ( )  remains pure virtual from CellIterator::Core
};


class CellIterator::Over::SegsOfChain
: public virtual CellIterator::Over::CellsOfOneDimMesh,
	public virtual CellIterator::Over::OrderedCells::Open

// iterates over all segments of a given open one-dimensional mesh
// not useful, at least not yet

{	public :

	inline SegsOfChain ( Mesh::OneDim::Positive * msh )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh )
	{	}

	inline SegsOfChain ( Mesh::OneDim::Positive * msh, const tag::DoNotOrderMesh & )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh, tag::do_not_order_mesh )
	{	}

	// void reset ( Cell::Core * cll )  remains pure virtual from CellIterator::Core,
	//    through IterOver::OrderedCells, IterOver::OrderedCells::Open
	// void reset ( )  remains pure virtual from CellIterator::Core,
	//    through IterOver::CellsOfOneDimMesh, IterOver::OrderedCells, IterOver::OrderedCells::Open
	// Cell::Core * deref ( )  defined by IterOver::OrderedCells
	// void advance ( )  remains pure virtual from CellIterator::Core,
	//   through IterOver::CellsOfOneDimMesh, IterOver::OrderedCells, IterOver::OrderedCells::Open
	// bool in_range ( )  defined by IterOver::OrderedCells::Open
};

	
class CellIterator::Over::SegsOfPosChain : public virtual CellIterator::Over::SegsOfChain

// iterates over all segments of a given open one-dimensional mesh

{	public :

	inline SegsOfPosChain ( Mesh::OneDim::Positive * msh )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh ), CellIterator::Over::SegsOfChain ( msh )
	{	}

	inline SegsOfPosChain ( Mesh::OneDim::Positive * msh, const tag::DoNotOrderMesh & )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfChain ( msh, tag::do_not_order_mesh )
	{	}

	void reset ( Cell::Core * cll );  // virtual from CellIterator::Core,
	//    through IterOver::OrderedCells, IterOver::OrderedCells::Open, CellIterator::Over::SegsOfChain

	void reset ( );  // virtual from CellIterator::Core, through IterOver::CellsOfOneDimMesh, 
	//    IterOver::OrderedCells, IterOver::OrderedCells::Open, CellIterator::Over::SegsOfChain
	
	void advance ( );  // virtual from CellIterator::Core, through IterOver::CellsOfOneDimMesh, 
	//    IterOver::OrderedCells, IterOver::OrderedCells::Open, CellIterator::Over::SegsOfChain
	
	// Cell::Core * deref ( )  defined by IterOver::OrderedCells
	// bool in_range ( )  defined by IterOver::OrderedCells::Open

	class Reverse;  class Positive;
};


class CellIterator::Over::SegsOfPosChain::Reverse :
public virtual CellIterator::Over::SegsOfChain

// iterates over all segments of a given positive one-dimensional mesh, in reverse order

{	public :

	inline Reverse ( Mesh::OneDim::Positive * msh )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh ), CellIterator::Over::SegsOfChain ( msh )
	{	}

	inline Reverse ( Mesh::OneDim::Positive * msh, const tag::DoNotOrderMesh & )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfChain ( msh, tag::do_not_order_mesh )
	{	}

	void reset ( Cell::Core * cll );  // virtual from CellIterator::Core,
	//    through IterOver::OrderedCells, IterOver::OrderedCells::Open, CellIterator::Over::SegsOfChain

	void reset ( );  // virtual from CellIterator::Core, through IterOver::CellsOfOneDimMesh, 
	//    IterOver::OrderedCells, IterOver::OrderedCells::Open, CellIterator::Over::SegsOfChain
	
	void advance ( );  // virtual from CellIterator::Core, through IterOver::CellsOfOneDimMesh, 
	//    IterOver::OrderedCells, IterOver::OrderedCells::Open, CellIterator::Over::SegsOfChain
	
	// Cell::Core * deref ( )  defined by IterOver::OrderedCells
	// bool in_range ( )  defined by IterOver::OrderedCells::Open

	class Positive;
};


class CellIterator::Over::SegsOfPosChain::Positive
: public virtual CellIterator::Over::SegsOfPosChain

// iterates over all segments of a positive open one-dimensional mesh
// returns positive segments

{	public :

	inline Positive ( Mesh::OneDim::Positive * msh )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh ),
		CellIterator::Over::SegsOfChain ( msh ),
		CellIterator::Over::SegsOfPosChain ( msh )
	{	}

	inline Positive ( Mesh::OneDim::Positive * msh, const tag::DoNotOrderMesh & )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfChain ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfPosChain ( msh, tag::do_not_order_mesh )
	{	}

	Cell::Core * deref ( ) override;
	// virtual, overrides definition by IterOver::OrderedCells

	// void reset ( Cell::Core * cll )  defined by IterOver::SegsOfPosChain
	// void reset ( )  defined by IterOver::SegsOfPosChain
	// bool in_range ( )  defined by IterOver::OrderedCells::Open
	// void advance ( )  defined by IterOver::SegsOfPosChain
};


class CellIterator::Over::SegsOfPosChain::Reverse::Positive
: public virtual CellIterator::Over::SegsOfPosChain::Reverse

// iterates over all segments of a positive open one-dimensional mesh, in reverse order
// returns positive segments

{	public :

	inline Positive ( Mesh::OneDim::Positive * msh )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh ),
		CellIterator::Over::SegsOfChain ( msh ),
		CellIterator::Over::SegsOfPosChain::Reverse ( msh )
	{	}

	inline Positive ( Mesh::OneDim::Positive * msh, const tag::DoNotOrderMesh & )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfChain ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfPosChain::Reverse ( msh, tag::do_not_order_mesh )
	{	}

	Cell::Core * deref ( ) override;
	// virtual, overrides definition by IterOver::OrderedCells

	// void reset ( Cell::Core * cll )  defined by IterOver::SegsOfPosChain
	// void reset ( )  defined by IterOver::SegsOfPosChain
	// bool in_range ( )  defined by IterOver::OrderedCells::Open
	// void advance ( )  defined by IterOver::SegsOfPosChain::Reverse
};


class CellIterator::Over::SegsOfNegChain :
public virtual CellIterator::Over::SegsOfPosChain::Reverse

// iterates over all segments of a negative open one-dimensional mesh

{	public :

	inline SegsOfNegChain ( Mesh::OneDim::Positive * msh )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh ),
		CellIterator::Over::SegsOfChain ( msh ),
		CellIterator::Over::SegsOfPosChain::Reverse ( msh )
	{	}

	inline SegsOfNegChain ( Mesh::OneDim::Positive * msh, const tag::DoNotOrderMesh & )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfChain ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfPosChain::Reverse ( msh, tag::do_not_order_mesh )
	{	}

	void reset ( Cell::Core * cll ) override;
	// virtual, overrides definition by IterOver::SegsOfPosChain::Reverse
	
	void reset ( ) override;  // virtual, overrides definition by IterOver::SegsOfPosChain::Reverse
	// the two functions are identical, but the compiler forces me to override both 'reset's

	Cell::Core * deref ( ) override;
	// virtual, overrides definition by IterOver::OrderedCells

	// bool in_range ( )  defined by IterOver::OrderedCells::Open
	// void advance ( )  defined by IterOver::SegsOfPosChain::Reverse

	class Reverse;  class Positive;
};


class CellIterator::Over::SegsOfNegChain::Reverse
:	public virtual CellIterator::Over::SegsOfPosChain

// iterates over all segments of a negative open one-dimensional mesh, in reverse order

{	public :

	inline Reverse ( Mesh::OneDim::Positive * msh )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh ),
		CellIterator::Over::SegsOfChain ( msh ),
		CellIterator::Over::SegsOfPosChain ( msh )
	{	}

	inline Reverse ( Mesh::OneDim::Positive * msh, const tag::DoNotOrderMesh & )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfChain ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfPosChain ( msh, tag::do_not_order_mesh )
	{	}

	void reset ( Cell::Core * cll ) override;
	// virtual, overrides definition by IterOver::SegsOfPosChain

	void reset ( ) override;  // virtual, overrides definition by IterOver::SegsOfPosChain
	// the two functions are identical, but the compiler forces me to override both 'reset's
	
	Cell::Core * deref ( ) override;
	// virtual, overrides definition by IterOver::OrderedCells

	// bool in_range ( )  defined by IterOver::OrderedCells::Open
	// void advance ( )  defined by IterOver::SegsOfPosChain

	class Positive;
};


class CellIterator::Over::SegsOfNegChain::Positive
: public virtual CellIterator::Over::SegsOfNegChain

// iterates over all segments of a negative open one-dimensional mesh
// returns positive segments

{	public :

	inline Positive ( Mesh::OneDim::Positive * msh )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh ),
		CellIterator::Over::SegsOfChain ( msh ),
		CellIterator::Over::SegsOfPosChain::Reverse ( msh ),
		CellIterator::Over::SegsOfNegChain ( msh )
	{	}

	inline Positive ( Mesh::OneDim::Positive * msh, const tag::DoNotOrderMesh & )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfChain ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfPosChain::Reverse ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfNegChain ( msh, tag::do_not_order_mesh )
	{	}

	Cell::Core * deref ( ) override;
	// virtual, overrides definition by IterOver::OrderedCells and by IterOver::SegsOfNegChain

	// void reset ( Cell::Core * cll )  defined by IterOver::SegsOfNegChain
	// void reset ( )  defined by IterOver::SegsOfNegChain
	// bool in_range ( )  defined by IterOver::OrderedCells::Open
	// void advance ( )  defined by IterOver::SegsOfPosChain::Reverse
};


class CellIterator::Over::SegsOfNegChain::Reverse::Positive
: public virtual CellIterator::Over::SegsOfNegChain::Reverse

// iterates over all segments of a negative open one-dimensional mesh, in reverse order
// returns positive segments

{	public :

	inline Positive ( Mesh::OneDim::Positive * msh )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh ),
		CellIterator::Over::SegsOfChain ( msh ),
		CellIterator::Over::SegsOfPosChain ( msh ),
		CellIterator::Over::SegsOfNegChain::Reverse ( msh )
	{	}

	inline Positive ( Mesh::OneDim::Positive * msh, const tag::DoNotOrderMesh & )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfChain ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfPosChain ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfNegChain::Reverse ( msh, tag::do_not_order_mesh )
	{	}

	Cell::Core * deref ( ) override;
	// virtual, overrides definition by IterOver::OrderedCells and by SegsOfNegChain::Reverse

	// void reset ( Cell::Core * cll )  defined by IterOver::SegsOfNegChain::Reverse
	// void reset ( )  defined by IterOver::SegsOfNegChain::Reverse
	// bool in_range ( )  defined by IterOver::OrderedCells::Open
	// void advance ( )  defined by IterOver::SegsOfPosChain
};


class CellIterator::Over::VerOfChain
: public virtual CellIterator::Over::CellsOfOneDimMesh,
	public virtual CellIterator::Over::OrderedCells::Open

// iterates over all vertices of a given open one-dimensional mesh

{	public :

	inline VerOfChain ( Mesh::OneDim::Positive * msh )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh )
	{	}

	inline VerOfChain ( Mesh::OneDim::Positive * msh, const tag::DoNotOrderMesh & )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh, tag::do_not_order_mesh )
	{	}

	void reset ( Cell::Core * cll );
	// virtual from CellIterator::Core, through IterOver::OrderedCells, IterOver::OrderedCells::Open

	// void advance ( )  remains pure virtual from CellIterator::Core,
	    // through IterOver::CellsOfOneDimMesh, IterOver::OrderedCells, IterOver::OrderedCells::Open

	// void reset ( )  remains pure virtual from CellIterator::Core,
	//    through IterOver::CellsOfOneDimMesh, IterOver::OrderedCells, IterOver::OrderedCells::Open
	// Cell::Core * deref ( )  defined by IterOver::OrderedCells
	// bool in_range ( )  defined by IterOver::OrderedCells::Open
};

	
class CellIterator::Over::VerOfPosChain : public virtual CellIterator::Over::VerOfChain
	
// iterates over all vertices of a positive open one-dimensional mesh

{	public :

	inline VerOfPosChain ( Mesh::OneDim::Positive * msh )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh ), CellIterator::Over::VerOfChain ( msh )
	{	}

	inline VerOfPosChain ( Mesh::OneDim::Positive * msh, const tag::DoNotOrderMesh & )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::VerOfChain ( msh, tag::do_not_order_mesh )
	{	}

	void reset ( );  // virtual from CellIterator::Core, through IterOver::CellsOfOneDimMesh, 
	//   IterOver::OrderedCells, IterOver::OrderedCells::Open, CellIterator::Over::VerOfChain
	
	void advance ( );  // virtual from CellIterator::Core, through IterOver::CellsOfOneDimMesh,
	//   IterOver::OrderedCells, IterOver::OrderedCells::Open, IterOver::VerOfChain

	// void reset ( Cell::Core * cll )  defined by CellIterator::Over::VerOfChain
	// Cell::Core * deref ( )  defined by IterOver::OrderedCells
	// bool in_range ( )  defined by IterOver::OrderedCells::Open

	class Reverse;
};


class CellIterator::Over::VerOfPosChain::Reverse : public virtual CellIterator::Over::VerOfChain

// iterates over all vertices of a positive open one-dimensional mesh, in reverse order

{	public :

	inline Reverse ( Mesh::OneDim::Positive * msh )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh ), CellIterator::Over::VerOfChain ( msh )
	{	}

	inline Reverse ( Mesh::OneDim::Positive * msh, const tag::DoNotOrderMesh & )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::VerOfChain ( msh, tag::do_not_order_mesh )
	{	}

	void reset ( );  // virtual from CellIterator::Core, through IterOver::CellsOfOneDimMesh, 
	//   IterOver::OrderedCells, IterOver::OrderedCells::Open, IterOver::VerOfChain
	
	void advance ( );  // virtual from CellIterator::Core, through IterOver::CellsOfOneDimMesh,
	//   IterOver::OrderedCells, IterOver::OrderedCells::Open, IterOver::VerOfChain

	// void reset ( Cell::Core * cll )  defined by CellIterator::Over::VerOfChain
	// Cell::Core * deref ( )  defined by IterOver::OrderedCells
	// bool in_range ( )  defined by IterOver::OrderedCells::Open
};


// no need for class IterOver::VerOfNegChain, just use IterOver::VerOfPosChain::Reverse
// no need for class IterOver::VerOfNegChain::Reverse, just use IterOver::VerOfPosChain


class CellIterator::Over::SegsOfLoop
: public virtual CellIterator::Over::CellsOfOneDimMesh,
	public virtual CellIterator::Over::OrderedCells::Loop

// iterates over all segments of a given closed one-dimensional mesh
// (for instance, the boundary of a triangle)
// not useful, at least not yet

{	public :

	inline SegsOfLoop ( Mesh::OneDim::Positive * msh )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh )
	{	}

	inline SegsOfLoop ( Mesh::OneDim::Positive * msh, const tag::DoNotOrderMesh & )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh, tag::do_not_order_mesh )
	{	}

	void reset ( Cell::Core * cll );  // virtual from CellIterator::Core,
	//    through IterOver::OrderedCells, IterOver::OrderedCells::Loop

	void reset ( );  // virtual from CellIterator::Core,
	//    through IterOver::CellsOfOneDimMesh, IterOver::OrderedCells, IterOver::OrderedCells::Loop
	
	// Cell::Core * deref ( )  defined by IterOver::OrderedCells
	// void advance ( )  remains pure virtual from CellIterator::Core,
	//   through IterOver::CellsOfOneDimMesh, IterOver::OrderedCells, IterOver::OrderedCells::Loop
	// bool in_range ( )  defined by IterOver::OrderedCells::Loop
};

	
class CellIterator::Over::SegsOfPosLoop : public virtual CellIterator::Over::SegsOfLoop

// iterates over all segments of a positive closed one-dimensional mesh
// (for instance, the boundary of a positive triangle)

{	public :

	inline SegsOfPosLoop ( Mesh::OneDim::Positive * msh )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh ), CellIterator::Over::SegsOfLoop ( msh )
	{	}

	inline SegsOfPosLoop ( Mesh::OneDim::Positive * msh, const tag::DoNotOrderMesh & )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfLoop ( msh, tag::do_not_order_mesh )
	{	}

	void advance ( );  // virtual from CellIterator::Core, through IterOver::CellsOfOneDimMesh, 
	//   IterOver::OrderedCells, IterOver::OrderedCells::Loop, CellIterator::Over::SegsOfLoop

	// void reset ( Cell::Core * cll )  defined by IterOver::SegsOfLoop
	// void reset ( )  defined by IterOver::SegsOfLoop
	// Cell::Core * deref ( )  defined by IterOver::OrderedCells
	// bool in_range ( )  defined by IterOver::OrderedCells::Loop

	class Reverse;  class Positive;
};


class CellIterator::Over::SegsOfPosLoop::Reverse
: public virtual CellIterator::Over::SegsOfLoop

// iterates over all segments of a positive closed one-dimensional mesh, in reverse order

{	public :

	inline Reverse ( Mesh::OneDim::Positive * msh )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh ), CellIterator::Over::SegsOfLoop ( msh )
	{	}

	inline Reverse ( Mesh::OneDim::Positive * msh, const tag::DoNotOrderMesh & )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfLoop ( msh, tag::do_not_order_mesh )
	{	}

	void advance ( );  // virtual from CellIterator::Core, through IterOver::CellsOfOneDimMesh, 
	//    IterOver::OrderedCells, IterOver::OrderedCells::Loop, CellIterator::Over::SegsOfLoop
	
	// void reset ( Cell::Core * cll )  defined by IterOver::SegsOfLoop
	// void reset ( )  defined by IterOver::SegsOfLoop
	// Cell::Core * deref ( )  defined by IterOver::OrderedCells
	// bool in_range ( )  defined by IterOver::OrderedCells::Loop

	class Positive;
};


class CellIterator::Over::SegsOfPosLoop::Positive : public virtual CellIterator::Over::SegsOfPosLoop

// iterates over all segments of a positive closed one-dimensional mesh
// returns positive segments

{	public :

	inline Positive ( Mesh::OneDim::Positive * msh )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh ),
		CellIterator::Over::SegsOfLoop ( msh ),
		CellIterator::Over::SegsOfPosLoop ( msh )
	{	}

	inline Positive ( Mesh::OneDim::Positive * msh, const tag::DoNotOrderMesh & )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfLoop ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfPosLoop ( msh, tag::do_not_order_mesh )
	{	}

	Cell::Core * deref ( ) override;
	// virtual, overrides definition by IterOver::OrderedCells

	// void reset ( Cell::Core * cll )  defined by IterOver::SegsOfPosLoop
	// void reset ( )  defined by IterOver::SegsOfPosLoop
	// bool in_range ( )  defined by IterOver::OrderedCells::Loop
	// void advance ( )  defined by IterOver::SegsOfPosLoop
};


class CellIterator::Over::SegsOfPosLoop::Reverse::Positive
: public virtual CellIterator::Over::SegsOfPosLoop::Reverse

// iterates over all segments of a positive closed one-dimensional mesh, in reverse order
// returns positive segments

{	public :

	inline Positive ( Mesh::OneDim::Positive * msh )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh ),
		CellIterator::Over::SegsOfLoop ( msh ),
		CellIterator::Over::SegsOfPosLoop::Reverse ( msh )
	{	}

	inline Positive ( Mesh::OneDim::Positive * msh, const tag::DoNotOrderMesh & )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfLoop ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfPosLoop::Reverse ( msh, tag::do_not_order_mesh )
	{	}

	Cell::Core * deref ( ) override;
	// virtual, overrides definition by IterOver::OrderedCells

	// void reset ( Cell::Core * cll )  defined by IterOver::SegsOfPosLoop
	// void reset ( )  defined by IterOver::SegsOfPosLoop
	// bool in_range ( )  defined by IterOver::OrderedCells::Loop
	// void advance ( )  defined by IterOver::SegsOfPosLoop::Reverse
};


class CellIterator::Over::SegsOfNegLoop
: public virtual CellIterator::Over::SegsOfPosLoop::Reverse

// iterates over all segments of a negative closed one-dimensional mesh
// (for instance, the boundary of a negative triangle)

{	public :

	inline SegsOfNegLoop ( Mesh::OneDim::Positive * msh )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh ),
		CellIterator::Over::SegsOfLoop ( msh ),
		CellIterator::Over::SegsOfPosLoop::Reverse ( msh )
	{	}

	inline SegsOfNegLoop ( Mesh::OneDim::Positive * msh, const tag::DoNotOrderMesh & )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfLoop ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfPosLoop::Reverse ( msh, tag::do_not_order_mesh )
	{	}

	void reset ( Cell::Core * cll ) override;
	// virtual, overrides definition by IterOver::SegsOfLoop
	
	void reset ( ) override;  // virtual, overrides definition by IterOver::SegsOfLoop
	// the two functions are identical, but the compiler forces me to override both 'reset's

	Cell::Core * deref ( ) override;
	// virtual, overrides definition by IterOver::OrderedCells

	// bool in_range ( )  defined by IterOver::OrderedCells::Loop
	// void advance ( )  defined by by IterOver::SegsOfPosLoop::Reverse

	class Reverse;  class Positive;
};


class CellIterator::Over::SegsOfNegLoop::Positive
: public virtual CellIterator::Over::SegsOfNegLoop

// iterates over all segments of a negative open one-dimensional mesh
// returns positive segments

{	public :

	inline Positive ( Mesh::OneDim::Positive * msh )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh ),
		CellIterator::Over::SegsOfLoop ( msh ),
		CellIterator::Over::SegsOfPosLoop::Reverse ( msh ),
		CellIterator::Over::SegsOfNegLoop ( msh )
	{	}

	inline Positive ( Mesh::OneDim::Positive * msh, const tag::DoNotOrderMesh & )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfLoop ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfPosLoop::Reverse ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfNegLoop ( msh, tag::do_not_order_mesh )
	{	}

	Cell::Core * deref ( ) override;
	// virtual, overrides definition by IterOver::OrderedCells and by IterOver::SegsOfNegLoop

	// void reset ( Cell::Core * cll )  defined by IterOver::SegsOfNegLoop
	// void reset ( )  defined by IterOver::SegsOfNegLoop
	// bool in_range ( )  defined by IterOver::OrderedCells::Loop
	// void advance ( )  defined by IterOver::SegsOfPosLoop::Reverse
};


class CellIterator::Over::SegsOfNegLoop::Reverse
: public virtual CellIterator::Over::SegsOfPosLoop

// iterates over all segments of a negative closed one-dimensional mesh
// (for instance, the boundary of a negative triangle)

{	public :

	inline Reverse ( Mesh::OneDim::Positive * msh )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh ),
		CellIterator::Over::SegsOfLoop ( msh ),
		CellIterator::Over::SegsOfPosLoop ( msh )
	{	}

	inline Reverse ( Mesh::OneDim::Positive * msh, const tag::DoNotOrderMesh & )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfLoop ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfPosLoop ( msh, tag::do_not_order_mesh )
	{	}

	void reset ( Cell::Core * cll ) override;
	// virtual, overrides definition by IterOver::SegsOfLoop
	
	void reset ( ) override;  // virtual, overrides definition by IterOver::SegsOfLoop
	// the two functions are identical, but the compiler forces me to override both 'reset's

	Cell::Core * deref ( ) override;
	// virtual, overrides definition by IterOver::OrderedCells

	// bool in_range ( )  defined by IterOver::OrderedCells::Loop
	// void advance ( )  defined by by IterOver::SegsOfPosLoop::Reverse

	class Positive;
};


class CellIterator::Over::SegsOfNegLoop::Reverse::Positive
: public virtual CellIterator::Over::SegsOfNegLoop::Reverse

// iterates over all segments of a negative closed one-dimensional mesh, in reverse order
// returns positive segments

{	public :

	inline Positive ( Mesh::OneDim::Positive * msh )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh ),
		CellIterator::Over::SegsOfLoop ( msh ),
		CellIterator::Over::SegsOfPosLoop ( msh ),
		CellIterator::Over::SegsOfNegLoop::Reverse ( msh )
	{	}

	inline Positive ( Mesh::OneDim::Positive * msh, const tag::DoNotOrderMesh & )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfLoop ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfPosLoop ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::SegsOfNegLoop::Reverse ( msh, tag::do_not_order_mesh )
	{	}

	Cell::Core * deref ( ) override;
	// virtual, overrides definition by IterOver::OrderedCells and by SegsOfNegLoop::Reverse

	// void reset ( Cell::Core * cll )  defined by IterOver::SegsOfNegLoop::Reverse
	// void reset ( )  defined by IterOver::SegsOfNegLoop::Reverse
	// bool in_range ( )  defined by IterOver::OrderedCells::Loop
	// void advance ( )  defined by IterOver::SegsOfPosLoop
};


class CellIterator::Over::VerOfLoop
: public virtual CellIterator::Over::CellsOfOneDimMesh,
	public virtual CellIterator::Over::OrderedCells::Loop

// iterates over all vertices of a given closed one-dimensional mesh
// (for instance, the boundary of a triangle)

{	public :

	inline VerOfLoop ( Mesh::OneDim::Positive * msh )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh )
	{	}

	inline VerOfLoop ( Mesh::OneDim::Positive * msh, const tag::DoNotOrderMesh & )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh, tag::do_not_order_mesh )
	{	}

	void reset ( Cell::Core * cll );  // virtual from CellIterator::Core,
	//    through IterOver::OrderedCells, IterOver::OrderedCells::Loop
		
	void reset ( );  // virtual from CellIterator::Core,
	//    through IterOver::CellsOfOneDimMesh, IterOver::OrderedCells, IterOver::OrderedCells::Loop

	// void advance ( )  remains pure virtual from CellIterator::Core,
	//    through IterOver::CellsOfOneDimMesh, IterOver::OrderedCells, IterOver::OrderedCells::Loop
	// Cell::Core * deref ( )  defined by IterOver::OrderedCells
	// bool in_range ( )  defined by IterOver::OrderedCells::Loop
};

	
class CellIterator::Over::VerOfPosLoop : public virtual CellIterator::Over::VerOfLoop

// iterates over all vertices of a given closed one-dimensional mesh
// (for instance, the boundary of a positive triangle)

{	public :

	inline VerOfPosLoop ( Mesh::OneDim::Positive * msh )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh ), CellIterator::Over::VerOfLoop ( msh )
	{	}

	inline VerOfPosLoop ( Mesh::OneDim::Positive * msh, const tag::DoNotOrderMesh & )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::VerOfLoop ( msh, tag::do_not_order_mesh )
	{	}

	void advance ( );  // virtual from CellIterator::Core, through IterOver::CellsOfOneDimMesh,
	//   IterOver::OrderedCells, IterOver::OrderedCells::Loop, IterOver::VerOfLoop

	// void reset ( Cell::Core * cll )  defined by CellIterator::Over::VerOfLoop
	// void reset ( )  defined by CellIterator::Over::VerOfLoop
	// Cell::Core * deref ( )  defined by IterOver::OrderedCells
	// bool in_range ( )  defined by IterOver::OrderedCells::Loop

	class Reverse;
};


class CellIterator::Over::VerOfPosLoop::Reverse : public virtual CellIterator::Over::VerOfLoop

// iterates over all vertices of a given closed one-dimensional mesh
// (for instance, the boundary of a negative triangle)

{	public :

	inline Reverse ( Mesh::OneDim::Positive * msh )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh ), CellIterator::Over::VerOfLoop ( msh )
	{	}

	inline Reverse ( Mesh::OneDim::Positive * msh, const tag::DoNotOrderMesh & )
	:	CellIterator::Over::CellsOfOneDimMesh ( msh, tag::do_not_order_mesh ),
		CellIterator::Over::VerOfLoop ( msh, tag::do_not_order_mesh )
	{	}

	void advance ( );  // virtual from CellIterator::Core, through IterOver::CellsOfOneDimMesh,
	//   IterOver::OrderedCells, IterOver::OrderedCells::Loop, IterOver::VerOfLoop, IterOver::VerOfPosLoop

	// void reset ( Cell::Core * cll )  defined by CellIterator::Over::VerOfLoop
	// void reset ( )  defined by CellIterator::Over::VerOfLoop
	// Cell::Core * deref ( )  defined by IterOver::OrderedCells
	// bool in_range ( )  defined by IterOver::OrderedCells::Loop
};


// no need for class IterOver::VerOfNegLoop, just use IterOver::VerOfPosLoop::Reverse
// no need for class IterOver::VerOfNegLoop::Reverse, just use IterOver::VerOfPosLoop

// ----------------------------------------------------------------------------------//


inline CellIterator::CellIterator
( const tag::OverCellsOf &, const Mesh & msh, const tag::OfDimension &, size_t dim )

{	assert ( dim <= msh.dim() );
	if ( msh.is_positive() )
		if ( msh.dim() > 1 )
			core = new CellIterator::Over::CellsOfPosMesh ( msh.core, dim );
		else
		{	assert ( msh.dim() == 1 );
			Mesh::OneDim::Positive * msh1 = (Mesh::OneDim::Positive*) msh.core;
			msh1->order();
			if ( dim == 1 )
				if ( msh1->first_ver == Cell::ghost )
					core = new CellIterator::Over::SegsOfPosLoop ( msh1, tag::do_not_order_mesh );
				else
					core = new CellIterator::Over::SegsOfPosChain ( msh1, tag::do_not_order_mesh );
			else
			{	assert ( dim == 0 );
				if ( msh1->first_ver == Cell::ghost )
					core = new CellIterator::Over::VerOfPosLoop ( msh1, tag::do_not_order_mesh );
				else
					core = new CellIterator::Over::VerOfPosChain ( msh1, tag::do_not_order_mesh );  }  }
	else  // negative mesh
		if ( msh.dim() > 1 )
			if ( msh.dim() == dim )
				core = new CellIterator::Over::CellsOfNegMesh::MaxDim ( msh.core, dim );
			else  core = new CellIterator::Over::CellsOfNegMesh ( msh.core, dim );
		else
		{	assert ( msh.dim() == 1 );
			Mesh::OneDim::Positive * msh1 = (Mesh::OneDim::Positive*) msh.core;
			msh1->order();
			if ( dim == 1 )
				if ( msh1->first_ver == Cell::ghost )
					core = new CellIterator::Over::SegsOfNegLoop ( msh1, tag::do_not_order_mesh );
				else
					core = new CellIterator::Over::SegsOfNegChain ( msh1, tag::do_not_order_mesh );
			else
			{	assert ( dim == 0 );
				if ( msh1->first_ver == Cell::ghost )
					core = new CellIterator::Over::VerOfPosLoop::Reverse ( msh1, tag::do_not_order_mesh );
				else
					core = new CellIterator::Over::VerOfPosChain::Reverse ( msh1, tag::do_not_order_mesh );  }  }
}


inline CellIterator::CellIterator ( const tag::OverCellsOf &, const Mesh & msh,
           const tag::OfDimension &, size_t dim, const tag::ForcePositive & )

{	assert ( dim <= msh.dim() );
	if ( msh.is_positive() )
		if ( msh.dim() > 1 )
			if ( msh.dim() == dim )
				core = new CellIterator::Over::CellsOfPosMesh::Positive ( msh.core, dim );
			else  core = new CellIterator::Over::CellsOfPosMesh ( msh.core, dim );
		else
		{	assert ( msh.dim() == 1 );
			Mesh::OneDim::Positive * msh1 = (Mesh::OneDim::Positive*) msh.core;
			msh1->order();
			if ( dim == 1 )
				if ( msh1->first_ver == Cell::ghost )
					core = new CellIterator::Over::SegsOfPosLoop::Positive ( msh1, tag::do_not_order_mesh );
				else
					core = new CellIterator::Over::SegsOfPosChain::Positive ( msh1, tag::do_not_order_mesh );
			else
			{	assert ( dim == 0 );
				if ( msh1->first_ver == Cell::ghost )
					core = new CellIterator::Over::VerOfPosLoop ( msh1, tag::do_not_order_mesh );
				else
					core = new CellIterator::Over::VerOfPosChain ( msh1, tag::do_not_order_mesh );  }  }
	else  // negative mesh
		if ( msh.dim() > 1 )
			if ( msh.dim() == dim )
				core = new CellIterator::Over::CellsOfNegMesh::Positive ( msh.core, dim );
			else  core = new CellIterator::Over::CellsOfNegMesh ( msh.core, dim );
		else
		{	assert ( msh.dim() == 1 );
			Mesh::OneDim::Positive * msh1 = (Mesh::OneDim::Positive*) msh.core;
			msh1->order();
			if ( dim == 1 )
				if ( msh1->first_ver == Cell::ghost )
					core = new CellIterator::Over::SegsOfNegLoop::Positive ( msh1, tag::do_not_order_mesh );
				else  core = new CellIterator::Over::SegsOfNegChain::Positive ( msh1, tag::do_not_order_mesh );
			else
			{	assert ( dim == 0 );
				if ( msh1->first_ver == Cell::ghost )
					core = new CellIterator::Over::VerOfPosLoop::Reverse ( msh1, tag::do_not_order_mesh );
				else  core = new CellIterator::Over::VerOfPosChain::Reverse ( msh1, tag::do_not_order_mesh );  }  }
}


inline CellIterator::CellIterator ( const tag::OverCellsOf &, const Mesh & msh,
           const tag::OfDimension &, size_t dim, const tag::Reverse & )

{	assert ( dim <= msh.dim() );
	assert ( msh.dim() == 1 );
	Mesh::OneDim::Positive * msh1 = (Mesh::OneDim::Positive*) msh.core;
	msh1->order();
	if ( msh.is_positive() )
		if ( dim == 1 )
			if ( msh1->first_ver == Cell::ghost )
				core = new CellIterator::Over::SegsOfPosLoop::Reverse ( msh1, tag::do_not_order_mesh );
			else
				core = new CellIterator::Over::SegsOfPosChain::Reverse ( msh1, tag::do_not_order_mesh );
		else
		{	assert ( dim == 0 );
			if ( msh1->first_ver == Cell::ghost )
				core = new CellIterator::Over::VerOfPosLoop::Reverse ( msh1, tag::do_not_order_mesh );
			else
				core = new CellIterator::Over::VerOfPosChain::Reverse ( msh1, tag::do_not_order_mesh );  }
	else  // negative mesh
		if ( dim == 1 )
			if ( msh1->first_ver == Cell::ghost )
				core = new CellIterator::Over::SegsOfNegLoop::Reverse ( msh1, tag::do_not_order_mesh );
			else
				core = new CellIterator::Over::SegsOfNegChain::Reverse ( msh1, tag::do_not_order_mesh );
		else
		{	assert ( dim == 0 );
			if ( msh1->first_ver == Cell::ghost )
				core = new CellIterator::Over::VerOfPosLoop ( msh1, tag::do_not_order_mesh );
			else
				core = new CellIterator::Over::VerOfPosChain ( msh1, tag::do_not_order_mesh );  }
}


inline CellIterator::CellIterator ( const tag::OverCellsOf &, const Mesh & msh,
	const tag::OfDimension &, size_t dim, const tag::Reverse &, const tag::ForcePositive & )

{	assert ( dim <= msh.dim() );
	assert ( msh.dim() == 1 );
	Mesh::OneDim::Positive * msh1 = (Mesh::OneDim::Positive*) msh.core;
	msh1->order();
	if ( msh.is_positive() )
		if ( dim == 1 )
			if ( msh1->first_ver == Cell::ghost )
				core = new CellIterator::Over::SegsOfPosLoop::Reverse::Positive ( msh1, tag::do_not_order_mesh );
			else
				core = new CellIterator::Over::SegsOfPosChain::Reverse::Positive ( msh1, tag::do_not_order_mesh );
		else
		{	assert ( dim == 0 );
			if ( msh1->first_ver == Cell::ghost )
				core = new CellIterator::Over::VerOfPosLoop::Reverse ( msh1, tag::do_not_order_mesh );
			else
				core = new CellIterator::Over::VerOfPosChain::Reverse ( msh1, tag::do_not_order_mesh );  }
	else  // negative mesh
		if ( dim == 1 )
			if ( msh1->first_ver == Cell::ghost )
				core = new CellIterator::Over::SegsOfNegLoop::Reverse::Positive ( msh1, tag::do_not_order_mesh );
			else
				core = new CellIterator::Over::SegsOfNegChain::Reverse::Positive ( msh1, tag::do_not_order_mesh );
		else
		{	assert ( dim == 0 );
			if ( msh1->first_ver == Cell::ghost )
				core = new CellIterator::Over::VerOfPosLoop ( msh1, tag::do_not_order_mesh );
			else
				core = new CellIterator::Over::VerOfPosChain ( msh1, tag::do_not_order_mesh );  }
}


inline CellIterator::CellIterator ( const tag::OverCellsOf &, const Mesh & msh,
	const tag::OfDimension &, size_t dim, const tag::ForcePositive &, const tag::Reverse & )
:	CellIterator ( tag::over_cells_of, msh, tag::of_dim, dim, tag::reverse, tag::force_positive )
{	}


inline CellIterator Mesh::iter_over	( const tag::CellsOfDim &, size_t d ) const
{	return CellIterator ( tag::over_cells_of, *this, tag::of_dim, d );  }

inline CellIterator Mesh::iter_over
( const tag::CellsOfDim &, size_t d, const tag::Reverse & ) const
{	return CellIterator ( tag::over_cells_of, *this, tag::of_dim, d, tag::reverse );  }

inline CellIterator Mesh::iter_over
( const tag::CellsOfDim &, size_t d, const tag::ForcePositive & ) const
{	return CellIterator ( tag::over_cells_of, *this, tag::of_dim, d, tag::force_positive );  }

inline CellIterator Mesh::iter_over
( const tag::CellsOfDim &, size_t d, const tag::ForcePositive &, const tag::Reverse & ) const
{	return CellIterator ( tag::over_cells_of, *this, tag::of_dim, d, tag::reverse, tag::force_positive );  }

inline CellIterator Mesh::iter_over
( const tag::CellsOfDim &, size_t d, const tag::Reverse &, const tag::ForcePositive & ) const
{	return CellIterator ( tag::over_cells_of, *this, tag::of_dim, d, tag::reverse, tag::force_positive );  }

inline CellIterator Mesh::iter_over ( const tag::Vertices & ) const
{	return iter_over ( tag::cells_of_dim, 0 );  }

inline CellIterator Mesh::iter_over ( const tag::Vertices &, const tag::ForcePositive & ) const
{	return iter_over ( tag::cells_of_dim, 0, tag::force_positive );  }

inline CellIterator Mesh::iter_over ( const tag::Vertices &, const tag::Reverse & ) const
{	return iter_over ( tag::cells_of_dim, 0, tag::reverse );  }

inline CellIterator Mesh::iter_over
( const tag::Vertices &, const tag::ForcePositive &, const tag::Reverse & ) const
{	return iter_over ( tag::cells_of_dim, 0, tag::force_positive, tag::reverse );  }

inline CellIterator Mesh::iter_over
( const tag::Vertices &, const tag::Reverse &, const tag::ForcePositive & ) const
{	return iter_over ( tag::cells_of_dim, 0, tag::force_positive, tag::reverse );  }

inline CellIterator Mesh::iter_over ( const tag::Segments & ) const
{	return iter_over ( tag::cells_of_dim, 1 );  }

inline CellIterator Mesh::iter_over ( const tag::Segments &, const tag::ForcePositive & ) const
{	return iter_over ( tag::cells_of_dim, 1, tag::force_positive );  }

inline CellIterator Mesh::iter_over ( const tag::Segments &, const tag::Reverse & ) const
{	return iter_over ( tag::cells_of_dim, 1, tag::reverse );  }

inline CellIterator Mesh::iter_over
( const tag::Segments &, const tag::ForcePositive &, const tag::Reverse & ) const
{	return iter_over ( tag::cells_of_dim, 1, tag::force_positive, tag::reverse );  }

inline CellIterator Mesh::iter_over
( const tag::Segments &, const tag::Reverse &, const tag::ForcePositive & ) const
{	return iter_over ( tag::cells_of_dim, 1, tag::force_positive, tag::reverse );  }


}  // namespace maniFEM
	
#endif
// ifndef MANIFEM_ITER_H
