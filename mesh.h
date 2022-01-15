
//   mesh.h  2022.01.15

//   This file is part of maniFEM, a C++ library for meshes and finite elements on manifolds.

//   Copyright 2019, 2020, 2021, 2022 Cristian Barbarosie cristian.barbarosie@gmail.com
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


#ifndef MANIFEM_MESH_H
#define MANIFEM_MESH_H

#include <iostream>
#include <vector>
#include <list>
#include <deque>
#include <map>
#include <memory>
#include "assert.h"

namespace maniFEM {

// tags are used to distinguish between different versions of the same
// function or method or constructor
// they create a sort of spoken language ...
	
namespace tag {  // see paragraph 11.3 in the manual

	// tags are listed in no particular order
	
	struct Positive { };  static const Positive positive;
	struct PositiveCell { };  static const PositiveCell positive_cell;
	struct IsPositive { };  static const IsPositive is_positive;
	struct IsNegative { };  static const IsNegative is_negative;
	struct ThisMeshIsPositive { };  static const ThisMeshIsPositive this_mesh_is_positive;
	struct ThisMeshIsNegative { };  static const ThisMeshIsNegative this_mesh_is_negative;
	struct RequireOrder { };  static const RequireOrder require_order;
	struct ReverseOrder { };  static const ReverseOrder reverse_order;
	struct ReverseOrderIfAny { };  static const ReverseOrderIfAny reverse_order_if_any;
	struct ReverseOf { };  static const ReverseOf reverse_of;
	struct Reversed { };  static const Reversed reversed;
	struct Backwards { };  static const Backwards backwards;
	                       static const Backwards backward;
  struct Around { };  static const Around around;
	struct ReverseEachCell { };  static const ReverseEachCell reverse_each_cell;
	struct HasNoReverse { };  static const HasNoReverse has_no_reverse;
	struct AsTheyAre { };  static const AsTheyAre as_they_are;
	struct NonExistent { };  static const NonExistent non_existent;
	struct BuildIfNotExists { };  static const BuildIfNotExists build_if_not_exists;
	struct SeenFrom { };  static const SeenFrom seen_from;
	struct Fuzzy { };  static const Fuzzy fuzzy;
	struct MayNotExist { };  static const MayNotExist may_not_exist;
	struct DoNotBuildCells { };  static const DoNotBuildCells do_not_build_cells;
	struct DoNotBother { };  static const DoNotBother do_not_bother;
	struct SurelyExists { };  static const SurelyExists surely_exists;
	struct CellsSurelyExist { };  static const CellsSurelyExist cells_surely_exist;
	struct OfDimension { };  static const OfDimension of_dim;
	                         static const OfDimension of_dimension;
	struct OfMaxDim { };  static const OfMaxDim of_max_dim;
	struct SameDim { };  static const SameDim same_dim;
	struct CellHasLowDim { }; static const CellHasLowDim cell_has_low_dim;
	struct MinusOne { };  static const MinusOne minus_one;
	struct GreaterThanOne { };  static const GreaterThanOne greater_than_one;
	struct MightBeOne { };  static const MightBeOne might_be_one;
	struct Oriented { };  static const Oriented oriented;
	struct NotOriented { };  static const NotOriented not_oriented;
	struct DeepCopyOf { };  static const DeepCopyOf deep_copy_of;
	struct BuildCellsIfNec { };  static const BuildCellsIfNec build_cells_if_necessary;
	struct BuildNewVertices { };  static const BuildNewVertices build_new_vertices;
	struct UseExistingVertices { };  static const UseExistingVertices use_existing_vertices;
	struct Progressive { };  static const Progressive progressive;
	struct StartAt { };  static const StartAt start_at;
	struct StopAt { };  static const StopAt stop_at;
	struct Towards { };  static const Towards towards;
	struct Boundary { };  static const Boundary boundary;
	struct BoundaryOf { };  static const BoundaryOf boundary_of;
	struct SizeMeshes { };  static const SizeMeshes size_meshes;
	struct BehindFace { };  static const BehindFace behind_face;
	struct InFrontOfFace { };  static const InFrontOfFace in_front_of_face;
	struct WithinMesh { };  static const WithinMesh within_mesh;
	struct Vertices { };  static const Vertices vertices;
	struct Segments { };  static const Segments segments;
	struct Over { };  static const Over over;
	struct OverVertices { };  static const OverVertices over_vertices;
	struct OverSegments { };  static const OverSegments over_segments;
	struct CellsOfDim { };  static const CellsOfDim cells_of_dim;
	struct CellsOfMaxDim { };  static const CellsOfMaxDim cells_of_max_dim;
	struct OverCells { };  static const OverCells over_cells;
	struct OverCellsOfDim { };  static const OverCellsOfDim over_cells_of_dim;
	struct OverCellsOfMaxDim { };  static const OverCellsOfMaxDim over_cells_of_max_dim;
	struct OverCellsOfReverseOf { }; static const OverCellsOfReverseOf over_cells_of_reverse_of;
	struct Vertex { };  static const Vertex vertex; static const Vertex point;
	struct Segment { };  static const Segment segment;
	struct DividedIn { };  static const DividedIn divided_in;
	struct Triangle { };  static const Triangle triangle;
	struct Quadrangle { };  static const Quadrangle quadrangle;
	                        static const Quadrangle quadrilateral;
	struct Parallelogram { };  static const Parallelogram parallelogram;
	struct Rectangle { };  static const Rectangle rectangle;
	struct Square { };  static const Square square;
	struct Pentagon { };  static const Pentagon pentagon;
	struct Hexagon { };  static const Hexagon hexagon;
	struct Tetrahedron { };  static const Tetrahedron tetrahedron;
	struct Hexahedron { };  static const Hexahedron hexahedron;
	struct WhoseBoundaryIs { };  static const WhoseBoundaryIs whose_bdry_is;
	                             static const WhoseBoundaryIs whose_boundary_is;
	struct WhoseCoreIs { };  static const WhoseCoreIs whose_core_is;
	struct ForcePositive { };  static const ForcePositive force_positive;
	struct HasSize { };  static const HasSize has_size;
	struct ReserveSize { };  static const ReserveSize reserve_size;
	struct Pretty { };  static const Pretty pretty;
	struct Adapt { };  static const Adapt adapt;
	struct Identify { };  static const Identify identify;
	struct With { };  static const With with;
	struct Within { };  static const Within within;
	struct Replace { };  static const Replace replace;
	struct By { };  static const By by;
	struct OfDegree { };  static const OfDegree of_degree;
	struct MeshIsBdry { };  static const MeshIsBdry mesh_is_bdry;
	struct MeshIsNotBdry { };  static const MeshIsNotBdry mesh_is_not_bdry;
	enum WithTriangles { with_triangles, not_with_triangles };
	struct Join { };  static const Join join;
	enum KeyForHook { tangent_vector, normal_vector, node_in_cloud };
	struct Onto { };  static const Onto onto;
	struct EntireManifold { };  static const EntireManifold entire_manifold;
	struct DesiredLength { };  static const DesiredLength desired_length;
	struct Orientation { };  static const Orientation orientation;
	enum OrientationChoice { random, intrinsic, inherent, not_provided };
	struct Import { };  static const Import import;
	struct Msh { };  static const Msh msh;
	struct Winding { };  static const Winding winding;
	struct Singular { };  static const Singular singular;
	struct Unfold { };  static const Unfold unfold;
	struct ReturnMapBetween { };  static const ReturnMapBetween return_map_between;
	struct OverRegion { };  static const OverRegion over_region;
	struct OneGenerator { };  static const OneGenerator one_generator;
	struct TwoGenerators { };  static const TwoGenerators two_generators;
	struct ShadowVertices { };  static const ShadowVertices shadow_vertices;

	struct Util
	{ template < class T > class Wrapper;
		class Core;
		class CellCore;  //  aka  class Cell::Core
		class MeshCore;  //  aka  class Mesh::Core
		class CellIterator;  // aka  class Mesh::Iterator
		class Action;  //  aka class Function::Action, aka class Manfold::Action
		// we define it in function.h because we need it for Function::MultiValued
		// but we prefer the user to see it as an attribute of class Manifold
		class InequalitySet;  //  aka  class Function::Inequality::Set
		// defined in function.f
		inline static size_t assert_diff ( const size_t a, const size_t b )
		{	assert ( a >= b );  return  a - b;  }
		template < typename X, typename Y > inline static Y assert_cast ( X x )
		#ifndef NDEBUG
		{	Y y = dynamic_cast < Y > (x);  assert (y);  return y;  }
		#else
		{	Y y = static_cast < Y > (x);  return y;  }
		#endif
		template < typename X, typename Y > inline static const Y assert_const_cast ( const X x )
		#ifndef NDEBUG
		{	const Y y = dynamic_cast < Y* > (x);  assert (y);  return y;  }
		#else
		{	const Y y = static_cast < Y* > (x);  return y;  }
		#endif
		static bool return_true ( );
		static bool return_false ( );
		static const double one_third, minus_one_third, one_sixth, minus_one_sixth,
		                    two_thirds, minus_two_thirds;
	};
	struct MayBeNull { };  static const MayBeNull may_be_null;
	struct SurelyNotNull { };  static const SurelyNotNull surely_not_null;
	struct FreshlyCreated { };  static const FreshlyCreated freshly_created;
	struct PreviouslyExisting { };  static const PreviouslyExisting previously_existing;
	struct ZeroWrappers { };  static const ZeroWrappers zero_wrappers;
	struct OldCore { };  static const OldCore old_core;
	struct NoNeedToDisposeOf { };  static const NoNeedToDisposeOf no_need_to_dispose_of;
	struct Empty { };  static const Empty empty;
	struct OneDummyWrapper { };  static const OneDummyWrapper one_dummy_wrapper;

	namespace local_functions { };

}  // end of namespace tag

//-----------------------------------------------------------------------------//


// wrappers act like shared pointers towards cores
// you may want to forbid copying for your cores

// cores should have virtual destructor

template < class core_type >	class tag::Util::Wrapper 

{	public :

	core_type * core;

	inline Wrapper ( const tag::Empty & ) : core { nullptr }  { }

	// this constructor is avoided (replaced by direct assignment of core) in iterator.cpp
	inline Wrapper ( const tag::WhoseCoreIs &, core_type * c, const tag::FreshlyCreated & )
	: core { c }  // c just constructed with nb_of_wrappers == 1, no need to increment
	{	assert ( c );  }
	
	inline Wrapper ( const tag::WhoseCoreIs &, core_type * c,
	                 const tag::PreviouslyExisting &, const tag::SurelyNotNull & )
	: core { c }
	{	assert ( c );  c->nb_of_wrappers ++;  }
	
	inline Wrapper ( const tag::WhoseCoreIs &, core_type * c,
	                 const tag::PreviouslyExisting &, const tag::MayBeNull )
	: core { c }
	{	if ( c ) c->nb_of_wrappers ++;  }
	
	inline ~Wrapper ( )
	{	this->dispose_core ( tag::may_be_null );  }

	inline Wrapper operator= ( const Wrapper & w )
	{	dispose_core ( tag::may_be_null );
		set_core ( w .core, tag::may_be_null );
		return *this;                          }
	
	inline void dispose_core ( const tag::MayBeNull & )
	{	if ( this->core )
			if ( this->core->dispose_query () )
				delete this->core;                 }
	
	inline void dispose_core ( const tag::SurelyNotNull & )
	{	assert ( this->core );
		if ( this->core->dispose_query () )
			delete this->core;                }
	
	class Inactive;

};  // end of class tag::Util::Wrapper


// class tag::Util::Core is intended to be used by constructing
// your own Core class which inherits from tag::Util::Core
// you will probably implement polymorphic cores
// we declare the destructors virtual, thus making the class polymorphic
// this is useful if you want to do things like
//   tag::Util::Core * c = pointer_to_your_own_type_of_core;
//   ... some code ...
//   delete c;
// also useful for dynamic_cast
// see e.g. item 7 in the book of Scott Meyers, Effective C++

// more often than not, cores are built with nb_of_wrappers = 1
// although they did not meet any wrapper yet (they are newborn)
// we assume that the 'new' command was issued by a wrapper
// who is waiting for them to see the light of day and is eager to embrace them

class tag::Util::Core

{	public :

	size_t nb_of_wrappers;

	inline Core ( const tag::OneDummyWrapper & ) : nb_of_wrappers { 1 } { }
	// a new core is often built through a wrapper
	// thus, we initialize with one instead of zero

	inline Core ( const tag::ZeroWrappers & ) : nb_of_wrappers { 0 } { }
	// in some cases, a core is built without a wrapper
	// like for negative cells
	
	virtual ~Core ( ) { }

  inline bool dispose_query ( );

	class DelegateDispose;
	class Inactive;
	static bool default_dispose_query ( tag::Util::Core::DelegateDispose * );
	static bool dispose_query_cell_with_reverse ( tag::Util::Core::DelegateDispose * );

};


class tag::Util::Core::DelegateDispose : public tag::Util::Core

// some classes need specialized 'dispose_query' method, e.g. cells having reverse

{	public :

	bool (*dispose_query_p) ( tag::Util::Core::DelegateDispose * )
	{ & tag::Util::Core::default_dispose_query };  // this is not the body of a function
	// it is an initialization of the pointer-to-function

	inline DelegateDispose ( bool (*disp) ( tag::Util::Core::DelegateDispose * ),
                           const tag::OneDummyWrapper &                        )
	// a new core is built through a wrapper; thus, we initialize with one instead of zero
	:	tag::Util::Core ( tag::one_dummy_wrapper ), dispose_query_p { disp } { }

	inline DelegateDispose ( bool (*disp) ( tag::Util::Core::DelegateDispose * ),
                           const tag::ZeroWrappers &                           )
	// in rare cases, a core is built without a wrapper, like for negative cells
	:	tag::Util::Core ( tag::zero_wrappers ), dispose_query_p { disp } { }
	
	virtual ~DelegateDispose ( ) { }

	inline bool dispose_query ( )
	{	return dispose_query_p ( this );  }

};  // end of class tag::Util::Core


template < class core_type >	class tag::Util::Wrapper<core_type>::Inactive
// does nothing
{	public :
	core_type * core;
	inline Inactive ( const tag::Empty & ) : core { nullptr } { }
	inline Inactive ( core_type * c ) : core { c } { }
	inline ~Inactive ( ) { }
};

class tag::Util::Core::Inactive
// does nothing
{	public :
	inline Inactive ( const tag::OneDummyWrapper & ) { }
	inline Inactive ( const tag::ZeroWrappers & ) { }
	inline ~Inactive ( ) { }
};

//-----------------------------------------------------------------------------//

class Cell;  class Mesh;
class Manifold;  class Function;

//-----------------------------------------------------------------------------//


class tag::Util::CellIterator  // aka  class Mesh::Iterator

// a thin wrapper around a CellIterator::Core, with most methods delegated to 'core'

// iterates over all cells of a given mesh (cells of given dimension)

// for maximum dimension (equal to the dimension of the mesh) returns oriented cells
// which may be positive or negative (tag::force_positive overrides this)
// for lower dimension, returns positive cells

// or, iterates over all cells around a given cell

// see paragraph 9.5 in the manual
	
{	public :

	class Core;
	
	std::unique_ptr < tag::Util::CellIterator::Core > core;

	inline CellIterator ( const tag::Util::CellIterator & it );
	// inline CellIterator ( tag::Util::CellIterator && it );

	inline CellIterator ( const tag::WhoseCoreIs &, tag::Util::CellIterator::Core * c )
	:	core { c }
	{	assert ( c );  }
	
	inline ~CellIterator ( ) { };
	
	inline tag::Util::CellIterator & operator= ( const tag::Util::CellIterator & it );
	// inline tag::Util::CellIterator & operator= ( tag::Util::CellIterator && it );

	inline void reset ( );	
	inline void reset ( const tag::StartAt &, const Cell & );	
	inline Cell operator* ( );
	inline tag::Util::CellIterator & operator++ ( );
	inline tag::Util::CellIterator & operator++ ( int );
	inline tag::Util::CellIterator & advance ( );
	inline bool in_range ( );

	struct Over
	{	class TwoVerticesOfSeg;
		class CellsOfConnectedOneDimMesh;
		class VerticesOfConnectedOneDimMesh;
		class SegmentsOfConnectedOneDimMesh;
		class CellsOfFuzzyMesh;
		class SegmentsAboveVertex;           };
	struct Around  {  class OneCell;  struct OneVertex { struct OfAnyCodim
	       {  class OverSegments;  class OverHDCells; };  };  };
	struct Adaptor  {  class ForcePositive;  };

};  // end of class CellIterator


//-----------------------------------------------------------------------------//
//-----------------  wrappers Cell and Mesh  ----------------------------------//
//-----------------------------------------------------------------------------//


// a cell of dimension zero is a point, see classes Cell::Positive::Vertex and Negative
// a cell of dimension one is a segment, see classes Cell::Positive::Segment and Negative
// for cells of dimension two or more, see classes Cell::Positive::HighDim and Negative
// a cell of dimension two may be a triangle, a quadrangle or some other polygon
// a cell of dimension three may be a tetrahedron, a cube or some other polyhedron
// cells of dimension four or higher may be constructed,
// but their usefulness is questionable

// cells may be positively or negatively oriented
// see class Cell::Positive and class Cell::Negative
// which is which depends only on the construction process
// the first one to be created is positive, the other one will be negative,
// created when we call the 'reverse' method

// a cell is mainly defined by its boundary (which is a mesh of lower dimension)
// the orientation of a cell is nothing more than an orientation of its boundary
// see the comments on orientation in class Mesh below

// see paragraphs 1.2, 9.1 and 11.4 in the manual

// I would very much prefer the name 'Cell::Core' instead of 'tag::Util::CellCore'
// but it was not possible ...

#ifdef MANIFEM_COLLECT_CM	

class Cell : public tag::Util::Wrapper < tag::Util::CellCore >

#else  // no MANIFEM_COLLECT_CM

class Cell : public tag::Util::Wrapper < tag::Util::CellCore > ::Inactive 

#endif  // MANIFEM_COLLECT_CM	

// a thin wrapper around a Cell::Core, with most methods delegated to 'core'

{	public :

	typedef tag::Util::CellCore Core;

	// static attributes at the end

	// Cell::Core * core  inherited from tag::Util::Wrapper < Cell::Core >

	inline Cell ( const tag::NonExistent & )
	#ifdef MANIFEM_COLLECT_CM	
	:	tag::Util::Wrapper < Cell::Core > ( tag::empty )
	#else  // no MANIFEM_COLLECT_CM	
	:	tag::Util::Wrapper < Cell::Core > ::Inactive ( tag::empty )
	#endif  // MANIFEM_COLLECT_CM	
	{ }

	inline Cell ( const tag::WhoseCoreIs &, Cell::Core * c, const tag::FreshlyCreated & )
	#ifdef MANIFEM_COLLECT_CM	
	:	tag::Util::Wrapper < Cell::Core > ( tag::whose_core_is, c, tag::freshly_created )
	#else  // no MANIFEM_COLLECT_CM	
	:	tag::Util::Wrapper < Cell::Core > ::Inactive ( c )
	#endif  // MANIFEM_COLLECT_CM	
	{ }

	inline Cell ( const tag::WhoseCoreIs &, Cell::Core * c,
	              const tag::PreviouslyExisting &, const tag::SurelyNotNull & )
	#ifdef MANIFEM_COLLECT_CM	
	:	tag::Util::Wrapper < Cell::Core > ( tag::whose_core_is, c,
		                                    tag::previously_existing, tag::surely_not_null )
	#else  // no MANIFEM_COLLECT_CM	
	:	tag::Util::Wrapper < Cell::Core > ::Inactive ( c )
	#endif  // MANIFEM_COLLECT_CM	
	{	}

	inline Cell ( const tag::WhoseCoreIs &, Cell::Core * c,
	              const tag::PreviouslyExisting &, const tag::MayBeNull & )
	#ifdef MANIFEM_COLLECT_CM	
	:	tag::Util::Wrapper < Cell::Core > ( tag::whose_core_is, c,
		                                    tag::previously_existing, tag::may_be_null )
	#else  // no MANIFEM_COLLECT_CM	
	:	tag::Util::Wrapper < Cell::Core > ::Inactive ( c )
	#endif  // MANIFEM_COLLECT_CM	
	{	}

	inline Cell ( ) : Cell ( tag::non_existent ) { assert ( false );  }
	
	inline Cell ( const Cell & c )
	: Cell ( tag::whose_core_is, c.core, tag::previously_existing, tag::may_be_null )
	{	}

	inline Cell ( Cell && c ) : Cell ( tag::whose_core_is, c.core, tag::freshly_created )
	// name of tag::freshly_created is not the best choice
	// could be tag::do_not_increment_nb_of_wrappers
	{	c .core = nullptr;  }

	inline Cell ( const tag::WhoseBoundaryIs &, Mesh & );
	inline Cell ( const tag::Vertex &, const tag::IsPositive & ispos = tag::is_positive );
	inline Cell ( const tag::Segment &, const Cell & A, const Cell & B );
	inline Cell ( const tag::Triangle &, const Cell & AB, const Cell & BC, const Cell & CA );
	inline Cell ( const tag::Quadrangle &, const Cell & AB, const Cell & BC,
                                         const Cell & CD, const Cell & DA );
	inline Cell ( const tag::Parallelogram &, const Cell & AB, const Cell & BC,
                                            const Cell & CD, const Cell & DA );
	inline Cell ( const tag::Rectangle &, const Cell & AB, const Cell & BC,
                                        const Cell & CD, const Cell & DA );
	inline Cell ( const tag::Square &, const Cell & AB, const Cell & BC,
                                     const Cell & CD, const Cell & DA );
	inline Cell ( const tag::Pentagon &, const Cell & AB, const Cell & BC,
                                       const Cell & CD, const Cell & DE, const Cell & EA );
	inline Cell ( const tag::Hexagon &, const Cell & AB, const Cell & BC,
                       const Cell & CD, const Cell & DE, const Cell & EF, const Cell & FA );

	inline Cell & operator= ( const Cell & c );
	inline Cell & operator= ( Cell && c );
	
	// we are still in class Cell
	
	// methods delegated to 'core'

	inline Cell reverse
	( const tag::BuildIfNotExists & build = tag::build_if_not_exists ) const;
	inline Cell reverse ( const tag::SurelyExists & ) const;
	inline Cell reverse ( const tag::MayNotExist & ) const;
	inline Mesh boundary () const;
	inline bool exists () const  { return this->core != nullptr;  }
	inline bool is_positive () const;
	inline Cell get_positive ();
	inline size_t dim () const;
	inline bool has_reverse () const;
	inline Cell tip () const;
	inline Cell base () const;

	inline bool belongs_to
		( const Mesh & msh, const tag::SameDim &, const tag::Oriented & o = tag::oriented ) const;
	inline bool belongs_to
		( const Mesh & msh, const tag::SameDim &, const tag::NotOriented & ) const;
	inline bool belongs_to
		( const Mesh & msh, const tag::CellHasLowDim &,
		  const tag::NotOriented & no = tag::not_oriented ) const;
	inline bool belongs_to ( const Mesh & msh, const tag::Oriented & ) const;
	inline bool belongs_to ( const Mesh & msh, const tag::NotOriented & ) const;
	inline bool belongs_to ( const Mesh & msh ) const;
	
	inline bool is_inner_to ( const Mesh & msh ) const;  // defined in iterator.h

	// method 'glue_on_bdry_of' is intensively used when building a mesh
	// it glues 'this' cell to the boundary of 'cll'
	inline void glue_on_bdry_of ( Cell & cll );
	inline void glue_on_bdry_of ( Cell & cll, const tag::DoNotBother & );
	
	// method 'cut_from_bdry_of' does the reverse : cuts 'this' cell from
	// the boundary of 'cll' - used mainly in remeshing
	inline void cut_from_bdry_of ( Cell & cll );
	inline void cut_from_bdry_of ( Cell & cll, const tag::DoNotBother & );
	
	// methods 'add_to_mesh' and 'remove_from_mesh' add/remove 'this' cell to/from the mesh 'msh'
	// if 'msh' is the boundary of some cell, methods 'glue_on_bdry_of'
	// and 'cut_from_bdry_of' (see above) should be used instead
	inline void add_to_mesh ( Mesh & msh );
	inline void add_to_mesh ( Mesh & msh, const tag::DoNotBother & );
	inline void remove_from_mesh ( Mesh & msh );
	inline void remove_from_mesh ( Mesh & msh, const tag::DoNotBother & );

	// we are still in class Cell

	inline void project ( ) const;  // both defined in manifold.h
	inline void project ( const tag::Onto &, const Manifold m ) const;

	class Winding;  // defined in function.h
	inline Winding winding ( ) const;  // defined in manifold.h
		
#ifndef NDEBUG
	inline void print_everything ( );
#endif

	// static int counter;

	// a list of functions to be called each time a new cell is created
	// two lists, in fact, one for positive cells the other for negative cells
	static std::vector < std::vector < void(*)(Cell::Core*,void*) > > init_pos_cell, init_neg_cell;
	// more data can be passed to the above functions by using
	static std::vector < std::vector < void* > > data_for_init_pos, data_for_init_neg;
	// probably through a 'this' pointer
	
	struct field_to_meshes
	{	short int counter_pos;
		short int counter_neg;
		std::list<Cell>::iterator where;
		inline field_to_meshes ( short int i, short int j )
		:	counter_pos {i}, counter_neg {j} { }
		inline field_to_meshes ( short int i, short int j,
		                         std::list < Cell > ::iterator w )
		:	counter_pos {i}, counter_neg {j}, where {w} { }              };

	struct field_to_meshes_same_dim
	{	short int sign;
		std::list < Cell > ::iterator where;
		inline field_to_meshes_same_dim ( short int i ) : sign {i} { }
		inline field_to_meshes_same_dim
		( short int i, std::list < Cell > ::iterator w )
		:	sign {i}, where {w} { }                                       };

  class Positive;  class Negative;
	
	// any way to rewrite the names below as Positive::Vertex and such ?
	// I get an error because Cell::Core::add_to_seg expects a PositiveSegment ...
	class PositiveVertex;  class PositiveNotVertex;
	class PositiveSegment;  class PositiveHighDim;
	class NegativeVertex;  class NegativeNotVertex;
	class NegativeSegment;  class NegativeHighDim;

	class Numbering;  // a callable object with a Cell as argument, returning a size_t
	
}; // end of  class Cell


inline bool operator== ( const Cell & c1, const Cell & c2 )
{	return c1.core == c2.core;  }

inline bool operator!= ( const Cell & c1, const Cell & c2 )
{	return c1.core != c2.core;  }

inline bool operator< ( const Cell & c1, const Cell & c2 )
{	return c1.core < c2.core;  }

//-----------------------------------------------------------------------------//
//-----------------------------------------------------------------------------//


// roughly speaking, a mesh is a collection of cells of the same dimension

// an orientation of the mesh
// is nothing more than an orientation of each of its cells (of maximum dimension)
// but these orientations cannot be arbitrary, they must be compatible
// in the sense that a face common to two cells must be seen as positive
// from one of the cells and as negative from the other cell

// as a consequence, for connected meshes there are only two possible orientations
// although nothing prevents a mesh to be disconnected, we only allow for two orientations
// which is which depends only on the construction process
// the first one to be created is positive, the other one will be negative,
// created when we call the 'reverse' method

// negative meshes will appear mostly as boundaries of negative cells
// that's why we do not store their core in the computer's memory
// wrappers for negative meshes are built on-the-fly, e.g. in method Cell::boundary
// their core points to the (positive) reverse, the only difference is in 'is_pos'

// see paragraphs 1.2, 9.1 and 11.4 in the manual

// I would very much prefer the name 'Mesh::Core' instead of 'tag::Util::MeshCore'
// but it was not possible ...

#ifdef MANIFEM_COLLECT_CM	

class Mesh : public tag::Util::Wrapper < tag::Util::MeshCore >

#else  // no MANIFEM_COLLECT_CM	

class Mesh : public tag::Util::Wrapper < tag::Util::MeshCore > ::Inactive

#endif  // MANIFEM_COLLECT_CM	

// a thin wrapper around a Mesh::Core, with most methods delegated to 'core'

{	public :
	
	typedef tag::Util::MeshCore Core;
	typedef tag::Util::CellIterator Iterator;

	// Mesh::Core * core  inherited from tag::Util::Wrapper < Mesh::Core >

	// there are no cores for negative meshes
	// instead, we keep here a pointer to the direct (positive) core
	// thus, for meshes, 'core' points always to a positive Mesh::Core

	// this is how we distinguish between a positive mesh and a negative one
	bool (*is_pos) ();

	// low-level constructors :

	inline Mesh ( const tag::NonExistent & )
	#ifdef MANIFEM_COLLECT_CM	
	:	tag::Util::Wrapper < Mesh::Core > ( tag::empty ),
	#else  // no MANIFEM_COLLECT_CM	
	:	tag::Util::Wrapper < Mesh::Core > ::Inactive ( tag::empty ),
	#endif  // MANIFEM_COLLECT_CM	
		is_pos { & tag::Util::return_true }
	{ }

	inline Mesh ( const tag::WhoseCoreIs &, Mesh::Core * c, const tag::FreshlyCreated &,
	              const tag::IsPositive & ispos = tag::is_positive                       )
	// c is freshly created with tag::one_dummy-wrapper
	#ifdef MANIFEM_COLLECT_CM	
	: tag::Util::Wrapper < Mesh::Core > ( tag::whose_core_is, c, tag::freshly_created ),
	#else  // no MANIFEM_COLLECT_CM	
	:	tag::Util::Wrapper < Mesh::Core > ::Inactive ( c ),
	#endif  // MANIFEM_COLLECT_CM	
		is_pos { & tag::Util::return_true }
	{	}
	
	inline Mesh ( const tag::WhoseCoreIs &, Mesh::Core * c, const tag::PreviouslyExisting &,
                const tag::IsPositive & ispos = tag::is_positive                           )
	// c is either pre-existent or is freshly created with tag::zero_wrappers
	#ifdef MANIFEM_COLLECT_CM	
	: tag::Util::Wrapper < Mesh::Core > ( tag::whose_core_is, c,
	                                      tag::previously_existing, tag::surely_not_null ),
	#else  // no MANIFEM_COLLECT_CM	
	:	tag::Util::Wrapper < Mesh::Core > ::Inactive ( c ),
	#endif  // MANIFEM_COLLECT_CM	
		is_pos { & tag::Util::return_true }
	{	}
	
	// build a negative mesh from a positive one
	// without worrying whether reverse cells exist or not
	inline Mesh ( const tag::WhoseCoreIs &, Mesh::Core * c, const tag::FreshlyCreated &,
	              const tag::IsNegative &, const tag::DoNotBuildCells &                 )
	#ifdef MANIFEM_COLLECT_CM	
	: tag::Util::Wrapper < Mesh::Core > ( tag::whose_core_is, c,
	                                      tag::previously_existing, tag::surely_not_null ),
	#else  // no MANIFEM_COLLECT_CM	
	:	tag::Util::Wrapper < Mesh::Core > ::Inactive ( c ),
	#endif  // MANIFEM_COLLECT_CM	
		is_pos { & tag::Util::return_false }
	{	}
	
	// build a negative mesh from a positive one
	// without worrying whether reverse cells exist or not
	inline Mesh ( const tag::WhoseCoreIs &, Mesh::Core * c, const tag::PreviouslyExisting &,
	              const tag::IsNegative &, const tag::DoNotBuildCells &                     )
	#ifdef MANIFEM_COLLECT_CM	
	: tag::Util::Wrapper < Mesh::Core > ( tag::whose_core_is, c,
                                        tag::previously_existing, tag::surely_not_null ),
	#else  // no MANIFEM_COLLECT_CM	
	:	tag::Util::Wrapper < Mesh::Core > ::Inactive ( c ),
	#endif  // MANIFEM_COLLECT_CM	
		is_pos { & tag::Util::return_false }
	{	}

	inline Mesh ( const Mesh & msh )
	: Mesh ( tag::whose_core_is, msh.core, tag::previously_existing )
	{	this->is_pos = msh .is_pos;  }
		
	inline Mesh ( Mesh && msh )
	: Mesh ( tag::whose_core_is, msh.core, tag::freshly_created )
	// name of tag::freshly_created is not the best choice
	// could be tag::do_not_increment_nb_of_wrappers
	{	this->is_pos = msh .is_pos;
		msh.core = nullptr;        }
		
	// more elaborate (high-level) constructors :
	// they call one of the above, then manipulate the mesh
	
	// build a negative mesh from a positive one, assuming all cells have reverse :
	inline Mesh ( const tag::WhoseCoreIs &, Mesh::Core *, const tag::PreviouslyExisting &,
	              const tag::IsNegative &, const tag::CellsSurelyExist &                   );
	
	// build a negative mesh from a positive one, creating reverse cells if necessary :
	inline Mesh ( const tag::WhoseCoreIs &, Mesh::Core *, const tag::PreviouslyExisting &,
	              const tag::IsNegative &, const tag::BuildCellsIfNec &                    );

	inline Mesh ( const tag::DeepCopyOf &, const Mesh & msh );
	inline Mesh ( const tag::DeepCopyOf &, const Mesh & msh, const tag::Fuzzy & );

	inline Mesh ( const tag::Fuzzy &, const tag::OfDimension &, const size_t dim,
                const tag::IsPositive & ispos = tag::is_positive                );

	inline Mesh ( const tag::Fuzzy &, const tag::OfDimension &, const size_t dim,
                const tag::MinusOne &, const tag::IsPositive & ispos = tag::is_positive );

	// we are still in class Mesh
	
	// constructors below build a chain of n segment cells
	inline Mesh ( const tag::Segment &, const Cell & A, const Cell & B,
	              const tag::DividedIn &, const size_t n               );
	inline Mesh ( const tag::Segment &, const Cell & A, const Cell & B,
	              const tag::DividedIn &, const size_t n,
	              const tag::Winding &, const tag::Util::Action & );
	
	// constructors below build a triangular mesh from three sides
	// the number of divisions defined by the divisions of the sides (must be the same)
	inline Mesh ( const tag::Triangle &, const Mesh & AB, const Mesh & BC, const Mesh & CA );
	inline Mesh ( const tag::Triangle &, const Mesh & AB, const Mesh & BC, const Mesh & CA,
	              const tag::Winding &                                                     );
	inline Mesh ( const tag::Triangle &, const Mesh & AB, const Mesh & BC, const Mesh & CA,
	              const tag::Winding &, const tag::Singular &, const Cell & S              );
	
	// constructors below build a rectangular mesh from four sides
	// the number of divisions are already in the sides (must be the same for opposite sides)
	// if last argument is 'with_triangles', each rectangular cell will be cut in two triangles
	inline Mesh ( const tag::Quadrangle &, const Mesh & south, const Mesh & east,
	                                       const Mesh & north, const Mesh & west,
	              const tag::WithTriangles & wt = tag::not_with_triangles         );
	inline Mesh ( const tag::Parallelogram &, const Mesh & south, const Mesh & east,
	                                          const Mesh & north, const Mesh & west,
	              const tag::WithTriangles & wt = tag::not_with_triangles         )
	:	Mesh ( tag::quadrangle, south, east, north, west, wt )  { }
	inline Mesh ( const tag::Rectangle &, const Mesh & south, const Mesh & east,
	                                      const Mesh & north, const Mesh & west,
	              const tag::WithTriangles & wt = tag::not_with_triangles         )
	:	Mesh ( tag::quadrangle, south, east, north, west, wt )  { }
	inline Mesh ( const tag::Square &, const Mesh & south, const Mesh & east,
	                                   const Mesh & north, const Mesh & west,
	              const tag::WithTriangles & wt = tag::not_with_triangles         )
	:	Mesh ( tag::quadrangle, south, east, north, west, wt )  { }

	// same as above, but take into account winding segments
	// specific information about widing numbers is included in the four segments
	inline Mesh ( const tag::Quadrangle &, const Mesh & south, const Mesh & east,
	                                       const Mesh & north, const Mesh & west,
	              const tag::Winding &, const tag::WithTriangles & wt = tag::not_with_triangles );
	inline Mesh ( const tag::Quadrangle &, const Mesh & south, const Mesh & east,
	                                       const Mesh & north, const Mesh & west,
	              const tag::WithTriangles &, const tag::Winding &                  );
	inline Mesh ( const tag::Parallelogram &, const Mesh & south, const Mesh & east,
	                                          const Mesh & north, const Mesh & west,
	              const tag::Winding &, const tag::WithTriangles & wt = tag::not_with_triangles )
	:	Mesh ( tag::quadrangle, south, east, north, west, tag::winding, wt )  { }
	inline Mesh ( const tag::Parallelogram &, const Mesh & south, const Mesh & east,
	                                          const Mesh & north, const Mesh & west,
	              const tag::WithTriangles &, const tag::Winding &                  )
	:	Mesh ( tag::quadrangle, south, east, north, west, tag::winding, tag::with_triangles ) { }
	inline Mesh ( const tag::Rectangle &, const Mesh & south, const Mesh & east,
	                                      const Mesh & north, const Mesh & west,
	              const tag::Winding &, const tag::WithTriangles & wt = tag::not_with_triangles )
	:	Mesh ( tag::quadrangle, south, east, north, west, tag::winding, wt )  { }
	inline Mesh ( const tag::Rectangle &, const Mesh & south, const Mesh & east,
	                                      const Mesh & north, const Mesh & west,
	              const tag::WithTriangles &, const tag::Winding &                  )
	:	Mesh ( tag::quadrangle, south, east, north, west, tag::winding, tag::with_triangles ) { }
	inline Mesh ( const tag::Square &, const Mesh & south, const Mesh & east,
	                                   const Mesh & north, const Mesh & west,
	              const tag::Winding &, const tag::WithTriangles & wt = tag::not_with_triangles )
	:	Mesh ( tag::quadrangle, south, east, north, west, tag::winding, wt )  { }
	inline Mesh ( const tag::Square &, const Mesh & south, const Mesh & east,
	                                   const Mesh & north, const Mesh & west,
	              const tag::WithTriangles &, const tag::Winding &                  )
	:	Mesh ( tag::quadrangle, south, east, north, west, tag::winding, tag::with_triangles ) { }

	inline Mesh
	( const tag::Quadrangle &, const Cell & SW, const Cell & SE,
	  const Cell & NE, const Cell & NW, const size_t m, const size_t n,
	  const tag::WithTriangles & wt = tag::not_with_triangles           );
	// builds a rectangular mesh from four vertices
	// if last argument is 'with_triangles', each rectangular cell will be cut in two triangles

	inline Mesh ( const tag::Join &, const Mesh &, const Mesh & );
	inline Mesh ( const tag::Join &, const Mesh &, const Mesh &, const Mesh & );
	inline Mesh ( const tag::Join &, const Mesh &, const Mesh &, const Mesh &, const Mesh & );
	inline Mesh ( const tag::Join &,
	              const Mesh &, const Mesh &, const Mesh &, const Mesh &, const Mesh & );

	template < typename container >
	inline Mesh ( const tag::Join &, const container & l );

	template < typename container >
	static inline size_t join ( Mesh * const that, const container & l );

	template < typename container >
	static inline void join_meshes ( Mesh * const that, const container & l );
		
	inline Mesh ( const tag::Import &, const tag::Msh &, const std::string );

	void import_msh ( Mesh * that, const std::string filename );  // defined in global.cpp

	// we are still in class Mesh
	// constructors with tag::Progressive are defined in progressive.cpp
	
	Mesh ( const tag::Progressive &, const tag::DesiredLength &, const Function & length );

	Mesh ( const tag::Progressive &, const tag::EntireManifold &, Manifold manif,
         const tag::DesiredLength &, const Function & length                   );

	Mesh ( const tag::Progressive &, const tag::DesiredLength &, const Function & length,
	       const tag::Orientation &, const tag::OrientationChoice &                      );

	Mesh ( const tag::Progressive &, const tag::EntireManifold &, Manifold manif,
	       const tag::DesiredLength &, const Function & length,
	       const tag::Orientation &, const tag::OrientationChoice &              );

	Mesh ( const tag::Progressive &, const tag::Boundary &, Mesh interface,
	       const tag::DesiredLength &, const Function & length             );

	Mesh ( const tag::Progressive &, const tag::Boundary &, Mesh interface,
	       const tag::DesiredLength &, const Function & length,
				 const tag::Orientation &, const tag::OrientationChoice &        );
 
	Mesh ( const tag::Progressive &, const tag::Boundary &, Mesh interface,
	       const tag::StartAt &, const Cell & start,
	       const tag::Towards &, std::vector < double > normal,
	       const tag::DesiredLength &, const Function & length               );

	Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	       const tag::Towards &, std::vector < double > tangent,
	       const tag::StopAt &, const Cell & stop,
	       const tag::DesiredLength &, const Function & length                 );

	Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	       const tag::Towards &, std::vector < double > tangent,
	       const tag::DesiredLength &, const Function & length                  );

	Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	       const tag::StopAt &, const Cell & stop,
	       const tag::DesiredLength &, const Function & length                  );

	Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	       const tag::DesiredLength &, const Function & length                  );

	Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	       const tag::DesiredLength &, const Function & length,
	       const tag::Orientation &, const tag::OrientationChoice &           );

	Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	       const tag::StopAt &, const Cell & stop,
	       const tag::DesiredLength &, const Function & length,
	       const tag::Orientation &, const tag::OrientationChoice &           );

	inline Mesh & operator= ( const Mesh & c );
	inline Mesh & operator= ( Mesh && c );

	// we are still in class Mesh

	// methods 'fold' defined in global.cpp -- 'wrap' is a synonym
	// some of them take a mesh having as external boundary a parallelogram
	// and identify one or two pairs of opposite sides
	// see paragraphs 7.16 - 7.18 in the manual
	
	Mesh fold ( const tag::BuildNewVertices & );

	Mesh fold ( const tag::BuildNewVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
										size_t dim, std::map < Cell, Cell > & m                 );

	inline Mesh wrap ( const tag::BuildNewVertices & )
	{	return this->fold ( tag::build_new_vertices );  }

	inline Mesh wrap ( const tag::BuildNewVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
										size_t dim, std::map < Cell, Cell > & m                 )
	{	return this->fold ( tag::build_new_vertices, tag::return_map_between,
												tag::cells_of_dim, dim, m                         );  }

	Mesh fold ( const tag::UseExistingVertices & );

	Mesh fold ( const tag::UseExistingVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
										size_t dim, std::map < Cell, Cell > & m                 );
	
	inline Mesh wrap ( const tag::UseExistingVertices & )
	{	return this->fold ( tag::use_existing_vertices );  }

	inline Mesh wrap ( const tag::UseExistingVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
										size_t dim, std::map < Cell, Cell > & m                 )
	{	return this->fold ( tag::use_existing_vertices, tag::return_map_between,
												tag::cells_of_dim, dim, m                         );  }

	Mesh fold ( const tag::Identify &, const Mesh & msh1, const tag::With &, const Mesh & msh2,
	            const tag::BuildNewVertices &                                                  );

	Mesh fold ( const tag::Identify &, const Mesh & msh1, const tag::With &, const Mesh & msh2,
	            const tag::BuildNewVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
										size_t dim, std::map < Cell, Cell > & m                 );

	inline Mesh wrap
	( const tag::Identify &, const Mesh & msh1, const tag::With &, const Mesh & msh2,
	  const tag::BuildNewVertices &                                                   )
	{	return this->fold ( tag::identify, msh1, tag::with, msh2, tag::build_new_vertices );  }

	inline Mesh wrap
		( const tag::Identify &, const Mesh & msh1, const tag::With &, const Mesh & msh2,
	            const tag::BuildNewVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
										size_t dim, std::map < Cell, Cell > & m                 )
	{	return this->fold ( tag::identify, msh1, tag::with, msh2,
		                    tag::build_new_vertices, tag::return_map_between,
		                    tag::cells_of_dim, dim, m  );  }

	Mesh fold ( const tag::Identify &, const Mesh & msh1, const tag::With &, const Mesh & msh2,
	            const tag::UseExistingVertices &                                               );

	Mesh fold ( const tag::Identify &, const Mesh & msh1, const tag::With &, const Mesh & msh2,
	            const tag::UseExistingVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
										size_t dim, std::map < Cell, Cell > & m                 );

	inline Mesh wrap
	( const tag::Identify &, const Mesh & msh1, const tag::With &, const Mesh & msh2,
	  const tag::UseExistingVertices &                                                )
	{	return this->fold ( tag::identify, msh1, tag::with, msh2, tag::use_existing_vertices );  }

	inline Mesh wrap ( const tag::Identify &, const Mesh & msh1,
	                   const tag::With &, const Mesh & msh2,
	                   const tag::UseExistingVertices &,
	                   const tag::ReturnMapBetween &, const tag::CellsOfDim &,
	                   size_t dim, std::map < Cell, Cell > & m                 )
	{	return this->fold ( tag::identify, msh1, tag::with, msh2,
		                    tag::use_existing_vertices, tag::return_map_between,
	                      tag::cells_of_dim, dim, m                           );  }

	Mesh fold ( const tag::Identify &, const Mesh & msh1, const tag::With &, const Mesh & msh2,
	            const tag::Identify &, const Mesh & msh3, const tag::With &, const Mesh & msh4,
	            const tag::BuildNewVertices &                                                  );

	Mesh fold ( const tag::Identify &, const Mesh & msh1, const tag::With &, const Mesh & msh2,
	            const tag::Identify &, const Mesh & msh3, const tag::With &, const Mesh & msh4,
	            const tag::BuildNewVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
										size_t dim, std::map < Cell, Cell > & m                 );

	inline Mesh wrap
	( const tag::Identify &, const Mesh & msh1, const tag::With &, const Mesh & msh2,
	  const tag::Identify &, const Mesh & msh3, const tag::With &, const Mesh & msh4,
	  const tag::BuildNewVertices &                                                  )
	{	return this->fold ( tag::identify, msh1, tag::with, msh2,
		                    tag::identify, msh3, tag::with, msh4, tag::build_new_vertices );  }

	inline Mesh wrap
	( const tag::Identify &, const Mesh & msh1, const tag::With &, const Mesh & msh2,
	  const tag::Identify &, const Mesh & msh3, const tag::With &, const Mesh & msh4,
	  const tag::BuildNewVertices &, const tag::ReturnMapBetween &,
	  const tag::CellsOfDim &, size_t dim, std::map < Cell, Cell > & m                 )
	{	return this->fold ( tag::identify, msh1, tag::with, msh2,
		                    tag::identify, msh3, tag::with, msh4, tag::build_new_vertices,
		                    tag::return_map_between, tag::cells_of_dim, dim, m             );  }
	
	Mesh fold ( const tag::Identify &, const Mesh & msh1, const tag::With &, const Mesh & msh2,
	            const tag::Identify &, const Mesh & msh3, const tag::With &, const Mesh & msh4,
	            const tag::UseExistingVertices &                                               );
	
	Mesh fold ( const tag::Identify &, const Mesh & msh1, const tag::With &, const Mesh & msh2,
	            const tag::Identify &, const Mesh & msh3, const tag::With &, const Mesh & msh4,
	            const tag::UseExistingVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
										size_t dim, std::map < Cell, Cell > & m                 );
	
	inline Mesh wrap
	( const tag::Identify &, const Mesh & msh1, const tag::With &, const Mesh & msh2,
	  const tag::Identify &, const Mesh & msh3, const tag::With &, const Mesh & msh4,
	  const tag::UseExistingVertices &                                               )
	{	return this->fold ( tag::identify, msh1, tag::with, msh2,
		                    tag::identify, msh3, tag::with, msh4, tag::use_existing_vertices );  }

	Mesh fold ( const tag::Identify &, const Mesh & msh1, const tag::With &, const Mesh & msh2,
	            const tag::Identify &, const Mesh & msh3, const tag::With &, const Mesh & msh4,
	            const tag::Identify &, const Mesh & msh5, const tag::With &, const Mesh & msh6,
	            const tag::BuildNewVertices &                                                  );

	Mesh fold ( const tag::Identify &, const Mesh & msh1, const tag::With &, const Mesh & msh2,
	            const tag::Identify &, const Mesh & msh3, const tag::With &, const Mesh & msh4,
	            const tag::Identify &, const Mesh & msh5, const tag::With &, const Mesh & msh6,
	            const tag::BuildNewVertices &, const tag::ReturnMapBetween &,
	            const tag::CellsOfDim &, size_t dim, std::map < Cell, Cell > & m                 );

	inline Mesh wrap
	( const tag::Identify &, const Mesh & msh1, const tag::With &, const Mesh & msh2,
	  const tag::Identify &, const Mesh & msh3, const tag::With &, const Mesh & msh4,
	  const tag::Identify &, const Mesh & msh5, const tag::With &, const Mesh & msh6,
	  const tag::BuildNewVertices &                                                  )
	{	return this->fold ( tag::identify, msh1, tag::with, msh2,
		                    tag::identify, msh3, tag::with, msh4,
		                    tag::identify, msh5, tag::with, msh6, tag::build_new_vertices );  }

	Mesh fold ( const tag::Identify &, const Mesh & msh1, const tag::With &, const Mesh & msh2,
	            const tag::Identify &, const Mesh & msh3, const tag::With &, const Mesh & msh4,
	            const tag::Identify &, const Mesh & msh5, const tag::With &, const Mesh & msh6,
	            const tag::UseExistingVertices &                                               );
	
	Mesh fold ( const tag::Identify &, const Mesh & msh1, const tag::With &, const Mesh & msh2,
	            const tag::Identify &, const Mesh & msh3, const tag::With &, const Mesh & msh4,
	            const tag::Identify &, const Mesh & msh5, const tag::With &, const Mesh & msh6,
	            const tag::UseExistingVertices &,
									const tag::ReturnMapBetween &, const tag::CellsOfDim &,
										size_t dim, std::map < Cell, Cell > & m                 );
	
	inline Mesh wrap
	( const tag::Identify &, const Mesh & msh1, const tag::With &, const Mesh & msh2,
	  const tag::Identify &, const Mesh & msh3, const tag::With &, const Mesh & msh4,
	  const tag::Identify &, const Mesh & msh5, const tag::With &, const Mesh & msh6,
	  const tag::UseExistingVertices &                                               )
	{	return this->fold ( tag::identify, msh1, tag::with, msh2,
		                    tag::identify, msh3, tag::with, msh4,
		                    tag::identify, msh5, tag::with, msh6, tag::use_existing_vertices );  }

	// 'unfold' does the opposite : takes a mesh and unfolds it
	// defined in global.cpp
	// over a given region of the plane so we can export the resulting mesh in msh format
	
	Mesh unfold ( const tag::OverRegion &, const tag::Util::InequalitySet & constraints,
	              const tag::ReturnMapBetween &, const tag::CellsOfDim &,
	              size_t dim, std::map < Cell, std::pair < Cell, tag::Util::Action > > & mapping )
	const;

	Mesh unfold ( const std::vector < tag::Util::Action > &,
	              const tag::OverRegion &, const tag::Util::InequalitySet & constraints,
	              const tag::ReturnMapBetween &, const tag::CellsOfDim &,
	              size_t dim, std::map < Cell, std::pair < Cell, tag::Util::Action > > & mapping )
	const;
	
	Mesh unfold ( const tag::OverRegion &, const tag::Util::InequalitySet & constraints ) const;

	Mesh unfold ( const std::vector < tag::Util::Action > &,
	              const tag::OverRegion &, const tag::Util::InequalitySet & constraints ) const;

	// inline versions of unfold defined in function.h
	
	inline Mesh unfold ( const tag::OverRegion &,
	                     const tag::Util::InequalitySet & c1,
	                     const tag::Util::InequalitySet & c2 ) const;
	inline Mesh unfold ( const tag::OverRegion &,
	                     const tag::Util::InequalitySet & c1,
	                     const tag::Util::InequalitySet & c2,
	                     const tag::Util::InequalitySet & c3 ) const;
	inline Mesh unfold ( const tag::OverRegion &,
	                     const tag::Util::InequalitySet & c1,
	                     const tag::Util::InequalitySet & c2,
	                     const tag::Util::InequalitySet & c3,
	                     const tag::Util::InequalitySet & c4 ) const;
	inline Mesh unfold ( const tag::OverRegion &,
	                     const tag::Util::InequalitySet & c1,
	                     const tag::Util::InequalitySet & c2,
	              const tag::ReturnMapBetween &, const tag::CellsOfDim &,
	              size_t dim, std::map < Cell, std::pair < Cell, tag::Util::Action > > & mapping )
	const;
	inline Mesh unfold ( const tag::OverRegion &,
	                     const tag::Util::InequalitySet & c1,
	                     const tag::Util::InequalitySet & c2,
	                     const tag::Util::InequalitySet & c3,
	              const tag::ReturnMapBetween &, const tag::CellsOfDim &,
	              size_t dim, std::map < Cell, std::pair < Cell, tag::Util::Action > > & mapping )
	const;
	inline Mesh unfold ( const tag::OverRegion &,
	                     const tag::Util::InequalitySet & c1,
	                     const tag::Util::InequalitySet & c2,
	                     const tag::Util::InequalitySet & c3,
	                     const tag::Util::InequalitySet & c4,
	              const tag::ReturnMapBetween &, const tag::CellsOfDim &,
	              size_t dim, std::map < Cell, std::pair < Cell, tag::Util::Action > > & mapping )
	const;

	inline Mesh unfold ( const std::vector < tag::Util::Action > &, const tag::OverRegion &,
	                     const tag::Util::InequalitySet & c1,
	                     const tag::Util::InequalitySet & c2           ) const;
	inline Mesh unfold ( const std::vector < tag::Util::Action > &, const tag::OverRegion &,
	                     const tag::Util::InequalitySet & c1,
	                     const tag::Util::InequalitySet & c2,
	                     const tag::Util::InequalitySet & c3           ) const;
	inline Mesh unfold ( const std::vector < tag::Util::Action > &, const tag::OverRegion &,
	                     const tag::Util::InequalitySet & c1,
	                     const tag::Util::InequalitySet & c2,
	                     const tag::Util::InequalitySet & c3,
	                     const tag::Util::InequalitySet & c4           ) const;
	inline Mesh unfold ( const std::vector < tag::Util::Action > &, const tag::OverRegion &,
	                     const tag::Util::InequalitySet & c1,
	                     const tag::Util::InequalitySet & c2,
	              const tag::ReturnMapBetween &, const tag::CellsOfDim &,
	              size_t dim, std::map < Cell, std::pair < Cell, tag::Util::Action > > & mapping )
	const;
	inline Mesh unfold ( const std::vector < tag::Util::Action > &, const tag::OverRegion &,
	                     const tag::Util::InequalitySet & c1,
	                     const tag::Util::InequalitySet & c2,
	                     const tag::Util::InequalitySet & c3,
	              const tag::ReturnMapBetween &, const tag::CellsOfDim &,
	              size_t dim, std::map < Cell, std::pair < Cell, tag::Util::Action > > & mapping )
	const;
	inline Mesh unfold ( const std::vector < tag::Util::Action > &, const tag::OverRegion &,
	                     const tag::Util::InequalitySet & c1,
	                     const tag::Util::InequalitySet & c2,
	                     const tag::Util::InequalitySet & c3,
	                     const tag::Util::InequalitySet & c4,
	              const tag::ReturnMapBetween &, const tag::CellsOfDim &,
	              size_t dim, std::map < Cell, std::pair < Cell, tag::Util::Action > > & mapping )
	const;

	// we are still in class Mesh

	inline void copy_all_cells_to ( Mesh & msh ) const;

	inline bool exists () const  { return this->core != nullptr;  }
	inline bool is_positive () const  {  assert ( this->exists() );  return this->is_pos();  }
	inline size_t dim () const;
	inline Mesh reverse () const;

	inline size_t number_of ( const tag::Vertices & ) const;
	inline size_t number_of ( const tag::Segments & ) const;
	inline size_t number_of ( const tag::CellsOfDim &, const size_t d ) const;
	inline size_t number_of ( const tag::CellsOfMaxDim & ) const;

	inline Cell first_vertex ( ) const;
	inline Cell last_vertex ( ) const;
	inline Cell first_segment ( ) const;
	inline Cell last_segment ( ) const;

	// we are still in class Mesh
	
	inline Cell cell_in_front_of
	( const Cell face, const tag::SurelyExists & se = tag::surely_exists ) const;

	inline Cell cell_behind
	( const Cell face, const tag::SurelyExists & se = tag::surely_exists ) const;

	inline Cell cell_in_front_of ( const Cell face, const tag::MayNotExist & ) const;

	inline Cell cell_behind ( const Cell face, const tag::MayNotExist & ) const;

	// we are still in class Mesh
	// the eight methods below are only relevant for STSI meshes

	inline Cell cell_in_front_of
	( const Cell face, const tag::SeenFrom &, const Cell neighbour,
	  const tag::SurelyExists & se = tag::surely_exists              ) const;

	inline Cell cell_in_front_of
	( const Cell face, const tag::SeenFrom &, const Cell neighbour,
	  const tag::MayNotExist &                                       ) const;

	inline Cell cell_behind
	( const Cell face, const tag::SeenFrom &, const Cell neighbour,
	  const tag::SurelyExists & se = tag::surely_exists              ) const;

	inline Cell cell_behind
	( const Cell face, const tag::SeenFrom &, const Cell neighbour,
	  const tag::MayNotExist &                                       ) const;

	inline void closed_loop ( const Cell & ver );
	// for connected one-dim meshes, set both first_ver and last_ver to 'ver'

	inline void closed_loop ( const Cell & ver, size_t );
	// for connected one-dim meshes, set both first_ver and last_ver to 'ver'
	// and number of segments

	// we are still in class Mesh
	
	void baricenter ( const Cell & ver );
	void baricenter ( const Cell & ver, const tag::Winding & );
	void baricenter ( const Cell & ver, const tag::Winding &,
                    const tag::ShadowVertices &, const std::vector < Cell > & vec_cll );

	// iterators defined in iterator.h
	inline Mesh::Iterator iterator ( const tag::OverVertices & ) const;
	inline Mesh::Iterator iterator ( const tag::OverVertices &, const tag::AsTheyAre & ) const;
	inline Mesh::Iterator iterator ( const tag::OverVertices &, const tag::RequireOrder & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverVertices &, const tag::AsTheyAre &, const tag::RequireOrder & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverVertices &, const tag::RequireOrder &, const tag::AsTheyAre & ) const;
	inline Mesh::Iterator iterator ( const tag::OverVertices &, const tag::Backwards & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverVertices &, const tag::AsTheyAre &, const tag::Backwards & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverVertices &, const tag::Backwards &, const tag::AsTheyAre & ) const;
	inline Mesh::Iterator iterator ( const tag::OverVertices &, const tag::ForcePositive & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverVertices &, const tag::ForcePositive &, const tag::RequireOrder & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverVertices &, const tag::RequireOrder &, const tag::ForcePositive & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverVertices &, const tag::ForcePositive &, const tag::Backwards & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverVertices &, const tag::Backwards &, const tag::ForcePositive & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::RequireOrder & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverVertices &, const tag::RequireOrder &,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Backwards &    ) const;
	inline Mesh::Iterator iterator
	( const tag::OverVertices &, const tag::Backwards &,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const;

	// we are still in class Mesh

	inline Mesh::Iterator iterator ( const tag::Over &, const tag::Vertices & ) const
	{	return this->iterator ( tag::over_vertices );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Vertices &, const tag::AsTheyAre & ) const
	{	return this->iterator ( tag::over_vertices, tag::as_they_are );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Vertices &, const tag::RequireOrder & ) const
	{	return this->iterator ( tag::over_vertices, tag::require_order );  }
	inline Mesh::Iterator iterator ( const tag::Over &, const tag::Vertices &,
                                 const tag::AsTheyAre &, const tag::RequireOrder & ) const
	{	return this->iterator ( tag::over_vertices, tag::as_they_are, tag::require_order );  }
	inline Mesh::Iterator iterator ( const tag::Over &, const tag::Vertices &,
                                 const tag::RequireOrder &, const tag::AsTheyAre & ) const
	{	return this->iterator ( tag::over_vertices, tag::require_order, tag::as_they_are );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Vertices &, const tag::Backwards & ) const
	{	return this->iterator ( tag::over_vertices, tag::backwards );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Vertices &, const tag::AsTheyAre &, const tag::Backwards & ) const
	{	return this->iterator ( tag::over_vertices, tag::as_they_are, tag::backwards );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Vertices &, const tag::Backwards &, const tag::AsTheyAre & ) const
	{	return this->iterator ( tag::over_vertices, tag::backwards, tag::as_they_are );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Vertices &, const tag::ForcePositive & ) const
	{	return this->iterator ( tag::over_vertices, tag::force_positive );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Vertices &,
	  const tag::ForcePositive &, const tag::RequireOrder & ) const
	{	return this->iterator ( tag::over_vertices, tag::force_positive, tag::require_order );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Vertices &,
	  const tag::RequireOrder &, const tag::ForcePositive & ) const
	{	return this->iterator ( tag::over_vertices, tag::require_order, tag::force_positive );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Vertices &,
	  const tag::ForcePositive &, const tag::Backwards & ) const
	{	return this->iterator ( tag::over_vertices, tag::force_positive, tag::backwards );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Vertices &,
	  const tag::Backwards &, const tag::ForcePositive & ) const
	{	return this->iterator ( tag::over_vertices, tag::backwards, tag::force_positive );  }
	inline Mesh::Iterator iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const;
	inline Mesh::Iterator iterator ( const tag::Over &, const tag::Vertices &,
                                 const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
	{	return this->iterator
			( tag::over_vertices, tag::reverse_each_cell, tag::do_not_build_cells );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Vertices &,
    const tag::ReverseEachCell &, const tag::DoNotBuildCells &, const tag::RequireOrder & ) const
	{	return this->iterator ( tag::over_vertices,
			tag::reverse_each_cell, tag::do_not_build_cells, tag::require_order );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Vertices &,
    const tag::RequireOrder &, const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
	{	return this->iterator ( tag::over_vertices,
			tag::require_order, tag::reverse_each_cell, tag::do_not_build_cells );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Vertices &,
    const tag::ReverseEachCell &, const tag::DoNotBuildCells &, const tag::Backwards & ) const
	{	return this->iterator
			( tag::over_vertices, tag::reverse_each_cell, tag::do_not_build_cells, tag::backwards );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Vertices &, const tag::Backwards &,
    const tag::ReverseEachCell &, const tag::DoNotBuildCells &       ) const
	{	return this->iterator
			( tag::over_vertices, tag::backwards, tag::reverse_each_cell, tag::do_not_build_cells );  }
	
	// we are still in class Mesh
	
	inline Mesh::Iterator iterator ( const tag::OverSegments & ) const;
	inline Mesh::Iterator iterator ( const tag::OverSegments &, const tag::AsTheyAre & ) const;
	inline Mesh::Iterator iterator ( const tag::OverSegments &, const tag::RequireOrder & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::RequireOrder & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverSegments &, const tag::RequireOrder &, const tag::AsTheyAre & ) const;
	inline Mesh::Iterator iterator ( const tag::OverSegments &, const tag::Backwards & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::Backwards & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverSegments &, const tag::Backwards &, const tag::AsTheyAre & ) const;
	inline Mesh::Iterator iterator ( const tag::OverSegments &, const tag::ForcePositive & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::RequireOrder & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverSegments &, const tag::RequireOrder &, const tag::ForcePositive & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::Backwards & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverSegments &, const tag::Backwards &, const tag::ForcePositive & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::RequireOrder & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverSegments &, const tag::RequireOrder &,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Backwards &    ) const;
	inline Mesh::Iterator iterator
	( const tag::OverSegments &, const tag::Backwards &,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const;

	// we are still in class Mesh

	inline Mesh::Iterator iterator ( const tag::Over &, const tag::Segments & ) const
	{	return this->iterator ( tag::over_segments );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Segments &, const tag::AsTheyAre & ) const
	{	return this->iterator ( tag::over_segments, tag::as_they_are );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Segments &, const tag::RequireOrder & ) const
	{	return this->iterator ( tag::over_segments, tag::require_order );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Segments &,
	  const tag::AsTheyAre &, const tag::RequireOrder & ) const
	{	return this->iterator ( tag::over_segments, tag::as_they_are, tag::require_order );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Segments &,
	  const tag::RequireOrder &, const tag::AsTheyAre & ) const
	{	return this->iterator ( tag::over_segments, tag::require_order, tag::as_they_are );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Segments &, const tag::Backwards & ) const
	{	return this->iterator ( tag::over_segments, tag::backwards );  }
	inline Mesh::Iterator iterator ( const tag::Over &, const tag::Segments &,
                                 const tag::AsTheyAre &, const tag::Backwards & ) const
	{	return this->iterator ( tag::over_segments, tag::as_they_are, tag::backwards );  }
	inline Mesh::Iterator iterator ( const tag::Over &, const tag::Segments &,
                                 const tag::Backwards &, const tag::AsTheyAre & ) const
	{	return this->iterator ( tag::over_segments, tag::backwards, tag::as_they_are );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Segments &, const tag::ForcePositive & ) const
	{	return this->iterator ( tag::over_segments, tag::force_positive );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Segments &,
	  const tag::ForcePositive &, const tag::RequireOrder & ) const
	{	return this->iterator ( tag::over_segments, tag::force_positive, tag::require_order );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Segments &,
	  const tag::RequireOrder &, const tag::ForcePositive & ) const
	{	return this->iterator ( tag::over_segments, tag::require_order, tag::force_positive );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Segments &,
	  const tag::ForcePositive &, const tag::Backwards & ) const
	{	return this->iterator ( tag::over_segments, tag::force_positive, tag::backwards );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Segments &,
	  const tag::Backwards &, const tag::ForcePositive & ) const
	{	return this->iterator ( tag::over_segments, tag::backwards, tag::force_positive );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Segments &,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
	{	return this->iterator
			( tag::over_segments, tag::reverse_each_cell, tag::do_not_build_cells );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Segments &,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells &, const tag::RequireOrder & ) const
	{	return this->iterator ( tag::over_segments, tag::reverse_each_cell,
	                          tag::do_not_build_cells, tag::require_order );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Segments &, const tag::RequireOrder &,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells &          ) const
	{	return this->iterator ( tag::over_segments, tag::require_order,
	                          tag::reverse_each_cell, tag::do_not_build_cells );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Segments &,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells &, const tag::Backwards & ) const
	{	return this->iterator
			( tag::over_segments, tag::reverse_each_cell, tag::do_not_build_cells, tag::backwards );  }
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::Segments &, const tag::Backwards &,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
	{	return this->iterator
			( tag::over_segments, tag::backwards, tag::reverse_each_cell, tag::do_not_build_cells );  }
	
	// we are still in class Mesh
	
	inline Mesh::Iterator iterator ( const tag::OverCellsOfDim &, const size_t ) const;
	inline Mesh::Iterator iterator ( const tag::Over &, const tag::CellsOfDim &, const size_t d ) const
	{	return this->iterator ( tag::over_cells_of_dim, d );  }
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre & ) const;
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::CellsOfDim &, const size_t d, const tag::AsTheyAre & ) const
	{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are );  }
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::RequireOrder & ) const;
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::CellsOfDim &, const size_t d, const tag::RequireOrder & ) const
	{	return this->iterator ( tag::over_cells_of_dim, d, tag::require_order );  }
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfDim &, const size_t,
	 const tag::AsTheyAre &, const tag::RequireOrder & ) const;
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::CellsOfDim &, const size_t d,
	  const tag::AsTheyAre &, const tag::RequireOrder &          ) const
	{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are, tag::require_order );  }
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::RequireOrder &, const tag::AsTheyAre & ) const;
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::CellsOfDim &, const size_t d,
	  const tag::RequireOrder &, const tag::AsTheyAre &          ) const
	{	return this->iterator ( tag::over_cells_of_dim, d, tag::require_order, tag::as_they_are );  }
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::Backwards & ) const;
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::CellsOfDim &, const size_t d, const tag::Backwards & ) const
	{	return this->iterator ( tag::over_cells_of_dim, d, tag::backwards );  }
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::AsTheyAre &, const tag::Backwards & ) const;
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::CellsOfDim &, const size_t d,
	  const tag::AsTheyAre &, const tag::Backwards &             ) const
	{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are, tag::backwards );  }
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::Backwards &, const tag::AsTheyAre & ) const;
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::CellsOfDim &, const size_t d,
	  const tag::Backwards &, const tag::AsTheyAre &             ) const
	{	return this->iterator ( tag::over_cells_of_dim, d, tag::backwards, tag::as_they_are );  }
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive & ) const;
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::CellsOfDim &, const size_t d, const tag::ForcePositive & ) const
	{	return this->iterator ( tag::over_cells_of_dim, d, tag::force_positive );  }
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::RequireOrder &                                             ) const;
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::CellsOfDim &, const size_t d,
	  const tag::ForcePositive &, const tag::RequireOrder &      ) const
	{	return this->iterator ( tag::over_cells_of_dim, d,
	                          tag::force_positive, tag::require_order );  }
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::RequireOrder &,
	  const tag::ForcePositive &                                           ) const;
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::CellsOfDim &, const size_t d,
	  const tag::RequireOrder &, const tag::ForcePositive &      ) const
	{	return this->iterator ( tag::over_cells_of_dim, d,
	                          tag::require_order, tag::force_positive );  }
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::Backwards &                                             ) const;
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::CellsOfDim &, const size_t d,
	  const tag::ForcePositive &, const tag::Backwards &         ) const
	{	return this->iterator ( tag::over_cells_of_dim, d, tag::force_positive, tag::backwards );  }
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::Backwards &,
	  const tag::ForcePositive &                                        ) const;
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::CellsOfDim &, const size_t d,
	  const tag::Backwards &, const tag::ForcePositive &         ) const
	{	return this->iterator ( tag::over_cells_of_dim, d, tag::backwards, tag::force_positive );  }
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfDim &, const size_t,
		const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const;
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::CellsOfDim &, const size_t d,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
	{	return this->iterator ( tag::over_cells_of_dim, d,
	                          tag::reverse_each_cell, tag::do_not_build_cells );  }
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::RequireOrder &                 ) const;
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::CellsOfDim &, const size_t d,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells &, const tag::RequireOrder & ) const
	{	return this->iterator ( tag::over_cells_of_dim, d,
	                    tag::reverse_each_cell, tag::do_not_build_cells, tag::require_order );  }
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::RequireOrder &,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells &           ) const;
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::CellsOfDim &, const size_t d, const tag::RequireOrder &,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells &                            ) const
	{	return this->iterator ( tag::over_cells_of_dim, d, tag::require_order,
		                        tag::reverse_each_cell, tag::do_not_build_cells );  }
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Backwards &                    ) const;
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::CellsOfDim &, const size_t d,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells &, const tag::Backwards & ) const
	{	return this->iterator ( tag::over_cells_of_dim, d,
		                        tag::reverse_each_cell, tag::do_not_build_cells, tag::backwards );  }
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::Backwards &,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells &        ) const;
	inline Mesh::Iterator iterator
	( const tag::Over &, const tag::CellsOfDim &, const size_t d,
	  const tag::Backwards &, const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
	{	return this->iterator ( tag::over_cells_of_dim, d,
	                    tag::backwards, tag::reverse_each_cell, tag::do_not_build_cells );  }

	// we are still in class Mesh
	
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t d ) const
	{	return this->iterator ( tag::over_cells_of_dim, d );  }
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t d,
	  const tag::AsTheyAre &                                           ) const
	{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are );  }
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t d,
	  const tag::RequireOrder &                                        ) const
	{	return this->iterator ( tag::over_cells_of_dim, d, tag::require_order );  }
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t d,
	  const tag::AsTheyAre &, const tag::RequireOrder &                ) const
	{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are, tag::require_order );  }
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t d,
	  const tag::RequireOrder &, const tag::AsTheyAre &                ) const
	{	return this->iterator ( tag::over_cells_of_dim, d, tag::require_order, tag::as_they_are );  }
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t d,
	  const tag::Backwards &                                           ) const
	{	return this->iterator ( tag::over_cells_of_dim, d, tag::backwards );  }
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t d,
	  const tag::AsTheyAre &, const tag::Backwards &                   ) const
	{	return this->iterator ( tag::over_cells_of_dim, d, tag::as_they_are, tag::backwards );  }
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t d,
	  const tag::Backwards &, const tag::AsTheyAre &                   ) const
	{	return this->iterator ( tag::over_cells_of_dim, d, tag::backwards, tag::as_they_are );  }
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t d,
	  const tag::ForcePositive &                                       ) const
	{	return this->iterator ( tag::over_cells_of_dim, d, tag::force_positive );  }
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t d,
	  const tag::ForcePositive &, const tag::RequireOrder &            ) const
	{	return this->iterator
			( tag::over_cells_of_dim, d, tag::force_positive, tag::require_order );  }
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t d,
	  const tag::RequireOrder &, const tag::ForcePositive &            ) const
	{	return this->iterator
			( tag::over_cells_of_dim, d, tag::require_order, tag::force_positive );  }
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t d,
	  const tag::ForcePositive &, const tag::Backwards &               ) const
	{	return this->iterator ( tag::over_cells_of_dim, d, tag::force_positive, tag::backwards );  }
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t d,
	  const tag::Backwards &, const tag::ForcePositive &               ) const
	{	return this->iterator ( tag::over_cells_of_dim, d, tag::backwards, tag::force_positive );  }
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t d,
		const tag::ReverseEachCell &, const tag::DoNotBuildCells &       ) const
	{	return this->iterator ( tag::over_cells_of_dim, d,
		                        tag::reverse_each_cell, tag::do_not_build_cells );  }
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t d,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells &, const tag::RequireOrder & ) const
	{	return this->iterator ( tag::over_cells_of_dim, d, tag::reverse_each_cell,
		                        tag::do_not_build_cells, tag::require_order       );  }
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t d,
	  const tag::RequireOrder &, const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
	{	return this->iterator ( tag::over_cells_of_dim, d,
	                    tag::require_order, tag::reverse_each_cell, tag::do_not_build_cells );  }
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t d,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells &, const tag::Backwards & ) const
	{	return this->iterator ( tag::over_cells_of_dim, d, tag::reverse_each_cell,
														tag::do_not_build_cells, tag::backwards           );  }
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t d,
	  const tag::Backwards &, const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const
	{	return this->iterator ( tag::over_cells_of_dim, d, tag::backwards,
														tag::reverse_each_cell, tag::do_not_build_cells );  }

	// we are still in class Mesh
	
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfMaxDim & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre & ) const;
	inline Mesh::Iterator iterator ( const tag::OverCellsOfMaxDim &, const tag::RequireOrder & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::RequireOrder & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::RequireOrder &, const tag::AsTheyAre & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::Backwards & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::Backwards & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::Backwards &, const tag::AsTheyAre & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseOrderIfAny & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::ReverseOrderIfAny & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseOrderIfAny &, const tag::AsTheyAre & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::RequireOrder &                                  ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::RequireOrder &,
	  const tag::ForcePositive &                                ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::Backwards &                                     ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::Backwards &,
	  const tag::ForcePositive &                             ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &                             ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseOrderIfAny &,
	  const tag::ForcePositive &                                     ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfMaxDim &,
		const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::RequireOrder &      ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::RequireOrder &,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Backwards &         ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::Backwards &,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::ReverseOrderIfAny & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseOrderIfAny &,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells &     ) const;

	// we are still in class Mesh
	
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfMaxDim & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &, const tag::AsTheyAre & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &, const tag::RequireOrder & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &,
	  const tag::AsTheyAre &, const tag::RequireOrder & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &,
	  const tag::RequireOrder &, const tag::AsTheyAre & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &, const tag::Backwards & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &,
	  const tag::AsTheyAre &, const tag::Backwards & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &,
	  const tag::Backwards &, const tag::AsTheyAre & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &, const tag::ReverseOrderIfAny & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &,
	  const tag::AsTheyAre &, const tag::ReverseOrderIfAny & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &,
	  const tag::ReverseOrderIfAny &, const tag::AsTheyAre & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &, const tag::ForcePositive & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &,
	  const tag::ForcePositive &, const tag::RequireOrder & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &,
	  const tag::RequireOrder &, const tag::ForcePositive & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &,
	  const tag::ForcePositive &, const tag::Backwards & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &,
	  const tag::Backwards &, const tag::ForcePositive & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &,
	  const tag::ForcePositive &, const tag::ReverseOrderIfAny & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &,
	  const tag::ReverseOrderIfAny &, const tag::ForcePositive & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &,
		const tag::ReverseEachCell &, const tag::DoNotBuildCells & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::RequireOrder &                     ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &, const tag::RequireOrder &,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells &               ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Backwards &                        ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &, const tag::Backwards &,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells &            ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::ReverseOrderIfAny &                ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &, const tag::ReverseOrderIfAny &,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells &                    ) const;

	// we are still in class Mesh

	inline Mesh::Iterator iterator
	( const tag::OverVertices &, const tag::Around &, const Cell & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverSegments &, const tag::Around &, const Cell & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::Around &, const Cell & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverSegments &, const tag::AsTheyAre &,
	  const tag::BuildCellsIfNec &, const tag::Around &, const Cell & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverSegments &, const tag::DoNotBuildCells &,
	  const tag::Around &, const Cell &                        ) const;
	inline Mesh::Iterator iterator
	( const tag::OverSegments &, const tag::AsTheyAre &,
	  const tag::DoNotBuildCells &, const tag::Around &, const Cell & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::Around &, const Cell &                              ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
		const tag::AsTheyAre &, const tag::Around &, const Cell &      ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, const Cell &                              ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::DoNotBuildCells &, const tag::Around &, const Cell & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, const Cell &                              ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::Around &, const Cell & ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::Around &, const Cell &                                 ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, const Cell &                    ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::DoNotBuildCells &,
	  const tag::Around &, const Cell &                                       ) const;
	inline Mesh::Iterator iterator
	( const tag::OverCellsOfDim, const size_t,
	  const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, const Cell &                    ) const;

	// we are still in class Mesh
	// methods  draw_ps  defined in global.cpp

	void draw_ps ( std::string file_name ) const;
	void draw_ps ( std::string file_name,
	               const tag::Unfold &, const std::vector < tag::Util::Action > &,
	               const tag::OverRegion &, const tag::Util::InequalitySet & constraints ) const;
	void draw_ps ( std::string file_name,
	               const tag::Unfold &, const tag::OverRegion &,
	               const tag::Util::InequalitySet & constraints ) const;
	void draw_ps ( std::string file_name,
	               const tag::Unfold &, const tag::OneGenerator &, const tag::OverRegion &,
	               const tag::Util::InequalitySet & constraints                            ) const;
	void draw_ps ( std::string file_name,
	               const tag::Unfold &, const tag::TwoGenerators &, const tag::OverRegion &,
	               const tag::Util::InequalitySet & constraints                            ) const;

	// inline versions of draw_ps defined in function.h

	inline void draw_ps ( std::string file_name,
	           const tag::Unfold &, const std::vector < tag::Util::Action > &,
	           const tag::OverRegion &, const tag::Util::InequalitySet & c1,
	                                    const tag::Util::InequalitySet & c2   ) const;
	inline void draw_ps ( std::string file_name,
	           const tag::Unfold &, const std::vector < tag::Util::Action > &,
	           const tag::OverRegion &, const tag::Util::InequalitySet & c1,
	                                    const tag::Util::InequalitySet & c2,
	                                    const tag::Util::InequalitySet & c3    ) const;
	inline void draw_ps ( std::string file_name,
	           const tag::Unfold &, const std::vector < tag::Util::Action > &,
	           const tag::OverRegion &, const tag::Util::InequalitySet & c1,
	                                    const tag::Util::InequalitySet & c2,
	                                    const tag::Util::InequalitySet & c3,
	                                    const tag::Util::InequalitySet & c4    ) const;
	inline void draw_ps ( std::string file_name,
	                      const tag::Unfold &, const tag::OverRegion &,
	                      const tag::Util::InequalitySet & c1,
	                      const tag::Util::InequalitySet & c2           ) const;
	inline void draw_ps ( std::string file_name,
	                      const tag::Unfold &, const tag::OverRegion &,
	                      const tag::Util::InequalitySet & c1,
	                      const tag::Util::InequalitySet & c2,
	                      const tag::Util::InequalitySet & c3           ) const;
	inline void draw_ps ( std::string file_name,
	                      const tag::Unfold &, const tag::OverRegion &,
	                      const tag::Util::InequalitySet & c1,
	                      const tag::Util::InequalitySet & c2,
	                      const tag::Util::InequalitySet & c3,
	                      const tag::Util::InequalitySet & c4           ) const;
	void draw_ps_3d ( std::string file_name ) const;
	
	// we are still in class Mesh
	// methods  export_msh  defined in global.cpp

	void export_msh ( std::string f, Cell::Numbering & ) const;
	void export_msh ( std::string f, std::map < Cell, size_t > & ) const;
	void export_msh ( std::string f ) const;
	
	// several versions of 'build' below are defined in global.cpp

	void build ( const tag::Segment &,  // builds a chain of n segment cells
	             const Cell & A, const Cell & B, const tag::DividedIn &, const size_t n );
	
	void build ( const tag::Segment &,  // builds a chain of n segment cells
	             const Cell & A, const Cell & B, const tag::DividedIn &, const size_t n,
	             const tag::Winding &, const tag::Util::Action &                        );

	void build ( const tag::Triangle &, const Mesh & AB, const Mesh & BC, const Mesh & CA );

	void build ( const tag::Triangle &, const Mesh & AB, const Mesh & BC, const Mesh & CA,
	             const tag::Winding &                                                     );

	void build ( const tag::Triangle &, const Mesh & AB, const Mesh & BC, const Mesh & CA,
	             const tag::Winding &, const tag::Singular &, const Cell & S              );
	
	void build ( const tag::Quadrangle &, const Mesh & south, const Mesh & east,
	             const Mesh & north, const Mesh & west, bool cut_rectangles_in_half );
	
	void build ( const tag::Quadrangle &, const Mesh & south, const Mesh & east,
	             const Mesh & north, const Mesh & west, bool cut_rectangles_in_half,
	             const tag::Winding &                                               );
	
	void build ( const tag::Quadrangle &, const Cell & SW, const Cell & SE,
	             const Cell & NE, const Cell & NW, const size_t m, const size_t n,
	             bool cut_rectangles_in_half                                      );
	
	#ifndef NDEBUG
	inline void print_everything ( );
	#endif
	
	// we keep here the topological dimension of the largest mesh we intend to build
	// see method 'set_max_dim' and paragraph 11.7 in the manual
	static size_t maximum_dimension_plus_one;

	inline static void set_max_dim ( const size_t d );
	// invoked at the beginning of mesh.cpp
	// see paragraph 11.7 in the manual
	
	struct Connected  {  class OneDim;  class HighDim;  };
	struct MultiplyConnected  {  class OneDim;  class HighDim; };
	class ZeroDim;  class NotZeroDim;  class Fuzzy;  class STSI;
	
}; // end of  class Mesh


inline bool operator== ( const Mesh & m1, const Mesh & m2 )
{	return m1 .core == m2 .core;  }


//-----------------------------------------------------------------------------//
//----------------      core Cells     ----------------------------------------//
//-----------------------------------------------------------------------------//


// I would very much prefer the name 'class Cell::Core' instead of 'tag::Util::CellCore'
// but it was not possible ...


#ifdef MANIFEM_COLLECT_CM	

class tag::Util::CellCore : public tag::Util::Core::DelegateDispose

#else  // no MANIFEM_COLLECT_CM	

class tag::Util::CellCore : public tag::Util::Core::Inactive

#endif  // MANIFEM_COLLECT_CM	

// abstract class
// specialized in Cell::Positive::{Vertex,Segment,HighDim} and
//                Cell::Negative::{Vertex,Segment,HighDim}

{	public :
	
	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM
	
	// we keep numeric values here :
	std::vector < double > double_heap;
	std::vector < size_t > size_t_heap;
	std::vector < short int > short_int_heap;
	// and heterogeneous information here :
	std::map < tag::KeyForHook, void * > hook;

	// if 'this' is a face of another cell and that cell belongs to some mesh msh,
	// cell_behind_within [msh] keeps that cell
	// see methods Mesh::cell_behind and Mesh::cell_in_front_of
	std::map < Mesh::Core *, Cell > cell_behind_within;
	// we use the Cell wrapper as pointer

	// we use a wrapper as pointer
	Cell reverse_attr;

	#ifndef NDEBUG
	std::string name;
	#endif

	inline CellCore ( const tag::OfDimension &, const size_t d,  // for positive cells
	                  const tag::HasNoReverse &, const tag::OneDummyWrapper & );

	inline CellCore ( const tag::OfDimension &, const size_t d,  // for negative cells
	           const tag::ReverseOf &, Cell::Core * direct_cell_p, const tag::OneDummyWrapper & );

	virtual ~CellCore ( )  // 
	{	}  // { Cell::counter--;  std::cout << Cell::counter << std::endl;  }

	// bool dispose_query ( )  defined by tag::Util::Core::DelegateDispose ifdef COLLECT_CM
	
	virtual bool is_positive ( ) const = 0;
	virtual Cell get_positive ( ) = 0;
	virtual size_t get_dim ( ) const = 0;

	virtual Cell tip ();
	virtual Cell base ();
	// execution forbidden
	// overridden by Cell::Positive::Segment and Cell::Negative::Segment

	virtual Mesh boundary ( ) = 0;

	virtual bool belongs_to
	( Mesh::Core *, const tag::CellHasLowDim &,
	  const tag::NotOriented & no = tag::not_oriented ) const = 0;
	virtual bool belongs_to
	( Mesh::Core *, const tag::SameDim &, const tag::Oriented & o = tag::oriented ) const = 0;
	bool belongs_to ( Mesh::Core *, const tag::SameDim &, const tag::NotOriented & ) const;

	// Method 'glue_on_bdry_of' is intensively used when building a mesh,
	// e.g. within factory functions in Cell class.
	// It glues 'this' cell to the boundary of 'cll'.
	inline void glue_on_bdry_of ( Cell::Core * cll )
	{	assert ( cll );
		cll->glue_on_my_bdry ( this );   }
	inline void glue_on_bdry_of ( Cell::Core * cll, const tag::DoNotBother & )
	{	assert ( cll );
		cll->glue_on_my_bdry ( this, tag::do_not_bother );   }
	// tag::do_not_bother is useful for a Mesh::Connected::OneDim

	// Method 'cut_from_bdry_of' does the reverse : cuts 'this' cell from
	// the boundary of 'cll'. Used mainly in remeshing.
	inline void cut_from_bdry_of ( Cell::Core * cll )
	{	cll->cut_from_my_bdry ( this );   }
	inline void cut_from_bdry_of ( Cell::Core * cll, const tag::DoNotBother & )
	{	cll->cut_from_my_bdry ( this, tag::do_not_bother );   }
	// tag::do_not_bother is useful for a Mesh::Connected::OneDim
	
	// the four methods below are only relevant for vertices
  virtual void add_to_seg ( Cell::PositiveSegment * seg ) = 0;
  virtual void add_to_seg ( Cell::PositiveSegment * seg, const tag::DoNotBother & ) = 0;
	virtual void remove_from_seg ( Cell::PositiveSegment * seg ) = 0;
	virtual void remove_from_seg ( Cell::PositiveSegment * seg, const tag::DoNotBother & ) = 0;

	// the eight methods below are not relevant for vertices
  virtual void add_to_mesh ( Mesh::Core * msh ) = 0;
  virtual void add_to_mesh ( Mesh::Core * msh, const tag::DoNotBother & ) = 0;
	virtual void remove_from_mesh ( Mesh::Core * msh ) = 0;
	virtual void remove_from_mesh ( Mesh::Core * msh, const tag::DoNotBother & ) = 0;
  virtual void add_to_bdry ( Mesh::Core * msh ) = 0;
  virtual void add_to_bdry ( Mesh::Core * msh, const tag::DoNotBother & ) = 0;
	virtual void remove_from_bdry ( Mesh::Core * msh ) = 0;
	virtual void remove_from_bdry ( Mesh::Core * msh, const tag::DoNotBother & ) = 0;
	// tag::do_not_bother is useful for a Mesh::Connected::OneDim

	virtual void glue_on_my_bdry ( Cell::Core * ) = 0;
	virtual void glue_on_my_bdry ( Cell::Core *, const tag::DoNotBother & ) = 0;
	virtual void cut_from_my_bdry ( Cell::Core * ) = 0;
	virtual void cut_from_my_bdry ( Cell::Core *, const tag::DoNotBother & ) = 0;
	// tag::do_not_bother is useful for a Mesh::Connected::OneDim
	// we feel that 'glue_on_bdry_of' and 'cut_from_bdry_of' are more readable
	// so we suggest to use those (see above)

	virtual void compute_sign ( short int & cp, short int & cn, Mesh::Core * const cell_bdry ) = 0;
		
	#ifndef NDEBUG
	virtual std::string get_name ( ) = 0;
	virtual void print_everything ( ) = 0;
	#endif

}; // end of class Cell::Core

//-----------------------------------------------------------------------------//


class Cell::Positive : public Cell::Core

// abstract class, introduces attribute  meshes
// and defines virtual methods  is_positive, get_positive

// specialized in Cell::Positive{Vertex,Segment,HighDim}

// I would very much prefer the name
// class 'Cell::Positive::Vertex' instead of 'Cell::PositiveVertex'
// but it was not possible ...

{	public :

	typedef Cell::PositiveVertex Vertex;
	typedef Cell::PositiveNotVertex NotVertex;
	typedef Cell::PositiveSegment Segment;
	typedef Cell::PositiveHighDim HighDim;
	
	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// double_heap, size_t_heap, short_int_heap, hook
	// and cell_behind_within inherited from Cell::Core

	// Cell reverse_attr  inherited from Cell::Core

	// the 'meshes' attribute keeps information about all meshes
	// "above" 'this' cell, that is, all meshes containing 'this' as a cell
	// it is indexed over the dimension of the mesh minus the dimension of 'this' 
	// for each mesh, it keeps a 'Cell::field_to_meshes' value, containing two counters
	// and an iterator into the 'cells' field of that particular mesh
	// the iterator is only meaningful for fuzzy or STSI meshes
	// of course this implies quite some amount of redundant information
	// but this redundancy makes the classes fast, especially for remeshing
	
	// since indices of vectors begin at zero
	// the keys of 'meshes[i]' will be meshes of dimension 'i + this->dim'
	// however, for meshes of the same dimension as 'this' cell,
	// we keep a different record in Cell::Positive::Vertex::segments
	// and Cell::Positive::NotVertex::meshes_same_dim
	// so meshes[0] will always be an empty map

	std::vector < std::map < Mesh::Core*, Cell::field_to_meshes > > meshes;

	static std::vector < size_t > double_heap_size, size_t_heap_size, short_int_heap_size;

	inline Positive ( const tag::OfDimension &, const size_t d,
	                  const tag::SizeMeshes &, const size_t sz, const tag::OneDummyWrapper & )
	:	Cell::Core ( tag::of_dim, d, tag::has_no_reverse, tag::one_dummy_wrapper ),
		meshes ( sz )
	{	}

	// bool dispose_query ( )  defined by tag::Util::Core::DelegateDispose ifdef COLLECT_CM
	
	bool is_positive ( ) const;  // virtual from Cell::Core
	Cell get_positive ( );  // virtual from Cell::Core

	virtual Cell::Negative * build_reverse ( const tag::OneDummyWrapper & ) = 0;

	bool belongs_to  // cell has dimension lower than the mesh, virtual from Cell::Core
	( Mesh::Core *, const tag::CellHasLowDim &,
	  const tag::NotOriented & no = tag::not_oriented ) const;
	// belongs_to (with tag::same_dim, tag::oriented )  stays pure virtual from Cell::Core
	// belongs_to (with tag::same_dim, tag::not_oriented )  defined by Cell::Core

	// two versions of glue_on_bdry_of and two of cut_from_bdry_of defined inline by Cell::Core

	// add_to_seg  and  remove_from_seg  stay pure virtual from Cell:Core

	// (two versions of each) add_to_mesh, remove_from_mesh, add_to_bdry, remove_from_bdry
	// stay pure virtual from Cell:Core

	// two versions of glue_on_my_bdry and two of cut_from_my_bdry
	// stay pure virtual from Cell:Core
	
	// compute_sign  stays pure virtual from Cell::Core
	
	#ifndef NDEBUG
	std::string get_name ( );  // virtual from Cell::Core
	// void print_everything ( );  stays pure virtual from Cell::Core
	#endif
	
}; // end of class Cell::Positive

//-----------------------------------------------------------------------------//


class Cell::Negative : public Cell::Core

// abstract class, defines virtual methods
// is_positive, get_positive, belongs_to,
// glue_on_my_bdry, cut_from_my_bdry

// specialized in Cell::Negative::{Vertex,Segment,HighDim}

{	public :

	typedef Cell::NegativeVertex Vertex;
	typedef Cell::NegativeNotVertex NotVertex;
	typedef Cell::NegativeSegment Segment;
	typedef Cell::NegativeHighDim HighDim;
	// I would very much prefer the name 'Cell::Negative::Vertex', 'Cell::Negative::Segment'
	// instead of 'Cell::NegativeVertex', 'Cell::NegativeSegment'
	// but it was not possible ...
	
	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// double_heap, size_t_heap, short_int_heap, hook
	// and cell_behind_within inherited from Cell::Core

	// Cell reverse_attr  inherited from Cell::Core

	static std::vector < size_t > double_heap_size, size_t_heap_size, short_int_heap_size;

	inline Negative ( const tag::OfDimension &, const size_t d,
	                  const tag::ReverseOf &, Cell::Positive * direct, const tag::OneDummyWrapper & )
	: Cell::Core ( tag::of_dim, d, tag::reverse_of, direct, tag::one_dummy_wrapper )
	{	}

	bool is_positive ( ) const;  // virtual from Cell::Core
	Cell get_positive ( );  // virtual from Cell::Core
	// virtual size_t get_dim ( ) const = 0;  // declared in Cell::Core

	bool belongs_to  // cell has dimension lower than the mesh, virtual from Cell::Core
	( Mesh::Core *, const tag::CellHasLowDim &,
	  const tag::NotOriented & no = tag::not_oriented ) const;
	// belongs_to (with tag::same_dim, tag::oriented )  stays pure virtual from Cell::Core
	// belongs_to (with tag::same_dim, tag::not_oriented )  defined by Cell::Core

	// two versions of glue_on_bdry_of and two of cut_from_bdry_of defined inline by Cell::Core

	// add_to_seg  and  remove_from_seg  stay pure virtual from Cell:Core

	// (two versions of each) add_to_mesh, remove_from_mesh, add_to_bdry, remove_from_bdry
	// stay pure virtual from Cell:Core

	// four methods below are virtual from Cell::Core
	void glue_on_my_bdry ( Cell::Core * );
	void glue_on_my_bdry ( Cell::Core *, const tag::DoNotBother & );
	void cut_from_my_bdry ( Cell::Core * );
	void cut_from_my_bdry ( Cell::Core *, const tag::DoNotBother & );
	// tag::do_not_bother is useful for a Mesh::Connected::OneDim
	// we feel that 'glue_on_bdry_of' and 'cut_from_bdry_of' are more readable
	// so we suggest to use those

	// method below is virtual from Cell::Core, here execution forbidden
	void compute_sign ( short int & cp, short int & cn, Mesh::Core * const cell_bdry );

#ifndef NDEBUG
	std::string get_name ( );  // virtual from Cell::Core
	// void print_everything ( );  stays pure virtual from Cell::Core
	#endif
	
}; // end of class Cell::Negative

//-----------------------------------------------------------------------------//


// I would very much prefer the name
// 'class Cell::Positive::Vertex' instead of 'Cell::PositiveVertex'
// but it was not possible ...

class Cell::PositiveVertex : public Cell::Positive

{	public :

	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// double_heap, size_t_heap, short_int_heap, hook
	// and cell_behind_within inherited from Cell::Core

	// Cell reverse_attr  inherited from Cell::Core

	// meshes  inherited from Cell::Positive

	// in 'segments' we keep record of meshes of dimension zero "above" 'this' vertex,
	// that is, of segments having 'this' extremity
	// the 'short int' is a sign, 1 or -1

	std::map < Cell::Positive::Segment*, short int > segments;
	
	inline PositiveVertex ( const tag::OneDummyWrapper & );

	PositiveVertex ( const Cell::Positive::Vertex & ) = delete;
	PositiveVertex ( const Cell::Positive::Vertex && ) = delete;
	Cell::Positive::Vertex & operator= ( const Cell::Positive::Vertex & ) = delete;
	Cell::Positive::Vertex & operator= ( const Cell::Positive::Vertex && ) = delete;

	// bool dispose_query ( )  defined by tag::Util::Core::DelegateDispose ifdef COLLECT_CM
	
	// is_positive  and  get_positive  defined by Cell::Positive
	size_t get_dim ( ) const; // virtual from Cell::Core

	// tip  and  base  defined by Cell::Core, execution forbidden
	
	Mesh boundary ( );  // virtual from Cell::Core, here execution forbidden

	Cell::Negative * build_reverse ( const tag::OneDummyWrapper & );
	// virtual from Cell::Positive

	// belongs_to (with tag::low_dim)  defined by Cell::Positive
	bool belongs_to  // virtual from Cell::Core
	( Mesh::Core *, const tag::SameDim &, const tag::Oriented & o = tag::oriented ) const;
	// belongs_to (with tag::same_dim, tag::not_oriented )  defined by Cell::Core

	// glue_on_bdry_of  and  cut_from_bdry_of  defined by Cell::Core

	// four methods below virtual from Cell::Core
	void add_to_seg ( Cell::Positive::Segment * seg );
	void add_to_seg ( Cell::Positive::Segment * seg, const tag::DoNotBother & );
	void remove_from_seg ( Cell::Positive::Segment * seg );
	void remove_from_seg ( Cell::Positive::Segment * seg, const tag::DoNotBother & );

	// the twelve methods below are virtual from Cell::Core, here execution forbidden
	void add_to_mesh ( Mesh::Core * msh );
	void add_to_mesh ( Mesh::Core * msh, const tag::DoNotBother & );
	void remove_from_mesh ( Mesh::Core * msh );
	void remove_from_mesh ( Mesh::Core * msh, const tag::DoNotBother & );
  void add_to_bdry ( Mesh::Core * msh );
  void add_to_bdry ( Mesh::Core * msh, const tag::DoNotBother & );
	void remove_from_bdry ( Mesh::Core * msh );
	void remove_from_bdry ( Mesh::Core * msh, const tag::DoNotBother & );
	void glue_on_my_bdry ( Cell::Core * );
	void glue_on_my_bdry ( Cell::Core *, const tag::DoNotBother & );
	void cut_from_my_bdry ( Cell::Core * );
	void cut_from_my_bdry ( Cell::Core *, const tag::DoNotBother & );

	// method below is virtual from Cell::Core
	void compute_sign ( short int & cp, short int & cn, Mesh::Core * const cell_bdry );
	
	#ifndef NDEBUG
	// std::string get_name ( )  defined by Cell::Positive
	void print_everything ( );  // virtual from Cell::Core
	#endif
	
}; // end of  class Cell::PositiveVertex, aka Cell::Positive::Vertex

//-----------------------------------------------------------------------------//


// I would very much prefer the name
// 'class Cell::Negative::Vertex' instead of 'Cell::NegativeVertex'
// but it was not possible ...

class Cell::NegativeVertex : public Cell::Negative

{	public :

	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// double_heap, size_t_heap, short_int_heap, hook
	// and cell_behind_within inherited from Cell::Core

	// Cell reverse_attr  inherited from Cell::Core

	inline NegativeVertex
	( const tag::ReverseOf &, Cell::Positive * direct_ver_p, const tag::OneDummyWrapper & );

	NegativeVertex ( const Cell::Negative::Vertex & ) = delete;
	NegativeVertex ( const Cell::Negative::Vertex && ) = delete;
	Cell::Negative::Vertex & operator= ( const Cell::Negative::Vertex & ) = delete;
	Cell::Negative::Vertex & operator= ( const Cell::Negative::Vertex && ) = delete;

	// bool dispose_query ( )  defined by tag::Util::Core::DelegateDispose ifdef COLLECT_CM
	
	// is_positive  and  get_positive  defined by Cell::Negative
	size_t get_dim ( ) const; // virtual from Cell::Core
	
	// tip  and  base  defined by Cell::Core, execution forbidden

	Mesh boundary ( );  // virtual from Cell::Core, here execution forbidden

	// belongs_to (with tag::low_dim)  defined by Cell::Negative
	bool belongs_to   // virtual from Cell::Core
	( Mesh::Core *, const tag::SameDim &, const tag::Oriented & o = tag::oriented ) const;
	// belongs_to (with tag::same_dim, tag::not_oriented )  defined by Cell::Core

	// glue_on_bdry_of  and  cut_from_bdry_of  defined by Cell::Core

	// methods 'add_to' and 'remove_from' add/remove 'this' cell to/from the mesh 'msh'
	// see paragraph 11.9 in the manual

	// four methods below virtual from Cell::Core
	void add_to_seg ( Cell::Positive::Segment * seg );
	void add_to_seg ( Cell::Positive::Segment * seg, const tag::DoNotBother & );
	void remove_from_seg ( Cell::Positive::Segment * seg );
	void remove_from_seg ( Cell::Positive::Segment * seg, const tag::DoNotBother & );

	// the eight methods below are virtual from Cell::Core, here execution forbidden
	void add_to_mesh ( Mesh::Core * msh );
	void add_to_mesh ( Mesh::Core * msh, const tag::DoNotBother & );
	void remove_from_mesh ( Mesh::Core * msh );
	void remove_from_mesh ( Mesh::Core * msh, const tag::DoNotBother & );
  void add_to_bdry ( Mesh::Core * msh );
  void add_to_bdry ( Mesh::Core * msh, const tag::DoNotBother & );
	void remove_from_bdry ( Mesh::Core * msh );
	void remove_from_bdry ( Mesh::Core * msh, const tag::DoNotBother & );

	// glue_on_my_bdry  and  cut_from_my_bdry  defined by Cell:Negative
	// compute_sign  defined by Cell::Negative, execution forbidden
	
	#ifndef NDEBUG
	// std::string get_name ( )  defined by Cell::Negative
	void print_everything ( );  // virtual from Cell::Core
	#endif
	
}; // end of  class Cell::NegativeVertex, aka Cell::Negative::Vertex

//-----------------------------------------------------------------------------//


// I would very much prefer the name
// 'class Cell::Positive::NotVertex' instead of 'Cell::PositiveNotVertex'
// but it was not possible ...

class Cell::PositiveNotVertex : public Cell::Positive

// abstract class, useful only for introducing the attribute  meshes_same_dim
// specialized in Cell::Positive::{Segment,HighDim}

{	public :

	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// double_heap, size_t_heap, short_int_heap, hook
	// and cell_behind_within inherited from Cell::Core

	// Cell reverse_attr  inherited from Cell::Core

	// meshes  inherited from Cell::Positive

	// in 'meshes_same_dim' we keep record of meshes "above" 'this' cell,
	//   of the same dimension
	// in Cell::field_to_meshes_same_dim, the 'short int sign' is a sign, 1 or -1
	// the iterator 'where' is only meaningful for fuzzy or STSI meshes

	std::map < Mesh::Core*, Cell::field_to_meshes_same_dim > meshes_same_dim;
	
	inline PositiveNotVertex ( const tag::OfDimension &, const size_t d,
                             const tag::SizeMeshes &, const size_t sz,
                             const tag::OneDummyWrapper &              )
	: Cell::Positive ( tag::of_dim, d, tag::size_meshes, sz, tag::one_dummy_wrapper )
	{	}
	
	// bool dispose_query ( )  defined by tag::Util::Core::DelegateDispose ifdef COLLECT_CM
	
	// Cell::Negative * build_reverse ( const tag::OneDummyWrapper & )
	//   stays pure virtual from Cell::Positive

	// belongs_to (with low_dim)  defined by Cell::Positive
	bool belongs_to  // virtual from Cell::Core
	( Mesh::Core *, const tag::SameDim &, const tag::Oriented & o = tag::oriented ) const;
	// belongs_to (with tag::same_dim, tag::not_oriented )  defined by Cell::Core

	// four methods below virtual from Cell::Core, here execution forbidden
	void add_to_seg ( Cell::Positive::Segment * seg );
	void add_to_seg ( Cell::Positive::Segment * seg, const tag::DoNotBother & );
	void remove_from_seg ( Cell::Positive::Segment * seg );
	void remove_from_seg ( Cell::Positive::Segment * seg, const tag::DoNotBother & );

	// add_to_mesh, remove_from_mesh, add_to_bdry, remove_from_bdry
	//   stay pure virtual from Cell::Core

	// method below is virtual from Cell::Core
	void compute_sign ( short int & cp, short int & cn, Mesh::Core * const cell_bdry );
		
	#ifndef NDEBUG
	// std::string get_name ( )  defined by Cell::Positive
	// void print_everything ( )  stays pure virtual from Cell::Core
	#endif

}; // end of  class Cell::PositiveNotVertex, aka Cell::Positive::NotVertex

//-----------------------------------------------------------------------------//


// I would very much prefer the name
// 'class Cell::Negative::NotVertex' instead of 'Cell::NegativeNotVertex'
// but it was not possible ...

class Cell::NegativeNotVertex : public Cell::Negative

// abstract class, useful only for introducing a method belongs_to
// specialized in Cell::Negative::{Segment,HighDim}

{	public :

	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// double_heap, size_t_heap, short_int_heap, hook
	// and cell_behind_within inherited from Cell::Core

	// Cell reverse_attr  inherited from Cell::Core

	inline NegativeNotVertex
	( const tag::OfDimension &, size_t d,
	  const tag::ReverseOf &, Cell::Positive * direct_cll_p, const tag::OneDummyWrapper & );

	// belongs_to (with tag::low_dim)  defined by Cell::Negative
	bool belongs_to   // virtual from Cell::Core
	( Mesh::Core *, const tag::SameDim &, const tag::Oriented & o = tag::oriented ) const;
	// belongs_to (with tag::same_dim, tag::not_oriented )  defined by Cell::Core

	// four methods below virtual from Cell::Core, here execution forbidden
	void add_to_seg ( Cell::Positive::Segment * seg );
	void add_to_seg ( Cell::Positive::Segment * seg, const tag::DoNotBother & );
	void remove_from_seg ( Cell::Positive::Segment * seg );
	void remove_from_seg ( Cell::Positive::Segment * seg, const tag::DoNotBother & );

	// add_to_mesh, remove_from_mesh, add_to_bdry, remove_from_bdry stay pure virtual from Cell::Core
	// compute_sign  defined by Cell::Negative, execution forbidden
	
	#ifndef NDEBUG
	// std::string get_name ( )  defined by Cell::Positive
	// void print_everything ( )  stays pure virtual from Cell::Core
	#endif

};  // end of  class Cell::NegativeNotVertex (aka class Cell::Negative::NotVertex)
	
//-----------------------------------------------------------------------------//


	
// I would very much prefer the name
// 'class Cell::Positive::Segment' instead of 'Cell::PositiveSegment'
// but it was not possible ...

class Cell::PositiveSegment : public Cell::Positive::NotVertex

{	public :

	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// double_heap, size_t_heap, short_int_heap, hook
	// and cell_behind_within inherited from Cell::Core

	// Cell reverse_attr  inherited from Cell::Core

	// meshes  inherited from Cell::Positive
	// meshes_same_dim  inherited from Cell::Positive::NotVertex

	// here we use Cell wrappers as pointers
	// because a segment should keep its extremities alive
	Cell base_attr, tip_attr;

	inline PositiveSegment ( Cell A, Cell B, const tag::OneDummyWrapper & );

	PositiveSegment ( const Cell::Positive::Segment & ) = delete;
	PositiveSegment ( const Cell::Positive::Segment && ) = delete;
	Cell::Positive::Segment & operator= ( const Cell::Positive::Segment & ) = delete;
	Cell::Positive::Segment & operator= ( const Cell::Positive::Segment && ) = delete;

	// is_positive  and  get_positive  defined by Cell::Positive
	size_t get_dim ( ) const; // virtual from Cell::Core

	Cell tip () override; // virtual, overrides definition by Cell::Core
	Cell base () override; // virtual, overrides definition by Cell::Core

	Mesh boundary ( );  // virtual from Cell::Core

	Cell::Negative * build_reverse ( const tag::OneDummyWrapper & );
  // virtual from Cell::Positive

	// belongs_to (with tag::low_dim)  defined by Cell::Positive
	// belongs_to (with tag::same_dim, tag::oriented )  defined by Cell::Positive::NotVertex
	// belongs_to (with tag::same_dim, tag::not_oriented )  defined by Cell::Core

	// glue_on_bdry_of  and  cut_from_bdry_of  defined by Cell::Core

	// methods 'add_to' and 'remove_from' add/remove 'this' cell to/from the mesh 'msh'
	// see paragraph 11.9 in the manual

	// add_to_seg and remove_from seg virtual from Cell::Core,
	// defined by Cell::Positive::NotVertex, execution forbidden

	// the twelve methods below are virtual from Cell::Core
	void add_to_mesh ( Mesh::Core * msh );
	void add_to_mesh ( Mesh::Core * msh, const tag::DoNotBother & );
	void remove_from_mesh ( Mesh::Core * msh );
	void remove_from_mesh ( Mesh::Core * msh, const tag::DoNotBother & );
  void add_to_bdry ( Mesh::Core * msh );
  void add_to_bdry ( Mesh::Core * msh, const tag::DoNotBother & );
	void remove_from_bdry ( Mesh::Core * msh );
	void remove_from_bdry ( Mesh::Core * msh, const tag::DoNotBother & );
	void glue_on_my_bdry ( Cell::Core * );
	void glue_on_my_bdry ( Cell::Core *, const tag::DoNotBother & );
	void cut_from_my_bdry ( Cell::Core * );
	void cut_from_my_bdry ( Cell::Core *, const tag::DoNotBother & );
	// tag::do_not_bother is useful for a Mesh::Connected::OneDim	
	// we feel that 'glue_on_bdry_of' and 'cut_from_bdry_of' are more readable
	// so we suggest to use those (see class Cell::Core)

	// compute_sign  defined by Cell::Positive::NotVertex

	#ifndef NDEBUG
	// std::string get_name ( )  defined by Cell::Positive
	void print_everything ( );  // virtual from Cell::Core
	#endif
	
}; // end of  class Cell::PositiveSegment, aka Cell::Positive::Segment

//-----------------------------------------------------------------------------//


// I would very much prefer the name
// 'class Cell::Negative::Segment' instead of 'Cell::NegativeSegment'
// but it was not possible ...

class Cell::NegativeSegment : public Cell::Negative::NotVertex

{	public :

	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// double_heap, size_t_heap, short_int_heap, hook
	// and cell_behind_within inherited from Cell::Core

	// Cell reverse_attr  inherited from Cell::Core

	inline NegativeSegment
		( const tag::ReverseOf &, Cell::Positive * direct_seg_p, const tag::OneDummyWrapper & );

	// bool dispose_query ( )  defined by tag::Util::Core::DelegateDispose ifdef COLLECT_CM
	
	NegativeSegment ( const Cell::Negative::Segment & ) = delete;
	NegativeSegment ( const Cell::Negative::Segment && ) = delete;
	Cell::Negative::Segment & operator= ( const Cell::Negative::Segment & ) = delete;
	Cell::Negative::Segment & operator= ( const Cell::Negative::Segment && ) = delete;

	// is_positive  and  get_positive  defined by Cell::Negative
	size_t get_dim ( ) const; // virtual from Cell::Core
	
	Cell tip () override;  // virtual, overrides definition by Cell::Core
	Cell base () override;  // virtual, overrides definition by Cell::Core

	Mesh boundary ( );  // virtual from Cell::Core

	// belongs_to (with tag::low_dim)  defined by Cell::Negative
	// belongs_to (with tag::same_dim, tag::oriented)  defined by Cell::Negative::NotVertex
	// belongs_to (with tag::same_dim, tag::not_oriented )  defined by Cell::Core

	// glue_on_bdry_of  and  cut_from_bdry_of  defined by Cell::Core

	// methods 'add_to' and 'remove_from' add/remove 'this' cell to/from the mesh 'msh'
	// see paragraph 11.9 in the manual

	// add_to_seg and remove_from seg virtual from Cell::Core,
	// defined by Cell::Negative::NotVertex, execution forbidden

	// the eight methods below are virtual from Cell::Core
	void add_to_mesh ( Mesh::Core * msh );
	void add_to_mesh ( Mesh::Core * msh, const tag::DoNotBother & );
	void remove_from_mesh ( Mesh::Core * msh );
	void remove_from_mesh ( Mesh::Core * msh, const tag::DoNotBother & );
  void add_to_bdry ( Mesh::Core * msh );
  void add_to_bdry ( Mesh::Core * msh, const tag::DoNotBother & );
	void remove_from_bdry ( Mesh::Core * msh );
	void remove_from_bdry ( Mesh::Core * msh, const tag::DoNotBother & );
	// tag::do_not_bother is useful for a Mesh::Connected::OneDim	

	// glue_on_my_bdry and cut_from_my_bdry  defined by Cell:Negative
	// compute_sign  defined by Cell::Negative, execution forbidden

	#ifndef NDEBUG
	// std::string get_name ( )  defined by Cell::Negative
	void print_everything ( );  // virtual from Cell::Core
	#endif
	
}; // end of  class Cell::NegativeSegment, aka Cell::Negative::Segment

//-----------------------------------------------------------------------------//


// I would very much prefer the name
// 'class Cell::Positive::HighDim' instead of 'Cell::PositiveHighDim'
// but it was not possible ...

class Cell::PositiveHighDim : public Cell::Positive::NotVertex

// a cell of dimension >= 2

{	public :

	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// double_heap, size_t_heap, short_int_heap, hook
	// and cell_behind_within inherited from Cell::Core

	// Cell reverse_attr  inherited from Cell::Core

	// meshes  inherited from Cell::Positive
	// meshes_same_dim  inherited from Cell::Positive::NotVertex

	// we use Mesh wrapper as pointer
	// because a Cell should keep its boundary alive
	Mesh boundary_attr;  // positive mesh

	inline PositiveHighDim ( const tag::OfDimension &, const size_t d,
	                         const tag::WhoseBoundaryIs &, const Mesh & msh,
	                         const tag::OneDummyWrapper &                   );
	inline PositiveHighDim ( const tag::WhoseBoundaryIs &, const Mesh & msh,
	                         const tag::OneDummyWrapper &                   );
	inline PositiveHighDim ( const tag::Triangle &,
	                         const Cell & AB, const Cell & BC, const Cell & CA,
	                         const tag::OneDummyWrapper &                       );
	inline PositiveHighDim ( const tag::Quadrangle &, const Cell & AB, const Cell & BC,
	                         const Cell & CD, const Cell & DA, const tag::OneDummyWrapper & );
	inline PositiveHighDim ( const tag::Pentagon &, const Cell & AB, const Cell & BC,
	                         const Cell & CD, const Cell & DE, const Cell & EA,
                           const tag::OneDummyWrapper &                             );
	inline PositiveHighDim ( const tag::Hexagon &, const Cell & AB, const Cell & BC,
	                         const Cell & CD, const Cell & DE, const Cell & EF, const Cell & FA,
                           const tag::OneDummyWrapper &                                        );

	PositiveHighDim ( const Cell::Positive::HighDim & ) = delete;
	PositiveHighDim ( const Cell::Positive::HighDim && ) = delete;
	Cell::Positive::HighDim & operator= ( const Cell::Positive::HighDim & ) = delete;
	Cell::Positive::HighDim & operator= ( const Cell::Positive::HighDim && ) = delete;

	// bool dispose_query ( )  defined by tag::Util::Core::DelegateDispose ifdef COLLECT_CM
	
	// is_positive  and  get_positive  defined by Cell::Positive
	size_t get_dim ( ) const; // virtual

	// tip  and  base  defined by Cell::Core, execution forbidden

	Mesh boundary ( );  // virtual from Cell::Core

	Cell::Negative * build_reverse ( const tag::OneDummyWrapper & );
  // virtual from Cell::Positive

	// belongs_to (with tag::low_dim)  defined by Cell::Positive
	// belongs_to (with tag::same_dim, tag::oriented )  defined by Cell::Positive::NotVertex
	// belongs_to (with tag::same_dim, tag::not_oriented )  defined by Cell::Core

	// glue_on_bdry_of  and  cut_from_bdry_of  defined by Cell::Core

	// methods 'add_to' and 'remove_from' add/remove 'this' cell to/from the mesh 'msh'
	// see paragraph 11.9 in the manual

	// add_to_seg and remove_from seg virtual from Cell::Core,
	// defined by Cell::Positive::NotVertex, execution forbidden

	// the twelve methods below are virtual from Cell::Core
	void add_to_mesh ( Mesh::Core * msh );
	void add_to_mesh ( Mesh::Core * msh, const tag::DoNotBother & );
	void remove_from_mesh ( Mesh::Core * msh );
	void remove_from_mesh ( Mesh::Core * msh, const tag::DoNotBother & );
  void add_to_bdry ( Mesh::Core * msh );
  void add_to_bdry ( Mesh::Core * msh, const tag::DoNotBother & );
	void remove_from_bdry ( Mesh::Core * msh );
	void remove_from_bdry ( Mesh::Core * msh, const tag::DoNotBother & );
	void glue_on_my_bdry ( Cell::Core * );
	void glue_on_my_bdry ( Cell::Core *, const tag::DoNotBother & );
	void cut_from_my_bdry ( Cell::Core * );
	void cut_from_my_bdry ( Cell::Core *, const tag::DoNotBother & );
	// tag::do_not_bother is useful for a Mesh::Connected::OneDim	
	// we feel that 'glue_on_bdry_of' and 'cut_from_bdry_of' are more readable
	// so we suggest to use those (see class Cell::Core)

	// compute_sign  defined by Cell::Positive::NotVertex
	
	#ifndef NDEBUG
	// std::string get_name ( )  defined by Cell::Positive
	void print_everything ( );  // virtual from Cell::Core
	#endif

}; // end of  class Cell::PositiveHighDim, aka Cell::Positive::HighDim

//-----------------------------------------------------------------------------//


// I would very much prefer the name
// 'class Cell::Negative::HighDim' instead of 'Cell::NegativeHighDim'
// but it was not possible ...

class Cell::NegativeHighDim : public Cell::Negative::NotVertex

// a negative cell of dimension >= 2

{	public :

	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// double_heap, size_t_heap, short_int_heap, hook
	// and cell_behind_within inherited from Cell::Core

	// Cell reverse_attr  inherited from Cell::Core

	inline NegativeHighDim ( const tag::OfDimension &, const size_t d,
	                         const tag::ReverseOf &, Cell::Positive * direct_cell_p,
                           const tag::OneDummyWrapper &                            );
	
	inline NegativeHighDim
		( const tag::ReverseOf &, Cell::Positive * direct_cell_p, const tag::OneDummyWrapper & );

	NegativeHighDim ( const Cell::Negative::HighDim & ) = delete;
	NegativeHighDim ( const Cell::Negative::HighDim && ) = delete;
	Cell::Negative::HighDim & operator= ( const Cell::Negative::HighDim & ) = delete;
	Cell::Negative::HighDim & operator= ( const Cell::Negative::HighDim && ) = delete;

	// bool dispose_query ( )  defined by tag::Util::Core::DelegateDispose ifdef COLLECT_CM
	
	// is_positive  and  get_positive  defined by Cell::Negative
	size_t get_dim ( ) const; // virtual from Cell::Core

	// tip  and  base  defined by Cell::Core, execution forbidden
	
	Mesh boundary ( );  // virtual from Cell::Core

	// belongs_to (with tag::low_dim)  defined by Cell::Negative
	// belongs_to (with tag::same_dim and tag::oriented)  defined by Cell::Negative::NotVertex
	// belongs_to (with tag::same_dim, tag::not_oriented )  defined by Cell::Core

	// glue_on_bdry_of  and  cut_from_bdry_of  defined by Cell::Core

	// methods 'add_to' and 'remove_from' add/remove 'this' cell to/from the mesh 'msh'
	// see paragraph 11.9 in the manual

	// add_to_seg and remove_from seg virtual from Cell::Core,
	// defined by Cell::Negative::NotVertex, execution forbidden

	// the eight methods below are virtual from Cell::Core
	void add_to_mesh ( Mesh::Core * msh );
	void add_to_mesh ( Mesh::Core * msh, const tag::DoNotBother & );
	void remove_from_mesh ( Mesh::Core * msh );
	void remove_from_mesh ( Mesh::Core * msh, const tag::DoNotBother & );
  void add_to_bdry ( Mesh::Core * msh );
  void add_to_bdry ( Mesh::Core * msh, const tag::DoNotBother & );
	void remove_from_bdry ( Mesh::Core * msh );
	void remove_from_bdry ( Mesh::Core * msh, const tag::DoNotBother & );
	// tag::do_not_bother is useful for a Mesh::Connected::OneDim	

	// glue_on_my_bdry  and  cut_from_my_bdry  defined by Cell:Negative
	// compute_sign  defined by Cell::Negative, execution forbidden

	#ifndef NDEBUG
	// std::string get_name ( )  defined by Cell::Negative
	void print_everything ( );  // virtual from Cell::Core
	#endif

}; // end of  class Cell::NegativeHighDim, aka Cell::Negative::HighDim


//-----------------------------------------------------------------------------//
//------------------       core Meshes     ------------------------------------//
//-----------------------------------------------------------------------------//


// I would very much prefer the name 'class Mesh::Core' instead of 'tag::Util::MeshCore'
// but it was not possible ...

#ifdef MANIFEM_COLLECT_CM	

class tag::Util::MeshCore : public tag::Util::Core

#else  // no MANIFEM_COLLECT_CM

class tag::Util::MeshCore

#endif  // MANIFEM_COLLECT_CM	

// represents a positive mesh
// negative meshes have no core (wrappers for negative meshes are built on-the-fly)

// abstract class, specialized in Mesh::ZeroDim, Mesh::Connected::{OneDim,HighDim},
// Mesh::MultiplyConnected::{OneDim,HighDim}, Mesh::Fuzzy, Mesh::STSI

{	public :

	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// we use an ordinary pointer, not a wrapper
	// because we do not want a cell to be kept alive by its boundary
	Cell::Positive * cell_enclosed { nullptr };
	// when 'this' is the boundary of a cell

	inline MeshCore ( const tag::OfDimension &, const size_t d , const tag::MinusOne &,
                    const tag::OneDummyWrapper &                                      )
	#ifdef MANIFEM_COLLECT_CM	
	:	tag::Util::Core ( tag::one_dummy_wrapper )
	#endif  // MANIFEM_COLLECT_CM	
	{	}

	inline MeshCore ( const tag::OfDimension &, const size_t d, const tag::MinusOne &,
	                  const tag::BoundaryOf &, const tag::PositiveCell &,
	                  Cell::Positive * cll, const tag::OneDummyWrapper &               )
	#ifdef MANIFEM_COLLECT_CM	
	:	tag::Util::Core ( tag::one_dummy_wrapper ),
		cell_enclosed { cll }
	#else  // no MANIFEM_COLLECT_CM
	:	cell_enclosed { cll }
	#endif  // MANIFEM_COLLECT_CM
	{ }

	// bool dispose_query ( )  defined by tag::Util::Core ifdef COLLECT_CM
	
	virtual size_t get_dim_plus_one ( ) = 0;

	virtual size_t number_of ( const tag::Vertices & ) = 0;
	virtual size_t number_of ( const tag::Segments & ) = 0;
	virtual size_t number_of ( const tag::CellsOfDim &, const size_t d ) = 0;
	virtual size_t number_of ( const tag::CellsOfMaxDim & ) = 0;

	// the four methods below are only relevant for connected one-dimensional meshes
	// so we forbid execution for now and then override them in Mesh::Connected::OneDim
	virtual Cell first_vertex ( );
	virtual Cell last_vertex ( );
	virtual Cell first_segment ( );
	virtual Cell last_segment ( );

	// the four methods below are only relevant for STSI meshes
	// so we forbid execution for now and then override them in Mesh::STSI
	virtual Cell::Core * cell_in_front_of
	( const Cell::Core * face_p, const tag::SeenFrom &, const Cell neighbour,
	  const tag::SurelyExists & se = tag::surely_exists                       ) const;
	virtual Cell::Core * cell_in_front_of
	( const Cell::Core * face_p, const tag::SeenFrom &, const Cell neighbour,
	  const tag::MayNotExist &                                                ) const;
	virtual Cell::Core * cell_behind
	( const Cell::Core * face_p, const tag::SeenFrom &, const Cell neighbour,
	  const tag::SurelyExists & se = tag::surely_exists                       ) const;
	virtual Cell::Core * cell_behind
	( const Cell::Core * face_p, const tag::SeenFrom &, const Cell neighbour,
	  const tag::MayNotExist &                                                ) const;

	virtual void add_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsNotBdry & ) = 0;
	virtual void add_pos_seg
		( Cell::Positive::Segment *, const tag::MeshIsNotBdry &, const tag::DoNotBother & ) = 0;
	virtual void add_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsBdry & ) = 0;
	virtual void add_pos_seg
		( Cell::Positive::Segment *, const tag::MeshIsBdry &, const tag::DoNotBother & ) = 0;
	virtual void remove_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsNotBdry & ) = 0;
	virtual void remove_pos_seg
		( Cell::Positive::Segment *, const tag::MeshIsNotBdry &, const tag::DoNotBother & ) = 0;
	virtual void remove_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsBdry & ) = 0;
	virtual void remove_pos_seg
		( Cell::Positive::Segment *, const tag::MeshIsBdry &, const tag::DoNotBother & ) = 0;
	virtual void add_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsNotBdry & ) = 0;
	virtual void add_neg_seg
		( Cell::Negative::Segment *, const tag::MeshIsNotBdry &, const tag::DoNotBother & ) = 0;
	virtual void add_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsBdry & ) = 0;
	virtual void add_neg_seg
		( Cell::Negative::Segment *, const tag::MeshIsBdry &, const tag::DoNotBother & ) = 0;
	virtual void remove_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsNotBdry & ) = 0;
	virtual void remove_neg_seg
		( Cell::Negative::Segment *, const tag::MeshIsNotBdry &, const tag::DoNotBother & ) = 0;
	virtual void remove_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsBdry & ) = 0;
	virtual void remove_neg_seg
		( Cell::Negative::Segment *, const tag::MeshIsBdry &, const tag::DoNotBother & ) = 0;
	virtual void add_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsNotBdry & ) = 0;
	virtual void add_pos_hd_cell
		( Cell::Positive::HighDim *, const tag::MeshIsNotBdry &, const tag::DoNotBother & ) = 0;
	virtual void add_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsBdry & ) = 0;
	virtual void add_pos_hd_cell
		( Cell::Positive::HighDim *, const tag::MeshIsBdry &, const tag::DoNotBother & ) = 0;
	virtual void remove_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsNotBdry & ) = 0;
	virtual void remove_pos_hd_cell
		( Cell::Positive::HighDim *, const tag::MeshIsNotBdry &, const tag::DoNotBother & ) = 0;
	virtual void remove_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsBdry & ) = 0;
	virtual void remove_pos_hd_cell
		( Cell::Positive::HighDim *, const tag::MeshIsBdry &, const tag::DoNotBother & ) = 0;
	virtual void add_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsNotBdry & ) = 0;
	virtual void add_neg_hd_cell
		( Cell::Negative::HighDim *, const tag::MeshIsNotBdry &, const tag::DoNotBother & ) = 0;
	virtual void add_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsBdry & ) = 0;
	virtual void add_neg_hd_cell
		( Cell::Negative::HighDim *, const tag::MeshIsBdry &, const tag::DoNotBother & ) = 0;
	virtual void remove_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsNotBdry & ) = 0;
	virtual void remove_neg_hd_cell
		( Cell::Negative::HighDim *, const tag::MeshIsNotBdry &, const tag::DoNotBother & ) = 0;
	virtual void remove_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsBdry & ) = 0;
	virtual void remove_neg_hd_cell
		( Cell::Negative::HighDim *, const tag::MeshIsBdry &, const tag::DoNotBother & ) = 0;

	virtual Mesh::Core * build_deep_copy ( ) = 0;
	
	virtual std::list<Cell>::iterator add_to_my_cells ( Cell::Core * const, const size_t ) = 0;
	virtual std::list<Cell>::iterator add_to_my_cells
		( Cell::Core * const, const size_t, const tag::DoNotBother & ) = 0;

	virtual void remove_from_my_cells
		( Cell::Core * const, const size_t, std::list<Cell>::iterator ) = 0;
	virtual void remove_from_my_cells
		( Cell::Core * const, const size_t, std::list<Cell>::iterator, const tag::DoNotBother & ) = 0;
		
	virtual void closed_loop ( const Cell & ver ) = 0;
	virtual void closed_loop ( const Cell & ver, size_t ) = 0;
	// for connected one-dim meshes, set both first_ver and last_ver to 'ver'
	// (and number of segments)
		
	// iterators defined in iterator.cpp
	// we are still in class Mesh::Core

	virtual Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                        ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                        ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &                   ) = 0;

	// we are still in class Mesh::Core

	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                        ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                        ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &                   ) = 0;

	// we are still in class Mesh::Core

	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::AsTheyAre &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &        ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &        ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &   ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ForcePositive &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &            ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &            ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &       ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive &               ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & ) = 0;

	// we are still in class Mesh::Core

	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
	  const tag::ThisMeshIsPositive &                        ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::ThisMeshIsPositive &                            ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::RequireOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                              ) = 0;

	// we are still in class Mesh::Core

	virtual Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &              ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &              ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                 ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                 ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                 ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                  ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                     ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                     ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                     ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                    ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                       ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                       ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                       ) = 0;

	// we are still in class Mesh::Core

	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &        ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &        ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &   ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &    ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &            ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &             ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &       ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ThisMeshIsPositive &                                         ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &              ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &              ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &         ) = 0;

	// we are still in class Mesh::Core

	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                   ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                      ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                      ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                      ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                              ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &      ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ThisMeshIsPositive &                                 ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &      ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &      ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & ) = 0;

	// we are still in class Mesh::Core

	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &              ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                 ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                 ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                 ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                  ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                     ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                     ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                     ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                    ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                       ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                       ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                       ) = 0;

	// we are still in class Mesh::Core

	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &        ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &        ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &   ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &    ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &            ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &             ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &       ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ThisMeshIsPositive &                                         ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &              ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &              ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &         ) = 0;

	// we are still in class Mesh::Core

	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                   ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                      ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                      ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                      ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                              ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &      ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ThisMeshIsPositive &                                 ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &      ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &      ) = 0;
	virtual Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & ) = 0;

	#ifndef NDEBUG
	std::string name;
	std::string get_name();
	virtual void print_everything () = 0;
	#endif

}; // end of  class Mesh::Core

//-----------------------------------------------------------------------------//


class Mesh::ZeroDim : public Mesh::Core

// represents two points, one negative one positive

// zero-dimensional meshes only exist as boundary of a segment,
// which can be retrieved through the cell_enclosed attribute
// they are built on-the-fly

{	public :

	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// Cell::Positive * cell_enclosed inherited from Mesh::Core
	
	inline ZeroDim
		( const tag::BoundaryOf &, const tag::Positive &, const tag::Segment &,
		  Cell::Positive * seg_p, const tag::OneDummyWrapper &                  )
	:	Mesh::Core ( tag::of_dimension, 1, tag::minus_one,
	               tag::boundary_of, tag::positive_cell, seg_p, tag::one_dummy_wrapper )
	{ }

	// bool dispose_query ( )  defined by tag::Util::Core ifdef COLLECT_CM
	
	size_t get_dim_plus_one ( );  // virtual from Mesh::Core

	size_t number_of ( const tag::Vertices & );  // virtual from Mesh::Core
	size_t number_of ( const tag::Segments & );
	// virtual from Mesh::Core, here execution forbidden
	size_t number_of ( const tag::CellsOfDim &, const size_t d );  // virtual from Mesh::Core
	size_t number_of ( const tag::CellsOfMaxDim & );  // virtual from Mesh::Core

	// first_vertex, last_vertex, first_segment and last_segment
	// are defined by Mesh::Core, execution forbidden

	// private:

	// cell_in_front_of  and  cell_behind  defined by Mesh::Core, execution forbidden

	// the thirty-two methods below are virtual from Mesh::Core, here execution forbidden
	void add_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsNotBdry & );
	void add_pos_seg
		( Cell::Positive::Segment *, const tag::MeshIsNotBdry &, const tag::DoNotBother & );
	void add_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsBdry & );
	void add_pos_seg
		( Cell::Positive::Segment *, const tag::MeshIsBdry &, const tag::DoNotBother & );
	void remove_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsNotBdry & );
	void remove_pos_seg
		( Cell::Positive::Segment *, const tag::MeshIsNotBdry &, const tag::DoNotBother & );
	void remove_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsBdry & );
	void remove_pos_seg
		( Cell::Positive::Segment *, const tag::MeshIsBdry &, const tag::DoNotBother & );
	void add_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsNotBdry & );
	void add_neg_seg
		( Cell::Negative::Segment *, const tag::MeshIsNotBdry &, const tag::DoNotBother & );
	void add_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsBdry & );
	void add_neg_seg
		( Cell::Negative::Segment *, const tag::MeshIsBdry &, const tag::DoNotBother & );
	void remove_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsNotBdry & );
	void remove_neg_seg
		( Cell::Negative::Segment *, const tag::MeshIsNotBdry &, const tag::DoNotBother & );
	void remove_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsBdry & );
	void remove_neg_seg
		( Cell::Negative::Segment *, const tag::MeshIsBdry &, const tag::DoNotBother & );
	void add_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsNotBdry & );
	void add_pos_hd_cell
		( Cell::Positive::HighDim *, const tag::MeshIsNotBdry &, const tag::DoNotBother & );
	void add_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsBdry & );
	void add_pos_hd_cell
		( Cell::Positive::HighDim *, const tag::MeshIsBdry &, const tag::DoNotBother & );
	void remove_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsNotBdry & );
	void remove_pos_hd_cell
		( Cell::Positive::HighDim *, const tag::MeshIsNotBdry &, const tag::DoNotBother & );
	void remove_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsBdry & );
	void remove_pos_hd_cell
		( Cell::Positive::HighDim *, const tag::MeshIsBdry &, const tag::DoNotBother & );
	void add_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsNotBdry & );
	void add_neg_hd_cell
		( Cell::Negative::HighDim *, const tag::MeshIsNotBdry &, const tag::DoNotBother & );
	void add_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsBdry & );
	void add_neg_hd_cell
		( Cell::Negative::HighDim *, const tag::MeshIsBdry &, const tag::DoNotBother & );
	void remove_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsNotBdry & );
	void remove_neg_hd_cell
		( Cell::Negative::HighDim *, const tag::MeshIsNotBdry &, const tag::DoNotBother & );
	void remove_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsBdry & );
	void remove_neg_hd_cell
		( Cell::Negative::HighDim *, const tag::MeshIsBdry &, const tag::DoNotBother & );
	
	Mesh::Core * build_deep_copy ( );  // virtual from Mesh::Core, here execution forbidden

	// we are still in class Mesh::ZeroDim
	
	std::list<Cell>::iterator add_to_my_cells ( Cell::Core * const, const size_t );
	std::list<Cell>::iterator add_to_my_cells
		( Cell::Core * const, const size_t, const tag::DoNotBother & );
	// virtual from Cell::Core, here execution forbidden

	void remove_from_my_cells ( Cell::Core * const, const size_t, std::list<Cell>::iterator );
	void remove_from_my_cells
		( Cell::Core * const, const size_t, std::list<Cell>::iterator, const tag::DoNotBother & );
	// virtual from Cell::Core, here execution forbidden
	
	void closed_loop ( const Cell & ver );
	void closed_loop ( const Cell & ver, size_t );
	// virtual from Mesh::Core, here execution forbidden
		
	// iterators are virtual from Mesh::Core and are defined in iterator.cpp

	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	// iterate over the two vertices, first base (negative) then tip (positive)

	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	// iterate over the two vertices, first tip (positive) then base (negative)

	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	// iterate over the two vertices, first base then tip (both positive)

	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	// we iterate over the two vertices, first tip then base (both positive)

	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                        );
	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                        );
	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &                   );

	// we are still in class Mesh::ZeroDim

	Mesh::Iterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                        );
	Mesh::Iterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                        );
	Mesh::Iterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &                   );

	// we are still in class Mesh::ZeroDim
	// twelve iterators below assert dim == 0 then call iterator over vertices

	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::AsTheyAre &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &        );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &        );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &   );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &            );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &            );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &       );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive &               );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );

	// we are still in class Mesh::ZeroDim
	// twelve iterators below simply produce iterator over vertices
	
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );

	// we are still in class Mesh::ZeroDim

	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &              );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &              );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                 );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                 );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                 );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                  );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                     );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                     );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                     );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                    );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                       );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                       );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                       );

	// we are still in class Mesh::ZeroDim

	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &        );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &        );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &   );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &    );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &            );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &             );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &       );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ThisMeshIsPositive &                                         );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &              );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &              );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &         );

	// we are still in class Mesh::ZeroDim

	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                   );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                              );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ThisMeshIsPositive &                                 );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );

	// we are still in class Mesh::ZeroDim

	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &              );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                 );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                 );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                 );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                  );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                     );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                     );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                     );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                    );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                       );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                       );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                       );

	// we are still in class Mesh::ZeroDim

	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &        );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &        );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &   );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &    );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &            );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &             );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &       );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ThisMeshIsPositive &                                         );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &              );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &              );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &         );

	// we are still in class Mesh::ZeroDim

	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                   );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                              );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ThisMeshIsPositive &                                 );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );

	#ifndef NDEBUG
	// attribute  name  inherited from Mesh::Core
	// std::string get_name()  defined by Mesh::Core
	void print_everything ();  // virtual from Mesh::Core
	#endif

};  // end of class Mesh::ZeroDim

//-----------------------------------------------------------------------------//


class Mesh::NotZeroDim : public Mesh::Core

// represents a positive mesh of dimension >= 1

// abstract class, specialized in Mesh::Connected::***Dim,
//   Mesh::MultiplyConnected::***Dim, Mesh::Fuzzy, Mesh::STSI
// see paragraph 11.4 in the manual

{	public :

	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// Cell::Positive * cell_enclosed inherited from Mesh::Core

	inline NotZeroDim ( const tag::OfDimension &, const size_t d , const tag::MinusOne &,
                      const tag::OneDummyWrapper &                                      )
	:	Mesh::Core ( tag::of_dim, d, tag::minus_one, tag::one_dummy_wrapper )  { }

	// bool dispose_query ( )  defined by tag::Util::Core ifdef COLLECT_CM
	
	// size_t get_dim_plus_one ( )  stays pure virtual from Mesh::Core

	// four versions of number_of  stay pure virtual from Mesh::Core

	// first_vertex, last_vertex, first_segment and last_segment
	// are defined by Mesh::Core, execution forbidden

	// cell_in_front_of and cell_behind ( tag::seen_from )
	// defined by Mesh::Core, execution forbidden

	// the thirty-two methods below are virtual from Mesh::Core
	// called from Cell::****tive::***::add_to_mesh and Cell::****tive::***::remove_from_mesh
	void add_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsNotBdry & );
	void add_pos_seg
		( Cell::Positive::Segment *, const tag::MeshIsNotBdry &, const tag::DoNotBother & );
	void add_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsBdry & );
	void add_pos_seg
		( Cell::Positive::Segment *, const tag::MeshIsBdry &, const tag::DoNotBother & );
	void remove_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsNotBdry & );
	void remove_pos_seg
		( Cell::Positive::Segment *, const tag::MeshIsNotBdry &, const tag::DoNotBother & );
	void remove_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsBdry & );
	void remove_pos_seg
		( Cell::Positive::Segment *, const tag::MeshIsBdry &, const tag::DoNotBother & );
	void add_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsNotBdry & );
	void add_neg_seg
		( Cell::Negative::Segment *, const tag::MeshIsNotBdry &, const tag::DoNotBother & );
	void add_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsBdry & );
	void add_neg_seg
		( Cell::Negative::Segment *, const tag::MeshIsBdry &, const tag::DoNotBother & );
	void remove_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsNotBdry & );
	void remove_neg_seg
		( Cell::Negative::Segment *, const tag::MeshIsNotBdry &, const tag::DoNotBother & );
	void remove_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsBdry & );
	void remove_neg_seg
		( Cell::Negative::Segment *, const tag::MeshIsBdry &, const tag::DoNotBother & );
	void add_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsNotBdry & );
	void add_pos_hd_cell
		( Cell::Positive::HighDim *, const tag::MeshIsNotBdry &, const tag::DoNotBother & );
	void add_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsBdry & );
	void add_pos_hd_cell
		( Cell::Positive::HighDim *, const tag::MeshIsBdry &, const tag::DoNotBother & );
	void remove_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsNotBdry & );
	void remove_pos_hd_cell
		( Cell::Positive::HighDim *, const tag::MeshIsNotBdry &, const tag::DoNotBother & );
	void remove_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsBdry & );
	void remove_pos_hd_cell
		( Cell::Positive::HighDim *, const tag::MeshIsBdry &, const tag::DoNotBother & );
	void add_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsNotBdry & );
	void add_neg_hd_cell
		( Cell::Negative::HighDim *, const tag::MeshIsNotBdry &, const tag::DoNotBother & );
	void add_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsBdry & );
	void add_neg_hd_cell
		( Cell::Negative::HighDim *, const tag::MeshIsBdry &, const tag::DoNotBother & );
	void remove_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsNotBdry & );
	void remove_neg_hd_cell
		( Cell::Negative::HighDim *, const tag::MeshIsNotBdry &, const tag::DoNotBother & );
	void remove_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsBdry & );
	void remove_neg_hd_cell
		( Cell::Negative::HighDim *, const tag::MeshIsBdry &, const tag::DoNotBother & );

	Mesh::Core * build_deep_copy ( );  // virtual from Mesh::Core, execution forbidden for now

	// two versions of add_to_my_cells  stay pure virtual from Cell::Core
	// two versions of remove_from_my_cells  stay pure virtual from Cell::Core

	// two versions of  closed_loop  stay pure virtual from Mesh::Core

	// iterators stay pure virtual from Mesh::Core

	#ifndef NDEBUG
	// attribute  name  inherited from Mesh::Core
	// std::string get_name()  defined by Mesh::Core
	// void print_everything ()  stays pure virtual from Mesh::Core
	#endif

};  // end of class Mesh::NotZeroDim

//-----------------------------------------------------------------------------//


class Mesh::Connected::OneDim : public Mesh::NotZeroDim

// represents a connected positive mesh of dimension 1
// a chain of segments, either open or closed

{	public :

	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// Cell::Positive * cell_enclosed inherited from Mesh::Core

	size_t nb_of_segs;
	// useful for quickly answering to 'number_of'

	// here we use Cell wrappers as pointers
	Cell first_ver, last_ver;
	// first_ver is negative, last_ver is positive
	// if last_ver == first_ver .reverse(), the mesh is a closed loop
	
	inline OneDim ( const tag::OneDummyWrapper & )
	:	Mesh::NotZeroDim ( tag::of_dimension, 2, tag::minus_one, tag::one_dummy_wrapper ),
		first_ver ( tag::non_existent ), last_ver ( tag::non_existent )
  { }

	inline OneDim
	( const tag::With &, size_t n, const tag::Segments &, const tag::OneDummyWrapper & )
	:	Mesh::NotZeroDim ( tag::of_dimension, 2, tag::minus_one, tag::one_dummy_wrapper ),
		nb_of_segs { n }, first_ver ( tag::non_existent ), last_ver ( tag::non_existent )
	{	}

	virtual ~OneDim ();
	
	// bool dispose_query ( )  defined by tag::Util::Core ifdef COLLECT_CM
	
	size_t get_dim_plus_one ( );  // virtual from Mesh::Core

	size_t number_of ( const tag::Vertices & );  // virtual from Mesh::Core
	size_t number_of ( const tag::Segments & );  // virtual from Mesh::Core
	size_t number_of ( const tag::CellsOfDim &, const size_t d );  // virtual from Mesh::Core
	size_t number_of ( const tag::CellsOfMaxDim & );  // virtual from Mesh::Core

	Cell first_vertex ( );  // virtual from Mesh::Core, here overridden
	Cell last_vertex ( );  // virtual from Mesh::Core, here overridden
	Cell first_segment ( );  // virtual from Mesh::Core, here overridden
	Cell last_segment ( );  // virtual from Mesh::Core, here overridden

	// cell_in_front_of and cell_behind ( tag::seen_from )
	// defined by Mesh::Core, execution forbidden

	// thirty-two methods add_*** and remove_***  defined by Mesh::NotZeroDim
	// called from Cell::****tive::***::add_to_mesh and Cell::****tive::***::remove_from_mesh

	// Mesh::Core * build_deep_copy ( )  defined by Mesh::NotZeroDim, execution forbidden for now

	std::list<Cell>::iterator add_to_my_cells ( Cell::Core * const, const size_t );
	// virtual from Cell::Core, here returns garbage, updates nb_of_segs, first_ver, last_ver

	std::list<Cell>::iterator add_to_my_cells
		( Cell::Core * const, const size_t, const tag::DoNotBother & );
	// virtual from Cell::Core, here returns garbage

	void remove_from_my_cells ( Cell::Core * const, const size_t, std::list<Cell>::iterator );
	// virtual from Cell::Core, here only updates nb_of_segs, first_ver, last_ver
	
	void remove_from_my_cells
		( Cell::Core * const, const size_t, std::list<Cell>::iterator, const tag::DoNotBother & );
	// virtual from Cell::Core, here does nothing
	
	void closed_loop ( const Cell & ver );  // virtual from Mesh::Core
	// sets both first_ver and last_ver to 'ver'

	void closed_loop ( const Cell & ver, size_t );  // virtual from Mesh::Core
	// sets both first_ver and last_ver to 'ver' and number of segments

	// iterators are virtual from Mesh::Core and are defined in iterator.cpp
	
	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	// iterate over the two vertices, first base (negative) then tip (positive)

	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	// iterate over the two vertices, first tip (positive) then base (negative)

	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	// iterate over the two vertices, first base then tip (both positive)

	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	// we iterate over the two vertices, first tip then base (both positive)

	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                        );
	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                        );
	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &                   );

	// we are still in class Mesh::Connected::OneDim

	Mesh::Iterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                        );
	Mesh::Iterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                        );
	Mesh::Iterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &                   );

	// we are still in class Mesh::Connected::OneDim

	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::AsTheyAre &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &        );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &        );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &   );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &            );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &            );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &       );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive &               );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );

	// we are still in class Mesh::Connected::OneDim
	// twelve iterators below simply produce iterator over segments
	
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );

	// we are still in class Mesh::Connected::OneDim

	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &              );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &              );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                 );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                 );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                 );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                  );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                     );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                     );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                     );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                    );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                       );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                       );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                       );

	// we are still in class Mesh::Connected::OneDim

	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &        );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &        );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &   );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &    );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &            );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &             );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &       );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ThisMeshIsPositive &                                         );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &              );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &              );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &         );

	// we are still in class Mesh::Connected::OneDim

	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                   );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                              );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ThisMeshIsPositive &                                 );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );

	// we are still in class Mesh::Connected::OneDim

	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &              );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                 );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                 );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                 );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                  );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                     );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                     );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                     );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                    );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                       );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                       );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                       );

	// we are still in class Mesh::Connected::OneDim

	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &        );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &        );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &   );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &    );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &            );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &             );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &       );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ThisMeshIsPositive &                                         );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &              );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &              );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &         );

	// we are still in class Mesh::Connected::OneDim

	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                   );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                              );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ThisMeshIsPositive &                                 );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	
	#ifndef NDEBUG
	// attribute  name  inherited from Mesh::Core
	// std::string get_name()  defined by Mesh::Core
	void print_everything ();  // virtual from Mesh::Core
	#endif

}; // end of  class Mesh::Connected::OneDim

//-----------------------------------------------------------------------------//


class Mesh::Connected::HighDim : public Mesh::NotZeroDim

// represents a connected positive mesh of dimension >= 2

{	public :

	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	std::vector < size_t > nb_of_cells;
	// useful for quickly answering to 'number_of'

	// Cell start;
	
	inline HighDim ( const tag::OfDimension &, const size_t dim_p1, const tag::MinusOne &,
	                 const tag::OneDummyWrapper &                                          )
	:	Mesh::NotZeroDim ( tag::of_dimension, dim_p1, tag::minus_one, tag::one_dummy_wrapper )
	{	}

  HighDim ( const tag::Triangle &, const Mesh & AB, const Mesh & BC, const Mesh & CA,
            const tag::OneDummyWrapper &                                              );
	// defined in global.cpp

	inline HighDim ( const tag::Quadrangle &,
		const Mesh & south, const Mesh & east, const Mesh & north, const Mesh & west,
	  const tag::OneDummyWrapper &, const tag::WithTriangles & wt = tag::not_with_triangles );

	inline HighDim ( const tag::Quadrangle &,
		const Cell & SW, const Cell & SE, const Cell & NE, const Cell & NW,
		const size_t m, const size_t n, const tag::OneDummyWrapper &,
									 const tag::WithTriangles & wt = tag::not_with_triangles );

	// bool dispose_query ( )  defined by tag::Util::Core ifdef COLLECT_CM
	
	size_t get_dim_plus_one ( );  // virtual from Mesh::Core

	size_t number_of ( const tag::Vertices & );  // virtual from Mesh::Core
	size_t number_of ( const tag::Segments & );  // virtual from Mesh::Core
	size_t number_of ( const tag::CellsOfDim &, const size_t d );  // virtual from Mesh::Core
	size_t number_of ( const tag::CellsOfMaxDim & );  // virtual from Mesh::Core

	// first_vertex, last_vertex, first_segment and last_segment
	// are defined by Mesh::Core, execution forbidden

	// private :
	
	// cell_in_front_of  and  cell_behind  defined by Mesh::Core, execution forbidden

	// add_to_my_cells ( Cell::Core *, size_t )  defined by Cell::Core, returns garbage
	// remove_from_my_cells defined by Cell::Core, does nothing

	inline Mesh::Iterator iterator ( const tag::OverCellsOfDim &, const size_t d );
	
	#ifndef NDEBUG
	// attribute  name  inherited from Mesh::Core
	// std::string get_name()  defined by Mesh::Core
	void print_everything ();  // virtual from Mesh::Core
	#endif
	
}; // end of  class Mesh::Connected::HighDim

//-----------------------------------------------------------------------------//


class Mesh::MultiplyConnected::OneDim : public Mesh::NotZeroDim

// represents a positive mesh of dimension 1 with several connected components

{	public :

	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// Cell::Positive * cell_enclosed inherited from Mesh::Core

	std::vector < size_t > nb_of_cells;
	// useful for quickly answering to 'number_of'

	// keep a list or a vector of starting/ending points

	inline OneDim ( const tag::OneDummyWrapper & )
		:	Mesh::NotZeroDim ( tag::of_dimension, 2, tag::minus_one, tag::one_dummy_wrapper )
	{ }

  OneDim ( const tag::Segment &, Cell::Negative::Vertex * A,
		Cell::Positive::Vertex * B, const tag::DividedIn &, const size_t n, const tag::OneDummyWrapper & );
	// defined in global.cpp
	
	// bool dispose_query ( )  defined by tag::Util::Core ifdef COLLECT_CM
	
	size_t get_dim_plus_one ( );  // virtual from Mesh::Core

	// size_t number_of ( const tag::Vertices & );  // virtual from Mesh::Core
	// size_t number_of ( const tag::Segments & );  // virtual from Mesh::Core
	// size_t number_of ( const tag::CellsOfDim &, const size_t d );  // virtual from Mesh::Core
	// size_t number_of ( const tag::CellsOfMaxDim & );  // virtual from Mesh::Core

	// first_vertex, last_vertex, first_segment and last_segment
	// are defined by Mesh::Core, execution forbidden

	// private:
	
	// cell_in_front_of  and  cell_behind  defined by Mesh::Core, execution forbidden

	// add_to_my_cells ( Cell::Core *, size_t )  defined by Cell::Core, returns garbage
	// remove_from_my_cells defined by Cell::Core, does nothing

	#ifndef NDEBUG
	// attribute  name  inherited from Mesh::Core
	// std::string get_name()  defined by Mesh::Core
	// void print_everything ();  // virtual from Mesh::Core
	#endif

}; // end of  class Mesh::MultiplyConnected::OneDim

//-----------------------------------------------------------------------------//


class Mesh::MultiplyConnected::HighDim : public Mesh::NotZeroDim

// represents a positive mesh of dimension >= 2 with several connected components

{	public :

	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// Cell::Positive * cell_enclosed inherited from Mesh::Core

	std::vector < size_t > nb_of_cells;
	// useful for quickly answering to 'number_of'

	// keep a list or a vector of starting cells
	
	inline HighDim ( const tag::OfDimension &, const size_t dim_p1, const tag::MinusOne &, const tag::OneDummyWrapper & )
		:	Mesh::NotZeroDim ( tag::of_dimension, dim_p1, tag::minus_one, tag::one_dummy_wrapper )
	{	}

	// bool dispose_query ( )  defined by tag::Util::Core ifdef COLLECT_CM
	
	size_t get_dim_plus_one ( );  // virtual from Mesh::Core

	// size_t number_of ( const tag::Vertices & );  // virtual from Mesh::Core
	// size_t number_of ( const tag::Segments & );  // virtual from Mesh::Core
	// size_t number_of ( const tag::CellsOfDim &, const size_t d );  // virtual from Mesh::Core
	// size_t number_of ( const tag::CellsOfMaxDim & );  // virtual from Mesh::Core

	// first_vertex, last_vertex, first_segment and last_segment
	// are defined by Mesh::Core, execution forbidden

	void build_rectangle ( const Mesh & south, const Mesh & east,
		const Mesh & north, const Mesh & west, bool cut_rectangles_in_half );
	// defined in global.cpp

	static bool is_positive ( );
	static Mesh reverse ( Mesh::Core * core );

	// private :
	
	// cell_in_front_of  and  cell_behind  defined by Mesh::Core, execution forbidden

	// add_to_my_cells ( Cell::Core *, size_t )  defined by Cell::Core, returns garbage
	// remove_from_my_cells defined by Cell::Core, does nothing

	inline Mesh::Iterator iterator ( const tag::OverCellsOfDim &, const size_t d );
	
	#ifndef NDEBUG
	// attribute  name  inherited from Mesh::Core
	// std::string get_name()  defined by Mesh::Core
	// void print_everything ();  // virtual from Mesh::Core
	#endif
	
}; // end of  class Mesh::MultiplyConnected::HighDim

//-----------------------------------------------------------------------------//


class Mesh::Fuzzy : public Mesh::NotZeroDim

// represents a positive mesh, unordered (includes 1D meshes)

// roughly speaking, a mesh is a collection of cells of the same dimension
// however, for efficiency purposes, we keep lists of cells of lower
// dimension as well; that's the purpose of the 'cells' vector :
// to keep lists of cells indexed by their dimension

{	public :
	
	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// Cell::Positive * cell_enclosed inherited from Mesh::Core

	// the 'cells' attribute holds lists of cells of 'this' mesh, indexed by their dimension
	// for maximum dimension, the cells are oriented
	// for lower dimension, the cells are always positive

	// here we use Cell wrappers as pointers because we want a mesh to keep its cells alive
	std::vector < std::list < Cell > > cells;

	inline Fuzzy ( const tag::OfDimension &, const size_t dim_p1, const tag::MinusOne &,
	               const tag::OneDummyWrapper &                                         )
	:	Mesh::NotZeroDim ( tag::of_dimension, dim_p1, tag::minus_one, tag::one_dummy_wrapper ),
		cells ( dim_p1 )
	{	assert ( dim_p1 > 1 );  }

	virtual ~Fuzzy ();
	
	// bool dispose_query ( )  defined by tag::Util::Core ifdef COLLECT_CM
	
	size_t get_dim_plus_one ( );  // virtual from Mesh::Core

	size_t number_of ( const tag::Vertices & );  // virtual from Mesh::Core
	size_t number_of ( const tag::Segments & );  // virtual from Mesh::Core
	size_t number_of ( const tag::CellsOfDim &, const size_t d );  // virtual from Mesh::Core
	size_t number_of ( const tag::CellsOfMaxDim & );  // virtual from Mesh::Core

	// first_vertex, last_vertex, first_segment and last_segment
	// are defined by Mesh::Core, execution forbidden

	// cell_in_front_of and cell_behind ( tag::seen_from )
	// defined by Mesh::Core, execution forbidden
	
	void build_rectangle ( const Mesh & south, const Mesh & east,  // defined in global.cpp
		const Mesh & north, const Mesh & west, bool cut_rectangles_in_half );
		
	// cell_in_front_of  and  cell_behind  defined by Mesh::Core, execution forbidden

	// thirty-two methods add_*** and remove_***  defined by Mesh::NotZeroDim
	// called from Cell::****tive::***::add_to_mesh and Cell::****tive::***::remove_from_mesh

	// Mesh::Core * build_deep_copy ( )  defined by Mesh::NotZeroDim, execution forbidden for now

	// add a cell to 'this->cells [d]' list, return iterator into that list
	std::list < Cell > ::iterator add_to_my_cells ( Cell::Core * const cll, const size_t d );
	std::list < Cell > ::iterator add_to_my_cells
		( Cell::Core * const cll, const size_t d, const tag::DoNotBother & );
	// virtual from Cell::Core, later overriden by Mesh::STSI
	
	// remove a cell from 'this->cells[d]' list using the provided iterator
	void remove_from_my_cells ( Cell::Core * const, const size_t d, std::list < Cell > ::iterator );
	void remove_from_my_cells
		( Cell::Core * const, const size_t d, std::list < Cell > ::iterator, const tag::DoNotBother & );
	// virtual from Cell::Core, later overriden by Mesh::STSI
	
	void closed_loop ( const Cell & ver );
	void closed_loop ( const Cell & ver, size_t );
	// virtual from Mesh::Core, here execution forbidden

	// we are still in class Mesh::Fuzzy
	// iterators are virtual from Mesh::Core and are defined in iterator.cpp

	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                        );
	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                        );
	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &                   );

	// we are still in class Mesh::Fuzzy

	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                        );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                        );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &                   );

	// we are still in class Mesh::Fuzzy

	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::AsTheyAre &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &        );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &        );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &   );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &            );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &            );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &       );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive &               );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );

	// we are still in class Mesh::Fuzzy
	
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	
	// we are still in class Mesh::Fuzzy

	Mesh::Iterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &              );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &              );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                 );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                 );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                 );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                  );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                     );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                     );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                     );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                    );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                       );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                       );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                       );

	// we are still in class Mesh::Fuzzy

	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &        );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &        );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &   );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &    );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &            );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &             );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &       );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ThisMeshIsPositive &                                         );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &              );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &              );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &         );

	// we are still in class Mesh::Fuzyy

	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                   );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                              );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ThisMeshIsPositive &                                 );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );

	// we are still in class Mesh::Fuzzy

	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &              );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                 );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                 );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                 );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                  );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                     );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                     );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                     );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                    );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                       );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                       );
	Mesh::Iterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                       );

	// we are still in class Mesh::Fuzzy

	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &        );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &        );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &   );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &    );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &            );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &             );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &       );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ThisMeshIsPositive &                                         );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &              );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &              );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &         );

	// we are still in class Mesh::Fuzzy

	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                   );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                                                      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                                      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
	  const tag::ThisMeshIsPositive &                                                      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &,
	  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
	  const tag::ThisMeshIsPositive &                              );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ThisMeshIsPositive &                                 );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &      );
	Mesh::Iterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	
	#ifndef NDEBUG
	// attribute  name  inherited from Mesh::Core
	// std::string get_name()  defined by Mesh::Core
	void print_everything ();  // virtual from Mesh::Core
	#endif
	
}; // end of  class Mesh::Fuzzy

//-----------------------------------------------------------------------------//


class Mesh::STSI : public Mesh::Fuzzy

// represents a positive mesh, maybe self-touching, maybe self-intersecting
// (includes 1D meshes)

// mainly used for iterators over connected meshes
// also used for progressive mesh generation

{	public :
	
	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// Cell::Positive * cell_enclosed  inherited from Mesh::Core

	// 'cells' inherited from Mesh::Fuzzy

	// in 'singular' we keep pair of adjacent cells
	// the common face of such a pair is a singular face
	// that is, a face where the mesh touches itself
	std::vector < std::pair < Cell, Cell > > singular;
	
	inline STSI ( const tag::OfDimension &, const size_t dim_p1, const tag::MinusOne &,
                const tag::OneDummyWrapper &                                         )
	:	Mesh::Fuzzy ( tag::of_dimension, dim_p1, tag::minus_one, tag::one_dummy_wrapper )
	{	}

	virtual ~STSI ();

	// bool dispose_query ( )  defined by tag::Util::Core ifdef COLLECT_CM
	
	// size_t get_dim_plus_one ( )  defined in Mesh::Fuzzy
	// size_t number_of ( const tag::Vertices & )  defined in Mesh::Fuzzy
	// size_t number_of ( const tag::Segments & )  defined in Mesh::Fuzzy
	// size_t number_of ( const tag::CellsOfDim &, const size_t d )  defined in Mesh::Fuzzy
	// size_t number_of ( const tag::CellsOfMaxDim & )  defined in Mesh::Fuzzy

	// first_vertex, last_vertex, first_segment and last_segment
	// are defined by Mesh::Core, execution forbidden

	void build_rectangle ( const Mesh & south, const Mesh & east,
		const Mesh & north, const Mesh & west, bool cut_rectangles_in_half );
	// defined in global.cpp
	
	static Mesh reverse ( Mesh::Core * core );

	// private :

	// the four methods below are defined in Cell::Core, here overriden
	virtual Cell::Core * cell_in_front_of
	( const Cell::Core * face_p, const tag::SeenFrom &, const Cell neighbour,
	  const tag::SurelyExists & se = tag::surely_exists                   ) const override;
	virtual Cell::Core * cell_in_front_of
	( const Cell::Core * face_p, const tag::SeenFrom &,
    const Cell neighbour, const tag::MayNotExist & ) const override;
	virtual Cell::Core * cell_behind
	( const Cell::Core * face_p, const tag::SeenFrom &, const Cell neighbour,
	  const tag::SurelyExists & se = tag::surely_exists                   ) const override;
	virtual Cell::Core * cell_behind
	( const Cell::Core * face_p, const tag::SeenFrom &,
    const Cell neighbour, const tag::MayNotExist & ) const override;

	// add a cell to 'this->cells [d]' list, return iterator into that list
	virtual std::list<Cell>::iterator add_to_my_cells
	( Cell::Core * const cll, const size_t d ) override;
	// virtual from Cell::Core, defined by Mesh::Fuzzy, here overriden

	// remove a cell from 'this->cells [d]' list using the provided iterator
	virtual void remove_from_my_cells
	( Cell::Core * const, const size_t d, std::list<Cell>::iterator );
	// virtual from Cell::Core, defined by Mesh::Fuzzy, here overriden

	// iterators are virtual from Mesh::Core and are defined in iterator.cpp

	inline Mesh::Iterator iterator ( const tag::OverCellsOfDim &, const size_t d );
	
	#ifndef NDEBUG
	// attribute  name  inherited from Mesh::Core
	// std::string get_name()  defined by Mesh::Core
	void print_everything ();  // virtual from Mesh::Core
	#endif
	
}; // end of  class Mesh::STSI

//-----------------------------------------------------------------------------//
//-----------------------------------------------------------------------------//


inline void Mesh::set_max_dim ( const size_t d )  // static
// invoked at the beginning of mesh.cpp
// see paragraph 11.7 in the manual
{	Mesh::maximum_dimension_plus_one = d + 1;
	Cell::Positive::double_heap_size .resize ( maximum_dimension_plus_one, 0 );
	Cell::Negative::double_heap_size .resize ( maximum_dimension_plus_one, 0 );
	Cell::Positive::size_t_heap_size .resize ( maximum_dimension_plus_one, 0 );
	Cell::Negative::size_t_heap_size .resize ( maximum_dimension_plus_one, 0 );
	Cell::Positive::short_int_heap_size .resize ( maximum_dimension_plus_one, 0 );
	Cell::Negative::short_int_heap_size .resize ( maximum_dimension_plus_one, 0 );  }

//-----------------------------------------------------------------------------//


inline tag::Util::CellCore::CellCore
( const tag::OfDimension &, const size_t d,  // for positive cells
  const tag::HasNoReverse &, const tag::OneDummyWrapper & )
#ifdef MANIFEM_COLLECT_CM	
:	tag::Util::Core::DelegateDispose
	( & tag::Util::Core::default_dispose_query, tag::one_dummy_wrapper ),
#else  // no MANIFEM_COLLECT_CM
:	tag::Util::Core::Inactive ( tag::one_dummy_wrapper ),
#endif  // MANIFEM_COLLECT_CM	
	double_heap ( Cell::Positive::double_heap_size [d] ),
	size_t_heap ( Cell::Positive::size_t_heap_size [d] ),
	short_int_heap ( Cell::Positive::short_int_heap_size [d], 0 ),
	reverse_attr ( tag::non_existent )
{	// Cell::counter++;
	std::vector < void(*)(Cell::Core*,void*) > & init = Cell::init_pos_cell [d];
	std::vector < void* > & data = Cell::data_for_init_pos [d];
	std::vector < void(*)(Cell::Core*,void*) > ::iterator it_f = init .begin();
	std::vector < void* > ::iterator it_d = data .begin();
	for ( ; it_f != init .end(); it_f++, it_d++ )
	{	assert ( it_d != data .end() );
		(*it_f) ( this, *it_d );       }
	assert ( it_d == data .end() );                                               }


inline tag::Util::CellCore::CellCore
( const tag::OfDimension &, const size_t d,  // for negative cells
  const tag::ReverseOf &, Cell::Core * direct_cell_p, const tag::OneDummyWrapper & )
#ifdef MANIFEM_COLLECT_CM	
:	tag::Util::Core::DelegateDispose
	( & tag::Util::Core::dispose_query_cell_with_reverse, tag::one_dummy_wrapper ),
#else  // no MANIFEM_COLLECT_CM
:	tag::Util::Core::Inactive ( tag::one_dummy_wrapper ),
#endif  // MANIFEM_COLLECT_CM	
	double_heap ( Cell::Negative::double_heap_size [d] ),
	size_t_heap ( Cell::Negative::size_t_heap_size [d] ),
	short_int_heap ( Cell::Negative::short_int_heap_size [d], 0 ),
	reverse_attr ( tag::whose_core_is, direct_cell_p, tag::previously_existing, tag::surely_not_null )
{	// Cell::counter++;
	std::vector < void(*)(Cell::Core*,void*) > & init = Cell::init_neg_cell [d];
	std::vector < void* > & data = Cell::data_for_init_neg [d];
	std::vector < void(*)(Cell::Core*,void*) > ::iterator it_f = init .begin();
	std::vector<void*>::iterator it_d = data .begin();
	for ( ; it_f != init .end(); it_f++, it_d++ )
	{	assert ( it_d != data .end() );
		( *it_f ) ( this, *it_d );       }
	assert ( it_d == data .end() );                                               }

//-----------------------------------------------------------------------------//


inline Cell & Cell::operator= ( const Cell & c )
{
	#ifdef MANIFEM_COLLECT_CM	
	this->dispose_core ( tag::may_be_null );  // old core may be null
	#endif
	this->core = c.core;
	#ifdef MANIFEM_COLLECT_CM	
	if ( c .core ) c .core->nb_of_wrappers++;
	#endif
	return *this;                            }
		

inline Cell & Cell::operator= ( Cell && c )
{
	#ifdef MANIFEM_COLLECT_CM	
	this->dispose_core ( tag::may_be_null );  // old core may be null
	#endif
	this->core = c .core;
	c .core = nullptr;
	return *this;                            }

	
inline Mesh & Mesh::operator= ( const Mesh & c )
{
	#ifdef MANIFEM_COLLECT_CM	
	this->dispose_core ( tag::may_be_null );  // old core may be null
	#endif
	this->core = c .core;
	this->is_pos = c .is_pos;
	#ifdef MANIFEM_COLLECT_CM	
	if ( c .core ) c .core->nb_of_wrappers++;
	#endif
	return *this;                            }


inline Mesh & Mesh::operator= ( Mesh && c )
{
	#ifdef MANIFEM_COLLECT_CM	
	this->dispose_core ( tag::may_be_null );  // old core may be null
	#endif
	this->core = c .core;
	this->is_pos = c .is_pos;
	c .core = nullptr;
	return *this;                            }
	
//------------------------------------------------------------------------------------------------//


inline void Mesh::copy_all_cells_to ( Mesh & msh ) const
{	assert ( this->exists() );
	assert ( msh.exists() );
	Mesh::Iterator it = this->iterator ( tag::over_cells_of_max_dim );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell cll = *it;  cll.add_to_mesh ( msh );  }                      }


inline Mesh::Mesh ( const tag::WhoseCoreIs &, Mesh::Core * c, const tag::PreviouslyExisting &,
                    const tag::IsNegative &, const tag::CellsSurelyExist &                    )
// builds a negative mesh from a positive one, assuming that reverse cells exist
// used in Cell::boundary and in Mesh::Mesh below

: Mesh ( tag::whose_core_is, c, tag::previously_existing,
         tag::is_negative, tag::do_not_build_cells        )

#ifndef NDEBUG
{	assert ( c );
	// check that all cells have reverse
	Mesh::Iterator it = this->iterator  // as they are : oriented
		( tag::over_cells_of_max_dim, tag::as_they_are );
	for ( it .reset() ; it .in_range(); it++ )
		assert ( ( *it ) .reverse ( tag::may_not_exist ) .exists() );  }
#else
{	}
#endif


inline Mesh::Mesh ( const tag::WhoseCoreIs &, Mesh::Core * c, const tag::PreviouslyExisting &,
                    const tag::IsNegative &, const tag::BuildCellsIfNec & b                   )
// builds a negative mesh from a positive one, creating reverse cells if necessary
// used in Mesh::Positive::reverse
	
: Mesh ( tag::whose_core_is, c, tag::previously_existing,
         tag::is_negative, tag::do_not_build_cells        )

{	Mesh::Iterator it ( tag::whose_core_is, c->iterator  // as they are : oriented
		( tag::over_cells_of_max_dim, tag::as_they_are, tag::this_mesh_is_positive ) );
	for ( it .reset() ; it .in_range(); it++ )
		( *it ) .reverse ( tag::build_if_not_exists );                                  }


inline Mesh::Mesh ( const tag::Fuzzy &, const tag::OfDimension &, const size_t d,
                    const tag::IsPositive & ispos                                 )
// by default, ispos = tag::is_positive, so may be called with only three arguments
:	Mesh ( tag::whose_core_is,
	 	     new Mesh::Fuzzy ( tag::of_dimension, d+1, tag::minus_one, tag::one_dummy_wrapper ),
	       tag::freshly_created, tag::is_positive                                             )
{	}


inline Mesh::Mesh ( const tag::Fuzzy &, const tag::OfDimension &, const size_t d,
                    const tag::MinusOne &, const tag::IsPositive & ispos          )
// by default, ispos = tag::is_positive, so may be called with only four arguments
:	Mesh ( tag::whose_core_is,
	 	     new Mesh::Fuzzy ( tag::of_dimension, d, tag::minus_one, tag::one_dummy_wrapper ),
	       tag::freshly_created, tag::is_positive                                           )
{	}
	

inline Mesh::Mesh ( const tag::DeepCopyOf &, const Mesh & msh )
:	Mesh ( tag::whose_core_is, msh .core->build_deep_copy(),
	       tag::freshly_created, tag::is_positive           )
// the core is built with tag::one_dummy_wrapper
// we must ensure that it does not meet any other wrapper before 'this' one
// in other words, we must ensure that the core is freshly created
// and eager to meet its parent ('this' wrapper)
{	}


inline Mesh::Mesh ( const tag::DeepCopyOf &, const Mesh & msh, const tag::Fuzzy & )
:	Mesh ( tag::fuzzy, tag::of_dim, msh.core->get_dim_plus_one(),
         tag::minus_one, tag::is_positive                       )
{	msh .copy_all_cells_to ( *this );  }

//------------------------------------------------------------------------------------------------//


inline Mesh::Mesh ( const tag::Segment &, const Cell & A, const Cell & B,
                    const tag::DividedIn &, const size_t n                )
: Mesh ( tag::whose_core_is,
         new Mesh::Connected::OneDim ( tag::with, n, tag::segments, tag::one_dummy_wrapper ),
         tag::freshly_created, tag::is_positive                                               )
{	assert ( not A .is_positive() );
	assert ( B .is_positive() );
	this->build ( tag::segment, A, B, tag::divided_in, n );
	Mesh::Connected::OneDim * this_core = tag::Util::assert_cast
		< Mesh::Core*, Mesh::Connected::OneDim* > ( this->core );
	this_core->first_ver = A;
	this_core->last_ver = B;                                      }


inline Mesh::Mesh ( const tag::Segment &, const Cell & A, const Cell & B,
                    const tag::DividedIn &, const size_t n,
                    const tag::Winding &, const tag::Util::Action & s )
// due to the winding, A.reverse may be equal to B (mesh may be a closed loop)
: Mesh ( tag::whose_core_is,
         new Mesh::Connected::OneDim ( tag::with, n, tag::segments, tag::one_dummy_wrapper ),
         tag::freshly_created, tag::is_positive                                               )
{	assert ( not A .is_positive() );
	assert ( B .is_positive() );
	this->build ( tag::segment, A, B, tag::divided_in, n, tag::winding, s );
	Mesh::Connected::OneDim * this_core = tag::Util::assert_cast
		< Mesh::Core*, Mesh::Connected::OneDim* > ( this->core );
	this_core->first_ver = A;
	this_core->last_ver = B;                                              }

//------------------------------------------------------------------------------------------------//


inline Mesh::Mesh
( const tag::Triangle &, const Mesh & AB, const Mesh & BC, const Mesh & CA )
: Mesh ( tag::whose_core_is,
         new Mesh::Fuzzy ( tag::of_dimension, 3, tag::minus_one, tag::one_dummy_wrapper ),
         tag::freshly_created, tag::is_positive                                           )
{	this->build ( tag::triangle, AB, BC, CA );  }


inline Mesh::Mesh
( const tag::Triangle &, const Mesh & AB, const Mesh & BC, const Mesh & CA, const tag::Winding & )

// the tag::winding provides no specific information,
// it just warns maniFEM that we are on a quotient manifold
// and that it must take winding segments into account
// specific information about winding numbers is included in the three segments

: Mesh ( tag::whose_core_is,
         new Mesh::Fuzzy ( tag::of_dimension, 3, tag::minus_one, tag::one_dummy_wrapper ),
         tag::freshly_created, tag::is_positive                                           )

{	this->build ( tag::triangle, AB, BC, CA, tag::winding );  }


inline Mesh::Mesh
( const tag::Triangle &, const Mesh & AB, const Mesh & BC, const Mesh & CA,
  const tag::Winding &, const tag::Singular &, const Cell & O               )

// the tag::winding provides no specific information,
// it just warns maniFEM that we are on a quotient manifold
// and that it must take winding segments into account
// specific information about winding numbers is included in the three segments

// cell O is special, it's like the vertex of a cone
// cell O must be the extremity of two of the three segments provided as arguments

: Mesh ( tag::whose_core_is,
         new Mesh::Fuzzy ( tag::of_dimension, 3, tag::minus_one, tag::one_dummy_wrapper ),
         tag::freshly_created, tag::is_positive                                           )

{	this->build ( tag::triangle, AB, BC, CA, tag::winding, tag::singular, O );  }

//------------------------------------------------------------------------------------------------//


inline Mesh::Mesh ( const tag::Quadrangle &, const Mesh & south,
                    const Mesh & east, const Mesh & north, const Mesh & west,
                    const tag::WithTriangles & wt                             )

// 'wt' defaults to 'tag::not_with_triangles',
// so constructor can be called with only five arguments

: Mesh ( tag::whose_core_is,
         new Mesh::Fuzzy ( tag::of_dimension, 3, tag::minus_one, tag::one_dummy_wrapper ),
         tag::freshly_created, tag::is_positive                                           )

{	this->build ( tag::quadrangle, south, east, north, west, wt != tag::not_with_triangles );  }


inline Mesh::Mesh ( const tag::Quadrangle &, const Mesh & south,
                    const Mesh & east, const Mesh & north, const Mesh & west,
                    const tag::Winding &,
                    const tag::WithTriangles & wt                             )

// 'wt' defaults to 'tag::not_with_triangles',
// so constructor can be called with only six arguments

// the tag::winding provides no specific information,
// it just warns maniFEM that we are on a quotient manifold
// and that it must take winding segments into account
// specific information about winding numbers is included in the four segments

: Mesh ( tag::whose_core_is,
         new Mesh::Fuzzy ( tag::of_dimension, 3, tag::minus_one, tag::one_dummy_wrapper ),
         tag::freshly_created, tag::is_positive                                           )

{	this->build ( tag::quadrangle, south, east, north, west,
                wt != tag::not_with_triangles, tag::winding  );  }


inline Mesh::Mesh ( const tag::Quadrangle &, const Mesh & south,
                    const Mesh & east, const Mesh & north, const Mesh & west,
                    const tag::WithTriangles & wt, const tag::Winding &         )

// the tag::winding provides no specific information,
// it just warns maniFEM that we are on a quotient manifold
// and that it must take winding segmtnts into account
// specific information about winding numbers is included in the four segments

: Mesh ( tag::whose_core_is,
         new Mesh::Fuzzy ( tag::of_dimension, 3, tag::minus_one, tag::one_dummy_wrapper ),
         tag::freshly_created, tag::is_positive                                           )

{	assert ( wt == tag::with_triangles );
	this->build ( tag::quadrangle, south, east, north, west, true, tag::winding );  }


inline Mesh::Mesh ( const tag::Quadrangle &, const Cell & SW, const Cell & SE,
                    const Cell & NE, const Cell & NW, const size_t m, const size_t n,
                    const tag::WithTriangles & wt                                     )

// 'wt' defaults to 'tag::not_with_triangles',
// so constructor can be called with only seven arguments

: Mesh ( tag::whose_core_is,
         new Mesh::Fuzzy ( tag::of_dimension, 3, tag::minus_one, tag::one_dummy_wrapper ),
         tag::freshly_created, tag::is_positive                                           )

{	Mesh south ( tag::segment, SW.reverse ( tag::build_if_not_exists ), SE, tag::divided_in, m );
	Mesh east  ( tag::segment, SE.reverse ( tag::build_if_not_exists ), NE, tag::divided_in, n );
	Mesh north ( tag::segment, NE.reverse ( tag::build_if_not_exists ), NW, tag::divided_in, m );
	Mesh west  ( tag::segment, NW.reverse ( tag::build_if_not_exists ), SW, tag::divided_in, n );
	this->build ( tag::quadrangle, south, east, north, west, wt != tag::not_with_triangles );     }


inline Mesh::Mesh ( const tag::Join &, const Mesh & m1, const Mesh & m2 )
: Mesh ( tag::join, std::vector { m1, m2 } )  { }

	
inline Mesh::Mesh ( const tag::Join &, const Mesh & m1, const Mesh & m2, const Mesh & m3 )
: Mesh ( tag::join, std::vector { m1, m2, m3 } )  { }


inline Mesh::Mesh
( const tag::Join &, const Mesh & m1, const Mesh & m2, const Mesh & m3, const Mesh & m4 )
: Mesh ( tag::join, std::vector { m1, m2, m3, m4 } )  { }


inline Mesh::Mesh ( const tag::Join &, const Mesh & m1, const Mesh & m2, const Mesh & m3,
                                       const Mesh & m4, const Mesh & m5 )
: Mesh ( tag::join, std::vector { m1, m2, m3, m4, m5 } )  { }


template < typename container >  // static
inline size_t Mesh::join ( Mesh * const that, const container & l )

{	// check the dimensions
	#ifndef NDEBUG
	typename container::const_iterator it0 = l .begin();
	assert ( it0 != l .end() );
	for ( ; it0 != l .end(); it0++ ) assert ( that->dim() == it0->dim() );
	#endif  // DEBUG
	size_t n = 0;
	// sweep over the list of meshes	
	typename container::const_iterator it = l .begin();
	for ( ; it != l .end(); it++ )
	{	Mesh m = *it;  // sweep all cells of m
		n += m.number_of ( tag::cells_of_max_dim );
		Mesh::Iterator itt = m.iterator ( tag::over_cells_of_max_dim );
		for ( itt .reset(); itt .in_range(); itt++ )
		{	Cell cll = *itt;
			cll .add_to_mesh ( *that, tag::do_not_bother );  }          }
	return n;                                                                 }


template < typename container >  // static
inline void Mesh::join_meshes ( Mesh * const that, const container & l )

// if any of the meshes is not a Mesh::Connected::OneDim, 'this' will be Mesh::Fuzzy
{	container ll = l;
	std::deque < Mesh > d;
	typename container::iterator it = ll .begin();
	assert ( it != ll .end() );
	Mesh m = *it;
	Mesh::Connected::OneDim * mm = dynamic_cast < Mesh::Connected::OneDim* > ( m .core );
	if ( mm == nullptr )  goto fuzzy;
	d .push_back ( m );  ll .erase ( it );
	again :
	for ( it = ll .begin(); it != ll .end(); it++ )
	{	m = *it;
		mm = dynamic_cast < Mesh::Connected::OneDim* > ( m .core );
		if ( mm == nullptr )  goto fuzzy;
		Mesh m_front = d .front();
		if ( m_front.first_vertex() .reverse() == m .last_vertex() )
		{  d .push_front ( m );  ll .erase ( it );  goto again;  }
		Mesh m_back = d .back();
		if ( m_back .last_vertex() == m .first_vertex() .reverse() )
		{  d .push_back ( m );  ll .erase ( it );  goto again;  }    }
	// if some meshes have not moved to d (are still in ll), the mesh will be disconnected
	if ( ll .size() ) goto fuzzy;
  mm = new Mesh::Connected::OneDim ( tag::one_dummy_wrapper );
  that->core = mm;
	mm->nb_of_segs = Mesh::join ( that, l );
  mm->first_ver = d .front() .first_vertex();
  mm->last_ver = d .back() .last_vertex();
	return;
	fuzzy :
	that->core = new Mesh::Fuzzy ( tag::of_dim, l .front() .core->get_dim_plus_one(),
	                               tag::minus_one, tag::one_dummy_wrapper             );
	Mesh::join ( that, l );
	return;                                                                                  }


template < typename container >
inline Mesh::Mesh ( const tag::Join &, const container & l )
: Mesh ( tag::non_existent )
{	Mesh::join_meshes ( this, l );  }
	
//-----------------------------------------------------------------------------//


inline Mesh::Mesh ( const tag::Import &, const tag::Msh &, const std::string filename )
:	Mesh ( tag::non_existent )
{	Mesh::import_msh ( this, filename );  }

//-----------------------------------------------------------------------------//


inline Cell::Cell ( const tag::WhoseBoundaryIs &, Mesh & msh )
:	Cell ( tag::whose_core_is,
	       new Cell::Positive::HighDim ( tag::whose_boundary_is, msh, tag::one_dummy_wrapper ),
				 tag::freshly_created                                         )
#ifndef NDEBUG
{	assert ( msh .is_positive() );
	assert ( msh .dim() >= 1 );
	Mesh::STSI * msh_stsi = dynamic_cast < Mesh::STSI* > ( msh .core );
	assert ( msh_stsi == nullptr );                                     }
#else
{	}
#endif


inline Cell::Cell ( const tag::Vertex &, const tag::IsPositive & ispos )
// by default, ispos = tag::is_positive, so may be called with only one argument
:	Cell ( tag::whose_core_is, new Cell::Positive::Vertex ( tag::one_dummy_wrapper ),
         tag::freshly_created                                                       )
{	}


inline Cell::Cell ( const tag::Segment &, const Cell & A, const Cell & B )
: Cell ( tag::whose_core_is, new Cell::Positive::Segment ( A, B, tag::one_dummy_wrapper ),
	       tag::freshly_created                                                             )
{	}


inline Cell::Cell ( const tag::Triangle &, const Cell & AB, const Cell & BC, const Cell & CA )
:	Cell ( tag::whose_core_is, new Cell::Positive::HighDim
	         ( tag::triangle, AB, BC, CA, tag::one_dummy_wrapper ),
	       tag::freshly_created                                     )
{	}


inline Cell::Cell ( const tag::Quadrangle &, const Cell & AB, const Cell & BC,
                                             const Cell & CD, const Cell & DA )
:	Cell ( tag::whose_core_is, new Cell::Positive::HighDim
	         ( tag::quadrangle, AB, BC, CD, DA, tag::one_dummy_wrapper ),
	       tag::freshly_created                                           )
{	}


inline Cell::Cell ( const tag::Parallelogram &, const Cell & AB, const Cell & BC,
                                                const Cell & CD, const Cell & DA )
:	Cell ( tag::whose_core_is, new Cell::Positive::HighDim
	         ( tag::quadrangle, AB, BC, CD, DA, tag::one_dummy_wrapper ),
	       tag::freshly_created                                           )
{	}


inline Cell::Cell ( const tag::Rectangle &, const Cell & AB, const Cell & BC,
                                            const Cell & CD, const Cell & DA )
:	Cell ( tag::whose_core_is, new Cell::Positive::HighDim
	         ( tag::quadrangle, AB, BC, CD, DA, tag::one_dummy_wrapper ),
	       tag::freshly_created                                           )
{	}


inline Cell::Cell ( const tag::Square &, const Cell & AB, const Cell & BC,
                                         const Cell & CD, const Cell & DA )
:	Cell ( tag::whose_core_is, new Cell::Positive::HighDim
	         ( tag::quadrangle, AB, BC, CD, DA, tag::one_dummy_wrapper ),
	       tag::freshly_created                                           )
{	}


inline Cell::Cell ( const tag::Pentagon &, const Cell & AB, const Cell & BC,
                         const Cell & CD, const Cell & DE, const Cell & EA )
:	Cell ( tag::whose_core_is, new Cell::Positive::HighDim
	       ( tag::pentagon, AB, BC, CD, DE, EA, tag::one_dummy_wrapper ),
	       tag::freshly_created                                           )
{	}


inline Cell::Cell ( const tag::Hexagon &, const Cell & AB, const Cell & BC,
                    const Cell & CD, const Cell & DE, const Cell & EF, const Cell & FA )
:	Cell ( tag::whose_core_is, new Cell::Positive::HighDim
	       ( tag::hexagon, AB, BC, CD, DE, EF, FA, tag::one_dummy_wrapper ),
	       tag::freshly_created                                             )
{	}

//-----------------------------------------------------------------------------//


inline size_t Mesh::number_of ( const tag::Vertices & ) const
{	return this->number_of ( tag::cells_of_dim, 0 );  }

inline size_t Mesh::number_of ( const tag::Segments & ) const
{	return this->number_of ( tag::cells_of_dim, 1 );  }

inline size_t Mesh::number_of ( const tag::CellsOfDim &, const size_t d ) const
{	return this->core->number_of ( tag::cells_of_dim, d );  }

inline size_t Mesh::number_of ( const tag::CellsOfMaxDim & ) const
{	return this->core->number_of ( tag::cells_of_max_dim );  }


inline Cell Mesh::first_vertex ( ) const
{	assert ( this->dim() == 1 );
	if ( this->is_positive() ) return this->core->first_vertex();
	// else
	return this->core->last_vertex() .reverse ( tag::surely_exists );  }

inline Cell Mesh::last_vertex ( ) const
{	assert ( this->dim() == 1 );
	if ( this->is_positive() ) return this->core->last_vertex();
	// else
	return this->core->first_vertex() .reverse ( tag::surely_exists );  }


inline Cell Mesh::first_segment ( ) const
{	assert ( this->dim() == 1 );
	if ( this->is_positive() ) return this->core->first_segment();
	// else
	return this->core->last_segment() .reverse ( tag::surely_exists );  }

inline Cell Mesh::last_segment ( ) const
{	assert ( this->dim() == 1 );
	if ( this->is_positive() ) return this->core->last_segment();
	// else
	return this->core->first_segment() .reverse ( tag::surely_exists );        }


inline bool Cell::belongs_to
( const Mesh & msh, const tag::SameDim &, const tag::Oriented & ) const
// third argument defaults to tag::not_oriented, so method can be called with two arguments only

{	if ( msh.is_positive() )
		return this->core->belongs_to ( msh .core, tag::same_dim, tag::oriented );
	// mesh is negative
	Cell rev_this = this->reverse ( tag::may_not_exist );
	if ( rev_this .exists() )
		return rev_this .core->belongs_to ( msh .core, tag::same_dim, tag::oriented );
	return false;                                                                   }


inline bool Cell::belongs_to
( const Mesh & msh, const tag::SameDim &, const tag::NotOriented & ) const

{	return this->core->belongs_to ( msh.core, tag::same_dim, tag::not_oriented );  }


inline bool Cell::belongs_to
( const Mesh & msh, const tag::CellHasLowDim &, const tag::NotOriented & ) const

{	return this->core->belongs_to ( msh.core, tag::cell_has_low_dim, tag::not_oriented );  }

	
inline bool Cell::belongs_to ( const Mesh & msh, const tag::Oriented & ) const

{	return this->belongs_to ( msh, tag::same_dim, tag::oriented );  }


inline bool Cell::belongs_to ( const Mesh & msh, const tag::NotOriented & ) const

{	if ( this->dim() == msh .dim() )
		return this->core->belongs_to ( msh .core, tag::same_dim, tag::not_oriented );
	// else --
	assert ( this->dim() < msh .dim() );
	return this->core->belongs_to ( msh .core, tag::cell_has_low_dim, tag::not_oriented );  }


inline bool Cell::belongs_to ( const Mesh & msh ) const

{	assert ( this->dim() < msh .dim() );
	// when the dimensions are equal, we require a more specific call
	return this->core->belongs_to ( msh .core, tag::cell_has_low_dim, tag::not_oriented );  }

//-----------------------------------------------------------------------------//


inline Cell Mesh::cell_in_front_of ( const Cell face, const tag::MayNotExist & ) const

// return the cell towards which 'face' is looking
// recall that the faces of a cell are looking outwards

{	Cell face_rev = face .reverse ( tag::may_not_exist );
	if ( not face_rev .exists() ) return Cell ( tag::non_existent );
	else return this->cell_behind ( face_rev, tag::may_not_exist );  }
	

inline Cell Mesh::cell_in_front_of ( const Cell face, const tag::SurelyExists & se ) const

// 'se' defaults to tag::surely_exists, so method may be called with only one argument

// return the cell towards which 'face' is looking
// recall that the faces of a cell are looking outwards

{	Cell face_rev = face .reverse ( tag::surely_exists );
	return this->cell_behind ( face_rev, tag::surely_exists );      }
	

inline Cell Mesh::cell_behind
( const Cell face, const tag::MayNotExist & ) const

// return the cell to which 'face' belongs, non-existent if we are facing the boundary

{	assert ( this->dim() == face .dim() + 1 );
	if ( this->is_positive() )
	{	std::map < Mesh::Core*, Cell > ::const_iterator
			it = face .core->cell_behind_within .find ( this->core );
		if ( it == face .core->cell_behind_within .end() ) return Cell ( tag::non_existent );
			// nothing behind us, we are touching the boundary
		assert ( it->second .exists() );
		return it->second;                                                                    }
		// face .core->cell_behind_within [ this->core ]
	else
	{	Cell face_rev = face .reverse ( tag::surely_exists );
		// we are in a negative mesh, all faces must have reverse
		std::map < Mesh::Core*, Cell >::const_iterator
			it = face_rev .core->cell_behind_within .find ( this->core );
		if ( it == face_rev .core->cell_behind_within .end() )
			return Cell ( tag::non_existent );  // we are facing the boundary
		Cell cll_rev = it->second;
		assert ( cll_rev .exists() );
		return cll_rev .reverse ( tag::surely_exists );                  }  }


inline Cell Mesh::cell_behind ( const Cell face, const tag::SurelyExists & se ) const

// 'se' defaults to tag::surely_exists, so method may be called with only one argument

// return the cell to which 'face' belongs

{	assert ( this->dim() == face .dim() + 1 );
	if ( this->is_positive() )
	{	std::map < Mesh::Core*, Cell > ::const_iterator
			it = face .core->cell_behind_within .find ( this->core );
		assert ( it != face .core->cell_behind_within .end() );
		assert ( it->second .exists() );
		return it->second;                                        }
		// face .core->cell_behind_within [ this->core ]
	else
	{	Cell face_rev = face .reverse ( tag::surely_exists );
		std::map <Mesh::Core*, Cell > ::const_iterator
			it = face_rev .core->cell_behind_within .find ( this->core );
		assert ( it != face_rev .core->cell_behind_within .end() );
		Cell cll_rev = it->second;
		assert ( cll_rev .exists() );
		return cll_rev .reverse ( tag::surely_exists );                  }  }


#ifndef NDEBUG
inline void Mesh::print_everything ( )
{	if ( not is_positive() ) std::cout << "(negative Mesh) ";
	core->print_everything ();                                }
#endif

//-----------------------------------------------------------------------------//


inline bool Cell::is_positive ( ) const
{	assert ( this->exists() );
	return this->core->is_positive ( );  }

inline Cell Cell::get_positive ( )
{	assert ( this->exists() );
	return this->core->get_positive();  }


inline size_t Cell::dim ( ) const
{	assert ( this->exists() );
	return this->core->get_dim ( );  }

inline size_t Mesh::dim ( ) const
{	return tag::Util::assert_diff ( this->core->get_dim_plus_one(), 1 );      }
// tag::Util::assert_diff  provides a safe way to substract two size_t numbers


inline Cell Cell::reverse ( const tag::BuildIfNotExists & build ) const
// 'build' defaults to tag::build_if_not_exists, so method may be called with no arguments
{	assert ( this->exists() );
	if ( not this->core->reverse_attr .exists() )
	{	assert ( this->is_positive() );
		Cell::Positive * pos_this_core = tag::Util::assert_cast
			< Cell::Core*, Cell::Positive* > ( this->core );
		this->core->reverse_attr = Cell ( tag::whose_core_is,
			pos_this_core->build_reverse ( tag::one_dummy_wrapper ), tag::freshly_created );  }
	return this->core->reverse_attr;                                                        }


inline Cell Cell::reverse ( const tag::MayNotExist & ) const
{	assert ( this->exists() );
	return this->core->reverse_attr;  }


inline Cell Cell::reverse ( const tag::SurelyExists & ) const
{	assert ( this->exists() );
	assert ( this->core->reverse_attr .exists() );
	return this->core->reverse_attr;              }


inline Mesh Mesh::reverse ( ) const
{	if ( this->is_positive() )
		return Mesh ( tag::whose_core_is, core, tag::previously_existing,
	                tag::is_negative, tag::build_cells_if_necessary     );
	// else
	return Mesh ( tag::whose_core_is, core, tag::previously_existing, tag::is_positive );  }
	

inline bool Cell::has_reverse ( ) const
{	assert ( this->exists() );
	return this->core->reverse_attr .exists();  }


inline Mesh Cell::boundary ( ) const

{	assert ( this->exists() );
	return this->core->boundary(); }

inline Cell Cell::tip () const
{	assert ( this->exists() );
	assert ( this->core->tip() .exists() );
	return this->core->tip();               }

inline Cell Cell::base () const
{	assert ( this->exists() );
	assert ( this->core->base() .exists() );
	return this->core->base();               }


#ifndef NDEBUG
inline void Cell::print_everything ( )
{	core->print_everything ( );  }
#endif


inline void Cell::glue_on_bdry_of ( Cell & cll )

// glue 'this' face on the boundary of cell 'cll'
// any of them may be negative

{	assert ( this->exists() );
	this->core->glue_on_bdry_of ( cll .core );  }


inline void Cell::glue_on_bdry_of ( Cell & cll, const tag::DoNotBother & )

// glue 'this' face on the boundary of cell 'cll'
// any of them may be negative

{	assert ( this->exists() );
	this->core->glue_on_bdry_of ( cll .core, tag::do_not_bother );  }


inline void Cell::cut_from_bdry_of ( Cell & cll )

// cut 'this' face from the boundary of cell 'cll'
// any of them may be negative

{	assert ( this->exists() );
	this->core->cut_from_bdry_of ( cll .core );  }


inline void Cell::cut_from_bdry_of ( Cell & cll, const tag::DoNotBother & )

// cut 'this' face from the boundary of cell 'cll'
// any of them may be negative

{	assert ( this->exists() );
	this->core->cut_from_bdry_of ( cll .core, tag::do_not_bother );  }


inline void Cell::add_to_mesh ( Mesh & msh )

// add 'this' cell to the mesh 'msh' by calling the virtual method add_to_mesh
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

{	assert ( this->exists() );
	assert ( msh .exists() );
	assert ( this->dim() == msh .dim() );
	assert ( this->dim() > 0 );
	if ( msh .is_positive() )  this->core->add_to_mesh ( msh .core );
	else
	{	assert( this->core->reverse_attr .exists() );
		this->core->reverse_attr .core->add_to_mesh ( msh .core );  }   }
// for negative Meshes, the core points towards the reverse, positive, Mesh::Core


inline void Cell::add_to_mesh ( Mesh & msh, const tag::DoNotBother & )

// add 'this' cell to the mesh 'msh' by calling the virtual method add_to_mesh
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

// tag::do_not_bother means : for Mesh::Connected::***,
// do not check for connectivity/openness/closedeness
// do not try to update first_vertex, last_vertex

{	assert ( this->exists() );
	assert ( msh .exists() );
	assert ( this->dim() == msh .dim() );
	assert ( this->dim() > 0 );
	if ( msh .is_positive() )  this->core->add_to_mesh ( msh .core, tag::do_not_bother );
	else
	{	assert( this->core->reverse_attr .exists() );
		this->core->reverse_attr .core->add_to_mesh ( msh .core, tag::do_not_bother );  }   }
// for negative Meshes, the core points towards the reverse, positive, Mesh::Core


inline void Cell::remove_from_mesh ( Mesh & msh )

// remove 'this' cell from the mesh 'msh' by calling the virtual method remove_from_mesh
// if 'msh' is the boundary of a cell, 'cut_from_bdry_of' should be used instead

{	assert ( this->exists() );
	assert ( this->dim() == msh .dim() );
	assert ( this->dim() > 0 );
	if ( msh .is_positive() )  this->core->remove_from_mesh ( msh .core );
	else
	{	assert( this->core->reverse_attr .exists() );
		this->core->reverse_attr .core->remove_from_mesh ( msh .core );  }   }
// for negative Meshes, the core points towards the reverse, positive, Mesh::Core


inline void Cell::remove_from_mesh ( Mesh & msh, const tag::DoNotBother & )

// remove 'this' cell from the mesh 'msh' by calling the virtual method remove_from_mesh
// if 'msh' is the boundary of a cell, 'cut_from_bdry_of' should be used instead

// tag::do_not_bother means : for Mesh::Connected::***,
// do not check for connectivity/openness/closedeness
// do not try to update first_vertex, last_vertex

{	assert ( this->exists() );
	assert ( this->dim() == msh .dim() );
	assert ( this->dim() > 0 );
	if ( msh .is_positive() )  this->core->remove_from_mesh ( msh .core, tag::do_not_bother );
	else
	{	assert ( this->core->reverse_attr .exists() );
		this->core->reverse_attr .core->remove_from_mesh ( msh .core, tag::do_not_bother );  }   }
// for negative Meshes, the core points towards the reverse, positive, Mesh::Core

//-----------------------------------------------------------------------------//


inline void Mesh::closed_loop ( const Cell & ver )
// for connected one-dim meshes, set both first_ver and last_ver to 'ver'
{	assert ( this->exists() );
	this->core->closed_loop ( ver );  }


inline void Mesh::closed_loop ( const Cell & ver, size_t n )
// for connected one-dim meshes, set both first_ver and last_ver to 'ver'
// and number of segments
{	assert ( this->exists() );
	this->core->closed_loop ( ver, n );  }

//-----------------------------------------------------------------------------//


inline Cell::PositiveVertex::PositiveVertex ( const tag::OneDummyWrapper & )
: Cell::Positive ( tag::of_dim, 0, tag::size_meshes, Mesh::maximum_dimension_plus_one,
                   tag::one_dummy_wrapper                                              )
{	}


inline Cell::NegativeVertex::NegativeVertex
( const tag::ReverseOf &, Cell::Positive * direct_ver_p, const tag::OneDummyWrapper & )
: Cell::Negative ( tag::of_dim, 0, tag::reverse_of, direct_ver_p, tag::one_dummy_wrapper )
{	assert ( direct_ver_p );
	assert ( direct_ver_p->get_dim() == 0 );    }


inline Cell::PositiveSegment::PositiveSegment ( Cell A, Cell B, const tag::OneDummyWrapper & )

:	Cell::Positive::NotVertex ( tag::of_dim, 1, tag::size_meshes,
	                            tag::Util::assert_diff ( Mesh::maximum_dimension_plus_one, 1 ),
	                            tag::one_dummy_wrapper                                         ),
	// tag::Util::assert_diff provides a safe way to substract two 'size_t' numbers
	base_attr { A }, tip_attr { B }

{	// below is a much simplified version of Cell::Negative::Vertex::add_to_seg
	// that's because 'this' segment has just been created, so it has no meshes above
	// also, the base and tip have already been correctly initialized
	assert ( not A .is_positive() );
	assert ( B .is_positive() );
	Cell pos_A =  A .reverse ( tag::surely_exists );
	assert ( pos_A .is_positive() );
	Cell::Positive::Vertex * pos_Aa = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( pos_A .core );
	assert ( pos_Aa->segments .find (this) == pos_Aa->segments.end() );
	// pos_Aa->meshes [0] [msh] = Cell::field_to_meshes { 0, 1 };
	// the third component 'where' is irrelevant here
	pos_Aa->segments.emplace ( std::piecewise_construct,
	      std::forward_as_tuple (this), std::forward_as_tuple (-1) );
	// below is a much simplified version of Cell::Positive::Vertex::add_to_seg
	// that's because 'this' segment has just been created, so it has no meshes above
	// also, the tip has already been correctly initialized
	Cell::Positive::Vertex * Bb = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( B .core );
	assert ( Bb->segments .find (this) == Bb->segments .end() );
	// Bb->meshes [0] [msh] = Cell::field_to_meshes { 1, 0 };
	// the third component 'where' is irrelevant here
	Bb->segments.emplace ( std::piecewise_construct,
	   std::forward_as_tuple (this), std::forward_as_tuple (1) );           }


inline Cell::NegativeNotVertex::NegativeNotVertex
( const tag::OfDimension &, size_t d,
  const tag::ReverseOf &, Cell::Positive * direct_cll_p, const tag::OneDummyWrapper & )

:	Cell::Negative ( tag::of_dim, d, tag::reverse_of, direct_cll_p, tag::one_dummy_wrapper )

{	}


inline Cell::NegativeSegment::NegativeSegment
( const tag::ReverseOf &, Cell::Positive * direct_seg_p, const tag::OneDummyWrapper & )

:	Cell::NegativeNotVertex ( tag::of_dim, 1, tag::reverse_of, direct_seg_p,
	                          tag::one_dummy_wrapper                        )

// we must make sure that both extremities of 'direct_seg_p' have a reverse
// well, the base surely has one since it's a NegativeVertex

{	assert ( direct_seg_p );
	Cell Ar = direct_seg_p->base(), B = direct_seg_p->tip();
	assert ( direct_seg_p->base() .exists() );
	assert ( direct_seg_p->tip() .exists() );
	assert ( direct_seg_p->base() .reverse ( tag::may_not_exist ) .exists() );
	Cell rev_tip = direct_seg_p->tip() .reverse ( tag::build_if_not_exists );  }	


inline Cell::PositiveHighDim::PositiveHighDim
( const tag::OfDimension &, const size_t d,
  const tag::WhoseBoundaryIs &, const Mesh & msh, const tag::OneDummyWrapper &    )

:	Cell::Positive::NotVertex ( tag::of_dim, d, tag::size_meshes,
	         tag::Util::assert_diff ( Mesh::maximum_dimension_plus_one, d ),
	                           tag::one_dummy_wrapper                       ),
	// tag::Util::assert_diff provides a safe way to substract two 'size_t' numbers
	boundary_attr ( msh )
	
{	assert ( msh .core );
	assert ( msh .core->get_dim_plus_one() == d );
	msh.core->cell_enclosed = this;                }


inline Cell::PositiveHighDim::PositiveHighDim
( const tag::WhoseBoundaryIs &, const Mesh & msh, const tag::OneDummyWrapper & )
:	Cell::Positive::HighDim ( tag::of_dim, msh.core->get_dim_plus_one(),
	                          tag::whose_bdry_is, msh, tag::one_dummy_wrapper )
{	}


inline Cell::NegativeHighDim::NegativeHighDim
( const tag::ReverseOf &, Cell::Positive * direct_cell_p, const tag::OneDummyWrapper & )
:	Cell::Negative::HighDim
	( tag::of_dim, direct_cell_p->get_dim(), tag::reverse_of, direct_cell_p,
	  tag::one_dummy_wrapper                                                )
{	}


inline Cell::NegativeHighDim::NegativeHighDim
(	const tag::OfDimension &, const size_t d, const tag::ReverseOf &,
	Cell::Positive * direct_cell_p, const tag::OneDummyWrapper &                    )
	
:	Cell::NegativeNotVertex
	( tag::of_dim, d, tag::reverse_of, direct_cell_p, tag::one_dummy_wrapper )

// we must make sure that all faces of 'direct_cell_p' have a reverse
// we build the reverse faces if necessary

{	assert ( direct_cell_p );
	assert ( direct_cell_p->get_dim() == d );
	assert ( direct_cell_p->boundary() .is_positive() );
	Cell direct_cell ( tag::whose_core_is, direct_cell_p,
                     tag::previously_existing, tag::surely_not_null );
	// Mesh::Iterator it ( tag::whose_core_is, direct_cell_p->iterator ( ...  !!
	Mesh::Iterator it = direct_cell_p->boundary() .iterator
		( tag::over_cells_of_max_dim, tag::as_they_are );
	for ( it .reset(); it .in_range(); it++ )
	{	Cell face = *it;
		assert ( face .exists() );
		Cell face_rev = face .reverse ( tag::build_if_not_exists );  }        }


inline Cell::PositiveHighDim::PositiveHighDim
( const tag::Triangle &, const Cell & AB, const Cell & BC, const Cell & CA,
  const tag::OneDummyWrapper &                                              )
	
:	Cell::Positive::HighDim
	( tag::whose_boundary_is,
	  Mesh ( tag::whose_core_is,
	         new Mesh::Connected::OneDim ( tag::with, 3, tag::segments,
//         new Mesh::Fuzzy ( tag::of_dimension, 2, tag::minus_one,
	                                       tag::one_dummy_wrapper       ),
	         tag::freshly_created                                         ),
	  tag::one_dummy_wrapper                                                )
	
{	assert ( AB .exists() );
	assert ( BC .exists() );
	assert ( CA .exists() );
	assert ( AB .dim() == 1 );
	assert ( BC .dim() == 1 );
	assert ( CA .dim() == 1 );
	assert ( AB .tip() == BC.base() .reverse() );
	assert ( BC .tip() == CA.base() .reverse() );
	assert ( CA .tip() == AB.base() .reverse() );
	// no need for glue_on_bdry_of : 'this' cell has just been created, it has no meshes above
	// we call add_to_mesh instead (as if the mesh were not a boundary)
	assert  ( this->boundary() .is_positive() );
	AB .core->add_to_mesh ( this->boundary() .core, tag::do_not_bother );
	BC .core->add_to_mesh ( this->boundary() .core, tag::do_not_bother );
	CA .core->add_to_mesh ( this->boundary() .core, tag::do_not_bother );
	Mesh::Connected::OneDim * bdry = tag::Util::assert_cast
		< Mesh::Core*, Mesh::Connected::OneDim* > ( this->boundary() .core );
	bdry->first_ver = AB .base();  // negative
	bdry->last_ver = CA .tip();                                              }
	

inline Cell::PositiveHighDim::PositiveHighDim
( const tag::Quadrangle &, const Cell & AB, const Cell & BC, const Cell & CD,
  const Cell & DA, const tag::OneDummyWrapper &                               )
	
:	Cell::Positive::HighDim
	( tag::whose_boundary_is,
	  Mesh ( tag::whose_core_is,
	         new Mesh::Connected::OneDim ( tag::with, 4, tag::segments,
//         new Mesh::Fuzzy ( tag::of_dimension, 2, tag::minus_one,
	                                       tag::one_dummy_wrapper       ),
	         tag::freshly_created                                         ),
	  tag::one_dummy_wrapper                                                )

{	assert ( AB .exists() );
	assert ( BC .exists() );
	assert ( CD .exists() );
	assert ( DA .exists() );
	assert ( AB .dim() == 1 );
	assert ( BC .dim() == 1 );
	assert ( CD .dim() == 1 );
	assert ( DA .dim() == 1 );
	assert ( AB .tip() == BC .base() .reverse() );
	assert ( BC .tip() == CD .base() .reverse() );
	assert ( CD .tip() == DA .base() .reverse() );
	assert ( DA .tip() == AB .base() .reverse() );
	// no need for glue_on_bdry_of : 'this' cell has just been created, it has no meshes above
	// we call add_to_mesh instead (as if the mesh were not a boundary)
	assert  ( this->boundary() .is_positive() );
	AB .core->add_to_mesh ( this->boundary() .core, tag::do_not_bother );
	BC .core->add_to_mesh ( this->boundary() .core, tag::do_not_bother );
	CD .core->add_to_mesh ( this->boundary() .core, tag::do_not_bother );
	DA .core->add_to_mesh ( this->boundary() .core, tag::do_not_bother );
	Mesh::Connected::OneDim * bdry = tag::Util::assert_cast
		< Mesh::Core*, Mesh::Connected::OneDim* > ( this->boundary().core );
	bdry->first_ver = AB .base();  // negative
	bdry->last_ver = DA .tip();                                            }


inline Cell::PositiveHighDim::PositiveHighDim
( const tag::Pentagon &, const Cell & AB, const Cell & BC, const Cell & CD,
  const Cell & DE, const Cell & EA, const tag::OneDummyWrapper &            )
	
:	Cell::Positive::HighDim
	( tag::whose_boundary_is,
	  Mesh ( tag::whose_core_is,
	         new Mesh::Connected::OneDim ( tag::with, 5, tag::segments,
	                                       tag::one_dummy_wrapper       ),
	         tag::freshly_created                                         ),
	  tag::one_dummy_wrapper                                                )

{	assert ( AB .exists() );
	assert ( BC .exists() );
	assert ( CD .exists() );
	assert ( DE .exists() );
	assert ( EA .exists() );
	assert ( AB .dim() == 1 );
	assert ( BC .dim() == 1 );
	assert ( CD .dim() == 1 );
	assert ( DE .dim() == 1 );
	assert ( EA .dim() == 1 );
	assert ( AB .tip() == BC .base() .reverse() );
	assert ( BC .tip() == CD .base() .reverse() );
	assert ( CD .tip() == DE .base() .reverse() );
	assert ( DE .tip() == EA .base() .reverse() );
	assert ( EA .tip() == AB .base() .reverse() );
	// no need for glue_on_bdry_of : 'this' cell has just been created, it has no meshes above
	// we call add_to_mesh instead (as if the mesh were not a boundary)
	assert  ( this->boundary() .is_positive() );
	AB .core->add_to_mesh ( this->boundary() .core, tag::do_not_bother );
	BC .core->add_to_mesh ( this->boundary() .core, tag::do_not_bother );
	CD .core->add_to_mesh ( this->boundary() .core, tag::do_not_bother );
	DE .core->add_to_mesh ( this->boundary() .core, tag::do_not_bother );
	EA .core->add_to_mesh ( this->boundary() .core, tag::do_not_bother );
	Mesh::Connected::OneDim * bdry = tag::Util::assert_cast
		< Mesh::Core*, Mesh::Connected::OneDim* > ( this->boundary() .core );
	bdry->first_ver = AB .base();  // negative
	bdry->last_ver = EA .tip();                                             }


inline Cell::PositiveHighDim::PositiveHighDim
( const tag::Hexagon &, const Cell & AB, const Cell & BC, const Cell & CD,
  const Cell & DE, const Cell & EF, const Cell & FA, const tag::OneDummyWrapper & )
	
:	Cell::Positive::HighDim
	( tag::whose_boundary_is,
	  Mesh ( tag::whose_core_is,
	         new Mesh::Connected::OneDim ( tag::with, 6, tag::segments,
	                                       tag::one_dummy_wrapper       ),
	         tag::freshly_created                                         ),
	  tag::one_dummy_wrapper                                                )

{	assert ( AB .exists() );
	assert ( BC .exists() );
	assert ( CD .exists() );
	assert ( DE .exists() );
	assert ( EF .exists() );
	assert ( FA .exists() );
	assert ( AB .dim() == 1 );
	assert ( BC .dim() == 1 );
	assert ( CD .dim() == 1 );
	assert ( DE .dim() == 1 );
	assert ( EF .dim() == 1 );
	assert ( FA .dim() == 1 );
	assert ( AB .tip() == BC .base() .reverse() );
	assert ( BC .tip() == CD .base() .reverse() );
	assert ( CD .tip() == DE .base() .reverse() );
	assert ( DE .tip() == EF .base() .reverse() );
	assert ( EF .tip() == FA .base() .reverse() );
	assert ( FA .tip() == AB .base() .reverse() );
	// no need for glue_on_bdry_of : 'this' cell has just been created, it has no meshes above
	// we call add_to_mesh instead (as if the mesh were not a boundary)
	assert  ( this->boundary() .is_positive() );
	AB .core->add_to_mesh ( this->boundary() .core, tag::do_not_bother );
	BC .core->add_to_mesh ( this->boundary() .core, tag::do_not_bother );
	CD .core->add_to_mesh ( this->boundary() .core, tag::do_not_bother );
	DE .core->add_to_mesh ( this->boundary() .core, tag::do_not_bother );
	EF .core->add_to_mesh ( this->boundary() .core, tag::do_not_bother );
	FA .core->add_to_mesh ( this->boundary() .core, tag::do_not_bother );
	Mesh::Connected::OneDim * bdry = tag::Util::assert_cast
		< Mesh::Core*, Mesh::Connected::OneDim* > ( this->boundary() .core );
	bdry->first_ver = AB .base();  // negative
	bdry->last_ver = FA .tip();                                             }

//---------------------------------------------------------------------------------


inline bool tag::Util::Core::dispose_query ( )
{	assert ( this->nb_of_wrappers > 0 );
	this->nb_of_wrappers--;
	return this->nb_of_wrappers == 0;     }


//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------

class Cell::Numbering

// abstract class, specialized in Cell::Numbering::Field and Cell::Numbering::Map

{	public :

	virtual ~Numbering ( ) { };

	virtual size_t size ( ) = 0;

	virtual size_t & operator[] ( const Cell ) = 0;

	class Map;  class Field;  // latter defined in field.h

};	


class Cell::Numbering::Map : public Cell::Numbering

{	public :

	std::map < Cell, size_t > * map;

	inline Map ( )
	: map { new std::map < Cell, size_t > }
	{	}

	inline Map ( std::map < Cell, size_t > * numb_map )
	:	map { numb_map }
	{	}

	size_t size ( )  // virtual from Cell::Numbering
	{	return this->map->size();  }
	
	size_t & operator[] ( const Cell );  // virtual from Cell::Numbering, defined in global.cpp

};	


}  // namespace maniFEM

#endif  // ifndef MANIFEM_MESH_H

