
// maniFEM mesh.h 2020.01.26

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

#ifndef MANIFEM_MESH_H
#define MANIFEM_MESH_H

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include "assert.h"

namespace maniFEM {

namespace tag {
	// see paragraph 8.2 in the manual
	struct IsNegative { };  static const IsNegative is_negative;
	struct IsPositive { };  static const IsPositive is_positive;
	struct Reverse { };  static const Reverse reverse;
	struct ReverseOf { };  static const ReverseOf reverse_of;
	struct NonExistent { };  static const NonExistent non_existent;
	struct BuildIfNotExists { };  static const BuildIfNotExists build_if_not_exists;
	struct MayNotExist { };  static const MayNotExist may_not_exist;
	struct SurelyExists { };  static const SurelyExists surely_exists;
	struct OnTheFly { };  static const OnTheFly on_the_fly;
	struct OfDimension { };  static const OfDimension of_dim;
	                         static const OfDimension of_dimension;
	struct OfDimensionOne { };  static const OfDimensionOne of_dim_one;
	                            static const OfDimensionOne of_dimension_one;
	struct MinusOne { };  static const MinusOne minus_one;
	struct GreaterThanOne { };  static const GreaterThanOne greater_than_one;
	struct MightBeOne { };  static const MightBeOne might_be_one;
	struct Oriented { };  static const Oriented oriented;
	struct NotOriented { };  static const NotOriented not_oriented;
	struct BuildCellsIfNecessary { };  static const BuildCellsIfNecessary build_cells_if_necessary;
	struct Progressive { };  static const Progressive progressive;
	struct StartAt { };  static const StartAt start_at;
	struct StopAt { };  static const StopAt stop_at;
	struct Towards { };  static const Towards towards;
	struct Boundary { };  static const Boundary boundary;
	struct Ghost { };  static const Ghost ghost;
	struct Bizarre { };  static const Bizarre bizarre;
	struct SizeMeshes { };  static const SizeMeshes size_meshes;
	struct BehindFace { };  static const BehindFace behind_face;
	struct InFrontOfFace { };  static const InFrontOfFace in_front_of_face;
	struct WithinMesh { };  static const WithinMesh within_mesh;
	struct CellsOfDim { };  static const CellsOfDim cells_of_dim;
	struct CellsOfReverseOf { }; static const CellsOfReverseOf cells_of_reverse_of;
	struct Vertices { };  static const Vertices vertices;
	struct Segments { };  static const Segments segments;
	struct Vertex { };  static const Vertex vertex; static const Vertex point;
	struct Segment { };  static const Segment segment;
	struct DividedIn { };  static const DividedIn divided_in;
	struct Triangle { };  static const Triangle triangle;
	struct Quadrangle { };  static const Quadrangle quadrangle;
	                        static const Quadrangle rectangle;
	                        static const Quadrangle square;
	                        static const Quadrangle quadrilateral;
	struct WhoseBoundaryIs { };  static const WhoseBoundaryIs whose_bdry_is;
	                             static const WhoseBoundaryIs whose_boundary_is;
	struct WhoseCoreIs { };  static const WhoseCoreIs whose_core_is;
	struct OverCellsOf { };  static const OverCellsOf over_cells_of;
	struct ForcePositive { };  static const ForcePositive force_positive;
	struct HasSize { };  static const HasSize has_size;
	struct ReserveSize { };  static const ReserveSize reserve_size;
	struct lagrange { };  static const lagrange Lagrange;
	struct Pretty { };  static const Pretty pretty;
	struct OfDegree { };  static const OfDegree of_degree;
	enum WithTriangles { with_triangles, not_with_triangles };
	struct Join { };  static const Join join;
	struct Onto { };  static const Onto onto;
	struct Around { };  static const Around around;
	struct Iff { };  static const Iff iff;
	struct LessThan { };  static const LessThan less_than;
	struct IfLessThan { };  static const IfLessThan if_less_than;
	struct Otherwise { };  static const Otherwise otherwise;
	struct EntireManifold { };  static const EntireManifold entire_manifold;
	struct DesiredLength { };  static const DesiredLength desired_length;
	struct IntrinsicOrientation { };  static const IntrinsicOrientation intrinsic_orientation;
	struct InherentOrientation { };  static const InherentOrientation inherent_orientation;
	struct RandomOrientation { };  static const RandomOrientation random_orientation;   }

class Cell;  class Mesh;
class CellIterator;  class MeshIterator;
class Manifold;

//-----------------------------------------------------------------------------//


// a cell of dimension zero is a point, see class Cell::Positive::Vertex and Negative::Vertex
// a cell of dimension one is a segment, see class Cell::Positive::Segment and Negative::Segment
// for cells of dimension two or more, see class Cell::Positive and Cell::Negative
// a cell of dimension two may be a triangle, a quadrangle or some other polygon
// a cell of dimension three may be a tetrahedron, a cube or some other polyhedron
// cells of dimension four or higher may be constructed,
// but their usefulness is questionable

// cells may be positively or negatively oriented
// see class Cell::Core::Positive and class Cell::Core::Negative
// which is which depends only on the construction process
// the first one to be created is positive, the other one will be negative,
// created when we call the 'reverse' method

// a cell is mainly defined by its boundary (which is a mesh of lower dimension)
// the orientation of a cell is nothing more than an orientation of its boundary
// see the comments on orientation in class Mesh below

class Cell

// a thin wrapper around a Cell::Core, with most methods delegated to 'core'

{	public :

	class Core;
	
	Cell::Core * core;  // use sharedptr ?

	// many constructors defined after class Mesh::Negative

	inline Cell ( const tag::NonExistent & )
	:	core { nullptr }
	{ }

	inline Cell ( const tag::WhoseCoreIs &, Cell::Core * c )
	:	core { c }
	{	assert ( c );  }

	inline Cell ( const tag::WhoseBoundaryIs &, Mesh & );
	inline Cell ( const tag::ReverseOf &, const Cell & direct_cell,
                const tag::BuildIfNotExists & build = tag::build_if_not_exists );
	inline Cell ( const tag::ReverseOf &, const Cell & direct_cell, const tag::SurelyExists & );
	inline Cell ( const tag::ReverseOf &, const Cell & direct_cell, const tag::MayNotExist & );
	// inline Cell ( const tag::WhoseBoundaryIs &, Mesh::Core * );
	inline Cell ( const tag::Vertex &, const tag::IsPositive & ispos = tag::is_positive );
	inline Cell ( const tag::Segment &, const Cell & A, const Cell & B );
	inline Cell ( const tag::Bizarre &, const Cell &, const Cell & );
	inline Cell ( const tag::Triangle &, const Cell & AB, const Cell & BC, const Cell & CA );
	inline Cell ( const tag::Quadrangle &, const Cell & AB, const Cell & BC,
                                         const Cell & CD, const Cell & DA );
	inline Cell ( const tag::BehindFace &, const Cell &, const tag::WithinMesh &,
	              const Mesh &, const tag::SurelyExists & se = tag::surely_exists );
	inline Cell ( const tag::BehindFace &, const Cell &,
                const tag::WithinMesh &, const Mesh &, const tag::MayNotExist & );
	inline Cell ( const tag::InFrontOfFace &, const Cell &, const tag::WithinMesh &,
                const Mesh &, const tag::SurelyExists & se = tag::surely_exists     );
	inline Cell ( const tag::InFrontOfFace &, const Cell &,
                const tag::WithinMesh &, const Mesh &, const tag::MayNotExist & );
	// use : Cell tri ( tag::in_front_of_face, f, tag::within_mesh, msh );

	// methods delegated to 'core', defined after class Mesh::Negative

	inline Cell reverse
	( const tag::BuildIfNotExists & build = tag::build_if_not_exists ) const;
	inline Cell reverse ( const tag::SurelyExists & ) const;
	inline Cell reverse ( const tag::MayNotExist & ) const;
	inline Mesh boundary () const;
	inline bool exists () const  { return core != nullptr;  }
	inline bool is_positive () const;
	inline Cell get_positive ();
	inline size_t dim () const;
	inline bool has_reverse () const;
	inline Cell tip () const;
	inline Cell base () const;

	inline bool belongs_to ( const Mesh & msh, const tag::Oriented & ) const;
	inline bool belongs_to ( const Mesh & msh, const tag::NotOriented & ) const;

	// method 'glue_on_bdry_of' is intensively used when building a mesh
	// it glues 'this' cell to the boundary of 'cll'
	inline void glue_on_bdry_of ( Cell & cll );
	
	// method 'cut_from_bdry_of' does the reverse : cuts 'this' cell from
	// the boundary of 'cll' - used mainly in remeshing
	inline void cut_from_bdry_of ( Cell & cll );
	
	// methods 'add_to' and 'remove_from' add/remove 'this' cell to/from the mesh 'msh'
	// if 'msh' is the boundary of some cell, methods 'glue_on_bdry_of'
	// and 'cut_from_bdry_of' (see above) should be used instead
	inline void add_to ( Mesh & msh );
	inline void remove_from ( Mesh & msh );

	inline void dispose ( ) { };  //  do intersting things !!
	
	inline void project ( ) const;
	
	inline void project ( const tag::Onto &, const Manifold m ) const;
		
#ifndef NDEBUG
	inline void print_everything ( );
#endif

	void print_coords();

	static std::vector < size_t > double_heap_size_pos, double_heap_size_neg,
		size_t_heap_size_pos, size_t_heap_size_neg, short_int_heap_size_pos, short_int_heap_size_neg;

	static Cell::Core * const ghost;
	// see paragraph 8.14 in the manual
	
	struct field_to_meshes
	{	short int counter_pos;
		short int counter_neg;
		std::list<Cell::Core*>::iterator where;
		inline field_to_meshes ( short int i, short int j )
		:	counter_pos {i}, counter_neg {j} { }
		inline field_to_meshes ( short int i, short int j,
		                         std::list<Cell::Core*>::iterator w )
		:	counter_pos {i}, counter_neg {j}, where {w} { }                    };

  class Positive;  class Negative;  class Vertex;  class Segment;
	
}; // end of  class Cell


inline bool operator== ( const Cell & c1, const Cell & c2 )
{	return c1.core == c2.core;  }

inline bool operator!= ( const Cell & c1, const Cell & c2 )
{	return c1.core != c2.core;  }

inline bool operator< ( const Cell & c1, const Cell & c2 )
{	return c1.core < c2.core;  }

//-----------------------------------------------------------------------------//


// roughly speaking, a mesh is a collection of cells of the same dimension
// however, for efficiency purposes, we keep lists of cells of lower
// dimension as well; that's the purpose of 'IndexedList <Cell*> cells' :
// to keep lists of cells indexed by their dimension
// this implies quite some amount of redundant information,
// but this redundancy makes the classes fast, especially for remeshing

// an orientation of the mesh
// is nothing more than an orientation of each of its cells (of maximum dimension).
// but these orientations cannot be arbitrary, they must be compatible
// in the sense that a face common to two cells must be seen as positive
// from one of the cells and as negative from the other cell

// as a consequence, for connected meshes there are only two possible orientations
// although nothing prevents a mesh to be disconnected, we only allow for two orientations
// which is which depends only on the construction process
// the first one to be created is positive, the other one will be negative,
// created when we call the 'reverse' method

// negative meshes will appear only as boundaries of negative cells
// that's why we do not store their core in the computer's memory
// wrappers for negative meshes are built on-the-fly, e.g. in method Cell::boundary
// their core points to the (positive) reverse, the only difference is in 'meth'

class Mesh

{	public :
	
	class Core;

	struct Methods
	{	bool (*is_positive) ( );
		Mesh (*reverse) ( Mesh::Core * );  };

	// use sharedptr ?
	Mesh::Core * core;  // there are no cores for negative meshes
	// instead, we keep here a pointer to the direct (positive) core
	// thus, for meshes, 'core' points always to a positive Mesh::Core

	// emulate virtual methods
	// this is how we distinguish between a positive mesh and a negative one
	const Methods * meth;

	// we keep here the topological dimension of the largest mesh we intend to build
	// see method 'set_max_dim' and paragraph 8.5 in the manual
	static size_t maximum_dimension_plus_one;

	// constructors :
	
	inline Mesh ( const tag::OfDimension &, size_t dim, const tag::GreaterThanOne &,
                const tag::IsPositive & ispos = tag::is_positive                   );
	inline Mesh ( const tag::OfDimension &, size_t dim, const tag::MightBeOne &,
                const tag::IsPositive & ispos = tag::is_positive               );
	inline Mesh ( const tag::OfDimensionOne &,
                const tag::IsPositive & ispos = tag::is_positive );
	inline Mesh ( const tag::WhoseCoreIs &, Mesh::Core *,
                const tag::IsPositive & ispos = tag::is_positive );
	// builds a negative mesh from a positive one, assuming all cells have reverse :
	inline Mesh ( const tag::WhoseCoreIs &, Mesh::Core *, const tag::IsNegative &,
	              const tag::SurelyExists &                                        );
	// builds a negative mesh from a positive one, creating reverse cells if necessary :
	inline Mesh ( const tag::WhoseCoreIs &, Mesh::Core *, const tag::IsNegative &,
                const tag::BuildCellsIfNecessary &                                  );

	// geometric constructors are defined in global.cpp
	// segment, triangle, quadrangle, join
	inline Mesh ( const tag::Segment &, const Cell & A, const Cell & B, const tag::DividedIn &, size_t n );
	inline Mesh ( const tag::Pretty &, const tag::Segment &,
                const Cell & A, const Cell & B, const tag::DividedIn &, size_t n );
	inline Mesh ( const tag::Triangle &, const Mesh & AB, const Mesh & BC, const Mesh & CA );
	inline Mesh ( const tag::Pretty &, const tag::Triangle &,
                const Mesh & AB, const Mesh & BC, const Mesh & CA );
	inline Mesh ( const tag::Quadrangle &, const Mesh & south, const Mesh & east,
	                                       const Mesh & north, const Mesh & west,
	              const tag::WithTriangles & wt = tag::not_with_triangles         );
	inline Mesh ( const tag::Quadrangle &, const Cell & SW, const Cell & SE,
	                                       const Cell & NE, const Cell & NW, size_t m, size_t n,
	              const tag::WithTriangles & wt = tag::not_with_triangles                        );
	inline Mesh ( const tag::Pretty &, const tag::Quadrangle &,
	              const Mesh & south, const Mesh & east, const Mesh & north, const Mesh & west,
	              const tag::WithTriangles & wt = tag::not_with_triangles );
	inline Mesh ( const tag::Join &, const std::list<Mesh> & l );
	inline Mesh ( const tag::Join &, const Mesh &, const Mesh & );
  inline Mesh ( const tag::Join &, const Mesh &, const Mesh &, const Mesh & );
  inline Mesh ( const tag::Join &, const Mesh &, const Mesh &, const Mesh &, const Mesh & );
  inline Mesh ( const tag::Join &, const Mesh &, const Mesh &, const Mesh &, const Mesh &, const Mesh & );

	// constructors with tag::Progressive are defined in progressive.cpp
	
	Mesh ( const tag::Progressive &, const tag::DesiredLength &, double desired_length );

	Mesh ( const tag::Progressive &, const tag::EntireManifold, Manifold manif,
         const tag::DesiredLength &, double desired_length                    );

	Mesh ( const tag::Progressive &, const tag::DesiredLength &, double desired_length,
	       const tag::RandomOrientation &                                               );

	Mesh ( const tag::Progressive &, const tag::EntireManifold, Manifold manif,
         const tag::DesiredLength &, double desired_length, const tag::RandomOrientation & );

	Mesh ( const tag::Progressive &, const tag::DesiredLength &, double desired_length,
	       const tag::InherentOrientation &                                             );

	Mesh ( const tag::Progressive &, const tag::EntireManifold, Manifold manif,
	       const tag::DesiredLength &, double desired_length, const tag::InherentOrientation & );

	Mesh ( const tag::Progressive &, const tag::Boundary &, Mesh interface,
	       const tag::DesiredLength &, double desired_length                );

	inline Mesh ( const tag::Progressive &, const tag::Boundary &, Mesh interface,
	              const tag::DesiredLength &, double desired_length, const tag::RandomOrientation & )
	:	Mesh ( tag::progressive, tag::boundary, interface, tag::desired_length, desired_length )  { }

	Mesh ( const tag::Progressive &, const tag::Boundary &, Mesh interface,
	       const tag::DesiredLength &, double desired_length, const tag::IntrinsicOrientation & );

	Mesh ( const tag::Progressive &, const tag::Boundary &, Mesh interface,
	       const tag::DesiredLength &, double desired_length, const tag::InherentOrientation & );

	Mesh ( const tag::Progressive &, const tag::Boundary &, Mesh interface,
	       const tag::StartAt &, const Cell & start,
	       const tag::Towards &, std::vector<double> & normal,
	       const tag::DesiredLength &, double desired_length                   );

	Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	       const tag::Towards &, std::vector<double> & tangent,
	       const tag::StopAt &, const Cell & stop,
	       const tag::DesiredLength &, double desired_length                   );

	Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	       const tag::Towards &, std::vector<double> & tangent,
	       const tag::DesiredLength &, double desired_length                   );

	Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	       const tag::StopAt &, const Cell & stop,
	       const tag::DesiredLength &, double desired_length                   );

	Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	       const tag::DesiredLength &, double desired_length                   );

	inline Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	              const tag::StopAt &, const Cell & stop,
	              const tag::DesiredLength &, double desired_length, const tag::RandomOrientation & )
	:	Mesh ( tag::progressive, tag::start_at, start, tag::stop_at, stop,
		       tag::desired_length, desired_length                          )  { }

	Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	       const tag::StopAt &, const Cell & stop,
	       const tag::DesiredLength &, double desired_length, const tag::InherentOrientation & );

	inline Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	              const tag::DesiredLength &, double desired_length, const tag::InherentOrientation & )
	:	Mesh ( tag::progressive, tag::start_at, start, tag::stop_at, start,
		       tag::desired_length, desired_length, tag::inherent_orientation )  { }

	void pretty_constructor ( const tag::Segment &, const Cell & A, const Cell & B,
                            const tag::DividedIn &, size_t n );
	void pretty_constructor ( const tag::Triangle &, const Mesh & AB, const Mesh & BC, const Mesh & CA );
	void pretty_constructor ( const tag::Quadrangle &, const Mesh & south, const Mesh & east,
	             const Mesh & north, const Mesh & west, const tag::WithTriangles & wt        );

	inline bool is_positive () const;
	inline size_t dim () const;
	inline Mesh reverse () const;

	inline size_t number_of ( const tag::CellsOfDim &, size_t d ) const;
	inline size_t number_of ( const tag::Vertices & ) const;
	inline size_t number_of ( const tag::Segments & ) const;

	inline Cell first_vertex ( ) const;
	inline Cell last_vertex ( ) const;
	inline Cell first_segment ( ) const;
	inline Cell last_segment ( ) const;

	inline Cell cell_in_front_of
	( const Cell face, const tag::SurelyExists & se = tag::surely_exists ) const;

	inline Cell cell_behind
	( const Cell face, const tag::SurelyExists & se = tag::surely_exists ) const;

	inline Cell cell_in_front_of ( const Cell face, const tag::MayNotExist & ) const;

	inline Cell cell_behind ( const Cell face, const tag::MayNotExist & ) const;

	inline Cell::Core * cell_in_front_of
	( const Cell::Core * face_p, const tag::SurelyExists & se = tag::surely_exists ) const;

	inline Cell::Core * cell_in_front_of
	( const Cell::Core * face_p, const tag::MayNotExist & ) const;

	inline Cell::Core * cell_behind
	( const Cell::Core * face_p, const tag::SurelyExists & se = tag::surely_exists ) const;

	inline Cell::Core * cell_behind
	( const Cell::Core * face_p, const tag::MayNotExist & ) const;
	
	inline void dispose ( ) { };  //  do intersting things !!
	
	void join_list ( const std::list<Mesh> & l );

	void inline baricenter ( const Cell & ver, const Cell & seg );
	// defined in manifold.h

	inline CellIterator iter_over ( const tag::CellsOfDim &, size_t d ) const;
	inline CellIterator iter_over
	( const tag::CellsOfDim &, size_t d, const tag::ForcePositive & ) const;
	inline CellIterator iter_over
	( const tag::CellsOfDim &, size_t d, const tag::Reverse & ) const;
	inline CellIterator iter_over
	( const tag::CellsOfDim &, size_t d, const tag::ForcePositive &, const tag::Reverse & ) const;
	inline CellIterator iter_over
	( const tag::CellsOfDim &, size_t d, const tag::Reverse &, const tag::ForcePositive & ) const;
	inline CellIterator iter_over ( const tag::Vertices & ) const;
	inline CellIterator iter_over ( const tag::Vertices &, const tag::ForcePositive & ) const;
	inline CellIterator iter_over ( const tag::Vertices &, const tag::Reverse & ) const;
	inline CellIterator iter_over
	( const tag::Vertices &, const tag::ForcePositive &, const tag::Reverse & ) const;
	inline CellIterator iter_over
	( const tag::Vertices &, const tag::Reverse &, const tag::ForcePositive & ) const;
	inline CellIterator iter_over ( const tag::Segments & ) const;
	inline CellIterator iter_over ( const tag::Segments &, const tag::ForcePositive & ) const;
	inline CellIterator iter_over ( const tag::Segments &, const tag::Reverse & ) const;
	inline CellIterator iter_over
	( const tag::Segments &, const tag::ForcePositive &, const tag::Reverse & ) const;
	inline CellIterator iter_over
	( const tag::Segments &, const tag::Reverse &, const tag::ForcePositive & ) const;

	// methods draw_ps and export_msh defined in global.cpp

	void draw_ps ( std::string file_name );
	void draw_ps_3d ( std::string file_name );
	void export_msh ( std::string f, std::map<Cell::Core*,size_t> & ver_numbering );
	void export_msh ( std::string f );
	
#ifndef NDEBUG
	inline void print_everything ( );
#endif
	
	inline static size_t diff ( size_t a, size_t b )
	{	assert ( a >= b );  return  a - b;  }

	inline static void set_max_dim ( size_t d )
	// see paragraph 8.5 in the manual
	{	maximum_dimension_plus_one = d + 1;
		Cell::double_heap_size_pos.resize ( maximum_dimension_plus_one, 0 );
		Cell::double_heap_size_neg.resize ( maximum_dimension_plus_one, 0 );
		Cell::size_t_heap_size_pos.resize ( maximum_dimension_plus_one, 0 );
		Cell::size_t_heap_size_neg.resize ( maximum_dimension_plus_one, 0 );
		Cell::short_int_heap_size_pos.resize ( maximum_dimension_plus_one, 0 );
		Cell::short_int_heap_size_neg.resize ( maximum_dimension_plus_one, 0 );  }
	
	static void action_add
		( Cell::Core *, Cell::Core *, Mesh::Core *, short int, short int );
	static void action_remove
		( Cell::Core *, Cell::Core *, Mesh::Core *, short int, short int );
	static void action_add_rev
		( Cell::Core *, Cell::Core *, Mesh::Core *, short int, short int );
	static void action_remove_rev
		( Cell::Core *, Cell::Core *, Mesh::Core *, short int, short int );
	// methods action_* are passed to 'deep_connections' from  Cell::add_to and Cell::remove_from

	class Positive;  class Negative;
	struct OneDim  {  class Positive;  };
	
}; // end of  class Mesh


inline bool operator== ( const Mesh & m1, const Mesh & m2 )
{	return m1.core == m2.core;  }

//-----------------------------------------------------------------------------//
//-----------------------------------------------------------------------------//


class Cell::Core

{	public :

	Cell::Core * reverse_p { nullptr };

	// we keep numeric values here :
	std::vector < double > double_heap;
	std::vector < size_t > size_t_heap;
	std::vector < short int > short_int_heap;

	// if 'this' is a face of another cell and that other cell belongs to some mesh msh,
	// cell_behind_within[msh] keeps that cell
	// see methods Mesh::cell_behind and Mesh::cell_in_front_of
	std::map < Mesh::Core*, Cell::Core* > cell_behind_within;

	inline Core ( const tag::IsPositive &, const tag::OfDimension &, size_t d )
	:	double_heap ( Cell::double_heap_size_pos [d] ),
		size_t_heap ( Cell::size_t_heap_size_pos [d] ),
		short_int_heap ( Cell::short_int_heap_size_pos [d] )
	{ }

	inline Core ( const tag::IsNegative& ,
		const tag::ReverseOf &, Cell::Core * rev, const tag::OfDimension &, size_t d )
	:	reverse_p { rev },
		double_heap ( Cell::double_heap_size_neg [d] ),
		size_t_heap ( Cell::size_t_heap_size_neg [d] ),
		short_int_heap ( Cell::short_int_heap_size_neg [d] )
	{	assert ( rev );  }

	inline Core ( const tag::Ghost & )  { }

	virtual ~Core ( ) { };

	class Positive;  class Negative;

	virtual bool is_positive ( ) const = 0;
	virtual Cell::Core::Positive * get_positive ( ) = 0;
	virtual size_t get_dim ( ) const = 0;
	virtual Cell::Core * reverse ( const tag::BuildIfNotExists & ) = 0;
	virtual Cell::Core * tip ();
	virtual Cell::Core * base ();

	virtual bool belongs_to ( Mesh::Core *, const tag::Oriented & ) const = 0;
	virtual bool belongs_to ( Mesh::Core *, const tag::NotOriented & ) const = 0;

	// Method 'glue_on_bdry_of' is intensively used when building a mesh,
	// e.g. within factory functions in Cell class.
	// It glues 'this' cell to the boundary of 'cll'.
	inline void glue_on_bdry_of ( Cell::Core * cll );

	// Method 'cut_from_bdry_of' does the reverse : cuts 'this' cell from
	// the boundary of 'cll'. Used mainly in remeshing.
	inline void cut_from_bdry_of ( Cell::Core * cll );
	
	// Methods 'add_to' and 'remove_from' add/remove 'this' cell to/from the mesh 'msh'.
	// If 'msh' is the boundary of some cell, methods 'glue_on_bdry_of'
	// and 'cut_from_bdry_of' (see above) should be used instead.
	virtual void add_to ( Mesh::Core * msh ) = 0;
	virtual void remove_from ( Mesh::Core * msh ) = 0;

	virtual void glue_on_my_bdry ( Cell::Core * ) = 0;
	virtual void cut_from_my_bdry ( Cell::Core * ) = 0;
	// we feel that 'glue_on_bdry_of' and 'cut_from_bdry_of' are more readable
	// so we suggest to use those (see above)

	virtual void forget () { };
	// assert there are no meshes above
	// for each face, cut_from_bdry
	// if the face has no other meshes above, forget it
	// forget the (now empty) boundary, then delete this cell

#ifndef NDEBUG
	std::string name;
	virtual std::string get_name ( ) = 0;
	virtual void print_everything ( ) = 0;
#endif

}; // end of class Cell::Core

//-----------------------------------------------------------------------------//

class Cell::Core::Positive : public Cell::Core

{	public :

	// the 'meshes' attribute keeps information about all meshes
	// "above" 'this' cell, that is, all meshes containing 'this' as a cell
	// it is indexed over the dimension of the mesh minus the dimension of 'this' 
	// for each mesh, it keeps a 'Cell::field_to_meshes' value, containing two counters
	// and an iterator into the 'cells' field of that particular mesh
	// of course this implies quite some amount of redundant information
	// but this redundancy makes the classes fast, especially for remeshing
	
	// since indices of vectors begin at zero
	// the keys of 'meshes[i]' will be meshes of dimension 'i + this->dim'
	// on the other hand, there are no zero-dimensional meshes,
	// so if this->dim == 0 (when 'this' is a vertex) the keys of meshes[0] are pointers
	// to Positive::Segment (they must be explicitly converted to Positive::Segment*)

	std::vector < std::map < Mesh::Core *, Cell::field_to_meshes > > meshes;

#ifndef NDEBUG
	std::string get_name();
#endif
	
	inline Positive ( const tag::OfDimension &, size_t d,
                            const tag::SizeMeshes &, size_t sz  )
	:	Cell::Core ( tag::is_positive, tag::of_dim, d ),
		meshes ( sz )
	{	}

	bool is_positive ( ) const;  // virtual from Cell::Core
	Cell::Core::Positive * get_positive ( );  // virtual from Cell::Core

	bool belongs_to ( Mesh::Core *, const tag::Oriented & ) const;  // virtual from Cell::Core
	bool belongs_to ( Mesh::Core *, const tag::NotOriented & ) const;  // virtual from Cell::Core

	inline void glue_common ( Cell::Core * face );
	inline void cut_common ( Cell::Core * face );
	// do not use directly; called from glue_on_my_bdry and cut_from_my_bdry
		
}; // end of class Cell::Core::Positive

//-----------------------------------------------------------------------------//

class Cell::Core::Negative : public Cell::Core

{	public :

	inline Negative ( const tag::OfDimension &, size_t d,
														const tag::ReverseOf &, Cell::Core::Positive * rev )
	: Cell::Core ( tag::is_negative, tag::reverse_of, rev, tag::of_dim, d )
	{	}

	inline Negative ( const tag::Ghost & )
	: Cell::Core ( tag::ghost )
	{	}

	bool is_positive ( ) const;  // virtual from Cell::Core
	Cell::Core::Positive * get_positive ( );  // virtual from Cell::Core

	Cell::Core * reverse ( const tag::BuildIfNotExists & );  // virtual from Cell::Core

	bool belongs_to ( Mesh::Core *, const tag::Oriented & ) const;  // virtual from Cell::Core
	bool belongs_to ( Mesh::Core *, const tag::NotOriented & ) const;  // virtual from Cell::Core

	void glue_on_my_bdry ( Cell::Core * ); // virtual from Cell::Core
	void cut_from_my_bdry ( Cell::Core * ); // virtual from Cell::Core
	// we feel that 'glue_on_bdry_of' and 'cut_from_bdry_of' are more readable
	// so we suggest to use those (see class Cell::Core)

#ifndef NDEBUG
	std::string get_name();
#endif
	
}; // end of class Cell::Core::Negative

//-----------------------------------------------------------------------------//


class Cell::Positive : public Cell::Core::Positive

// a cell of dimension >= 2

{	public :

	Mesh::Core * boundary_p;

	inline Positive ( const tag::OfDimension &, size_t d,
                        const tag::WhoseBoundaryIs &, Mesh::Core * msh );
	inline Positive ( const tag::WhoseBoundaryIs &, Mesh::Core * msh );
	inline Positive ( const tag::Triangle &, Cell::Core * AB,
	                        Cell::Core * BC, Cell::Core * CA );
	inline Positive ( const tag::Quadrangle &, Cell::Core * AB, Cell::Core * BC,
	                                           Cell::Core * CD, Cell::Core * DA );

	Positive ( const Cell::Positive & ) = delete;
	Positive ( const Cell::Positive && ) = delete;
	Cell::Positive & operator= ( const Cell::Positive & ) = delete;
	Cell::Positive & operator= ( const Cell::Positive && ) = delete;

	size_t get_dim ( ) const; // virtual
	Cell::Core * reverse ( const tag::BuildIfNotExists & ); // virtual from Cell::Core

	// methods 'add_to' and 'remove_from' add/remove 'this' cell to/from the mesh 'msh'
	// if 'msh' is the boundary of some cell, methods 'glue_on_bdry_of'
	// and 'cut_from_bdry_of' should be used instead (see class Cell::Core)
	void add_to ( Mesh::Core * msh ); // virtual from Cell::Core
	void remove_from ( Mesh::Core * msh ); // virtual from Cell::Core

	void glue_on_my_bdry ( Cell::Core * ); // virtual from Cell::Core
	void cut_from_my_bdry ( Cell::Core * ); // virtual from Cell::Core
	// we feel that 'glue_on_bdry_of' and 'cut_from_bdry_of' are more readable
	// so we suggest to use those (see class Cell::Core)

#ifndef NDEBUG
	void print_everything ( ); // virtual from Cell::Core
#endif

	class Vertex;  class Segment;
}; // end of  class Cell::Positive

//-----------------------------------------------------------------------------//


class Cell::Negative : public Cell::Core::Negative

// a cell of dimension >= 2

{	public :

	inline Negative ( const tag::OfDimension, size_t d,
                        const tag::ReverseOf &, Cell::Positive * direct_cell_p );
	inline Negative ( const tag::ReverseOf &, Cell::Positive * direct_cell_p );

	Negative ( const Cell::Negative & ) = delete;
	Negative ( const Cell::Negative && ) = delete;
	Cell::Negative & operator= ( const Cell::Negative & ) = delete;
	Cell::Negative & operator= ( const Cell::Negative && ) = delete;

	size_t get_dim ( ) const; // virtual from Cell::Core

	// methods 'add_to' and 'remove_from' add/remove 'this' cell to/from the mesh 'msh'
	// if 'msh' is the boundary of some cell, methods 'glue_on_bdry_of'
	// and 'cut_from_bdry_of' should be used instead (see class Cell::Core)
	void add_to ( Mesh::Core * msh ); // virtual
	void remove_from ( Mesh::Core * msh ); // virtual

#ifndef NDEBUG
	void print_everything ( ); // virtual from Cell::Core
#endif

	class Vertex;  class Segment;
	
}; // end of  class Cell::Negative

//-----------------------------------------------------------------------------//


class Cell::Positive::Vertex : public Cell::Core::Positive

{	public :

	inline Vertex ( );

	Vertex ( const Cell::Positive::Vertex & ) = delete;
	Vertex ( const Cell::Positive::Vertex && ) = delete;
	Cell::Positive::Vertex & operator= ( const Cell::Positive::Vertex & ) = delete;
	Cell::Positive::Vertex & operator= ( const Cell::Positive::Vertex && ) = delete;

	size_t get_dim ( ) const; // virtual from Cell::Core
	Cell::Core * reverse ( const tag::BuildIfNotExists & ); // virtual from Cell::Core
	
	// Methods 'add_to' and 'remove_from' add/remove 'this' cell to/from the mesh 'msh'.
	// If 'msh' is the boundary of some cell, methods 'glue_on_bdry_of'
	// and 'cut_from_bdry_of' should be used instead.
	void add_to ( Mesh::Core * msh ); // virtual from Cell::Core
	void remove_from ( Mesh::Core * msh ); // virtual from Cell::Core

	void glue_on_my_bdry ( Cell::Core * ); // virtual from Cell::Core
	void cut_from_my_bdry ( Cell::Core * ); // virtual from Cell::Core
	// we feel that 'glue_on_bdry_of' and 'cut_from_bdry_of' are more readable
	// so we suggest to use those (see class Cell::Core)

#ifndef NDEBUG
	void print_everything ( ); // virtual from Cell::Core
#endif
	
}; // end of  class Cell::Positive::Vertex

//-----------------------------------------------------------------------------//

class Cell::Negative::Vertex : public Cell::Core::Negative

{	public :

	inline Vertex ( const tag::ReverseOf &, Cell::Positive::Vertex * direct_ver_p );

	inline Vertex ( const tag::Ghost & )
	:	Cell::Core::Negative ( tag::ghost )
	{	}

	Vertex ( const Cell::Negative::Vertex & ) = delete;
	Vertex ( const Cell::Negative::Vertex && ) = delete;
	Cell::Negative::Vertex & operator= ( const Cell::Negative::Vertex & ) = delete;
	Cell::Negative::Vertex & operator= ( const Cell::Negative::Vertex && ) = delete;

	size_t get_dim ( ) const; // virtual from Cell::Core
	
	// methods 'add_to' and 'remove_from' add/remove 'this' cell to/from the mesh 'msh'
	// if 'msh' is the boundary of some cell, methods 'glue_on_bdry_of'
	// and 'cut_from_bdry_of' should be used instead (see class Cell::Core)
	void add_to ( Mesh::Core * msh ); // virtual from Cell::Core
	void remove_from ( Mesh::Core * msh ); // virtual from Cell::Core

#ifndef NDEBUG
	void print_everything ( ); // virtual from Cell::Core
#endif
	
}; // end of  class Cell::Negative::Vertex

//-----------------------------------------------------------------------------//

class Cell::Positive::Segment : public Cell::Core::Positive

{	public :

	Cell::Negative::Vertex * base_p;
	Cell::Positive::Vertex * tip_p;

	inline Segment ( Cell::Negative::Vertex * a, Cell::Positive::Vertex * b );

	Segment ( const Cell::Positive::Segment & ) = delete;
	Segment ( const Cell::Positive::Segment && ) = delete;
	Cell::Segment & operator= ( const Cell::Positive::Segment & ) = delete;
	Cell::Segment & operator= ( const Cell::Positive::Segment && ) = delete;

	Cell::Core * tip () override; // virtual, overrides definition by Cell::Core
	Cell::Core * base () override; // virtual, overrides definition by Cell::Core

	size_t get_dim ( ) const; // virtual from Cell::Core
	Cell::Core * reverse ( const tag::BuildIfNotExists & ); // virtual from Cell::Core

	// methods 'add_to' and 'remove_from' add/remove 'this' cell to/from the mesh 'msh'
	// if 'msh' is the boundary of some cell, methods 'glue_on_bdry_of'
	// and 'cut_from_bdry_of' should be used instead (see class Cell::Core)
	void add_to ( Mesh::Core * msh ); // virtual from Cell::Core
	void remove_from ( Mesh::Core * msh ); // virtual from Cell::Core

	void glue_on_my_bdry ( Cell::Core * ); // virtual from Cell::Core
	void cut_from_my_bdry ( Cell::Core * ); // virtual from Cell::Core
	// we feel that 'glue_on_bdry_of' and 'cut_from_bdry_of' are more readable
	// so we suggest to use those

	// here is the core linking between cells and meshes
	// do not use directly; this is called from add_to and remove_from
	void deep_connections
	(	Cell::Core::Positive * cll,
		void (*action) ( Cell::Core *, Cell::Core *,
		                 Mesh::Core *, short int, short int )  );

#ifndef NDEBUG
	void print_everything ( ); // virtual from Cell::Core
#endif
	
}; // end of  class Cell::Positive::Segment

//-----------------------------------------------------------------------------//

class Cell::Negative::Segment : public Cell::Core::Negative

{	public :

	inline Segment ( const tag::ReverseOf &, Cell::Positive::Segment * direct_seg_p );

	Segment ( const Cell::Negative::Segment & ) = delete;
	Segment ( const Cell::Negative::Segment && ) = delete;
	Cell::Negative::Segment & operator= ( const Cell::Negative::Segment & ) = delete;
	Cell::Negative::Segment & operator= ( const Cell::Negative::Segment && ) = delete;

	Cell::Core * tip () override;  // virtual, overrides definition by Cell::Core
	Cell::Core * base () override;  // virtual, overrides definition by Cell::Core

	size_t get_dim ( ) const; // virtual from Cell::Core
	
	// methods 'add_to' and 'remove_from' add/remove 'this' cell to/from the mesh 'msh'
	// if 'msh' is the boundary of some cell, methods 'glue_on_bdry_of'
	// and 'cut_from_bdry_of' should be used instead (see class Cell::Core)
	void add_to ( Mesh::Core * msh ); // virtual
	void remove_from ( Mesh::Core * msh ); // virtual from Cell::Core

#ifndef NDEBUG
	void print_everything ( ); // virtual from Cell::Core
#endif
	
}; // end of  class Cell::Negative::Segment

//-----------------------------------------------------------------------------//
//-----------------------------------------------------------------------------//


class Mesh::Core

// represents a positive mesh
// negative meshes have no core (wrappers for negative meshes are built on-the-fly)

{	public :

	Cell::Positive * cell_enclosed { nullptr };

	// the 'cells' attribute holds lists of cells of 'this' mesh, indexed by their dimension
	// for maximum dimension, the cells are oriented
	// for lower dimension, the cells are always positive

	// in the future, 'cells' will belong to Mesh::Positive
	// thus, Mesh::OneDim::Positive will have no 'cells' attribute
	// vamos ter que rever muita coisa em add_to, remove_from etc

	std::vector < std::list < Cell::Core* > > cells;

	inline Core ( const tag::OfDimension &, size_t d , const tag::MinusOne & )
	:	cells ( d )
	{	}
	
	virtual size_t get_dim_plus_one ( ) = 0;

	virtual size_t number_of ( const tag::CellsOfDim &, size_t d ) = 0;

	virtual Cell::Core * first_vertex ( ) = 0;
	virtual Cell::Core * last_vertex ( ) = 0;
	virtual Cell::Core * first_segment ( ) = 0;
	virtual Cell::Core * last_segment ( ) = 0;

	//private :

	virtual const Mesh::Methods & get_meth_pos() = 0;
	virtual const Mesh::Methods & get_meth_neg() = 0;
	
	// here is where the low-level linking between cells and meshes happens
	// do not use directly
	// deep_connections is called from add_to and remove_from
	void deep_connections
	( Cell::Core::Positive * cll,  Cell::Core * o_cll,
	  void (*action) ( Cell::Core *, Cell::Core *,
	                   Mesh::Core *, short int, short int ) );

#ifndef NDEBUG
	std::string name;
  virtual std::string get_name() = 0;
 	virtual void print_everything () = 0;
#endif
	
}; // end of  class Mesh::Core


class Mesh::Positive : public Mesh::Core

// represents a positive mesh of dimension > 1
// negative meshes have no core (wrappers for negative meshes are built on-the-fly)

{	public :
	
	static const Mesh::Methods methods_pos;
	static const Mesh::Methods methods_neg;
	
	inline Positive ( const tag::OfDimension &, size_t dim_p1, const tag::MinusOne & )
	:	Mesh::Core ( tag::of_dimension, dim_p1, tag::minus_one )
	{	}

  Positive ( const tag::Triangle &, const Mesh & AB, const Mesh & BC, const Mesh & CA );
	// defined in global.cpp

	inline Positive ( const tag::Quadrangle &, const Mesh & south, const Mesh & east,
	  const Mesh & north, const Mesh & west, const tag::WithTriangles & wt = tag::not_with_triangles );

	inline Positive ( const tag::Quadrangle &, const Cell & SW, const Cell & SE, const Cell & NE,
		const Cell & NW, size_t m, size_t n, const tag::WithTriangles & wt = tag::not_with_triangles );

	size_t get_dim_plus_one ( );  // virtual from Mesh::Core

	size_t number_of ( const tag::CellsOfDim &, size_t d );  // virtual from Mesh::Core

	Cell::Core * first_vertex ( );  // virtual from Mesh::Core
	// returns a negative vertex
	Cell::Core * last_vertex ( );  // virtual from Mesh::Core
	// returns a positive vertex
	Cell::Core * first_segment ( );  // virtual from Mesh::Core
	Cell::Core * last_segment ( );  // virtual from Mesh::Core

	void build_rectangle ( const Mesh & south, const Mesh & east,
		const Mesh & north, const Mesh & west, bool cut_rectangles_in_half );
	// defined in global.cpp

	void forget () { };
	// remove_from this mesh each cell of maximum dimension
	// if the cell has no other meshes above, forget it
	// delete this (now empty) mesh
	
	static bool is_positive ( );
	static Mesh reverse ( Mesh::Core * core );

	// private :
	
	const Mesh::Methods & get_meth_pos(); // virtual
	const Mesh::Methods & get_meth_neg(); // virtual
	
	inline CellIterator iter_over ( const tag::CellsOfDim &, size_t d );
	
#ifndef NDEBUG
  std::string get_name();  // virtual from Mesh::Core
 	void print_everything ();  // virtual from Mesh::Core
#endif
	
}; // end of  class Mesh::Positive


class Mesh::OneDim::Positive : public Mesh::Core

{	public :

	Cell::Negative::Vertex * first_ver { nullptr };
	Cell::Positive::Vertex * last_ver;

	// we will have here a list of segments (no vertices kept)

	static const Mesh::Methods methods_pos;
	static const Mesh::Methods methods_neg;
	
	inline Positive ( )
	:	Mesh::Core ( tag::of_dimension, 2, tag::minus_one )
	{ }

	Positive ( const tag::Segment &, Cell::Positive::Vertex * A, Cell::Positive::Vertex * B,
             const tag::DividedIn &, size_t n );
	// defined in global.cpp
	
	size_t get_dim_plus_one ( );  // virtual from Mesh::Core

	size_t number_of ( const tag::CellsOfDim &, size_t d );  // virtual from Mesh::Core

	Cell::Core * first_vertex ( );  // virtual from Mesh::Core
	Cell::Core * last_vertex ( );  // virtual from Mesh::Core
	Cell::Core * first_segment ( );  // virtual from Mesh::Core
	Cell::Core * last_segment ( );  // virtual from Mesh::Core

	// private:
	
	const Mesh::Methods & get_meth_pos();  // virtual from Mesh::Core
	const Mesh::Methods & get_meth_neg();  // virtual from Mesh::Core

	inline void order ( );
	// run over all segments, order them linearly, check that the mesh is connected
	// if it is a loop, set first_ver = Cell::ghost
	// if it is an open chain, set first_ver and last_ver
	// see paragraph 8.14 in the manual [?]
	
#ifndef NDEBUG
	std::string get_name();  // virtual from Mesh::Core
	void print_everything ();  // virtual from Mesh::Core
#endif

}; // end of  class Mesh::OneDim::Positive

//-----------------------------------------------------------------------------//


// negative meshes aren't kept in the computer, they are just abstract concepts
// to be more precise: there are Mesh objects (wrappers) refering to negative meshes
// but their core points to a positive Mesh::Core; there are no negative Mesh::Corees
// we declare the Mesh::Negative namespace however to mimick static methods of a class

struct Mesh::Negative

{	// mimick virtual methods
	static bool is_positive ( );
	static Mesh reverse ( Mesh::Core * core );
	static inline CellIterator iter_over ( const tag::CellsOfReverseOf &,
		Mesh::Core * that, const tag::OfDimension &, size_t     );           };

//-----------------------------------------------------------------------------//


inline Mesh::Mesh ( const tag::OfDimension &, size_t d, const tag::GreaterThanOne &,
                    const tag::IsPositive & ispos                                     )
// by default, ispos = tag::is_positive, so may be called with only three arguments
:	core { new Mesh::Positive ( tag::of_dimension, d+1, tag::minus_one ) },
	meth { & Mesh::Positive::methods_pos }
{	assert ( d > 1 );  }
	

inline Mesh::Mesh ( const tag::OfDimension &, size_t d, const tag::MightBeOne &,
                    const tag::IsPositive & ispos                                     )
// by default, ispos = tag::is_positive, so may be called with only three arguments
:	core { nullptr },
	meth { & Mesh::Positive::methods_pos }
{	if ( d > 1 )
		this->core = new Mesh::Positive ( tag::of_dimension, d+1, tag::minus_one );
	else
	{	assert ( d == 1 );
		this->core = new Mesh::OneDim::Positive ( );  }                               }
	

inline Mesh::Mesh ( const tag::OfDimensionOne &, const tag::IsPositive & ispos )
// by default, ispos = tag::is_positive, so may be called with only one argument
:	core { new Mesh::OneDim::Positive ( ) }, meth { & Mesh::OneDim::Positive::methods_pos }
{ }
	

inline Mesh::Mesh ( const tag::WhoseCoreIs &, Mesh::Core * c, const tag::IsPositive & ispos )
// by default, ispos = tag::is_positive, so may be called with only two arguments
// used in Mesh::Negative::reverse and Cell::boundary
: core { c }, meth { & c->get_meth_pos() }
{	assert ( c );
	assert ( ( meth == & Mesh::Positive::methods_pos ) or
	         ( meth == & Mesh::OneDim::Positive::methods_pos ) );  }


inline Mesh::Mesh ( const tag::WhoseCoreIs &, Mesh::Core * c, const tag::IsNegative &,
                    const tag::SurelyExists &                                          )
// builds a negative mesh from a positive one, assuming that reverse cells exist
// used in Cell::boundary and Mesh::Mesh below
: core { c }, meth { & c->get_meth_neg() }
{	assert ( c );
	assert ( ( meth == & Mesh::Positive::methods_neg ) or
	         ( meth == & Mesh::OneDim::Positive::methods_neg ) );  }


inline Mesh::Mesh ( const tag::WhoseCoreIs &, Mesh::Core * c, const tag::IsNegative &,
                    const tag::BuildCellsIfNecessary & b                                )
// builds a negative mesh from a positive one, creating reverse cells if necessary
// used in Mesh::Positive::reverse
:	Mesh ( tag::whose_core_is, c, tag::is_negative, tag::surely_exists )
{	std::list < Cell::Core * > & cls = c->cells.back();  // maximum dimension
	std::list<Cell::Core*>::iterator it = cls.begin(), it_e = cls.end();
	// cell[0] means cells of same dimension as the mesh
	for ( ; it != it_e; it++ )
	{	Cell::Core * cll_p = *it;
		assert ( cll_p );
		Cell::Core * cll_rev_p = cll_p->reverse ( tag::build_if_not_exists );
		assert ( cll_rev_p );  assert ( cll_rev_p == cll_p->reverse_p );             }     }


// In Mesh ( tag::segment... ), the first two arguments are positive points.
// This is not consistent with Cell ( tag::segment... ), where the user
// must provide a negative point (the base) then a positive point (the tip).
// It is also inconsistent with other constructors like Mesh ( tag::quadrangle... )
// where we provide faces with orientation compatible with the orientation
// of the future mesh.
// However, we think it is easier for the user to build chains of segments like
//   Cell A ( tag::vertex ); Cell B ( tag::point );
//   Mesh AB ( tag::segment, A, B, tag::divided_in, 10 );
// rather than
//   Cell A ( tag::vertex ); Cell B ( tag::point );
//   Mesh AB ( tag::segment, A.reverse(), B, tag::divided_in, 10 );

inline Mesh::Mesh ( const tag::Segment &, const Cell & A, const Cell & B, const tag::DividedIn &, size_t n )
: Mesh ( tag::whose_core_is, new Mesh::OneDim::Positive
  ( tag::segment, (Cell::Positive::Vertex*) A.core, (Cell::Positive::Vertex*) B.core, tag::divided_in, n ) )
{	}


inline Mesh::Mesh
( const tag::Pretty &, const tag::Segment &, const Cell & A, const Cell & B,
  const tag::DividedIn &, size_t n )
:	Mesh ( tag::of_dimension_one )
{	this->pretty_constructor ( tag::segment, A, B, tag::divided_in, n );  }


inline Mesh::Mesh
(	const tag::Triangle &, const Mesh & AB, const Mesh & BC, const Mesh & CA )
: Mesh ( tag::whose_core_is, new Mesh::Positive ( tag::triangle, AB, BC, CA ) )
{	}


inline Mesh::Mesh
(	const tag::Pretty &, const tag::Triangle &, const Mesh & AB, const Mesh & BC, const Mesh & CA )
:	Mesh ( tag::of_dimension, 2, tag::greater_than_one )
{	this->pretty_constructor ( tag::triangle, AB, BC, CA );  }


inline Mesh::Mesh
(	const tag::Pretty &, const tag::Quadrangle &,
	const Mesh & south, const Mesh & east, const Mesh & north, const Mesh & west,
	const tag::WithTriangles & wt                                                  )
// 'wt' defaults to 'tag::not_with_triangles',
:	Mesh ( tag::of_dimension, 2, tag::greater_than_one )
{	this->pretty_constructor ( tag::quadrangle, south, east, north, west, wt );  }


inline Mesh::Mesh
(	const tag::Quadrangle &,
	const Mesh & south, const Mesh & east, const Mesh & north, const Mesh & west,
	const tag::WithTriangles & wt                                                  )
// 'wt' defaults to 'tag::not_with_triangles',
:	Mesh ( tag::whose_core_is, new Mesh::Positive
	        ( tag::quadrangle, south, east, north, west, wt ) )
{	}

inline Mesh::Mesh ( const tag::Quadrangle &, const Cell & SW, const Cell & SE,
                    const Cell & NE, const Cell & NW, size_t m, size_t n,
                    const tag::WithTriangles & wt                               )
// 'wt' defaults to 'tag::not_with_triangles',
:	Mesh ( tag::whose_core_is, new Mesh::Positive ( tag::quadrangle, SW, SE, NE, NW, m ,n, wt ) )
{	}


inline Mesh::Mesh ( const tag::Join &, const Mesh & m1, const Mesh & m2 )
:	Mesh ( tag::of_dim, m1.dim(), tag::might_be_one )
{	std::list < Mesh > l;
	l.push_back ( m1 );  l.push_back ( m2 );
  this->join_list ( l );                   }

inline Mesh::Mesh ( const tag::Join &, const Mesh & m1, const Mesh & m2, const Mesh & m3 )
:	Mesh ( tag::of_dim, m1.dim(), tag::might_be_one )
{	std::list < Mesh > l;
	l.push_back ( m1 );  l.push_back ( m2 );  l.push_back ( m3 );
  this->join_list ( l );                                         }

inline Mesh::Mesh
( const tag::Join &, const Mesh & m1, const Mesh & m2, const Mesh & m3, const Mesh & m4 )
:	Mesh ( tag::of_dim, m1.dim(), tag::might_be_one )
{	std::list < Mesh > l;
	l.push_back ( m1 );  l.push_back ( m2 );
	l.push_back ( m3 );  l.push_back ( m4 );
  this->join_list ( l );                   }

inline Mesh::Mesh ( const tag::Join &, const Mesh & m1, const Mesh & m2, const Mesh & m3,
                                       const Mesh & m4, const Mesh & m5 )
:	Mesh ( tag::of_dim, m1.dim(), tag::might_be_one )
{	std::list < Mesh > l;
	l.push_back ( m1 );  l.push_back ( m2 );  l.push_back ( m3 );
	l.push_back ( m4 );  l.push_back ( m5 );
  this->join_list ( l );                                         }

inline Mesh::Mesh ( const tag::Join &, const std::list<Mesh> & l )
:	Mesh ( tag::of_dim, l.front().dim(), tag::might_be_one )
{	this->join_list ( l );  }


//-----------------------------------------------------------------------------//


inline Cell::Cell ( const tag::WhoseBoundaryIs &, Mesh & msh )
:	Cell ( tag::whose_core_is, new Cell::Positive ( tag::whose_boundary_is, msh.core ) )
{	assert ( msh.is_positive() );  }


// inline Cell::Cell ( const tag::WhoseBoundaryIs &, Mesh::Core * msh )
// :	Cell ( tag::whose_core_is, new Cell::Positive ( tag::whose_boundary_is, msh ) )
// { }

inline Cell::Cell
( const tag::ReverseOf &, const Cell & direct_cell, const tag::BuildIfNotExists & build )
// 'build' defaults to tag::build_if_not_exists
// so constructor may be called with only two arguments
:	Cell ( tag::whose_core_is, direct_cell.core->reverse ( tag::build_if_not_exists ) )
{	}

inline Cell::Cell ( const tag::ReverseOf &, const Cell & direct_cell, const tag::SurelyExists & )
:	Cell ( tag::whose_core_is, direct_cell.core->reverse_p )
{	}

inline Cell::Cell ( const tag::ReverseOf &, const Cell & direct_cell, const tag::MayNotExist & )
:	core ( direct_cell.core->reverse_p )
{	}


inline Cell::Cell ( const tag::BehindFace &, const Cell & face,
	const tag::WithinMesh &, const Mesh & msh, const tag::SurelyExists & se )
// 'se' defaults to tag::surely_exists, so constructor may be called with only four arguments
:	core { msh.cell_behind ( face.core, tag::may_not_exist ) }
{	assert ( this->core );   }
	

inline Cell::Cell ( const tag::BehindFace &, const Cell & face,
	const tag::WithinMesh &, const Mesh & msh, const tag::MayNotExist & )
:	core { msh.cell_behind ( face.core, tag::may_not_exist ) }
{	}
	

inline Cell::Cell ( const tag::InFrontOfFace &, const Cell & face,
	const tag::WithinMesh &, const Mesh & msh, const tag::SurelyExists & se )
// 'se' defaults to tag::surely_exists, so constructor may be called with only four arguments
:	core { msh.cell_in_front_of ( face.core, tag::may_not_exist ) }
{	assert ( this->core );   }
		

inline Cell::Cell ( const tag::InFrontOfFace &, const Cell & face,
	const tag::WithinMesh &, const Mesh & msh, const tag::MayNotExist & )
:	core { msh.cell_in_front_of ( face.core, tag::may_not_exist ) }
{	}
		

inline Cell::Cell ( const tag::Vertex &, const tag::IsPositive & ispos )
// by default, ispos = tag::is_positive, so may be called with only one argument
:	Cell ( tag::whose_core_is, new Cell::Positive::Vertex )
{	}


inline Cell::Cell ( const tag::Segment &, const Cell & A, const Cell & B )
: Cell ( tag::whose_core_is, new Cell::Positive::Segment
( (Cell::Negative::Vertex*) A.core, (Cell::Positive::Vertex*) B.core ) )

{	assert ( not A.is_positive() );
	assert ( B.is_positive() );       }


inline Cell::Cell ( const tag::Triangle &, const Cell & AB, const Cell & BC, const Cell & CA )
:	Cell ( tag::whose_core_is, new Cell::Positive
            ( tag::triangle, AB.core, BC.core, CA.core ) )
{	}

inline Cell::Cell ( const tag::Quadrangle &, const Cell & AB, const Cell & BC,
                                             const Cell & CD, const Cell & DA )
:	Cell ( tag::whose_core_is, new Cell::Positive
            ( tag::quadrangle, AB.core, BC.core, CD.core, DA.core ) )
{	}

//-----------------------------------------------------------------------------//


inline size_t Mesh::number_of ( const tag::CellsOfDim &, size_t d ) const
{	return this->core->number_of ( tag::cells_of_dim, d );  }

inline size_t Mesh::number_of ( const tag::Segments & ) const
{	return this->number_of ( tag::cells_of_dim, 1 );  }

inline size_t Mesh::number_of ( const tag::Vertices & ) const
{	return this->number_of ( tag::cells_of_dim, 0 );  }


// perhaps include the four methods below in Mesh::meth in order to speed up the code

inline Cell Mesh::first_vertex ( ) const
{	assert ( this->dim() == 1 );
	if ( this->is_positive() )
		return Cell ( tag::whose_core_is, this->core->first_vertex() );
	else
		return Cell ( tag::whose_core_is, this->core->last_vertex() );  }

inline Cell Mesh::last_vertex ( ) const
{	assert ( this->dim() == 1 );
	if ( this->is_positive() )
		return Cell ( tag::whose_core_is, this->core->last_vertex() );
	else
		return Cell ( tag::whose_core_is, this->core->first_vertex() );  }

inline Cell Mesh::first_segment ( ) const
{	assert ( this->dim() == 1 );
	if ( this->is_positive() )
		return Cell ( tag::whose_core_is, this->core->first_segment() );
	else
		return Cell ( tag::whose_core_is, this->core->last_segment()->reverse_p );  }

inline Cell Mesh::last_segment ( ) const
{	assert ( this->dim() == 1 );
	if ( this->is_positive() )
		return Cell ( tag::whose_core_is, this->core->last_segment() );
	else
		return Cell ( tag::whose_core_is, this->core->first_segment()->reverse_p );  }


inline Cell Mesh::cell_in_front_of ( const Cell face, const tag::SurelyExists & se ) const
// 'se' defaults to tag::surely_exists, so method may be called with only one argument
{	return Cell ( tag::in_front_of_face, face, tag::within_mesh, *this, tag::surely_exists );  }

inline Cell Mesh::cell_in_front_of ( const Cell face, const tag::MayNotExist & ) const
{	return Cell ( tag::in_front_of_face, face, tag::within_mesh, *this, tag::may_not_exist );  }

inline Cell Mesh::cell_behind ( const Cell face, const tag::SurelyExists & se ) const
// 'se' defaults to tag::surely_exists, so method may be called with only one argument
{	return Cell ( tag::behind_face, face, tag::within_mesh, *this, tag::surely_exists );  }

inline Cell Mesh::cell_behind ( const Cell face, const tag::MayNotExist & ) const
{	return Cell ( tag::behind_face, face, tag::within_mesh, *this, tag::may_not_exist );  }


inline bool Cell::belongs_to ( const Mesh & msh, const tag::Oriented & ) const
{	return this->core->belongs_to ( msh.core, tag::oriented );  }

inline bool Cell::belongs_to ( const Mesh & msh, const tag::NotOriented & ) const
{	return this->core->belongs_to ( msh.core, tag::not_oriented );  }

//-----------------------------------------------------------------------------//


inline Mesh::Positive::Positive ( const tag::Quadrangle &, const Mesh & south,
  const Mesh & east, const Mesh & north, const Mesh & west, const tag::WithTriangles & wt )

// 'wt' defaults to 'tag::not_with_triangles',
// which means 'cut_rectangles_in_half' defaults to 'false'

:	Mesh::Core ( tag::of_dimension, 3, tag::minus_one )
	
{	bool cut_rectangles_in_half = ( wt == tag::with_triangles );

	this->build_rectangle ( south, east, north, west, cut_rectangles_in_half );  }


inline Mesh::Positive::Positive ( const tag::Quadrangle &, const Cell & SW, const Cell & SE,
	const Cell & NE, const Cell & NW, size_t m, size_t n, const tag::WithTriangles & wt )

// 'wt' defaults to 'tag::not_with_triangles',
// which means 'cut_rectangles_in_half' defaults to 'false'

:	Mesh::Core ( tag::of_dimension, 3, tag::minus_one )
	
{	bool cut_rectangles_in_half = ( wt == tag::with_triangles );

	Mesh south ( tag::segment, SW, SE, tag::divided_in, m );
	Mesh east  ( tag::segment, SE, NE, tag::divided_in, n );
	Mesh north ( tag::segment, NE, NW, tag::divided_in, m );
	Mesh west  ( tag::segment, NW, SW, tag::divided_in, n );

	// when these four meshes go out of scope, their core should be disposed of

	this->build_rectangle ( south, east, north, west, cut_rectangles_in_half );  }

//-----------------------------------------------------------------------------//


inline Cell::Core * Mesh::cell_behind
( const Cell::Core * face_p, const tag::MayNotExist & ) const

// return the cell to which 'face' belongs

{	assert ( this->dim() == face_p->get_dim() + 1 );
	if ( this->is_positive() )
	{	std::map<Mesh::Core*,Cell::Core*>::const_iterator
			it = face_p->cell_behind_within.find ( this->core );
		if ( it == face_p->cell_behind_within.end() ) return nullptr;
		return it->second;                                                }
		// face_p->cell_behind_within[this->core]
	else
	{	Cell::Core * face_rev_p = face_p->reverse_p;
		if ( face_rev_p == nullptr ) return nullptr;
		std::map<Mesh::Core*,Cell::Core*>::const_iterator
			it = face_rev_p->cell_behind_within.find ( this->core );
		if ( it == face_rev_p->cell_behind_within.end() ) return nullptr;
		Cell::Core * cll_rev_p = it->second;
		assert ( cll_rev_p );
		return cll_rev_p->reverse_p;                                        }  }


inline Cell::Core * Mesh::cell_behind
( const Cell::Core * face_p, const tag::SurelyExists & se ) const

// 'se' defaults to tag::surely_exists, so method may be called with only one argument

// return the cell to which 'face' belongs

{	assert ( this->dim() == face_p->get_dim() + 1 );
	if ( this->is_positive() )
	{	std::map<Mesh::Core*,Cell::Core*>::const_iterator
			it = face_p->cell_behind_within.find ( this->core );
		assert ( it != face_p->cell_behind_within.end() );
		assert ( it->second );
		return it->second;                                                }
		// face_p->cell_behind_within[this->core]
	else
	{	Cell::Core * face_rev_p = face_p->reverse_p;
		assert ( face_rev_p );
		std::map<Mesh::Core*,Cell::Core*>::const_iterator
			it = face_rev_p->cell_behind_within.find ( this->core );
		assert ( it != face_rev_p->cell_behind_within.end() );
		Cell::Core * cll_rev_p = it->second;
		assert ( cll_rev_p );  assert ( cll_rev_p->reverse_p );
		return cll_rev_p->reverse_p;                                        }  }


inline Cell::Core * Mesh::cell_in_front_of
( const Cell::Core * face_p, const tag::MayNotExist & ) const

// return the cell towards which 'face' is looking
// recall that the faces of a cell are looking outwards

{	Cell::Core * face_rev = face_p->reverse_p;
	if ( face_rev == nullptr ) return nullptr;
	else return this->cell_behind ( face_rev, tag::may_not_exist );  }
	

inline Cell::Core * Mesh::cell_in_front_of
( const Cell::Core * face_p, const tag::SurelyExists & se ) const

// 'se' defaults to tag::surely_exists, so method may be called with only one argument

// return the cell towards which 'face' is looking
// recall that the faces of a cell are looking outwards

{	Cell::Core * face_rev = face_p->reverse_p;
	assert ( face_rev );
	return this->cell_behind ( face_rev, tag::surely_exists );  }
	
// use :
// Cell f (...);
// Cell A ( tag::in_front_of_face, f, tag::within_mesh, msh );
// Cell A = msh.cell_in_front_of ( f.core );
// Cell::Core * A_p = msh.cell_in_front_of ( f.core, tag::may_not_exist );


#ifndef NDEBUG
inline void Mesh::print_everything ( )
{	if ( not is_positive() ) std::cout << "(negative Mesh) ";
	core->print_everything ();                                }
#endif

//-----------------------------------------------------------------------------//


inline bool Cell::is_positive ( ) const
{	return this->core->is_positive ( );  }

inline bool Mesh::is_positive ( ) const
{	return this->meth->is_positive ( );  }

inline Cell Cell::get_positive ( )
{	return Cell ( tag::whose_core_is, this->core->get_positive() );  }


inline size_t Cell::dim ( ) const
{	return this->core->get_dim ( );  }

inline size_t Mesh::dim ( ) const
{	return Mesh::diff ( this->core->get_dim_plus_one(), 1 );      }
// Mesh::diff  provides a safe way to substract two size_t numbers


inline Cell Cell::reverse ( const tag::BuildIfNotExists & build ) const
// 'build' defaults to tag::build_if_not_exists, so method may be called with no arguments
{	return Cell ( tag::reverse_of, *this, tag::build_if_not_exists );  }

inline Cell Cell::reverse ( const tag::MayNotExist & ) const
{	return Cell ( tag::reverse_of, *this, tag::may_not_exist );  }

inline Cell Cell::reverse ( const tag::SurelyExists & ) const
{	return Cell ( tag::reverse_of, *this, tag::surely_exists );  }

inline Mesh Mesh::reverse ( ) const
{	return this->meth->reverse ( this->core );  }

inline bool Cell::has_reverse ( ) const
{	return this->core->reverse_p;  }


inline Mesh Cell::boundary ( ) const

{	assert ( this->core );
	assert ( this->dim() > 1 );  // there are no zero-dimensional meshes
	if ( this->is_positive() )
	{	Cell::Positive * cll = (Cell::Positive*) this->core;
		return Mesh ( tag::whose_core_is, cll->boundary_p, tag::is_positive );  }
	else // negative mesh, boundary of negative cell, faces already have reverse
	{	assert ( this->core->reverse_p );
		Cell::Positive * cll = (Cell::Positive*) this->core->reverse_p;
		return Mesh ( tag::whose_core_is, cll->boundary_p, tag::is_negative, tag::surely_exists );  }  }


inline Cell Cell::tip () const
{	// assert ( this->core->tip() );
	return Cell ( tag::whose_core_is, this->core->tip() );  }

inline Cell Cell::base () const
{	// assert ( this->core->base() );
	return Cell ( tag::whose_core_is, this->core->base() );  }


#ifndef NDEBUG
inline void Cell::print_everything ( )
{	core->print_everything ( );  }
#endif


inline void Cell::glue_on_bdry_of ( Cell & cll )

// glue 'this' face on the boundary of cell 'cll'
// any of them may be negative

{	this->core->glue_on_bdry_of ( cll.core );  }


inline void Cell::cut_from_bdry_of ( Cell & cll )

// cut 'this' face from the boundary of cell 'cll'
// any of them may be negative

{	this->core->cut_from_bdry_of ( cll.core );  }


inline void Cell::add_to ( Mesh & msh )

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

{	assert ( this->dim() == msh.dim() );
	if ( msh.is_positive() )  this->core->add_to ( msh.core );
	else
	{	assert ( this->core->reverse_p );
		this->core->reverse_p->add_to ( msh.core );  }             }
// for negative Meshes, the core points towards the reverse Mesh::Core


inline void Cell::remove_from ( Mesh & msh )

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

{	assert ( this->dim() == msh.dim() );
	if ( msh.is_positive() )  this->core->remove_from ( msh.core );
	else
	{	assert ( this->core->reverse_p );
		this->core->reverse_p->remove_from ( msh.core );  }             }
// for negative Meshes, the core points towards the reverse Mesh::Core


//-----------------------------------------------------------------------------//


inline Cell::Positive::Vertex::Vertex ( )
: Cell::Core::Positive ( tag::of_dim, 0, tag::size_meshes, Mesh::maximum_dimension_plus_one )
{	}


inline Cell::Negative::Vertex::Vertex
( const tag::ReverseOf &, Cell::Positive::Vertex * direct_ver_p )
: Cell::Core::Negative ( tag::of_dim, 0, tag::reverse_of, direct_ver_p )
{	}


inline Cell::Positive::Segment::Segment
( Cell::Negative::Vertex * Aa, Cell::Positive::Vertex * Bb )

: Cell::Core::Positive
		( tag::of_dim, 1, tag::size_meshes, Mesh::diff ( Mesh::maximum_dimension_plus_one, 1 ) ),
	// Mesh::diff provides a safe way to substract two 'size_t' numbers
	base_p { Aa }, tip_p { Bb }

{	Mesh::Core * msh = (Mesh::Core*) this;
	// a Positive::Segment disguised as a zero-dimensional Mesh::Core
	// below is a much simplified version of Negative::Vertex::add_to
	// that's because 'this' segment has just been created, so it has no meshes above
	// also, the base has already been correctly initialized
	assert ( Aa->reverse_p );
	Cell::Positive::Vertex * pos_Aa { (Cell::Positive::Vertex*) Aa->reverse_p };
	assert ( pos_Aa->meshes.size() > 0 );
	assert ( pos_Aa->meshes[0].find(msh) == pos_Aa->meshes[0].end() );
	// pos_Aa->meshes[0][msh] = Cell::field_to_meshes { 0, 1 };
	// the third component 'where' is irrelevant here
	pos_Aa->meshes[0].emplace ( std::piecewise_construct,
	  std::forward_as_tuple(msh), std::forward_as_tuple(0,1) );
	// below is a much simplified version of Positive::Vertex::add_to
	// that's because 'this' segment has just been created, so it has no meshes above
	// also, the tip has already been correctly initialized
	assert ( Bb->meshes.size() > 0 );
	assert ( Bb->meshes[0].find(msh) == Bb->meshes[0].end() );
	// Bb->meshes[0][msh] = Cell::field_to_meshes { 1, 0 };
	// the third component 'where' is irrelevant here
	Bb->meshes[0].emplace ( std::piecewise_construct,
	  std::forward_as_tuple(msh), std::forward_as_tuple(1,0) );                         }


inline Cell::Negative::Segment::Segment
( const tag::ReverseOf &, Cell::Positive::Segment * direct_seg_p )

: Cell::Core::Negative ( tag::of_dim, 1, tag::reverse_of, direct_seg_p )

// we must make sure that both extremities of 'direct_seg_p' have a reverse
// well, the base surely has one since it's a Negative::Vertex

{	assert ( direct_seg_p->base_p );
	assert ( direct_seg_p->tip_p );
	assert ( direct_seg_p->base_p->reverse_p );
	direct_seg_p->tip_p->reverse ( tag::build_if_not_exists );
	assert ( direct_seg_p->tip_p->reverse_p );                  }

		
inline Cell::Positive::Positive
( const tag::OfDimension &, size_t d, const tag::WhoseBoundaryIs &, Mesh::Core * msh )

:	Cell::Core::Positive ( tag::of_dim, d, tag::size_meshes,
	                           Mesh::diff ( Mesh::maximum_dimension_plus_one, d ) ),
	// Mesh::diff provides a safe way to substract two 'size_t' numbers
	boundary_p ( msh )

{	assert ( msh );
	assert ( msh->get_dim_plus_one() == d );
	msh->cell_enclosed = this;                }


inline Cell::Positive::Positive
( const tag::WhoseBoundaryIs &, Mesh::Core * msh )
:	Cell::Positive ( tag::of_dim, msh->get_dim_plus_one(), tag::whose_bdry_is, msh )
{	}


inline Cell::Negative::Negative
( const tag::OfDimension, size_t d, const tag::ReverseOf &, Cell::Positive * direct_cell_p )
	
: Cell::Core::Negative ( tag::of_dim, d, tag::reverse_of, direct_cell_p )

// we must make sure that all faces of 'direct_cell_p' have a reverse

{	assert ( direct_cell_p );
	assert ( direct_cell_p->get_dim() == d );
	assert ( direct_cell_p->boundary_p );
	std::list < Cell::Core * > & cls = direct_cell_p->boundary_p->cells.back();
	// cells of same dimension as the mesh (in this case, faces of 'this')
	std::list<Cell::Core*>::iterator it = cls.begin(), it_e = cls.end();
	for ( ; it != it_e; it++ )
	{	Cell::Core * cll_p = *it;
		assert ( cll_p );
		Cell::Core * cll_rev_p = cll_p->reverse ( tag::build_if_not_exists );
		assert ( cll_rev_p == cll_p->reverse_p );  assert ( cll_rev_p );             }    }


inline Cell::Negative::Negative ( const tag::ReverseOf &, Cell::Positive * direct_cell_p )
: Cell::Negative ( tag::of_dim, direct_cell_p->get_dim(), tag::reverse_of, direct_cell_p )
{	}


inline Cell::Positive::Positive ( const tag::Triangle &,
	Cell::Core * AB, Cell::Core * BC, Cell::Core * CA )
	
: Cell::Positive ( tag::whose_boundary_is, new Mesh::OneDim::Positive )

{	assert ( AB->get_dim() == 1 );
	assert ( BC->get_dim() == 1 );
	assert ( CA->get_dim() == 1 );
	assert ( AB->tip() == BC->base()->reverse_p );
	assert ( BC->tip() == CA->base()->reverse_p );
	assert ( CA->tip() == AB->base()->reverse_p );
	// no need for glue_on_bdry_of : 'this' cell has just been created, it has no meshes above
	assert  ( this->boundary_p );
	AB->add_to ( this->boundary_p );
	BC->add_to ( this->boundary_p );
	CA->add_to ( this->boundary_p );                }
	

inline Cell::Positive::Positive ( const tag::Quadrangle &,
	Cell::Core * AB, Cell::Core * BC, Cell::Core * CD, Cell::Core * DA )
	
: Cell::Positive ( tag::whose_boundary_is, new Mesh::OneDim::Positive )

{	assert ( AB->get_dim() == 1 );
	assert ( BC->get_dim() == 1 );
	assert ( CD->get_dim() == 1 );
	assert ( DA->get_dim() == 1 );
	assert ( AB->tip() == BC->base()->reverse_p );
	assert ( BC->tip() == CD->base()->reverse_p );
	assert ( CD->tip() == DA->base()->reverse_p );
	assert ( DA->tip() == AB->base()->reverse_p );
	// no need for glue_on_bdry_of : 'this' cell has just been created, it has no meshes above
	assert  ( this->boundary_p );
	AB->add_to ( this->boundary_p );
	BC->add_to ( this->boundary_p );
	CD->add_to ( this->boundary_p );
	DA->add_to ( this->boundary_p );                }


inline void Cell::Core::Positive::glue_common ( Cell::Core * face )
	
{	if ( this->meshes.size() == 0 ) return;
	std::map < Mesh::Core*, Cell::field_to_meshes > & tm0 = this->meshes[0];
	// '0' means the same dimension as the cell
	std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
		it = tm0.begin(), it_e = tm0.end();
	for ( ; it != it_e; ++it )
	{	Mesh::Core * msh = it->first;
		std::list<Cell::Core*>::iterator wh = it->second.where;
		Cell::Core * other_cell = *wh;  assert ( other_cell );
		if ( other_cell == this )  // orientations match
			face->cell_behind_within[msh] = this;
		else  // mismatched orientations
		{	assert ( other_cell == this->reverse_p );
			Cell::Core * rev_face { face->reverse_p };
			assert ( rev_face );
			rev_face->cell_behind_within[msh] = other_cell;   }           }                 }
	

inline void Cell::Core::cut_from_bdry_of ( Cell::Core * cll )
{	cll->cut_from_my_bdry ( this );   }


inline void Cell::Core::Positive::cut_common ( Cell::Core * face )
	
{	if ( this->meshes.size() == 0 ) return;
	std::map < Mesh::Core *, Cell::field_to_meshes > & tm0 = this->meshes[0];
	// '0' means the same dimension as the cell
	std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
		it = tm0.begin(), it_e = tm0.end();
	for ( ; it != it_e; ++it )
	{	Mesh::Core * msh = it->first;
		std::list<Cell::Core*>::iterator wh = it->second.where;
		Cell::Core * other_cell = *wh;  assert ( other_cell );
		if ( other_cell == this )  // orientations match
		{	assert ( face->cell_behind_within.find(msh) !=
		           face->cell_behind_within.end()        );
			assert ( face->cell_behind_within[msh] == this );
			face->cell_behind_within.erase(msh);                  }
		else  // mismatched orientations
		{	assert ( other_cell == this->reverse_p );
			Cell::Core * rev_face { face->reverse_p };
			assert ( rev_face );
			assert ( rev_face->cell_behind_within.find(msh) !=
		           rev_face->cell_behind_within.end()        );
			assert ( rev_face->cell_behind_within[msh] == other_cell );
			rev_face->cell_behind_within.erase(msh);                     }  }                  }
	

inline void Cell::Core::glue_on_bdry_of ( Cell::Core * cll )
{	cll->glue_on_my_bdry ( this );   }

//-----------------------------------------------------------------------------//


inline void Mesh::OneDim::Positive::order ( )

// run over all segments, order them linearly, check that the mesh is connected
// if it is a loop, set first_ver = Cell::ghost
// if it is an open chain, set first_ver and last_ver
// see paragraph 8.14 in the manual [?]

{	if ( this->first_ver ) return;
	// if first_vertex is not null, the mesh is already ordered
	
	std::list<Cell::Core*>::iterator it0 = this->cells[1].begin();  // will change !
	if ( it0 == this->cells[1].end() )  // empty mesh
	{	assert ( this->cells[1].size() == 0 );  // will change !
		this->first_ver = ( Cell::Negative::Vertex * ) Cell::ghost;
		return;                                                       }
	Cell::Core * seg = *it0;
	assert ( seg );
	assert ( seg->get_dim() == 1 );
	size_t counter = 1;
	Cell::Core * ver, * neg_ver;
	Cell::Core * seg_ini = seg;
	while ( true )
	{	ver = seg->tip();
		assert ( ver );
		assert ( ver->is_positive() );
		neg_ver = ver->reverse_p;
		if ( neg_ver == nullptr )  goto open_chain;
		assert ( not neg_ver->is_positive() );
		std::map<Mesh::Core*,Cell::Core*>::const_iterator
			it = neg_ver->cell_behind_within.find ( this );
		if ( it == neg_ver->cell_behind_within.end() )  goto open_chain;
		seg = it->second;
		assert ( seg );
		if ( seg == seg_ini )  break;
		counter ++;                                              }
	// we assume there is no such loop having only one segment, OK ?
	if ( seg == seg_ini )  // we are dealing with a loop
	{	assert ( counter == this->cells[1].size() );  // may change !
		this->first_ver = ( Cell::Negative::Vertex * ) Cell::ghost;
		return;                                                      }
  open_chain :
	this->last_ver = ( Cell::Positive::Vertex* ) ver;
	seg = seg_ini;
	while ( true )
	{	neg_ver = seg->base();
		assert ( neg_ver );
		assert ( not neg_ver->is_positive() );
		ver = neg_ver->reverse_p;
		assert ( ver->is_positive() );
		std::map<Mesh::Core*,Cell::Core*>::const_iterator
			it = ver->cell_behind_within.find ( this );
		if ( it == ver->cell_behind_within.end() )  break;
		seg = it->second;
		assert ( seg );
		if ( seg == seg_ini )  break;
		counter ++;                                                       }
	assert ( counter == this->cells[1].size() );  // may change !
	this->first_ver = ( Cell::Negative::Vertex * ) neg_ver;                   }
	
//-----------------------------------------------------------------------------//


}  // namespace maniFEM

#endif
// ifndef MANIFEM_MESH_H
