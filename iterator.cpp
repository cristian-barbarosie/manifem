
// iterator.cpp 2021.07.23

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

#include "iterator.h"

using namespace maniFEM;


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base (negative) then tip (positive)
{	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
	  < Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::NormalOrder::AsTheyAre ( seg );  }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::AsTheyAre &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base (negative) then tip (positive)
{	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::NormalOrder::AsTheyAre ( seg );  }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::AsTheyAre &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip (positive) then base (negative)
{	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::AsTheyAre ( seg );  }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::AsTheyAre &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip (positive) then base (negative)
{	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::AsTheyAre ( seg );  }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base then tip (both positive)
{	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
	  < Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::NormalOrder::ForcePositive ( seg );  }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ForcePositive &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base then tip (both positive)
{	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::NormalOrder::ForcePositive ( seg );  }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ForcePositive &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip then base (both positive)
{	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::ForcePositive ( seg );  }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ForcePositive &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip then base (both positive)
{	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::ForcePositive ( seg );  }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base (positive) then tip (negative)
// do not bother whether reverse cells exist or not
{	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg
	             ::NormalOrder::ReverseEachCell::AssumeCellsExist ( seg );   }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base (positive) then tip (negative)
// do not bother whether reverse cells exist or not
{	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg
	             ::NormalOrder::ReverseEachCell::AssumeCellsExist ( seg );   }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip (negative) then base (positive)
// do not bother whether reverse cells exist or not
{	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg
	             ::ReverseOrder::ReverseEachCell::AssumeCellsExist ( seg );  }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip (negative) then base (positive)
// do not bother whether reverse cells exist or not
{	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg
	             ::ReverseOrder::ReverseEachCell::AssumeCellsExist ( seg );  }

//--------------------------------------------------------------------------------------------//


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                      }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                              }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                              }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                              }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                              }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                             }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                             }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                             }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                             }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                             }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                             }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                             }

//--------------------------------------------------------------------------------------------//


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d,
  const tag::AsTheyAre &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base (negative) then tip (positive)
{	assert ( d == 0 );
	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::NormalOrder::AsTheyAre ( seg );  }
	

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &          )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base (negative) then tip (positive)
{	assert ( d == 0 );
	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::NormalOrder::AsTheyAre ( seg );  }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &          )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip (positive) then base (negative)
{	assert ( d == 0 );
	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::AsTheyAre ( seg );  }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &          )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip (positive) then base (negative)
{	assert ( d == 0 );
	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::AsTheyAre ( seg );  }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d,
  const tag::ForcePositive &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base then tip (both positive)
{	assert ( d == 0 );
	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::NormalOrder::ForcePositive ( seg );  }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &              )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base then tip (both positive)
{	assert ( d == 0 );
	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::NormalOrder::ForcePositive ( seg );  }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &              )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip then base (both positive)
{	assert ( d == 0 );
	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::ForcePositive ( seg );  }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &         )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip then base (both positive)
{	assert ( d == 0 );
	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::ForcePositive ( seg );  }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive &                 )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base (positive) then tip (negative)
// do not bother whether reverse cells exist or not
{	assert ( d == 0 );
	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg
	             ::NormalOrder::ReverseEachCell::AssumeCellsExist ( seg );   }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base (positive) then tip (negative)
// do not bother whether reverse cells exist or not
{	assert ( d == 0 );
	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg
	             ::NormalOrder::ReverseEachCell::AssumeCellsExist ( seg );    }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip (negative) then base (positive)
// do not bother whether reverse cells exist or not
{	assert ( d == 0 );
	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg
	             ::ReverseOrder::ReverseEachCell::AssumeCellsExist ( seg );  }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip (negative) then base (positive)
// do not bother whether reverse cells exist or not
{	assert ( d == 0 );
	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg
	             ::ReverseOrder::ReverseEachCell::AssumeCellsExist ( seg );  }

//--------------------------------------------------------------------------------------------//


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base (negative) then tip (positive)
{	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::NormalOrder::AsTheyAre ( seg );  }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base (negative) then tip (positive)
{	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::NormalOrder::AsTheyAre ( seg );  }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip (positive) then base (negative)
{	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::AsTheyAre ( seg );  }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip (positive) then base (negative)
{	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::AsTheyAre ( seg );  }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base then tip (both positive)
{	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::NormalOrder::ForcePositive ( seg );  }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip then base (both positive)
{	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::NormalOrder::ForcePositive ( seg );  }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip then base (both positive)
{	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::ForcePositive ( seg );  }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip then base (both positive)
{	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::ForcePositive ( seg );  }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base (positive) then tip (negative)
// do not bother whether reverse cells exist or not
{	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg
	             ::NormalOrder::ReverseEachCell::AssumeCellsExist ( seg );    }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base (positive) then tip (negative)
// do not bother whether reverse cells exist or not
{	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg
	             ::NormalOrder::ReverseEachCell::AssumeCellsExist ( seg );    }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip (negative) then base (positive)
// do not bother whether reverse cells exist or not
{	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg
	             ::ReverseOrder::ReverseEachCell::AssumeCellsExist ( seg );   }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip (negative) then base (positive)
// do not bother whether reverse cells exist or not
{	assert ( this->cell_enclosed );
	Cell::PositiveSegment * seg = tag::Util::assert_cast
		< Cell::Positive*, Cell::PositiveSegment* > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg
	             ::ReverseOrder::ReverseEachCell::AssumeCellsExist ( seg );   }

//--------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------//


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &              )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &              )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &              )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
  const tag::ThisMeshIsPositive &                                                  )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
  const tag::ThisMeshIsPositive &                                                  )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core *, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core *, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                  )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                  )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core *, const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core *, const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &,  const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core *, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core *, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
  const tag::ThisMeshIsPositive &                                                      )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
  const tag::ThisMeshIsPositive &                                                      )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                  )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                  )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core *, const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core *, const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &,  const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core *, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core *, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
  const tag::ThisMeshIsPositive &                                                      )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
  const tag::ThisMeshIsPositive &                                                      )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }

//--------------------------------------------------------------------------------------------//


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ThisMeshIsPositive &                                     )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ThisMeshIsPositive &                                     )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &          )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &          )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &          )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &          )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &     )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &     )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &              )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &              )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &               )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &              )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &         )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &         )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                 )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &           )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &           )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }

//--------------------------------------------------------------------------------------------//


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ThisMeshIsPositive &                                     )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ThisMeshIsPositive &                                     )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &          )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &          )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &          )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &          )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &     )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &     )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &              )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &              )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &               )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &              )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &         )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &         )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                 )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &           )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &           )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a zero-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }

//--------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------//


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder ( this );  }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::AsTheyAre &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder ( this );  }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::AsTheyAre &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::VerticesOfConnectedOneDimMesh::ReverseOrder ( this );  }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::AsTheyAre &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::VerticesOfConnectedOneDimMesh::ReverseOrder ( this );  }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & )
// in one-dimensional meshes, vertices are positive anyway
{	return new CellIterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder ( this );  }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ForcePositive &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
// in one-dimensional meshes, vertices are positive anyway
{	return new CellIterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder ( this );  }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ForcePositive &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
// in one-dimensional meshes, vertices are positive anyway
{	return new CellIterator::Over::VerticesOfConnectedOneDimMesh::ReverseOrder ( this );  }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ForcePositive &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
// in one-dimensional meshes, vertices are positive anyway
{	return new CellIterator::Over::VerticesOfConnectedOneDimMesh::ReverseOrder ( this );  }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to require reversed vertices for a 1D mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                         )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to require reversed vertices for a 1D mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                         )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to require reversed vertices for a 1D mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &                   )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to require reversed vertices for a 1D mesh." << std::endl;
	exit ( 1 );                                                                                   }

//--------------------------------------------------------------------------------------------//


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::
		SegmentsOfConnectedOneDimMesh::NormalOrder::AsTheyAre ( this );  }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::
		SegmentsOfConnectedOneDimMesh::NormalOrder::AsTheyAre ( this );  }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::
		SegmentsOfConnectedOneDimMesh::ReverseOrder::AsTheyAre ( this );  }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::
		SegmentsOfConnectedOneDimMesh::ReverseOrder::AsTheyAre ( this );  }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::
		SegmentsOfConnectedOneDimMesh::NormalOrder::ForcePositive ( this );  }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::
		SegmentsOfConnectedOneDimMesh::NormalOrder::ForcePositive ( this );  }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::
		SegmentsOfConnectedOneDimMesh::ReverseOrder::ForcePositive ( this );  }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::
		SegmentsOfConnectedOneDimMesh::ReverseOrder::ForcePositive ( this );  }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::SegmentsOfConnectedOneDimMesh::
		NormalOrder::ReverseEachCell::AssumeCellsExist ( this );      }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                         )
{	return new CellIterator::Over::SegmentsOfConnectedOneDimMesh::
		NormalOrder::ReverseEachCell::AssumeCellsExist ( this );      }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                         )
{	return new CellIterator::Over::SegmentsOfConnectedOneDimMesh::
		ReverseOrder::ReverseEachCell::AssumeCellsExist ( this );      }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &                   )
{	return new CellIterator::Over::SegmentsOfConnectedOneDimMesh::
		ReverseOrder::ReverseEachCell::AssumeCellsExist ( this );      }

//--------------------------------------------------------------------------------------------//


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d,
  const tag::AsTheyAre &, const tag::ThisMeshIsPositive & )
{	if ( d == 0 )
		return new CellIterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder ( this );
	if ( d == 1 )
		return new CellIterator::Over::
			SegmentsOfConnectedOneDimMesh::NormalOrder::AsTheyAre ( this );
	assert ( false );                                                                      }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	if ( d == 0 )
		return new CellIterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder ( this );
	if ( d == 1 )
		return new CellIterator::Over::
			SegmentsOfConnectedOneDimMesh::NormalOrder::AsTheyAre ( this );
	assert ( false );                                                                      }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &           )
{	if ( d == 0 )
		return new CellIterator::Over::VerticesOfConnectedOneDimMesh::ReverseOrder ( this );
	if ( d == 1 )
		return new CellIterator::Over::
			SegmentsOfConnectedOneDimMesh::ReverseOrder::AsTheyAre ( this );
	assert ( false );                                                                      }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &      )
{	if ( d == 0 )
		return new CellIterator::Over::VerticesOfConnectedOneDimMesh::ReverseOrder ( this );
	if ( d == 1 )
		return new CellIterator::Over::
			SegmentsOfConnectedOneDimMesh::ReverseOrder::AsTheyAre ( this );
	assert ( false );                                                                      }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d,
  const tag::ForcePositive &, const tag::ThisMeshIsPositive & )
{	if ( d == 0 )  // in one-dimensional meshes, vertices are positive anyway
		return new CellIterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder ( this );
	if ( d == 1 )
		return new CellIterator::Over::SegmentsOfConnectedOneDimMesh::
			NormalOrder::ForcePositive ( this );
	assert ( false );                                                                      }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &              )
{	if ( d == 0 )  // in one-dimensional meshes, vertices are positive anyway
		return new CellIterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder ( this );
	if ( d == 1 )
		return new CellIterator::Over::SegmentsOfConnectedOneDimMesh::
			NormalOrder::ForcePositive ( this );
	assert ( false );                                                                      }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &              )
{	if ( d == 0 )  // in one-dimensional meshes, vertices are positive anyway
		return new CellIterator::Over::VerticesOfConnectedOneDimMesh::ReverseOrder ( this );
	if ( d == 1 )
		return new CellIterator::Over::SegmentsOfConnectedOneDimMesh::
			ReverseOrder::ForcePositive ( this );
	assert ( false );                                                                      }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &          )
{	if ( d == 0 )  // in one-dimensional meshes, vertices are positive anyway
		return new CellIterator::Over::VerticesOfConnectedOneDimMesh::ReverseOrder ( this );
	if ( d == 1 )
		return new CellIterator::Over::SegmentsOfConnectedOneDimMesh::
			ReverseOrder::ForcePositive ( this );
	assert ( false );                                                                      }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive &                    )
{	if ( d == 0 )  // in one-dimensional meshes, vertices are positive anyway
	{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
							<< __extension__ __PRETTY_FUNCTION__ << ": ";
		std::cout << "It makes no sense to require reversed vertices for a 1D mesh." << std::endl;
		exit ( 1 );                                                                                }
	if ( d == 1 )
		return new CellIterator::Over::SegmentsOfConnectedOneDimMesh::
			NormalOrder::ReverseEachCell::AssumeCellsExist ( this );
	assert ( false );                                                                              }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	if ( d == 0 )  // in one-dimensional meshes, vertices are positive anyway
	{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
							<< __extension__ __PRETTY_FUNCTION__ << ": ";
		std::cout << "It makes no sense to require reversed vertices for a 1D mesh." << std::endl;
		exit ( 1 );                                                                                }
	if ( d == 1 )
		return new CellIterator::Over::SegmentsOfConnectedOneDimMesh::
			NormalOrder::ReverseEachCell::AssumeCellsExist ( this );
	assert ( false );                                                                              }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	if ( d == 0 )  // in one-dimensional meshes, vertices are positive anyway
	{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
							<< __extension__ __PRETTY_FUNCTION__ << ": ";
		std::cout << "It makes no sense to require reversed vertices for a 1D mesh." << std::endl;
		exit ( 1 );                                                                                }
	if ( d == 1 )
		return new CellIterator::Over::SegmentsOfConnectedOneDimMesh::
			ReverseOrder::ReverseEachCell::AssumeCellsExist ( this );
	assert ( false );                                                                              }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
{	if ( d == 0 )  // in one-dimensional meshes, vertices are positive anyway
	{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
							<< __extension__ __PRETTY_FUNCTION__ << ": ";
		std::cout << "It makes no sense to require reversed vertices for a 1D mesh." << std::endl;
		exit ( 1 );                                                                                }
	if ( d == 1 )
		return new CellIterator::Over::SegmentsOfConnectedOneDimMesh::
			ReverseOrder::ReverseEachCell::AssumeCellsExist ( this );
	assert ( false );                                                                              }

//--------------------------------------------------------------------------------------------//


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::
		SegmentsOfConnectedOneDimMesh::NormalOrder::AsTheyAre ( this );  }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::
		SegmentsOfConnectedOneDimMesh::NormalOrder::AsTheyAre ( this );  }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::
		SegmentsOfConnectedOneDimMesh::ReverseOrder::AsTheyAre ( this );  }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::
		SegmentsOfConnectedOneDimMesh::ReverseOrder::AsTheyAre ( this );  }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::
		SegmentsOfConnectedOneDimMesh::NormalOrder::ForcePositive ( this );  }
	

CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::
		SegmentsOfConnectedOneDimMesh::NormalOrder::ForcePositive ( this );  }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::
		SegmentsOfConnectedOneDimMesh::ReverseOrder::ForcePositive ( this );  }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::
		SegmentsOfConnectedOneDimMesh::ReverseOrder::ForcePositive ( this );  }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive &       )
{	return new CellIterator::Over::SegmentsOfConnectedOneDimMesh::
		NormalOrder::ReverseEachCell::AssumeCellsExist ( this );      }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                             )
{	return new CellIterator::Over::SegmentsOfConnectedOneDimMesh::
		NormalOrder::ReverseEachCell::AssumeCellsExist ( this );      }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                             )
{	return new CellIterator::Over::SegmentsOfConnectedOneDimMesh::
		ReverseOrder::ReverseEachCell::AssumeCellsExist ( this );      }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &                        )
{	return new CellIterator::Over::SegmentsOfConnectedOneDimMesh::
		ReverseOrder::ReverseEachCell::AssumeCellsExist ( this );      }

//-----------------------------------------------------------------------------------------


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &              )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &              )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &              )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
  const tag::ThisMeshIsPositive &                                                  )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core *, const tag::RequireOrder &,
  const tag::ThisMeshIsPositive &                                                  )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core *, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core *, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                  )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                  )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core *, const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core *, const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &,  const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core *, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core *, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
  const tag::ThisMeshIsPositive &                                                      )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
  const tag::ThisMeshIsPositive &                                                      )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                  )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive &                  )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core *, const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core *, const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &,  const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core *, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core *, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
  const tag::ThisMeshIsPositive &                                                      )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core *, const tag::ReverseOrderIfAny &,
  const tag::ThisMeshIsPositive &                                                      )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }

//--------------------------------------------------------------------------------------------//


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ThisMeshIsPositive &                                     )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ThisMeshIsPositive &                                     )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &          )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &          )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &          )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &          )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &     )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &     )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &              )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &              )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &               )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &              )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &         )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &         )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                 )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &           )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &           )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }

//--------------------------------------------------------------------------------------------//


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ThisMeshIsPositive &                                     )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ThisMeshIsPositive &                                     )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &          )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &          )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &          )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &          )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &     )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &     )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &              )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &              )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &               )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &              )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &         )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &         )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                 )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &           )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core *,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &           )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot iterate around a cell in a one-dimensional mesh." << std::endl;
	exit ( 1 );                                                                           }

//--------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------//


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::CellsOfFuzzyMesh::AsTheyAre ( this, tag::cells_of_dim, 0 );  }
// perhaps provide a different constructor for calling directly msh->cells[0]
// or msh->cells.front(), it should be slightly faster when compiled with optimization


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::AsTheyAre &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to require order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::AsTheyAre &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to reverse order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::AsTheyAre &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::CellsOfFuzzyMesh::AsTheyAre ( this, tag::cells_of_dim, 0 );  }
// perhaps provide a different constructor for calling directly msh->cells[0]
// or msh->cells.front(), it should be slightly faster when compiled with optimization


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & )
// fuzzy meshes have dimension at least one; thus, vertices are positive anyway
{	return new CellIterator::Over::CellsOfFuzzyMesh::AsTheyAre ( this, tag::cells_of_dim, 0 );  }
// perhaps provide a different constructor for calling directly msh->cells[0]
// or msh->cells.front(), it should be slightly faster when compiled with optimization


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ForcePositive &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to require order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ForcePositive &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to reverse order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ForcePositive &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
// fuzzy meshes have dimension at least one; thus, vertices are positive anyway
{	return new CellIterator::Over::CellsOfFuzzyMesh::AsTheyAre ( this, tag::cells_of_dim, 0 );  }
// perhaps provide a different constructor for calling directly msh->cells[0]
// or msh->cells.front(), it should be slightly faster when compiled with optimization


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to require reversed vertices for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                                   }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                         )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "This makes no sense for two different reasons." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                         )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "This makes no sense for two different reasons." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &                   )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to require reversed vertices for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                                   }

//--------------------------------------------------------------------------------------------//


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::CellsOfFuzzyMesh::AsTheyAre ( this, tag::cells_of_dim, 1 );  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to require order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to reverse order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::CellsOfFuzzyMesh::AsTheyAre ( this, tag::cells_of_dim, 1 );  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & )
{	if ( this->get_dim_plus_one() == 2 )
		return new CellIterator::Over::CellsOfFuzzyMesh::ForcePositive ( this, tag::cells_of_dim, 1 );
	// else : dim >= 2, segments are positive anyway
	assert ( this->get_dim_plus_one() > 2 );
	return new CellIterator::Over::CellsOfFuzzyMesh::AsTheyAre ( this, tag::cells_of_dim, 1 );       }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to require order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to reverse order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
{	if ( this->get_dim_plus_one() == 2 )
		return new CellIterator::Over::CellsOfFuzzyMesh::ForcePositive
			( this, tag::cells_of_dim, 1 );
	// else : dim >= 2, segments are positive anyway
	assert ( this->get_dim_plus_one() > 2 );
	return new CellIterator::Over::CellsOfFuzzyMesh::AsTheyAre ( this, tag::cells_of_dim, 1 );  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive & )
{	assert ( this->get_dim_plus_one() == 2 );
	return new CellIterator::Over::CellsOfFuzzyMesh::ReverseEachCell::AssumeCellsExist
		( this, tag::cells_of_dim, 1 );                                                  }
// perhaps provide a different constructor for calling directly msh->cells.back(),
// it should be slightly faster when compiled with optimization


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                         )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to require order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                         )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to reverse order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &                   )
{	assert ( this->get_dim_plus_one() == 2 );
	return new CellIterator::Over::CellsOfFuzzyMesh::ReverseEachCell::AssumeCellsExist
		( this, tag::cells_of_dim, 1 );                                                  }
// perhaps provide a different constructor for calling directly msh->cells.back(),
// it should be slightly faster when compiled with optimization

//--------------------------------------------------------------------------------------------//


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d,
  const tag::AsTheyAre &, const tag::ThisMeshIsPositive & )
{	assert ( this->get_dim_plus_one() > d );
	return new CellIterator::Over::CellsOfFuzzyMesh::AsTheyAre ( this, tag::cells_of_dim, d );  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to require order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to reverse order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
{	assert ( this->get_dim_plus_one() > d );
	return new CellIterator::Over::CellsOfFuzzyMesh::AsTheyAre ( this, tag::cells_of_dim, d );  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d,
  const tag::ForcePositive &, const tag::ThisMeshIsPositive & )
{	assert ( this->get_dim_plus_one() > d );
	if ( this->get_dim_plus_one() == d+1 )
		return new CellIterator::Over::CellsOfFuzzyMesh::ForcePositive
			( this, tag::cells_of_dim, d );
	// else : cell dim < mesh dim, cells are positive anyway
	return new CellIterator::Over::CellsOfFuzzyMesh::AsTheyAre ( this, tag::cells_of_dim, d );  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &              )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to require order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &              )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to reverse order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &          )
{	assert ( this->get_dim_plus_one() > d );
	if ( this->get_dim_plus_one() == d+1 )
		return new CellIterator::Over::CellsOfFuzzyMesh::ForcePositive
			( this, tag::cells_of_dim, d );
	// else : cell dim < mesh dim, cells are positive anyway
	return new CellIterator::Over::CellsOfFuzzyMesh::AsTheyAre ( this, tag::cells_of_dim, d );  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive &                    )
{	assert ( this->get_dim_plus_one() == d+1 );
	return new CellIterator::Over::CellsOfFuzzyMesh::ReverseEachCell::AssumeCellsExist
		( this, tag::cells_of_dim, d );                                                  }
// perhaps provide a different constructor for calling directly msh->cells.back(),
// it should be slightly faster when compiled with optimization


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to require order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to reverse order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d,
  const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
{	assert ( this->get_dim_plus_one() == d+1 );
	return new CellIterator::Over::CellsOfFuzzyMesh::ReverseEachCell::AssumeCellsExist
		( this, tag::cells_of_dim, d );                                                  }
// perhaps provide a different constructor for calling directly msh->cells.back(),
// it should be slightly faster when compiled with optimization

//--------------------------------------------------------------------------------------------//


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & )
{ return new CellIterator::Over::CellsOfFuzzyMesh::AsTheyAre
		( this, tag::cells_of_dim, tag::Util::assert_diff ( this->get_dim_plus_one(), 1 ) );  }
// perhaps provide a different constructor for calling directly msh->cells.back(),
// it should be slightly faster when compiled with optimization


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to require order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to reverse order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::CellsOfFuzzyMesh::AsTheyAre
		( this, tag::cells_of_dim, tag::Util::assert_diff ( this->get_dim_plus_one(), 1 ) );  }
// perhaps provide a different constructor for calling directly msh->cells.back(),
// it should be slightly faster when compiled with optimization


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & )
{ return new CellIterator::Over::CellsOfFuzzyMesh::ForcePositive
		( this, tag::cells_of_dim, tag::Util::assert_diff ( this->get_dim_plus_one(), 1 ) );  }
// perhaps provide a different constructor for calling directly msh->cells.back(),
// it should be slightly faster when compiled with optimization
	

CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to require order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to reverse order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
{ return new CellIterator::Over::CellsOfFuzzyMesh::ForcePositive
		( this, tag::cells_of_dim, tag::Util::assert_diff ( this->get_dim_plus_one(), 1 ) );  }
// perhaps provide a different constructor for calling directly msh->cells.back(),
// it should be slightly faster when compiled with optimization


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::ThisMeshIsPositive &       )
{	return new CellIterator::Over::CellsOfFuzzyMesh::ReverseEachCell::AssumeCellsExist
		( this, tag::cells_of_dim, tag::Util::assert_diff ( this->get_dim_plus_one(), 1 ) );  }
// perhaps provide a different constructor for calling directly msh->cells.back(),
// it should be slightly faster when compiled with optimization


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                             )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to require order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                             )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to reverse order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &                        )
{	return new CellIterator::Over::CellsOfFuzzyMesh::ReverseEachCell::AssumeCellsExist
		( this, tag::cells_of_dim, tag::Util::assert_diff ( this->get_dim_plus_one(), 1 ) );  }
// perhaps provide a different constructor for calling directly msh->cells.back(),
// it should be slightly faster when compiled with optimization

//-----------------------------------------------------------------------------------------



CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core * c, const tag::ThisMeshIsPositive &              )
	
{	assert ( c->get_dim() == 0 );
	assert ( this->get_dim_plus_one() == 3 );  // 2D mesh
	return new CellIterator::Around::OneCell::OfCodimTwo
		::WorkAround2D::NormalOrder::BuildReverseCells ( this,
			tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core * c, const tag::ThisMeshIsPositive &              )
// in the future, change to CellIterator::Around::OneVertex::OfAnyCodim::OverSegments
	
{	assert ( c->get_dim() == 0 );
	assert ( this->get_dim_plus_one() == 3 );  // 2D mesh
	if ( c->is_positive() )
		return new CellIterator::Around::OneCell::OfCodimTwo
			::OverVertices::NormalOrder::BuildReverseCells::AsTheyAre ( this,
			  tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );
	// else : c is negative
	return new CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::NormalOrder::BuildReverseCells::ReverseEachCell ( this,
		  tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const>
	 	    ( c->reverse_attr.core )                                            );      }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core * c, const tag::ThisMeshIsPositive &              )
// in the future, change to CellIterator::Around::OneVertex::OfAnyCodim::OverSegments
	
{	assert ( c->get_dim() == 0 );
	assert ( this->get_dim_plus_one() == 3 );  // 2D mesh
	if ( c->is_positive() )
		return new CellIterator::Around::OneCell::OfCodimTwo
			::OverVertices::NormalOrder::AssumeCellsExist::AsTheyAre ( this,
			  tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );
	// else : c is negative
	return new CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::NormalOrder::AssumeCellsExist::ReverseEachCell ( this,
		  tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const>
	 	    ( c->reverse_attr.core )                                            );      }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core * c, const tag::RequireOrder &,
  const tag::ThisMeshIsPositive &                                                  )

{	assert ( c->get_dim() == 0 );
	assert ( this->get_dim_plus_one() == 3 );  // 2D mesh
	if ( c->is_positive() )
		return new CellIterator::Around::OneCell::OfCodimTwo
			::OverVertices::NormalOrder::BuildReverseCells::AsTheyAre ( this,
			  tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );
	// else : c is negative
	return new CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::NormalOrder::BuildReverseCells::ReverseEachCell ( this,
		  tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const>
	 	    ( c->reverse_attr.core )                                            );      }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core * c, const tag::RequireOrder &,
  const tag::ThisMeshIsPositive &                                                  )

{	assert ( c->get_dim() == 0 );
	assert ( this->get_dim_plus_one() == 3 );  // 2D mesh
	if ( c->is_positive() )
		return new CellIterator::Around::OneCell::OfCodimTwo
			::OverVertices::NormalOrder::AssumeCellsExist::AsTheyAre ( this,
			  tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );
	// else : c is negative
	return new CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::NormalOrder::AssumeCellsExist::ReverseEachCell ( this,
		  tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const>
	 	    ( c->reverse_attr.core )                                            );      }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core * c, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )

{	assert ( c->get_dim() == 0 );
	assert ( this->get_dim_plus_one() == 3 );  // 2D mesh
	if ( c->is_positive() )
		return new CellIterator::Around::OneCell::OfCodimTwo
			::OverVertices::ReverseOrder::BuildReverseCells::AsTheyAre ( this,
			  tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );
	// else : c is negative
	return new CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::ReverseOrder::BuildReverseCells::ReverseEachCell ( this,
		  tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const>
	 	    ( c->reverse_attr.core )                                            );      }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core * c, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )

{	assert ( c->get_dim() == 0 );
	assert ( this->get_dim_plus_one() == 3 );  // 2D mesh
	if ( c->is_positive() )
		return new CellIterator::Around::OneCell::OfCodimTwo
			::OverVertices::ReverseOrder::AssumeCellsExist::AsTheyAre ( this,
			  tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );
	// else : c is negative
	return new CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::ReverseOrder::AssumeCellsExist::ReverseEachCell ( this,
		  tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const>
	 	    ( c->reverse_attr.core )                                            );      }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core * c, const tag::ReverseOrderIfAny &,
  const tag::ThisMeshIsPositive &                                                 )
// in the future, change to CellIterator::Around::OneVertex::OfAnyCodim::OverSegments
	
{	assert ( c->get_dim() == 0 );
	assert ( this->get_dim_plus_one() == 3 );  // 2D mesh
	if ( c->is_positive() )
		return new CellIterator::Around::OneCell::OfCodimTwo
			::OverVertices::ReverseOrder::BuildReverseCells::AsTheyAre ( this,
			  tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );
	// else : c is negative
	return new CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::ReverseOrder::BuildReverseCells::ReverseEachCell ( this,
		  tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const>
	 	    ( c->reverse_attr.core )                                            );      }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core * c, const tag::ReverseOrderIfAny &,
  const tag::ThisMeshIsPositive &                                                 )
// in the future, change to CellIterator::Around::OneVertex::OfAnyCodim::OverSegments
	
{	assert ( c->get_dim() == 0 );
	assert ( this->get_dim_plus_one() == 3 );  // 2D mesh
	if ( c->is_positive() )
		return new CellIterator::Around::OneCell::OfCodimTwo
			::OverVertices::ReverseOrder::AssumeCellsExist::AsTheyAre ( this,
			  tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );
	// else : c is negative
	return new CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::ReverseOrder::AssumeCellsExist::ReverseEachCell ( this,
		  tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const>
	 	    ( c->reverse_attr.core )                                            );      }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core * c, const tag::ThisMeshIsPositive &                  )
// in the future, change to CellIterator::Around::OneVertex::OfAnyCodim::OverSegments
	
{	assert ( c->get_dim() == 0 );
	assert ( this->get_dim_plus_one() == 3 );  // 2D mesh
	return new CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::NormalOrder::BuildReverseCells::ForcePositive ( this,
		  tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const>
	       ( c->get_positive().core )                                        );  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core * c, const tag::ThisMeshIsPositive &                  )
// in the future, change to CellIterator::Around::OneVertex::OfAnyCodim::OverSegments
	
{	assert ( c->get_dim() == 0 );
	assert ( this->get_dim_plus_one() == 3 );  // 2D mesh
	return new CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::NormalOrder::AssumeCellsExist::ForcePositive ( this,
		  tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const>
	       ( c->get_positive().core )                                        );  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core * c, const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
	
{	assert ( c->get_dim() == 0 );
	assert ( this->get_dim_plus_one() == 3 );  // 2D mesh
	return new CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::NormalOrder::BuildReverseCells::ForcePositive ( this,
		  tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const>
	       ( c->get_positive().core )                                        );  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core * c, const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
	
{	assert ( c->get_dim() == 0 );
	assert ( this->get_dim_plus_one() == 3 );  // 2D mesh
	return new CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::NormalOrder::AssumeCellsExist::ForcePositive ( this,
		  tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const>
	       ( c->get_positive().core )                                        );  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &,  const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core * c, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )

{	assert ( c->get_dim() == 0 );
	assert ( this->get_dim_plus_one() == 3 );  // 2D mesh
	return new CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::ReverseOrder::BuildReverseCells::ForcePositive ( this,
		  tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const>
	       ( c->get_positive().core )                                        );  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core * c, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )

{	assert ( c->get_dim() == 0 );
	assert ( this->get_dim_plus_one() == 3 );  // 2D mesh
	return new CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::ReverseOrder::AssumeCellsExist::ForcePositive ( this,
		  tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const>
	       ( c->get_positive().core )                                        );  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core * c, const tag::ReverseOrderIfAny &,
  const tag::ThisMeshIsPositive &                                                      )
// in the future, change to CellIterator::Around::OneVertex::OfAnyCodim::OverSegments
	
{	assert ( c->get_dim() == 0 );
	assert ( this->get_dim_plus_one() == 3 );  // 2D mesh
	return new CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::ReverseOrder::BuildReverseCells::ForcePositive ( this,
		  tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const>
	       ( c->get_positive().core )                                        );  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core * c, const tag::ReverseOrderIfAny &,
  const tag::ThisMeshIsPositive &                                                      )
// in the future, change to CellIterator::Around::OneVertex::OfAnyCodim::OverSegments
	
{	assert ( c->get_dim() == 0 );
	assert ( this->get_dim_plus_one() == 3 );  // 2D mesh
	return new CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::ReverseOrder::AssumeCellsExist::ForcePositive ( this,
		  tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const>
	       ( c->get_positive().core )                                        );  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core * c, const tag::ThisMeshIsPositive &                  )
// in the future, change to CellIterator::Around::OneVertex::OfAnyCodim::OverSegments
	
{	assert ( c->get_dim() == 0 );
	assert ( this->get_dim_plus_one() == 3 );  // 2D mesh
	return new CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::NormalOrder::BuildReverseCells::ReverseEachCell ( this,
			tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core * c, const tag::ThisMeshIsPositive &                  )
// in the future, change to CellIterator::Around::OneVertex::OfAnyCodim::OverSegments
	
{	assert ( c->get_dim() == 0 );
	assert ( this->get_dim_plus_one() == 3 );  // 2D mesh
	return new CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::NormalOrder::AssumeCellsExist::ReverseEachCell ( this,
			tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core * c, const tag::RequireOrder &, const tag::ThisMeshIsPositive & )

{	assert ( c->get_dim() == 0 );
	assert ( this->get_dim_plus_one() == 3 );  // 2D mesh
	return new CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::NormalOrder::BuildReverseCells::ReverseEachCell ( this,
			tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core * c, const tag::RequireOrder &, const tag::ThisMeshIsPositive & )

{	assert ( c->get_dim() == 0 );
	assert ( this->get_dim_plus_one() == 3 );  // 2D mesh
	return new CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::NormalOrder::AssumeCellsExist::ReverseEachCell ( this,
			tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &,  const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core * c, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
	
{	assert ( c->get_dim() == 0 );
	assert ( this->get_dim_plus_one() == 3 );  // 2D mesh
	return new CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::ReverseOrder::BuildReverseCells::ReverseEachCell ( this,
			tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core * c, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
	
{	assert ( c->get_dim() == 0 );
	assert ( this->get_dim_plus_one() == 3 );  // 2D mesh
	return new CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::ReverseOrder::AssumeCellsExist::ReverseEachCell ( this,
			tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::BuildCellsIfNec &,
  const tag::Around &, Cell::Core * c, const tag::ReverseOrderIfAny &,
  const tag::ThisMeshIsPositive &                                                      )
// in the future, change to CellIterator::Around::OneVertex::OfAnyCodim::OverSegments
	
{	assert ( c->get_dim() == 0 );
	assert ( this->get_dim_plus_one() == 3 );  // 2D mesh
	return new CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::ReverseOrder::BuildReverseCells::ReverseEachCell ( this,
			tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBuildCells &,
  const tag::Around &, Cell::Core * c, const tag::ReverseOrderIfAny &,
  const tag::ThisMeshIsPositive &                                                      )
// in the future, change to CellIterator::Around::OneVertex::OfAnyCodim::OverSegments
	
{	assert ( c->get_dim() == 0 );
	assert ( this->get_dim_plus_one() == 3 );  // 2D mesh
	return new CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::ReverseOrder::AssumeCellsExist::ReverseEachCell ( this,
			tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );  }

//--------------------------------------------------------------------------------------------//


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core * c,
  const tag::ThisMeshIsPositive &                                     )

{	assert ( c->get_dim() + 2 < this->get_dim_plus_one() );  // co-dimension >= 2
	assert ( c->get_dim() < d );
	assert ( d < this->get_dim_plus_one() );
	// CellIterators OfAnyCodim are not implemented yet, so ...
	assert ( c->get_dim() + 3 == this->get_dim_plus_one() );  // co-dimension == 2
	if ( d == c->get_dim() + 1 )
		return new CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::NormalOrder::BuildReverseCells::AsTheyAre ( this,
			tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );
	else
	{	assert ( d == c->get_dim() + 2 );
		return new CellIterator::Around::OneCell::OfCodimTwo
        ::OverSegments::NormalOrder::AsTheyAre ( this,
        tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );  }  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core * c,
  const tag::ThisMeshIsPositive &                                     )

{	assert ( c->get_dim() + 2 < this->get_dim_plus_one() );  // co-dimension >= 2
	assert ( c->get_dim() < d );
	assert ( d < this->get_dim_plus_one() );
	// CellIterators OfAnyCodim are not implemented yet, so ...
	assert ( c->get_dim() + 3 == this->get_dim_plus_one() );  // co-dimension == 2
	if ( d == c->get_dim() + 1 )
		return new CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::NormalOrder::AssumeCellsExist::AsTheyAre ( this,
			tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );
	else
	{	assert ( d == c->get_dim() + 2 );
		return new CellIterator::Around::OneCell::OfCodimTwo
		    ::OverSegments::NormalOrder::AsTheyAre ( this,
		  tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );  }  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core * c,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &          )

{	assert ( c->get_dim() + 3 == this->get_dim_plus_one() );  // co-dimension == 2
	assert ( c->get_dim() < d );
	assert ( d < this->get_dim_plus_one() );
	if ( d == c->get_dim() + 1 )
		return new CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::NormalOrder::BuildReverseCells::AsTheyAre ( this,
			tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );
	else
	{	assert ( d == c->get_dim() + 2 );
		return new CellIterator::Around::OneCell::OfCodimTwo
		     ::OverSegments::NormalOrder::AsTheyAre ( this,
      tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );  }  }
 

CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core * c,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &          )

{	assert ( c->get_dim() + 3 == this->get_dim_plus_one() );  // co-dimension == 2
	assert ( c->get_dim() < d );
	assert ( d < this->get_dim_plus_one() );
	if ( d == c->get_dim() + 1 )
		return new CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::NormalOrder::AssumeCellsExist::AsTheyAre ( this,
			tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );
	else
	{	assert ( d == c->get_dim() + 2 );
		return new CellIterator::Around::OneCell::OfCodimTwo
		     ::OverSegments::NormalOrder::AsTheyAre ( this,
      tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );  }  }
 

CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core * c,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &          )

{	assert ( c->get_dim() + 3 == this->get_dim_plus_one() );  // co-dimension == 2
	assert ( c->get_dim() < d );
	assert ( d < this->get_dim_plus_one() );
	if ( d == c->get_dim() + 1 )
		return new CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::ReverseOrder::BuildReverseCells::AsTheyAre ( this,
			tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );
	else
	{	assert ( d == c->get_dim() + 2 );
		return new CellIterator::Around::OneCell::OfCodimTwo
			::OverSegments::ReverseOrder::AsTheyAre ( this,
        tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );  }  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core * c,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &          )

{	assert ( c->get_dim() + 3 == this->get_dim_plus_one() );  // co-dimension == 2
	assert ( c->get_dim() < d );
	assert ( d < this->get_dim_plus_one() );
	if ( d == c->get_dim() + 1 )
		return new CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::ReverseOrder::AssumeCellsExist::AsTheyAre ( this,
			tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );
	else
	{	assert ( d == c->get_dim() + 2 );
		return new CellIterator::Around::OneCell::OfCodimTwo
			::OverSegments::ReverseOrder::AsTheyAre ( this,
        tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );  }  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core * c,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &     )

{	assert ( c->get_dim() + 2 < this->get_dim_plus_one() );  // co-dimension >= 2
	assert ( c->get_dim() < d );
	assert ( d < this->get_dim_plus_one() );
	// CellIterators OfAnyCodim are not implemented yet, so ...
	assert ( c->get_dim() + 3 == this->get_dim_plus_one() );  // co-dimension == 2
	if ( d == c->get_dim() + 1 )
		return new CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::ReverseOrder::BuildReverseCells::AsTheyAre ( this,
			tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );
	else
	{	assert ( d == c->get_dim() + 2 );
		return new CellIterator::Around::OneCell::OfCodimTwo
        ::OverSegments::ReverseOrder::AsTheyAre ( this,
        tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );  }  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core * c,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &     )

{	assert ( c->get_dim() + 2 < this->get_dim_plus_one() );  // co-dimension >= 2
	assert ( c->get_dim() < d );
	assert ( d < this->get_dim_plus_one() );
	// CellIterators OfAnyCodim are not implemented yet, so ...
	assert ( c->get_dim() + 3 == this->get_dim_plus_one() );  // co-dimension == 2
	if ( d == c->get_dim() + 1 )
		return new CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::ReverseOrder::AssumeCellsExist::AsTheyAre ( this,
			tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );
	else
	{	assert ( d == c->get_dim() + 2 );
		return new CellIterator::Around::OneCell::OfCodimTwo
        ::OverSegments::ReverseOrder::AsTheyAre ( this,
        tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );  }  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core * c,
  const tag::ThisMeshIsPositive &                                         )

{	assert ( c->get_dim() + 2 < this->get_dim_plus_one() );  // co-dimension >= 2
	assert ( c->get_dim() < d );
	assert ( d < this->get_dim_plus_one() );
	// CellIterators OfAnyCodim are not implemented yet, so ...
	assert ( c->get_dim() + 3 == this->get_dim_plus_one() );  // co-dimension == 2
	if ( d == c->get_dim() + 1 )
		return new CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::NormalOrder::BuildReverseCells::ForcePositive ( this,
			tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );
	else
	{	assert ( d == c->get_dim() + 2 );
		return new CellIterator::Around::OneCell::OfCodimTwo
		//    ::OverSegments::NormalOrder::ForcePositive ( this,
        ::OverSegments::NormalOrder::AsTheyAre ( this,
        tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );  }  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core * c,
  const tag::ThisMeshIsPositive &                                         )

{	assert ( c->get_dim() + 2 < this->get_dim_plus_one() );  // co-dimension >= 2
	assert ( c->get_dim() < d );
	assert ( d < this->get_dim_plus_one() );
	// CellIterators OfAnyCodim are not implemented yet, so ...
	assert ( c->get_dim() + 3 == this->get_dim_plus_one() );  // co-dimension == 2
	if ( d == c->get_dim() + 1 )
		return new CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::NormalOrder::AssumeCellsExist::ForcePositive ( this,
			tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );
	else
	{	assert ( d == c->get_dim() + 2 );
		return new CellIterator::Around::OneCell::OfCodimTwo
		//    ::OverSegments::NormalOrder::ForcePositive ( this,
        ::OverSegments::NormalOrder::AsTheyAre ( this,
        tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );  }  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core * c,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &              )

{	assert ( c->get_dim() + 3 == this->get_dim_plus_one() );  // co-dimension == 2
	assert ( c->get_dim() < d );
	assert ( d < this->get_dim_plus_one() );
	if ( d == c->get_dim() + 1 )
		return new CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::NormalOrder::BuildReverseCells::ForcePositive ( this,
			tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );
	else
	{	assert ( d == c->get_dim() + 2 );
		return new CellIterator::Around::OneCell::OfCodimTwo
		//    ::OverSegments::NormalOrder::ForcePositive ( this,
        ::OverSegments::NormalOrder::AsTheyAre ( this,
        tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );  }  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core * c,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &              )

{	assert ( c->get_dim() + 3 == this->get_dim_plus_one() );  // co-dimension == 2
	assert ( c->get_dim() < d );
	assert ( d < this->get_dim_plus_one() );
	if ( d == c->get_dim() + 1 )
		return new CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::NormalOrder::AssumeCellsExist::ForcePositive ( this,
			tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );
	else
	{	assert ( d == c->get_dim() + 2 );
		return new CellIterator::Around::OneCell::OfCodimTwo
		//    ::OverSegments::NormalOrder::ForcePositive ( this,
        ::OverSegments::NormalOrder::AsTheyAre ( this,
        tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );  }  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core * c,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &               )

{	assert ( c->get_dim() + 3 == this->get_dim_plus_one() );  // co-dimension == 2
	assert ( c->get_dim() < d );
	assert ( d < this->get_dim_plus_one() );
	if ( d == c->get_dim() + 1 )
		return new CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::ReverseOrder::BuildReverseCells::ForcePositive ( this,
			tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );
	else
	{	assert ( d == c->get_dim() + 2 );
		return new CellIterator::Around::OneCell::OfCodimTwo
		//    ::OverSegments::ReverseOrder::ForcePositive ( this,
        ::OverSegments::ReverseOrder::AsTheyAre ( this,
        tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );  }  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core * c,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &              )

{	assert ( c->get_dim() + 3 == this->get_dim_plus_one() );  // co-dimension == 2
	assert ( c->get_dim() < d );
	assert ( d < this->get_dim_plus_one() );
	if ( d == c->get_dim() + 1 )
		return new CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::ReverseOrder::AssumeCellsExist::ForcePositive ( this,
			tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );
	else
	{	assert ( d == c->get_dim() + 2 );
		return new CellIterator::Around::OneCell::OfCodimTwo
		//    ::OverSegments::ReverseOrder::ForcePositive ( this,
        ::OverSegments::ReverseOrder::AsTheyAre ( this,
        tag::Util::assert_cast < Cell::Core * const, Cell::Positive * const> ( c ) );  }  }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core * c,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &         )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                             }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core * c,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &         )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                             }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core * c,
  const tag::ThisMeshIsPositive &                                            )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                            }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core * c, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                             }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core * c,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                             }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core * c,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                             }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core * c,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                 )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                            }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core * c,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                             }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core * c,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &           )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                             }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core * c,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &           )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                             }

//--------------------------------------------------------------------------------------------//


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core * c,
  const tag::ThisMeshIsPositive &                                   )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                             }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core * c,
  const tag::ThisMeshIsPositive &                                   )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                             }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core * c,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &        )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                             }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core * c,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &        )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                             }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core * c,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &        )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                             }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core * c,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &        )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                             }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core * c,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &   )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                             }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core * c,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &   )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                            }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core * c,
  const tag::ThisMeshIsPositive &                                   )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                            }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core * c,
  const tag::ThisMeshIsPositive &                                   )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                             }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core * c,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &        )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                             }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core * c,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &        )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                              }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core * c,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &        )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                             }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core * c,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &        )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                             }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core * c,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &   )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                              }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core * c,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &   )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                             }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core * c,
  const tag::ThisMeshIsPositive &                                   )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                            }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core * c,
  const tag::ThisMeshIsPositive &                                   )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                             }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core * c,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &        )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                             }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core * c,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &        )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                            }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core * c,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &        )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                             }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core * c,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &        )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                             }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::BuildCellsIfNec &, const tag::Around &, Cell::Core * c,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &   )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                              }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBuildCells &, const tag::Around &, Cell::Core * c,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &   )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "not yet implemented" << std::endl;
	exit ( 1 );                                             }

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


void CellIterator::Over::TwoVerticesOfSeg::reset ( ) // virtual from CellIterator::Core
{	assert ( this->seg_p );
	this->passage = 0;      }

void CellIterator::Over::TwoVerticesOfSeg::reset ( const tag::StartAt &, Cell::Core * cll )
// virtual from CellIterator::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot reset an iterator this way." << std::endl;
	exit ( 1 );                                                      }

void CellIterator::Over::TwoVerticesOfSeg::advance ( )
// virtual from CellIterator::Core
{	assert ( this->passage < 2 );
	this->passage++;                                                 }

bool CellIterator::Over::TwoVerticesOfSeg::in_range ( )
// virtual from CellIterator::Core
{	return this->passage < 2;   }

Cell CellIterator::Over::TwoVerticesOfSeg::NormalOrder::AsTheyAre::deref ( )
// virtual from CellIterator::Core
// iterate over the two vertices, first base (negative) then tip (positive)
{	assert ( this->seg_p );
	assert ( this->passage < 2 );
	if ( this->passage == 0 )
	{	assert ( not this->seg_p->base_attr.is_positive() );
		return this->seg_p->base_attr;                        }  // negative vertex
	else
	{	assert ( this->seg_p->tip_attr.is_positive() );
		return this->seg_p->tip_attr;                    } }  // positive vertex


Cell CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::AsTheyAre::deref ( )
// virtual from CellIterator::Core
// iterate over the two vertices, first tip (positive) then base (negative)
{	assert ( this->seg_p );
	assert ( this->passage < 2 );
	if ( this->passage == 0 )
  {	assert ( this->seg_p->tip_attr.is_positive() );
		return this->seg_p->tip_attr;                   }    // positive vertex
	else
	{	assert ( not this->seg_p->base_attr.is_positive() );
		return this->seg_p->base_attr;                        }  } // negative vertex


Cell CellIterator::Over::TwoVerticesOfSeg::NormalOrder::ForcePositive::deref ( )
// virtual from CellIterator::Core
// iterate over the two vertices, first base then tip (both positive)
{	assert ( this->seg_p );
	assert ( this->passage < 2 );
	if ( this->passage == 0 )
		return this->seg_p->base_attr.reverse ( tag::surely_exists );      // positive vertex
	else
		return this->seg_p->tip_attr;                                   }  // positive vertex


Cell CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::ForcePositive::deref ( )
// virtual from CellIterator::Core
// iterate over the two vertices, first tip then base (both positive)
{	assert ( this->seg_p );
	assert ( this->passage < 2 );
	if ( this->passage == 0 )
		return this->seg_p->tip_attr;                      // positive vertex
	else
		return this->seg_p->base_attr.reverse ( tag::surely_exists );  }   // positive vertex


Cell CellIterator::Over::TwoVerticesOfSeg
               ::NormalOrder::ReverseEachCell::AssumeCellsExist::deref ( )
// virtual from CellIterator::Core
// iterate over the two vertices, first base (negative) then tip (positive)
{	assert ( this->seg_p );
	assert ( this->passage < 2 );
	if ( this->passage == 0 )
		return this->seg_p->base_attr.reverse ( tag::surely_exists );     // positive vertex
	else
		return this->seg_p->tip_attr.reverse ( tag::surely_exists );   }  // negative vertex
	

Cell CellIterator::Over::TwoVerticesOfSeg
               ::NormalOrder::ReverseEachCell::BuildCells::deref ( )
// virtual from CellIterator::Core
// iterate over the two vertices, first base (negative) then tip (positive)
{	assert ( this->seg_p );
	assert ( this->passage < 2 );
	if ( this->passage == 0 )
		return this->seg_p->base_attr.reverse ( tag::surely_exists );     // positive vertex
	else
		return this->seg_p->tip_attr.reverse ( tag::build_if_not_exists );  }  // negative vertex
	

Cell CellIterator::Over::TwoVerticesOfSeg
               ::ReverseOrder::ReverseEachCell::AssumeCellsExist::deref ( )
// virtual from CellIterator::Core
// iterate over the two vertices, first tip (positive) then base (negative)
{	assert ( this->seg_p );
	assert ( this->passage < 2 );
	if ( this->passage == 0 )
		return this->seg_p->tip_attr.reverse ( tag::surely_exists );      // negative vertex
	else
		return this->seg_p->base_attr.reverse ( tag::surely_exists );  }  // positive vertex


Cell CellIterator::Over::TwoVerticesOfSeg
               ::ReverseOrder::ReverseEachCell::BuildCells::deref ( )
// virtual from CellIterator::Core
// iterate over the two vertices, first base (negative) then tip (positive)
{	assert ( this->seg_p );
	assert ( this->passage < 2 );
	if ( this->passage == 0 )
		return this->seg_p->base_attr.reverse ( tag::surely_exists );     // positive vertex
	else
		return this->seg_p->tip_attr.reverse ( tag::build_if_not_exists );  }  // negative vertex
	

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


void CellIterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder::reset ( )
// virtual from CellIterator::Core
{	assert ( this->msh );
	if ( this->msh->nb_of_segs == 0 )
	{	this->current_vertex = nullptr;  return;  }
	this->current_vertex = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* >
		( this->msh->first_ver.reverse(tag::surely_exists).core );
	this->last_vertex = this->msh->last_ver.core;
	if ( this->msh->first_ver.reverse(tag::surely_exists) == this->msh->last_ver ) // loop
	{	Cell rev_ver = this->current_vertex->reverse_attr;
		assert ( rev_ver.exists() );
		Cell seg = rev_ver.core->cell_behind_within[msh];
		this->current_vertex = tag::Util::assert_cast
			< Cell::Core*, Cell::Positive::Vertex* > ( seg.tip().core );  }
	assert ( this->current_vertex->is_positive() );
	assert ( this->last_vertex->is_positive() );                         }


void CellIterator::Over::VerticesOfConnectedOneDimMesh::ReverseOrder::reset ( )
// virtual from CellIterator::Core
{	assert ( this->msh );
	if ( this->msh->nb_of_segs == 0 )
	{	this->current_vertex = nullptr;  return;  }
	this->last_vertex = this->msh->first_ver.reverse(tag::surely_exists).core;
	this->current_vertex = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( this->msh->last_ver.core );
	if ( this->msh->first_ver.reverse(tag::surely_exists) == this->msh->last_ver ) // loop
	{	Cell seg = this->current_vertex->cell_behind_within[msh];
		this->current_vertex = tag::Util::assert_cast
			< Cell::Core*, Cell::Positive::Vertex* >
			( seg.base().reverse(tag::surely_exists).core );         }
	assert ( this->current_vertex->is_positive() );
	assert ( this->last_vertex->is_positive() );                                   }


void CellIterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder::reset
( const tag::StartAt &, Cell::Core * ver )  // virtual from CellIterator::Core
{	assert ( this->msh );
	assert ( this->msh->nb_of_segs > 0 );
	assert ( ver );
	assert ( ver->get_dim() == 0 );
	this->current_vertex = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( ver );
	if ( this->msh->first_ver.reverse(tag::surely_exists) == this->msh->last_ver ) // loop
	{	Cell seg = ver->cell_behind_within[msh];
		this->last_vertex = tag::Util::assert_cast
			< Cell::Core*, Cell::Positive::Vertex* >
			( seg.base().reverse(tag::surely_exists).core );  }
	else // open chain
		this->last_vertex = this->msh->last_ver.core;
	assert ( this->current_vertex->is_positive() );
	assert ( this->last_vertex->is_positive() );                                  }


void CellIterator::Over::VerticesOfConnectedOneDimMesh::ReverseOrder::reset
( const tag::StartAt &, Cell::Core * ver )  // virtual from CellIterator::Core
{	assert ( this->msh );
	assert ( this->msh->nb_of_segs > 0 );
	assert ( ver );
	assert ( ver->get_dim() == 0 );
	this->current_vertex = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( ver );
	if ( this->msh->first_ver.reverse(tag::surely_exists) == this->msh->last_ver ) // loop
	{	Cell rev_ver = ver->reverse_attr;
		assert ( rev_ver.exists() );
		Cell seg = rev_ver.core->cell_behind_within[msh];
		this->last_vertex = tag::Util::assert_cast
			< Cell::Core*, Cell::Positive::Vertex* > ( seg.tip().core );  }
	else // open chain
		this->last_vertex = this->msh->first_ver.reverse(tag::surely_exists).core;
	assert ( this->current_vertex->is_positive() );
	assert ( this->last_vertex->is_positive() );                                 }


Cell CellIterator::Over::VerticesOfConnectedOneDimMesh::deref ( )
// virtual from CellIterator::Core
{	return Cell ( tag::whose_core_is, this->current_vertex,
                tag::previously_existing, tag::surely_not_null );  }


void CellIterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder::advance ( )
// virtual from CellIterator::Core
{	assert ( this->msh );
	assert ( this->current_vertex );
	if ( this->current_vertex == this->last_vertex )              
	{	this->current_vertex = nullptr;  return;  }  // becomes out of range
	Cell rev_ver = this->current_vertex->reverse_attr;
	assert ( rev_ver.exists() );
	Cell rev_seg = rev_ver.core->cell_behind_within[msh];
	this->current_vertex = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( rev_seg.tip().core );
	assert ( this->current_vertex->is_positive() );                    }


void CellIterator::Over::VerticesOfConnectedOneDimMesh::ReverseOrder::advance ( )
// virtual from CellIterator::Core
{	assert ( this->msh );
	assert ( this->current_vertex );
	if ( this->current_vertex == this->last_vertex )              
	{	this->current_vertex = nullptr;  return;  }  // becomes out of range
	Cell seg = this->current_vertex->cell_behind_within[msh];
	this->current_vertex = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* >
		( seg.base().reverse(tag::surely_exists).core );
	assert ( this->current_vertex->is_positive() );            }


bool CellIterator::Over::VerticesOfConnectedOneDimMesh::in_range ( )
// virtual from CellIterator::Core
{	return this->current_vertex;  }

//-----------------------------------------------------------------------------------------


void CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder::reset ( )
// virtual from CellIterator::Core
{	assert ( this->msh );
	if ( this->msh->nb_of_segs == 0 )
	{ this->current_segment = nullptr;  return;  }
	Cell neg_ver = this->msh->first_ver;
	assert ( neg_ver.exists() );  assert ( not neg_ver.is_positive() );
	this->current_segment = neg_ver.core->cell_behind_within[this->msh].core;
	assert ( this->current_segment );
	this->last_vertex = this->msh->last_ver.core;
	assert ( this->last_vertex->is_positive() );                               }


void CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder::reset ( )
// virtual from CellIterator::Core
{	assert ( this->msh );
	if ( this->msh->nb_of_segs == 0 )
	{ this->current_segment = nullptr;  return;  }
	Cell pos_ver = this->msh->last_ver;
	assert ( pos_ver.exists() );  assert ( pos_ver.is_positive() );
	this->current_segment = pos_ver.core->cell_behind_within[this->msh].core;
	assert ( this->current_segment );
	this->last_vertex = this->msh->first_ver.core;
	assert ( not this->last_vertex->is_positive() );                           }


void CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder::reset
( const tag::StartAt &, Cell::Core * seg )  // virtual from CellIterator::Core
{	assert ( this->msh );  assert ( this->msh->nb_of_segs > 0 );
	assert ( seg );  assert ( seg->get_dim() == 1 );
	this->current_segment = seg;
	if ( this->msh->first_ver.reverse(tag::surely_exists) == this->msh->last_ver ) // loop
		this->last_vertex = seg->base().reverse(tag::surely_exists).core;
	else  // open chain
		this->last_vertex = this->msh->last_ver.core;
	assert ( this->last_vertex->is_positive() );                         }


void CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder::reset
( const tag::StartAt &, Cell::Core * seg )  // virtual from CellIterator::Core
{	assert ( this->msh );  assert ( this->msh->nb_of_segs > 0 );
	assert ( seg );  assert ( seg->get_dim() == 1 );
	this->current_segment = seg;
	if ( this->msh->first_ver.reverse(tag::surely_exists) == this->msh->last_ver ) // loop
		this->last_vertex = seg->tip().reverse(tag::surely_exists).core;
	else  // open chain
		this->last_vertex = this->msh->first_ver.core;
	assert ( not this->last_vertex->is_positive() );                    }


Cell CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder::AsTheyAre::deref ( )
// virtual from CellIterator::Core
{	return Cell ( tag::whose_core_is, this->current_segment,
                tag::previously_existing, tag::surely_not_null );  }


Cell CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder::AsTheyAre::deref ( )
// virtual from CellIterator::Core
{	return Cell ( tag::whose_core_is, this->current_segment,
                tag::previously_existing, tag::surely_not_null );  }


Cell CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder::ForcePositive::deref ( )
// virtual from CellIterator::Core
{	assert ( this->current_segment );
	return this->current_segment->get_positive();  }


Cell CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder::ForcePositive::deref ( )
// virtual from CellIterator::Core
{	assert ( this->current_segment );
	return this->current_segment->get_positive();  }


Cell CellIterator::Over::SegmentsOfConnectedOneDimMesh::
	NormalOrder::ReverseEachCell::AssumeCellsExist::deref ( )
// virtual from CellIterator::Core
{	assert ( this->current_segment );
	assert ( this->current_segment->reverse_attr.exists() );
	return this->current_segment->reverse_attr;              }


Cell CellIterator::Over::SegmentsOfConnectedOneDimMesh::
	NormalOrder::ReverseEachCell::BuildCells::deref ( )
// virtual from CellIterator::Core
{	assert ( this->current_segment );
	if ( not this->current_segment->reverse_attr.exists() )
		//	this->current_segment->reverse_attr =
		//		Cell ( tag::whose_core_is, this->current_segment->build_reverse
		//	         ( tag::one_dummy_wrapper ), tag::freshly_created         );
	{	Cell::Positive * tcs = tag::Util::assert_cast
			< Cell::Core*, Cell::Positive* > ( this->current_segment );
		tcs->reverse_attr.core = tcs->build_reverse ( tag::one_dummy_wrapper );  }
	assert ( this->current_segment->reverse_attr.exists() );
	return this->current_segment->reverse_attr;                                  }


Cell CellIterator::Over::SegmentsOfConnectedOneDimMesh::
	ReverseOrder::ReverseEachCell::AssumeCellsExist::deref ( )
// virtual from CellIterator::Core
{	assert ( this->current_segment );
	assert ( this->current_segment->reverse_attr.exists() );
	return this->current_segment->reverse_attr;              }


Cell CellIterator::Over::SegmentsOfConnectedOneDimMesh::
	ReverseOrder::ReverseEachCell::BuildCells::deref ( )
// virtual from CellIterator::Core
{	assert ( this->current_segment );
	if ( not this->current_segment->reverse_attr.exists() )
		//	this->current_segment->reverse_attr =
		//		Cell ( tag::whose_core_is, this->current_segment->build_reverse
		//	         ( tag::one_dummy_wrapper ), tag::freshly_created         );
	{	Cell::Positive * tcs = tag::Util::assert_cast
			< Cell::Core*, Cell::Positive* > ( this->current_segment );
		tcs->reverse_attr.core = tcs->build_reverse ( tag::one_dummy_wrapper );  }
	assert ( this->current_segment->reverse_attr.exists() );
	return this->current_segment->reverse_attr;                                  }


void CellIterator::Over::SegmentsOfConnectedOneDimMesh::NormalOrder::advance ( )
// virtual from CellIterator::Core
{	assert ( this->current_segment );
	Cell::Core * ver = this->current_segment->tip().core;
	assert ( ver );  assert ( ver->is_positive() );
	if ( ver == this->last_vertex )  // becomes out of range
	{	this->current_segment = nullptr;  return;  }
	Cell::Core * neg_ver = ver->reverse_attr.core;
	assert ( neg_ver );
	this->current_segment = neg_ver->cell_behind_within[msh].core;
	assert ( this->current_segment );                               }


void CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder::advance ( )
// virtual from CellIterator::Core
{	assert ( this->current_segment );
	Cell::Core * ver = this->current_segment->base().core;
	assert ( ver );  assert ( not ver->is_positive() );
	if ( ver == this->last_vertex )  // becomes out of range
	{	this->current_segment = nullptr;  return;  }
	Cell::Core * pos_ver = ver->reverse_attr.core;
	assert ( pos_ver );  assert ( pos_ver->is_positive() );
	this->current_segment = pos_ver->cell_behind_within[msh].core;
	assert ( this->current_segment );                               }


bool CellIterator::Over::SegmentsOfConnectedOneDimMesh::in_range ( )
// virtual from CellIterator::Core
{	return this->current_segment;  }

//-----------------------------------------------------------------------------------------


void CellIterator::Over::CellsOfFuzzyMesh::reset ( ) // virtual from CellIterator::Core
{	this->iter = this->list.begin();  }


void CellIterator::Over::CellsOfFuzzyMesh::reset ( const tag::StartAt &, Cell::Core * cll )
// virtual from CellIterator::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot reset an iterator this way." << std::endl;
	exit ( 1 );                                                      }


void CellIterator::Over::CellsOfFuzzyMesh::advance ( )
// virtual from CellIterator::Core
{	this->iter++;  }


bool CellIterator::Over::CellsOfFuzzyMesh::in_range ( )
// virtual from CellIterator::Core
{	return this->iter != this->list.end();  }


Cell CellIterator::Over::CellsOfFuzzyMesh::AsTheyAre::deref ( )
// virtual from CellIterator::Core
{	return * this->iter;  }


Cell CellIterator::Over::CellsOfFuzzyMesh::ReverseEachCell::AssumeCellsExist::deref ( )
// virtual from CellIterator::Core
{	return ( * this->iter ).reverse ( tag::surely_exists );  }


Cell CellIterator::Over::CellsOfFuzzyMesh::ForcePositive::deref ( )
// virtual from CellIterator::Core
{	return ( * this->iter ).get_positive();  }

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


void CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder::reset ( )
// virtual from CellIterator::Core

{	this->current_vertex = this->first_vertex;
	Mesh m ( tag::whose_core_is, this->msh_p, tag::previously_existing, tag::is_positive );
	Cell new_ver ( tag::whose_core_is, this->current_vertex,
                 tag::previously_existing, tag::surely_not_null );
	assert ( new_ver.tip().core == this->center );
	this->current_segment = m.cell_in_front_of ( new_ver, tag::may_not_exist ).core;        }

// the above could get simpler and faster by avoiding the use of wrappers


void CellIterator::Around::OneCell::OfCodimTwo::OverVertices::ReverseOrder::reset ( )
// virtual from CellIterator::Core

{	this->current_vertex = this->first_vertex;
	Mesh m ( tag::whose_core_is, this->msh_p, tag::previously_existing, tag::is_positive );
	Cell new_ver ( tag::whose_core_is, this->current_vertex,
                 tag::previously_existing, tag::surely_not_null );
	assert ( new_ver.tip().core == this->center );
	this->current_segment = m.cell_behind ( new_ver, tag::may_not_exist ).core;             }

// the above could get simpler and faster by avoiding the use of wrappers


void CellIterator::Around::OneCell::OfCodimTwo
	::OverVertices::NormalOrder::AssumeCellsExist::AsTheyAre::reset
( const tag::StartAt &, Cell::Core * cll )
// virtual from CellIterator::Core

{ // closed loop or open chain ?
	// in other words : central cell is interior or on the boundary ?
	Mesh m ( tag::whose_core_is, this->msh_p, tag::previously_existing, tag::is_positive );
	Cell fv ( tag::whose_core_is, this->first_vertex,
            tag::previously_existing, tag::surely_not_null );
	this->current_vertex = cll;
	Cell new_ver ( tag::whose_core_is, cll, tag::previously_existing, tag::surely_not_null );
	// center may have dimension different from zero, so we should not use 'tip' !!
	if ( this->center->is_positive() ) assert ( new_ver.tip().core == this->center );
	else assert ( new_ver.base().core == this->center );
	this->current_segment = m.cell_in_front_of ( new_ver, tag::may_not_exist ).core;
	if ( this->last_segment )
	if ( this->last_segment == m.cell_behind ( fv, tag::may_not_exist ).core )
	// closed loop, i.e., center is interior to the mesh
	{	this->first_vertex = cll;
		this->last_segment = m.cell_behind ( new_ver, tag::may_not_exist ).core;  }
}


void CellIterator::Around::OneCell::OfCodimTwo
	::OverVertices::NormalOrder::AssumeCellsExist::ForcePositive::reset
( const tag::StartAt &, Cell::Core * cll )
// virtual from CellIterator::Core

{ // closed loop or open chain ?
	// in other words : central cell is interior or on the boundary ?
	Mesh m ( tag::whose_core_is, this->msh_p, tag::previously_existing, tag::is_positive );
	Cell fv ( tag::whose_core_is, this->first_vertex,
            tag::previously_existing, tag::surely_not_null );
	this->current_vertex = cll;
	Cell new_ver ( tag::whose_core_is, cll, tag::previously_existing, tag::surely_not_null );
	// center may have dimension different from zero, so we should not use 'tip' !!
	if ( this->center->is_positive() ) assert ( new_ver.tip().core == this->center );
	else assert ( new_ver.base().core == this->center );
	this->current_segment = m.cell_in_front_of ( new_ver, tag::may_not_exist ).core;
	if ( this->last_segment )
	if ( this->last_segment == m.cell_behind ( fv, tag::may_not_exist ).core )
	// closed loop, i.e., center is interior to the mesh
	{	this->first_vertex = cll;
		this->last_segment = m.cell_behind ( new_ver, tag::may_not_exist ).core;  }
}


void CellIterator::Around::OneCell::OfCodimTwo
	::OverVertices::NormalOrder::AssumeCellsExist::ReverseEachCell::reset
( const tag::StartAt &, Cell::Core * cll )
// virtual from CellIterator::Core

{ // closed loop or open chain ?
	// in other words : central cell is interior or on the boundary ?
	Mesh m ( tag::whose_core_is, this->msh_p, tag::previously_existing, tag::is_positive );
	Cell fv ( tag::whose_core_is, this->first_vertex,
            tag::previously_existing, tag::surely_not_null );
	this->current_vertex = cll;
	Cell new_ver ( tag::whose_core_is, cll, tag::previously_existing, tag::surely_not_null );
	// center may have dimension different from zero, so we should not use 'tip' !!
	if ( this->center->is_positive() ) assert ( new_ver.tip().core == this->center );
	else assert ( new_ver.base().core == this->center );
	this->current_segment = m.cell_in_front_of ( new_ver, tag::may_not_exist ).core;
	if ( this->last_segment )
	if ( this->last_segment == m.cell_behind ( fv, tag::may_not_exist ).core )
	// closed loop, i.e., center is interior to the mesh
	{	this->first_vertex = cll;
		this->last_segment = m.cell_behind ( new_ver, tag::may_not_exist ).core;  }
}


void CellIterator::Around::OneCell::OfCodimTwo
	::OverVertices::NormalOrder::BuildReverseCells::AsTheyAre::reset
( const tag::StartAt &, Cell::Core * cll )
// virtual from CellIterator::Core

{ // closed loop or open chain ?
	// in other words : central cell is interior or on the boundary ?
	Mesh m ( tag::whose_core_is, this->msh_p, tag::previously_existing, tag::is_positive );
	Cell fv ( tag::whose_core_is, this->first_vertex,
            tag::previously_existing, tag::surely_not_null );
	this->current_vertex = cll;
	Cell new_ver ( tag::whose_core_is, cll, tag::previously_existing, tag::surely_not_null );
	// center may have dimension different from zero, so we should not use 'tip' !!
	if ( this->center->is_positive() ) assert ( new_ver.tip().core == this->center );
	else assert ( new_ver.base().core == this->center );
	this->current_segment = m.cell_in_front_of ( new_ver, tag::may_not_exist ).core;
	if ( this->last_segment )
	if ( this->last_segment == m.cell_behind ( fv, tag::may_not_exist ).core )
	// closed loop, i.e., center is interior to the mesh
	{	this->first_vertex = cll;
		this->last_segment = m.cell_behind ( new_ver, tag::may_not_exist ).core;  }
}


void CellIterator::Around::OneCell::OfCodimTwo
	::OverVertices::NormalOrder::BuildReverseCells::ForcePositive::reset
( const tag::StartAt &, Cell::Core * cll )
// virtual from CellIterator::Core

{ // closed loop or open chain ?
	// in other words : central cell is interior or on the boundary ?
	Mesh m ( tag::whose_core_is, this->msh_p, tag::previously_existing, tag::is_positive );
	Cell fv ( tag::whose_core_is, this->first_vertex,
            tag::previously_existing, tag::surely_not_null );
	this->current_vertex = cll;
	Cell new_ver ( tag::whose_core_is, cll, tag::previously_existing, tag::surely_not_null );
	// center may have dimension different from zero, so we should not use 'tip' !!
	if ( this->center->is_positive() ) assert ( new_ver.tip().core == this->center );
	else assert ( new_ver.base().core == this->center );
	this->current_segment = m.cell_in_front_of ( new_ver, tag::may_not_exist ).core;
	if ( this->last_segment )
	if ( this->last_segment == m.cell_behind ( fv, tag::may_not_exist ).core )
	// closed loop, i.e., center is interior to the mesh
	{	this->first_vertex = cll;
		this->last_segment = m.cell_behind ( new_ver, tag::may_not_exist ).core;  }
}


void CellIterator::Around::OneCell::OfCodimTwo
	::OverVertices::NormalOrder::BuildReverseCells::ReverseEachCell::reset
( const tag::StartAt &, Cell::Core * cll )
// virtual from CellIterator::Core

{ // closed loop or open chain ?
	// in other words : central cell is interior or on the boundary ?
	Mesh m ( tag::whose_core_is, this->msh_p, tag::previously_existing, tag::is_positive );
	Cell fv ( tag::whose_core_is, this->first_vertex,
            tag::previously_existing, tag::surely_not_null );
	this->current_vertex = cll;
	Cell new_ver ( tag::whose_core_is, cll, tag::previously_existing, tag::surely_not_null );
	// center may have dimension different from zero, so we should not use 'tip' !!
	if ( this->center->is_positive() ) assert ( new_ver.tip().core == this->center );
	else assert ( new_ver.base().core == this->center );
	this->current_segment = m.cell_in_front_of ( new_ver, tag::may_not_exist ).core;
	if ( this->last_segment )
	if ( this->last_segment == m.cell_behind ( fv, tag::may_not_exist ).core )
	// closed loop, i.e., center is interior to the mesh
	{	this->first_vertex = cll;
		this->last_segment = m.cell_behind ( new_ver, tag::may_not_exist ).core;  }
}


void CellIterator::Around::OneCell::OfCodimTwo
	::OverVertices::ReverseOrder::AssumeCellsExist::AsTheyAre::reset
( const tag::StartAt &, Cell::Core * cll )
// virtual from CellIterator::Core

{ // closed loop or open chain ?
	// in other words : central cell is interior or on the boundary ?
	Mesh m ( tag::whose_core_is, this->msh_p, tag::previously_existing, tag::is_positive );
	Cell fv ( tag::whose_core_is, this->first_vertex,
            tag::previously_existing, tag::surely_not_null );
	this->current_vertex = cll;
	Cell new_ver ( tag::whose_core_is, cll, tag::previously_existing, tag::surely_not_null );
	assert ( new_ver.tip().core == this->center );
	this->current_segment = m.cell_behind ( new_ver, tag::may_not_exist ).core;
	std::cout << "reset ::OverVertices::ReverseOrder";
	if ( this->last_segment )
	if ( this->last_segment == m.cell_in_front_of ( fv, tag::may_not_exist ).core )
	// closed loop, i.e., center is interior to the mesh
	{	this->first_vertex = cll;
		std::cout << ", center is interior to the mesh " << this->last_segment;
		this->last_segment = m.cell_in_front_of ( new_ver, tag::may_not_exist ).core;  }
	std::cout << std::endl;
}
	

void CellIterator::Around::OneCell::OfCodimTwo
	::OverVertices::ReverseOrder::AssumeCellsExist::ForcePositive::reset
( const tag::StartAt &, Cell::Core * cll )
// virtual from CellIterator::Core

{ // closed loop or open chain ?
	// in other words : central cell is interior or on the boundary ?
	Mesh m ( tag::whose_core_is, this->msh_p, tag::previously_existing, tag::is_positive );
	Cell fv ( tag::whose_core_is, this->first_vertex,
            tag::previously_existing, tag::surely_not_null );
	this->current_vertex = cll;
	Cell new_ver ( tag::whose_core_is, cll, tag::previously_existing, tag::surely_not_null );
	assert ( new_ver.tip().core == this->center );
	this->current_segment = m.cell_behind ( new_ver, tag::may_not_exist ).core;
	std::cout << "reset ::OverVertices::ReverseOrder";
	if ( this->last_segment )
	if ( this->last_segment == m.cell_in_front_of ( fv, tag::may_not_exist ).core )
	// closed loop, i.e., center is interior to the mesh
	{	this->first_vertex = cll;
		std::cout << ", center is interior to the mesh " << this->last_segment;
		this->last_segment = m.cell_in_front_of ( new_ver, tag::may_not_exist ).core;  }
	std::cout << std::endl;
}
	

void CellIterator::Around::OneCell::OfCodimTwo
	::OverVertices::ReverseOrder::AssumeCellsExist::ReverseEachCell::reset
( const tag::StartAt &, Cell::Core * cll )
// virtual from CellIterator::Core

{ // closed loop or open chain ?
	// in other words : central cell is interior or on the boundary ?
	Mesh m ( tag::whose_core_is, this->msh_p, tag::previously_existing, tag::is_positive );
	Cell fv ( tag::whose_core_is, this->first_vertex,
            tag::previously_existing, tag::surely_not_null );
	this->current_vertex = cll;
	Cell new_ver ( tag::whose_core_is, cll, tag::previously_existing, tag::surely_not_null );
	assert ( new_ver.tip().core == this->center );
	this->current_segment = m.cell_behind ( new_ver, tag::may_not_exist ).core;
	std::cout << "reset ::OverVertices::ReverseOrder";
	if ( this->last_segment )
	if ( this->last_segment == m.cell_in_front_of ( fv, tag::may_not_exist ).core )
	// closed loop, i.e., center is interior to the mesh
	{	this->first_vertex = cll;
		std::cout << ", center is interior to the mesh " << this->last_segment;
		this->last_segment = m.cell_in_front_of ( new_ver, tag::may_not_exist ).core;  }
	std::cout << std::endl;
}
	

void CellIterator::Around::OneCell::OfCodimTwo
	::OverVertices::ReverseOrder::BuildReverseCells::AsTheyAre::reset
( const tag::StartAt &, Cell::Core * cll )
// virtual from CellIterator::Core

{ // closed loop or open chain ?
	// in other words : central cell is interior or on the boundary ?
	Mesh m ( tag::whose_core_is, this->msh_p, tag::previously_existing, tag::is_positive );
	Cell fv ( tag::whose_core_is, this->first_vertex,
            tag::previously_existing, tag::surely_not_null );
	this->current_vertex = cll;
	Cell new_ver ( tag::whose_core_is, cll, tag::previously_existing, tag::surely_not_null );
	assert ( new_ver.tip().core == this->center );
	this->current_segment = m.cell_behind ( new_ver, tag::may_not_exist ).core;
	std::cout << "reset ::OverVertices::ReverseOrder";
	if ( this->last_segment )
	if ( this->last_segment == m.cell_in_front_of ( fv, tag::may_not_exist ).core )
	// closed loop, i.e., center is interior to the mesh
	{	this->first_vertex = cll;
		std::cout << ", center is interior to the mesh " << this->last_segment;
		this->last_segment = m.cell_in_front_of ( new_ver, tag::may_not_exist ).core;  }
	std::cout << std::endl;
}
	

void CellIterator::Around::OneCell::OfCodimTwo
	::OverVertices::ReverseOrder::BuildReverseCells::ForcePositive::reset
( const tag::StartAt &, Cell::Core * cll )
// virtual from CellIterator::Core

{ // closed loop or open chain ?
	// in other words : central cell is interior or on the boundary ?
	Mesh m ( tag::whose_core_is, this->msh_p, tag::previously_existing, tag::is_positive );
	Cell fv ( tag::whose_core_is, this->first_vertex,
            tag::previously_existing, tag::surely_not_null );
	this->current_vertex = cll;
	Cell new_ver ( tag::whose_core_is, cll, tag::previously_existing, tag::surely_not_null );
	assert ( new_ver.tip().core == this->center );
	this->current_segment = m.cell_behind ( new_ver, tag::may_not_exist ).core;
	std::cout << "reset ::OverVertices::ReverseOrder";
	if ( this->last_segment )
	if ( this->last_segment == m.cell_in_front_of ( fv, tag::may_not_exist ).core )
	// closed loop, i.e., center is interior to the mesh
	{	this->first_vertex = cll;
		std::cout << ", center is interior to the mesh " << this->last_segment;
		this->last_segment = m.cell_in_front_of ( new_ver, tag::may_not_exist ).core;  }
	std::cout << std::endl;
}
	

void CellIterator::Around::OneCell::OfCodimTwo
	::OverVertices::ReverseOrder::BuildReverseCells::ReverseEachCell::reset
( const tag::StartAt &, Cell::Core * cll )
// virtual from CellIterator::Core

{ // closed loop or open chain ?
	// in other words : central cell is interior or on the boundary ?
	Mesh m ( tag::whose_core_is, this->msh_p, tag::previously_existing, tag::is_positive );
	Cell fv ( tag::whose_core_is, this->first_vertex,
            tag::previously_existing, tag::surely_not_null );
	this->current_vertex = cll;
	Cell new_ver ( tag::whose_core_is, cll, tag::previously_existing, tag::surely_not_null );
	assert ( new_ver.tip().core == this->center );
	this->current_segment = m.cell_behind ( new_ver, tag::may_not_exist ).core;
	std::cout << "reset ::OverVertices::ReverseOrder";
	if ( this->last_segment )
	if ( this->last_segment == m.cell_in_front_of ( fv, tag::may_not_exist ).core )
	// closed loop, i.e., center is interior to the mesh
	{	this->first_vertex = cll;
		std::cout << ", center is interior to the mesh " << this->last_segment;
		this->last_segment = m.cell_in_front_of ( new_ver, tag::may_not_exist ).core;  }
	std::cout << std::endl;
}
	

Cell CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::NormalOrder::AssumeCellsExist::AsTheyAre::deref ( )
// virtual from CellIterator::Core
{	return Cell ( tag::whose_core_is, this->current_vertex,
                tag::previously_existing, tag::surely_not_null );  }


Cell CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::NormalOrder::BuildReverseCells::AsTheyAre::deref ( )
// virtual from CellIterator::Core
{	return Cell ( tag::whose_core_is, this->current_vertex,
                tag::previously_existing, tag::surely_not_null );  }


Cell CellIterator::Around::OneCell
         ::OfCodimTwo::OverVertices::NormalOrder::AssumeCellsExist::ForcePositive::deref ( )
// virtual from CellIterator::Core
{	assert ( this->current_vertex );
	return this->current_vertex->get_positive();  }


Cell CellIterator::Around::OneCell
         ::OfCodimTwo::OverVertices::NormalOrder::BuildReverseCells::ForcePositive::deref ( )
// virtual from CellIterator::Core
{	assert ( this->current_vertex );
	return this->current_vertex->get_positive();  }


Cell CellIterator::Around::OneCell::OfCodimTwo::OverVertices::
	NormalOrder::AssumeCellsExist::ReverseEachCell::deref ( )
// virtual from CellIterator::Core
{	assert ( this->current_vertex );
	assert ( this->current_vertex->reverse_attr.exists() );
	return this->current_vertex->reverse_attr;              }


Cell CellIterator::Around::OneCell::OfCodimTwo::OverVertices::
	NormalOrder::BuildReverseCells::ReverseEachCell::deref ( )
// virtual from CellIterator::Core

{	assert ( this->current_vertex );
	if ( not this->current_vertex->reverse_attr.exists() )
	//	this->current_vertex->reverse_attr =
	//		Cell ( tag::whose_core_is, this->current_vertex->build_reverse
	//	         ( tag::one_dummy_wrapper ), tag::freshly_created         );
	{	Cell::Positive * tcv = tag::Util::assert_cast
			< Cell::Core*, Cell::Positive* > ( this->current_vertex );
		tcv->reverse_attr.core = tcv->build_reverse ( tag::one_dummy_wrapper );  }
	assert ( this->current_vertex->reverse_attr.exists() );
	return this->current_vertex->reverse_attr;              }


Cell CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::ReverseOrder::AssumeCellsExist::AsTheyAre::deref ( )
// virtual from CellIterator::Core
{	return Cell ( tag::whose_core_is, this->current_vertex,
                tag::previously_existing, tag::surely_not_null );  }


Cell CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::ReverseOrder::BuildReverseCells::AsTheyAre::deref ( )
// virtual from CellIterator::Core
{	return Cell ( tag::whose_core_is, this->current_vertex,
                tag::previously_existing, tag::surely_not_null );  }


Cell CellIterator::Around::OneCell
         ::OfCodimTwo::OverVertices::ReverseOrder::AssumeCellsExist::ForcePositive::deref ( )
// virtual from CellIterator::Core
{	assert ( this->current_vertex );
	return this->current_vertex->get_positive();  }


Cell CellIterator::Around::OneCell
         ::OfCodimTwo::OverVertices::ReverseOrder::BuildReverseCells::ForcePositive::deref ( )
// virtual from CellIterator::Core
{	assert ( this->current_vertex );
	return this->current_vertex->get_positive();  }


Cell CellIterator::Around::OneCell::OfCodimTwo::OverVertices::
	ReverseOrder::AssumeCellsExist::ReverseEachCell::deref ( )
// virtual from CellIterator::Core
{	assert ( this->current_vertex );
	assert ( this->current_vertex->reverse_attr.exists() );
	return this->current_vertex->reverse_attr;              }


Cell CellIterator::Around::OneCell::OfCodimTwo::OverVertices::
	ReverseOrder::BuildReverseCells::ReverseEachCell::deref ( )
// virtual from CellIterator::Core

{	assert ( this->current_vertex );
	if ( not this->current_vertex->reverse_attr.exists() )
	//	this->current_vertex->reverse_attr =
	//		Cell ( tag::whose_core_is, this->current_vertex->build_reverse
	//	         ( tag::one_dummy_wrapper ), tag::freshly_created         );
	{	Cell::Positive * tcv = tag::Util::assert_cast
			< Cell::Core*, Cell::Positive* > ( this->current_vertex );
		tcv->reverse_attr.core = tcv->build_reverse ( tag::one_dummy_wrapper );  }
	assert ( this->current_vertex->reverse_attr.exists() );
	return this->current_vertex->reverse_attr;              }


Cell CellIterator::Around::OneCell
         ::OfCodimTwo::WorkAround2D::NormalOrder::BuildReverseCells::deref ( )
// virtual from CellIterator::Core,
// defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder::
// ::BuildReverseCells::AsTheyAre, here overriden to return the reversed base of the segment
{	Cell seg = CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::NormalOrder::BuildReverseCells::AsTheyAre::deref();
	return seg.base().reverse();                                          }


Cell CellIterator::Around::OneCell
         ::OfCodimTwo::WorkAround2D::ReverseOrder::BuildReverseCells::deref ( )
// virtual from CellIterator::Core,
// defined by CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder::
// ::BuildReverseCells::AsTheyAre, here overriden to return the reversed base of the segment
{	Cell seg = CellIterator::Around::OneCell::OfCodimTwo
		::OverVertices::ReverseOrder::BuildReverseCells::AsTheyAre::deref();
	return seg.base().reverse();                                          }


void CellIterator::Around::OneCell::OfCodimTwo::OverVertices::NormalOrder::advance ( )
// virtual from CellIterator::Core

{	Cell this_center ( tag::whose_core_is, this->center,
	                   tag::previously_existing, tag::surely_not_null );
	Mesh m ( tag::whose_core_is, this->msh_p, tag::previously_existing, tag::is_positive );
	if ( this->current_segment == nullptr )  // out of range
	{	this->current_vertex = nullptr; return;  }
	Cell new_ver = this->current_segment->boundary().cell_behind
		( this_center, tag::surely_exists );
	this->current_vertex = new_ver.core;
	this->current_segment = m.cell_in_front_of ( new_ver, tag::may_not_exist ).core;
	if ( this->current_segment == this->last_segment ) this->current_segment = nullptr;  }
	// if we are near the boundary, current_segment may be nullptr
	// the calling program must check this using the 'in_range' method

// the above could get simpler and faster by avoiding the use of wrappers


void CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::ReverseOrder::AssumeCellsExist::advance ( )
// virtual from CellIterator::Core

{	Cell this_center ( tag::whose_core_is, this->center,
	                   tag::previously_existing, tag::surely_not_null );
	Mesh m ( tag::whose_core_is, this->msh_p, tag::previously_existing, tag::is_positive );
	if ( this->current_segment == nullptr )  // out of range
	{	this->current_vertex = nullptr; return;  }
	Cell new_ver = this->current_segment->boundary().cell_in_front_of
		( this_center, tag::surely_exists ).reverse ( tag::surely_exists );
	this->current_vertex = new_ver.core;
	this->current_segment = m.cell_behind ( new_ver, tag::may_not_exist ).core;
	if ( this->current_segment == this->last_segment ) this->current_segment = nullptr;  }
	// if we are near the boundary, current_segment may be nullptr
	// the calling program must check this using the 'in_range' method

// the above could get simpler and faster by avoiding the use of wrappers


void CellIterator::Around::OneCell::OfCodimTwo
         ::OverVertices::ReverseOrder::BuildReverseCells::advance ( )
// virtual from CellIterator::Core

{	Cell this_center ( tag::whose_core_is, this->center,
	                   tag::previously_existing, tag::surely_not_null );
	Mesh m ( tag::whose_core_is, this->msh_p, tag::previously_existing, tag::is_positive );
	if ( this->current_segment == nullptr )  // out of range
	{	this->current_vertex = nullptr; return;  }
	Cell new_ver = this->current_segment->boundary().cell_in_front_of
		( this_center, tag::surely_exists ).reverse ( tag::build_if_not_exists );
	this->current_vertex = new_ver.core;
	this->current_segment = m.cell_behind ( new_ver, tag::may_not_exist ).core;
	if ( this->current_segment == this->last_segment ) this->current_segment = nullptr;  }
	// if we are near the boundary, current_segment may be nullptr
	// the calling program must check this using the 'in_range' method

// the above could get simpler and faster by avoiding the use of wrappers


bool CellIterator::Around::OneCell::OfCodimTwo::OverVertices::in_range ( )
// virtual from CellIterator::Core
{	return this->current_vertex;   }

//-----------------------------------------------------------------------------------------



void CellIterator::Around::OneCell::OfCodimTwo::OverSegments::reset ( )
// virtual from CellIterator::Core

{	this->current_segment = this->first_segment;
	assert ( this->center->is_positive() );
	assert ( this->current_segment );             }	

// the above could get simpler and faster by using a method like
// inline Cell::Core * Mesh::Core::cell_behind_ptr
// ( const Cell::Core * face, const tag::MayNotExist & ) const
// which is not difficult to implement


void CellIterator::Around::OneCell::OfCodimTwo::OverSegments::NormalOrder::reset
( const tag::StartAt &, Cell::Core * cll_p )  // virtual from CellIterator::Core

{	assert ( this->center->is_positive() );
	Cell cll ( tag::whose_core_is, cll_p,
             tag::previously_existing, tag::surely_not_null );
	Mesh m ( tag::whose_core_is, this->msh_p, tag::previously_existing, tag::is_positive );
	assert ( cll.belongs_to ( m, tag::oriented ) );
	this->current_segment = cll_p;
	if ( this->last_vertex )  // closed loop
	{	this->first_segment = cll_p;
		Cell this_center ( tag::whose_core_is, this->center,
		                   tag::previously_existing, tag::surely_not_null );
		this->last_vertex = cll.boundary().cell_in_front_of
			( this_center, tag::surely_exists ). reverse ( tag::surely_exists ). core;  }
	assert ( this->current_segment );                                                       }

// the above could get simpler and faster by using a method like
// inline Cell::Core * Mesh::Core::cell_behind_ptr
// ( const Cell::Core * face, const tag::MayNotExist & ) const
// which is not difficult to implement


void CellIterator::Around::OneCell::OfCodimTwo::OverSegments::ReverseOrder::reset
( const tag::StartAt &, Cell::Core * cll_p )  // virtual from CellIterator::Core

{	assert ( this->center->is_positive() );
	Cell cll ( tag::whose_core_is, cll_p,
             tag::previously_existing, tag::surely_not_null );
	Mesh m ( tag::whose_core_is, this->msh_p, tag::previously_existing, tag::is_positive );
	assert ( cll.belongs_to ( m, tag::oriented ) );
	this->current_segment = cll_p;
	if ( this->last_vertex )  // closed loop
	{	this->first_segment = cll_p;
		Cell this_center ( tag::whose_core_is, this->center,
		                   tag::previously_existing, tag::surely_not_null );
		this->last_vertex = cll.boundary().cell_behind
			( this_center, tag::surely_exists ). reverse ( tag::surely_exists ). core;  }
	assert ( this->current_segment );                                                       }

// the above could get simpler and faster by using a method like
// inline Cell::Core * Mesh::Core::cell_behind_ptr
// ( const Cell::Core * face, const tag::MayNotExist & ) const
// which is not difficult to implement


void CellIterator::Around::OneCell::OfCodimTwo::OverSegments::NormalOrder::advance ( )
// virtual from CellIterator::Core

{	Cell this_center ( tag::whose_core_is, this->center,
	                   tag::previously_existing, tag::surely_not_null );
	Mesh m ( tag::whose_core_is, this->msh_p, tag::previously_existing, tag::is_positive );
	assert ( this->current_segment );
	Cell new_ver = this->current_segment->boundary().cell_behind
		( this_center, tag::surely_exists );
	if ( new_ver.core == this->last_vertex )
		this->current_segment = nullptr;
	else
		this->current_segment = m.cell_in_front_of ( new_ver, tag::may_not_exist ).core;     }

// the above could get simpler and faster by avoiding the use of wrappers


void CellIterator::Around::OneCell::OfCodimTwo::OverSegments::ReverseOrder::advance ( )
// virtual from CellIterator::Core

{	Cell this_center ( tag::whose_core_is, this->center,
	                   tag::previously_existing, tag::surely_not_null );
	Mesh m ( tag::whose_core_is, this->msh_p, tag::previously_existing, tag::is_positive );
	assert ( this->current_segment );
	Cell new_ver = this->current_segment->boundary().cell_in_front_of
		( this_center, tag::surely_exists );
	if ( new_ver.core == this->last_vertex )
		this->current_segment = nullptr;
	else
		this->current_segment = m.cell_in_front_of ( new_ver, tag::may_not_exist ).core;     }

// the above could get simpler and faster by avoiding the use of wrappers


Cell CellIterator::Around::OneCell::OfCodimTwo::OverSegments::NormalOrder::AsTheyAre::deref ( )
// virtual from CellIterator::Core

{ return Cell ( tag::whose_core_is, this->current_segment,
                tag::previously_existing, tag::surely_not_null );  }


Cell CellIterator::Around::OneCell::OfCodimTwo::OverSegments::ReverseOrder::AsTheyAre::deref ( )
// virtual from CellIterator::Core

{ return Cell ( tag::whose_core_is, this->current_segment,
                tag::previously_existing, tag::surely_not_null );  }


bool CellIterator::Around::OneCell::OfCodimTwo::OverSegments::in_range ( )
// virtual from CellIterator::Core

{	return this->current_segment;  }

	// for an open chain, last_vertex stays nullptr
	// for a closed loop, it is used as stopping criterion

//-----------------------------------------------------------------------------------------




