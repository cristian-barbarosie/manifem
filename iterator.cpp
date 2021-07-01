
// iterator.cpp 2021.06.12

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
  const tag::DoNotBother &, const tag::ThisMeshIsPositive & )
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
  const tag::DoNotBother &, const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
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
  const tag::DoNotBother &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
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
  const tag::DoNotBother &, const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
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
  const tag::DoNotBother &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                             }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                             }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                             }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
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
  const tag::DoNotBother &, const tag::ThisMeshIsPositive &                 )
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
  const tag::DoNotBother &, const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
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
  const tag::DoNotBother &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
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
  const tag::DoNotBother &, const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
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
  const tag::DoNotBother &, const tag::ThisMeshIsPositive & )
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
  const tag::DoNotBother &, const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
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
  const tag::DoNotBother &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
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
  const tag::DoNotBother &, const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
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
  const tag::DoNotBother &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to require reversed vertices for a 1D mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBother &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                         )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to require reversed vertices for a 1D mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBother &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                         )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to require reversed vertices for a 1D mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBother &,
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
  const tag::DoNotBother &, const tag::ThisMeshIsPositive & )
{	return new CellIterator::Over::SegmentsOfConnectedOneDimMesh::
		NormalOrder::ReverseEachCell::AssumeCellsExist ( this );      }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBother &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                         )
{	return new CellIterator::Over::SegmentsOfConnectedOneDimMesh::
		NormalOrder::ReverseEachCell::AssumeCellsExist ( this );      }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBother &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                         )
{	return new CellIterator::Over::SegmentsOfConnectedOneDimMesh::
		ReverseOrder::ReverseEachCell::AssumeCellsExist ( this );      }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBother &,
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
  const tag::DoNotBother &, const tag::ThisMeshIsPositive &                    )
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
  const tag::ReverseEachCell &, const tag::DoNotBother &,
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
  const tag::ReverseEachCell &, const tag::DoNotBother &,
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
  const tag::ReverseEachCell &, const tag::DoNotBother &,
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
  const tag::DoNotBother &, const tag::ThisMeshIsPositive &       )
{	return new CellIterator::Over::SegmentsOfConnectedOneDimMesh::
		NormalOrder::ReverseEachCell::AssumeCellsExist ( this );      }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &, const tag::DoNotBother &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                             )
{	return new CellIterator::Over::SegmentsOfConnectedOneDimMesh::
		NormalOrder::ReverseEachCell::AssumeCellsExist ( this );      }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &, const tag::DoNotBother &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                             )
{	return new CellIterator::Over::SegmentsOfConnectedOneDimMesh::
		ReverseOrder::ReverseEachCell::AssumeCellsExist ( this );      }


CellIterator::Core * Mesh::Connected::OneDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &, const tag::DoNotBother &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &                        )
{	return new CellIterator::Over::SegmentsOfConnectedOneDimMesh::
		ReverseOrder::ReverseEachCell::AssumeCellsExist ( this );      }

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


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
  const tag::DoNotBother &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to require reversed vertices for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                                   }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBother &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                         )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "This makes no sense for two different reasons." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBother &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                         )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "This makes no sense for two different reasons." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBother &,
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
  const tag::DoNotBother &, const tag::ThisMeshIsPositive & )
{	assert ( this->get_dim_plus_one() == 2 );
	return new CellIterator::Over::CellsOfFuzzyMesh::ReverseEachCell::AssumeCellsExist
		( this, tag::cells_of_dim, 1 );                                                  }
// perhaps provide a different constructor for calling directly msh->cells.back(),
// it should be slightly faster when compiled with optimization


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBother &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                         )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to require order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBother &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                         )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to reverse order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBother &,
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
  const tag::DoNotBother &, const tag::ThisMeshIsPositive &                    )
{	assert ( this->get_dim_plus_one() == d+1 );
	return new CellIterator::Over::CellsOfFuzzyMesh::ReverseEachCell::AssumeCellsExist
		( this, tag::cells_of_dim, d );                                                  }
// perhaps provide a different constructor for calling directly msh->cells.back(),
// it should be slightly faster when compiled with optimization


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d,
  const tag::ReverseEachCell &, const tag::DoNotBother &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to require order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d,
  const tag::ReverseEachCell &, const tag::DoNotBother &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to reverse order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d,
  const tag::ReverseEachCell &, const tag::DoNotBother &,
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
  const tag::DoNotBother &, const tag::ThisMeshIsPositive &       )
{	return new CellIterator::Over::CellsOfFuzzyMesh::ReverseEachCell::AssumeCellsExist
		( this, tag::cells_of_dim, tag::Util::assert_diff ( this->get_dim_plus_one(), 1 ) );  }
// perhaps provide a different constructor for calling directly msh->cells.back(),
// it should be slightly faster when compiled with optimization


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &, const tag::DoNotBother &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                             )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to require order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &, const tag::DoNotBother &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                             )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to reverse order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &, const tag::DoNotBother &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &                        )
{	return new CellIterator::Over::CellsOfFuzzyMesh::ReverseEachCell::AssumeCellsExist
		( this, tag::cells_of_dim, tag::Util::assert_diff ( this->get_dim_plus_one(), 1 ) );  }
// perhaps provide a different constructor for calling directly msh->cells.back(),
// it should be slightly faster when compiled with optimization

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
               ::ReverseOrder::ReverseEachCell::AssumeCellsExist::deref ( )
// virtual from CellIterator::Core
// iterate over the two vertices, first tip (positive) then base (negative)
{	assert ( this->seg_p );
	assert ( this->passage < 2 );
	if ( this->passage == 0 )
		return this->seg_p->tip_attr.reverse ( tag::surely_exists );      // negative vertex
	else
		return this->seg_p->base_attr.reverse ( tag::surely_exists );  }  // positive vertex

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
	if ( this->msh->first_ver ==
	     this->msh->last_ver.reverse(tag::surely_exists) ) // loop
	{	Cell seg = this->current_vertex->cell_behind_within[msh];
		this->current_vertex = tag::Util::assert_cast
			< Cell::Core*, Cell::Positive::Vertex* >
			( seg.base().reverse(tag::surely_exists).core );             	}
	assert ( this->current_vertex->is_positive() );
	assert ( this->last_vertex->is_positive() );                                  }


void CellIterator::Over::VerticesOfConnectedOneDimMesh::NormalOrder::reset
( const tag::StartAt &, Cell::Core * ver )  // virtual from CellIterator::Core
{	assert ( this->msh );
	assert ( this->msh->nb_of_segs > 0 );
	assert ( ver );
	assert ( ver->get_dim() == 0 );
	this->current_vertex = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( ver );
	if ( this->msh->first_ver ==
       this->msh->last_ver.reverse(tag::surely_exists) ) // loop
	{	Cell seg = ver->cell_behind_within[msh];
		this->last_vertex = tag::Util::assert_cast
			< Cell::Core*, Cell::Positive::Vertex* >
			( seg.base().reverse(tag::surely_exists).core );    }
	else // open chain
		this->last_vertex = this->msh->last_ver.core;
	assert ( this->current_vertex->is_positive() );
	assert ( this->last_vertex->is_positive() );                         }


void CellIterator::Over::VerticesOfConnectedOneDimMesh::ReverseOrder::reset
( const tag::StartAt &, Cell::Core * ver )  // virtual from CellIterator::Core
{	assert ( this->msh );
	assert ( this->msh->nb_of_segs > 0 );
	assert ( ver );
	assert ( ver->get_dim() == 0 );
	this->current_vertex = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( ver );
	if ( this->msh->first_ver ==
			 this->msh->last_ver.reverse(tag::surely_exists) ) // loop
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
	if ( this->msh->first_ver ==
			 this->msh->last_ver.reverse(tag::surely_exists) ) // loop
		this->last_vertex = seg->base().reverse(tag::surely_exists).core;
	else  // open chain
		this->last_vertex = this->msh->last_ver.core;
	assert ( this->last_vertex->is_positive() );                         }


void CellIterator::Over::SegmentsOfConnectedOneDimMesh::ReverseOrder::reset
( const tag::StartAt &, Cell::Core * seg )  // virtual from CellIterator::Core
{	assert ( this->msh );  assert ( this->msh->nb_of_segs > 0 );
	assert ( seg );  assert ( seg->get_dim() == 1 );
	this->current_segment = seg;
	if ( this->msh->first_ver ==
			 this->msh->last_ver.reverse(tag::surely_exists) ) // loop
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
	ReverseOrder::ReverseEachCell::AssumeCellsExist::deref ( )
// virtual from CellIterator::Core
{	assert ( this->current_segment );
	assert ( this->current_segment->reverse_attr.exists() );
	return this->current_segment->reverse_attr;              }


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


	
