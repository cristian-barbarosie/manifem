
// metric-tree.h 2021.06.18

//   This is MetricTree, a tiny C++ library for hierarchical organization
//   of a cloud of points in a metric space.

//   Copyright 2020, 2021 Cristian Barbarosie cristian.barbarosie@gmail.com
//   https://github.com/cristian-barbarosie/MetricTree

//   MetricTree is free software: you can redistribute it and/or modify it
//   under the terms of the GNU Lesser General Public License as published
//   by the Free Software Foundation, either version 3 of the License
//   or (at your option) any later version.

//   MetricTree is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//   See the GNU Lesser General Public License for more details.

//   You should have received a copy of the GNU Lesser General Public License
//   along with MetricTree.  If not, see <https://www.gnu.org/licenses/>.


// a cloud, i.e. a set of points in a metric space, organized as a tree
// similar to m-tree, just not balanced
// triangular inequality is assumed
// no geometric or topologic dimension assumed
// so it's a sort of generalization of quad-trees and oct-trees

// each node of the tree corresponds to a point, leaves have no special status

// the tree is not balanced (just like a quad-tree isn't)

// each node has a 'rank' associated to it (an integer, possibly zero, possibly negative)
// the rank has no special meaning except that
//   the children of a node N have rank one unit less
//   the children of a node N are no farther than dist[rank[N]] from N
// rank zero nodes have no special status
// leaves may have any rank (positive, zero or negative)

// as a consequence, indirect descendants of a node N will be no farther than 'range[r]' from N
// where range[r] = dist[r] + dist[r-1] + dist[r-2] + ...  (an infinite sum)  where r = rank[N]

// if N is a node and P is another point in the cloud with dist(N,P) < dist[rank[N]]
//   this does not imply that P is a subaltern of N  (just that it could be)
// if N is a node and P is another point in the cloud with dist(N,P) < range[rank[N]]
//   this does not imply that P is an indirect subaltern of N  (just that it could be)
// in other words: domains overlap

// 'dist' should be a geometric sequence, that is,  dist[r+1] = ratio * dist[r]
// of course all 'dist' are positive
// 'ratio' must be greater than 2 (we recommend a value between 5 and 10)
// so  range[r] = dist[r] / ( 1 - 1/ratio )  (the infinite sum above)

// we prefer to work with squared distance (thus avoiding computing square roots)

// see paragraph 12.10 in the manual of maniFEM
// https://webpages.ciencias.ulisboa.pt/~cabarbarosie/manifem/manual-manifem.pdf

#include <iostream>
#include <fstream>
#include <list>
#include <set>
#include <vector>
#include "math.h"
#include "assert.h"


template < typename Point, typename SqDist >
class MetricTree

{	public:
	
	class Node;

	SqDist squared_distance;
	// callable object returning the square of the distance between two points

	const double ratio;
	const double sq_ratio { ratio * ratio };
	const double range_factor { 0.999 - 1./ratio };
	// indirect range = dist[rank] / range_factor
	// should be  1 - 1/ratio  but we decrease it a little to compensate for numerical errors
	const double sq_range_factor { range_factor * range_factor };

	const double dist_rank_zero;
	const double sq_dist_rank_zero;
	
	Node * root { nullptr };

	// two vectors holding distances, one for nodes of positive rank and one for negative rank
	// to obtain ranges (distance to indirect descendants) just divide by 'range_factor'
	std::vector < double > dist_pos_rank, dist_neg_rank;
	std::vector < double > sq_dist_pos_rank, sq_dist_neg_rank;

	inline MetricTree ( SqDist sd, double d0, double r );

	Node * add ( const Point & );
	void add ( Node * );

	inline void remove ( Node * );

	inline std::list < Point > find_close_neighbours_of ( const Point & P, double dd );
	// return all points in the cloud whose distance to P is less than or equal to dd

	inline double get_dist ( int r );
	inline double get_sq_dist ( int r );

	inline void register_rank ( int r );

	inline void promote_children_of ( Node * nod );

	inline size_t nb_of_nodes ()
	{	if ( root ) return root->nb_of_nodes();
		return 0;                               }

};  // end of  class MetricTree

//-----------------------------------------------------------------------------------------------//

template < typename Point, typename SqDist >
inline MetricTree<Point,SqDist>::MetricTree ( SqDist sd, double d0, double r )
:	squared_distance ( sd ), ratio { r }, dist_rank_zero { d0 }, sq_dist_rank_zero { d0 * d0 },
	dist_pos_rank ( 1, dist_rank_zero ), dist_neg_rank ( 1, dist_rank_zero ),
	sq_dist_pos_rank ( 1, sq_dist_rank_zero ), sq_dist_neg_rank ( 1, sq_dist_rank_zero )
{	}

//-----------------------------------------------------------------------------------------------//

template < typename Point, typename SqDist >
inline double MetricTree<Point,SqDist>::get_dist ( int r )
{	if ( r == 0 ) return this->dist_rank_zero;
	if ( r > 0 )
	{	size_t rr = r;
		assert ( this->dist_pos_rank.size() > rr );
		return this->dist_pos_rank[rr];             }
	size_t rr = -r;
	assert ( this->dist_neg_rank.size() > rr );
	return this->dist_neg_rank[rr];                    }

//-----------------------------------------------------------------------------------------------//

template < typename Point, typename SqDist >
inline double MetricTree<Point,SqDist>::get_sq_dist ( int r )
{	if ( r == 0 ) return this->sq_dist_rank_zero;
	if ( r > 0 )
	{	size_t rr = r;
		assert ( this->sq_dist_pos_rank.size() > rr );
		return this->sq_dist_pos_rank[rr];             }
	size_t rr = -r;
	assert ( this->sq_dist_neg_rank.size() > rr );
	return this->sq_dist_neg_rank[rr];                    }

//-----------------------------------------------------------------------------------------------//

template < typename Point, typename SqDist >
inline void MetricTree<Point,SqDist>::register_rank ( int r )
{	if ( r == 0 ) return;
	if ( r > 0 )
	{	size_t rr = r;
		while ( this->sq_dist_pos_rank.size() <= rr )
		{	this->dist_pos_rank.push_back ( this->dist_pos_rank.back() * this->ratio );
			this->sq_dist_pos_rank.push_back ( this->sq_dist_pos_rank.back() * this->sq_ratio );  }  }
	else
	{	size_t rr = -r;
		while ( this->sq_dist_neg_rank.size() <= rr )
		{	this->dist_neg_rank.push_back ( this->dist_neg_rank.back() / this->ratio );
			this->sq_dist_neg_rank.push_back ( this->sq_dist_neg_rank.back() / this->sq_ratio );  }  }   }

//-----------------------------------------------------------------------------------------------//

template < typename Point, typename SqDist >
class MetricTree<Point,SqDist>::Node

{	public:

	Point point;  // point in the metric space

	int rank;

	Node * parent { nullptr };
	typename std::list<Node*>::iterator loc_in_parents_list;
	std::list < Node * > children;

	inline Node ( Point P, int r )
	:	point ( P ), rank { r }, children { }
	{ }

	void remove_from ( MetricTree * cloud );
	
	void get_close_neighbours_of
	( const Point & P, double dd, std::list < Point > & ll, MetricTree<Point,SqDist> * cloud );
	// return all points in the cloud whose distance to P is less than or equal to dd
	// the cloud is used as source of information (ratio etc)
	
	void promote ( MetricTree * cloud );
	// increases rank by one
	// the cloud is used as source of information (ratio etc) and is modified
	
	void adopt ( Node * nod, MetricTree * cloud );
	inline void raw_adopt ( Node * nod, MetricTree * cloud );
	void adopt_children_of ( Node * nod, MetricTree * cloud );
	// the cloud is used as source of information (ratio etc) and is modified

	size_t nb_of_nodes ();

};  // end of  class MetricTree::Node

//-----------------------------------------------------------------------------------------------//

namespace {
double cloud_power ( double x, int exp )
{	if ( exp == 0 )  return 1.;
	if ( exp < 0 )  return 1. / cloud_power ( x, -exp );
	if ( exp == 1 )  return x;
	int e = exp / 2;
	return cloud_power ( x, e ) * cloud_power ( x, exp - e );   }
}  //  end of anonymous namespace

//-----------------------------------------------------------------------------------------------//

template < typename Point, typename SqDist >
void MetricTree<Point,SqDist>::Node::promote ( MetricTree<Point,SqDist> * cloud )

// increase the rank by one

{	// promote children first
	for ( typename std::list<Node*>::const_iterator it_a = this->children.begin();
	      it_a != this->children.end(); it_a++                                       )
	{	typename MetricTree<Point,SqDist>::Node * ambitious = *it_a;
		if ( ambitious->rank == this->rank ) continue;  // already promoted
		ambitious->promote ( cloud );
		// transfer ownership of some brothers from 'this' to 'ambitious'
		for ( typename std::list<Node*>::const_iterator it_b = this->children.begin();
		      it_b != this->children.end();                                              )
		{	typename MetricTree<Point,SqDist>::Node * brother = *it_b;
			if ( brother->rank < this->rank )
			if ( cloud->squared_distance ( ambitious->point, brother->point ) <=
		       cloud->get_sq_dist ( ambitious->rank )                           )
			{	// transfer 'brother' from 'this' to 'ambitious'
				assert  ( it_b == brother->loc_in_parents_list );
				it_b = this->children.erase ( it_b );
				brother->parent = ambitious;
				ambitious->children.push_front ( brother );
				brother->loc_in_parents_list = ambitious->children.begin();
				continue;                                                    }
			it_b++;                                                                  }  }
	this->rank++;
}

//-----------------------------------------------------------------------------------------------//

template < typename Point, typename SqDist >
typename MetricTree<Point,SqDist>::Node * MetricTree<Point,SqDist>::add ( const Point & P )

// returns the newly created node

{	MetricTree<Point,SqDist>::Node * N = new MetricTree<Point,SqDist>::Node ( P, 0 );
	// rank is irrelevant, we give zero, will be set correctly soon
	this->add(N);  return N;                                                           }

//-----------------------------------------------------------------------------------------------//

template < typename Point, typename SqDist >
void MetricTree<Point,SqDist>::add ( MetricTree<Point,SqDist>::Node * nod )

{	if ( this->root == nullptr )
	{	this->root = nod;  return;  }
	this->root->adopt ( nod, this );
	if ( nod->parent ) return;  // has parent means has been adopted
	// nod cannot be adopted by this->root unless we increase its rank
	this->root->promote ( this );
	this->register_rank ( this->root->rank );
	this->add ( nod );                             }

//-----------------------------------------------------------------------------------------------//

template < typename Point, typename SqDist >
void MetricTree<Point,SqDist>::Node::adopt
(	MetricTree<Point,SqDist>::Node * nod, MetricTree<Point,SqDist> * cloud )

// 'this' tries to adopt 'nod' (first ask children of 'this' to adopt 'nod')
	
{	double sq_dist = cloud->get_sq_dist ( this->rank );
	double sq_range = sq_dist / cloud->sq_range_factor;
	double sq_d = cloud->squared_distance ( this->point, nod->point );
	if ( sq_d > sq_range )  return;  // 'nod' is out of the range of 'this'
	for ( typename std::list<Node*>::const_iterator it = this->children.begin();
	      it != this->children.end(); it++                                       )
	{	typename MetricTree<Point,SqDist>::Node * child = *it;
		assert ( child );
		child->adopt ( nod, cloud );
		if ( nod->parent ) return;                              }
		// has parent means has been adopted
	if ( sq_d <= sq_dist )  // yes, 'this' will adopt 'nod'
		this->raw_adopt ( nod, cloud );                                                 }
	// else ... sorry, nobody wants noddy :-(

//-----------------------------------------------------------------------------------------------//

template < typename Point, typename SqDist >
inline void MetricTree<Point,SqDist>::Node::raw_adopt
(	typename MetricTree<Point,SqDist>::Node * nod, MetricTree<Point,SqDist> * cloud )

{	cloud->register_rank ( this->rank - 1 );
	nod->rank = this->rank - 1;
	this->children.push_front ( nod );
	nod->parent = this;
	nod->loc_in_parents_list = this->children.begin();  }

//-----------------------------------------------------------------------------------------------//

template < typename Point, typename SqDist >
inline void MetricTree<Point,SqDist>::promote_children_of
( typename MetricTree<Point,SqDist>::Node * nod )

// 'nod' is in the process of being removed from the cloud so its rank is irrelevant
// and not necessarily correctly related to the rank of its children
// (children may have been promoted in the process)

{	for ( typename std::list<typename MetricTree<Point,SqDist>::Node*>::
	      const_iterator it = nod->children.begin(); it != nod->children.end(); )
	{	typename MetricTree<Point,SqDist>::Node * child = *it;
		child->promote ( this );
		// 'child' will attempt to adopt its own brothers now
		typename std::list<typename MetricTree<Point,SqDist>::Node*>::
			const_iterator itt = it;
		for ( itt++; itt != nod->children.end(); )
		{	typename MetricTree<Point,SqDist>::Node * little_brother = *itt;
			double sq_dist = this->get_sq_dist ( child->rank );
			double d = this->squared_distance ( child->point, little_brother->point );
			if ( d <= sq_dist )
			{	child->raw_adopt ( little_brother, this );
				itt = nod->children.erase ( itt );          }
			else itt++;                                                                    }
		it++;                                                                               }  }

//-----------------------------------------------------------------------------------------------//

template < typename Point, typename SqDist >
void MetricTree<Point,SqDist>::Node::adopt_children_of
( MetricTree<Point,SqDist>::Node * nod, MetricTree<Point,SqDist> * cloud )

// well, try to adopt some of them
// unlike in method 'adopt', do not ask children of 'this' to adopt

// 'nod' is in the process of being removed from the cloud so it rank is irrelevant
// and not necessarily correctly related to the rank of its children
// (children may have been promoted in the process)

{	for ( typename std::list<Node*>::const_iterator it = nod->children.begin();
	      it != nod->children.end();                                            )
	{	MetricTree<Point,SqDist>::Node * child = *it;
		assert ( child->loc_in_parents_list == it );
		double sq_dist = cloud->get_sq_dist ( this->rank );
		double d = cloud->squared_distance ( this->point, child->point );
		if ( d <= sq_dist )
		{	this->raw_adopt ( child, cloud );
			it = nod->children.erase ( it );  }
		else it++;                                                                     }  }

//-----------------------------------------------------------------------------------------------//

template < typename Point, typename SqDist >
inline void MetricTree<Point,SqDist>::remove
( MetricTree<Point,SqDist>::Node * nod )

{	nod->remove_from ( this );  }

//-----------------------------------------------------------------------------------------------//

template < typename Point, typename SqDist >
void MetricTree<Point,SqDist>::Node::remove_from
( MetricTree<Point,SqDist> * cloud )

{	MetricTree<Point,SqDist>::Node * p = this->parent;
	if ( p )
	{	typename std::list<Node*>::const_iterator it_p = this->loc_in_parents_list;
		assert ( it_p != p->children.end() );
		p->children.erase ( it_p );
		if ( this->children.empty() )  {  delete this;  return;  }
		while ( true )
		// p is available for adopting the children of 'this'
		// rank[p] == rank[this] + 1 == rank[children] + 2
		// but before that, other children of p may accept to adopt these children
		{	for ( typename std::list<Node*>::const_reverse_iterator it = p->children.rbegin();
		        it != p->children.rend(); it++                                               )
			// we use reverse iterator because raw_adopt inserts new children at the beginning
			{	(*it)->adopt_children_of ( this, cloud );
				if ( this->children.empty() )  {  delete this;  return;  }  }
			cloud->promote_children_of ( this );
			if ( p->parent == nullptr ) break;
			p = p->parent;                                                                        }  }
	else  // p == nullptr
	{	assert ( this == cloud->root );
		typename std::list<Node*>::const_iterator it = this->children.begin();
		if ( it == this->children.end() )  // cloud becomes empty
		{	cloud->root = nullptr;  delete this;  return;  }
		MetricTree<Point,SqDist>::Node * child = *it;
		child->parent = nullptr;
		cloud->root = child;
		assert ( child->loc_in_parents_list == it );
		it = this->children.erase ( it );
		if ( it == this->children.end() )  // we're good, 'child' will be the new root
		{	delete this;  return;  }
		child->promote ( cloud );  p = child;                                      }

	// p is available for adopting the children of 'this'
	// rank[p] == rank[children] + 1

	// basta promover p e os children um passo de cada vez
	assert ( p == cloud->root );
	if ( this->children.empty() )  {  delete this;  return;  }
	while ( true )
	{	p->adopt_children_of ( this, cloud );
		if ( this->children.empty() )  {  delete this;  return;  }
		p->promote ( cloud );  // register new rank !
		cloud->register_rank ( p->rank );
		cloud->promote_children_of ( this );                         }                                }

//-----------------------------------------------------------------------------------------------//

template < typename Point, typename SqDist >
inline std::list < Point > MetricTree<Point,SqDist>::find_close_neighbours_of
( const Point & P, double dd )

// return all points in the cloud whose distance to P is less than or equal to dd

{	std::list < Point > ll;
	this->root->get_close_neighbours_of ( P, dd, ll, this );
	return ll;                                                 }

//-----------------------------------------------------------------------------------------------//


template < typename Point, typename SqDist >
void MetricTree<Point,SqDist>::Node::get_close_neighbours_of
( const Point & P, double dd, std::list < Point > & ll, MetricTree<Point,SqDist> * cloud )

// return all points in the cloud whose distance to P is less than or equal to dd

{	double dist = cloud->get_dist ( this->rank );
	double range = dist / cloud->range_factor;
	double sq_d = cloud->squared_distance ( P, this->point );
	double sum = range + dd;
	if ( sq_d > sum * sum ) return;  // P is too far from 'this'
	if ( sq_d <= dd * dd ) ll.push_back ( this->point );
	for ( typename std::list<Node*>::const_iterator it = this->children.begin();
	      it != this->children.end(); it++                                       )
		(*it)->get_close_neighbours_of ( P, dd, ll, cloud );                         }

//-----------------------------------------------------------------------------------------------//

template < typename Point, typename SqDist >
	size_t MetricTree<Point,SqDist>::Node::nb_of_nodes ( )

{	size_t res = 1;
	for ( typename std::list<Node*>::const_iterator it = this->children.begin();
	      it != this->children.end(); it++                                       )
		res += (*it)->nb_of_nodes();
	return res;                                                                     }

//-----------------------------------------------------------------------------------------------//

