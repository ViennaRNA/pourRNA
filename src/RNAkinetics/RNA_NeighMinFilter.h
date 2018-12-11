/*
 * RNA_NeighMinFilter.h
 *
 *  Created on: 09.08.2014
 *      Author: Hamdi
 *        This version is copied and adjusted by Gregor Entzian
 */

#ifndef RNA_NEIGHMINFILTER_H_
#define RNA_NEIGHMINFILTER_H_

#include "../SC_PartitionFunction.h"
#include <unordered_map>
extern "C" {
#include <ViennaRNA/move_set.h>
}
#include <algorithm>    // std::sort
#include <map>
#include <vector>
#include "../StateCollector.h"
#include "iostream"
#include "../PairHashTable.h"
//##############################################################################
//##############################################################################ظ

/*!
 * index pair
 */
//typedef std::pair<size_t, size_t> PairID;

/*!
 * This is an abstract interface implemented by Neighbor minimum Filtering
 * classes
 *
 * @author Hamdi Aldaoos
 */
class NeighMinFilter
{
protected:
//! set of local minima found so far, which  the indices to filter is subset of it
//  const std::map<size_t, MyState> & minimaSet;

public:
//! defines a list of local minima indices
typedef std::vector<size_t> IndexList;

/*!
 * sets the central minimum information
 * @param minimaSet the set of minima where the indices to filter belong
 *           to; ie. all indices to filter have to less than
 *           minimaSet.getMinCount();
 */
NeighMinFilter ()
{
}


/*!
 * destructor
 */
virtual
~NeighMinFilter ()
{
}


/*!
 * Filters a given list of minima indices and removes filtered indices
 * from the given container
 *
 * @param minimaSet, set of Local minimas found so far.
 * @param Z, container of partition functions, that already generated, which
 *  Z[i,j] holds partition function of saddle points beteween minima i and j
 *  Z[i,i] holds basins partition function.
 * @param curMinIndex the index of the current minimum the others are
 *         neighbored to
 * @param indicesToFilter the list of neighbored indices to prune
 *
 */
virtual
void
filter(const PairHashTable::HashTable &       minimaSet,
       const SC_PartitionFunction::Z_Matrix & Z,
       const size_t                           curMinIndex,
       const MyState                          currentMin,
       IndexList &                            indicesToFilter) = 0;
};
// end class

//##############################################################################

/*!
 * NeighMinCombinator Class implements NeighMinFilter, which takes
 * two NeighMinFilters and applys them one after the another.
 */

class NeighMinCombinator : public NeighMinFilter
{
protected:

//! represents first filtering technique to be applied
NeighMinFilter &  firstFilter;

//! represents second filtering technique to be applied
NeighMinFilter &  secondFilter;

public:

/*! constructor of NeighborCombinator Object
 * @param first is an Object of NeighMinFilter as first filter
 *  initialize the minimaSet of the Super Class  NeighMinFilter
 * @param second is an Object of NeighMinFilter as second filter
 *
 */
NeighMinCombinator (NeighMinFilter &  first,
                    NeighMinFilter &  second) :
  firstFilter(first), secondFilter(second)
{
}


/*!
 * destruction
 */
virtual
~NeighMinCombinator ()
{
}


/*!
 * Filters a given list of minima indices by a successive application of
 * two internally accessible filters
 *
 * @param minimaSet, set of Local minimas found so far.
 * @param Z, container of partition functions, that already generated, which
 *  Z[i,j] holds partition function of saddle points beteween minima i and j
 *  Z[i,i] holds basins partition function.
 * @param curMinIndex the index of the current minimum the others are
 *         neighbored to
 * @param indicesToFilter the list of neighbored indices to prune
 *
 */
virtual
void
filter(const PairHashTable::HashTable &       minimaSet,
       const SC_PartitionFunction::Z_Matrix & Z,
       const size_t                           curMinIndex,
       const MyState                          currentMin,
       IndexList &                            indicesToFilter);
};
// end of class

//##############################################################################

/*! This Class implements  criteria of sort the set of
 * Neighbors indices of the given State, through creating an Object of type
 * GreatByZMatrix, which will call third parameter Operator () of std::sort  .
 *
 */
class GreaterByZMatrix
{
private:
//! the Z matrix to use for order decisions
const SC_PartitionFunction::Z_Matrix &  Z;
//! index of current row in Z to access
const std::uint32_t                            curRow;

public:

/*! constructor
 * @param Zmatrix the Z values to use for comparison
 * @param row the row in the matrix to use
 */
GreaterByZMatrix (const SC_PartitionFunction::Z_Matrix &Zmatrix,
                  const std::uint32_t                          row) :
  Z(Zmatrix), curRow(row)
{
}


/*!
 * less comparison of two minima indices.
 * @param x the index to be checked if Z[curRow,x] is smaller
 * @param y the  to compare to (Z[curRow,y])
 * @return true if Z[curRow,x] > Z[curRow,y];
 *         false otherwise
 */
bool
operator ()(const std::uint32_t  x,
            const std::uint32_t  y)
{
  auto c_to_x = SC_PartitionFunction::PairID((std::uint32_t)curRow, (std::uint32_t)x);
  auto c_to_y = SC_PartitionFunction::PairID((std::uint32_t)curRow, (std::uint32_t)y);
  auto pf_cx_it = ((SC_PartitionFunction::Z_Matrix &)Z).find(c_to_x);
  auto pf_cy_it = ((SC_PartitionFunction::Z_Matrix &)Z).find(c_to_y);
  double  Zx  = pf_cx_it->second.getZ();
  double  Zy  = pf_cy_it->second.getZ();

  return Zx > Zy;
}
};
// ends the Class

/*!
 * This class uses  NeighMinFilter class to make filtering for the Neighbor
 * minima set to return the K highest Energy states
 *
 * @author Hamdi Aldaoos
 */
class NeighMin_K_Filter : public NeighMinFilter

{
protected:

//! access to the minima set stored in the super class
//using NeighMinFilter::minimaSet;

//! transition Z matrix
//  const SC_PartitionFunction::Z_Matrix & Z;

//! maximal number of best neighbors to keep
const size_t k_Neigh;

public:

/*!
 * constructs a neighbored minima filter that prunes a list to the K
 * minima indices with highest transition Z value.
 *
 *  @param K, given integer value, indicated Kth best Neighbors of current state, it will be used
 *  in Kth best Neighbors filter method to prune the unrelevant Neighbors.
 */
NeighMin_K_Filter (const size_t K) :
  NeighMinFilter(), k_Neigh(K)
{
}


/*!
 * destruction
 */
virtual
~NeighMin_K_Filter ();

/*!
 * Filters a given list of minima indices and removes filtered indices
 * from the given container.
 *
 * This is done by reducing the list to the indices with the highest
 * K transition Z-values
 *
 * @param minimaSet, set of Local minimas found so far.
 * @param Z, container of partition functions, that already generated, which
 *  Z[i,j] holds partition function of saddle points beteween minima i and j
 *  Z[i,i] holds basins partition function.
 * @param curMinIndex the index of the current minimum the others are
 *         neighbored to
 * @param indicesToFilter the list of indices to prune
 *
 */
virtual
void
filter(const PairHashTable::HashTable &       minimaSet,
       const SC_PartitionFunction::Z_Matrix & Z,
       const size_t                           curMinIndex,
       const MyState                          currentMin,
       IndexList &                            indicesToFilter);
};
// end of the class

//##############################################################################

/*!
 * This class uses  NeighMinFilter class to make filtering for the Neighbor
 * minima set to return the Neighbor Minima, which have lower Energy then
 * current Minimum
 *
 * @author Hamdi Aldaoos
 */
class NeighMin_E_Filter : public NeighMinFilter

{
protected:

//! access to the minima set stored in the super class
// using NeighMinFilter::minimaSet;

//! transition Z matrix
// const SC_PartitionFunction::Z_Matrix & Z;

//! Given param to be added to Neighbors Energy
const double deltaE;

public:

/*!
 * constructs a neighbored minima filter that prunes a list to the sublist
 * minima indices with higher Energy value then current minima.
 *
 *  @param E, given double value, it will be used
 *  in Energy filter method to prune the unrelevant Neighbors.
 */
NeighMin_E_Filter (const double E) :
  NeighMinFilter(), deltaE(E)
{
}


/*!
 * destruction
 */
virtual
~NeighMin_E_Filter ();

/*!
 * Filters a given list of minima indices and removes filtered indices
 * from the given container.
 *
 * This is done by reducing the list to the indices with Energy (Z-Values) higher
 * then the given minima
 *
 * @param minimaSet, set of Local minimas found so far.
 * @param Z, container of partition functions, that already generated, which
 *  Z[i,j] holds partition function of saddle points beteween minima i and j
 *  Z[i,i] holds basins partition function.
 * @param curMinIndex the index of the current minimum the others are
 *         neighbored to
 * @param indicesToFilter the list of indices to prune
 *
 */
virtual
void
filter(const PairHashTable::HashTable &       minimaSet,
       const SC_PartitionFunction::Z_Matrix & Z,
       const size_t                           curMinIndex,
       const MyState                          currentMin,
       IndexList &                            indicesToFilter);
};
// ends of Class

//##############################################################################

#endif
