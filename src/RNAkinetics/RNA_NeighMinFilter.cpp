/*
 * RNA_NeighMinFilter.cpp
 *
 *  Created on: 09.08.2014
 *      Author: Hamdi
 *        This version is copied and adjusted by Gregor Entzian
 */

#include "RNA_NeighMinFilter.h"
#include <unordered_map>
#include "../StructureUtils.h"

//################################################################################################################################

/*
 * Destruction
 */
NeighMin_K_Filter::~NeighMin_K_Filter ()
{
}


/*Filter the set of Neighbors Indices of the current State
 * @param curMinIndex, the index of the current State in the MinimaSet
 * @param indicesToFilter, the set of Neighbors indices of current State
 *   This Method Will reduce the Neighbor Indices set.
 * */
void
NeighMin_K_Filter::filter(const PairHashTable::HashTable &        minimaSet,
                          const SC_PartitionFunction::Z_Matrix &  Z,
                          const size_t                            curMinIndex,
                          const MyState                           currentMin,
                          IndexList &                             indicesToFilter)
{
  // Sort by decreasing Z (= decreasing rate). Neighbors with equal Z are ordered by their
  // minimum (lower energy, then lexicographically smaller structure). Sorting by Z alone left
  // the order of ties to the order of indicesToFilter, which comes from a hash map and from
  // thread scheduling, so which of several equally good neighbors was kept was arbitrary.
  std::unordered_map<size_t, const MyState *> stateOfIndex;
  for (auto it = minimaSet.begin(); it != minimaSet.end(); it++)
    stateOfIndex[it->second] = &it->first;

  auto zOf = [&Z, curMinIndex](size_t x) {
               return Z.find(SC_PartitionFunction::PairID(curMinIndex, x))->second.getZ();
             };
  std::sort(indicesToFilter.begin(), indicesToFilter.end(),
            [&](size_t x, size_t y) {
    double zx = zOf(x), zy = zOf(y);
    if (zx != zy)
      return zx > zy;

    return *stateOfIndex.at(x) < *stateOfIndex.at(y);
  });

  //consider two cases for K_Neigh compare to number of Neighbors  of the given state,
  //remove from list all that are not among the best K
  if (k_Neigh <= indicesToFilter.size())
    // resize to the first k_Neigh elements and remove all others
    indicesToFilter.resize(k_Neigh);
} // ends of Filter Method


//################################################################################################################################

/*
 * Destruction
 * */

NeighMin_E_Filter::~NeighMin_E_Filter ()
{
}


/* ! Filter the set of Neighbors Indices of the current State
 * ! @ param curMinIndex, the index of the current State in the MinimaSet
 * ! @ param indicesToFilter, the set of Neighbors indices of current State
 *   This Method Will reduce the Neighbor Indices set, depends in the Z energy
 *   of the Neighbors of the current State using deltaE value.
 */
void
NeighMin_E_Filter::filter(const PairHashTable::HashTable &        minimaSet,
                          const SC_PartitionFunction::Z_Matrix &  Z,
                          const size_t                            curMinIndex,
                          const MyState                           currentMin,
                          IndexList &                             indicesToFilter)
{
  // create the Temporary Neighbors indices, which will store the Neighbors
  // that pass the filter.
  IndexList tempList;

  // index -> minimum (the former linear search over minimaSet for every index was quadratic)
  std::unordered_map<size_t, const MyState *> stateOfIndex;
  for (auto it = minimaSet.begin(); it != minimaSet.end(); it++)
    stateOfIndex[it->second] = &it->first;

  // Keep a neighbor if it is energetically smaller than the current minimum raised by the
  // filter value: E(neigh) < E(cur) + filterValue, and at equality the lexicographically
  // smaller structure wins, as in the order used for gradient basins. With filterValue 0 this
  // keeps exactly the more stable neighbors; the former test E(neigh) - E(cur) <= filterValue
  // also kept neighbors of equal energy, so a macroscopic gradient walk could step sideways.
  // Energies are compared in integer 10 cal/mol.
  const long long threshold = StructureUtils::KcalToDcal(deltaE);
  for (IndexList::const_iterator i = indicesToFilter.begin(); i != indicesToFilter.end(); i++) {
    auto other = stateOfIndex.find(*i);
    if (other == stateOfIndex.end())
      continue;

    const MyState   *otherMin = other->second;
    const long long diff      = (long long)otherMin->energy - currentMin.energy;
    if (diff < threshold
        || (diff == threshold
            && StructureUtils::IsSmaller(otherMin->structure, currentMin.structure)))
      tempList.push_back(*i);
  }

  //resize the size of indicesToFilter List by take just the relevant elements, which already stored
  // at tempList.
  indicesToFilter.clear();
  indicesToFilter.assign(tempList.begin(), tempList.end());
}     // ends Filter method


//################################################################################################################################

/* ! Filter the set of Neighbors Indices of the current State
 * ! @ param curMinIndex, the index of the current State in the MinimaSet
 * ! @ param indicesToFilter, the set of Neighbors indices of current State
 *   This Method Will reduce the Neighbor Indices set, depends in two Filter Techniques
 *   which represented as firstfilter, secondfilter parameters
 */
void
NeighMinCombinator::filter(const PairHashTable::HashTable &       minimaSet,
                           const SC_PartitionFunction::Z_Matrix & Z,
                           const size_t                           curMinIndex,
                           const MyState                          currentMin,
                           IndexList &                            indicesToFilter)
{
  // apply both filters and be done

  firstFilter.filter(minimaSet, Z, curMinIndex, currentMin, indicesToFilter);
  secondFilter.filter(minimaSet, Z, curMinIndex, currentMin, indicesToFilter);
}
