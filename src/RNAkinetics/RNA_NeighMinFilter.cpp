/*
 * RNA_NeighMinFilter.cpp
 *
 *  Created on: 09.08.2014
 *      Author: Hamdi
 *      	This version is copied and adjusted by Gregor Entzian
 */

#include "RNA_NeighMinFilter.h"

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
NeighMin_K_Filter::filter (const PairHashTable::HashTable & minimaSet, const SC_PartitionFunction::Z_Matrix & Z,
			   const size_t curMinIndex, const MyState currentMin, IndexList & indicesToFilter)
{
  // creating an greaterByZ Object as third parameter of std:: sort method
  GreaterByZMatrix greaterByZ (Z, curMinIndex);

  // sort bz Z value
  std::sort (indicesToFilter.begin (), indicesToFilter.end (), greaterByZ);

  //consider two cases for K_Neigh compare to number of Neighbors  of the given state,
  //remove from list all that are not among the best K
  if (k_Neigh <= indicesToFilter.size ())
    {

      // resize to the first k_Neigh elements and remove all others
      indicesToFilter.resize (k_Neigh);
    }

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
NeighMin_E_Filter::filter (const PairHashTable::HashTable & minimaSet, const SC_PartitionFunction::Z_Matrix & Z,
			   const size_t curMinIndex, const MyState currentMin, IndexList & indicesToFilter)
{
  // create the Temporary Neighbors indices, which will store the Neighbors
  // with energy higher than current state.
  IndexList tempList;

  // compare the Energy of current minima to all its Neighbors
  for (IndexList::const_iterator i = indicesToFilter.begin (); i != indicesToFilter.end (); i++)
    {

      //find minimum with index i.
      MyState* otherMin = NULL;
      for (auto it = minimaSet.begin (); it != minimaSet.end (); it++)
	{
	  if (*i == it->second)
	    {
	      otherMin = (MyState*) &it->first;
	      break;
	    }
	}

      if (otherMin != NULL)
	{
	  // keep : E(curMin) >= (E(neighMin)-deltaE)
	  // remove E(curMin) < (E(neighMin)-deltaE)
	  if ((currentMin.energy / 100.0 >= (otherMin->energy / 100.0 - deltaE)))
	    {
	      tempList.push_back (*i);
	    }
	}
    }

  //resize the size of indicesToFilter List by take just the relevant elements, which already stored
  // at tempList.
  indicesToFilter.clear ();
  indicesToFilter.assign (tempList.begin (), tempList.end ());

}	    // ends Filter method

//################################################################################################################################

/* ! Filter the set of Neighbors Indices of the current State
 * ! @ param curMinIndex, the index of the current State in the MinimaSet
 * ! @ param indicesToFilter, the set of Neighbors indices of current State
 *   This Method Will reduce the Neighbor Indices set, depends in two Filter Techniques
 *   which represented as firstfilter, secondfilter parameters
 */
void
NeighMinCombinator::filter (const PairHashTable::HashTable & minimaSet, const SC_PartitionFunction::Z_Matrix & Z,
			    const size_t curMinIndex, const MyState currentMin, IndexList & indicesToFilter)
{
  // apply both filters and be done

  firstFilter.filter (minimaSet, Z, curMinIndex, currentMin, indicesToFilter);
  secondFilter.filter (minimaSet, Z, curMinIndex, currentMin, indicesToFilter);
}

