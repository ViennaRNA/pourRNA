/*
 * StatePairCollector.h
 *
 *  Created on: 26.07.2014
 *      Author: Gregor Entzian
 */

#ifndef STATEPAIRCOLLECTOR_H_
#define STATEPAIRCOLLECTOR_H_

#include <stdio.h>
extern "C"
{
#include <ViennaRNA/fold_vars.h>
}
#include "PairHashTable.h"
#include "unordered_map"
#include <vector>
#include "SC_PartitionFunction.h"
#include <algorithm>
#include <math.h>
#include "StructureUtils.h"
#include "WalkGradientHashed.h"
#include <list>

/**
 * ! This class calculates the partition function for the contact surface.
 *   The partition function will be updated, when a statepair is added.
 */
class StatePairCollector
{
public:
  /**
   * ! Collects all pairs of neighbors between the currentMinID and other states.
   * Calculates the partition function.
   * @param minima set of all minima.
   * @param currentMinID the id of the current minimum which represents the current gradient basin.
   * @param z a container for the partition functions.
   * @param discoveredMinima queue to report new discovered minima.
   */
  StatePairCollector (size_t currentMinID,
		      PairHashTable::HashTable& minima, SC_PartitionFunction::Z_Matrix& z,
		      const size_t maxGradWalkHashed,std::list<MyState> *discoveredMinima);
  virtual
  ~StatePairCollector ();
  /**
   * ! add the statepair to the partition function.
   * ! Assign a ID for the local minimum if it is not in the minima set, and add it.
   */
  virtual void
  add (vrna_fold_compound_t *vc, const MyState* const state1, const MyState* const state2);
  /**
   * ! do a gradient walk for the given structure.
   *   @structure a rna-structure as pair table.
   *   @return the minimal structure and energy.
   */
  static MyState *
  doGradientWalk (const std::string rnaSequence, short * structure);

  size_t
  getNumberOfOuterStates () const
  {
    return NumberOfOuterStates;
  }

private:
  // ! identifier of the current minimum.
  const size_t CurMinID;
  // ! minima map: structure, index.
  PairHashTable::HashTable& Minima;
  // ! the partition funcion matrix.
  SC_PartitionFunction::Z_Matrix& Z;
  // ! minimum identification for this state pair
  WalkGradientHashed GradWalk;
  // ! handle minima id for outer states to save gradient walks.
  PairHashTable::HashTable HandledOuterStates;
  // ! number of states which are not in basin, but in the contactsurface.
  size_t NumberOfOuterStates;
  // ! unique minima queue (ready to flood)
  std::list<MyState> *DiscoveredMinima;

};

#endif /* STATEPAIRCOLLECTOR_H_ */
