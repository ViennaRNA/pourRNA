/*
 * WalkGradientHashed.h
 *
 *  Created on: 14.08.2014
 *      Author: Gregor Entzian
 */

#ifndef WALKGRADIENTHASHED_H_
#define WALKGRADIENTHASHED_H_

#include "StructureUtils.h"
#include "GlobalParameter.h"
#include <stdlib.h>
#include <cstring>
extern "C" {
#include <ViennaRNA/move_set.h>
#include <ViennaRNA/pair_mat.h>
#include <ViennaRNA/fold.h>
}
#include "MyState.h"
#include "SpookyHash/SpookyHashMap.h"

/**
 * ! A gradient walk based on the ViennaRNA browseNeighbor function,
 * which is faster than the ell-functions.
 */
class WalkGradientHashed
{
public:
  /**
   * ! Initialize the gradient walk.
   * @param sequence the rna-sequence (ACGT).
   * @param maxNeighbors the maximal number of neighbors which can occur (neighbors of open chain).
   * @param maxHashSize is the maximal number of states which will be stored.
   */
  WalkGradientHashed (std::string sequence, size_t maxNeighbors, const size_t maxHashSize = 10000);
  /**
   * ! Clean up.
   */
  virtual
  ~WalkGradientHashed ();
  /**
   * ! Start a walk.
   * @param rnaSequence the rna sequence with characters ACGT.
   * @param maxNeighbors the maximal number of neighbors which can occur (neighbors of open chain).
   * @param startState with the initial structure and energy.
   * @param state2min a limited hash for the states.
   */
  MyState*
  walkGradient (char * rnaSequence, const MyState& startState, SpookyHashMap::HashTable & state2min);
  /**
   * ! Start a walk by calling the walkGradient function.
   * @param rnaSequence the rna sequence with characters ACGT.
   * @param startState with the initial structure and energy.
   */
  MyState *
  walk (char * rnaSequence, const MyState& startState);
private:
  // ! the maximal number of neighbors which can occur (neighbors of open chain).
  size_t MaxNeighbors;
  // ! to store the minimum for each visited state.
  SpookyHashMap::HashTable State2min;
  // ! energy parameter for the ViennaRNA functions (energy_of_..., browseNeighs)
  paramT & EnergyParameter;

  // ! the list of neighbors for one thread.
  neighborList NeighborList;
  // ! the encoded sequence. (used by some energy-functions)
  short *S0;
  // ! another encoded sequence. (used by some energy-functions)
  short *S1;

};

#endif /* WALKGRADIENTHASHED_H_ */
