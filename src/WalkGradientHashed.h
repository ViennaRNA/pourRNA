/*
 * WalkGradientHashed.h
 *
 *  Created on: 14.08.2014
 *      Author: Gregor Entzian
 */

#ifndef WALKGRADIENTHASHED_H_
#define WALKGRADIENTHASHED_H_

#include "StructureUtils.h"
#include <stdlib.h>
#include <cstring>
extern "C" {
#include <ViennaRNA/neighbor.h>
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/eval.h>
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
  WalkGradientHashed (std::string sequence, const size_t maxHashSize = 10000);
  /**
   * ! Clean up.
   */
  virtual
  ~WalkGradientHashed ();

  struct_en*
  get_Neighbors_pt (vrna_fold_compound_t *vc, struct_en* structureEnergy);
  /**
   * ! Start a walk.
   * @param rnaSequence the rna sequence with characters ACGT.
   * @param maxNeighbors the maximal number of neighbors which can occur (neighbors of open chain).
   * @param startState with the initial structure and energy.
   * @param state2min a limited hash for the states.
   */
  MyState*
  walkGradient (vrna_fold_compound_t *vc,char * rnaSequence, const MyState& startState, SpookyHashMap::HashTable & state2min);
  /**
   * ! Start a walk by calling the walkGradient function.
   * @param rnaSequence the rna sequence with characters ACGT.
   * @param startState with the initial structure and energy.
   */
  MyState *
  walk (vrna_fold_compound_t *vc,char * rnaSequence, const MyState& startState);
private:
  // ! to store the minimum for each visited state.
  SpookyHashMap::HashTable State2min;
  // ! the new foldcompound object==> replaces Energy parameters TODO
  vrna_fold_compound_t *VC;

};

#endif /* WALKGRADIENTHASHED_H_ */
