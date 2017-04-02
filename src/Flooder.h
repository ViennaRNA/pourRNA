/*
 * flooder.h
 *
 *  Created on: 26.07.2014
 *      Author: Gregor Entzian
 */

#ifndef FLOODER_H_
#define FLOODER_H_

#include "MyState.h"
#include <iostream>
#include <limits.h>
extern "C" {
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/move_set.h>
}
#include <list>
#include <unordered_map>
#include "StateCollector.h"
#include "StatePairCollector.h"
#include "PairHashTable.h"
#include "PriorityQueue.h"
#include "StructureUtils.h"

class Flooder
{
public:

  /**
   * ! local basin flooder.
   * @param sequence the rna-sequence as "acgu"-string.
   * @param maxEnergy the maximal energy all states should be below in kcal/mol
   * @param maxToQueue the maximal number of elements the underlying queue is allowed to contain
   */
  Flooder (std::string sequence, double maxEnergy, size_t maxToQueue);

  /**
   * ! flood the basin of the given local minimum. Compute the partition function sum.
   *  @param localMinState the minimum of the basin which should be flooded.
   *  @param scBasin a state collector which counts the state of the basin.
   *  @param scSurface a collector which calculates the contact surface.
   */
  int
  floodBasin (vrna_fold_compound_t *vc, const MyState& localMinState, StateCollector& scBasin, StatePairCollector& scSurface);
  virtual
  ~Flooder ();

private:
  // ! the rna-sequence in "acgu"-format.
  char* Sequence;
  //! the maximal energy all states should be below in 10cal/mol
  int MaxEnergy;
  //! the maximal number of elements the underlying queue is allowed to contain
  size_t MaxStatesToQueue;
  //! states that where considered during the flooding.
  size_t ProcessedStates;

public:
  /**
   * ! returns the finial energy threshold that has been changed if
   * states from the priority queue are removed (in kcal/mol).
   */
  double
  getFinalEnergyThreshold ()
  {
    return (double) MaxEnergy / 100.0;
  }
  ;
  /**
   * ! get the neighbored structures which differ in one base-pair.
   *   The neighbors will be stored in the NeighborList variable.
   *  @param structureEnergy the rna-structure in pair-table format and energy (see ViennaRNA-package).
   */
  struct_en*
  get_Neighbors_pt (vrna_fold_compound_t *vc, struct_en* structureEnergy);

  /**
   * ! states that where considered during the flooding.
   */
  size_t
  getProcessedStates () const
  {
    return ProcessedStates;
  }
};
#endif /* FLOODER_H_ */
