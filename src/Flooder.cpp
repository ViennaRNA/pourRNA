/*
 * flooder.cpp
 *
 *  Created on: 26.07.2014
 *      Author: Gregor Entzian
 */

#include "Flooder.h"
extern "C" {
#include <ViennaRNA/neighbor.h>
#include <ViennaRNA/eval.h>
#include <ViennaRNA/structure_utils.h>
}

//energy*100 --> convert from kcal/mol to 10cal/mol (the ViennaRNA internal integer format for energies)
Flooder::Flooder(double       maxEnergy,
                 size_t       maxToQueue,
                 unsigned int move_set) :
  MaxEnergy(maxEnergy * 100.0), MaxStatesToQueue(
    maxToQueue), ProcessedStates(0), Move_set(move_set)
{
}


Flooder::~Flooder()
{
}


struct_en *
Flooder::get_Neighbors_pt(vrna_fold_compound_t  *vc,
                          const MyState         *structureEnergy,
                          size_t                *numberOfNeighbors)
{
  //browse_neighs_pt_par_list_alloc_energy (Sequence, structureEnergy, S0, S1, 0, 0, 0, &NeighborList, &EnergyParameter, MaxNeighbors);
  vrna_move_t *tmp_neighbors = vrna_neighbors(vc, structureEnergy->structure,
                                              this->Move_set);
  size_t      count = 0;

  for (vrna_move_t *m = tmp_neighbors; m->pos_5 != 0; m++)
    count++;
  *numberOfNeighbors = count;

  struct_en *neighbors  = (struct_en *)malloc(sizeof(struct_en) * (count + 1));
  size_t       i           = 0;
  int       energy;
  short     *pt_neighbor;
  struct_en *neighbor;
  for (vrna_move_t *m = tmp_neighbors; m->pos_5 != 0; m++, i++) {
    energy      = vrna_eval_move_shift_pt(vc, m, structureEnergy->structure); // vrna_eval_move_pt(vc, structureEnergy->structure, m->pos_5,m->pos_3);
    pt_neighbor = vrna_ptable_copy(structureEnergy->structure);
    vrna_move_apply(pt_neighbor, m);
    neighbor            = &neighbors[i];
    neighbor->energy    = structureEnergy->energy + energy;
    neighbor->structure = pt_neighbor;
  }
  neighbors[i].structure = NULL;
  free(tmp_neighbors);
  return neighbors;
}


int
Flooder::floodBasin(vrna_fold_compound_t  *vc,
                    const MyState&        localMinState,
                    StateCollector        *scBasin,
                    StatePairCollector    *scSurface)
{
  if (localMinState.structure == NULL)
    std::cout << "No local Min given! (structure==NULL)" << std::endl;

  if (scBasin == NULL)
    std::cout << "No StateCollector given! (scBasin==NULL)" << std::endl;

  if (localMinState.energy >= MaxEnergy)
    std::cout << "E(localMin) >= maxE" << std::endl;

  ProcessedStates = 0;

  // the priority queue to store neighbors to handle in
  PriorityQueue<QueueValue>               pq;
  // the hash to store seen elements of the basin
  //PairHashTable::HashTable handled;
  HashSet::UnorderedHashSet               handled;
  // init using the given local minimum
  PriorityQueue<QueueValue>::InsertResult it = pq.insert(localMinState);
  it.first->second.QueueState.structure = NULL; // has no minimal neighbor (is min)
  it.first->second.QueueState.energy    = it.first->first.QueueState.energy;
  it.first->second.StateID              = 0;

  typedef std::list<PriorityQueue<QueueValue>::QueueKey> NeighListPQ;

  // create temporary objects that are altered during the flooding
  MyState     *topState = NULL;
  MyState     *curNeigh = NULL;

  // contains all neighbors with E(top) < E <= maxE
  NeighListPQ toStore;

  // contains all neighbors with E < E(top)
  NeighListPQ toCheck;

  while (!pq.empty()) {
    // get top element
    PriorityQueue<QueueValue>::const_iterator top = pq.top();
    // reference to top State
    topState = (MyState *)&top->first.QueueState;

    // get neighbor list
    size_t                                    numberOfNeighbors = 0;
    struct_en                                 *neighbors        = get_Neighbors_pt(
      vc,
      topState,
      &
      numberOfNeighbors);

    short *minNeigh = topState->structure;  //init neighbor with minimal energy
    int   minNeighE = topState->energy;     //init energy of the minimal neighbor

    // iterate through all neighbors
    for (struct_en *neighbor = neighbors; neighbor->structure != NULL;
         neighbor++) {
      // check for lowest neighbor
      if ((neighbor->energy < minNeighE)
          || (neighbor->energy == minNeighE
              && StructureUtils::IsSmaller(neighbor->structure,
                                           minNeigh))) {
        minNeigh  = neighbor->structure;
        minNeighE = neighbor->energy;
      }

      // check if neighbor is of interest and not rejected
      if (neighbor->energy <= MaxEnergy) {
        // energy does not exceed maxE
        // check if neighbor is of higher order --> to PQ
        // energy higher than top or
        if (neighbor->energy > top->first.QueueState.energy
            || (neighbor->energy == top->first.QueueState.energy
                && StructureUtils::IsGreater(neighbor->structure,
                                             top->first.QueueState.structure))) {
          // greater in sequence order

          // store neighbor with higher energy that has to be added to pq
          toStore.push_front(MyState(neighbor->energy, neighbor->structure));
        } else if (scSurface != NULL) {
          // store all neighbors with lower energy
          toCheck.push_back(MyState(neighbor->energy, neighbor->structure));
        }
      }
    }

    // handle top depending on its basin membership
    if (top->second.QueueState.structure == NULL || (StructureUtils::IsEqual(minNeigh,
                                                                             top->second.QueueState
                                                                             .structure)) // belongs to the basin
        ) {
      // or top is the local minimum
      // store because top belongs to the basin
      scBasin->add(*topState);
      // add to hashed elements and mark that part of the basin
      handled.insert(top->first.QueueState);

      // put all neighbors with higher energy to PQ
      PriorityQueue<QueueValue>::iterator toFill;
      // sort NeighborList with increasing energy
      toStore.sort();

      for (NeighListPQ::const_iterator n = toStore.begin();
           n != toStore.end() && n->QueueState.energy <= MaxEnergy;
           n++) {
        toFill = pq.find(*n);
        if (toFill == pq.end()) {
          // insert empty value object
          toFill = pq.insert(n->QueueState).first;
          // feed data
          toFill->second.QueueState.energy =
            top->first.QueueState.energy;
          toFill->second.QueueState.structure = vrna_ptable_copy(
            top->first.QueueState.structure);
          toFill->second.StateID = top->second.StateID;
          // check queue size and reduce if necessary
          if (pq.size() >= MaxStatesToQueue) {
            // queue to big

            // delete last element (i.e. the element with the highest energy).
            PriorityQueue<QueueValue>::iterator toErase =
              --pq.end();
            pq.erase(toErase->first);

            // adjust flood level to current maximal E in queue.
            MaxEnergy = pq.getMaxE();
          }
        } else  // update of entry if this one is smaller
        if (toFill->second.QueueState.energy
            > top->first.QueueState.energy
            || (toFill->second.QueueState.energy
                == top->first.QueueState.energy
                && StructureUtils::IsGreater(
                  toFill->second.QueueState.structure,
                  top->first.QueueState.structure))) {
          // feed data
          toFill->second.QueueState.energy =
            top->first.QueueState.energy;
          toFill->second.QueueState.structure = vrna_ptable_copy(
            top->first.QueueState.structure);
          toFill->second.StateID = top->second.StateID;
        }
      }

      // check if surface / contact plane is of interest
      if (scSurface != NULL) {
        // check all neighbors with lower energy if NOT part of the
        // basin to enumerate the surface of the basin
        for (NeighListPQ::const_iterator n = toCheck.begin(); n != toCheck.end(); n++) {
          // check if neighbor is NOT part of basin
          if (handled.find(n->QueueState) == handled.end()) {
            // get neighbor
            curNeigh = (MyState *)&n->QueueState;
            // add to surface reporter
            scSurface->add(vc, topState, curNeigh, false);
          }
        }
      }
    } else {
      // belongs NOT to the basin but in surface
      // check if surface / contact plane is of interest
      if (scSurface != NULL) {
        // check all neighbors with lower energy IF part of the
        // basin to enumerate the surface of the basin
        for (NeighListPQ::const_iterator n = toCheck.begin(); n != toCheck.end(); n++) {
          // check if neighbor IS part of basin
          if (handled.find(n->QueueState) != handled.end()) {
            // get neighbor
            curNeigh = (MyState *)&n->QueueState;
            // add to surface reporter
            scSurface->add(vc, curNeigh, topState, true);
          }
        }
      }
    }

    // remove top element
    pq.pop();
    ProcessedStates++;
    toStore.clear();
    toCheck.clear();

    for (struct_en *neighbor = neighbors; neighbor->structure != NULL;
         neighbor++)
      free(neighbor->structure);
    free(neighbors);
  }

  handled.clear();

  return MaxEnergy;
}
