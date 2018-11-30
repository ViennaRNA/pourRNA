/*
 * WalkGradientHashed.cpp
 *
 *  Created on: 14.08.2014
 *      Author: Quin
 */

#include "WalkGradientHashed.h"

WalkGradientHashed::WalkGradientHashed(unsigned int move_set,
                                       const size_t maxHashSize) :
  Move_set(move_set), State2min(maxHashSize)
{
}


WalkGradientHashed::~WalkGradientHashed()
{
}


struct_en *
WalkGradientHashed::convert_moves_to_neighbors(vrna_fold_compound_t *vc,
                                               struct_en            *structureEnergy,
                                               vrna_move_t          *moves,
                                               int                  moves_count)
{
  struct_en *neighbors = (struct_en *)malloc(
    sizeof(struct_en) * (moves_count + 1));
  int       i = 0;
  int       energy;
  short     *pt_neighbor;
  struct_en *neighbor;

  for (vrna_move_t *m = moves; m->pos_5 != 0; m++, i++) {
    energy      = vrna_eval_move_shift_pt(vc, m, structureEnergy->structure); //vrna_eval_move_pt(vc, structureEnergy->structure, m->pos_5, m->pos_3);
    pt_neighbor = vrna_ptable_copy(structureEnergy->structure);
    vrna_move_apply(pt_neighbor, m);
    neighbor            = &neighbors[i];
    neighbor->energy    = structureEnergy->energy + energy;
    neighbor->structure = pt_neighbor;
  }
  neighbors[i].structure = NULL;
  return neighbors;
}


MyState *
WalkGradientHashed::walkGradient(vrna_fold_compound_t       *vc,
                                 const MyState&             startState,
                                 SpookyHashMap::HashTable&  state2min_)
{
  MyState                             *minimalNeighbor = new MyState(startState);
  struct_en                           currentStructure;

  currentStructure.structure = vrna_ptable_copy(minimalNeighbor->structure);

  SpookyHashMap::HashTable::iterator  hashed;

  bool                                foundBetterNeighbor = false;
  size_t                              minimal_move_index  = -1;
  vrna_move_t                         *tmp_neighbors      = NULL;
  int                                 size_neighbors      = 0;
  short                               *prev_pt            = vrna_ptable_copy(
    minimalNeighbor->structure);
  struct_en                           currentNeighbor;
  currentNeighbor.structure = prev_pt;
  do {
    foundBetterNeighbor = false;

    // check if we already know the local minimum for this State
    //  uint64_t hash=SpookyHash::Hash64(minimalNeighbor->structure,minimalNeighbor->structure[0],0);
    hashed = state2min_.find(*minimalNeighbor);

    if (hashed != state2min_.end()) {
      delete minimalNeighbor;
      free(prev_pt);
      free(tmp_neighbors);
      free(currentStructure.structure);
      minimalNeighbor = hashed->second.clone();
      return minimalNeighbor;
    }

    currentStructure.energy = minimalNeighbor->energy;
    copy_arr(currentStructure.structure, minimalNeighbor->structure);

    free(tmp_neighbors);
    tmp_neighbors = vrna_neighbors(vc, currentStructure.structure,
                                   this->Move_set);
    size_neighbors = 0;
    for (vrna_move_t *m = tmp_neighbors; m->pos_5 != 0; m++)
      size_neighbors++;

    size_t move_index = 0;
    int move_energy;
    for (vrna_move_t *m = tmp_neighbors; m->pos_5 != 0; m++, move_index++) {
      move_energy             = vrna_eval_move_shift_pt(vc, m, currentStructure.structure); //vrna_eval_move_pt(vc, currentStructure.structure, m->pos_5, m->pos_3);
      currentNeighbor.energy  = currentStructure.energy + move_energy;
      if (currentNeighbor.energy <= minimalNeighbor->energy) {
        //...compare also the structure
        copy_arr(currentNeighbor.structure, currentStructure.structure);
        vrna_move_apply(currentNeighbor.structure, m);
        if (currentNeighbor.energy < minimalNeighbor->energy
            || (currentNeighbor.energy == minimalNeighbor->energy
                && StructureUtils::IsSmaller(
                  currentNeighbor.structure,
                  minimalNeighbor->structure))) {
          minimalNeighbor->energy = currentNeighbor.energy;
          copy_arr(minimalNeighbor->structure, currentNeighbor.structure);
          foundBetterNeighbor = true;
          minimal_move_index  = move_index;
        }
      }
    }
    copy_arr(currentNeighbor.structure, currentStructure.structure);
  } while (foundBetterNeighbor);

  // store mapping of start to local minimum in hash
  state2min_.insert(
    std::pair<MyState, MyState>(startState, *minimalNeighbor));
  free(prev_pt);
  free(tmp_neighbors);
  free(currentStructure.structure);
  return minimalNeighbor;
}


MyState *
WalkGradientHashed::walk(vrna_fold_compound_t *vc,
                         const MyState&       startState)
{
  MyState *result = walkGradient(vc, startState, State2min);

  return result;
}
