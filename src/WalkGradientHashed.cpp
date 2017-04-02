/*
 * WalkGradientHashed.cpp
 *
 *  Created on: 14.08.2014
 *      Author: Quin
 */

#include "WalkGradientHashed.h"

WalkGradientHashed::WalkGradientHashed (std::string sequence, const size_t maxHashSize) :
    State2min (maxHashSize)
{
}

WalkGradientHashed::~WalkGradientHashed ()
{
}

struct_en*
WalkGradientHashed::get_Neighbors_pt (vrna_fold_compound_t *vc, struct_en* structureEnergy)
{
	  //browse_neighs_pt_par_list_alloc_energy (Sequence, structureEnergy, S0, S1, 0, 0, 0, &NeighborList, &EnergyParameter, MaxNeighbors);
	  vrna_move_t *tmp_neighbors = vrna_neighbors (vc, structureEnergy->structure,  VRNA_MOVESET_DEFAULT);
	  size_t count = 0;
	  for(vrna_move_t *m = tmp_neighbors; m->pos_5 != 0; m++)
		  count++;

	  struct_en* neighbors = (struct_en*)malloc(sizeof(struct_en)*(count+1));
	  int i = 0;
	  for(vrna_move_t *m = tmp_neighbors; m->pos_5 != 0; m++, i++){
		  int energy = vrna_eval_move_pt(vc,structureEnergy->structure,m->pos_5,m->pos_3);
		  short *pt_neighbor = vrna_ptable_copy(structureEnergy->structure);
		  vrna_move_apply(pt_neighbor, m);
		  struct_en * neighbor = &neighbors[i];
		  neighbor->energy = energy;
		  neighbor->structure = pt_neighbor;
	  }
	  neighbors[i].structure = NULL;
	  free(tmp_neighbors);
	  return neighbors;
}

MyState*
WalkGradientHashed::walkGradient (vrna_fold_compound_t *vc, char * rnaSequence, const MyState& startState,
				  SpookyHashMap::HashTable& state2min_)
{
  MyState* minimalNeighbor = (MyState*) startState.clone ();
  struct_en tmpMinimalNeighbor;

  SpookyHashMap::HashTable::iterator hashed;

  bool foundBetterNeighbor = false;
  do
    {
      foundBetterNeighbor = false;

      // check if we already know the local minimum for this State
      //  uint64_t hash=SpookyHash::Hash64(minimalNeighbor->structure,minimalNeighbor->structure[0],0);
      hashed = state2min_.find (*minimalNeighbor);

      if (hashed != state2min_.end ())
	{
	  delete minimalNeighbor;
	  //free (minimalNeighbor->structure);
	  //minimalNeighbor->structure = NULL;
	  minimalNeighbor = hashed->second.clone ();
	  return minimalNeighbor;
	}

      tmpMinimalNeighbor.energy = minimalNeighbor->energy;
      tmpMinimalNeighbor.structure = minimalNeighbor->structure;


      struct_en *neighbors = get_Neighbors_pt (vc, &tmpMinimalNeighbor);
      //browse_neighs_pt_par_list_alloc_energy (rnaSequence, &tmpMinimalNeighbor, S0, S1, 0, 0, 0, &NeighborList, &EnergyParameter, MaxNeighbors);

      for (struct_en *neighbor = neighbors; neighbor->structure != NULL; neighbor++)
	{
	  if (neighbor->energy < minimalNeighbor->energy
	      || (neighbor->energy == minimalNeighbor->energy
		  && StructureUtils::IsSmaller (neighbor->structure, minimalNeighbor->structure)))
	    {
	      minimalNeighbor->energy = neighbor->energy;
	      copy_arr (minimalNeighbor->structure, neighbor->structure);
	      foundBetterNeighbor = true;
	    }

	  free(neighbor->structure);
	}

      free(neighbors);

    }
  while (foundBetterNeighbor);

  // store mapping of start to local minimum in hash
  state2min_.insert (std::pair<MyState, MyState> (startState, *minimalNeighbor));

  return minimalNeighbor;
}

MyState*
WalkGradientHashed::walk (vrna_fold_compound_t *vc,char * rnaSequence, const MyState& startState)
{
  MyState* result = walkGradient (vc, rnaSequence, startState, State2min);
  return result;
}

