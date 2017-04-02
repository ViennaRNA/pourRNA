/*
 * WalkGradientHashed.cpp
 *
 *  Created on: 14.08.2014
 *      Author: Quin
 */

#include "WalkGradientHashed.h"

WalkGradientHashed::WalkGradientHashed (std::string sequence, size_t maxNeighbors, const size_t maxHashSize) :
    State2min (maxHashSize), EnergyParameter (*GlobalParameter::getInstance ()->getEnergyParameter ())
{
  NeighborList.list = NULL;
  NeighborList.list_length = 0;

  MaxNeighbors = maxNeighbors;

  //init for browse_neighs function.

  make_pair_matrix (); //is important before encode_sequence.
  S0 = encode_sequence (sequence.c_str (), 0);
  S1 = encode_sequence (sequence.c_str (), 1);
}

WalkGradientHashed::~WalkGradientHashed ()
{
  freeNeighborList ();

  free (S0);
  S0 = NULL;
  free (S1);
  S1 = NULL;
}

MyState*
WalkGradientHashed::walkGradient (char * rnaSequence, const MyState& startState,
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
      browse_neighs_pt_par_list_alloc_energy (rnaSequence, &tmpMinimalNeighbor, S0, S1, 0, 0, 0, &NeighborList,
					      &EnergyParameter, MaxNeighbors);

      for (size_t i = 0; i < NeighborList.list_length; i++)
	{
	  struct_en* neighbor = &NeighborList.list[i];
	  if (neighbor->energy < minimalNeighbor->energy
	      || (neighbor->energy == minimalNeighbor->energy
		  && StructureUtils::IsSmaller (neighbor->structure, minimalNeighbor->structure)))
	    {
	      minimalNeighbor->energy = neighbor->energy;
	      copy_arr (minimalNeighbor->structure, neighbor->structure);
	      foundBetterNeighbor = true;
	    }
	}

      clearNeighborList ();

    }
  while (foundBetterNeighbor);

  // store mapping of start to local minimum in hash
  state2min_.insert (std::pair<MyState, MyState> (startState, *minimalNeighbor));

  return minimalNeighbor;
}

MyState*
WalkGradientHashed::walk (char * rnaSequence, const MyState& startState)
{
  MyState* result = walkGradient (rnaSequence, startState, State2min);
  return result;
}

