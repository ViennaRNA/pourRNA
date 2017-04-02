/*
 * StatePairCollector.cpp
 *
 *  Created on: 26.07.2014
 *      Author: Gregor Entzian
 */

#include "StatePairCollector.h"

StatePairCollector::StatePairCollector (std::string rnaSequence, size_t maxNeighbors, size_t currentMinID,
					PairHashTable::HashTable& minima, SC_PartitionFunction::Z_Matrix& z,
					const size_t maxGradWalkHashed) :
    RNAsequence ((char*) rnaSequence.c_str ()), CurMinID (currentMinID), Minima (minima), Z (z), GradWalk (
	rnaSequence, maxNeighbors, maxGradWalkHashed), NumberOfOuterStates (0)
{

}

StatePairCollector::~StatePairCollector ()
{
}

void
StatePairCollector::add (const MyState* const state1, const MyState* const state2)
{
  // int minimumEnergy = move_gradient (RNAsequence, structurePairTable, s0, s1, 0, 0, 0);
  size_t neighborMinID = -1;
  MyState* newMin = GradWalk.walk ((char *)RNAsequence.c_str(), *state2);
  // search for newMin
  PairHashTable::HashTable::const_iterator minEntry = Minima.find (*newMin);
  // check if newMin is known in Minima
  if (minEntry == Minima.end ())
    {
      //If the min does not exist, create a new index and add it.
      neighborMinID = Minima.size ();
      Minima.insert (
	{ MyState (*newMin), neighborMinID });

    }
  else
    {
      // get index of newMin in Minima
      neighborMinID = minEntry->second;
    }
  //now we should free the memory of newMin, because it already exists in the minima list.
  delete newMin;

  SC_PartitionFunction::PairID pairID (CurMinID, neighborMinID);
  SC_PartitionFunction::Z_Matrix::iterator zIt = Z.find(pairID);
  if(zIt == Z.end())
    {
      double temp = GlobalParameter::getInstance ()->getBoltzmannWeightTemperature();
      Z[pairID].initialize(temp);
    }
  //  identify the higher energy state of the current state pair == saddle point
  if ((state1)->energy < (state2)->energy
      || ((state1)->energy == (state2)->energy && StructureUtils::IsSmaller (state1->structure, state2->structure))) //compare structure to break ties.
    {
      // update Z matrix with basin state
      Z[pairID].add (*state2);
    }
  else
    {
      // update Z matrix with non-basin state
      Z[pairID].add (*state1);
      NumberOfOuterStates++;
    }
}
