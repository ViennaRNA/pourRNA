/*
 * StatePairCollector.cpp
 *
 *  Created on: 26.07.2014
 *      Author: Gregor Entzian
 */

#include "StatePairCollector.h"

StatePairCollector::StatePairCollector(size_t currentMinID,
		PairHashTable::HashTable& minima, SC_PartitionFunction::Z_Matrix& z,
		const size_t maxGradWalkHashed, Concurrent_Queue<MyState> *discoveredMinima,
		double boltzmannWeightTemperature, unsigned int move_set, PairHashMap::HashMap& all_saddles) :
		CurMinID(currentMinID), Minima(minima), Z(z), GradWalk(
		    move_set, maxGradWalkHashed), NumberOfOuterStates(0), DiscoveredMinima(
				discoveredMinima), BoltzmannWeightTemperature(
				boltzmannWeightTemperature), All_Saddles(all_saddles) {

  current_min = Minima.begin()->first.clone();
}

StatePairCollector::~StatePairCollector() {
  delete current_min;
}

void StatePairCollector::add(vrna_fold_compound_t *vc,
		const MyState* const state1, const MyState* const state2,
		bool firstIsSmaller) {
	// int minimumEnergy = move_gradient (RNAsequence, structurePairTable, s0, s1, 0, 0, 0);
	size_t neighborMinID = -1;
	MyState* newMin = GradWalk.walk(vc, *state2);

	// search for newMin
	PairHashTable::HashTable::const_iterator minEntry = Minima.find(*newMin);
	// check if newMin is known in Minima
	if (minEntry == Minima.end()) {
		//If the min does not exist, create a new index and add it.
		neighborMinID = Minima.size();
		Minima.insert( { MyState(*newMin), neighborMinID });
		if (DiscoveredMinima != NULL)
			DiscoveredMinima->push((*newMin));

	} else {
		// get index of newMin in Minima
		neighborMinID = minEntry->second;

		//add the first saddle for this new minimum
	  MyState* a = current_min->clone();
	  MyState* b = newMin->clone();
	  std::pair<MyState,MyState> state_pair = std::pair<MyState,MyState>({*a,*b});

	  std::unique_lock<std::mutex> mlock(mutex_);
	  auto pair_it =  All_Saddles.find(state_pair);
	  if(pair_it == All_Saddles.end()){
	    std::pair<MyState,MyState> state_pair2 = std::pair<MyState,MyState>({*b,*a});
	    auto pair_it2 =  All_Saddles.find(state_pair2);
	    if(pair_it2 == All_Saddles.end()){
	      // add higher state as saddle
	      if(firstIsSmaller)
	        All_Saddles[state_pair] = MyState(*state2);
	      else
	        All_Saddles[state_pair] = MyState(*state1);
	    }
	  }

	  mlock.unlock();
	  cond_.notify_one();

	  delete a;
	  delete b;
	}



	//now we should free the memory of newMin, because it already exists in the minima list.
	delete newMin;

	SC_PartitionFunction::PairID pairID = { CurMinID, neighborMinID };
	SC_PartitionFunction::Z_Matrix::iterator zIt = Z.find(pairID);
	if (zIt == Z.end()) {
		// vc->params->temperature is in Celsius.
		Z[pairID].initialize(vc->params->temperature);
	}

	if (firstIsSmaller) {
		// update Z matrix with basin state
		Z[pairID].add(*state2);
	} else {
		// update Z matrix with non-basin state
		//TODO: compute the contact surface only once!
		Z[pairID].add(*state1);
		NumberOfOuterStates++;
	}

	/*
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
	 */
}
