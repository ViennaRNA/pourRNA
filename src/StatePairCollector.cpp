/*
 * StatePairCollector.cpp
 *
 *  Created on: 26.07.2014
 *      Author: Gregor Entzian
 */

#include "StatePairCollector.h"
#include "SC_PartitionFunction.h"

StatePairCollector::StatePairCollector(size_t                           currentMinID,
                                       PairHashTable::HashTable&        minima,
                                       SC_PartitionFunction::Z_Matrix&  z,
                                       const size_t                     maxGradWalkHashed,
                                       Concurrent_Queue<MyState>        *discoveredMinima,
                                       double                           boltzmannWeightTemperature,
                                       double                           gas_constant,
                                       double                           mfe,
                                       unsigned int                     move_set,
                                       MapOfMaps*                       all_saddles,
                                       const char                       *sourceStructure,
                                       const char                       *targetStructure,
                                       int                              maxBPdist) :
  CurMinID(currentMinID), Minima(minima), Z(z), GradWalk(
    move_set, maxGradWalkHashed), NumberOfOuterStates(0), DiscoveredMinima(
    discoveredMinima), BoltzmannWeightTemperature(
    boltzmannWeightTemperature), GasConstant(gas_constant), MFE(mfe), All_Saddles(all_saddles),
  SourceStructure(sourceStructure), TargetStructure(targetStructure), MaxBPdist(maxBPdist),
  MininaToIgnore()
{
  current_min = Minima.begin()->first.clone();
  CurrentStructure = vrna_db_from_ptable(current_min->getStructure());
}


StatePairCollector::~StatePairCollector()
{
  delete current_min;
  free(CurrentStructure);
  MininaToIgnore.clear();
}


void
StatePairCollector::add(vrna_fold_compound_t  *vc,
                        const MyState *const  state1,
                        const MyState *const  state2,
                        bool                  firstIsSmaller)
{
  // int minimumEnergy = move_gradient (RNAsequence, structurePairTable, s0, s1, 0, 0, 0);
  size_t                                    neighborMinID = -1;
  MyState                                   *newMin       = GradWalk.walk(vc, *state2);

  // search for newMin
  PairHashTable::HashTable::const_iterator  minEntry = Minima.find(*newMin);

  // search in black list
  HashSet::UnorderedHashSet::const_iterator badMin = MininaToIgnore.find(*newMin);

  if (badMin != MininaToIgnore.end()) {
    // bad minimum: dont add to contact surface, dont add to all minima, dont report as discovered minimum.
    delete newMin;
    return;
  }

  // check if newMin is known in Minima
  if (minEntry == Minima.end()) {
    //filter discovered minima at first
    // if the min is in the already filtered list, abort this function!!!

    //max base pair distance filter:
    if (MaxBPdist >= 0 && SourceStructure != NULL && TargetStructure != NULL) {
      int dist_from_source  = vrna_bp_distance((const char *)CurrentStructure, SourceStructure);
      int dist_from_target  = vrna_bp_distance((const char *)CurrentStructure, TargetStructure);
      if (dist_from_source + dist_from_target > MaxBPdist) {
        MininaToIgnore.insert(MyState(*newMin));
        delete newMin;
        return;
      }
    }

    //If the min does not exist, create a new index and add it.
    neighborMinID = Minima.size();
    Minima.insert({ MyState(*newMin), neighborMinID });

    //report discovered min.
    if (DiscoveredMinima != NULL)
      DiscoveredMinima->push((*newMin));

    //add the first saddle for this new minimum
    if(All_Saddles != NULL){
      MyState                       *a          = current_min->clone();
      MyState                       *b          = newMin->clone();
      std::pair<MyState, MyState>   state_pair  = std::pair<MyState, MyState>({ *a, *b });

      std::unique_lock<std::mutex>  mlock(mutex_);

      auto saddles_from_a = (*All_Saddles)[*a];
      auto saddle_a_b = saddles_from_a.find(*b);
      if (saddle_a_b == saddles_from_a.end()){
        if (firstIsSmaller){
           (*All_Saddles)[*a][*b] = MyState(*state2);
        }
        else{
           (*All_Saddles)[*a][*b] = MyState(*state1);
        }
      }
      mlock.unlock();
      cond_.notify_one();
      delete a;
      delete b;
    }
  } else {
    // get index of newMin in Minima
    neighborMinID = minEntry->second;
  }

  //now we should free the memory of newMin, because it already exists in the minima list.
  delete newMin;

  SC_PartitionFunction::PairID              pairID = {
    CurMinID, neighborMinID
  };
  SC_PartitionFunction::Z_Matrix::iterator  zIt = Z.find(pairID);
  if (zIt == Z.end())
    // vc->params->temperature is in Celsius.
    Z[pairID].initialize(BoltzmannWeightTemperature, GasConstant, MFE);

  if (firstIsSmaller) {
    // update Z matrix with basin state
    Z[pairID].add(*state2);
  } else {
    // update Z matrix with non-basin state
    //TODO: compute the contact surface only once!
    Z[pairID].add(*state1);
    NumberOfOuterStates++;
  }
}
