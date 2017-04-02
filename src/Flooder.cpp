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
Flooder::Flooder (std::string sequence, double maxEnergy, size_t maxToQueue) :
    Sequence ((char*) sequence.c_str ()), MaxEnergy (maxEnergy * 100.0), MaxStatesToQueue (
	maxToQueue), ProcessedStates (0)
{

}

Flooder::~Flooder ()
{
}

struct_en*
Flooder::get_Neighbors_pt (vrna_fold_compound_t *vc, struct_en* structureEnergy)
{
  //browse_neighs_pt_par_list_alloc_energy (Sequence, structureEnergy, S0, S1, 0, 0, 0, &NeighborList, &EnergyParameter, MaxNeighbors);
  vrna_move_t *tmp_neighbors = vrna_neighbors (vc, structureEnergy->structure,  VRNA_MOVESET_DEFAULT);
  size_t count = 0;
  for(vrna_move_t *m = tmp_neighbors; m->pos_5 != 0; m++)
	  count++;

  struct_en* neighbors = (struct_en*)malloc(sizeof(struct_en)*(count+1));
  int i = 0;
  for(vrna_move_t *m = tmp_neighbors; m->pos_5 != 0; m++, i++){
	  double energy = vrna_eval_move_pt(vc,structureEnergy->structure,m->pos_5,m->pos_3)/100.0f;
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

int
Flooder::floodBasin (vrna_fold_compound_t *vc, const MyState& localMinState, StateCollector& scBasin, StatePairCollector& scSurface)
{
  if (localMinState.structure == NULL)
    std::cout << "No local Min given! (structure==NULL)" << std::endl;
  if (&scBasin == NULL)
    std::cout << "No StateCollector given! (scBasin==NULL)" << std::endl;
  if (localMinState.energy >= MaxEnergy)
    std::cout << "E(localMin) >= maxE" << std::endl;
  ProcessedStates = 0;

  //construct first element of pq.
  MyState localMin = MyState (localMinState);

  // the priority queue to store neighbors to handle in
  PriorityQueue<QueueValue> pq;
  // the hash to store seen elements of the basin
  PairHashTable::HashTable handled;
  // init using the given local minimum
  PriorityQueue<QueueValue>::InsertResult it = pq.insert (localMin);
  it.first->second.QueueState.structure = NULL; // has no minimal neighbor (is min)
  it.first->second.QueueState.energy = it.first->first.QueueState.energy;
  it.first->second.StateID = 0;

  typedef std::list<PriorityQueue<QueueValue>::QueueKey> NeighListPQ;

  // create temporary objects that are altered during the flooding
  MyState* topState = NULL;
  MyState* curNeigh = NULL;

  while (!pq.empty ())
    {
      // get top element
      PriorityQueue<QueueValue>::const_iterator top = pq.top ();
      // reference to top State
      topState = (MyState*) &top->first.QueueState;

      // get neighbor list
      struct_en tmpTopState (
	{ topState->energy, vrna_ptable_copy(topState->structure) });
      struct_en *neighbors = get_Neighbors_pt (vc, &tmpTopState);

      // contains all neighbors with E(top) < E <= maxE
      NeighListPQ toStore;
      // contains all neighbors with E < E(top)
      NeighListPQ toCheck;

      short * minNeigh = tmpTopState.structure; //init neighbor with minimal energy
      int minNeighE = tmpTopState.energy; //init energy of the minimal neighbor
      // iterate through all neighbors
      for (struct_en *neighbor = neighbors; neighbor->structure != NULL; neighbor++)
	{
	  // get neighbor data
	  int neighborE = neighbor->energy;
	  // the compressed sequence of the current neighbor checked
	  short * neighborStructure = neighbor->structure;

	  // check for lowest neighbor
	  if ((neighborE < minNeighE)
	      || (neighborE == minNeighE && StructureUtils::IsSmaller (neighborStructure, minNeigh)))
	    {
	      minNeigh = neighborStructure;
	      minNeighE = neighborE;
	    }
	  // check if neighbor is of interest and not rejected
	  if (neighborE <= MaxEnergy) // energy does not exceed maxE
	    {
	      // check if neighbor is of higher order --> to PQ
	      // energy higher than top or
	      if (neighborE > top->first.QueueState.energy
		  || (neighborE == top->first.QueueState.energy
		      && StructureUtils::IsGreater (neighborStructure, top->first.QueueState.structure)))
		{ // greater in sequence order

		  // store neighbor with higher energy that has to be added to pq
		  toStore.push_front (MyState (neighborE, neighborStructure));

		}
	      else if (&scSurface != NULL)
		{
		  // store all neighbors with lower energy
		  toCheck.push_back (MyState (neighborE, neighborStructure));
		}
	    }
	}

      // handle top depending on its basin membership
      if ((minNeighE == top->second.QueueState.energy
	  && StructureUtils::IsEqual (minNeigh, top->second.QueueState.structure)) // belongs to the basin
      || top->second.QueueState.structure == NULL) // top is the local minimum
	{
	  // store because top belongs to the basin
	  scBasin.add (*topState);
	  // add to hashed elements and mark that part of the basin
	  handled.insert (std::pair<MyState, size_t> (top->first.QueueState, 0));

	  // put all neighbors with higher energy to PQ
	  PriorityQueue<QueueValue>::iterator toFill;
	  // sort NeighborList with increasing energy
	  toStore.sort ();

	  for (NeighListPQ::const_iterator n = toStore.begin ();
	      n != toStore.end () && n->QueueState.energy <= MaxEnergy; n++)
	    {
	      toFill = pq.find (*n);
	      if (toFill == pq.end ())
		{
		  // insert empty value object
		  toFill = pq.insert (n->QueueState).first;
		  // feed data
		  toFill->second.QueueState.energy = top->first.QueueState.energy;
		  toFill->second.QueueState.structure = vrna_ptable_copy (top->first.QueueState.structure);
		  toFill->second.StateID = top->second.StateID;
		  // check queue size and reduce if necessary
		  if (pq.size () >= MaxStatesToQueue)
		    { // queue to big

		      // delete last element (i.e. the element with the highest energy).
		      PriorityQueue<QueueValue>::iterator toErase = --pq.end ();
		      pq.erase (toErase->first);

		      // adjust flood level to current maximal E in queue.
		      MaxEnergy = pq.getMaxE ();
		    }
		}
	      else  // update of entry if this one is smaller
	      if (toFill->second.QueueState.energy > top->first.QueueState.energy
		  || (toFill->second.QueueState.energy == top->first.QueueState.energy
		      && StructureUtils::IsGreater (toFill->second.QueueState.structure,
						    top->first.QueueState.structure)))
		{
		  // feed data
		  toFill->second.QueueState.energy = top->first.QueueState.energy;
		  toFill->second.QueueState.structure = vrna_ptable_copy (top->first.QueueState.structure);
		  toFill->second.StateID = top->second.StateID;
		}
	    }

	  // check if surface / contact plane is of interest
	  if (&scSurface != NULL)
	    {
	      // check all neighbors with lower energy if NOT part of the
	      // basin to enumerate the surface of the basin
	      for (NeighListPQ::const_iterator n = toCheck.begin (); n != toCheck.end (); n++)
		{
		  // check if neighbor is NOT part of basin
		  if (handled.find (n->QueueState) == handled.end ())
		    {
		      // get neighbor
		      curNeigh = (MyState*) &n->QueueState;
		      // add to surface reporter
		      scSurface.add (vc, topState, curNeigh);
		    }
		}
	    }
	}
      else
	{ // belongs NOT to the basin but in surface
	  // check if surface / contact plane is of interest
	  if (&scSurface != NULL)
	    {
	      // check all neighbors with lower energy IF part of the
	      // basin to enumerate the surface of the basin
	      for (NeighListPQ::const_iterator n = toCheck.begin (); n != toCheck.end (); n++)
		{
		  // check if neighbor IS part of basin
		  if (handled.find (n->QueueState) != handled.end ())
		    {
		      // get neighbor
		      curNeigh = (MyState*) &n->QueueState;
		      // add to surface reporter
		      scSurface.add (vc, curNeigh, topState);
		    }
		}
	    }
	}

      // remove top element
      pq.pop ();
      free(tmpTopState.structure);
      tmpTopState.structure = NULL;
      ProcessedStates++;
      toStore.clear ();
      toCheck.clear ();

      for (struct_en *neighbor = neighbors; neighbor->structure != NULL; neighbor++)
    	  free(neighbor->structure);
      free(neighbors);
    }

  handled.clear ();

  return MaxEnergy;
}

