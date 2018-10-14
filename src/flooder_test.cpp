/*
 * parKin.cpp
 *
 * Main file for ParKin 1.0
 *
 *  Created on: 21.07.2014
 *      Author: Gregor Entzian
 * Purpose: Test the basin flooder.
 */

#include "Flooder.h"
#include <stdio.h>
#include <chrono>
#include <time.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////
int MaximalNeighbors = 0;
int countNeighbor(struct_en* a, struct_en* b) {
	MaximalNeighbors++;
	return 0;//return zero=do not break up.
}
/* getMaxNeighborNumber function takes a RNA Sequence as parameter
 *  and return the maximum number of Neighbors that current State has
 *  @param rnaSeq, represents RNA sequence.
 */
size_t getMaxNeighborNumber(std::string rnaSequence) {
// create open chain State
	std::string openChain(rnaSequence.size(), '.');
	browse_neighs((char*) rnaSequence.c_str(), (char*) openChain.c_str(), 0, 0,
			0, countNeighbor);
	return MaximalNeighbors;
} // end function

//////////////////////////////////////////////////////////////////////////


/**
 * program main entry
 *
 * @param argc number of program arguments
 * @param argv array of program arguments of length argc
 */
int main( int argc, char* argv[] ){
  try
    {
      std::cout << "This is ParKin, the basin flooder test.\n";

      size_t maxToStore=50000000;
      size_t maxToHash=50000000;

      std::string rnaSequence = "GGGAAUUAUUGUUCCCUGAGAGCGGUAGUUCUC";//d33   //boris"GACCGGAAGGUCCGCCUUCC";
      //"ACUGUAUGCGCGU";
      /* start structure is open chain */
      string openChain (rnaSequence.length (), '.');
      string rnaStructure = openChain; //"....(((.(.....).).))";
     rnaStructure = "((((((((((((((.....))))))))))))))"; //mfe

      float energy = round (energy_of_structure (rnaSequence.c_str (), rnaStructure.c_str (), 0) * 100);
      const size_t maxNeighbors = getMaxNeighborNumber(rnaSequence);
      //flooder test:
      Flooder * myFlooder = new Flooder (rnaSequence, maxNeighbors, 9999, maxToStore);
      short * rnaStructurePT = make_pair_table (rnaStructure.c_str ());
      //openChain="........((.......)).";
      MyState startState (energy, rnaStructurePT);
      // start walking to find gradient basin minimum for start state
      MyState* currentMin = WalkGradientHashed (rnaSequence, maxToHash).walk ((char*) rnaSequence.c_str (),
									   startState);
      size_t currentMinID = 0;
      PairHashTable::HashTable Minima;
      Minima.insert (
	{ *currentMin, currentMinID });
      vector<size_t> indicesOfNeighbors;
      SC_PartitionFunction::Z_Matrix z;
      SC_PartitionFunction scBasin = z[SC_PartitionFunction::PairID (currentMinID, currentMinID)];
      StatePairCollector* spc = new StatePairCollector (rnaSequence, maxNeighbors, currentMinID, Minima, z,
                                                        maxToHash);


      //init stopwatch:
      std::chrono::time_point<std::chrono::system_clock> start, end;
      start = std::chrono::system_clock::now ();

      myFlooder->floodBasin (startState, scBasin, *spc);

      end = std::chrono::system_clock::now ();
      std::chrono::duration<double> elapsed_seconds = end - start;
      std::time_t end_time = std::chrono::system_clock::to_time_t (end);


      size_t Nneighbors = indicesOfNeighbors.size ();
      size_t Nmin = Minima.size ();
      std::cout << "BasinSize: " << scBasin.size () << "\n";
      std::cout << "BasinPartitionFunction: " << scBasin.getZ () << "\n";
      std::cout << "Neighbors " << Nneighbors << "\n";
      std::cout << "Minima " << Nmin << "\n";
      PairHashTable::HashTable::iterator it = Minima.begin ();
      while (it != Minima.end ())
	{
	  struct_en str_en = struct_en({it->first.energy,it->first.structure});
	  print_stren (stdout, &str_en);
	  ++it;
	}

      //clean up
      std::free (rnaStructurePT);

      //print stopwatch:
        std::cout << "finished computation at " << std::ctime (&end_time) << "elapsed time: " << elapsed_seconds.count ()
            << "s\n";

    }
  catch (std::exception & e)
    {
      std::cerr << "\n Exception raised : " << e.what () << "\n";
      return (-1);
    }

  return (0);
}

