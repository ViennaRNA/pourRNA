/*
 * SC_DotPlot.cpp
 *
 *  Created on: 02.04.2015
 *      Author: Martin Mann
 */

#include "SC_DotPlot.h"

#include <sstream>

SC_DotPlot::SC_DotPlot (const double temperature, const bool storeEnergies) :
    SC_PartitionFunction (), bpWeightSum ()
{
  // initialize data structures
  initialize (temperature, storeEnergies);
}

void
SC_DotPlot::add (const MyState& state)
{
  // super class function : update partition function
  SC_PartitionFunction::add (state);

  // for all base pairs : update base pair weight sum
  const MyState::Structure& stateStructure = state.getStructure ();
  PairID bp;
  const double stateWeight = getBoltzmannWeight (state);
  for (int i = 1; i <= stateStructure[0]; i++)
    {
      // check if paired to higher index = base pair
      if (stateStructure[i] > i)
	{
	  // correct indices to start counting at 0
	  bp.first = i - 1;
	  bp.second = stateStructure[i] - 1;
	  // check if this base pair is not known
	  if (bpWeightSum.find (bp) == bpWeightSum.end ())
	    {
	      // initialize sum with current weight
	      bpWeightSum[bp] = stateWeight;
	    }
	  else
	    {
	      // update weight sum
	      bpWeightSum[bp] += stateWeight;
	    }
	}
    }
}

SC_DotPlot::DotPlot
SC_DotPlot::getBasePairProbabilities () const
{
  // use static function for computation
  return SC_DotPlot::getBasePairProbabilities (bpWeightSum, getZ ());
}

SC_DotPlot::~SC_DotPlot ()
{
}

void
SC_DotPlot::initialize (const double temperature, const bool storeEnergies)
{
  // super class function
  SC_PartitionFunction::initialize (temperature, storeEnergies);

  // init local data structures
  bpWeightSum.clear ();
}

SC_DotPlot::DotPlot
SC_DotPlot::getBasePairProbabilities (const DotPlot& bpWeightSums, const double Z)
{
  // copy weight sums
  DotPlot bpProbs (bpWeightSums);
  // scale weight sums with partition functions to probabilities
  for (DotPlot::iterator bp2weight = bpProbs.begin (); bp2weight != bpProbs.end (); bp2weight++)
    {
      bp2weight->second /= Z;
    }
  // return probabilities
  return bpProbs;
}

bool
SC_DotPlot::writeDotPlot_PS (const std::string & absoluteFileName, const std::string & sequence,
			     const DotPlot & dotPlot)
{
  // check if filename given
  if (absoluteFileName.empty())
    {
      return false;
    }
  // check if sequence given
  if (sequence.empty())
    {
      return false;
    }

  // allocate data structure for vienna package
  plist* bpProbVienna = (plist *) space ((dotPlot.size () + 1) * sizeof(plist));
  plist mfeViennaDummy;

  // fill dummy content = end of list marking entry
  mfeViennaDummy.i = 0;
  mfeViennaDummy.j = 0;
  mfeViennaDummy.p = 0;
  mfeViennaDummy.type = 0;

  // mark end of base pair prob list
  size_t i = dotPlot.size ();
  bpProbVienna[i] = mfeViennaDummy;

  // copy base pair probabilities
  bool allSane = true;
  DotPlot::const_iterator bp2weight;
  for (i = 0, bp2weight = dotPlot.begin (); allSane && bp2weight != dotPlot.end (); bp2weight++)
    {
      // check if entry is sane given the sequence
      allSane = bp2weight->first.first < sequence.size () && bp2weight->first.second < sequence.size ()
	  && bp2weight->second >= 0.0 && bp2weight->second <= 1.0;
      if (allSane)
	{
    	  // NOTE: need to shift the indices by 1 for Vienna package data
	  bpProbVienna[i].i = 1 + bp2weight->first.first;
	  bpProbVienna[i].j = 1 + bp2weight->first.second;
	  bpProbVienna[i].p = bp2weight->second;
	  bpProbVienna[i].type = 0;
	  i++;
	} else {
		return false;
	}
    }
  if (allSane)
    {
      // get char array handles for sequence and filename
      char * seqChar = (char*) sequence.c_str ();
      char * fileChar = (char*) absoluteFileName.c_str ();
      char * emptyChar = (char*) std::string().c_str();

      // write dot plot postscript file and check for success
      allSane = (1 == PS_dot_plot_list (seqChar, fileChar, bpProbVienna, &mfeViennaDummy, emptyChar));
    }

  // free memory
  free (bpProbVienna);

  // return whether or not all went fine
  return allSane;
}

SC_DotPlot::DotPlot
SC_DotPlot::readDotPlot_PS (std::istream & input)
{
  // container to fill
  DotPlot bpProb;

  // line of the file to be filled
  std::string line;

  size_t i, j;
  double prob;
  std::stringstream lineStream;

  // linewise read from the stream
  while (input.good ())
    {
      // read next line from stream
      std::getline (input, line);
      // check if base pair probability line
      if (line.find ("ubox") != std::string::npos)
	{
	  // initialize parsing
	  lineStream.str (line);
	  i = 0;
	  j = 0;
	  prob = -1.0;
	  // read data
	  lineStream >> i;
	  lineStream >> j;
	  lineStream >> prob;
	  // check parsing success
	  if (i != 0 && j != 0 && prob >= 0.0)
	    {
	      // store base pair probability data
	      // NOTE: store prob^2 when read from vienna dot plot
	      bpProb[SC_PartitionFunction::PairID(i, j)] = prob * prob;
	    }
	}
    }

  return bpProb;
}
