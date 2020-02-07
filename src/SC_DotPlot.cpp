/*
 * SC_DotPlot.cpp
 *
 *  Created on: 02.04.2015
 *      Author: Martin Mann
 */

#include "SC_DotPlot.h"

#include <sstream>
#include <fstream>
extern "C" {
#include <ViennaRNA/utils.h>
#include <ViennaRNA/plotting/probabilities.h>

#include <ViennaRNA/utils/structures.h> //plist types
#include <ViennaRNA/MEA.h>
}

SC_DotPlot::SC_DotPlot (const double  temperature,
                        const double  gas_constant,
                        double mfe,
                        const bool    storeEnergies) :
  SC_PartitionFunction(temperature, gas_constant, mfe, storeEnergies), bpWeightSum()
{
  // initialize data structures
  // initialize(temperature, storeEnergies);
}


void
SC_DotPlot::add(const MyState& state)
{
  // super class function : update partition function
  SC_PartitionFunction::add(state);

  // for all base pairs : update base pair weight sum
  const MyState::Structure& stateStructure = state.getStructure();
  PairID                    bp;
  const double              stateWeight = getBoltzmannWeight(state);
  for (int i = 1; i <= stateStructure[0]; i++) {
    // check if paired to higher index = base pair
    if (stateStructure[i] > i) {
      // correct indices to start counting at 0
      bp.first  = i - 1;
      bp.second = stateStructure[i] - 1;
      // check if this base pair is not known
      if (bpWeightSum.find(bp) == bpWeightSum.end())
        // initialize sum with current weight
        bpWeightSum[bp] = stateWeight;
      else
        // update weight sum
        bpWeightSum[bp] += stateWeight;
    }
  }
}


SC_DotPlot::DotPlot
SC_DotPlot::getBasePairProbabilities() const
{
  // use static function for computation
  return SC_DotPlot::getBasePairProbabilities(bpWeightSum, getZ());
}


SC_DotPlot::~SC_DotPlot ()
{
}


void
SC_DotPlot::initialize(const double temperature,
                       const double gas_constant,
                       double mfe,
                       const bool   storeEnergies)
{
  // super class function
  SC_PartitionFunction::initialize(temperature, gas_constant, mfe, storeEnergies);

  // init local data structures
  bpWeightSum.clear();
}


SC_DotPlot::DotPlot
SC_DotPlot::getBasePairProbabilities(const DotPlot& bpWeightSums,
                                     const double   Z)
{
  // copy weight sums
  DotPlot bpProbs(bpWeightSums);

  // scale weight sums with partition functions to probabilities
  for (DotPlot::iterator bp2weight = bpProbs.begin(); bp2weight != bpProbs.end(); bp2weight++)
    bp2weight->second /= Z;
  // return probabilities
  return bpProbs;
}


bool
SC_DotPlot::writeDotPlot_PS(const std::string & absoluteFileName,
                            const std::string & sequence,
                            const DotPlot &     dotPlot)
{
  // check if filename given
  if (absoluteFileName.empty())
    return false;

  // check if sequence given
  if (sequence.empty())
    return false;

  // allocate data structure for vienna package
  plist *bpProbVienna = (plist *)space((dotPlot.size() + 1) * sizeof(plist));
  plist mfeViennaDummy;

  // fill dummy content = end of list marking entry
  mfeViennaDummy.i    = 0;
  mfeViennaDummy.j    = 0;
  mfeViennaDummy.p    = 0;
  mfeViennaDummy.type = 0;

  // mark end of base pair prob list
  size_t                  i = dotPlot.size();
  bpProbVienna[i] = mfeViennaDummy;

  // copy base pair probabilities
  bool                    allSane = true;
  DotPlot::const_iterator bp2weight;
  for (i = 0, bp2weight = dotPlot.begin(); allSane && bp2weight != dotPlot.end(); bp2weight++) {
    // check if entry is sane given the sequence
    allSane = bp2weight->first.first < sequence.size() && bp2weight->first.second < sequence.size()
              && bp2weight->second >= 0.0 && bp2weight->second <= 1.0;
    if (allSane) {
      // NOTE: need to shift the indices by 1 for Vienna package data
      bpProbVienna[i].i     = 1 + bp2weight->first.first;
      bpProbVienna[i].j     = 1 + bp2weight->first.second;
      bpProbVienna[i].p     = bp2weight->second;
      bpProbVienna[i].type  = 0;
      i++;
    } else {
      return false;
    }
  }
  if (allSane) {
    // get char array handles for sequence and filename
    char  *seqChar    = (char *)sequence.c_str();
    char  *fileChar   = (char *)absoluteFileName.c_str();
    char  *emptyChar  = NULL; //(char *)std::string().c_str();

    // write dot plot postscript file and check for success
    allSane = (1 == PS_dot_plot_list(seqChar, fileChar, bpProbVienna, &mfeViennaDummy, emptyChar));
  }

  // free memory
  free(bpProbVienna);

  // return whether or not all went fine
  return allSane;
}

bool
SC_DotPlot::writeDotPlot_PS_with_mfe_and_mea(const std::string & absoluteFileName,
                            const std::string & sequence,
                            const DotPlot &     dotPlot,
                            const std::string & mfe_structure,
                            const std::string & mea_structure)
{
  // check if filename given
  if (absoluteFileName.empty())
    return false;

  // check if sequence given
  if (sequence.empty())
    return false;

  // allocate data structure for vienna package
  plist *bpProbVienna = (plist *)space((dotPlot.size() + 1) * sizeof(plist));
  plist mfeViennaDummy;

  // fill dummy content = end of list marking entry
  mfeViennaDummy.i    = 0;
  mfeViennaDummy.j    = 0;
  mfeViennaDummy.p    = 0;
  mfeViennaDummy.type = 0;

  // mark end of base pair prob list
  size_t                  i = dotPlot.size();
  bpProbVienna[i] = mfeViennaDummy;

  // copy base pair probabilities
  bool                    allSane = true;
  DotPlot::const_iterator bp2weight;
  for (i = 0, bp2weight = dotPlot.begin(); allSane && bp2weight != dotPlot.end(); bp2weight++) {
    // check if entry is sane given the sequence
    allSane = bp2weight->first.first < sequence.size() && bp2weight->first.second < sequence.size()
              && bp2weight->second >= 0.0 && bp2weight->second <= 1.0;
    if (allSane) {
      // NOTE: need to shift the indices by 1 for Vienna package data
      bpProbVienna[i].i     = 1 + bp2weight->first.first;
      bpProbVienna[i].j     = 1 + bp2weight->first.second;
      bpProbVienna[i].p     = bp2weight->second;
      bpProbVienna[i].type  = 0;
      i++;
    } else {
      return false;
    }
  }
  if (allSane) {
    // get char array handles for sequence and filename
    char  *seqChar    = (char *)sequence.c_str();
    char  *fileChar   = (char *)absoluteFileName.c_str();
    char  *emptyChar  = NULL; //(char *)std::string().c_str();

    // write dot plot postscript file and check for success
    allSane = (1 == PS_dot_plot_list(seqChar, fileChar, bpProbVienna, &mfeViennaDummy, emptyChar));

    if (allSane){
      std::string line;
      std::ifstream infile_s(fileChar);
      std::stringstream infile;
      infile << infile_s.rdbuf();
      infile_s.close();

      std::fstream outfile (fileChar, std::fstream::binary | std::fstream::in | std::fstream::out);
      long position = 0;
      while (std::getline(infile, line))
      {
        // find box definition and skip 4 lines
        if (line.rfind("/box {", 0) == 0) {
          for(int i = 0; i < 4; i++){
            std::getline(infile, line);
          }
          position = infile.tellg();
          outfile.seekp(position, std::fstream::beg);

          std::string rect_definitions;
          rect_definitions.append("          /rect_flat { %size x y box - draws box centered on x,y\n");
          rect_definitions.append("             2 index 0.5 mul sub            % x -= 0.5\n");
          rect_definitions.append("             exch 2 index 0.5 mul sub exch  % y -= 0.5\n");
          rect_definitions.append("             3 -1 roll dup 0.5 mul rectfill\n");
          rect_definitions.append("          } bind def\n");

          rect_definitions.append("          /lboxl {\n");
          rect_definitions.append("             0 0 1 setrgbcolor\n");
          rect_definitions.append("             3 1 roll\n");
          rect_definitions.append("             len exch sub 1 add rect_flat\n");
          rect_definitions.append("             0 0 0 setrgbcolor\n");
          rect_definitions.append( "          } bind def\n");

          rect_definitions.append("          /lboxu {\n");
          rect_definitions.append("             1 0 0 setrgbcolor\n");
          rect_definitions.append("             3 1 roll\n");
          rect_definitions.append("             len exch sub 1.5 add rect_flat\n");
          rect_definitions.append("             0 0 0 setrgbcolor\n");
          rect_definitions.append("          } bind def\n");

          outfile << rect_definitions;

          short *mfe_pt = vrna_ptable(mfe_structure.c_str());
          short *mea_pt = vrna_ptable(mea_structure.c_str());
          std::string one_sf = "1.0000000";
          std::stringstream out_lines;

          while (std::getline(infile, line)){
            if (line.rfind("showpage",0) == 0) {
              position = infile.tellg();
              outfile.seekp(position + rect_definitions.length() -9, std::fstream::beg);
              for (int i = 1; i <= mfe_pt[0]; i++) {
                if (mfe_pt[i] > i) {
                  out_lines  << i << " " << mfe_pt[i] << " " << one_sf << " " << "lboxu\n";
                }
              }

              for (int i = 1; i <= mea_pt[0]; i++) {
                if (mea_pt[i] > i) {
                  out_lines << i << " " << mea_pt[i] << " " << one_sf << " " << "lboxl\n";
                }
              }
              std::string oline = out_lines.str();
              outfile << oline;
            }
            outfile << line + "\n"; // append overwritten lines.
          }
          free(mfe_pt);
          free(mea_pt);
          break;
        }
      }
      outfile.close();
    }
  }
  // free memory
  free(bpProbVienna);

  // return whether or not all went fine
  return allSane;
}

char * SC_DotPlot::mea_from_dotplot(const std::string & sequence, const DotPlot &     dotPlot, vrna_exp_param_t  *params){
  double weight;
  size_t nbp, n;
  nbp = dotPlot.size();
  vrna_ep_t *plist = (vrna_ep_t*)calloc(nbp+1, sizeof(vrna_ep_t));
  n = 0;
  for (auto bp2weight = dotPlot.begin(); bp2weight != dotPlot.end(); bp2weight++){
      plist[n].i = bp2weight->first.first +1;
      plist[n].j = bp2weight->first.second +1;
      plist[n].p = bp2weight->second;
      plist[n].type = VRNA_PLIST_TYPE_BASEPAIR;
      n++;
  }
  // list terminator
  plist[n].i = 0;
  plist[n].j = 0;
  plist[n].p = 0;
  char *mea_structure = (char *)calloc(sequence.length()+1, sizeof(char));
  //vrna_exp_param_t  *params = vrna_exp_params(NULL);
  memset(mea_structure, '.', sequence.length() * sizeof(char));
  mea_structure[sequence.length()] = '\0';

  float mea = MEA_seq(plist, sequence.c_str(), mea_structure, 2.0, params);
  free(plist);
  return mea_structure;
}


SC_DotPlot::DotPlot
SC_DotPlot::readDotPlot_PS(std::istream & input)
{
  // container to fill
  DotPlot           bpProb;

  // line of the file to be filled
  std::string       line;

  size_t            i, j;
  double            prob;
  std::stringstream lineStream;

  // linewise read from the stream
  while (input.good()) {
    // read next line from stream
    std::getline(input, line);
    // check if base pair probability line
    if (line.find("ubox") != std::string::npos) {
      // initialize parsing
      lineStream.str(line);
      i     = 0;
      j     = 0;
      prob  = -1.0;
      // read data
      lineStream >> i;
      lineStream >> j;
      lineStream >> prob;
      // check parsing success
      if (i != 0 && j != 0 && prob >= 0.0) {
        // store base pair probability data
        // NOTE: store prob^2 when read from vienna dot plot
        bpProb[SC_PartitionFunction::PairID(i, j)] = prob * prob;
      }
    }
  }

  return bpProb;
}
