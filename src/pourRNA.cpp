/*
 * pourRNA.cpp
 *
 * Main file for pourRNA 1.0
 *
 *  Created on: 09.08.2014
 *      Author: Gregor Entzian
 *
 */

#include <fstream>
#include <sstream>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <queue>
#include <list>
#include <chrono>
#include <time.h>
#include <vector>
#include <future>
#include <thread>
//#include <omp.h>
//basic file operations
#include <ostream>
#include <string>
#include <iostream>
#include <unistd.h>
#include <algorithm>
extern "C" {
#include "pourRNA_cmdl.h"
#include <ViennaRNA/move_set.h>
#include <ViennaRNA/fold.h>
//#include <ViennaRNA/read_epars.h>
#include <ViennaRNA/params/io.h>
#include <ViennaRNA/eval.h>
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/utils/structures.h>
#include <ViennaRNA/neighbor.h>
#include <ViennaRNA/file_formats.h>
}
#include "BIUlibPart/MatrixSparse.hh"
#include "RNAkinetics/RateMatrixUtil.h"
#include "RNAkinetics/RNA_NeighMinFilter.h"
#include "Flooder.h"
#include "TypeID.h"
#include "SC_DotPlot.h"
#include "StatePairCollector.h"

#include "Concurrent_Queue.h"
#include "PairHashMap.h"

#include "BarriersTree.h"

/*! Exception class for exceptions thrown during argument and input parsing.
 */
class ArgException : public std::exception {
public:
//! the error message
std::string errorMsg;
ArgException(std::string errorMsg_) :
  errorMsg(errorMsg_)
{
}


virtual
~ArgException() throw ()
{
}


virtual const char *
what() const throw ()
{
  return errorMsg.c_str();
}
};
// end class

/*
 * ! the maximum number of macrostates is used to normalize the rates.
 * @param z the matrix with transition partition functions.
 * @param done_List the list of handled minima.
 * @param minimaForReverseSearch the minima pairs with minID and structure&energy.
 * @return the maximal number of neighbored macrostates.
 */
size_t
getMaximalNeighborsOfAMacroState(SC_PartitionFunction::Z_Matrix& z,
                                 std::unordered_set<size_t>& done_List,
                                 std::unordered_map<size_t, MyState>& minimaForReverseSearch)
{
  size_t  maxNeighbors            = 0;
  size_t  countNeighborsForState  = 0;

  for (std::unordered_set<size_t>::iterator from = done_List.begin();
       from != done_List.end(); from++) {
    SC_PartitionFunction::PairID basinID = SC_PartitionFunction::PairID(
      *from, *from);
    if (z.find(basinID) != z.end()) {
      countNeighborsForState = 0;
      bool transitionPartitionsumExists = false;
      for (std::unordered_set<size_t>::iterator to = std::next(from, 1);
           to != done_List.end(); to++) {
        SC_PartitionFunction::PairID  transitionID =
          SC_PartitionFunction::PairID(*from, *to);
        SC_PartitionFunction::PairID  reverseTransitionID =
          SC_PartitionFunction::PairID(*to, *from);

        double                        z_transition        = 0;
        double                        z_reverseTransition = 0;
        if (z.find(transitionID) != z.end()) {
          z_transition = z.at(transitionID).getZ();
          if (z_transition > 0.0)
            transitionPartitionsumExists = true;
        }

        if (z.find(reverseTransitionID) != z.end()) {
          z_reverseTransition = z.at(reverseTransitionID).getZ();
          if (z_reverseTransition > 0.0)
            transitionPartitionsumExists = true;
        }

        if (transitionPartitionsumExists)
          countNeighborsForState++;
      }
      if (maxNeighbors < countNeighborsForState)
        maxNeighbors = countNeighborsForState;
    }
  }
  return maxNeighbors;
}


/*
 * Calculate the rate matrix by dividing the transition partition sum by
 * (the maximum number of neighbors times the basin partition sum).
 * @param z the matrix with transition partition functions.
 * @param done_List the list of handled minima.
 * @param minimaForReverseSearch the minima pairs with minID and structure&energy.
 * @return the biu::MatrixSparseC<double> RateMatrix.
 */
biu::MatrixSparseC<double> &
calculateRateMatrix(SC_PartitionFunction::Z_Matrix& z,
                    std::vector<std::pair<size_t, MyState *> >& done_List,
                    bool computeDiagonal = false)
{
  size_t max_size = 0;
  size_t done_id;
  for (auto done_it = done_List.begin(); done_it != done_List.end(); done_it++) {
    done_id = done_it->first;
    if(done_id > max_size)
      max_size = done_id;
  }
  max_size++;
  biu::MatrixSparseC<double>&                 final_Rate = *new biu::MatrixSparseC<double>(
      max_size, max_size, 0.0);
  // iterate through all elements of done List to produce the final rate matrix
  auto  fromIt  = done_List.begin();
  for (std::uint32_t from = 0; from < done_List.size(); from++, fromIt++ ) {
    // copy according rates
    std::uint32_t  fromOrig = fromIt->first;
    auto toIt = std::next(fromIt, 1); //done_List.begin();

    SC_PartitionFunction::PairID  basinID =
      SC_PartitionFunction::PairID(fromOrig, fromOrig);

    double  z_basin     = 0;
    auto    it_z_basin  = z.find(basinID);
    if (it_z_basin != z.end())
      z_basin = it_z_basin->second.getZ();

    if (z_basin > 0.0) {
      for (std::uint32_t to = from + 1; to < done_List.size(); to++ , toIt++) {
        std::uint32_t toOrig = toIt->first;
        // copy non-diagonal entries
        //if (toOrig != fromOrig) {
          //compute the rate from row-state two column-state
          //in Z-Matrix: from i to j.
          //in final_Rate: from j to i !
          SC_PartitionFunction::PairID  transitionID =
            SC_PartitionFunction::PairID(fromOrig, toOrig);
          SC_PartitionFunction::PairID  reverseTransitionID =
            SC_PartitionFunction::PairID(toOrig, fromOrig);


          double                        z_transition        = 0.0;
          double                        z_reverseTransition = 0.0;
          auto it_transition = z.find(transitionID);
          auto it_reverse_transition = z.find(reverseTransitionID);
          if (it_reverse_transition != z.end())
            z_reverseTransition = it_reverse_transition->second.getZ();
          if (it_transition != z.end())
            z_transition = it_transition->second.getZ();
          //max is important if filters are applied! Otherwise it should be equal.
          z_transition = max(z_transition, z_reverseTransition);

          if(z_transition != 0.0){
            double rate_Val = z_transition / (/*maxNeighNum * */ z_basin);
            if (rate_Val > 0.0) {
              final_Rate.at(fromOrig, toOrig) = rate_Val;
            }

            SC_PartitionFunction::PairID  basinID2    = SC_PartitionFunction::PairID(toOrig, toOrig);
            double                        z_basin2    = 0.0;
            auto                          it_z_basin2 = z.find(basinID2);
            if (it_z_basin2 != z.end()){
              z_basin2 = it_z_basin2->second.getZ();

              if (z_basin2 > 0.0) {
                double rate_Val = z_transition / (/*maxNeighNum * */ z_basin2);
                if (rate_Val > 0.0) {
                  final_Rate.at(toOrig, fromOrig) = rate_Val;
                }
              }
            }
          }
        //}
      }
    }
  }   // end for
  // compute diagonal
  if(computeDiagonal){
    double sum;
    for (auto from_it = done_List.begin(); from_it != done_List.end(); from_it++) {
      sum = 0.0;
      for (auto to_it = done_List.begin(); to_it != done_List.end(); to_it++) {
        if(from_it->first != to_it->first)
          sum += final_Rate.at(from_it->first, to_it->first);
      }
      final_Rate.at(from_it->first, from_it->first) = -sum;
    }
  }
  return final_Rate;
}


struct flooderInputParameter {
  size_t                    MaxToQueue;
  size_t                    MaxToHash;
  double                    MaxEnergy;
  double                    DeltaE;
  size_t                    BasinID;
  MyState                   *CurrentMinimum;
  NeighMinFilter            *Filter;
  Concurrent_Queue<MyState> *DiscoveredMinima;
  double                    TemperatureForBoltzmannWeight;
  double                    GasConstant;
  unsigned int              Move_set;
  StatePairCollector::MapOfMaps                 *All_Saddles;
  int                       MaxBPdist;
  char                      *SourceStructure;
  char                      *TargetStructure;
  double                    MFE; //used for pf scaling
  ~flooderInputParameter()
  {
    if (CurrentMinimum != NULL)
      delete CurrentMinimum;
  }
};

struct flooderOutputAnylse {
  size_t  NumberOfStatesInBasin;
  size_t  NumberOfStatesOutOfBasin; //outer contact surface.
  size_t  TimeForStatesInBasin;
  size_t  TimeForStatesOutOfBasin;
  double  MaxEnergyThreshold;
  size_t  ProcessedStates;
};

struct flooderOutputParameter {
  // ! All neighbored minima.
  PairHashTable::HashTable        LocalThreadMinima;
  // ! Only the indices are filtered.
  std::vector<size_t>             IndicesOfFilteredMinima;
  SC_PartitionFunction::Z_Matrix  PartitionFunctions;
  SC_PartitionFunction            *ScBasin;
  flooderOutputAnylse             Analysis;
  ~flooderOutputParameter()
  {
    if (ScBasin != NULL)
      delete ScBasin;
  }
};

/**
 * ! flood the basin with the given structure "flooderInputParameter->CurrentMinimum".
 *   The transition-partition function and the partition function of the basin will be written to
 *   "flooderOutputParameter->PartitionFunctions".
 *   @param inParameter contains a structure which identifies the basin and filter parameters.
 *   @param outParameter contains the partition function and neigbored minima.
 */
int
floodBasin(vrna_fold_compound_t   *vc,
           flooderInputParameter  *inParameter,
           flooderOutputParameter *outParameter)
{
  size_t                    loopCurrentMinID = inParameter->BasinID;

  PairHashTable::HashTable& localThreadMinima =
    outParameter->LocalThreadMinima;
  // initializing set, which contains the NEighbors indices  of current state
  size_t                    loopCurrentLocalMinID = 0;

  MyState                   *minState = new MyState(*inParameter->CurrentMinimum);

  //set local id to zero because of the statePairCollector id-assign-method.

  localThreadMinima.insert({ *minState, loopCurrentLocalMinID });

  // start walking and computing the Partition Function

  // creating the StatePaircollector to manage the state of Pairs between the neighbor basins
  StatePairCollector spc(loopCurrentLocalMinID, localThreadMinima,
                         outParameter->PartitionFunctions, inParameter->MaxToHash,
                         inParameter->DiscoveredMinima,
                         inParameter->TemperatureForBoltzmannWeight,
                         inParameter->GasConstant,
                         inParameter->MFE,
                         inParameter->Move_set,
                         inParameter->All_Saddles,
                         inParameter->SourceStructure,
                         inParameter->TargetStructure,
                         inParameter->MaxBPdist
                         );

  // set the maximal energy as upper floodlevel.
  double  maxEnergy = min(inParameter->MaxEnergy,
                          (inParameter->DeltaE + minState->energy / 100.0));

  Flooder myFlooder(maxEnergy, inParameter->MaxToQueue, inParameter->Move_set);

  // perform local basin flooding
  // SC_PartitionFunction scBasin;
  myFlooder.floodBasin(vc, *minState, outParameter->ScBasin, &spc);

  outParameter->PartitionFunctions.insert(
    std::pair<SC_PartitionFunction::PairID, SC_PartitionFunction>(
      SC_PartitionFunction::PairID(loopCurrentLocalMinID,
                                   loopCurrentLocalMinID), *outParameter->ScBasin));

  for (auto it = localThreadMinima.begin(); it != localThreadMinima.end();
       it++) {
    if (it->second != loopCurrentLocalMinID) // add only the neighbors.

      outParameter->IndicesOfFilteredMinima.push_back(it->second);
  }

  if (inParameter->Filter != NULL) {
    inParameter->Filter->filter(localThreadMinima,
                                outParameter->PartitionFunctions, loopCurrentLocalMinID,
                                *minState, outParameter->IndicesOfFilteredMinima);
  }

  delete minState;
  //add features for analysis.
  outParameter->Analysis.NumberOfStatesInBasin =
    outParameter->ScBasin->size();
  outParameter->Analysis.NumberOfStatesOutOfBasin =
    spc.getNumberOfOuterStates();
  outParameter->Analysis.TimeForStatesInBasin     = 0;
  outParameter->Analysis.TimeForStatesOutOfBasin  = 0;
  outParameter->Analysis.MaxEnergyThreshold       =
    myFlooder.getFinalEnergyThreshold();
  outParameter->Analysis.ProcessedStates = myFlooder.getProcessedStates();

  return 0;
}


//######################################################## -- MAIN BODY -- ########################################################

/*
 * return the maximal available physical memory.
 */
size_t
getTotalSystemMemoryInBytes()
{
  long  pages     = sysconf(_SC_PHYS_PAGES);
  long  page_size = sysconf(_SC_PAGE_SIZE);

  return pages * page_size;
}


/*
 * The hashmap for the gradientwalk consumes the most memory.
 * To prevent the system from swapping (and in the end from crashing), the maximal size of the hashmap should be estimated.
 */
size_t
estimateMaxToHash(size_t  numberOfHashMaps,
                  size_t  sequenceLength)
{
  //determine estimated maximum maxToHash value without swapping.
  size_t  memInBytes      = getTotalSystemMemoryInBytes();
  size_t  sizeOfHashEntry = sizeof(size_t) + sizeof(MyState)
                            + sizeof(short) * (sequenceLength + 1);
  size_t  numberOfEntries   = memInBytes / sizeOfHashEntry;
  size_t  entriesPerHashMap = numberOfEntries / numberOfHashMaps;

  return entriesPerHashMap;
}

int
merge_results(std::vector<std::pair<flooderInputParameter *,
                                    flooderOutputParameter *> > *threadParameter,
              bool logEnergies,
              ofstream &energyFile,
              std::unordered_set<size_t> &addedMinIDs,
              bool writeDotplot,
              SC_DotPlot::DotPlot &dotplot,
              PairHashTable::HashTable &Minima,
              std::unordered_map<size_t, MyState> &MinimaForReverseSearch,
              std::unordered_map<size_t, size_t> &Minimum_index_and_basin_size,
              SC_PartitionFunction::Z_Matrix &z,
              bool dynamicBestK,
              std::unordered_map<size_t, std::vector<size_t> > &dynamicBestKFilterNeighborList,
              std::list<size_t> &toDo_List,
              std::unordered_set<size_t> &done_List,
              std::unordered_map<MyState,  SC_DotPlot::DotPlot, PairHashTable::PairTableHash, PairHashTable::PairTableEqual> &dot_plot_per_basin,
              bool writeDotplotPerBasin,
              bool verbose)
{
  for (int k = 0; k < threadParameter->size(); ++k) {
    std::pair<flooderInputParameter *, flooderOutputParameter *>  *bothParameter;
    bothParameter = &threadParameter->at(k);

    flooderInputParameter                                         *inParameter =
      bothParameter->first;
    flooderOutputParameter                                        *outParameter =
      bothParameter->second;

    if (verbose) {
      std::cout << "States out of basin: "
                << outParameter->Analysis.NumberOfStatesOutOfBasin
                << std::endl;
      std::cout << "States in basin: "
                << outParameter->Analysis.NumberOfStatesInBasin
                << std::endl;
      std::cout << "Energy Threshold: "
                << outParameter->Analysis.MaxEnergyThreshold << std::endl;
      std::cout << "Max Energy in Basin: "
                << outParameter->ScBasin->getMaxEnergy() << std::endl;
      std::cout << "Processed States by flooding: "
                << outParameter->Analysis.ProcessedStates << std::endl;
    }

    if (logEnergies) {
      if (addedMinIDs.find(inParameter->BasinID) == addedMinIDs.end()) {
        addedMinIDs.insert(inParameter->BasinID);

        const std::vector<int>& energies =
          outParameter->ScBasin->getEnergies();
        int                     energy;
        for (size_t m = 0; m < energies.size(); m++) {
          energy = energies[m];
          energyFile << energy << " ";
        }
      }
    }

    //sum up dotplot basepair partition sums.
    if (writeDotplot) {
      SC_DotPlot::DotPlot tmpDotPlot =
        ((SC_DotPlot *)outParameter->ScBasin)->getBasePairWeightSum();
      for (SC_DotPlot::DotPlot::iterator it = tmpDotPlot.begin();
           it != tmpDotPlot.end(); ++it) {
        SC_DotPlot::DotPlot::iterator entryIt = dotplot.find(it->first);
        if (entryIt != dotplot.end())
          dotplot[it->first] += it->second;
        else
          dotplot[it->first] = it->second;
      }
    }

    if (writeDotplotPerBasin){
      SC_DotPlot::DotPlot tmpDotPlot = ((SC_DotPlot *)outParameter->ScBasin)->getBasePairWeightSum();
      dot_plot_per_basin[*inParameter->CurrentMinimum] = tmpDotPlot;
    }

    PairHashTable::HashTable&           localThreadMinima =
      outParameter->LocalThreadMinima;
    size_t                              loopCurrentMinID = inParameter->BasinID;

    //add localMinima to total minima and adjust indices in localPartitionfunctions
    //and in indicesOfNeighbors if needed.

    std::unordered_map<size_t, size_t>  mapOldIndexToNewIndex;

    for (auto minIt = localThreadMinima.begin();
         minIt != localThreadMinima.end(); minIt++) {
      size_t          oldIndex = minIt->second;

      const MyState&  localThreadMin = minIt->first;
      if (Minima.find(localThreadMin) == Minima.end()) {
        //if local min is not in total list.
        //add new min and set lowest maximal index.
        size_t lowestMaxIndex = TypeID::value<size_t>();

        Minima.insert({ MyState(localThreadMin), lowestMaxIndex });
        MinimaForReverseSearch.insert({ lowestMaxIndex,
                                        MyState(localThreadMin) });
      }

      size_t newIndex = Minima.at(localThreadMin);

      //assingn total id to localID.
      mapOldIndexToNewIndex.insert({ oldIndex, newIndex });
    }

    size_t count =outParameter->ScBasin->get_state_count();
    Minimum_index_and_basin_size[inParameter->BasinID] = count;

    SC_PartitionFunction::Z_Matrix& localZ =
      outParameter->PartitionFunctions;
    //adjust index in partitionfunction z.
    for (auto it = localZ.begin(); it != localZ.end(); it++) {
      auto new_first_index_it = mapOldIndexToNewIndex.find(it->first.first);
      auto new_second_index_it = mapOldIndexToNewIndex.find(it->first.second);
      if (new_first_index_it
          != mapOldIndexToNewIndex.end()
          && new_second_index_it
          != mapOldIndexToNewIndex.end()) {
        size_t                        newFirstIndex = new_first_index_it->second;
        size_t                        newSecondIndex = new_second_index_it->second;
        SC_PartitionFunction::PairID  localPairID =
          SC_PartitionFunction::PairID(newFirstIndex,
                                       newSecondIndex);
        z.insert({ localPairID, it->second });
      }
    }

    //fill dynamic Best-k filter list.
    if (dynamicBestK) {
      if (dynamicBestKFilterNeighborList.find(loopCurrentMinID)
          == dynamicBestKFilterNeighborList.end()) {
        vector<size_t>  allNeighboredMinima;
        size_t          currentMinID = localThreadMinima.at(
          *inParameter->CurrentMinimum);
        for (auto localMinIt = localThreadMinima.begin();
             localMinIt != localThreadMinima.end(); ++localMinIt)
          if (localMinIt->second != currentMinID)
            allNeighboredMinima.push_back(localMinIt->second);

        // creating an greaterByZ Object as third parameter of std:: sort method
        GreaterByZMatrix greaterByZ(localZ, currentMinID);

        // sort bz Z value
        std::sort(allNeighboredMinima.begin(),
                  allNeighboredMinima.end(), greaterByZ);

        for (int nm = 0; nm < allNeighboredMinima.size(); ++nm)
          allNeighboredMinima[nm] =
            mapOldIndexToNewIndex[allNeighboredMinima.at(nm)];

        dynamicBestKFilterNeighborList.insert({ loopCurrentMinID,
                                                allNeighboredMinima });
      }
    }

    //  add all neighbored minima not part of DONE list to toDO_list
    for (size_t i = 0; i < outParameter->IndicesOfFilteredMinima.size(); i++) {
      size_t index = mapOldIndexToNewIndex.at(
        outParameter->IndicesOfFilteredMinima.at(i));
      // if the element at the index position is not already in done List, add it to toDo_list
      if (done_List.find(index) == done_List.end()
          && std::find(toDo_List.begin(), toDo_List.end(), index)
          == toDo_List.end())
        toDo_List.push_back(index);
    } // end for

    //clean up (thread specific minima lists)
    localThreadMinima.clear();
    mapOldIndexToNewIndex.clear();
    localZ.clear();

    delete inParameter;
    delete outParameter;
  } //end merge results.
  return 0;
}


int
threadFunction(vrna_fold_compound_t                               *vc,
               std::vector<std::pair<flooderInputParameter *,
                                     flooderOutputParameter *> >  *threadParameter)
{
  for (int k = 0; k < threadParameter->size(); ++k) {
    flooderInputParameter   *inParameter  = threadParameter->at(k).first;
    flooderOutputParameter  *outParameter = threadParameter->at(k).second;
    floodBasin(vc, inParameter, outParameter);
  }
  return 0;
}


typedef struct less_second_ {
  typedef std::pair<size_t, MyState *> type;
  bool
  operator ()(type const& a,
              type const& b) const
  {
    // check if
    // - smaller energy or
    // - equal energy and smaller string representation
    return (a.second->getEnergy() < b.second->getEnergy())
           || (a.second->getEnergy() == b.second->getEnergy()
               && StructureUtils::IsSmaller(a.second->structure,
                                            b.second->structure));
  }
} less_second;


std::string read_sequence_from_stdin(){
   char          *rec_sequence, *rec_id, **rec_rest;
   unsigned int  rec_type;
   unsigned int  read_opt = 0;

   rec_id          = NULL;
   rec_rest        = NULL;

   rec_type = vrna_file_fasta_read_record(&rec_id,
                                          &rec_sequence,
                                          &rec_rest,
                                          stdin,
                                          read_opt);
   if (rec_type & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT)){
        free(rec_sequence);
        free(rec_id);
        free(rec_rest);
        return "";
   }
   std::string result_sequence = "";
   if(rec_type & VRNA_INPUT_SEQUENCE)
     result_sequence = std::string(rec_sequence);
   free(rec_sequence);
   free(rec_id);
   free(rec_rest);
  return result_sequence;
}

int
main(int  argc,
     char **argv)
{
  StatePairCollector::MapOfMaps                                all_saddles;

  //init stopwatch:
  std::chrono::time_point<std::chrono::system_clock>  start, end;
  start = std::chrono::system_clock::now();

  int                                                 return_Value = 0;
  //fill the data

  // output stream to display invalidity of some given input
  std::ostringstream                                  error_Mas;
  // output container
  std::ostream                                        *out = &std::cout;
  // file to store the transition matrix.
  std::ofstream                                       *transOutFile = NULL;
  // output container for the transition matrix.
  std::ostream                                        *transOut = &std::cout;
  // the transition matrix file name.
  std::string                                         rateMatrixFileName;
  // the output saddle file
  std::string                                         saddleFileName;
  // the maximal number of elements the underlying queue of the Flooder is allowed to contain.
  size_t                                              maxToQueue = UINT32_MAX;
  // maximum number of states to be hashed for each gradient walk.
  size_t                                              maxToHash = UINT32_MAX;
  // the maximum energy that a state is allowed to have to be considered by the flooder.
  double                                              maxEnergy = 5;
  // the maximum energy difference that states in a basin can have w.r.t. the local minimum.
  double                                              deltaE = 65536;
  // value for best K filter. -1 = do not filter.
  size_t                                              filterValueK      = 0;
  bool                                                enableBestKFilter = false;
  // increases K if the MFE structure is not in the done list and start the exploration again.
  bool                                                dynamicBestK = false;

  // value for max-neigh-e filter.
  double                                              filterValueE          = 0;
  bool                                                enableMax_Neigh_E_Filter = false;

  //max base pair distance filter
  int                                                 maxBPdist = -1;

  // file to store all energies.
  std::string                                         energyFileName = "";
  // binary rates file
  std::string                                         binary_rates_file = "";
  // this bool tells the flooder if it should store energies. (not recommended for large sequences!)
  bool                                                logEnergies = false;
  /* parameter for writing a postscript-file with a DotPlot
   * (it contains the base pair probabilities for all structures in the (filtered) energylandscape).*/
  bool                                                writeDotplot = false;
  //parameter for writing a postscript-file with a DotPlot.
  std::string                                         dotPlotFileName = "";
  // parameter for writing a postscript dotplot per basin.
  bool                                                writeDotplotPerBasin = false;
  // file name prefix for dotplots per basin.
  std::string                                         dotPlotPerBasinFileName = "";

  std::string                                         barrierTreeFileName = "";
  // output file that contains the minh representative and the elementary gradient basins which are mapped to them.
  std::string                                         minh_mapping_file = "";
  // file name for barriers-like output
  std::string                                         barriers_prefix;
  // file to store all partition functions.
  std::string                                         partitionFunctionFileName = "";
  // write partitionfunctions to a file if the name is given as parameter.
  bool                                                writePartitionFunctions = false;

  // parameter for temperature in celsius.
  double                                              temperatureForEnergy = 37;
  // paramter for the temperature for the Boltzmann weight in celsius.
  double                                              temperatureForBoltzmannWeight = 37;
  
  double                                              gas_constant = 0.00198717; /* in [kcal/(K*mol)] */

  // compute the diagonal of the rate matrix (can be skipped because some post-processing tools like treekin recompute it)
  bool                                                computeDiagonal = true;

  // How to treat \"dangling end\" energies for bases adjacent to helices in free ends and multi-loops
  int                                                 danglingEnd = 2;
  // The move set type (see neighbor.h vrna_package)
  unsigned int                                        move_set  = 0;
  bool                                                basin_size = false;
  bool                                                verbose   = false;

  // parameter: path to file for energy model that is placed in the share/misc folder of this build.
  // (set by th "M" parameter: 0=Turner2004,1=Turner1999,2=Andronescu2007)
  std::string                                         energyModelFile = "";
  /*
   * estimate the maxToHash parameter. The estimate depends on the available memory and the number of threads.
   * This estimation may slow down the computation, but it will be faster than with a swapping system.
   * The main purpose for this parameter is, that it reduces the probability that the system will crash if too much
   * memory is used with large sequences.
   */
  bool                    dynamicMaxToHash = false;

  double minh = -1;
  std::vector<size_t> dynamic_minh_max_states;


  //---------------------------- Parsing the Parameters---------------------------------

  struct pourRNA_args_info args_info;
  // check if there is any unparseable argument
  if (pourRNA_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  try {
    // parse sequence
    std::string rnaSequence = "";
    if (!args_info.sequence_given) {

      rnaSequence = read_sequence_from_stdin();
      if(rnaSequence == "")
        throw ArgException("No RNA sequence given");
    } else {
      rnaSequence = std::string(args_info.sequence_arg);
      if (!StructureUtils::IsValidSequence(rnaSequence)) {
        error_Mas << "Sequence " << rnaSequence << " is no valid Sequence ";
        throw ArgException(error_Mas.str());
      }
    }

    // parse structure
    std::string rna_start_str = "";
    if (!args_info.start_structure_given) {
      // initialize open chain structure
      rna_start_str = std::string(rnaSequence.size(), '.');
    } else {
      rna_start_str = std::string(args_info.start_structure_arg);
      if (!StructureUtils::IsValidStructure(rna_start_str)) {
        error_Mas << " start Structure " << rna_start_str
                  << " is no valid structure";
        throw ArgException(error_Mas.str());
      }
    }

    std::string rna_final_str = "";
    if (args_info.final_structure_given) {
      rna_final_str = args_info.final_structure_arg;
      if (!StructureUtils::IsValidStructure(rna_final_str)
          || rna_final_str.size() == 0) {
        error_Mas << " final Structure " << rna_final_str
                  << " is no valid structure";
        throw ArgException(error_Mas.str());
      }
      if (args_info.max_threads_given && args_info.max_threads_arg > 1){
        error_Mas << "the final structure is not allowed with multiple threads. Please set --max-threads=1";
        throw ArgException(error_Mas.str());
      }
    }

    int maxBP_add = 0;
    if (args_info.max_bp_dist_add_given) {
      maxBP_add = args_info.max_bp_dist_add_arg;
      //if (maxBPdist < 0)
      //  throw ArgException("The base pair distance has to be positive!");
      if(!args_info.final_structure_given || !args_info.start_structure_given)
        throw ArgException("Error: max-bp-dist-add expects that you also set a start structure and a final structure!");
    }

    std::list<std::string> start_structure_list;
    if (args_info.start_structure_file_given) {
      std::string start_structure_file = args_info.start_structure_file_arg;
      if (start_structure_file.size() != 0) {
        std::ifstream infile(start_structure_file);
        std::string   line;
        while (std::getline(infile, line))
          if (StructureUtils::IsValidStructure(line)){
              start_structure_list.push_back(line);
          }
      }
    }

    // container for all partition functions to generate
    // indices are according to minima container
    // Z(i,i) will hold the basins partition function
    // Z(i,j) will hold the partition function of the saddle points between minimum i and j
    SC_PartitionFunction::Z_Matrix      z;

    // set of all minima found so far
    PairHashTable::HashTable            Minima;
    std::unordered_map<size_t, MyState> MinimaForReverseSearch;

    // set best K filter if required
    if (args_info.filter_best_k_given) {
      filterValueK      = args_info.filter_best_k_arg;
      enableBestKFilter = true;
    }

    // set dynamic best K filter if required
    if (args_info.dynamic_best_k_given)
      dynamicBestK = true;

    // set Energy filter if required
    if (args_info.max_neigh_e_given) {
      filterValueE          = args_info.max_neigh_e_arg;
      enableMax_Neigh_E_Filter = true;
    }

    // set maxToStore variable for the states in each gradient basin.
    if (args_info.max_to_queue_given)
      maxToQueue = args_info.max_to_queue_arg;

    // set maxToHash variable for the states to be hashed in a gradient walk.
    if (args_info.max_to_hash_given)
      maxToHash = args_info.max_to_hash_arg;

    // set dynamicMaxToHash variable for estimating the maximal number of states to be hashed in a gradient walk,
    // by considering the maximal available physical memory.
    if (args_info.dynamic_max_to_hash_given)
      dynamicMaxToHash = true;

    // set the maximum energy that a state is allowed to have to be considered by the flooder.
    if (args_info.max_energy_given)
      maxEnergy = args_info.max_energy_arg;

    // set the maximum energy difference that states in a basin can have w.r.t. the local minimum.
    if (args_info.delta_e_given)
      deltaE = args_info.delta_e_arg;

    // set max_threads for parallelized computation.
    int maxThreads = 1;
    if (args_info.max_threads_given) {
      maxThreads = args_info.max_threads_arg;
      if (maxThreads < 1)
        throw ArgException("you should use at least one thread.");
    }

    if (args_info.energy_file_given) {
      // set output file name for energies.
      energyFileName  = args_info.energy_file_arg;
      logEnergies     = true;
    }

    if (args_info.binary_rates_file_given)
      // set output file name for energies.
      binary_rates_file = args_info.binary_rates_file_arg;

    //partitionFunctions
    if (args_info.partition_functions_given) {
      // set output file name for partitionfunctions.
      partitionFunctionFileName = args_info.partition_functions_arg;
      writePartitionFunctions   = true;
    }

    //dotPlot
    if (args_info.dot_plot_given) {
      // set output file name for dotPlot.
      dotPlotFileName = args_info.dot_plot_arg;
      writeDotplot    = true;
    }

    if (args_info.dot_plot_per_basin_given) {
      dotPlotPerBasinFileName = args_info.dot_plot_per_basin_arg;
      writeDotplotPerBasin    = true;
    }

    if(args_info.barrier_tree_file_arg)
      barrierTreeFileName = std::string(args_info.barrier_tree_file_arg);

    if(args_info.barriers_like_output_given)
      barriers_prefix = std::string(args_info.barriers_like_output_arg);

    if(args_info.saddle_file_given)
      saddleFileName = std::string(args_info.saddle_file_arg);

    if (args_info.transition_prob_given) {
      // set output stream
      rateMatrixFileName = std::string(args_info.transition_prob_arg);
      if (rateMatrixFileName.compare("STDOUT") == 0) {
        transOut = out;
      } else {
        transOutFile = new std::ofstream(
          rateMatrixFileName.c_str(),
          std::ofstream::out);
        if (!transOutFile->is_open()) {
          std::ostringstream oss;
          oss << "cannot open output file '"
              << args_info.transition_prob_arg << "'";
          throw ArgException(oss.str());
        }

        transOut = transOutFile;
      }
    }

    basin_size = args_info.basin_size_flag;

    if (args_info.dangling_end_given) {
      // set dangling ends.
      danglingEnd = args_info.dangling_end_arg;
      if (danglingEnd < 0) {
        std::ostringstream oss;
        oss << "Error: dangling end size must be >= 0 !";
        throw ArgException(oss.str());
      }
    }

    if (args_info.skip_diagonal_given)
      computeDiagonal = false;
    else
      computeDiagonal = true;

    if (args_info.temperature_given) {
      // set temperatureForEnergy.
      temperatureForEnergy = args_info.temperature_arg;
      if (!(args_info.boltzmann_temp_given))
        temperatureForBoltzmannWeight = temperatureForEnergy;
    }
    
    if (args_info.gas_constant_given) {
      gas_constant = args_info.gas_constant_arg;
    }

    if (args_info.boltzmann_temp_given)
      // set temperature for the boltzmann weight.
      temperatureForBoltzmannWeight = args_info.boltzmann_temp_arg;

    int tmp_move_set = 0;
    if (args_info.move_set_given)
      // set the move set integer
      tmp_move_set = args_info.move_set_arg;

    // convert the integer to the vrna_type
    switch (tmp_move_set) {
      case 0:
        move_set = VRNA_MOVESET_DEFAULT;
        break;
      case 1:
        move_set = VRNA_MOVESET_DEFAULT | VRNA_MOVESET_SHIFT;
        break;
      case 2:
        move_set = VRNA_MOVESET_NO_LP;
        break;
      default:
        move_set = VRNA_MOVESET_DEFAULT;
    }


    if (args_info.energy_model_given) {
      // set energyModel.
      int energyModelNumber = args_info.energy_model_arg;
      if (energyModelNumber >= 0 || energyModelNumber <= 2) {
        // Get the last position of '/'
        char        buf[UINT16_MAX] = "";
        ssize_t     res             = readlink("/proc/self/exe", buf, UINT16_MAX);
        if (res == -1)
          fprintf(stderr, "Error: could not read /proc/self/exe");

        std::string aux(buf); //(argv[0]);
        // get '/' or '\\' depending on unix/mac or windows.
#if defined(_WIN32) || defined(WIN32)
        int         pos = aux.rfind('\\');
#else
        int         pos = aux.rfind('/');
#endif
        std::string shareFolder = "misc";
#ifdef DEBUG
        shareFolder = "misc";
#else
        shareFolder = "share/pourRNA";
#endif

        // Get the path of the executable of this program.
        int         progPosition  = aux.substr(0, pos).rfind('/');
        std::string path          = aux.substr(0, progPosition + 1);

        switch (energyModelNumber) {
          case 0:
            energyModelFile = path + shareFolder
                              + "/rna_turner2004.par";
            //printf("%s \n",energyModelFile.c_str());
            read_parameter_file(energyModelFile.c_str());
            /*TODO: find out why the default setting of viennaRNA 2.3.5
             *      produces different results than the results after
             *      reading the default file
             */
            break;
          case 1:
            energyModelFile = path + shareFolder
                              + "/rna_turner1999.par";
            read_parameter_file(energyModelFile.c_str());
            break;
          case 2:
            energyModelFile = path + shareFolder
                              + "/rna_andronescu2007.par";
            read_parameter_file(energyModelFile.c_str());
            break;
        }
      } else {
        throw ArgException(
                "Error: Only values between 0 and 2 are allowed for parameter M (energy model)!");
      }
    }

    std::string mapped_structures_filename = "";
    if(args_info.map_structures_given)
      mapped_structures_filename = std::string(args_info.map_structures_arg);

    if (args_info.minh_given){
      minh = args_info.minh_arg;
    }

    if (args_info.dynamic_minh_given) {
      for (size_t i = 0; i < args_info.dynamic_minh_given; ++i)
        dynamic_minh_max_states.push_back((size_t)args_info.dynamic_minh_arg[i]);
      std::sort(dynamic_minh_max_states.begin(), dynamic_minh_max_states.end(), std::greater<size_t>());
      if (minh < 0)
        minh = 0.0;
    }

    if (args_info.verbose_given)
      verbose = true;

    vrna_md_t md;
    vrna_md_set_default(&md);
    md.circ = 0;
    //md.uniq_ML = 1; /* in case we need M1 arrays */
    md.noLP         = 0;
    md.gquad        = 0;
    md.noGU         = 0;
    md.compute_bpp  = 0;
    md.temperature  = temperatureForEnergy;
    md.dangles      = danglingEnd;
    char                  *sequence = (char *)rnaSequence.c_str();
    vrna_fold_compound_t  *vc       = vrna_fold_compound(sequence, &md,
                                                         VRNA_OPTION_MFE); // | VRNA_OPTION_PF);

    // set best K filter if required
    // create an filter Object type NeighMinFilter as first filtering Technique
    NeighMinFilter        *neighbor_K_Filter = NULL;

    // create an filter Object type NeighMinFilter as first filtering Technique
    NeighMinFilter        *neighbor_E_Filter = NULL;

    // create an Combinator object to apply both filtering Techniques in specific order
    NeighMinFilter        *neighborCombinator = NULL;

    // final overall filter to be used
    NeighMinFilter        *neighborFilter = NULL;
    if (enableBestKFilter) {
      neighbor_K_Filter = new NeighMin_K_Filter(filterValueK);
      neighborFilter    = neighbor_K_Filter;
    }

    // set Energy filter if required
    if (enableMax_Neigh_E_Filter) {
      neighbor_E_Filter = new NeighMin_E_Filter(filterValueE);
      neighborFilter    = neighbor_E_Filter;
    }

    // set final Filter to be combination of both Filters
    if (neighbor_K_Filter != NULL && neighbor_E_Filter != NULL) {
      neighborCombinator = new NeighMinCombinator(*neighbor_E_Filter,
                                                  *neighbor_K_Filter);
      neighborFilter = neighborCombinator;
    }

    //do flooding

    // start walking to find gradient basin minimum for start state
    short                                               *rnaStructurePT = vrna_ptable(
      rna_start_str.c_str());
    //printf("%s\n", rna_start_str.c_str());
    int                                                 energy = vrna_eval_structure_pt(
      vc,
      rnaStructurePT);
    MyState                                             startState(energy, rnaStructurePT);

    //init stopwatch for initial walk:
    std::chrono::time_point<std::chrono::system_clock>  startInitialWalk,
                                                        endInitialWalk;
    if (verbose)
      startInitialWalk = std::chrono::system_clock::now();

    MyState                                             *startStateMinimum = WalkGradientHashed(
      move_set,
      maxToHash).walk(vc,
                      startState);

    if (verbose) {
      //print stopwatch for initial walk:
      endInitialWalk = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds_InitialWalk =
        endInitialWalk - startInitialWalk;
      std::cout << "Initial Walk elapsed time: "
                << elapsed_seconds_InitialWalk.count() << "s\n";
    }

    char    *startStructureMinimum  = vrna_db_from_ptable(startState.getStructure());
    short   *rnaFinalStructurePT    = NULL;
    MyState *finalStructureMinimum  = NULL;
    char    *endStructureMinimum    = NULL;
    if (rna_final_str != "") {
      rnaFinalStructurePT = vrna_ptable(rna_final_str.c_str());
      int     energy_final = vrna_eval_structure_pt(vc, rnaFinalStructurePT);
      MyState finalState(energy_final, rnaFinalStructurePT);
      finalStructureMinimum = WalkGradientHashed(move_set, maxToHash).walk(vc, finalState);
      endStructureMinimum   = vrna_db_from_ptable(finalStructureMinimum->getStructure());
    }

    // add current minimum to set of minima and store its index
    size_t                  currentMinID = TypeID::value<size_t>();
    Minima.insert({ MyState(*startStateMinimum), currentMinID });
    MinimaForReverseSearch.insert({ currentMinID, MyState(*startStateMinimum) });

    std::unordered_map<size_t, size_t> minimum_index_and_basin_size;

    ///// for tests only -> write file with all Energies ///
    std::unordered_set<size_t> addedMinIDs;
    ofstream                energyFile;
    if ((energyFileName.size() > 0) & logEnergies) {
      energyFile.open(energyFileName);
      energyFile << "/* energy in 10kcal/mol*/.\n";
    }

    //initialize DotPlot file:
    SC_DotPlot::DotPlot dotplot;
    ////////

    std::unordered_map<MyState,  SC_DotPlot::DotPlot, PairHashTable::PairTableHash, PairHashTable::PairTableEqual> dot_plot_per_basin;

    ///////////   ITERATE EXPLORATIVE LOCAL FLOODING  /////////////

    // initializing the DO list as queue, which contain all states that should be processed
    std::list<size_t>           toDo_List;

    // initializing the List of all states, which had already visited
    std::unordered_set<size_t>  done_List;

    // insert the index of current minima to toDo_list
    toDo_List.push_back(currentMinID);

    if(args_info.max_bp_dist_add_given){
      maxBPdist = vrna_bp_distance(startStructureMinimum, endStructureMinimum);
      maxBPdist += maxBP_add;
    }

    ofstream mapped_structures_file;
    if (start_structure_list.size() > 0  && mapped_structures_filename.compare("") != 0){
      mapped_structures_file.open(mapped_structures_filename);
      mapped_structures_file << "from input, to basin" << std::endl;
    }

    for (auto it = start_structure_list.begin(); it != start_structure_list.end(); it++) {
      short   *tmpMinPairTable  = vrna_ptable(it->c_str());
      int     energy            = vrna_eval_structure_pt(vc, tmpMinPairTable);
      MyState start_struct_i(energy, tmpMinPairTable);
      free(tmpMinPairTable);
      MyState *start_struct_min = WalkGradientHashed(move_set, maxToHash).walk(vc, start_struct_i);
      auto    it_min            = Minima.find(*start_struct_min);
      if (it_min == Minima.end()) {
        bool allow_min = true;
        if(args_info.max_bp_dist_add_given){
            int start_dist = vrna_bp_distance(startStructureMinimum, start_struct_min->toString().c_str());
            int final_dist = vrna_bp_distance(endStructureMinimum, start_struct_min->toString().c_str());
            if(start_dist + final_dist > maxBPdist){
              allow_min = false;
            }
        }
        if(allow_min){
          size_t min_id_from_list = TypeID::value<size_t>();
          Minima.insert({ MyState(*start_struct_min), min_id_from_list });
          MinimaForReverseSearch.insert({ min_id_from_list, MyState(*start_struct_min) });
          toDo_List.push_back(min_id_from_list);
        }
      }

      char * mapped_to_basin = vrna_db_from_ptable(start_struct_min->structure);
      mapped_structures_file << it->c_str() << ", " << mapped_to_basin << std::endl;
      free(mapped_to_basin);

      delete start_struct_min;
    }

    if (start_structure_list.size() > 0 && mapped_structures_filename.compare("") != 0){
      mapped_structures_file.close();
    }

    ////// init dynamic-Best-K feature /////
    bool                                              mfeFound            = false;
    bool                                              finalStructureFound = false;
    char                                              *mfeStructure       = (char *)vrna_alloc(
      (vc->length + 1) * sizeof(char));
    double                                            mfeEnergy = (double)vrna_mfe(vc, mfeStructure);

    //init dynamic k-best filter list.
    std::unordered_map<size_t, std::vector<size_t> >  dynamicBestKFilterNeighborList;

    do {
      //start dynamic best k filter loop.
      // go through all the elements of toDo_list to get the Neighbors of that
      while (!toDo_List.empty()) {
        int listSize    = toDo_List.size();
        int numThreads  = maxThreads;  //std::min(listSize, maxThreads);

        if (verbose)
          std::cout << "Todo list size: " << listSize << std::endl;

        //set dynamic maxToHash
        if (dynamicMaxToHash)
          maxToHash = estimateMaxToHash(numThreads,
                                        rnaSequence.length());

        /////// asynchronous method
        /*
         * while todolist not empty
         * while not all threads finished
         *   for x in threads
         *      if x finished
         *        merge results and update todolist
         *        if todolist not empty
         *          start new thread with basin from todolist
         *
         *
         */

        //create threads for local flooding:
        std::vector<
          std::pair<flooderInputParameter *,
                    flooderOutputParameter *> > threadParameterLists[maxThreads];
        Concurrent_Queue<MyState>               discoveredMinimaForEachThread[maxThreads];

        std::future<int>                        v[maxThreads];

        //std::cout << "Waiting..." << std::flush;
        bool                                    finish = false;
        //v[0].wait();
        while (!finish) {
          finish = true;
          for (int index = 0; index < maxThreads; index++) {
            bool res = true;
            try {
              if (v[index].valid()) {
                res = v[index].wait_for(
                  std::chrono::nanoseconds(0))
                      == std::future_status::ready;
                //v[index].get();
              }
            } catch (std::future_error e) {
              //std::cout << "thread not ready!\n" << e.what();
              finish  = false;
              res     = false;
            }
            //std::cout << "" << index << ' ' << res << '\n' << std::flush;
            if (!res) {
              finish = false;
              if(!finalStructureFound){
                //test if new discovered minima are on the stack (pull all available minima)
                //TODO: warning! This kind of asynchronous flooding disables some Filter Options.
                if (!enableBestKFilter && !enableMax_Neigh_E_Filter) {
                  while (!discoveredMinimaForEachThread[index].empty()) {
                    MyState newMin =
                      discoveredMinimaForEachThread[index].pop();

                    //start flooding.
                    if (Minima.find(newMin) == Minima.end()) {
                      //if local min is not in total list.
                      //add new min and set lowest maximal index.
                      size_t lowestMaxIndex = TypeID::value<
                        size_t>();

                      Minima.insert(
                        { newMin, lowestMaxIndex });
                      MinimaForReverseSearch.insert({
                                                      lowestMaxIndex, newMin
                                                    });

                      // if the element at the index position is not already in done List, add it to toDo_list
                      if (done_List.find(lowestMaxIndex)
                          == done_List.end()
                          && std::find(toDo_List.begin(),
                                       toDo_List.end(),
                                       lowestMaxIndex)
                          == toDo_List.end())
                        toDo_List.push_back(lowestMaxIndex);
                    }

                    //discoveredMinimaForEachThread[index].pop_front();
                  }
                }
              }
            } else {
              //finished -> merge result and create a new thread
              std::vector<
                std::pair<flooderInputParameter *,
                          flooderOutputParameter *> > *threadParameter =
                &threadParameterLists[index];
              if (threadParameter->size() > 0) {
                merge_results(threadParameter, logEnergies,
                              energyFile, addedMinIDs, writeDotplot,
                              dotplot, Minima, MinimaForReverseSearch,
                              minimum_index_and_basin_size,
                              z, dynamicBestK,
                              dynamicBestKFilterNeighborList,
                              toDo_List, done_List, dot_plot_per_basin,
                              writeDotplotPerBasin,
                              verbose);

                threadParameter->clear();
              }

              // stop exploration at final structure only if no maximal base pair distance is given!
              if(!args_info.max_bp_dist_add_given){
                //if final structure is in done list --> break;
                if (!finalStructureFound && finalStructureMinimum != NULL) {
                  auto final_min_it = Minima.find(*finalStructureMinimum);
                  if (final_min_it != Minima.end()) {
                    auto done_final_it = done_List.find(final_min_it->second);
                    if (done_final_it != done_List.end()) {
                      finalStructureFound = true;
                      break;
                    }
                  }
                }
              }

              if (finalStructureFound) {
                // do not create new threads if the final structure has been found.
                toDo_List.clear();
                break;
              }

              if (!toDo_List.empty()) {
                finish = false;
                //create a new thread
                int listSize          = toDo_List.size();
                int fractionPerThread = int(
                  ceil(listSize / float(numThreads)));
                fractionPerThread = min(fractionPerThread,
                                        listSize);

                //create parameter lists for each thread
                threadParameter->clear();
                discoveredMinimaForEachThread[index].clear();
                for (int k = 0; k < fractionPerThread; ++k) {
                  size_t minID = toDo_List.front();
                  done_List.insert(minID);
                  toDo_List.pop_front();
                  // create inParameter before starting the thread!
                  // make a list as big as the number of threads and clean the entry before the replacement and in the end.

                  flooderInputParameter *inParameter =
                    new flooderInputParameter();
                  inParameter->BasinID        = minID;
                  inParameter->CurrentMinimum =
                    new MyState(MinimaForReverseSearch.at(minID));
                  inParameter->MaxToHash        = maxToHash;
                  inParameter->MaxToQueue       = maxToQueue;
                  inParameter->Filter           = neighborFilter;
                  inParameter->MaxEnergy        = maxEnergy;
                  inParameter->DeltaE           = deltaE;
                  inParameter->MFE              = mfeEnergy;

                  // check if we can report the minima on-the-fly, or if we have to buffer them for macro-state filtering
                  if (inParameter->Filter == NULL) {
                    inParameter->DiscoveredMinima =
                                          &discoveredMinimaForEachThread[index];
                  }
                  else{
                    /* don't report minima in the concurrent queue, but add them to the todo-list after
                     * the basin has been flooded (see merge()).
                     */
                    inParameter->DiscoveredMinima = NULL;
                  }
                  inParameter->TemperatureForBoltzmannWeight =
                    temperatureForBoltzmannWeight;
                  inParameter->GasConstant = gas_constant;
                  inParameter->Move_set         = move_set;

                  if(args_info.saddle_file_given || args_info.barrier_tree_file_given || args_info.minh_given || args_info.dynamic_minh_given)
                    inParameter->All_Saddles      = &all_saddles;
                  else
                    inParameter->All_Saddles      = NULL;
                  inParameter->MaxBPdist        = maxBPdist;
                  inParameter->SourceStructure  = startStructureMinimum;
                  inParameter->TargetStructure  = endStructureMinimum;

                  flooderOutputParameter *outParameter =
                    new flooderOutputParameter();

                  // check what state collector to instantiate
                  if (writeDotplot || writeDotplotPerBasin) {
                    outParameter->ScBasin = new SC_DotPlot(
                      temperatureForBoltzmannWeight,
                      gas_constant,
                      mfeEnergy,
                      logEnergies);
                  } else {
                    outParameter->ScBasin =
                      new SC_PartitionFunction(
                        temperatureForBoltzmannWeight,
                        gas_constant,
                        mfeEnergy,
                        logEnergies);
                  }

                  threadParameter->push_back(
                    std::pair<flooderInputParameter *,
                              flooderOutputParameter *>(
                      inParameter, outParameter));
                }
                //threadFunction( vc, threadParameter);
                v[index] = std::async(std::launch::async,
                                      threadFunction, vc, threadParameter);
              }
            }
          }
        }
        //clean remaining threadparameters and minima lists
        for (int index = 0; index < maxThreads; index++) {
          std::vector<
            std::pair<flooderInputParameter *,
                      flooderOutputParameter *> > *threadParameter =
            &threadParameterLists[index];
          if (threadParameter->size() > 0)
            threadParameter->clear();

          discoveredMinimaForEachThread[index].clear();
        }
      } // end while

      if (logEnergies)
        energyFile.close();

      /////// Search the MFE structure in the done_list ////
      short *mfeStructurePT = vrna_ptable(mfeStructure);
      for (auto doneIt = done_List.begin(); doneIt != done_List.end();
           ++doneIt) {
        if (StructureUtils::IsEqual(
              MinimaForReverseSearch[*doneIt].structure,
              mfeStructurePT)) {
          mfeFound = true;
          break;
        }
      }
      free(mfeStructurePT);
      mfeStructurePT = NULL;
      /////// K-Filter-MFE-Check ////
      if (dynamicBestK) {
        if (!mfeFound) {
          std::cout
          << "Warning: The mfe structure could not be found for K="
          << filterValueK << std::endl;
          //mfe is not in done_list --> increase K-Best filter.
          ++filterValueK;

          if (neighborCombinator != NULL)
            delete neighborCombinator;

          if (neighbor_K_Filter != NULL)
            delete neighbor_K_Filter;

          neighbor_K_Filter = new NeighMin_K_Filter(filterValueK);
          neighborFilter    = neighbor_K_Filter;
          // set final Filter to be combination of both Filters
          if (neighbor_K_Filter != NULL && neighbor_E_Filter != NULL) {
            neighborCombinator = new NeighMinCombinator(
              *neighbor_E_Filter, *neighbor_K_Filter);
            neighborFilter = neighborCombinator;
          }

          //add the k-th neighbor for each macrostate to the todolist.
          for (auto dynBKit = dynamicBestKFilterNeighborList.begin();
               dynBKit != dynamicBestKFilterNeighborList.end();
               ++dynBKit) {
            if (dynBKit->second.size() >= filterValueK) {
              size_t kNeighbor = dynBKit->second.at(
                filterValueK - 1);
              if (done_List.find(kNeighbor) == done_List.end()
                  && std::find(toDo_List.begin(),
                               toDo_List.end(), kNeighbor)
                  == toDo_List.end())
                toDo_List.push_back(kNeighbor);
            }
          }
        }
      }

      /////////////////////////////
    } while (dynamicBestK && !mfeFound && !toDo_List.empty()); //end dynamic best k filter loop.

    if (dynamicBestK) {
      if (mfeFound) {
        std::cout
        << "DynamicBestK: The mfe structure could be found for K="
        << filterValueK << std::endl;
      } else {
        std::cout
        << "Warning: The mfe structure could not be found by increasing K. Last K="
        << filterValueK << std::endl;
      }
    } else {
      if (!mfeFound)
        std::cout << "Warning: The mfe structure could not be found!"
                  << std::endl;
    }

    // sortedMinimaIDs will be sorted by energy, it contains the global id from the done list and the minimum pointer.
    std::vector<std::pair<size_t, MyState *> > sortedMinimaIDs;
    for (std::unordered_set<size_t>::const_iterator fromIt = done_List.begin(); fromIt != done_List.end(); fromIt++)
      sortedMinimaIDs.push_back(
        std::pair<size_t, MyState *>(*fromIt, (MyState *)&MinimaForReverseSearch.at(*fromIt)));

    std::sort(sortedMinimaIDs.begin(), sortedMinimaIDs.end(), less_second());

    PairHashTable::HashTable *sorted_min_and_output_ids = new PairHashTable::HashTable();
    for (size_t c = 0; c < sortedMinimaIDs.size(); c++) {
      (*sorted_min_and_output_ids)[*(sortedMinimaIDs[c].second)] = c;
    }

    // filter minh -- requires saddles.
    BarriersTree bt;
    std::vector<saddle_t> minimal_saddle_list;
    if (args_info.barrier_tree_file_given || args_info.minh_given || args_info.dynamic_minh_given)
      minimal_saddle_list = bt.create_minimal_saddle_list(sortedMinimaIDs, *sorted_min_and_output_ids, all_saddles);

    SC_PartitionFunction::Z_Matrix z_minh;

    bool change_output_filenames = false;
    std::string tmp_dotPlotFileName,
                tmp_dotPlotPerBasinFileName,
                tmp_barrierTreeFileName,
                tmp_rateMatrixFileName,
                tmp_partitionFunctionFileName,
                tmp_saddleFileName,
                tmp_barriers_prefix,
                tmp_binary_rates_file,
                tmp_minh_mapping_file;

    SC_PartitionFunction::Z_Matrix tmp_z;
    std::vector<saddle_t> tmp_minimal_saddle_list;
    std::vector<std::pair<size_t, MyState *>> tmp_sortedMinimaIDs;
    PairHashTable::HashTable* tmp_sorted_min_and_output_ids = NULL;
    std::unordered_map<MyState, SC_DotPlot::DotPlot, PairHashTable::PairTableHash, PairHashTable::PairTableEqual> tmp_dot_plot_per_basin;


    if (dynamic_minh_max_states.size() > 1){
      change_output_filenames = true;
      if(args_info.dot_plot_given)
        tmp_dotPlotFileName = dotPlotFileName;
      if(args_info.dot_plot_per_basin_arg)
        tmp_dotPlotPerBasinFileName = dotPlotPerBasinFileName;
      if(args_info.barrier_tree_file_given)
        tmp_barrierTreeFileName = barrierTreeFileName;
      if(args_info.transition_prob_given && transOutFile != NULL)
        tmp_rateMatrixFileName = rateMatrixFileName;
      if(args_info.partition_functions_given)
        tmp_partitionFunctionFileName = partitionFunctionFileName;
      if(args_info.saddle_file_given)
        tmp_saddleFileName = saddleFileName;
      if(args_info.barriers_like_output_given)
        tmp_barriers_prefix = barriers_prefix;
      if (args_info.binary_rates_file_given)
        tmp_binary_rates_file = binary_rates_file;
      if(args_info.minh_given || args_info.dynamic_minh_given)
        tmp_minh_mapping_file = minh_mapping_file;

      tmp_z.insert(z.begin(), z.end());
      tmp_minimal_saddle_list.assign(minimal_saddle_list.begin(), minimal_saddle_list.end());
      tmp_sortedMinimaIDs.assign(sortedMinimaIDs.begin(), sortedMinimaIDs.end());
      tmp_sorted_min_and_output_ids = new PairHashTable::HashTable();
      tmp_sorted_min_and_output_ids->insert(sorted_min_and_output_ids->begin(), sorted_min_and_output_ids->end());
      tmp_dot_plot_per_basin.insert(dot_plot_per_basin.begin(),dot_plot_per_basin.end());
    }

    std::vector<std::pair<MyState*, std::vector<MyState*>>> minh_representatives_and_basins;
    do{
      size_t minh_max_states;
      if (dynamic_minh_max_states.size() > 0){
        minh_max_states  = dynamic_minh_max_states[0];
        dynamic_minh_max_states.erase(dynamic_minh_max_states.begin());
      }
      if (args_info.minh_given || args_info.dynamic_minh_given){
        double adjusted_minh = minh;
        if (args_info.dynamic_minh_given){
          //Additionally change output file names if several dynamic minh thresholds will be applied.
          if (change_output_filenames){
            z.clear();
            z.insert(tmp_z.begin(), tmp_z.end());
            minimal_saddle_list.clear();
            minimal_saddle_list = std::vector<saddle_t>(tmp_minimal_saddle_list);
            sortedMinimaIDs.clear();
            sortedMinimaIDs = std::vector<std::pair<size_t, MyState *>>(tmp_sortedMinimaIDs);
            sorted_min_and_output_ids->clear();
            sorted_min_and_output_ids->insert(tmp_sorted_min_and_output_ids->begin(), tmp_sorted_min_and_output_ids->end());
            dot_plot_per_basin.clear();
            dot_plot_per_basin.insert(tmp_dot_plot_per_basin.begin(), tmp_dot_plot_per_basin.end());

            adjusted_minh = bt.determin_optimal_min_h(minh_max_states, minimal_saddle_list, sortedMinimaIDs);
            if (verbose){
              fprintf(stdout, "The dynamic minh for maximal %ld states is: %.2f\n", minh_max_states, adjusted_minh);
            }

            char minh_string[50];
            sprintf(minh_string,"%.2f",adjusted_minh);
            std::string appendix = "_dynamic_minh_max_"+std::to_string(minh_max_states)+"_minh_"+minh_string;
            if(args_info.dot_plot_given){
              dotPlotFileName = tmp_dotPlotFileName;
              dotPlotFileName.append(appendix);
            }
            if(args_info.dot_plot_per_basin_arg){
              dotPlotPerBasinFileName = tmp_dotPlotPerBasinFileName;
              dotPlotPerBasinFileName.append(appendix);
            }
            if(args_info.barrier_tree_file_given){
              barrierTreeFileName = tmp_barrierTreeFileName;
              barrierTreeFileName.append(appendix);
            }
            if(args_info.transition_prob_given && transOutFile != NULL){
              rateMatrixFileName = tmp_rateMatrixFileName;
              rateMatrixFileName.append(appendix);
              transOutFile->close();
              transOutFile = new std::ofstream(
                        rateMatrixFileName.c_str(),
                        std::ofstream::out);
                      if (!transOutFile->is_open()) {
                        std::ostringstream oss;
                        oss << "cannot open output file '"
                            << args_info.transition_prob_arg << "'";
                        throw ArgException(oss.str());
                      }
                      transOut = transOutFile;
            }
            if(args_info.partition_functions_given){
              partitionFunctionFileName = tmp_partitionFunctionFileName;
              partitionFunctionFileName.append(appendix);
            }
            if(args_info.saddle_file_given){
              saddleFileName = tmp_saddleFileName;
              saddleFileName.append(appendix);
            }
            if(args_info.barriers_like_output_given){
              barriers_prefix = tmp_barriers_prefix;
              barriers_prefix.append(appendix);
            }
            if (args_info.binary_rates_file_given){
              binary_rates_file = tmp_binary_rates_file;
              binary_rates_file.append(appendix);
            }
            if(args_info.minh_given || args_info.dynamic_minh_given){
              minh_mapping_file = tmp_minh_mapping_file;
              minh_mapping_file.append(appendix);
            }

            /**
             * TODO: write it not only once for the landscape.
            if(args_info.energy_file_given){
              energyFileName.append(appendix);
            }
            */
            /**
             * Is only written once for the whole landscape.
             * However, it the gradient basins can be mapped
             * with an additional basin mapping file for the
             * minh output.
            if(args_info.map_structures_given){
              mapped_structures_filename.append(appendix);
            }
            */
          }
        }

        std::unordered_map<size_t, size_t> merge_state_map = bt.filter_minh(minimal_saddle_list, sortedMinimaIDs, adjusted_minh);
        if (args_info.minh_given || args_info.dynamic_minh_given){
          // merge partition functions and dotplots and assign new IDs!
          std::unordered_map<size_t, std::unordered_set<size_t>> representatives_and_clustered_ids;
          std::unordered_map<size_t, size_t> merged_min_to_representative;
          for(size_t i =0; i < sortedMinimaIDs.size(); i++){
            size_t state_index = i;
            auto it_merged = merge_state_map.find(state_index);
            if(it_merged != merge_state_map.end()){
              // merge this state into its deeper neighbors.
              size_t cluster_root_state_index;
              std::unordered_set<size_t> states_to_cluster;
              size_t max_rounds = merge_state_map.size(); //prevent endless loops.
              while(it_merged != merge_state_map.end() && max_rounds > 0){
                  cluster_root_state_index = it_merged->second;
                  it_merged = merge_state_map.find(cluster_root_state_index);
                  states_to_cluster.insert(cluster_root_state_index);
                  max_rounds--;
              }
              merged_min_to_representative[state_index] = cluster_root_state_index;
              representatives_and_clustered_ids[cluster_root_state_index].insert(state_index);
            }
          }

          // update start and final state
          auto start_min_it = sorted_min_and_output_ids->find(*startStateMinimum);
          if (start_min_it != sorted_min_and_output_ids->end()){
            auto merged_it = merged_min_to_representative.find(start_min_it->second);
            if (merged_it != merged_min_to_representative.end()){
              startStructureMinimum = strcpy(startStructureMinimum, sortedMinimaIDs[merged_it->second].second->toString().c_str());
            }
          }
          if (finalStructureMinimum != NULL){
            auto final_min_it = sorted_min_and_output_ids->find(*finalStructureMinimum);
            if (final_min_it != sorted_min_and_output_ids->end()){
              auto merged_it = merged_min_to_representative.find(final_min_it->second);
              if (merged_it != merged_min_to_representative.end()){
                endStructureMinimum = strcpy(endStructureMinimum, sortedMinimaIDs[merged_it->second].second->toString().c_str());
              }
            }
          }

          // merge partition functions
          for(size_t i = 0; i < sortedMinimaIDs.size(); i++){
            size_t id_i = sortedMinimaIDs[i].first;
            for(size_t j = 0; j < sortedMinimaIDs.size(); j++){
              size_t id_j = sortedMinimaIDs[j].first;
              auto trans = z.find(SC_PartitionFunction::PairID(id_i, id_j));
              if (trans != z.end()){
                size_t cluster_i = id_i;
                size_t cluster_j = id_j;
                double pf_from =  trans->second.getZ();
                auto i_rep = merged_min_to_representative.find(i);
                auto j_rep = merged_min_to_representative.find(j);
                if (i_rep != merged_min_to_representative.end()){
                  cluster_i = sortedMinimaIDs[i_rep->second].first;
                }
                if (j_rep != merged_min_to_representative.end()){
                  cluster_j = sortedMinimaIDs[j_rep->second].first;
                }
                if ((cluster_i == cluster_j) && (j != i))
                   continue; // don't add inner cluster transitions
                auto trans_c = z_minh.find(SC_PartitionFunction::PairID(cluster_i, cluster_j));
                if (trans_c != z_minh.end()){
                  double pf_to = trans_c->second.getZ();
                  trans_c->second.setZ(pf_to + pf_from);
                }
                else{
                  z_minh[SC_PartitionFunction::PairID(cluster_i, cluster_j)] = trans->second;
                  //z_minh[SC_PartitionFunction::PairID(cluster_i, cluster_j)].setZ(pf_from);
                }
              }
            }
          }
          // replace z with new pf.
          z.clear();
          z = z_minh;

          // merge dotplots
          if(writeDotplotPerBasin){
            for(auto it = representatives_and_clustered_ids.begin(); it != representatives_and_clustered_ids.end(); it++){
              size_t representative = it->first;
              MyState *rep_minimum = sortedMinimaIDs[representative].second;
              // TODO: sum up all other properties of SC_PartitionFunction (for example the Energies) this is currently not supported for minh-clustering.
              SC_DotPlot::DotPlot* dp_rep = &dot_plot_per_basin[*rep_minimum];
              for(auto it_cluster = it->second.begin(); it_cluster != it->second.end(); it_cluster++){
                if (it->first != *it_cluster){
                  // add all clustered dotplots to the representative.
                  MyState *state = sortedMinimaIDs[*it_cluster].second;
                  //*dp_rep += dot_plot_per_basin[*state];
                  SC_DotPlot::DotPlot* dp_state = &dot_plot_per_basin[*state];
                  for(auto it_dp = dp_state->begin(); it_dp != dp_state->end(); it_dp++){
                    (*dp_rep)[it_dp->first] += it_dp->second;
                  }
                }
              }
            }
          }

          // merge saddles
          if (args_info.saddle_file_arg){
            for(auto it = representatives_and_clustered_ids.begin(); it != representatives_and_clustered_ids.end(); it++){
              size_t representative = it->first;
              MyState *rep_minimum = sortedMinimaIDs[representative].second;
              for(auto it_cluster = it->second.begin(); it_cluster != it->second.end(); it_cluster++){
                // merge all saddles to the representative.
                if (args_info.saddle_file_arg){
                  MyState *state = sortedMinimaIDs[*it_cluster].second;
                  auto it_saddles = all_saddles.find(*state);
                  // insert minimal saddle to neighbor which is not in the same cluster
                  for(auto it_saddle = it_saddles->second.begin(); it_saddle != it_saddles->second.end(); it_saddle++){
                    auto it_output_id = sorted_min_and_output_ids->find(it_saddle->first);
                    if (it_output_id != sorted_min_and_output_ids->end()){
                      size_t output_id = it_output_id->second;
                      if (it->second.find(output_id) == it->second.end() && output_id != representative){
                        // if not in the same cluster --> add transition
                        auto rep_neighbor_it = all_saddles[*rep_minimum].find(it_saddle->first);
                        if (rep_neighbor_it != all_saddles[*rep_minimum].end()){
                          // replace saddle if smaller
                          if(it_saddle->second.energy < rep_neighbor_it->second.energy){
                            all_saddles[*rep_minimum][it_saddle->first] = it_saddle->second;
                            all_saddles[it_saddle->first][*rep_minimum] = it_saddle->second;
                          }
                        }
                        else{
                          // add new saddle
                          all_saddles[*rep_minimum][it_saddle->first] = it_saddle->second;
                          all_saddles[it_saddle->first][*rep_minimum] = it_saddle->second;
                        }
                      }
                    }
                  }
                  // remove merged saddles
                  all_saddles.erase(*state);
                }
              }
            }
          }

          if (basin_size){
            // merge basin sizes.
            for(auto it = representatives_and_clustered_ids.begin(); it != representatives_and_clustered_ids.end(); it++){
              size_t basin_rep_global_id = sortedMinimaIDs[it->first].first;
              for(auto itc = it->second.begin(); itc != it->second.end(); itc++){
                if (*itc != it->first)
                  minimum_index_and_basin_size[basin_rep_global_id] += minimum_index_and_basin_size[sortedMinimaIDs[*itc].first];
              }
            }
          }

          //store sorted mapped minima in list
          MyState *min, *rep;
          for(size_t i = 0; i < sortedMinimaIDs.size(); i++){
            auto itrep = representatives_and_clustered_ids.find(i);
            std::vector<MyState*> basins;
            if(itrep != representatives_and_clustered_ids.end()){
              rep = sortedMinimaIDs[itrep->first].second;
              for(auto it = itrep->second.begin(); it != itrep->second.end(); it++){
                min = sortedMinimaIDs[*it].second;
                basins.push_back(min);
              }
              std::sort(basins.begin(), basins.end(), [](const MyState *a, const MyState *b) -> bool {return *a < *b;});
            }
            else{
              rep = sortedMinimaIDs[i].second;
              basins.push_back(rep);
            }
            minh_representatives_and_basins.push_back(std::pair<MyState*,std::vector<MyState*>>(rep,basins));
          }

          // then create new sorted_min list and hash map for output --> simply filter the original list and assign new IDs.
          for(size_t min_index = sortedMinimaIDs.size(); min_index > 0; min_index--){
            if(merge_state_map.find(min_index-1) != merge_state_map.end()){
              sorted_min_and_output_ids->erase(*(sortedMinimaIDs[min_index-1].second));
              sortedMinimaIDs.erase(sortedMinimaIDs.begin() + min_index-1);
            }
          }
          // assign new ids (still sorted).
          std::unordered_map<size_t, size_t> map_ids;
          for(size_t i = 0; i < sortedMinimaIDs.size(); i++){
            map_ids[(*sorted_min_and_output_ids)[*(sortedMinimaIDs[i].second)]] = i;
            (*sorted_min_and_output_ids)[*(sortedMinimaIDs[i].second)] = i;
          }
          // adjust saddle list to new ids.
          for(auto it = minimal_saddle_list.begin(); it != minimal_saddle_list.end(); it++){
            it->minimum_from = map_ids[it->minimum_from];
            it->minimum_to = map_ids[it->minimum_to];
          }

          if (basin_size){
            // filter and adjust ids.
            std::unordered_map<size_t, size_t> minimum_index_and_basin_size_filtered;
            for(auto it = sortedMinimaIDs.begin(); it != sortedMinimaIDs.end(); it++){
              minimum_index_and_basin_size_filtered[it->first] = minimum_index_and_basin_size[it->first];
            }
            minimum_index_and_basin_size.clear();
            minimum_index_and_basin_size.insert(minimum_index_and_basin_size_filtered.begin(), minimum_index_and_basin_size_filtered.end());
          }
        }
      }

      std::cout << std::endl;
      std::string out               = "Sequence: ";
      printf("%s", out.c_str());
      int         offset          = 31;
      int         padding_length  = offset - out.length(); //structure_length - out.length();
      std::string pad             = std::string(padding_length, ' ');
      printf("%s%s\n", pad.c_str(), rnaSequence.c_str());

      out = "MFE structure: ";
      printf("%s", out.c_str());
      padding_length  = offset - out.length();
      pad             = std::string(padding_length, ' ');
      printf("%s%s\n", pad.c_str(), mfeStructure);

      out = "The start state is: ";
      printf("%s", out.c_str());
      padding_length  = offset - out.length();
      pad             = std::string(padding_length, ' ');
      printf("%s%s\n", pad.c_str(), rna_start_str.c_str());

      out = "The start state ends in basin: ";
      printf("%s", out.c_str());
      padding_length  = offset - out.length();
      pad             = std::string(padding_length, ' ');
      printf("%s%s\n", pad.c_str(), startStructureMinimum);
      out = "The final state is: ";
      printf("%s", out.c_str());
      padding_length  = offset - out.length();
      pad             = std::string(padding_length, ' ');
      printf("%s%s\n", pad.c_str(), rna_final_str.c_str());
      out = "The final minimum is: ";
      printf("%s", out.c_str());
      if (endStructureMinimum != NULL) {
        padding_length  = offset - out.length();
        pad             = std::string(padding_length, ' ');
        printf("%s%s\n", pad.c_str(), endStructureMinimum);
      } else {
        printf("\n");
      }

      printf("Number of minima: %ld\n", done_List.size());

      std::cout << std::endl;

      // Calculate the final Rate Matrix:
      *transOut
      << "              --------------THE FINAL RATE MATRIX----------------------- "
      << std::endl;

      biu::MatrixSparseC<double>&  final_Rate = calculateRateMatrix(z, sortedMinimaIDs, computeDiagonal);

      // Print the Final rate matrix of the States in Final minima set.
      // printRateMatrix (*final_Rate, final_minima, *transOut, true);
      if(writeDotplotPerBasin){
        size_t min_id = 1;
        double pf_basin;
        std::vector<std::string> mea_structures_per_basin;
        for(auto it = sortedMinimaIDs.begin(); it != sortedMinimaIDs.end(); it++, min_id++){
          size_t minIndex = it->first;
          SC_PartitionFunction* sc = &z.at(SC_PartitionFunction::PairID(minIndex, minIndex));
          pf_basin = z.at(SC_PartitionFunction::PairID(minIndex, minIndex)).getZ();
          SC_DotPlot::DotPlot dp_tmp = dot_plot_per_basin[*it->second];
          SC_DotPlot::DotPlot dp = SC_DotPlot::getBasePairProbabilities(dp_tmp, pf_basin);

          char * dp_file_name = (char *)dotPlotPerBasinFileName.c_str();
          std::string dp_file_pattern = dotPlotPerBasinFileName +  "_%u.ps";
          int res = asprintf(&dp_file_name, dp_file_pattern.c_str(), min_id);
          if(res < 0) fprintf(stderr, "Error: failed to concatenate dotplot name!");
          std::string basin_representative = it->second->toString();
          char * mea_structure = SC_DotPlot::mea_from_dotplot(sequence, dp, vc->exp_params);
          std::string mea(mea_structure);
          mea_structures_per_basin.push_back(mea);
          free(mea_structure);
          bool sane = SC_DotPlot::writeDotPlot_PS_with_mfe_and_mea(dp_file_name, rnaSequence, dp, basin_representative, mea);
          if (!sane){
            fprintf(stderr,"Error: could not write the dot plot! %s\n", dp_file_name);
            fprintf(stderr, "%s %s %s %s\n", dp_file_name, rnaSequence.c_str(), basin_representative.c_str(), mea.c_str());
          }
          free(dp_file_name);
        }

        printRateMatrixSortedWithMEA(final_Rate, sortedMinimaIDs, mea_structures_per_basin, *transOut);
      }
      else{
        printRateMatrixSorted(final_Rate, sortedMinimaIDs, *transOut);
      }
      std::cout << std::endl;
      printEquilibriumDensities(z, sortedMinimaIDs, *transOut);
      std::cout << std::endl;

      //if parameter given: print Z-matrix to file.
      if ((partitionFunctionFileName.size() > 0) & writePartitionFunctions) {
        ofstream ptFile;
        ptFile.open(partitionFunctionFileName);
        ptFile
        << "              --------------THE PARTITIONFUNCTION MATRIX----------------------- "
        << std::endl << std::endl;
        printZMatrixSorted(z, sortedMinimaIDs, Minima, ptFile);
        ptFile.close();
      }

      double partitionFunctionLandscape = 0;
      //get partition Sum for the energy landscape.
      for (auto it = sortedMinimaIDs.begin(); it != sortedMinimaIDs.end(); it++) {
        size_t  minIndex  = it->first;
        partitionFunctionLandscape +=
          z.at(SC_PartitionFunction::PairID(minIndex, minIndex)).get_unscaled_Z();
      }
      std::cout << "The overall partition function is: "
                << partitionFunctionLandscape << std::endl;

      if(args_info.minh_given || args_info.dynamic_minh_given){
        // write mapping file
        std::ofstream mapping_file_minh(minh_mapping_file);
        const size_t LEAD = 6;
        for(size_t i = 0; i < minh_representatives_and_basins.size(); i++){
          mapping_file_minh << std::setw(LEAD) << std::to_string(i) << " [" << minh_representatives_and_basins[i].first->toString() << "]:" << std::endl;
          std::vector<MyState*>& basin_vec = minh_representatives_and_basins[i].second;
          for(size_t j = 0; j < basin_vec.size(); j++){
            mapping_file_minh << basin_vec[j]->toString() << std::endl;
          }
        }
      }
      //dotplot
      if (writeDotplot) {
        SC_DotPlot::DotPlot normalizedDotplot =
          SC_DotPlot::getBasePairProbabilities(dotplot,
                                               partitionFunctionLandscape);
        std::string mfe(mfeStructure);
        char * mea_structure = SC_DotPlot::mea_from_dotplot(sequence, normalizedDotplot, vc->exp_params);
        std::string mea(mea_structure);
        free(mea_structure);
        bool                dotplotWritten = SC_DotPlot::writeDotPlot_PS_with_mfe_and_mea(dotPlotFileName,
                                                                                          rnaSequence,
                                                                                          normalizedDotplot, mfe, mea);
        if (!dotplotWritten) {
          std::cerr
          << "\nWarning: error during the writing of the dot plot file "
          << dotPlotFileName << "\n";
        }
      }

      print_number_of_rates(final_Rate, sortedMinimaIDs, std::cout);

      if (basin_size){
        std::cout << std::endl;
        std::cout << "Basin size:" << std::endl;
        printf("(");
        size_t n = sortedMinimaIDs.size()-1;
        for (size_t i = 0; i <= n; i++){
          if (i < n)
            fprintf(stdout, "%ld, ", minimum_index_and_basin_size[sortedMinimaIDs[i].first]);
          else
            fprintf(stdout, "%ld", minimum_index_and_basin_size[sortedMinimaIDs[i].first]);
        }
        printf(")\n");
      }


      if (!binary_rates_file.empty())
        write_binary_rates_file(binary_rates_file, final_Rate, sortedMinimaIDs);

      /**
       * write saddle file
       * */
      if(args_info.saddle_file_given){
        FILE* saddle_file = fopen(saddleFileName.c_str(), "w");
        if (!saddle_file) {
            fprintf(stderr, "Error: could not open saddle file.\n");
            exit(101);
          }
        std::string header = std::string("id_from, loc_min_from, loc_min_from_energy, id_to, loc_min_to, ")+
            std::string("loc_min_to_energy, saddle, saddle_energy\n");
        fwrite(header.c_str(), sizeof(char), header.length(), saddle_file);

        for (auto it = all_saddles.begin(); it != all_saddles.end(); it++) {
          //const std::pair<MyState, MyState>&  from_to       = it->first;
          const MyState *state_from = &it->first;
          auto id_from_it = sorted_min_and_output_ids->find(*state_from);
          if (id_from_it != sorted_min_and_output_ids->end()){
            for (auto it_to = it->second.begin(); it_to != it->second.end(); it_to++) {
              const MyState *state_to = &it_to->first;
              auto id_to_it   = sorted_min_and_output_ids->find(*state_to);
              if ( id_to_it != sorted_min_and_output_ids->end()){
                const MyState *state_saddle = &it_to->second;
                double                              saddle_height = state_saddle->energy / 100.0;
                std::string                         s1            = state_from->toString();
                std::string                         s2            = state_to->toString();
                std::string                         saddle        = state_saddle->toString();
                  int id_from = id_from_it->second;
                  int id_to = id_to_it->second;
                  std::fprintf(saddle_file,"%d, %s, %.2f, %d, %s, %.2f, %s, %.2f\n", id_from,
                                        s1.c_str(), state_from->energy / 100.0,
                                        id_to, s2.c_str(), state_to->energy / 100.0, saddle.c_str(), saddle_height);
              }
            }
          }
        }
        fclose(saddle_file);
      }

      if(args_info.barriers_like_output_given){
        write_barriers_like_output(barriers_prefix,
                                    final_Rate,
                                    sortedMinimaIDs,
                                    rnaSequence,
                                    basin_size,
                                    minimum_index_and_basin_size);
      }

      if(args_info.binary_rates_file_sparse_given){
        std::string sparse_matrix_file(args_info.binary_rates_file_sparse_arg);
        write_binary_rates_file_sparse(sparse_matrix_file,
                                      final_Rate,
                                      sortedMinimaIDs);
      }

      // write barriers tree
      if (args_info.barrier_tree_file_given){
        std::unordered_map<size_t, int> structure_index_to_energy;
        for(auto it = sorted_min_and_output_ids->begin(); it != sorted_min_and_output_ids->end(); it++){
          structure_index_to_energy[it->second] = it->first.energy;
        }

        size_t output_mfe_id = sortedMinimaIDs[0].first;
        int output_mfe = sortedMinimaIDs[0].second->energy;
        structure_index_to_energy[output_mfe_id] = output_mfe;

        std::vector<node_t*> forest = bt.create_barrier_tree(minimal_saddle_list, structure_index_to_energy);
        //node_t* tree = forest[0]; //TODO: search mfe tree.
        std::cout << std::endl;
        //bt.free_tree(tree);
        for(size_t t = 0; t < forest.size(); t++){
          std::string newick_tree = bt.newick_string_builder(forest[t]);
          std::ofstream tree_file(barrierTreeFileName);
          tree_file << newick_tree << std::endl;
          tree_file.close();
          //std::cout << newick_tree << std::endl;
        }
        if(forest.size() > 0){
          size_t inorder_index = 0;
          std::string svg_tree = bt.svg_string_builder(forest[0], output_mfe, &inorder_index);
          std::ofstream tree_file(barrierTreeFileName + "_tree.svg");
          tree_file << svg_tree << std::endl;
          tree_file.close();
        }
        for(size_t t = 0; t < forest.size(); t++){
          bt.free_tree(forest[t]);
        }
      }

      delete &final_Rate;

    } while(dynamic_minh_max_states.size() > 0);

    if(writeDotplotPerBasin){
      size_t min_id = 1;
      double pf_basin;
      for(auto it = sortedMinimaIDs.begin(); it != sortedMinimaIDs.end(); it++, min_id++){
        size_t minIndex = it->first;
        pf_basin = z.at(SC_PartitionFunction::PairID(minIndex, minIndex)).getZ();
        SC_DotPlot::DotPlot dp_tmp = dot_plot_per_basin[*it->second];
        SC_DotPlot::DotPlot dp = SC_DotPlot::getBasePairProbabilities(dp_tmp, pf_basin);

        char * dp_file_name = (char *)dotPlotPerBasinFileName.c_str();
        std::string dp_file_pattern = dotPlotPerBasinFileName +  "_%u.ps";
        int res = asprintf(&dp_file_name, dp_file_pattern.c_str(), min_id);
        if(res < 0) fprintf(stderr, "Error: failed to concatenate dotplot name!");
        std::string basin_representative = it->second->toString();
        char * mea_structure = SC_DotPlot::mea_from_dotplot(sequence, dp, vc->exp_params);
        std::string mea(mea_structure);
        free(mea_structure);
        bool sane = SC_DotPlot::writeDotPlot_PS_with_mfe_and_mea(dp_file_name, rnaSequence, dp, basin_representative, mea);
        if (!sane)
          fprintf(stderr,"Error: could not write the dot plot!\n");
        free(dp_file_name);
      }
    }

    /***********Garbage collection ************/
    if (startStructureMinimum != NULL)
      free(startStructureMinimum);

    if (rnaFinalStructurePT != NULL)
      free(rnaFinalStructurePT);

    if (finalStructureMinimum != NULL)
      delete  finalStructureMinimum;

    free(endStructureMinimum);
    delete sorted_min_and_output_ids;
    //	delete ScMinimum;
    vrna_fold_compound_free(vc);
    if (mfeStructure != NULL)
      free(mfeStructure);

    neighborFilter = NULL;
    if (neighborCombinator != NULL)
      delete neighborCombinator;

    if (neighbor_K_Filter != NULL)
      delete neighbor_K_Filter;

    if (neighbor_E_Filter != NULL)
      delete neighbor_E_Filter;

    free(rnaStructurePT);
    Minima.clear();
    MinimaForReverseSearch.clear();
    //final_minima.clear();
    z.clear();
    if (startStateMinimum != NULL)
      delete startStateMinimum;

    if (tmp_sorted_min_and_output_ids != NULL)
      delete tmp_sorted_min_and_output_ids;

    pourRNA_cmdline_parser_free(&args_info);

  } catch (std::exception & ex) {
    std::cerr << "\n\n ERORR : " << ex.what() << "\n" << std::endl;
    return_Value = -1;
  }


  //print stopwatch:
  end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::time_t                   end_time        = std::chrono::system_clock::to_time_t(end);
  std::cout << "finished computation at " << std::ctime(&end_time)
            << "elapsed time: " << elapsed_seconds.count() << "s\n";
  return return_Value;
}   // end main body.
