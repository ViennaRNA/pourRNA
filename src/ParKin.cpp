/*
 * ParKin_Explore.cpp
 *
 * Main file for ParKin_Explore 1.0
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
#include <omp.h>
//basic file operations
#include <ostream>
#include <string>
#include <iostream>
#include <fstream>
#include <unistd.h>
extern "C" {
#include <ViennaRNA/move_set.h>
#include <ViennaRNA/fold.h>
//#include <ViennaRNA/read_epars.h>
#include <ViennaRNA/params/io.h>
#include <ViennaRNA/eval.h>
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/utils/structures.h>
#include <ViennaRNA/neighbor.h>
}
#include "BIUlibPart/OptionParser.h"
#include "BIUlibPart/MatrixSparse.hh"
#include "RNAkinetics/RateMatrixUtil.h"
#include "RNAkinetics/RNA_NeighMinFilter.h"
#include "Flooder.h"
#include "TypeID.h"
#include "SC_DotPlot.h"

#include "Concurrent_Queue.h"
#include "PairHashMap.h"

/*! Exception class for exceptions thrown during argument and input parsing.
 */
class ArgException: public std::exception {
public:
//! the error message
	std::string errorMsg;
	ArgException(std::string errorMsg_) :
			errorMsg(errorMsg_) {
	}
	virtual ~ArgException() throw () {
	}
	virtual const char*
	what() const throw () {
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
size_t getMaximalNeighborsOfAMacroState(SC_PartitionFunction::Z_Matrix& z,
		std::unordered_set<size_t>& done_List,
		std::unordered_map<size_t, MyState>& minimaForReverseSearch) {

	size_t maxNeighbors = 0;
	size_t countNeighborsForState = 0;
	for (std::unordered_set<size_t>::iterator from = done_List.begin();
			from != done_List.end(); from++) {
		SC_PartitionFunction::PairID basinID = SC_PartitionFunction::PairID(
				*from, *from);
		if (z.find(basinID) != z.end()) {
			countNeighborsForState = 0;
			bool transitionPartitionsumExists = false;
			for (std::unordered_set<size_t>::iterator to = std::next(from, 1);
					to != done_List.end(); to++) {
				SC_PartitionFunction::PairID transitionID =
						SC_PartitionFunction::PairID(*from, *to);
				SC_PartitionFunction::PairID reverseTransitionID =
						SC_PartitionFunction::PairID(*to, *from);

				double z_transition = 0;
				double z_reverseTransition = 0;
				if (z.find(transitionID) != z.end()) {
					z_transition = z.at(transitionID).getZ();
					if (z_transition > 0.0) {
						transitionPartitionsumExists = true;
					}
				}
				if (z.find(reverseTransitionID) != z.end()) {
					z_reverseTransition = z.at(reverseTransitionID).getZ();
					if (z_reverseTransition > 0.0) {
						transitionPartitionsumExists = true;
					}
				}
				if (transitionPartitionsumExists) {
					countNeighborsForState++;
				}
			}
			if (maxNeighbors < countNeighborsForState) {
				maxNeighbors = countNeighborsForState;
			}
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
 * @param maxNeighNum the maximum number of macrostates is used to normalize the rates.
 * @return the biu::MatrixSparseC<double> RateMatrix.
 */
biu::MatrixSparseC<double> &
calculateRateMatrix(SC_PartitionFunction::Z_Matrix& z,
		std::unordered_set<size_t>& done_List,
		std::unordered_map<size_t, MyState>& minimaForReverseSearch,
		size_t maxNeighNum) {
	biu::MatrixSparseC<double>& final_Rate = *new biu::MatrixSparseC<double>(
			done_List.size(), done_List.size(), 0.0);
	// iterate through all elements of done List to produce the final rate matrix
	std::unordered_set<size_t>::const_iterator toIt = done_List.begin();
	std::unordered_set<size_t>::const_iterator fromIt = done_List.begin();
	for (size_t from = 0; from < done_List.size(); from++, fromIt++) {
		// copy according rates
    size_t fromOrig = *fromIt;
		toIt = done_List.begin();
		double rateSum = 0.0;
		for (size_t to = 0; to < done_List.size(); to++, toIt++) {
			size_t toOrig = *toIt;
			// copy non-diagonal entries
			if (toOrig != fromOrig) {
				//compute the rate from row-state two column-state
				//in Z-Matrix: from i to j.
				//in final_Rate: from j to i !
				SC_PartitionFunction::PairID transitionID =
						SC_PartitionFunction::PairID(fromOrig, toOrig);
				SC_PartitionFunction::PairID reverseTransitionID =
						SC_PartitionFunction::PairID(toOrig, fromOrig);
				SC_PartitionFunction::PairID basinID =
						SC_PartitionFunction::PairID(fromOrig, fromOrig);

				double z_transition = 0;
				double z_reverseTransition = 0;
				if (z.find(transitionID) != z.end()) {
					z_transition = z.at(transitionID).getZ();
				}
				if (z.find(reverseTransitionID) != z.end()) {
					z_reverseTransition = z.at(reverseTransitionID).getZ();
				}
				//max is important if filters are applied! Otherwise it should be equal.
				z_transition = max(z_transition, z_reverseTransition);

				double z_basin = 0;
				if (z.find(basinID) != z.end()) {
					z_basin = z.at(basinID).getZ();
				}
				if (z_basin > 0.0) {
					double rate_Val = z_transition
							/ (/*maxNeighNum * */z_basin);
					if (rate_Val > 0.0) {
						final_Rate.at(to, from) = rate_Val;
						rateSum += rate_Val;
					}
				}
			}
		}
		final_Rate.at(from, from) = 1.0 - rateSum;
	}		// end for
	return final_Rate;
}

/* !AllowedArgument function to check the given input arguments
 * @param allowed, type OptionMap, which is list of possible options for the parameter list.
 *  allowed parameter is container of given input, any given input will be pushed back
 *  to allowed param, each input of type biu::COption will contain parameters indicated type
 *  of the input value. e.g "seq" as Sequence or "str" as Structure etc. and if it is optional
 *  or not .In addition description of each input of type biu::COption. all the input parameters
 *  of type biu::COption that allowed param contains, will be extracted during Parsing the parameters
 *  phase using biu::COption_Parser object.
 *
 * @param infos, contains some description about the possible options.
 */
void
writeDescription(biu::OptionMap & allowed, std::string & infos);

struct flooderInputParameter {
	size_t MaxToQueue;
	size_t MaxToHash;
	double MaxEnergy;
	double DeltaE;
	size_t BasinID;
	MyState* CurrentMinimum;
	NeighMinFilter* Filter;
	Concurrent_Queue<MyState> *DiscoveredMinima;
	double TemperatureForBoltzmannWeight;
	unsigned int Move_set;
	PairHashMap::HashMap* All_Saddles;
	int MaxBPdist;
	std::string SourceStructure;
	std::string TargetStructure;
  ~flooderInputParameter() {
    if (CurrentMinimum != NULL)
      delete CurrentMinimum;
  }
};

struct flooderOutputAnylse {
	size_t NumberOfStatesInBasin;
	size_t NumberOfStatesOutOfBasin; //outer contact surface.
	size_t TimeForStatesInBasin;
	size_t TimeForStatesOutOfBasin;
	double MaxEnergyThreshold;
	size_t ProcessedStates;
};

struct flooderOutputParameter {
	// ! All neighbored minima.
	PairHashTable::HashTable LocalThreadMinima;
	// ! Only the indices are filtered.
	std::vector<size_t> IndicesOfFilteredMinima;
	SC_PartitionFunction::Z_Matrix PartitionFunctions;
	SC_PartitionFunction* ScBasin;
	flooderOutputAnylse Analysis;
	~flooderOutputParameter() {
		if (ScBasin != NULL)
			delete ScBasin;
	}
};

/**
 * ! flood the basin with the given structure "flooderInputParameter->CurrentMinimum".
 *   The transition-partitionsums and the partitionsum of the basin will be written to
 *   "flooderOutputParameter->PartitionFunctions".
 *   @param inParameter contains a structure which identifies the basin and filter parameters.
 *   @param outParameter contains the partitionfunctions and neigbored minima.
 */
int floodBasin(vrna_fold_compound_t *vc, flooderInputParameter* inParameter,
		flooderOutputParameter* outParameter) {
	size_t loopCurrentMinID = inParameter->BasinID;

	PairHashTable::HashTable& localThreadMinima =
			outParameter->LocalThreadMinima;
	// initializing set, which contains the NEighbors indices  of current state
	size_t loopCurrentLocalMinID = 0;

	MyState* minState = new MyState(*inParameter->CurrentMinimum); //TODO: create new instance! (cleaned if localThreadMinima is cleaned)
	//set local id to zero because of the statePairCollector id-assign-method.

	localThreadMinima.insert( { *minState, loopCurrentLocalMinID });

	// start walking and computing the Partition Function

	// creating the StatePaircollector to manage the state of Pairs between the neighbor basins
	StatePairCollector spc(loopCurrentLocalMinID, localThreadMinima,
			outParameter->PartitionFunctions, inParameter->MaxToHash,
			inParameter->DiscoveredMinima,
			inParameter->TemperatureForBoltzmannWeight,
			inParameter->Move_set,
			*inParameter->All_Saddles,
			inParameter->SourceStructure,
			inParameter->TargetStructure,
			inParameter->MaxBPdist
	);

	// set the maximal energy as upper floodlevel.
	double maxEnergy = min(inParameter->MaxEnergy,
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
				{
        outParameter->IndicesOfFilteredMinima.push_back(it->second);
		}
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
	outParameter->Analysis.TimeForStatesInBasin = 0;
	outParameter->Analysis.TimeForStatesOutOfBasin = 0;
	outParameter->Analysis.MaxEnergyThreshold =
			myFlooder.getFinalEnergyThreshold();
	outParameter->Analysis.ProcessedStates = myFlooder.getProcessedStates();

	return 0;
}

//######################################################## -- MAIN BODY -- ########################################################

/*
 * return the maximal available physical memory.
 */
size_t getTotalSystemMemoryInBytes() {
	long pages = sysconf(_SC_PHYS_PAGES);
	long page_size = sysconf(_SC_PAGE_SIZE);
	return pages * page_size;
}

/*
 * The hashmap for the gradientwalk consumes the most memory.
 * To prevent the system from swapping (and in the end from crashing), the maximal size of the hashmap should be estimated.
 */
size_t estimateMaxToHash(size_t numberOfHashMaps, size_t sequenceLength) {
	//determine estimated maximum maxToHash value without swapping.
	size_t memInBytes = getTotalSystemMemoryInBytes();
	size_t sizeOfHashEntry = sizeof(size_t) + sizeof(MyState)
			+ sizeof(short) * (sequenceLength + 1);
	size_t numberOfEntries = memInBytes / sizeOfHashEntry;
	size_t entriesPerHashMap = numberOfEntries / numberOfHashMaps;
	return entriesPerHashMap;
}

int merge_results(
		std::vector<std::pair<flooderInputParameter*, flooderOutputParameter*>> *threadParameter,
		bool logEnergies, ofstream &energyFile,
		std::unordered_set<int> &addedMinIDs, bool writeDotplot,
		SC_DotPlot::DotPlot &dotplot, PairHashTable::HashTable &Minima,
		std::unordered_map<size_t, MyState> &MinimaForReverseSearch,
		SC_PartitionFunction::Z_Matrix &z, bool dynamicBestK,
		std::unordered_map<size_t, std::vector<size_t>> &dynamicBestKFilterNeighborList,
		std::list<size_t> &toDo_List, std::unordered_set<size_t> &done_List, bool verbose) {
	for (int k = 0; k < threadParameter->size(); ++k) {
		std::pair<flooderInputParameter*, flooderOutputParameter*>* bothParameter;
		bothParameter = &threadParameter->at(k);

		flooderInputParameter* inParameter = bothParameter->first;
		flooderOutputParameter* outParameter = bothParameter->second;

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
				int energy;
				for (int m = 0; m < energies.size(); m++) {
					energy = energies[m];
					energyFile << energy << " ";
				}
			}
		}

		//sum up dotplot basepair partition sums.
		if (writeDotplot) {
			SC_DotPlot::DotPlot tmpDotPlot =
					((SC_DotPlot*) outParameter->ScBasin)->getBasePairWeightSum();
			for (SC_DotPlot::DotPlot::iterator it = tmpDotPlot.begin();
					it != tmpDotPlot.end(); ++it) {
				SC_DotPlot::DotPlot::iterator entryIt = dotplot.find(it->first);
				if (entryIt != dotplot.end()) {
					dotplot[it->first] += it->second;
				} else {
					dotplot[it->first] = it->second;
				}
			}
		}
		PairHashTable::HashTable& localThreadMinima =
				outParameter->LocalThreadMinima;
		size_t loopCurrentMinID = inParameter->BasinID;

		//add localMinima to total minima and adjust indices in localPartitionfunctions
		//and in indicesOfNeighbors if needed.

		std::unordered_map<size_t, size_t> mapOldIndexToNewIndex;

		for (auto minIt = localThreadMinima.begin();
				minIt != localThreadMinima.end(); minIt++) {
			size_t oldIndex = minIt->second;

			const MyState& localThreadMin = minIt->first;
			if (Minima.find(localThreadMin) == Minima.end()) { //if local min is not in total list.
															   //add new min and set lowest maximal index.
				size_t lowestMaxIndex = TypeID::value<size_t>();

				Minima.insert( { MyState(localThreadMin), lowestMaxIndex });
				MinimaForReverseSearch.insert( { lowestMaxIndex,
				  MyState(localThreadMin) });
			}

			size_t newIndex = Minima.at(localThreadMin);

			//assingn total id to localID.
			mapOldIndexToNewIndex.insert( { oldIndex, newIndex });
		}

		SC_PartitionFunction::Z_Matrix& localZ =
				outParameter->PartitionFunctions;
		//adjust index in partitionfunction z.
		for (auto it = localZ.begin(); it != localZ.end(); it++) {
			if (mapOldIndexToNewIndex.find(it->first.first)
					!= mapOldIndexToNewIndex.end()
					&& mapOldIndexToNewIndex.find(it->first.second)
							!= mapOldIndexToNewIndex.end()) {
				size_t newFirstIndex = mapOldIndexToNewIndex.at(
						it->first.first);
				size_t newSecondIndex = mapOldIndexToNewIndex.at(
						it->first.second);
				SC_PartitionFunction::PairID localPairID =
						SC_PartitionFunction::PairID(newFirstIndex,
								newSecondIndex);
				z.insert( { localPairID, it->second });
			}
		}

		//fill dynamic Best-k filter list.
		if (dynamicBestK) {
			if (dynamicBestKFilterNeighborList.find(loopCurrentMinID)
					== dynamicBestKFilterNeighborList.end()) {
				vector<size_t> allNeighboredMinima;
				size_t currentMinID = localThreadMinima.at(
						*inParameter->CurrentMinimum);
				for (auto localMinIt = localThreadMinima.begin();
						localMinIt != localThreadMinima.end(); ++localMinIt) {
					if (localMinIt->second != currentMinID) {
						allNeighboredMinima.push_back(localMinIt->second);
					}
				}

				// creating an greaterByZ Object as third parameter of std:: sort method
				GreaterByZMatrix greaterByZ(localZ, currentMinID);

				// sort bz Z value
				std::sort(allNeighboredMinima.begin(),
						allNeighboredMinima.end(), greaterByZ);

				for (int nm = 0; nm < allNeighboredMinima.size(); ++nm) {
					allNeighboredMinima[nm] =
							mapOldIndexToNewIndex[allNeighboredMinima.at(nm)];
				}

				dynamicBestKFilterNeighborList.insert( { loopCurrentMinID,
						allNeighboredMinima });
			}
		}

		//  add all neighbored minima not part of DONE list to toDO_list
		for (int i = 0; i < outParameter->IndicesOfFilteredMinima.size(); i++) {
			size_t index = mapOldIndexToNewIndex.at(
					outParameter->IndicesOfFilteredMinima.at(i));
			// if the element at the index position is not already in done List, add it to toDo_list
			if (done_List.find(index) == done_List.end()
					&& std::find(toDo_List.begin(), toDo_List.end(), index)
							== toDo_List.end()) {
				toDo_List.push_back(index);
			}
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

int threadFunction(vrna_fold_compound_t *vc,
		std::vector<std::pair<flooderInputParameter*, flooderOutputParameter*>> *threadParameter) {
	for (int k = 0; k < threadParameter->size(); ++k) {
		flooderInputParameter* inParameter = threadParameter->at(k).first;
		flooderOutputParameter* outParameter = threadParameter->at(k).second;
		floodBasin(vc, inParameter, outParameter);
	}
	return 0;
}

int main(int argc, char** argv) {

  PairHashMap::HashMap all_saddles;

	//init stopwatch:
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	int return_Value = 0;
//fill the data

// output stream to display invalidity of some given input
	std::ostringstream error_Mas;
// output container
	std::ostream* out = &std::cout;
// file to store the transition matrix.
	std::ofstream* transOutFile = NULL;
// output container for the transition matrix.
	std::ostream* transOut = &std::cout;
// the maximal number of elements the underlying queue of the Flooder is allowed to contain.
	size_t maxToQueue = UINT64_MAX;
// maximum number of states to be hashed for each gradient walk.
	size_t maxToHash = UINT64_MAX;
// the maximum energy that a state is allowed to have to be considered by the flooder.
	double maxEnergy = 9999;
// the maximum energy difference that states in a basin can have w.r.t. the local minimum.
	double deltaE = 9999;
// value for best K filter. -1 = do not filter.
	size_t filterValueK = 0;
	bool enableBestKFilter = false;
// increases K if the MFE structure is not in the done list and start the exploration again.
	bool dynamicBestK = false;

// value for min. E filter. -1 = do not filter.
	double filterValueE = 0;
	bool enableDeltaMinEFilter = false;

	//max base pair distance filter
	int maxBPdist = 65536;

	// file to store all energies.
	std::string energyFileName = "";
	// binary rates file
	std::string binary_rates_file = "";
	// this bool tells the flooder if it should store energies. (not recommended for large sequences!)
	bool logEnergies = false;
	/* parameter for writing a postscript-file with a DotPlot
	 (it contains the base pair probabilities for all structures in the (filtered) energylandscape).*/
	bool writeDotplot = false;
  //parameter for writing a postscript-file with a DotPlot.
  std::string dotPlotFileName = "";
	// parameter for writing a postscript dotplot per basin.
	bool writeDotplotPerBasin = false;
	// file name prefix for dotplots per basin.
	//TODO: implement this!
	std::string dotPlotPerBasinFileName = "";

	// file to store all partition functions.
	std::string partitionFunctionFileName = "";
	// write partitionfunctions to a file if the name is given as parameter.
	bool writePartitionFunctions = false;

// parameter for temperature in celsius.
	double temperatureForEnergy = 37;
// paramter for the temperature for the Boltzmann weight in celsius.
	double temperatureForBoltzmannWeight = 37;

	// How to treat \"dangling end\" energies for bases adjacent to helices in free ends and multi-loops
	int danglingEnd = 2;
	// The move set type (see neighbor.h vrna_package)
	unsigned int move_set = 0;
	bool verbose = false;

// parameter: path to file for energy model that is placed in the share/misc folder of this build.
// (set by th "M" parameter: 0=Turner2004,1=Turner1999,2=Andronescu2007)
	std::string energyModelFile = "";
	/*
	 * estimate the maxToHash parameter. The estimate depends on the available memory and the number of threads.
	 * This estimation may slow down the computation, but it will be faster than with a swapping system.
	 * The main purpose for this parameter is, that it reduces the probability that the system will crash if too much
	 * memory is used with large sequences.
	 */
	bool dynamicMaxToHash = false;

//---------------------------- Parsing the Parameters---------------------------------

// to map the allowed arguments.
	biu::OptionMap allowed_Arg;

// information string
	std::string info = "";

	writeDescription(allowed_Arg, info);

// To parse the Program Arguments

	biu::COptionParser out_Parser = biu::COptionParser(allowed_Arg, argc, argv,
			info);

// check if there is any unparseable argument
	if (out_Parser.noErrors()) {

		// show the help output
		if (out_Parser.getBoolVal("help")) {
			out_Parser.coutUsage();
			return 0;
		}
	} else {
		return -1;
	}
	try {

		// parse sequence
		std::string rnaSequence = out_Parser.getStrVal("seq");
		if (rnaSequence.size() == 0) {
			throw ArgException("No RNA sequence given");

		} else if (!StructureUtils::IsValidSequence(rnaSequence)) {
			error_Mas << "Sequence " << rnaSequence << " is no valid Sequence ";
			throw ArgException(error_Mas.str());
		}
		// parse structure
		std::string rna_start_str = out_Parser.getStrVal("startStr");
		if (rna_start_str.size() == 0) {
			// initialize open chain structure
			rna_start_str = std::string(rnaSequence.size(), '.');
		} else if (!StructureUtils::IsValidStructure(rna_start_str)) {
			error_Mas << " start Structure " << rna_start_str
					<< " is no valid structure";
			throw ArgException(error_Mas.str());
		}

		std::string rna_final_str = "";
		if (out_Parser.argExist("finalStr")) {
			rna_final_str = out_Parser.getStrVal("finalStr");
			if (!StructureUtils::IsValidStructure(rna_final_str)
					|| rna_final_str.size() == 0) {
				error_Mas << " final Structure " << rna_final_str
						<< " is no valid structure";
				throw ArgException(error_Mas.str());
			}
		}

		// container for all partition functions to generate
		// indices are according to minima container
		// Z(i,i) will hold the basins partition function
		// Z(i,j) will hold the partition function of the saddle points between minimum i and j
		SC_PartitionFunction::Z_Matrix z;

		// set of all minima found so far
		PairHashTable::HashTable Minima;
		std::unordered_map<size_t, MyState> MinimaForReverseSearch;

		// set best K filter if required
		if (out_Parser.argExist("filterBestK")) {
			filterValueK = out_Parser.getIntVal("filterBestK");
			enableBestKFilter = true;
		}

		// set dynamic best K filter if required
		if (out_Parser.argExist("dynamicBestK")) {
			dynamicBestK = out_Parser.getBoolVal("dynamicBestK");
		}

		// set Energy filter if required
		if (out_Parser.argExist("filterMinE")) {
			filterValueE = out_Parser.getDoubleVal("filterMinE");
			enableDeltaMinEFilter = true;
		}

		// set maxToStore variable for the states in each gradient basin.
		if (out_Parser.argExist("maxToQueue")) {
			maxToQueue = out_Parser.getIntVal("maxToQueue");
		}

		// set maxToHash variable for the states to be hashed in a gradient walk.
		if (out_Parser.argExist("maxToHash")) {
			maxToHash = out_Parser.getIntVal("maxToHash");
		}

		// set dynamicMaxToHash variable for estimating the maximal number of states to be hashed in a gradient walk,
		// by considering the maximal available physical memory.
		if (out_Parser.argExist("dynamicMaxToHash")) {
			dynamicMaxToHash = out_Parser.getBoolVal("dynamicMaxToHash");
		}

		// set the maximum energy that a state is allowed to have to be considered by the flooder.
		if (out_Parser.argExist("maxEnergy")) {
			maxEnergy = out_Parser.getDoubleVal("maxEnergy");
		}

		// set the maximum energy difference that states in a basin can have w.r.t. the local minimum.
		if (out_Parser.argExist("deltaE")) {
			deltaE = out_Parser.getDoubleVal("deltaE");
		}

		// set maxThreads for parallelized computation.
		int maxThreads = 1;
		if (out_Parser.argExist("maxThreads")) {
			maxThreads = out_Parser.getIntVal("maxThreads");
			if (maxThreads < 1) {
				throw ArgException("you should use at least one thread.");
			}
		}


    if (out_Parser.argExist("maxBPdist")) {
      maxBPdist = out_Parser.getIntVal("maxBPdist");
      if (maxBPdist < 0) {
        throw ArgException("The base pair distance has to be positive!");
      }
    }

		if (out_Parser.argExist("energyFile")) {
			// set output file name for energies.
			if (out_Parser.getStrVal("energyFile").size() != 0) {
				energyFileName = out_Parser.getStrVal("energyFile");
				logEnergies = true;
			}
		}

    if (out_Parser.argExist("binary_rates_file")) {
      // set output file name for energies.
      if (out_Parser.getStrVal("binary_rates_file").size() != 0) {
        binary_rates_file = out_Parser.getStrVal("binary_rates_file");
      }
    }

		//partitionFunctions
		if (out_Parser.argExist("partitionFunctions")) {
			// set output file name for partitionfunctions.
			if (out_Parser.getStrVal("partitionFunctions").size() != 0) {
				partitionFunctionFileName = out_Parser.getStrVal(
						"partitionFunctions");
				writePartitionFunctions = true;
			}
		}

		//dotPlot
		if (out_Parser.argExist("dotPlot")) {
			// set output file name for dotPlot.
			if (out_Parser.getStrVal("dotPlot").size() != 0) {
				dotPlotFileName = out_Parser.getStrVal("dotPlot");
				writeDotplot = true;
			}
		}

		if(out_Parser.argExist("dotPlot_per_basin")){
		  if(out_Parser.getStrVal("dotPlot_per_basin").size() != 0){
		    dotPlotPerBasinFileName = out_Parser.getStrVal("dotPlot_per_basin");
		    writeDotplotPerBasin = true;
		  }
		}

		if (out_Parser.argExist("transProb")) {
			// set output stream
			if (out_Parser.getStrVal("transProb").size() == 0) {
				throw ArgException(
						"no transition probability output file given but argument present");
			} else if (out_Parser.getStrVal("transProb").compare("STDOUT")
					!= 0) {
				if (out_Parser.getStrVal("out").compare(
						out_Parser.getStrVal("transProb")) == 0) {
					transOut = out;
				} else {
					transOutFile = new std::ofstream(
							out_Parser.getStrVal("transProb").c_str(),
							std::ofstream::out);
					if (!transOutFile->is_open()) {
						std::ostringstream oss;
						oss << "cannot open output file '"
								<< out_Parser.getStrVal("transProb") << "'";
						throw ArgException(oss.str());
					}
					transOut = transOutFile;
				}
			}
		}

		if (out_Parser.argExist("d")) {
			// set dangling ends.
			danglingEnd = out_Parser.getIntVal("d");
			if (danglingEnd < 0) {
				std::ostringstream oss;
				oss << "Error: dangling end size must be >= 0 !";
				throw ArgException(oss.str());
			}
		}

		if (out_Parser.argExist("T")) {
			// set temperatureForEnergy.
			temperatureForEnergy = out_Parser.getIntVal("T");
			if (!out_Parser.argExist("B")) {
				temperatureForBoltzmannWeight = temperatureForEnergy;
			}
		}

		if (out_Parser.argExist("B")) {
			// set temperature for the boltzmann weight.
			temperatureForBoltzmannWeight = out_Parser.getIntVal("B");
		}

		int tmp_move_set = 0;
    if (out_Parser.argExist("move_set")) {
      // set the move set integer
      tmp_move_set = out_Parser.getIntVal("move_set");
    }
    // convert the integer to the vrna_type
    switch(tmp_move_set){
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


		if (out_Parser.argExist("M")) {
			// set energyModel.
			int energyModelNumber = out_Parser.getIntVal("M");
			if (energyModelNumber >= 0 || energyModelNumber <= 2) {
				// Get the last position of '/'
				char buf[UINT16_MAX] = "";
				ssize_t res = readlink("/proc/self/exe", buf, UINT16_MAX);
				if(res == -1){
				  fprintf(stderr,"Error: could not read /proc/self/exe");
				}
				std::string aux(buf); //(argv[0]);
				// get '/' or '\\' depending on unix/mac or windows.
#if defined(_WIN32) || defined(WIN32)
				int pos = aux.rfind ('\\');
#else
				int pos = aux.rfind('/');
#endif
				std::string shareFolder = "misc";
#ifdef DEBUG
				shareFolder = "misc";
#else
				shareFolder = "share/ParKin";
#endif

				// Get the path of the executable of this program.
				int progPosition = aux.substr(0, pos).rfind('/');
				std::string path = aux.substr(0, progPosition + 1);

				switch (energyModelNumber) {
				case 0:
					energyModelFile = path + shareFolder
							+ "/rna_turner2004.par";
					//printf("%s \n",energyModelFile.c_str());
					read_parameter_file(energyModelFile.c_str());
					/*TODO: find out why the default setting of viennaRNA 2.3.5
					        produces different results than the results after
					        reading the default file
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

		if (out_Parser.argExist("v")) {
			verbose = out_Parser.getBoolVal("v");
		}

		vrna_md_t md;
		vrna_md_set_default(&md);
		md.circ = 0;
		//md.uniq_ML = 1; /* in case we need M1 arrays */
		md.noLP = 0;
		md.gquad = 0;
		md.noGU = 0;
		md.compute_bpp = 0;
		md.temperature = temperatureForEnergy;
		md.dangles = danglingEnd;
		char * sequence = (char *) rnaSequence.c_str();
		vrna_fold_compound_t *vc = vrna_fold_compound(sequence, &md,
		VRNA_OPTION_MFE); // | VRNA_OPTION_PF);

		// set best K filter if required
		// create an filter Object type NeighMinFilter as first filtering Technique
		NeighMinFilter* neighbor_K_Filter = NULL;

		// create an filter Object type NeighMinFilter as first filtering Technique
		NeighMinFilter* neighbor_E_Filter = NULL;

		// create an Combinator object to apply both filtering Techniques in specific order
		NeighMinFilter* neighborCombinator = NULL;

		// final overall filter to be used
		NeighMinFilter * neighborFilter = NULL;
		if (enableBestKFilter) {
			neighbor_K_Filter = new NeighMin_K_Filter(filterValueK);
			neighborFilter = neighbor_K_Filter;
		}

		// set Energy filter if required
		if (enableDeltaMinEFilter) {
			neighbor_E_Filter = new NeighMin_E_Filter(filterValueE);
			neighborFilter = neighbor_E_Filter;
		}

		// set final Filter to be combination of both Filters
		if (neighbor_K_Filter != NULL && neighbor_E_Filter != NULL) {
			neighborCombinator = new NeighMinCombinator(*neighbor_E_Filter,
					*neighbor_K_Filter);
			neighborFilter = neighborCombinator;
		}

		//do flooding

		// start walking to find gradient basin minimum for start state
		short * rnaStructurePT = vrna_ptable(rna_start_str.c_str());
		//printf("%s\n", rna_start_str.c_str());
		int energy = vrna_eval_structure_pt(vc, rnaStructurePT);
		MyState startState(energy, rnaStructurePT);

		//init stopwatch for initial walk:
		std::chrono::time_point<std::chrono::system_clock> startInitialWalk,
				endInitialWalk;
		if (verbose)
			startInitialWalk = std::chrono::system_clock::now();

		MyState* startStateMinimum = WalkGradientHashed(move_set, maxToHash).walk(vc,
				startState);

		if (verbose) {
			//print stopwatch for initial walk:
			endInitialWalk = std::chrono::system_clock::now();
			std::chrono::duration<double> elapsed_seconds_InitialWalk =
					endInitialWalk - startInitialWalk;
			std::cout << "Initial Walk elapsed time: "
					<< elapsed_seconds_InitialWalk.count() << "s\n";
		}

		// add current minimum to set of minima and store its index
		size_t currentMinID = 0;
		Minima.insert( { *startStateMinimum, currentMinID });
		MinimaForReverseSearch.insert( { currentMinID, *startStateMinimum });

		///// for tests only -> write file with all Energies ///
		std::unordered_set<int> addedMinIDs;
		ofstream energyFile;
		if (energyFileName.size() > 0 & logEnergies) {
			energyFile.open(energyFileName);
			energyFile << "/* energy in 10kcal/mol*/.\n";
		}

		//initialize DotPlot file:
		SC_DotPlot::DotPlot dotplot;
		////////

		///////////   ITERATE EXPLORATIVE LOCAL FLOODING  /////////////

		// initializing the DO list as queue, which contain all states that should be processed
		std::list<size_t> toDo_List;

		// initializing the List of all states, which had already visited
		std::unordered_set<size_t> done_List;

		// insert the index of current minima to toDo_list
		toDo_List.push_back(currentMinID);

		////// init dynamic-Best-K feature /////
		bool mfeFound = false;
		char * mfeStructure = (char*) vrna_alloc(
				(vc->length + 1) * sizeof(char));
		float mfeEnergy = vrna_mfe(vc, mfeStructure);

		//init dynamic k-best filter list.
		std::unordered_map<size_t, std::vector<size_t>> dynamicBestKFilterNeighborList;

		do		//start dynamic best k filter loop.
		{
			// go through all the elements of toDo_list to get the Neighbors of that
			while (!toDo_List.empty()) {
				int listSize = toDo_List.size();
				int numThreads = maxThreads;   //std::min(listSize, maxThreads);

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
						std::pair<flooderInputParameter*,
								flooderOutputParameter*>> threadParameterLists[maxThreads];
				Concurrent_Queue<MyState> discoveredMinimaForEachThread[maxThreads];

				std::future<int> v[maxThreads];

				//std::cout << "Waiting..." << std::flush;
				bool finish = false;
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
							finish = false;
							res = false;
						}
						//std::cout << "" << index << ' ' << res << '\n' << std::flush;
						if (!res) {
							finish = false;
							//test if new discovered minima are on the stack (pull all available minima)
							//TODO: warning! This kind of asynchronous flooding disables some Filter Options.
							if (!enableBestKFilter && !enableDeltaMinEFilter) {
								while (!discoveredMinimaForEachThread[index].empty()) {
									MyState newMin =
											discoveredMinimaForEachThread[index].pop();

                  //start flooding.
                  if (Minima.find(newMin) == Minima.end()) { //if local min is not in total list.
                                         //add new min and set lowest maximal index.
                    size_t lowestMaxIndex = TypeID::value<
                        size_t>();

                    Minima.insert(
                        { newMin, lowestMaxIndex });
                    MinimaForReverseSearch.insert( {
                        lowestMaxIndex, newMin });

                    // if the element at the index position is not already in done List, add it to toDo_list
                    if (done_List.find(lowestMaxIndex)
                        == done_List.end()
                        && std::find(toDo_List.begin(),
                            toDo_List.end(),
                            lowestMaxIndex)
                            == toDo_List.end()) {
                      toDo_List.push_back(lowestMaxIndex);
                    }
                  }
                  //discoveredMinimaForEachThread[index].pop_front();
								}
							}
						} else {
							//finished -> merge result and create a new thread
							std::vector<
									std::pair<flooderInputParameter*,
											flooderOutputParameter*>>* threadParameter =
									&threadParameterLists[index];
							if (threadParameter->size() > 0) {
								merge_results(threadParameter, logEnergies,
										energyFile, addedMinIDs, writeDotplot,
										dotplot, Minima, MinimaForReverseSearch,
										z, dynamicBestK,
										dynamicBestKFilterNeighborList,
										toDo_List, done_List, verbose);

								threadParameter->clear();
							}
							if (!toDo_List.empty()) {
								finish = false;
								//create a new thread
								int listSize = toDo_List.size();
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

									flooderInputParameter* inParameter =
											new flooderInputParameter();
									inParameter->BasinID = minID;
									inParameter->CurrentMinimum =
											new MyState(MinimaForReverseSearch.at(minID));
									inParameter->MaxToHash = maxToHash;
									inParameter->MaxToQueue = maxToQueue;
									inParameter->Filter = neighborFilter;
									inParameter->MaxEnergy = maxEnergy;
									inParameter->DeltaE = deltaE;
									inParameter->DiscoveredMinima =
											&discoveredMinimaForEachThread[index];
									inParameter->TemperatureForBoltzmannWeight =
											temperatureForBoltzmannWeight;
									inParameter->Move_set = move_set;
									inParameter->All_Saddles = &all_saddles;
									inParameter->MaxBPdist = maxBPdist;
									inParameter->SourceStructure = rna_start_str;
									inParameter->TargetStructure = rna_final_str;

									flooderOutputParameter* outParameter =
											new flooderOutputParameter();

									// check what state collector to instantiate
									if (writeDotplot) {
										outParameter->ScBasin = new SC_DotPlot(
												temperatureForBoltzmannWeight,
												logEnergies);
									} else {
										outParameter->ScBasin =
												new SC_PartitionFunction(
														temperatureForBoltzmannWeight,
														logEnergies);
									}
									threadParameter->push_back(
											std::pair<flooderInputParameter*,
													flooderOutputParameter*>(
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
							std::pair<flooderInputParameter*,
									flooderOutputParameter*>>* threadParameter =
							&threadParameterLists[index];
					if (threadParameter->size() > 0) {
						threadParameter->clear();
					}
					discoveredMinimaForEachThread[index].clear();
				}

			} // end while

			if (logEnergies) {
				energyFile.close();
			}

			/////// Search the MFE structure in the done_list ////
			short * mfeStructurePT = vrna_ptable(mfeStructure);
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
					neighborFilter = neighbor_K_Filter;
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
			if (!mfeFound) {
				std::cout << "Warning: The mfe structure could not be found!"
						<< std::endl;
			}
		}

		std::cout << std::endl;
		std::cout << "MFE structure: " << mfeStructure << std::endl;
		std::cout << "The start state is: " << startState.toString()
				<< std::endl;
		std::cout << "The start state ends in basin: "
				<< startStateMinimum->toString() << std::endl;
		std::cout << std::endl;

		// Calculate the final Rate Matrix:
		*transOut
				<< "              --------------THE FINAL RATE MATRIX----------------------- "
				<< std::endl;
		// the maximum number of macrostates is used to normalize the rates.
		size_t maxNeighNum = getMaximalNeighborsOfAMacroState(z, done_List,
				MinimaForReverseSearch);
		biu::MatrixSparseC<double>* final_Rate = &calculateRateMatrix(z,
				done_List, MinimaForReverseSearch, maxNeighNum);

		// To store the States after applying all the Filtering Techniques
		std::unordered_map<size_t, MyState> final_minima;
		std::unordered_set<size_t>::const_iterator fromIt = done_List.begin();
		for (size_t from = 0; from < done_List.size(); from++, fromIt++) {
			// store minimum with correct index for final rate matrix
			// (from column-state to row-state)
			final_minima.insert( { from, MinimaForReverseSearch.at(*fromIt) });
		}

		// Print the Final rate matrix of the States in Final minima set.
		// printRateMatrix (*final_Rate, final_minima, *transOut, true);
		PairHashTable::HashTable* sorted_min_and_output_ids = printRateMatrixSorted(*final_Rate, final_minima, *transOut);
		std::cout << std::endl;
		printEquilibriumDensities(z, final_minima, Minima, *transOut);
		std::cout << std::endl;

		//if parameter given: print Z-matrix to file.
		if (partitionFunctionFileName.size() > 0 & writePartitionFunctions) {
			ofstream ptFile;
			ptFile.open(partitionFunctionFileName);
			ptFile
					<< "              --------------THE PARTITIONFUNCTION MATRIX----------------------- "
					<< std::endl << std::endl;
			printZMatrixSorted(z, maxNeighNum, final_minima, Minima, ptFile);
			ptFile.close();
		}

    double partitionFunctionLandscape = 0;
    //get partition Sum for the energy landscape.
    for (size_t i = 0; i < final_minima.size(); i++) {
      MyState tmpMin = final_minima.at(i);
      size_t minIndex = Minima.at(tmpMin);
      partitionFunctionLandscape +=
          z.at(SC_PartitionFunction::PairID(minIndex, minIndex)).getZ();
    }
    std::cout << "The overall partition function is: "
        << partitionFunctionLandscape << std::endl;

		//dotplot
		if (writeDotplot) {
			SC_DotPlot::DotPlot normalizedDotplot =
					SC_DotPlot::getBasePairProbabilities(dotplot,
							partitionFunctionLandscape);
			bool dotplotWritten = SC_DotPlot::writeDotPlot_PS(dotPlotFileName,
					rnaSequence, normalizedDotplot);
			if (!dotplotWritten) {
				std::cerr
						<< "\nWarning: error during the writing of the dot plot file "
						<< dotPlotFileName << "\n";
			}
		}

		print_number_of_rates (*final_Rate,final_minima, Minima, std::cout);

		if(!binary_rates_file.empty()){
		  write_binary_rates_file(binary_rates_file, *final_Rate,final_minima, Minima);
		}

		/**
		 * print saddles
		 * */
	  std::printf("id_from, loc_min_from, loc_min_from_energy, id_to, loc_min_to, loc_min_to_energy, saddle, saddle_energy\n");
	  for(auto it = all_saddles.begin(); it != all_saddles.end(); it++){
	    const std::pair<MyState,MyState>& from_to = it->first;
	    double saddle_height = it->second.energy/100.0;
	    std::string s1 = from_to.first.toString();
	    std::string s2 = from_to.second.toString();
	    std::string saddle = it->second.toString();
	    size_t id_from = sorted_min_and_output_ids->at(from_to.first);
	    size_t id_to = sorted_min_and_output_ids->at(from_to.second);
	    std::printf("%ld, %s, %.2f, %ld, %s, %.2f, %s, %.2f\n", id_from, s1.c_str(), from_to.first.energy/100.0,
	        id_to, s2.c_str(), from_to.second.energy/100.0, saddle.c_str(), saddle_height);
	  }


		/***********Garbage collection ************/
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
		final_minima.clear();
		z.clear();
		if (startStateMinimum != NULL)
			delete startStateMinimum;
		if (final_Rate != NULL)
			delete final_Rate;
	} catch (std::exception & ex) {
		std::cerr << "\n\n ERORR : " << ex.what() << "\n" << std::endl;
		return_Value = -1;
	}


//print stopwatch:
	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	std::cout << "finished computation at " << std::ctime(&end_time)
			<< "elapsed time: " << elapsed_seconds.count() << "s\n";
	return return_Value;

}	  // end main body.

/*!writes the description of the program to the "info" string and
 * writes the allowed parameters in the OptionMap "allowed".
 * @param allowed, list of possible options.
 * @param infos contains description of the options.
 */
void writeDescription(biu::OptionMap & allowed, std::string & info) {

	std::stringstream info_string;

//description
	info_string << "\n"
			<< "RNA_Expore takes an RNA sequence as input and explores the landscape topology locally. "
					"This means the flooding algorithm will be applied for each gradient basin. "
					"The partition function for the basin and also for the transitions to neighbored minima will "
					"be calculated during the flooding. This approach consumes more computations than "
					"global flooding, because the contact surface for two neighbored states is calculated twice. "
					"The advantage of RNA_Explore is, that the filtering techniques can be calculated locally. "
					"The filters are used to prune non-relevant transitions directly after flooding a gradient basin. "
					"As a result, the transition rates for the filtered landscape topology can be calculated faster "
					"than with global approaches. The advantage increases with increasing size of the energy landscape."
					"      " << "\n";

	info = info_string.str();

	allowed.push_back(
			biu::COption("seq", false, biu::COption::STRING,
					"The RNA sequence of the molecule", "ACUGUAUGCGCGU"));
	allowed.push_back(
			biu::COption("startStr", true, biu::COption::STRING,
					"the start structure of the exploration defining the first gradient basin; defaults"
							" to the open chain"));
	allowed.push_back(
			biu::COption("finalStr", true, biu::COption::STRING,
					"the final structure of the exploration defining the last gradient basin"));
	allowed.push_back(
			biu::COption("filterBestK", true, biu::COption::INT,
					"reduces outgoing transitions to the best K for each gradient basin"));
	allowed.push_back(
			biu::COption("dynamicBestK", true, biu::COption::BOOL,
					"Increases K if the MFE structure is not explored."));
	allowed.push_back(
			biu::COption("filterMinE", true, biu::COption::DOUBLE,
					"reduces outgoing transitions to the Neighbor minima, which have  "
							" Energy lower than current (E(minima)+filterValue) for each gradient basin"));
	allowed.push_back(
			biu::COption("maxToQueue", true, biu::COption::INT,
					"Sets the maximum number of states to be stored in the priority queue of the flooder.",
					std::to_string(SIZE_MAX)));
	allowed.push_back(
			biu::COption("maxToHash", true, biu::COption::INT,
					"Sets the maximum number of states to be hashed for a gradient walk.",
					std::to_string(SIZE_MAX)));
	allowed.push_back(
			biu::COption("dynamicMaxToHash", true, biu::COption::BOOL,
					"Sets the dynamicMaxToHash variable for estimating the maximal number of states to be hashed in a gradient walk, "
							"by considering the maximal available physical memory and the number of threads. This reduces the probability of swapping."));

	allowed.push_back(
			biu::COption("maxEnergy", true, biu::COption::DOUBLE,
					"Sets the maximum energy that a state is allowed to have to be considered by the flooder (in kcal/mol).",
					"5"));
	allowed.push_back(
			biu::COption("deltaE", true, biu::COption::DOUBLE,
					"Set the maximum energy difference that states in a basin can have w.r.t. the local minimum (in kcal/mol).",
					"9999"));
	allowed.push_back(
			biu::COption("T", true, biu::COption::INT,
					"Set the temperature for the free energy calculation (in °C). (If \"T\" is set and \"B\" not, \"B\" is equals \"T\").",
					"37"));
	allowed.push_back(
			biu::COption("d", true, biu::COption::INT,
					"How to treat \"dangling end\" energies for bases adjacent to helices in free ends and multi-loops.",
					"2"));
	allowed.push_back(
			biu::COption("B", true, biu::COption::INT,
					"Set the temperature for the Boltzmann weight (in °C).",
					"37"));
	allowed.push_back(
			biu::COption("M", true, biu::COption::INT,
					"Set the energy model. 0=Turner model 2004, 1=Turner model 1999, 2=Andronescu model, 2007",
					"0"));
  allowed.push_back(
      biu::COption("move_set", true, biu::COption::INT,
          "Move set: 0 = insertion and deletion, 1 = shift moves, 2 = no lonely pair moves.",
          "0"));
	allowed.push_back(
			biu::COption("transProb", true, biu::COption::STRING,
					"If provided, the transition probability matrix will be written to the given file name or 'STDOUT' when to write to standard output"));
	allowed.push_back(
			biu::COption("energyFile", true, biu::COption::STRING,
					"File to store all energies."));
  allowed.push_back(
      biu::COption("binary_rates_file", true, biu::COption::STRING,
          "File to store all rates in a treekin readable format."));
	allowed.push_back(
			biu::COption("partitionFunctions", true, biu::COption::STRING,
					"If provided, the partition function matrix will be written to the given file name."));
	allowed.push_back(
			biu::COption("dotPlot", true, biu::COption::STRING,
					"If provided, the dotPlot will be written to the given file name. \
	  The dotPlot contains the base pair probabilities for all structures in the (filtered) energy landscape."));
  allowed.push_back(
      biu::COption("dotPlot_per_basin", true, biu::COption::STRING,
          "Creates a dotplot for each gradient basin in the enrgy landscape. It shows the Maximum Expected Accuracy (MEA) structure \
           in the upper right triangle and the basin representative in the lower left triangle"));
	allowed.push_back(
			biu::COption("maxThreads", true, biu::COption::INT,
					"Sets the maximum number of threads for parallelized computation.",
					"1"));
	allowed.push_back(
	      biu::COption("maxBPdist", true, biu::COption::INT,
	          "Sets the maximum base pair distance for direct neighbor minima to be explored.",
	          "65536"));
	allowed.push_back(
			biu::COption("v", true, biu::COption::BOOL, "Verbose"));

}		// end function

