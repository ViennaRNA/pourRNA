/*
 * RateMatrixUtil.h
 *
 *  Created on: 09.08.2014
 *      Author: From RNAKinetics copied and adjusted by Gregor Entzian
 */

#ifndef RATEMATRIXUTIL_H_
#define RATEMATRIXUTIL_H_

#include <iostream>
#include <iomanip>
#include <math.h>
#include "../BIUlibPart/MatrixSparse.hh"
#include <unordered_map>
extern "C" {
#include <ViennaRNA/move_set.h>
}
#include "../StructureUtils.h"
#include "../MyState.h"
#include <algorithm> //for sort, max function
#include "../SC_PartitionFunction.h"
#include <sstream>
#include <cmath>
#include "../PairHashTable.h"

/**
 * ! print the rate matrix.
 * @param R the rateMatrix.
 * @param minima the indices and structures of the minima.
 * @param out the outputstream, default is std::cout.
 * @param noZeros defines if zero-rates should be printed or not.
 */
void
printRateMatrix (const biu::MatrixSparseC<double>& R, const std::unordered_map<size_t, MyState> & minima,
		 std::ostream & out = std::cout, const bool noZeros = true);

/**
 * ! print the rate matrix without zero rates and sorted.
 * @param R the rateMatrix.
 * @param minima the indices and structures of the minima.
 * @param out the outputstream, default is std::cout.
 */
void
printRateMatrixSorted (const biu::MatrixSparseC<double>& R, const std::unordered_map<size_t, MyState>& minima,
		       std::ostream & out = std::cout);

/**
 * ! print the partitionfunction sparse matrix, sorted.
 * @param z the partitionfunctions.
 * @param maxNeighbors the maximal number of neighbors is used to calculate the correct transition partition sums.
 * @param minima the indices and structures of the minima.
 * @param originalMinima all unfiltered discovered minima with indices that are consistent with the z matrix.
 * @param out the outputstream, default is std::cout.
 */
void
printZMatrixSorted (const SC_PartitionFunction::Z_Matrix& z, size_t maxNeighbors, const std::unordered_map<size_t, MyState>& minimaMap,
		    const PairHashTable::HashTable& originalMinima, std::ostream& out = std::cout);

/**
 * ! print a vector with the equilibrium densities (Zb/sum(Zb')).
 * @param z the matrix with all transition partition sums and indices that are consistent with original minima.
 * @param finalMinima all (maybe filtered) minima that are in the rate matrix.
 * @param originalMinima all unfiltered discovered minima with indices that are consistent with the z matrix.
 * @param out the outputstream, default is std::cout.
 */
void
printEquilibriumDensities (SC_PartitionFunction::Z_Matrix& z,
			   const std::unordered_map<size_t, MyState>& finalMinima,
			   const PairHashTable::HashTable& originalMinima, std::ostream& out);

#endif /*RATEMATRIXUTIL_H_*/
