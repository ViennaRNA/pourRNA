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
#include <vector>

/**
 * ! print the rate matrix.
 * @param R the rateMatrix.
 * @param minima the indices and structures of the minima.
 * @param out the outputstream, default is std::cout.
 * @param noZeros defines if zero-rates should be printed or not.
 */
void
printRateMatrix(const biu::MatrixSparseC<double>& R,
                const std::unordered_map<size_t, MyState> & minima,
                std::ostream & out  = std::cout,
                const bool noZeros  = true);


/**
 * ! print the rate matrix without zero rates and sorted.
 * @param R the rateMatrix.
 * @param minima the indices and structures of the minima.
 * @param mea a vector with the corresponding mea structures for each local minimum.
 * @param out the outputstream, default is std::cout.
 */
void
printRateMatrixSorted(const biu::MatrixSparseC<double>& R,
                      const std::vector<std::pair<size_t, MyState *> >& sortedMinimaIDs,
                      std::ostream & out = std::cout);


void
printRateMatrixSortedWithMEA(const biu::MatrixSparseC<double>& R,
                             const std::vector<std::pair<size_t, MyState *> >& sortedMinimaIDs,
                             const std::vector<std::string>& mea,
                             std::ostream & out = std::cout);


/**
 *  ! write binary rates file.
 *  first value is the number of states, then the all outgoing rates for state 1, then all for state 2 etc.
 * @param rates_file the name of the rate matrix file.
 * @param R the rateMatrix.
 * @param minima the indices and structures of the minima.
 */
void
write_binary_rates_file(std::string rates_file,
                        const biu::MatrixSparseC<double>& R,
                        const std::vector<std::pair<size_t, MyState *> >& sortedMinimaIDs);


/**
 *  ! write sparse binary rates file.
 *  first value is the number of states, then <integer from>, <integer number of how many value pairs to>,
 *  <value pair <integer to, double rate from, to>> etc.
 * @param rates_file the name of the rate matrix file.
 * @param R the rateMatrix.
 * @param minima the indices and structures of the minima.
 */
void
write_binary_rates_file_sparse(std::string rates_file,
                               const biu::MatrixSparseC<double>& R,
                               const std::vector<std::pair<size_t, MyState *> >& sortedMinimaIDs);


/**
 *  ! write a text file that contains the rate matrix and a text file that contains the structures.
 */
void
write_barriers_like_output(std::string file_prefix,
                           const biu::MatrixSparseC<double>& R,
                           const std::vector<std::pair<size_t, MyState *> >& sortedMinimaIDs,
                           std::string sequence,
                           bool basin_size,
                           std::unordered_map<size_t, size_t> minimum_index_and_basin_size);


/**
 * ! print the number of non-zero rates.
 * @param R the rateMatrix.
 * @param minima the indices and structures of the minima.
 * @param out the outputstream, default is std::cout.
 */
void
print_number_of_rates(const biu::MatrixSparseC<double>& R,
                      const std::vector<std::pair<size_t, MyState *> >& minimaMap,
                      std::ostream& out);


/**
 * ! print the partitionfunction sparse matrix, sorted.
 * @param z the partitionfunctions.
 * @param minima the indices and structures of the minima.
 * @param originalMinima all unfiltered discovered minima with indices that are consistent with the z matrix.
 * @param out the outputstream, default is std::cout.
 */
void
printZMatrixSorted(const SC_PartitionFunction::Z_Matrix& z,
                   const std::vector<std::pair<size_t, MyState *> >& sortedMinimaIDs,
                   const PairHashTable::HashTable& originalMinima,
                   std::ostream& out = std::cout);


/**
 * ! print a vector with the equilibrium densities (Zb/sum(Zb')).
 * @param z the matrix with all transition partition sums and indices that are consistent with original minima.
 * @param finalMinima all (maybe filtered) minima that are in the rate matrix.
 * @param out the outputstream, default is std::cout.
 */
void
printEquilibriumDensities(SC_PartitionFunction::Z_Matrix& z,
                          const std::vector<std::pair<size_t, MyState *> >& sortedMinimaIDs,
                          std::ostream& out);


#endif /*RATEMATRIXUTIL_H_*/
