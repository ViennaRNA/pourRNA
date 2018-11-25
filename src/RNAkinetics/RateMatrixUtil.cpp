/*
 * RateMatrixUtil.cpp
 *
 *  Created on: 09.08.2014
 *      Author: From RNAKinetics copied and adjusted by Gregor Entzian
 */

#include "RateMatrixUtil.h"

void printRateMatrix(const biu::MatrixSparseC<double>& R,
    const std::unordered_map<size_t, MyState> & minimaMap, std::ostream & out,
    const bool noZeros) {
  assertbiu(R.numColumns() == R.numRows(), "R is no square matrix");
  assertbiu(R.numRows() <= minimaMap.size(), "less minima than rates");

  const size_t LEAD = 6;

  if (noZeros) {

    out << std::scientific << "\n from : to\n";

    // print only non-empty rates
    for (size_t c = 0; c < R.numColumns(); c++) {
      out << "\n" << std::setw(LEAD) << c << " [" << minimaMap.at(c).toString()
          << "] :";

      const biu::MatrixSparseC<double>::EntryMap col = R.columnValues(c);
      biu::MatrixSparseC<double>::EntryMap::const_iterator row;
      for (row = col.begin(); row != col.end(); row++) {
        if (row != col.begin()) {
          out << ",";
        }
        out << " " << row->first << " = " << row->second;
      }
    }
    out << std::endl;

  } else {
    const size_t PREC = 6;
    double delta = 0.1;
    for (size_t i = 0; i < PREC; i++)
      delta = delta / 10.0;

    out << "\n" << std::setw(LEAD) << "from |";
    for (size_t i = 0; i < R.numColumns(); i++) {
      out << std::setw(PREC + 3) << i;
    }
    out << "\n" << std::setfill('-') << std::setw(LEAD) << "-";
    for (size_t i = 0; i < R.numColumns(); i++) {
      out << std::setfill('-') << std::setw(PREC + 3) << "-";
    }
    out << "\n" << std::setfill(' ') << "to 0 |";
    for (size_t from = 0; from < R.numColumns(); from++) {
      if (noZeros && std::abs(R.at(0, from)) < delta) {
        out << std::setw(PREC + 3) << " ";
      } else {
        out << std::fixed << std::setprecision(PREC) << std::setw(PREC + 3)
            << R.at(0, from);
      }
    }
    for (size_t to = 1; to < R.numRows(); to++) {
      out << "\n" << std::setw(LEAD - 2) << to << " |";
      for (size_t from = 0; from < R.numColumns(); from++) {
        if (noZeros && std::abs(R.at(to, from)) < delta) {
          out << std::setw(PREC + 3) << " ";
        } else {
          out << std::fixed << std::setprecision(PREC) << std::setw(PREC + 3)
              << R.at(to, from);
        }
      }
    }
    out << std::endl;
  }
}

struct less_second {
  typedef std::pair<size_t, MyState*> type;
  bool operator ()(type const& a, type const& b) const {
    // check if
    // - smaller energy or
    // - equal energy and smaller string representation
    return (a.second->getEnergy() < b.second->getEnergy())
        || (a.second->getEnergy() == b.second->getEnergy()
            && StructureUtils::IsSmaller(a.second->structure,
                b.second->structure));
  }
};

PairHashTable::HashTable* printRateMatrixSorted(const biu::MatrixSparseC<double>& R,
    const std::unordered_map<size_t, MyState>& minimaMap, std::ostream& out) {
  assertbiu(R.numColumns() == R.numRows(), "R is no square matrix");
  assertbiu(R.numRows() <= minimaMap.size(), "less minima than rates");

  const size_t LEAD = 6;

  std::vector<std::pair<size_t, MyState*>> sortedMinimaIDs;
  for (auto it = minimaMap.begin(); it != minimaMap.end(); it++) {
    sortedMinimaIDs.push_back(
        std::pair<size_t, MyState*>(it->first, (MyState*) &(it->second)));
  }

  std::sort(sortedMinimaIDs.begin(), sortedMinimaIDs.end(), less_second());
  PairHashTable::HashTable* states_and_output_ids = new PairHashTable::HashTable();
  out << std::scientific << "\n from : to\n";
  size_t nextMinID;
  // print only non-empty rates
  for (size_t c = 0; c < sortedMinimaIDs.size(); c++) {
    nextMinID = sortedMinimaIDs[c].first;
    out << "\n" << std::setw(LEAD) << c << " ["
        << minimaMap.at(nextMinID).toString() << "] :";
    states_and_output_ids->insert({MyState(minimaMap.at(nextMinID)), c});
    std::vector<double> columnVector = R.columnVec(nextMinID);
    bool bPrinted = false;
    size_t rowMinID;
    std::stringstream sstmp;
    sstmp << std::scientific;
    for (size_t r = 0; r < sortedMinimaIDs.size(); r++) {
      rowMinID = sortedMinimaIDs[r].first;
      if (columnVector[rowMinID] != 0.0) {
        sstmp << " " << r << " = " << columnVector[rowMinID] << ",";
        bPrinted = true;
      }
    }
    if (bPrinted) { //remove last comma.
      sstmp.seekp(-1, sstmp.cur);
      sstmp << " ";
      out << sstmp.str();
    }
  }
  out << std::endl;

  return states_and_output_ids;
}

void write_binary_rates_file(std::string rates_file,
    const biu::MatrixSparseC<double>& R,
    const std::unordered_map<size_t, MyState>& minimaMap,
    const PairHashTable::HashTable& originalMinima) {
  assertbiu(R.numColumns() == R.numRows(), "R is no square matrix");
  assertbiu(R.numRows() <= minimaMap.size(), "less minima than rates");

  FILE *BINOUT;
  const char *binfile = rates_file.c_str(); //"rates.bin";
  double tmprate;
  BINOUT = fopen(binfile, "w");
  if (!BINOUT) {
    fprintf(stderr, "could not open file pointer 4 binary outfile\n");
    exit(101);
  }

  std::vector<std::pair<size_t, MyState*>> sortedMinimaIDs;
  for (auto it = minimaMap.begin(); it != minimaMap.end(); it++) {
    sortedMinimaIDs.push_back(
        std::pair<size_t, MyState*>(it->first, (MyState*) &(it->second)));
  }

  std::sort(sortedMinimaIDs.begin(), sortedMinimaIDs.end(), less_second());

  size_t n = sortedMinimaIDs.size();
  /* first write dim to file */
  fwrite(&n, sizeof(int), 1, BINOUT);

  size_t nextMinID;
  // print only non-empty rates
  double rate = 0.0;
  for (size_t c = 0; c < sortedMinimaIDs.size(); c++) {
//MyState min = *sortedMinimaIDs[c].second;
    nextMinID = sortedMinimaIDs[c].first; //originalMinima.at(min);

    std::vector<double> columnVector = R.rowVec(nextMinID);
    bool bPrinted = false;
    size_t rowMinID;
    for (size_t r = 0; r < sortedMinimaIDs.size(); r++) {
      rowMinID = sortedMinimaIDs[r].first;
      rate = columnVector[rowMinID];
      fwrite(&rate, sizeof(double), 1, BINOUT);
    }
  }

  fprintf(stderr, "rate matrix written to binfile\n");
  fclose(BINOUT);
}

void print_number_of_rates(const biu::MatrixSparseC<double>& R,
    const std::unordered_map<size_t, MyState>& minimaMap,
    const PairHashTable::HashTable& originalMinima, std::ostream& out) {
  assertbiu(R.numColumns() == R.numRows(), "R is no square matrix");
  assertbiu(R.numRows() <= minimaMap.size(), "less minima than rates");

  size_t count_rates = 0;
  size_t nextMinID;
  size_t rowMinID;
  double rate;
  // count only non-empty rates
  for (auto it = minimaMap.begin(); it != minimaMap.end(); it++){
    nextMinID = it->first;
    std::vector<double> columnVector = R.columnVec(nextMinID);
    for (auto it_r = minimaMap.begin(); it_r != minimaMap.end(); it_r++){
      rowMinID = it_r->first;
      rate = columnVector[rowMinID];
      if (rate != 0.0) {
        count_rates += 1;
      }
    }
  }
  out << "number of rates: " << count_rates << std::endl;
  out << std::endl;
}

void printZMatrixSorted(const SC_PartitionFunction::Z_Matrix& z,
    size_t maxNeighbors, const std::unordered_map<size_t, MyState>& minimaMap,
    const PairHashTable::HashTable& originalMinima, std::ostream& out) {
  const size_t LEAD = 6;

  std::vector<std::pair<size_t, MyState*>> sortedMinimaIDs;
  for (auto it = minimaMap.begin(); it != minimaMap.end(); it++) {
    sortedMinimaIDs.push_back(
        std::pair<size_t, MyState*>(it->first, (MyState*) &(it->second)));
  }

  std::sort(sortedMinimaIDs.begin(), sortedMinimaIDs.end(), less_second());

  out << std::scientific << "# max. neighbors: " << maxNeighbors
  << "\n from : to\n";
  size_t nextMinID;
  // print only non-empty rates
  for (size_t c = 0; c < sortedMinimaIDs.size(); c++) {
    MyState min = *sortedMinimaIDs[c].second;
    nextMinID = originalMinima.at(min);
    out << "\n" << std::setw(LEAD) << c << " [" << min.toString() << "] :";

    bool bPrinted = false;
    size_t rowMinID;
    std::stringstream sstmp;
    sstmp << std::scientific;
    for (size_t r = 0; r < sortedMinimaIDs.size(); r++) {
      rowMinID = originalMinima.at(*sortedMinimaIDs[r].second);
      //in Z-Matrix: from i to j.
      //in final_Rate: from j to i !
      SC_PartitionFunction::PairID transitionID = SC_PartitionFunction::PairID(
          nextMinID, rowMinID);
      SC_PartitionFunction::PairID reverseTransitionID =
      SC_PartitionFunction::PairID(rowMinID, nextMinID);

      double z_transition = 0;
      double z_reverseTransition = 0;
      auto tmpIt = z.find(transitionID);
      bool bExists = false;
      if (tmpIt != z.end()) {
        z_transition = tmpIt->second.getZ();
        bExists = true;
      }
      tmpIt = z.find(reverseTransitionID);
      if (tmpIt != z.end()) {
        z_reverseTransition = tmpIt->second.getZ();
        bExists = true;
      }
      //max is important if filters are applied! Otherwise it should be equal.
      z_transition = std::max(z_transition, z_reverseTransition);
      // if (nextMinID != rowMinID)
      //   { // compute correct transition to another state.
      //     z_transition = z_transition / maxNeighbors;
      //   }

      if (bExists) {
        sstmp << " " << r << " = " << z_transition << ",";
        bPrinted = true;
      }
    }
    if (bPrinted) { //remove last comma.
      sstmp.seekp(-1, sstmp.cur);
      sstmp << " ";
      out << sstmp.str();
    }
  }
  out << std::endl;

}

void printEquilibriumDensities(SC_PartitionFunction::Z_Matrix& z,
    const std::unordered_map<size_t, MyState>& finalMinima,
    const PairHashTable::HashTable& originalMinima, std::ostream& out) {
  out << "Equilibrium Densities:" << std::endl;

  std::vector<std::pair<size_t, MyState*>> sortedMinimaIDs;
  for (auto it = finalMinima.begin(); it != finalMinima.end(); it++) {
    sortedMinimaIDs.push_back(
        std::pair<size_t, MyState*>(it->first, (MyState*) &(it->second)));
  }

  std::sort(sortedMinimaIDs.begin(), sortedMinimaIDs.end(), less_second());

  size_t nextMinID;
  double sumZb = 0;
  // calc. sum of all basin partition functions.
  for (size_t c = 0; c < sortedMinimaIDs.size(); c++) {
    nextMinID = originalMinima.at(*sortedMinimaIDs[c].second);
    sumZb += z[SC_PartitionFunction::PairID(nextMinID, nextMinID)].getZ();
  }
  double equilibriumDensity;
  out << "(";
  std::stringstream sstmp;
  sstmp << std::scientific;
  for (size_t c = 0; c < sortedMinimaIDs.size(); c++) {
    nextMinID = originalMinima.at(*sortedMinimaIDs[c].second);
    equilibriumDensity =
    z[SC_PartitionFunction::PairID(nextMinID, nextMinID)].getZ() / sumZb;
    // print probability and state
    sstmp << equilibriumDensity;
    // print spacer if needed

    if (c < sortedMinimaIDs.size() - 1) {
      sstmp << ", ";
    }
  }
  out << sstmp.str() << ")" << std::endl;
}
