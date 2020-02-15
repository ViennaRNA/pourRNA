/*
 * BarriersTree.h
 *
 *  Created on: Jan 26, 2020
 *      Author: Gregor Entzian
 */

#ifndef SRC_BARRIERSTREE_H_
#define SRC_BARRIERSTREE_H_

#include <vector>
#include <unordered_map>
#include <string>
#include "MyState.h"
#include "StatePairCollector.h"
#include "PairHashTable.h"

typedef struct node_{
  node_* parent;
  node_* left_child;
  node_* right_child;
  size_t minimum_id;
  int minimum_energy;
  int saddle;
  int branch_length;
  size_t node_index;
} node_t;

typedef struct saddle_{
  size_t minimum_from;
  size_t minimum_to;
  int saddle;
} saddle_t;

class BarriersTree {
public:
  BarriersTree();
  virtual ~BarriersTree();
  std::vector<node_t*> create_barrier_tree(std::vector<saddle_t> minimal_saddle_list, std::unordered_map<size_t, int> structure_index_to_energy);
  std::string newick_string_builder(node_t* root);
  std::string svg_string_builder(node_t *tree, int mfe, size_t *inorder_index);
  std::vector<saddle_t> create_minimal_saddle_list(const std::vector<std::pair<size_t, MyState *> > &sortedMinimaIDs,
      const PairHashTable::HashTable &sorted_min_and_output_ids, const StatePairCollector::MapOfMaps &all_saddles);
  double determin_optimal_min_h(const size_t maximal_number_of_states, const std::vector<saddle_t> &minimal_saddle_list, const std::vector<std::pair<size_t, MyState *>> &sortedMinimaIDs);
  std::unordered_map<size_t, size_t> filter_minh(std::vector<saddle_t> &minimal_saddle_list, const std::vector<std::pair<size_t, MyState *>> &sortedMinimaIDs, const double min_h);
  void free_tree(node_t* root);
};

#endif /* SRC_BARRIERSTREE_H_ */
