/*
 * BarriersTree.cpp
 *
 *  Created on: Jan 26, 2020
 *      Author: Gregor Entzian
 */

#include "BarriersTree.h"
#include <unordered_set>
#include <sstream>
#include <algorithm>
#include <iostream>

BarriersTree::BarriersTree() {
  // TODO Auto-generated constructor stub

}

BarriersTree::~BarriersTree() {
  // TODO Auto-generated destructor stub
}

struct pair_hash {
  template<class T1, class T2>
  std::size_t operator ()(std::pair<T1, T2> const &pair) const {
    std::size_t h1 = std::hash<T1>()(pair.first);
    std::size_t h2 = std::hash<T2>()(pair.second);

    return h1 ^ h2;
  }
};

void print_set(std::unordered_set<size_t> & set){
  printf("{");
  for(auto it = set.begin(); it != set.end(); it++)
    printf("%ld,",*it);
  printf("}");
}

void print_sets(std::vector<std::unordered_set<size_t>*>& clusters_id_sets){
  for(size_t t= 0; t < clusters_id_sets.size(); t++){
    print_set(*clusters_id_sets[t]);
    printf(", ");
  }
}

/**
 * ! create a minimal saddle list for each state. For each state it outputs the lowest saddle to a deeper minimum.
 *   The input is a hash table that contains all saddles to direct neighbors.
 *   Since direct neighbors do not contain all saddles, this procedure computes the transitive hull,
 *   in order to get the lowest saddle for all paths between minima.
 *   @param sortedMinimaIDs - a vector with ids and states as for all minima in the output.
 *   @param sorted_min_and_output_ids - a hash map, which maps the minima to ids.
 *   @param all_saddles - a hash map, that contains all saddle states.
 *   @return a vector of minimal outgoing saddles to a deeper minimum for each state (indices according to output IDs.).
 */
std::vector<saddle_t> BarriersTree::create_minimal_saddle_list(const std::vector<std::pair<size_t, MyState *> > &sortedMinimaIDs,
    const PairHashTable::HashTable &sorted_min_and_output_ids, const StatePairCollector::MapOfMaps &all_saddles){
  // init a saddle matrix. Indices from 0 to n according to output IDs.
  int** all_to_all_saddles = new int*[sorted_min_and_output_ids.size()];
  for(auto it = sorted_min_and_output_ids.begin(); it != sorted_min_and_output_ids.end(); it++){
    all_to_all_saddles[it->second] = new int[sorted_min_and_output_ids.size()];
    // init with inf distance
    std::fill_n(all_to_all_saddles[it->second], sorted_min_and_output_ids.size(), INT32_MAX);
  }

  // insert all saddles
  for (auto it = all_saddles.begin(); it != all_saddles.end(); it++) {
    auto it_from = sorted_min_and_output_ids.find(it->first);
    if (it_from != sorted_min_and_output_ids.end()){
      size_t state_from = it_from->second;
      for (auto it_to = it->second.begin(); it_to != it->second.end(); it_to++){
        auto it_to_min = sorted_min_and_output_ids.find(it_to->first);
        if (it_to_min != sorted_min_and_output_ids.end()){
          int saddle_energy = it_to->second.energy;
          size_t state_to = it_to_min->second;
          all_to_all_saddles[state_from][state_to] = saddle_energy;
          all_to_all_saddles[state_to][state_from] = saddle_energy;
        }
      }
    }
  }

  /* now insert all saddles (inclusive indirect connections).
   * warshall algorithm to compute the transitive hull (O(n^3) time):
   * for k = 1 to n
   *   for i = 1 to n
   *     if d[i,k] == 1
   *       for j = 1 to n
   *         if d[k,j] == 1
   *           d[i,j] = 1
   *
   * however we adjust d[i,j] = min(d[i,j], max(d[i,k], d[k,j]))
   * i.e. the minimal maximal saddle on the path.
   * and if d[i,k] < INF
   */
  size_t n = sorted_min_and_output_ids.size();
  for(size_t k = 0; k < n; k++){
    for(size_t i = 0; i < n; i++){
      if (all_to_all_saddles[i][k] < INT32_MAX){
        for(size_t j = 0; j < n; j++){
          if (all_to_all_saddles[k][j] < INT32_MAX){
            //if (all_to_all_saddles[i][j] == INT32_MAX){
              all_to_all_saddles[i][j] = std::min(all_to_all_saddles[i][j], std::max(all_to_all_saddles[i][k], all_to_all_saddles[k][j]));
              all_to_all_saddles[j][i] = std::min(all_to_all_saddles[i][j], std::max(all_to_all_saddles[i][k], all_to_all_saddles[k][j]));
           // }
          }
        }
      }
    }
  }


  std::vector<saddle_t> minimal_saddle_list;
  // for each min insert the minimal saddle to a deeper min (they are already sorted by energy --> only from higher to lower.
  for(size_t t = n; t > 0; t--){
    size_t i = t-1;
    saddle_t min_saddle;
    min_saddle.saddle = INT32_MAX;
    for(size_t j = 0; j < i; j++){
      if (all_to_all_saddles[i][j] < INT32_MAX){
        if (sortedMinimaIDs[j].second->energy < sortedMinimaIDs[i].second->energy){
          if (all_to_all_saddles[i][j] < min_saddle.saddle){
            min_saddle.minimum_from = i;
            min_saddle.minimum_to = j;
            min_saddle.saddle = all_to_all_saddles[i][j];
          }
        }
      }
    }
    if (min_saddle.saddle < INT32_MAX){
      minimal_saddle_list.push_back(min_saddle);
    }
  }

  for(size_t i = 0; i < n; i++){
    delete[] all_to_all_saddles[i];
  }
  delete[] all_to_all_saddles;
  return minimal_saddle_list;
}

bool compare_saddle(saddle_t& a, saddle_t& b)
{
    return (a.saddle < b.saddle);
}

std::vector<node_t*> BarriersTree::create_barrier_tree(
    std::vector<saddle_t> minimal_saddle_list,
    std::unordered_map<size_t, int> structure_index_to_energy) {
  std::sort(minimal_saddle_list.begin(), minimal_saddle_list.end(), compare_saddle);

  std::vector<node_t*> clusters;
  std::vector<std::unordered_set<size_t>*> clusters_id_sets; // list which contains as many sets as clusters, which contain the leaf ids.
  std::unordered_set<std::pair<size_t, size_t>, pair_hash> has_updated_length_for_ids; // = []
  std::unordered_set<size_t> clustered_ids; // = set()

  for (size_t i = 0; i < minimal_saddle_list.size(); i++) {
    saddle_t s = minimal_saddle_list[i];
    size_t id_from = s.minimum_from;
    size_t id_to = s.minimum_to;
    int saddle_energy = s.saddle;

    auto it_id_from = clustered_ids.find(id_from);
    auto it_id_to = clustered_ids.find(id_to);
    if ((it_id_from != clustered_ids.end())
        && (it_id_to != clustered_ids.end())) {
      //update lengths of two branches, if not already done (use lowest saddle of sorted saddle list)
      std::pair<size_t, size_t> id_pair = { id_from, id_to };
      auto it_id_pair = has_updated_length_for_ids.find(id_pair);
      if (it_id_pair == has_updated_length_for_ids.end()) {
        //get branches where i is in one and j in the other
        size_t c_id_to = UINT64_MAX;
        size_t c_id_from = UINT64_MAX;
        size_t j = 0;
        for (auto it_clust_id = clusters_id_sets.begin();
            it_clust_id != clusters_id_sets.end(); it_clust_id++, j++) {
          auto it_contains_from = (*it_clust_id)->find(id_from);
          if (it_contains_from != (*it_clust_id)->end()) {
            c_id_from = j;
            break;
          }
        }
        size_t k = 0;
        for (auto it_clust_id = clusters_id_sets.begin();
            it_clust_id != clusters_id_sets.end(); it_clust_id++, k++) {
          auto it_contains_from = (*it_clust_id)->find(id_to);
          if (it_contains_from != (*it_clust_id)->end()) {
            c_id_to = k;
            break;
          }
        }
        if (c_id_from == UINT64_MAX || c_id_to == UINT64_MAX) {
          if (c_id_from == UINT64_MAX) {
            fprintf(stderr, "Error: handled node is not in tree. %ld\n",
                id_from);
          }
          if (c_id_to == UINT64_MAX) {
            fprintf(stderr, "Error: handled node is not in tree. %ld\n", id_to);
          }

          std::exit(EXIT_FAILURE);
        }
        if (c_id_from != c_id_to) {
          clusters[c_id_from]->branch_length = saddle_energy
              - clusters[c_id_from]->minimum_energy;
          clusters[c_id_to]->branch_length = saddle_energy
              - clusters[c_id_to]->minimum_energy;

          node_t *parent = new node_t();
          parent->branch_length = 0;
          parent->minimum_energy = saddle_energy;
          parent->saddle = saddle_energy;
          parent->minimum_id = 0;
          // parent -> children
          parent->left_child = clusters[c_id_from];
          parent->right_child = clusters[c_id_to];
          clusters[c_id_from]->parent = parent;
          clusters[c_id_to]->parent = parent;

          clusters.erase(clusters.begin() + std::max(c_id_from, c_id_to));
          clusters.erase(clusters.begin() + std::min(c_id_from, c_id_to));
          clusters.push_back(parent);
          std::unordered_set<size_t> *set_ids_a = clusters_id_sets[std::max(
              c_id_from, c_id_to)];
          std::unordered_set<size_t> *set_ids_b = clusters_id_sets[std::min(
              c_id_from, c_id_to)];
          clusters_id_sets.erase(
              clusters_id_sets.begin() + std::max(c_id_from, c_id_to));
          clusters_id_sets.erase(
              clusters_id_sets.begin() + std::min(c_id_from, c_id_to));
          std::unordered_set<size_t> *union_a_b =
              new std::unordered_set<size_t>();
          union_a_b->insert(set_ids_a->begin(), set_ids_a->end());
          union_a_b->insert(set_ids_b->begin(), set_ids_b->end());
          clusters_id_sets.push_back(union_a_b);
          delete set_ids_a;
          delete set_ids_b;
        } else {
          // update lengths in same cluster
          clusters[c_id_from]->branch_length = saddle_energy
              - clusters[c_id_from]->minimum_energy;
        }
        has_updated_length_for_ids.insert(id_pair);
      }
      //std::cout << "both in cluster " << id_from << " " << id_to << std::endl;
      continue;
    }

    if ((it_id_from == clustered_ids.end())
        && (it_id_to == clustered_ids.end())) {
      // create new cluster with both ids
      int energy_from = structure_index_to_energy[id_from];
      int energy_to = structure_index_to_energy[id_to];
      node_t *node_from = new node_t();
      node_t *node_to = new node_t();
      node_from->minimum_energy = energy_from;
      node_from->saddle = saddle_energy;
      node_from->branch_length = saddle_energy - energy_from;
      node_from->minimum_id = id_from;

      node_to->minimum_energy = energy_to;
      node_to->saddle = saddle_energy;
      node_to->branch_length = saddle_energy - energy_to;
      node_to->minimum_id = id_to;

      node_t *parent = new node_t();
      parent->minimum_id = 0;
      parent->minimum_energy = saddle_energy;
      parent->saddle = saddle_energy;
      parent->parent = NULL;
      parent->left_child = node_from;
      parent->right_child = node_to;

      node_from->parent = parent;
      node_to->parent = parent;

      clusters.push_back(parent);
      std::unordered_set<size_t> *id_set_cluster =
          new std::unordered_set<size_t>();
      id_set_cluster->insert(id_from);
      id_set_cluster->insert(id_to);
      clusters_id_sets.push_back(id_set_cluster);

      clustered_ids.insert(id_from);
      clustered_ids.insert(id_to);

      //std::cout << "none in cluster " << id_from << " " << id_to << std::endl;
      continue;
    }
    if ((it_id_from != clustered_ids.end())
        && (it_id_to == clustered_ids.end())) {
      // create cluster for id to and merge it to cluster with id_from
      size_t c_from_id = 0; //= None
      size_t j = 0;
      for (auto it_clust_id = clusters_id_sets.begin();
          it_clust_id != clusters_id_sets.end(); it_clust_id++, j++) {
        auto it_contains_from = (*it_clust_id)->find(id_from);
        if (it_contains_from != (*it_clust_id)->end()) {
          c_from_id = j;
          break;
        }
      }
      node_t *node_from = clusters[c_from_id];
      clusters.erase(clusters.begin() + c_from_id);

      int energy_to = structure_index_to_energy[id_to];
      node_t *node_to = new node_t();
      node_to->minimum_id = id_to;
      node_to->minimum_energy = energy_to;
      node_to->saddle = saddle_energy;

      node_from->branch_length = saddle_energy - node_from->minimum_energy;
      node_to->branch_length = saddle_energy - energy_to;

      //saddle_energy = None
      node_t *parent = new node_t();
      parent->minimum_id = 0;
      parent->minimum_energy = saddle_energy;
      parent->saddle = saddle_energy;
      parent->parent = NULL;
      parent->left_child = node_from;
      parent->right_child = node_to;

      node_from->parent = parent;
      node_to->parent = parent;

      clusters.push_back(parent);
      std::unordered_set<size_t> *c_set_from_union_to = new std::unordered_set<size_t>();
      c_set_from_union_to->insert(clusters_id_sets[c_from_id]->begin(), clusters_id_sets[c_from_id]->end());
      delete clusters_id_sets[c_from_id];
      clusters_id_sets.erase(clusters_id_sets.begin() + c_from_id);
      c_set_from_union_to->insert(id_to);
      clusters_id_sets.push_back(c_set_from_union_to);

      clustered_ids.insert(id_to);

      //std::cout << "from in cluster " << id_from << " " << id_to << std::endl;
      continue;
    }
    if ((it_id_from == clustered_ids.end())
        && (it_id_to != clustered_ids.end())) {
      // create cluster for id_from and merge it to id_to
      size_t c_to_id;
      size_t k = 0;
      for (auto it_clust_id = clusters_id_sets.begin();
          it_clust_id != clusters_id_sets.end(); it_clust_id++, k++) {
        auto it_contains_from = (*it_clust_id)->find(id_to);
        if (it_contains_from != (*it_clust_id)->end()) {
          c_to_id = k;
          break;
        }
      }
      node_t *node_to = clusters[c_to_id];
      clusters.erase(clusters.begin() + c_to_id);

      int energy_from = structure_index_to_energy[id_from];
      node_t *node_from = new node_t();
      node_from->minimum_id = id_from;
      node_from->minimum_energy = energy_from;
      node_from->saddle = saddle_energy;

      node_from->branch_length = saddle_energy - energy_from;
      node_to->branch_length = saddle_energy - node_to->minimum_energy;

      node_t *parent = new node_t();
      parent->minimum_id = 0;
      parent->minimum_energy = saddle_energy;
      parent->saddle = saddle_energy;
      parent->parent = NULL;
      parent->left_child = node_from;
      parent->right_child = node_to;

      node_from->parent = parent;
      node_to->parent = parent;

      clusters.push_back(parent);
      std::unordered_set<size_t> *c_set_to_union_from = new std::unordered_set<size_t>();
      c_set_to_union_from->insert(clusters_id_sets[c_to_id]->begin(), clusters_id_sets[c_to_id]->end());
      delete clusters_id_sets[c_to_id];
      clusters_id_sets.erase(clusters_id_sets.begin() + c_to_id);
      c_set_to_union_from->insert(id_from);
      clusters_id_sets.push_back(c_set_to_union_from);

      clustered_ids.insert(id_from);

      //std::cout << "to in cluster " << id_from << " " << id_to << std::endl;
      continue;
    }
  }
  for(size_t i = 0; i < clusters_id_sets.size(); i++){
    delete clusters_id_sets[i];
  }
  clusters_id_sets.clear();

  return clusters;
}

/**
 * ! determine the min_h criterion.
 * @param maximal_number_of_states - filter criterion of maximal output states.
 * @param minimal_saddle_list - the list of minimal saddles for each state.
 * @return the minimal min_h that fulfills the maximal state criterion.
 */
double BarriersTree::determin_optimal_min_h(const size_t maximal_number_of_states, const std::vector<saddle_t> &minimal_saddle_list,
    const std::vector<std::pair<size_t, MyState *>> &sortedMinimaIDs){
  double min_h = 0;
  int* barrier_list = new int[minimal_saddle_list.size()];
  size_t n =  minimal_saddle_list.size();
  for(size_t i = 0; i < n; i++){
    barrier_list[i] = minimal_saddle_list[i].saddle - sortedMinimaIDs[minimal_saddle_list[i].minimum_from].second->energy;
  }
  // sort from highest to lowest barrier.
  std::sort(barrier_list, barrier_list + n, std::greater<int>());

  if(maximal_number_of_states > n){
    min_h = (double)barrier_list[n-1] / 100.0;
  }
  else if (maximal_number_of_states <= 2){
    min_h = (double)barrier_list[0] / 100.0;
    if (maximal_number_of_states <= 1)
      min_h += 0.01;
    /* The minimal saddle list contains all states, except the mfe (because there is no deeper basin).
       For the barrier at 0 we would connect with the mfe and would end up with 2 states.
       Thus we choose the barrier a little bit higher. The minimum number of states is 1, because there
       is not edge from the mfe.
    */
  }
  else{ // if (maximal_number_of_states >= 2 && maximal_number_of_states < n)
    // if the minh is equal to the next state, we have to determine the next higher minh
    // (because we will remove all states < minh, states >= minh will remain).
    size_t max_index = maximal_number_of_states-2; // -2 because: 1 for the mfe (because it is not in the minimal saddle list size) and 1 for the zero based index.
    if (barrier_list[max_index] == barrier_list[max_index+1]){
      while(max_index > 0 && barrier_list[max_index] == barrier_list[max_index-1]){
        max_index--;
      }
    }
    min_h = (double)barrier_list[max_index-1] / 100.0;
  }
  delete[] barrier_list;
  return min_h;
}

/**
 * ! minimum saddle height filter.
 * @param minimal_saddle_list - sorted list of minimal output transitions. Output is the filtered list. Used for input and output!
 * @param sortedMinimaIDs - map state id to local minimum.
 * @param min_h - the minimum height, all states with a lower barrier will be merged into a deeper basin according to the barrier tree.
 * @return map of minima that are merged together.
 */
std::unordered_map<size_t, size_t> BarriersTree::filter_minh(std::vector<saddle_t> &minimal_saddle_list, const std::vector<std::pair<size_t, MyState *>> &sortedMinimaIDs, const double min_h){
  std::unordered_map<size_t, size_t> merge_state_map;
  std::unordered_map<size_t, size_t>index_map;
  size_t number_result_states = minimal_saddle_list.size();
  double barrier_height;
  size_t state_index, neighbor_index;
  size_t i;
  for(size_t ii = minimal_saddle_list.size(); ii > 0; ii--){
    i = ii-1;
    state_index = minimal_saddle_list[i].minimum_from;
    barrier_height = (double)(minimal_saddle_list[i].saddle - sortedMinimaIDs[state_index].second->energy)/100.0;
    if (barrier_height < min_h){
      neighbor_index = minimal_saddle_list[i].minimum_to;
      merge_state_map[state_index] = neighbor_index;
      minimal_saddle_list.erase(minimal_saddle_list.begin() + i);
      number_result_states--;
      continue;
    }
    //index_map[state_index] = i; // index according to sorting by saddles.
  }
  number_result_states = minimal_saddle_list.size();

  for(size_t i =0; i < minimal_saddle_list.size(); i++){
    size_t state_index = minimal_saddle_list[i].minimum_from;
    index_map[state_index] = i; // index according to sorting by saddles.
  }
  return merge_state_map;
}


std::string BarriersTree::newick_string_builder(node_t *tree) {
  std::string sub_tree_string = "(";
  if (tree->left_child != NULL || tree->right_child != NULL) {
    sub_tree_string += "(";
    if (tree->left_child != NULL) {
      sub_tree_string += this->newick_string_builder(tree->left_child);
      sub_tree_string += ",";
    }
    if (tree->right_child != NULL) {
      sub_tree_string += this->newick_string_builder(tree->right_child);
      sub_tree_string += ",";
    }
    if (sub_tree_string[sub_tree_string.length() - 1] == ',') {
      sub_tree_string = sub_tree_string.substr(0, sub_tree_string.length() - 1);
    }
    sub_tree_string += ")";
  } else {
    // it is a leaf node --> add index
    sub_tree_string += std::to_string(tree->minimum_id);
  }
  if (tree->parent == NULL) {
    sub_tree_string.append(":").append(std::to_string(0)).append(")");
  } else {
    sub_tree_string.append(":").append(std::to_string(tree->branch_length)).append(
        ")");
  }
  return sub_tree_string;
}

/**
 * Create a vector graphic of the tree.
 * It does an inorder traversal of the tree in order to determine the x-coordinates.
 * @param tree - the barriers tree.
 * @param mfe - the minimum free energy in dcal/mol
 * @param inorder_index - a pointer to an index variable.
 * @return - returns the svg as string.
 */
std::string BarriersTree::svg_string_builder(node_t *tree, int mfe, size_t* inorder_index) {
  // traverse tree from lowest left to right.
  float node_scale = 20;
  float legend_margin = 50;
  std::string sub_tree_string = ""; //"(";
  if (tree->left_child != NULL || tree->right_child != NULL) {
    size_t left_idx, current_idx, right_idx;
    //sub_tree_string += "(";
    if (tree->left_child != NULL) {
      sub_tree_string += this->svg_string_builder(tree->left_child, mfe, inorder_index);
      sub_tree_string += "";
    }
    left_idx = tree->left_child->node_index;
    *inorder_index += 1;
    current_idx = *inorder_index;
    tree->node_index = current_idx;
    // draw root
    if (tree->right_child != NULL) {
      sub_tree_string += this->svg_string_builder(tree->right_child, mfe, inorder_index);
      sub_tree_string += "";
    }
    right_idx = tree->right_child->node_index;

    // line from current to parent
    sub_tree_string += "<path class=\"barrier\" d=\"M " + std::to_string(legend_margin + (float)tree->node_index * node_scale) + " "+ std::to_string(-tree->saddle);
    sub_tree_string += " l "+ std::to_string(0) + " "+ std::to_string(-tree->branch_length) + "\" stroke-width=\"3\" fill=\"none\"/>\n";

    // line from left to right.
    sub_tree_string += "<path class=\"saddle\" d=\"M " + std::to_string(legend_margin +(float)tree->left_child->node_index * node_scale) + " "+ std::to_string(-tree->saddle);
    sub_tree_string += " l "+ std::to_string(((float)tree->right_child->node_index - (float)tree->left_child->node_index) * node_scale) + " "+ std::to_string(0) + " \" stroke-width=\"3\" fill=\"none\"/>\n";

    // branch label
    char branch_label[50];
    sprintf(branch_label,"%.2f",tree->branch_length/100.0);
    float x = legend_margin + (float)tree->node_index * node_scale;
    float y = (-tree->saddle -tree->branch_length/2.0);
    sub_tree_string += "<text x=\""+std::to_string(x - 15)+"\" y=\""+std::to_string(y - 5)+
                  "\" transform=\"rotate(-90 "+std::to_string(x)+","+std::to_string(y)+")\">"+branch_label+"</text>\n";
  } else {
    // it is a leaf node --> add index
    *inorder_index += 1;
    tree->node_index = *inorder_index;
    sub_tree_string += "<text class=\"leaf_label\" x=\""+std::to_string(legend_margin+ (float)*inorder_index * node_scale - 5)+"\" y=\""+std::to_string(-tree->minimum_energy + 20)+
        "\" fill=\"black\">"+std::to_string(tree->minimum_id)+"</text>\n";

    // line from node to saddle
    sub_tree_string += "<path class=\"barrier\" d=\"M " + std::to_string(legend_margin + *inorder_index * node_scale) + " "+ std::to_string(-tree->minimum_energy)+
                                     " l "+ std::to_string(0) + " "+ std::to_string(-tree->branch_length) + "\" stroke-width=\"3\" fill=\"none\"/>\n";
    // branch label
    char branch_label[50];
    sprintf(branch_label,"%.2f",tree->branch_length/100.0);
    float x = legend_margin + (float)tree->node_index * node_scale;
    float y = (-tree->minimum_energy  - tree->branch_length/2.0);
    sub_tree_string += "<text class=\"barrier_label\" x=\""+std::to_string(x - 15)+"\" y=\""+std::to_string(y - 5)+
                  "\" transform=\"rotate(-90 "+std::to_string(x)+","+std::to_string(y)+")\">"+branch_label+"</text>\n";
  }
  if (tree->parent == NULL) {
    // create a viewbox, that allows negative indices.
    int upper_tree_bound = std::ceil((tree->saddle)/100.0) *100.0;
    int lower_tree_bound = std::floor((mfe)/100.0)*100.0;
    int tree_height = upper_tree_bound - lower_tree_bound;
    double double_border = 100.0;
    sub_tree_string = "<svg height=\""+std::to_string(tree_height+double_border)+"\" width=\""+std::to_string(legend_margin+ (float)*inorder_index * node_scale + double_border)+ "\"" +
        "  viewBox=\""+std::to_string(-double_border/2)+" "+std::to_string(-upper_tree_bound - double_border/2)+" "+std::to_string(legend_margin +(float)*inorder_index * node_scale+ double_border)+" "+std::to_string(tree_height+ double_border)+"\">\n"
        + "<g class=\"tree\" stroke=\"black\">\n" +sub_tree_string + "</g>\n";
    //add ruler.
    double ruler_offset = -upper_tree_bound;
    sub_tree_string +=  "<g class=\"ruler\">\n";
    sub_tree_string += "<text x=\""+std::to_string(0)+"\" y=\""+std::to_string(ruler_offset - 20)+
                  "\" fill=\"black\">Energy [kcal/mol]</text>\n";
    sub_tree_string += "<line x1=\""+std::to_string(0)+"\" y1=\""+std::to_string(ruler_offset)+
        "\" x2=\""+std::to_string(0)+"\" y2=\""+std::to_string(-lower_tree_bound)+"\" style=\"stroke:rgb(0,0,0);stroke-width:2\"/>\n";
    // add ticks
    for (int i = -upper_tree_bound; i <= -lower_tree_bound; i+= 100){
      sub_tree_string += "<line x1=\""+std::to_string(0)+"\" y1=\""+std::to_string(i)+
              "\" x2=\""+std::to_string(1*(legend_margin/10))+"\" y2=\""+std::to_string(i)+"\" style=\"stroke:rgb(0,0,0);stroke-width:2\"/>\n";
      // add tick label
      char tick_label[50];
      sprintf(tick_label,"%.2f",-i/100.0);
      sub_tree_string += "<text x=\""+std::to_string(1*(legend_margin/10) + 10)+"\" y=\""+std::to_string(i + 3)+
              "\" fill=\"black\">"+tick_label+"</text>\n";
    }
    sub_tree_string += "</g>\n";
    sub_tree_string += " </svg>\n";
  }
  return sub_tree_string;
}

void BarriersTree::free_tree(node_t* tree){
  if (tree->left_child != NULL || tree->right_child != NULL) {
    if (tree->left_child != NULL) {
      this->free_tree(tree->left_child);
    }
    if (tree->right_child != NULL) {
      this->free_tree(tree->right_child);
    }
  }
  delete tree;
}


