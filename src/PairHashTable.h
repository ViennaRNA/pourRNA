/*
 * PairHashTable.h
 *
 *  Created on: 30.07.2014
 *      Author: Gregor Entzian
 */

#ifndef PAIRHASHTABLE_H_
#define PAIRHASHTABLE_H_

#include <ViennaRNA/move_set.h>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <functional>
#include "StructureUtils.h"
#include "MyState.h"
#include "SpookyHash/SpookyV2.h"
#include "Concurrent_Pair_Hash_Map.h"
/**
 * This class maps a rna-structure (in pairTableFormat (see ViennaRNA-Package)) to an integer.
 */
class PairHashTable
{
public:

struct PairTableHash {
  std::size_t
  operator()(const MyState& k) const
  {
    // std::string keyStructure = k.toString();
    // return std::hash<std::string> () (keyStructure);
    return SpookyHash::Hash64(k.structure, sizeof(short) * (k.structure[0] + 1), 0);
  }
};

struct PairTableEqual {
  bool
  operator()(const MyState& lhs,
             const MyState& rhs) const
  {
    return StructureUtils::IsEqual(lhs.structure, rhs.structure);
  }
};
typedef Concurrent_Pair_Hash_Map<MyState, std::uint32_t, PairTableHash, PairTableEqual> HashTable;
typedef Concurrent_Pair_Hash_Map<std::uint32_t, MyState, std::hash<std::uint32_t>, std::equal_to<std::uint32_t>> HashTable_Reverse;
};

class HashSet
{
public:

struct SetHash {
  std::size_t
  operator()(const MyState& k) const
  {
    // std::string keyStructure = k.toString();
    // return std::hash<std::string> () (keyStructure);
    return SpookyHash::Hash64(k.structure, sizeof(short) * (k.structure[0] + 1), 0);
  }
};

struct HashSetEqual {
  bool
  operator()(const MyState& lhs,
             const MyState& rhs) const
  {
    return StructureUtils::IsEqual(lhs.structure, rhs.structure);
  }
};
typedef std::unordered_set<MyState, SetHash, HashSetEqual> UnorderedHashSet;
};


#endif /* PAIRHASHTABLE_H_ */
