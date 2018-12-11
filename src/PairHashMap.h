/*
 * PairHashMap.h
 *
 *  Created on: 24.06.2017
 *      Author: Gregor Entzian
 */

#ifndef PAIRHASHMAP_H_
#define PAIRHASHMAP_H_

#include <ViennaRNA/move_set.h>
#include <unordered_map>
#include <unordered_set>

#include <thread>
#include <mutex>
#include <condition_variable>

#include <string>
#include <functional>
#include "StructureUtils.h"
#include "MyState.h"
#include "SpookyHash/SpookyV2.h"
#include "Concurrent_Pair_Hash_Map.h"


/**
 * This class maps a rna-structure (in pairTableFormat (see ViennaRNA-Package)) to an integer.
 */
class PairHashMap
{
public:

 /*
struct PairMapHash {
  std::uint64_t
  operator()(const std::pair<MyState, MyState>& k) const
  {
    // std::string keyStructure = k.toString();
    // return std::hash<std::string> () (keyStructure);
    std::uint64_t first =
      SpookyHash::Hash32(k.first.structure, sizeof(short) * (k.first.structure[0] + 1), 0);
    std::uint64_t second =
      SpookyHash::Hash32(k.second.structure, sizeof(short) * (k.second.structure[0] + 1), 0);

    first << 32;
    first += second;
    return first;
  }
};

struct PairMapEqual {
  bool
  operator()(const std::pair<MyState, MyState>& lhs,
             const std::pair<MyState, MyState>& rhs) const
  {
    return StructureUtils::IsEqual(lhs.first.structure, rhs.first.structure) &&
           StructureUtils::IsEqual(lhs.second.structure, rhs.second.structure);
  }
};
typedef Concurrent_Pair_Hash_Map<std::pair<MyState, MyState>, MyState, PairMapHash, PairMapEqual> HashMapStates;
*/


struct PairMapHash_uint32 {
  std::uint64_t
  operator()(const std::pair<std::uint32_t, std::uint32_t>& k) const
  {
    // std::string keyStructure = k.toString();
    // return std::hash<std::string> () (keyStructure);
    std::uint64_t first = k.first;
    std::uint64_t second = k.second;

    first << 32;
    first += second;
    return first;
  }
};

struct PairMapEqual_uint32 {
  bool
  operator()(const std::pair<std::uint32_t, std::uint32_t>& lhs,
             const std::pair<std::uint32_t, std::uint32_t>& rhs) const
  {
    return (lhs.first == rhs.first) &&
           (lhs.second == rhs.second);
  }
};

typedef Concurrent_Pair_Hash_Map<std::pair<std::uint32_t,
    std::uint32_t>, MyState, PairMapHash_uint32, PairMapEqual_uint32> HashMap;
};


#endif /* PAIRHASHMAP_H_ */
