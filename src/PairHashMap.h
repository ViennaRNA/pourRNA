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
typedef Concurrent_Pair_Hash_Map<std::pair<MyState,
                                           MyState>, MyState, PairMapHash, PairMapEqual> HashMap;
};


#endif /* PAIRHASHMAP_H_ */
