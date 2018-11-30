/*
 * SpookyHashMap.h
 *
 *  Created on: Nov 7, 2014
 *      Author: Gregor Entzian
 */

#ifndef SPOOKYHASHMAP_H_
#define SPOOKYHASHMAP_H_

#include "../StructureUtils.h"
#include <unordered_map>
#include "SpookyV2.h"
#include "../BIUlibPart/LimitedHash.hh"

class SpookyHashMap
{
public:

struct SpookyHashMapHash {
  std::size_t
  operator()(const MyState& k) const
  {
    return SpookyHash::Hash64(k.structure, sizeof(short) * (k.structure[0] + 1), 0);
  }
};

struct SpookyHashMapEqual {
  bool
  operator()(const MyState& lhs,
             const MyState& rhs) const
  {
    return StructureUtils::IsEqual(lhs.structure, rhs.structure);
  }
};
typedef biu::LimitedHash<MyState, MyState, SpookyHashMapHash, SpookyHashMapEqual> HashTable;
};

#endif /* SPOOKYHASHMAP_H_ */
