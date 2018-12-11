/*
 * TypeID.h
 *
 *  Created on: 14.09.2014
 *      Author: Gregor Entzian
 */

#ifndef TYPEID_H_
#define TYPEID_H_

#include <cstddef>
#include <cstdint>
using namespace std;
/**
 * ! Assign unique integer IDs.
 */
class TypeID
{
public:
static std::uint32_t counter;

public:
static std::uint32_t
value()
{
  return ++counter;
}
};

#endif /* TYPEID_H_ */
