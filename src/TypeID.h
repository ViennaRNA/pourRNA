/*
 * TypeID.h
 *
 *  Created on: 14.09.2014
 *      Author: Gregor Entzian
 */

#ifndef TYPEID_H_
#define TYPEID_H_

#include <cstddef>
using namespace std;
/**
 * ! Assign unique integer IDs.
 */
class TypeID
{
public:
static size_t counter;

public:
template<typename T>
static size_t
value()
{
  return ++counter;
}
};

#endif /* TYPEID_H_ */
