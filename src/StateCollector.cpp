/*
 * StateCollector.cpp
 *
 *  Created on: 26.07.2014
 *      Author: Gregor Entzian
 */

#include "StateCollector.h"

StateCollector::StateCollector () :
    StateCount (0), MaxEnergy (INT32_MIN)
{
}

StateCollector::~StateCollector ()
{
}

void
StateCollector::add (const MyState& state)
{
  StateCount++;
  if (state.energy > MaxEnergy)
    {
      MaxEnergy = state.energy;
    }
}

size_t
StateCollector::size () const
{
  return StateCount;
}

