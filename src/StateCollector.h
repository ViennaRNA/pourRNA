/*
 * StateCollector.h
 *
 *  Created on: 26.07.2014
 *      Author: Gregor Entzian
 */

#ifndef STATECOLLECTOR_H_
#define STATECOLLECTOR_H_

#include "MyState.h"
#include <stdint.h>

/**
 * ! This is a state collector which only counts the added states.
 */
class StateCollector
{
protected:
//! the number of added States
size_t  StateCount;
//! the energy of the state with maximum energy in the basin.
int     MaxEnergy;

public:
StateCollector ();
virtual
~StateCollector ();
/**
 * ! count this state.
 * @param state the state which consists of structure and energy.
 */
virtual void
add(const MyState& state);


//! Returns number of added States.
//! @return the number of added states
virtual size_t
size() const;


/**
 * ! return the energy of the state with maximum energy in the basin in kcal/mol.
 */
double
getMaxEnergy()
{
  return (double)MaxEnergy / 100.0;
}
};

#endif /* STATECOLLECTOR_H_ */
