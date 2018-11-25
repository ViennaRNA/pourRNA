/*
 * MyState.h
 *
 *  Created on: 18.08.2014
 *      Author: Gregor Entzian
 */

#ifndef MYSTATE_H_
#define MYSTATE_H_

#include <cstdio>
#include <stdlib.h>
#include <iostream>
extern "C" {
#include <ViennaRNA/move_set.h>
}
#include <string>
#include "StructureUtils.h"

/**
 * ! wrapper class for struct_en
 * with clone function and
 * destructor to clean up.
 */
class MyState {
public:

	//! type of structure representation
	typedef short* Structure;

	MyState();
	MyState(int energy, short* structure);
	MyState(struct_en* state);
	MyState(const MyState& state);
	//! less comparison based on E and using the lexicographical order
	//! of s as tiebreaker
	bool
	operator<(const MyState& toCompare) const;
	MyState
	operator=(MyState toCopy);
	MyState*
	clone() const;
	std::string
	toString() const;
	~MyState();

	int getEnergy() const {
		return energy;
	}

	void setEnergy(int energy) {
		this->energy = energy;
	}

	const Structure &
	getStructure() const {
		return structure;
	}

	void setStructure(Structure structure) {
		this->structure = structure;
	}

public:
	// ! energy in 10kcal/mol
	int energy;
	// ! structure in pairTable format
	Structure structure;

};

#endif /* MYSTATE_H_ */
