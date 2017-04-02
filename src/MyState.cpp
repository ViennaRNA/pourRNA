/*
 * MyState.cpp
 *
 *  Created on: 18.08.2014
 *      Author: Gregor Entzian
 */

#include "MyState.h"

MyState::MyState () :
energy (0), structure (NULL)
{
}

MyState::MyState (int energy, Structure structure) :
    		energy(energy)
{
	this->structure=NULL;
	if (structure != NULL)
	{
		this->structure = (allocopy (structure));
	}
}

MyState::MyState (struct_en* structureAndEnergy) :
    		energy (structureAndEnergy->energy)
{
	this->structure=NULL;
	if (structureAndEnergy->structure != NULL)
	{
		this->structure = (allocopy (structureAndEnergy->structure));
	}
}

MyState::MyState (const MyState& state) :
    		energy (state.energy)
{
	this->structure=NULL;
	if (state.structure != NULL)
	{
		this->structure = (allocopy (state.structure));
	}
}

MyState*
MyState::clone () const
{
	MyState* newState = new MyState ();
	newState->energy = this->energy;
	newState->structure = allocopy (this->structure);
	return newState;
}

std::string
MyState::toString () const
{
	if (structure == NULL)
	{
		// std::cout << "pair_table is empty!" << std::endl;
		return "";
	}
	std::string s (structure[0], ' ');
	for (int i = 1; i <= structure[0]; i++)
	{
		if (structure[i] == 0)
			s[i-1] = '.';
		else if (structure[i] < i)
			s[i-1] = ')';
		else
			s[i-1] = '(';
	}
	return s;
}

MyState::~MyState ()
{
	energy = 0;
	free (structure);
	structure = NULL;
}

