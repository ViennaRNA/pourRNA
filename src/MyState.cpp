/*
 * MyState.cpp
 *
 *  Created on: 18.08.2014
 *      Author: Gregor Entzian
 */


extern "C" {
//#include <ViennaRNA/utils/structures.h>
#include <ViennaRNA/data_structures.h>
}
#include "MyState.h"

MyState::MyState () :
  energy(0), structure(NULL)
{
}


MyState::MyState (int       energy,
                  Structure structure) :
  energy(energy)
{
  this->structure = NULL;
  if (structure != NULL)
    this->structure = vrna_ptable_copy(structure);
}


MyState::MyState (struct_en *structureAndEnergy) :
  energy(structureAndEnergy->energy)
{
  this->structure = NULL;
  if (structureAndEnergy->structure != NULL)
    this->structure = (vrna_ptable_copy(structureAndEnergy->structure));
}


MyState::MyState (const MyState& state) :
  energy(state.energy)
{
  this->structure = NULL;
  if (state.structure != NULL)
    this->structure = (vrna_ptable_copy(state.structure));
}


bool
MyState::operator <(const MyState& s) const
{
  if (this->energy < s.energy)
    return true;

  if (this->energy == s.energy
      && StructureUtils::IsSmaller(this->structure, s.structure))  // use lexicographic order as tiebreaker
    return true;

  return false;
}


MyState
MyState::operator=(MyState toCopy)
{
  this->structure = NULL;
  if (toCopy.structure != NULL)
    this->structure = (vrna_ptable_copy(toCopy.structure));

  this->energy = toCopy.energy;
  return *this;
}


MyState *
MyState::clone() const
{
  MyState *newState = new MyState();

  newState->energy    = this->energy;
  newState->structure = vrna_ptable_copy(this->structure);
  return newState;
}


std::string
MyState::toString() const
{
  if (structure == NULL)
    // std::cout << "pair_table is empty!" << std::endl;
    return "";

  std::string s(structure[0], ' ');
  for (int i = 1; i <= structure[0]; i++) {
    if (structure[i] == 0)
      s[i - 1] = '.';
    else if (structure[i] < i)
      s[i - 1] = ')';
    else
      s[i - 1] = '(';
  }
  return s;
}


MyState::~MyState ()
{
  energy = 0;
  free(structure);
  structure = NULL;
}
