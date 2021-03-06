/*
 * SC_PartitionFunction.cpp
 *
 *  Created on: 01.08.2014
 *      Author: Gregor Entzian
 */

#include "SC_PartitionFunction.h"

SC_PartitionFunction::SC_PartitionFunction (const double  temperature,
                                            const double  gas_constant,
                                            const bool    storeEnergies)
  :
  StateCollector(),
  StoreEnergies(storeEnergies),
  Energies(),
  Z(0.0),
  GAS_CONSTANT_KCAL(gas_constant),
  kT(GAS_CONSTANT_KCAL * (temperature + 273.15))
{
  // initialize data structures
  //initialize(temperature, storeEnergies);
}


void
SC_PartitionFunction::add(const MyState& state)
{
  StateCollector::add(state);
  Z += getBoltzmannWeight(state);
  if (StoreEnergies)
    Energies.push_back(state.energy);
}


SC_PartitionFunction::~SC_PartitionFunction ()
{
}


void
SC_PartitionFunction::initialize(const double temperature,
                                 const double gas_constant,
                                 const bool   storeEnergies)
{
  // update constant parameters
  StoreEnergies = storeEnergies;
  double tempInKelvin = temperature + 273.15;
  kT = GAS_CONSTANT_KCAL * tempInKelvin;

  // init data structures
  Z = 0.0;
  Energies.clear();
}
