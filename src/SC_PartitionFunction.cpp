/*
 * SC_PartitionFunction.cpp
 *
 *  Created on: 01.08.2014
 *      Author: Gregor Entzian
 */

#include "SC_PartitionFunction.h"

SC_PartitionFunction::SC_PartitionFunction (const double  temperature,
                                            const double  gas_constant,
                                            double        mfe,
                                            const bool    storeEnergies)
  :
  StateCollector(),
  StoreEnergies(storeEnergies),
  Energies(),
  Z(0.0),
  GAS_CONSTANT_KCAL(gas_constant),
  MFE(mfe),
  kT(GAS_CONSTANT_KCAL * (temperature + 273.15)),
  Temperature(temperature)
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
                                 double mfe,
                                 const bool   storeEnergies)
{
  // update constant parameters
  StoreEnergies = storeEnergies;
  double tempInKelvin = temperature + 273.15;
  kT = GAS_CONSTANT_KCAL * tempInKelvin;

  // init data structures
  Z = 0.0;
  Energies.clear();
  MFE = mfe;
  Temperature = temperature;
}

SC_PartitionFunction&
SC_PartitionFunction::operator=(const SC_PartitionFunction& toCopy){
  this->initialize(toCopy.getTemperature(), toCopy.getGasConstant(), toCopy.getMFE(), toCopy.getStoreEnergies());
  const std::vector<int>& energies = toCopy.getEnergies();
  this->Energies.assign(energies.begin(), energies.end());
  this->setZ(toCopy.getZ());
  return *this;
}



