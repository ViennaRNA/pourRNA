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
                                            bool   storeStructures,
                                            const bool    storeEnergies)
  :
  StateCollector(),
  StoreEnergies(storeEnergies),
  StoreStructures(storeStructures),
  Energies(),
  Structures(),
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
  if (StoreStructures)
    Structures.push_back(*state.clone());
}


SC_PartitionFunction::~SC_PartitionFunction ()
{
}


void
SC_PartitionFunction::initialize(const double temperature,
                                 const double gas_constant,
                                 double mfe,
                                 bool   storeStructures,
                                 const bool   storeEnergies)
{
  // update constant parameters
  StoreEnergies = storeEnergies;
  StoreStructures = storeStructures;
  double tempInKelvin = temperature + 273.15;
  kT = GAS_CONSTANT_KCAL * tempInKelvin;

  // init data structures
  Z = 0.0;
  Energies.clear();
  Structures.clear();
  MFE = mfe;
  Temperature = temperature;
}

SC_PartitionFunction&
SC_PartitionFunction::operator=(const SC_PartitionFunction& toCopy){
  this->initialize(toCopy.getTemperature(), toCopy.getGasConstant(), toCopy.getMFE(), toCopy.getStoreStructures(), toCopy.getStoreEnergies());
  const std::vector<int>& energies = toCopy.getEnergies();
  this->Energies.assign(energies.begin(), energies.end());
  const std::vector<MyState>& structures = toCopy.getStructures();
  this->Structures.assign(structures.begin(), structures.end());
  this->setZ(toCopy.getZ());
  return *this;
}

SC_PartitionFunction&
SC_PartitionFunction::operator+=(const SC_PartitionFunction& toAdd){
  bool allSane = true;
  if (this->getTemperature() != toAdd.getTemperature()){
    fprintf(stderr, "Error: could not add partition function! Different temperature!");
    allSane = false;
  }
  if (this->getGasConstant() != toAdd.getGasConstant()){
    fprintf(stderr, "Error: could not add partition function! Different gas constant!");
    allSane = false;
  }
  if (this->getMFE() != toAdd.getMFE()){
    fprintf(stderr, "Error: could not add partition function! Different mfe!");
    allSane = false;
  }
  if (this->getStoreStructures() != toAdd.getStoreStructures()){
    fprintf(stderr, "Error: could not add partition function! Different store structure options!");
    allSane = false;
  }
  if (this->getStoreEnergies() != toAdd.getStoreEnergies()){
    fprintf(stderr, "Error: could not add partition function! Different store energy options!");
    allSane = false;
  }
  if (!allSane)
    std::exit(EXIT_FAILURE);
  const std::vector<int>& energies = toAdd.getEnergies();
  this->Energies.assign(energies.begin(), energies.end());
  const std::vector<MyState>& structures = toAdd.getStructures();
  this->Structures.assign(structures.begin(), structures.end());
  double to_Z = this->getZ();
  to_Z += toAdd.getZ();
  this->setZ(to_Z);
  return *this;
}



