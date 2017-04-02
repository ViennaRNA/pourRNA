/*
 * SC_PartitionFunction.cpp
 *
 *  Created on: 01.08.2014
 *      Author: Gregor Entzian
 */

#include "SC_PartitionFunction.h"

<<<<<<< SC_PartitionFunction.cpp
SC_PartitionFunction::SC_PartitionFunction () :
    StateCollector (), StoreEnergies (false), Z (0.0), GAS_CONSTANT_RT (GAS_CONSTANT_R * JOULE_TO_KCAL * 310.15), kT (
	GAS_CONSTANT_RT)
=======
SC_PartitionFunction::SC_PartitionFunction (const double temperature, const bool storeEnergies)
 :
	StateCollector(),
	StoreEnergies(false),
	Energies(),
	Z (0.0),
	JOULE_TO_KCAL(0.000239),
	GAS_CONSTANT_R(8.314472),
	GAS_CONSTANT_RT(GAS_CONSTANT_R * JOULE_TO_KCAL * 310.15),
	kT (GAS_CONSTANT_RT)
>>>>>>> 1.9
{
	// initialize data structures
	initialize(temperature, storeEnergies);
}

void
SC_PartitionFunction::add (const MyState& state)
{
<<<<<<< SC_PartitionFunction.cpp
  StateCollector::add (state);
  //(calculate with float energies to be consistent with RNA_Explore from RNAKinetics-Package.
  Z += exp (-(float) (state.energy / 100.0) / kT);
  if (StoreEnergies)
    {
      Energies.push_back (state.energy);
    }
=======
	StateCollector::add (state);
	Z += getBoltzmannWeight( state );
	if(StoreEnergies)
	{
		Energies.push_back (state.energy);
	}
>>>>>>> 1.9
}

SC_PartitionFunction::~SC_PartitionFunction ()
{
}

void
SC_PartitionFunction::initialize (const double temperature, const bool storeEnergies)
{
<<<<<<< SC_PartitionFunction.cpp
  StoreEnergies = storeEnergies;
  double tempInKelvin = temperature + 273.15;
  kT = GAS_CONSTANT_R * JOULE_TO_KCAL * tempInKelvin;
=======
	// update constant parameters
	StoreEnergies=storeEnergies;
	double tempInKelvin = temperature + 273.15;
	kT = GAS_CONSTANT_R * JOULE_TO_KCAL * tempInKelvin;

	// init data structures
	Z = 0.0;
	Energies.clear();
>>>>>>> 1.9
}
