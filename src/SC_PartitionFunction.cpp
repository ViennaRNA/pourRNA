/*
 * SC_PartitionFunction.cpp
 *
 *  Created on: 01.08.2014
 *      Author: Gregor Entzian
 */

#include "SC_PartitionFunction.h"

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
{
	// initialize data structures
	initialize(temperature, storeEnergies);
}

void
SC_PartitionFunction::add (const MyState& state)
{
	StateCollector::add (state);
	Z += getBoltzmannWeight( state );
	if(StoreEnergies)
	{
		Energies.push_back (state.energy);
	}
}

SC_PartitionFunction::~SC_PartitionFunction ()
{
}

void
SC_PartitionFunction::initialize (const double temperature, const bool storeEnergies)
{
	// update constant parameters
	StoreEnergies=storeEnergies;
	double tempInKelvin = temperature + 273.15;
	kT = GAS_CONSTANT_R * JOULE_TO_KCAL * tempInKelvin;

	// init data structures
	Z = 0.0;
	Energies.clear();
}
