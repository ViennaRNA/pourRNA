/*
 * SC_PartitionFunction.h
 *
 *  Created on: 01.08.2014
 *      Author: Gregor Entzian
 */

#ifndef SC_PARTITIONFUNCTION_H_
#define SC_PARTITIONFUNCTION_H_

#include <cmath>
#include <unordered_map>
#include <vector>
#include "StateCollector.h"

/**
 * ! This state collector updates the partition function if a state is added.
 */
class SC_PartitionFunction : public StateCollector
{
public:
/*! constructor
 * @param temperature the temperature to be used to compute Boltzmann weights
 * @param storeEnergies whether or not to store a list of all energies added
 * to the container
 */
SC_PartitionFunction(double temperature = 37.0,
                     double gas_constant = 0.00198717,
                     double mfe = 0,
                     bool   storeEnergies = false);

/*! sets the partition function sum to 0 and calculate the gasconstant (RT)
 * with the given temperature in degrees Celsius.
 * @param temperature the temperature in Celsius.
 * @param storeEnergies decides if the energies of all states will be stored in a vector.
 */
virtual
void initialize(const double  temperature,
                const double  gas_constant,
                double mfe,
                const bool    storeEnergies = false);


/*!
 * Adds the Boltzmann weight of the given state to the partition function.
 * @param state the state to add
 */
virtual void
add(const MyState& state);


virtual
~SC_PartitionFunction ();

/*! Access to the partition function.
 * @return the partition function over all states reported to this collector
 */
virtual double
getZ() const
{
  return Z;
}


void
setZ(double z = 0)
{
  Z = z;
}


/*!
 * Computes the Boltzmann weight of the given state based on its energy and
 * the currently used temperature. It is scaled by the minimum free energy!
 * @param state the state of interest
 * @return the Boltzmann weight of the state
 */
inline
double
getBoltzmannWeight(const MyState & state) const
{
  // NOTE: need to divide energy from vienna package by 100 to get kcal/mol
  //(calculate with float energies to be consistent with RNA_Explore from RNAKinetics-Package.
  return std::exp((MFE -(double)(state.energy) / 100.0) / kT);
}


const std::vector<int>&
getEnergies() const
{
  return Energies;
}


/*!
 * index pair
 */
typedef std::pair<size_t, size_t> PairID;


class UInt_Pair_Hash {
public:
  std::uint64_t
  operator()(const std::pair<std::uint32_t, std::uint32_t>& k) const
  {
    std::uint64_t first = k.first;
    std::uint64_t second = k.second;

    first = first << 32;
    first += second;
    return first;
  }
};
/*!
 * container for all partition functions to generate
 * indices are according to minima container
 * Z(i,i) will hold the basins partition function
 * Z(i,j) will hold the partition function of the saddle points between minimum i and j
 */
typedef std::unordered_map<PairID, SC_PartitionFunction, UInt_Pair_Hash> Z_Matrix;

protected:
//! variable decides if energies should be stored.
bool              StoreEnergies;
//! stores the energies of the basin if requested.
std::vector<int>  Energies;
//! The partition function sum to fill;
double            Z;
//! Gas constant unit : kcal / (Kelvin * mol)
const double      GAS_CONSTANT_KCAL;
//! The temperature scaled "Boltzmann" constant used to compute Boltzmann weights.
double            kT;
//! The global minimum free energy -> used for scaling
double            MFE;
};


#endif /* SC_PARTITIONFUNCTION_H_ */
