/*
 * SC_DotPlot.h
 *
 *  Created on: 02.04.2015
 *      Author: Martin Mann
 */

#ifndef SC_DOTPLOT_H_
#define SC_DOTPLOT_H_

#include "SC_PartitionFunction.h"

extern "C" {
#include <ViennaRNA/utils.h>
#include <ViennaRNA/PS_dot.h>
}

/**
 * ! This state collector updates the base pair probabilities and the partition
 * function for the given ensemble if a state is added.
 */
class SC_DotPlot : public SC_PartitionFunction
{
public:

	using SC_PartitionFunction::PairID;

	//! container to store base pair weight sums or probabilities
	typedef std::map< PairID, double > DotPlot;

public:

	/*! construction
	 * @param temperature the temperature to be used to compute Boltzmann weights
	 * @param storeEnergies whether or not to store a list of all energies added
	 * to the container
	 */
	SC_DotPlot (const double temperature=37.0, const bool storeEnergies=false);

	/*! resets the data structure and calculate the gasconstant (RT)
	 * with the given temperature in degrees Celsius.
	 * @param temperature the temperature in Celsius.
	 * @param storeEnergies decides if the energies of all states will be stored in a vector.
	 */
	virtual
	void
	initialize(const double temperature, const bool storeEnergies=false);

	/*!
	 * updates the overall partition function as well as the probabilities of
	 * all base pair part of the current state
	 * @param state the state to add
	 */
	virtual
	void
	add (const MyState& state);

	/*! destruction
	 */
	virtual
	~SC_DotPlot ();

	/*! Access to the partition functions of all structures per base pair. To
	 * derive base pair probabilities, the values have to be divided by the
	 * overall partition function accessible via getZ().
	 * @return the base pair specific partition functions
	 */
	const DotPlot&
	getBasePairWeightSum() const {
		return bpWeightSum;
	}

	/*! Computes the base pair probabilities.
	 * @return the probability for each base pair part of a structure reported
	 * to the container.
	 */
	DotPlot
	getBasePairProbabilities() const;


	/*! Computes the base pair probabilities based on a given set of base pair
	 * partition functions and the overall partition function. Each base pair
	 * specific partition function is divided by the overall partition function
	 * to derive the base pair probability.
	 *
	 * @param bpWeightSum the partition function for each base pair
	 * @param Z the overall partition function to derive base pair probabilities
	 *
	 * @return the probability for each base pair part of the input
	 */
	static
	DotPlot
	getBasePairProbabilities(const DotPlot& bpWeightSums, const double Z);


	/*! Writes the dot plot in vienna postscript format to file. Note, no mfe
	 * structure plot is made, only base pair probabilities (upper triangle)
	 * are plotted.
	 * @param absoluteFileName the absolute file name to create and write the
	 *           dot plot to
	 * @param sequence the sequence of the RNA molecule the dot plot is about
	 * @param dotPlot the dot plot base pair probabilities to print
	 * @return true if the write was successful; false otherwise
	 */
	static
	bool
	writeDotPlot_PS(
					const std::string & absoluteFileName,
					const std::string & sequence,
					const DotPlot & dotPlot
					);

	/*! Reads the dot plot information from a postscript file produced by the
	 * vienna rna package.
	 * @param input stream of the dot plot file
	 * @return the base pair probabilities read from the dot plot file (ubox)
	 */
	static
	DotPlot
	readDotPlot_PS( std::istream & input );

protected:

	//! the partition function of all structures that contain a specific base pair
	DotPlot bpWeightSum;

};

#endif /* SC_DOTPLOT_H_ */

