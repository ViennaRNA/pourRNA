/*
 * StructureUtils.h
 *
 *  Created on: 05.08.2014
 *      Author: Gregor Entzian
 */

#ifndef STRUCTUREUTILS_H_
#define STRUCTUREUTILS_H_

#include <string>
#include <iostream>
#include <regex>
/**
 * ! Purpose: convert structures or compare them.
 */
class StructureUtils
{
public:
StructureUtils ();
virtual
~StructureUtils ();
/**
 * ! read a pair table and convert it to dot-bracket notation.
 *  @param pair_table the encoded rna-structure.
 */
static std::string
getStructure(short *pair_table);


/**
 * ! returns true if a < b, w.r.t. lexicographical order.
 * @param a the first rna-structure as pairtable.
 * @param b the second rna-structure as pairtable.
 */
static bool
IsSmaller(short *a,
          short *b);


/**
 * ! returns true if a > b, w.r.t. lexicographical order.
 * @param a the first rna-structure as pairtable.
 * @param b the second rna-structure as pairtable.
 */
static bool
IsGreater(short *a,
          short *b);


/**
 * ! returns true if a == b, w.r.t. the dot-bracket structure.
 * @param a the first rna-structure as pairtable.
 * @param b the second rna-structure as pairtable.
 */
static bool
IsEqual(const short *a,
        const short *b);


/**
 * ! check if the sequence is a valid RNA sequence (only a|A|c|C|g|G|u|U is allowed)
 * @param sequence the RNA sequence.
 */
static bool
IsValidSequence(std::string sequence);


/**
 * ! check if the structure is in valid dot-bracket notation (only .|(|) is allowed)
 * @param structure the dot-bracket structure.
 */
static bool
IsValidStructure(std::string structure);


/**
 * ! Full validation of an RNA sequence: non-empty and only a|A|c|C|g|G|u|U.
 * @param sequence the RNA sequence.
 * @param error set to a human-readable reason if the sequence is invalid.
 * @return true if the sequence is valid.
 */
static bool
ValidateSequence(const std::string& sequence,
                 std::string&       error);


/**
 * ! Full validation of a dot-bracket structure against its sequence.
 * Checks the alphabet, that the length equals the sequence length, that the
 * brackets are balanced, that every base pair is canonical (AU, GC, GU) and
 * that every base pair encloses at least minLoopSize unpaired bases.
 * @param structure the dot-bracket structure.
 * @param sequence the RNA sequence the structure belongs to.
 * @param minLoopSize the minimal hairpin loop size of the energy model.
 * @param error set to a human-readable reason if the structure is invalid.
 * @return true if the structure is valid.
 */
static bool
ValidateStructure(const std::string& structure,
                  const std::string& sequence,
                  int                minLoopSize,
                  std::string&       error);
};

#endif /* STRUCTUREUTILS_H_ */
