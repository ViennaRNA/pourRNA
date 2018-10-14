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
  getStructure (short* pair_table);
  /**
   * ! returns true if a < b, w.r.t. lexicographical order.
   * @param a the first rna-structure as pairtable.
   * @param b the second rna-structure as pairtable.
   */
  static bool
  IsSmaller (short * a, short * b);
  /**
   * ! returns true if a > b, w.r.t. lexicographical order.
   * @param a the first rna-structure as pairtable.
   * @param b the second rna-structure as pairtable.
   */
  static bool
  IsGreater (short * a, short * b);
  /**
   * ! returns true if a == b, w.r.t. the dot-bracket structure.
   * @param a the first rna-structure as pairtable.
   * @param b the second rna-structure as pairtable.
   */
  static bool
  IsEqual (const short * a, const short * b);
  /**
   * ! check if the sequence is a valid RNA sequence (only a|A|c|C|g|G|u|U is allowed)
   * @param sequence the RNA sequence.
   */
  static bool IsValidSequence(std::string sequence);
  /**
   * ! check if the structure is in valid dot-bracket notation (only .|(|) is allowed)
   * @param structure the dot-bracket structure.
   */
  static bool IsValidStructure(std::string structure);

};

#endif /* STRUCTUREUTILS_H_ */
