/*
 * StructureUtils.cpp
 *
 *  Created on: 05.08.2014
 *      Author: Gregor Entzian
 */

#include <cstring>
#include <cctype>
#include <sstream>
#include <vector>
#include "StructureUtils.h"
#include <cstdlib>
#include <cmath>
#include <climits>
extern "C" {
#include <ViennaRNA/eval.h>
#include <ViennaRNA/utils/structures.h>
}

StructureUtils::StructureUtils()
{
}


StructureUtils::~StructureUtils()
{
}


std::string
StructureUtils::getStructure(short *pair_table)
{
  if (pair_table == NULL) {
    std::cout << "pair_table is empty!" << std::endl;
    return "";
  }

  std::string s(pair_table[0], ' ');
  for (int i = 1; i <= pair_table[0]; i++) {
    if (pair_table[i] == 0)
      s[i - 1] = '.';
    else if (pair_table[i] < i)
      s[i - 1] = ')';
    else
      s[i - 1] = '(';
  }
  return s;
}


bool
StructureUtils::IsSmaller(short *a,
                          short *b)
{
  //rule: "." < "(" < ")".
  if (a != NULL && b != NULL) {
    short encodedCharacterInA;
    short encodedCharacterInB;
    //ascii-lexicographic sort ( < ) < . (like in rna-subopt)
    for (int i = 1; i <= a[0]; i++) {
      if (a[i] != b[i]) {
        if (a[i] == 0)
          encodedCharacterInA = 2;  //="."
        else if (a[i] < i)
          encodedCharacterInA = 1;  //=")"
        else
          encodedCharacterInA = 0;  //="("

        if (b[i] == 0)
          encodedCharacterInB = 2;  //="."
        else if (b[i] < i)
          encodedCharacterInB = 1;  //=")"
        else
          encodedCharacterInB = 0;  //="("

        if (encodedCharacterInA != encodedCharacterInB)
          return encodedCharacterInA < encodedCharacterInB;
      }
    }
  }

  return false;
}


bool
StructureUtils::IsGreater(short *a,
                          short *b)
{
  if (a != NULL && b != NULL) {
    //encoding: "."=0, ")"=2, "(" = 1.
    for (int i = 1; i <= a[0]; i++) {
      if (a[i] != b[i]) {
        short encodedCharacterInA;
        short encodedCharacterInB;
        if (a[i] == 0)
          encodedCharacterInA = 2;  //="."
        else if (a[i] < i)
          encodedCharacterInA = 1;  //=")"
        else
          encodedCharacterInA = 0;  //="("

        if (b[i] == 0)
          encodedCharacterInB = 2;  //="."
        else if (b[i] < i)
          encodedCharacterInB = 1;  //=")"
        else
          encodedCharacterInB = 0;  //="("

        if (encodedCharacterInA != encodedCharacterInB)
          return encodedCharacterInA > encodedCharacterInB;
      }
    }
  }

  return false;
}


bool
StructureUtils::IsEqual(const short *a,
                        const short *b)
{
  if (a != NULL && b != NULL)
    return memcmp(a, b, (a[0] + 1) * sizeof(short)) == 0;

  return false;
}


bool
StructureUtils::IsValidSequence(std::string sequence)
{
  for (size_t i = 0; i < sequence.size(); i++) {
    char letter = sequence.at(i);
    if (!((letter == 'a') | (letter == 'A') | (letter == 'g')
          | (letter == 'G') | (letter == 'c') | (letter == 'C')
          | (letter == 'u') | (letter == 'U')))
      return false;
  }
  return true;
}


bool
StructureUtils::IsValidStructure(std::string structure)
{
  for (size_t i = 0; i < structure.size(); i++) {
    char letter = structure.at(i);
    if (!((letter == '.') | (letter == '(') | (letter == ')')))
      return false;
  }
  return true;
}


bool
StructureUtils::ValidateSequence(const std::string& sequence,
                                 std::string&       error)
{
  if (sequence.empty()) {
    error = "the sequence is empty";
    return false;
  }

  for (size_t i = 0; i < sequence.size(); i++) {
    char c = std::toupper(static_cast<unsigned char>(sequence[i]));
    if (c != 'A' && c != 'C' && c != 'G' && c != 'U') {
      std::ostringstream msg;
      msg << "invalid character '" << sequence[i] << "' at position " << (i + 1)
          << " (only A, C, G, U are allowed";
      if (c == 'T')
        msg << "; replace T by U";

      msg << ")";
      error = msg.str();
      return false;
    }
  }
  return true;
}


static bool
isCanonicalPair(char a,
                char b)
{
  a = std::toupper(static_cast<unsigned char>(a));
  b = std::toupper(static_cast<unsigned char>(b));
  return (a == 'A' && b == 'U') || (a == 'U' && b == 'A')
         || (a == 'G' && b == 'C') || (a == 'C' && b == 'G')
         || (a == 'G' && b == 'U') || (a == 'U' && b == 'G');
}


bool
StructureUtils::ValidateStructure(const std::string& structure,
                                  const std::string& sequence,
                                  int                minLoopSize,
                                  std::string&       error)
{
  std::ostringstream msg;

  if (structure.size() != sequence.size()) {
    msg << "structure length " << structure.size()
        << " does not match sequence length " << sequence.size();
    error = msg.str();
    return false;
  }

  std::vector<size_t> open;  // 0-based positions of unmatched '('
  for (size_t j = 0; j < structure.size(); j++) {
    char c = structure[j];
    if (c == '.') {
      continue;
    } else if (c == '(') {
      open.push_back(j);
    } else if (c == ')') {
      if (open.empty()) {
        msg << "unbalanced brackets: ')' at position " << (j + 1)
            << " has no matching '('";
        error = msg.str();
        return false;
      }

      size_t i = open.back();
      open.pop_back();
      if (!isCanonicalPair(sequence[i], sequence[j])) {
        msg << "non-canonical base pair (" << (i + 1) << "," << (j + 1) << ") "
            << sequence[i] << "-" << sequence[j]
            << " (only AU, GC and GU pairs are allowed)";
        error = msg.str();
        return false;
      }

      if (static_cast<int>(j - i - 1) < minLoopSize) {
        msg << "base pair (" << (i + 1) << "," << (j + 1) << ") encloses only "
            << (j - i - 1) << " unpaired bases (minimum hairpin size is "
            << minLoopSize << ")";
        error = msg.str();
        return false;
      }
    } else {
      msg << "invalid character '" << c << "' at position " << (j + 1)
          << " (only '.', '(' and ')' are allowed)";
      error = msg.str();
      return false;
    }
  }

  if (!open.empty()) {
    msg << "unbalanced brackets: '(' at position " << (open.back() + 1)
        << " has no matching ')'";
    error = msg.str();
    return false;
  }

  return true;
}


bool
StructureUtils::HasLonelyPair(const std::string& structure,
                              size_t&            i,
                              size_t&            j)
{
  short *pt = vrna_ptable(structure.c_str());

  if (pt == NULL)
    return false;

  bool  found = false;
  for (int p = 1; p <= pt[0] && !found; p++) {
    int q = pt[p];
    if (q <= p)
      continue;  // unpaired or closing bracket

    bool outer_stack  = (p > 1 && q < pt[0] && pt[p - 1] == q + 1);
    bool inner_stack  = (pt[p + 1] == q - 1);
    if (!outer_stack && !inner_stack) {
      i     = p;
      j     = q;
      found = true;
    }
  }
  free(pt);
  return found;
}


int
StructureUtils::EvalMove(vrna_fold_compound_t *vc,
                         const vrna_move_t    *m,
                         const short          *pt)
{
  vrna_move_t first = *m;

  first.next = NULL;

  if (m->next == NULL || m->next->pos_5 == 0)
    return vrna_eval_move_shift_pt(vc, &first, (short *)pt);

  // evaluate the parts one after another on a copy that gets each part applied
  short *tmp    = vrna_ptable_copy(pt);
  int   energy  = vrna_eval_move_shift_pt(vc, &first, tmp);
  vrna_move_apply(tmp, &first);
  for (const vrna_move_t *n = m->next; n->pos_5 != 0; n++) {
    energy += StructureUtils::EvalMove(vc, n, tmp);
    vrna_move_apply(tmp, n);
  }
  free(tmp);
  return energy;
}


int
StructureUtils::KcalToDcal(double kcal)
{
  double dcal = std::round(kcal * 100.0);

  if (dcal >= (double)INT_MAX)
    return INT_MAX;

  if (dcal <= (double)INT_MIN)
    return INT_MIN;

  return (int)dcal;
}
