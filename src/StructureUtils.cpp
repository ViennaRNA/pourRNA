/*
 * StructureUtils.cpp
 *
 *  Created on: 05.08.2014
 *      Author: Gregor Entzian
 */

#include "StructureUtils.h"

StructureUtils::StructureUtils ()
{
}

StructureUtils::~StructureUtils ()
{
}

std::string
StructureUtils::getStructure (short * pair_table)
{
  if (pair_table == NULL)
    {
      std::cout << "pair_table is empty!" << std::endl;
      return "";
    }
  std::string s (pair_table[0], ' ');
  for (int i = 1; i <= pair_table[0]; i++)
    {
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
StructureUtils::IsSmaller (short * a, short * b)
{
  //rule: "." > ")" > "(".
  if (a != NULL && b != NULL)
    {
      //encoding: "."=2, ")"=1, "(" = 0.
      for (int i = 1; i <= a[0]; i++)
	{
	  if (a[i] != b[i])
	    {
	      short encodedCharacterInA;
	      short encodedCharacterInB;
	      if (a[i] == 0)
		{
		  encodedCharacterInA = 2; //="."
		}
	      else if (a[i] < i)
		{
		  encodedCharacterInA = 1; //=")"
		}
	      else
		{
		  encodedCharacterInA = 0; //="("
		}

	      if (b[i] == 0)
		{
		  encodedCharacterInB = 2; //="."
		}
	      else if (b[i] < i)
		{
		  encodedCharacterInB = 1; //=")"
		}
	      else
		{
		  encodedCharacterInB = 0; //="("
		}

	      if (encodedCharacterInA != encodedCharacterInB)
		{
		  return encodedCharacterInA < encodedCharacterInB;
		}
	    }
	}
    }

  return 0;

}

bool
StructureUtils::IsGreater (short * a, short * b)
{
  if (a != NULL && b != NULL)
    {
      //encoding: "."=2, ")"=1, "(" = 0.
      for (int i = 1; i <= a[0]; i++)
	{
	  if (a[i] != b[i])
	    {
	      short encodedCharacterInA;
	      short encodedCharacterInB;
	      if (a[i] == 0)
		{
		  encodedCharacterInA = 2; //="."
		}
	      else if (a[i] < i)
		{
		  encodedCharacterInA = 1; //=")"
		}
	      else
		{
		  encodedCharacterInA = 0; //="("
		}

	      if (b[i] == 0)
		{
		  encodedCharacterInB = 2; //="."
		}
	      else if (b[i] < i)
		{
		  encodedCharacterInB = 1; //=")"
		}
	      else
		{
		  encodedCharacterInB = 0; //="("
		}
	      if (encodedCharacterInA != encodedCharacterInB)
		{
		  return encodedCharacterInA > encodedCharacterInB;
		}
	    }
	}
    }
  return 0;
}

bool
StructureUtils::IsEqual (const short * a, const short * b)
{
  if (a != NULL && b != NULL)
    {
      bool bIsEqual = true;
      for (int i = 1; i <= a[0]; i++)
	{
	  if (a[i] != b[i])
	    {
	      bIsEqual = false;
	      break;
	    }
	}
      return bIsEqual;
    }
  else
    {
      return 0;
    }
}

bool
StructureUtils::IsValidSequence (std::string sequence)
{
  for (size_t i = 0; i < sequence.size (); i++)
    {
      char letter = sequence.at (i);
      if (!((letter == 'a') | (letter == 'A') | (letter == 'g') | (letter == 'G') | (letter == 'c')
	  | (letter == 'C') | (letter == 'u') | (letter == 'U')))
	{
	  return false;
	}
    }
  return true;
}

bool
StructureUtils::IsValidStructure (std::string structure)
{
  for (size_t i = 0; i < structure.size (); i++)
    {
      char letter = structure.at (i);
      if (!((letter == '.') | (letter == '(') | (letter == ')')))
	{
	  return false;
	}
    }
  return true;
}
