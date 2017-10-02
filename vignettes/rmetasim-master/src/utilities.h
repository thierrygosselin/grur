/*

$Modified: astrand $

Copyright (C) 1999 Allan E. Strand

This file is part of Metasim

This file contains the interfaces of utility functions for Metasim 

*/

#ifndef UTILITIES_H
#define UTILITIES_H


#include "const.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>

#include <stdlib.h>

/*
extern "C" {
#include <randlib.h>
};
*/


//extern void SetAll(long l1, long l2);

extern int Error(int errnum);



//comparison operator for using strings in STL maps

using namespace std;
struct ltstr
{
  inline bool operator()(const string s1, const string s2) const
  {
    return (s1 < s2);
  }
};


#endif /* UTILITIES */
/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; minor-mode: font-lock ***
;;; End: ***
*/
