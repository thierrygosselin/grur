/*

$Modified: astrand $

Copyright (C) 1999 Allan E. Strand

This file is part of Metasim

This is the declaration of the main header file that includes other
files.

*/

#ifndef METASIM_H
#define METASIM_H

/*
  #define __USE_MALLOC
*/

#include "const.h"
#include "utilities.h"

#include <assert.h>



//#include <math.h>
#include <cmath>
#include <time.h>

#include <vector>
#include <map>
#include <algorithm>
/**
 #include <ext/algorithm>
**/
#include <string>

/**
UTILITY MACROS
 */

///#define RDEBUG
///#define DEBUG

//// #define ROUND_2_INT(f) ((int)(f >= 0.0 ? (f + 0.5) : (f - 0.5))) 
#define HAVE_LAPACK

#endif /*METASIM_H*/
/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; minor-mode: font-lock ***
;;; End: ***
*/

























