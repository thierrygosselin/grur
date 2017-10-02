/**

   Constants for the metasim system

const.h
 */

#ifndef CONST_H
#define CONST_H


/**
   This is a flag for compilation of debugging code inserted in production code.
 */
//#define DEBUG


/**

Constants

*/

/**
Largest random number possible given the stdlib random number generator
*/
#define MAXRANDOM             2147483647.0
#define MAXIDS                1000000000

#define HABDELTA              0.01
#define RANDOMSTATE           100
#define MAXSTAGES             500
#define MAXLOCI               20000
#define MAXALLELES            2000
#define MAXPLOIDY             2
#define MAXSEQLEN             1000
#define GENLEN                MAXLOCI*MAXPLOIDY
#define TOKENLEN              10
#define RESERVE_CLASS_SIZE    6000
#define SORT_I_THRESHOLD      200


#define ONELESS       0.999999
#define ONEMORE       1.000001




/*
  Object identifiers
 */

#define BASEOBJ  1

#define SITEOBJ  100

#define ALLELE 200
#define SEQALLELE 201


#define ALLELETBL     250
#define INFALLELETBL  251
#define STEPALLELETBL 252
#define SEQALLELETBL  253


#endif /*CONST_H*/

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; End: ***
*/

