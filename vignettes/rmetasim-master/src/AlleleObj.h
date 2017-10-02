/*

$Modified: astrand $

Copyright (C) 1999 Allan E. Strand

This file is part of Metasim

This is the declaration of the Allele object types.

*/

#ifndef ALLELEOBJ_H
#define ALLELEOBJ_H

/*includes
*/

#include <metasim.h>
#include <BaseObj.h>
#include <RandLib.h>
#include <iostream>
#include <R.h>

using namespace std;

///Allele Class

/**
   An allele with an integer-based state.  This allele is maintained
   in a lookup table and possesses a frequency, a time of birth a time
   it is lost from the system and a proportion in the system overall
   
 */
class Allele: public BaseObj {
protected: 

  int state;
  int birth;
  int freq;
  double prop;

public:

  Allele (double p=0, int s=0, int b=0, int d=0);
  virtual ~Allele ();

  inline void SetFreq (int cp)
    {
      freq=cp;
    }

  inline int GetFreq ()
    {
      return freq;
    }

  inline void SetProp (double pr)
    {
      prop=pr;
    }

  inline double GetProp ()
    {
      return prop;
    }

  inline void SetState (int st)
    {
      state=st;
    }

  inline int GetState ()
    {
      return state;
    }


  inline void SetBirth (int b)
    {
      birth = b;
    }

  inline int GetBirth ()
    {
      return birth;
    }

    inline void SetDeath (int /*d*/)
    {
      //      death = d;
    }

  inline int GetDeath ()
    {
      return 0;
    }
  
  virtual void WriteState(ostream &stream);
  virtual void Write(ostream & stream);
  virtual void Scan(istream & stream);
  ///compares two alleleles returns true or false
  //  virtual bool Compare(Allele obj);		

  friend ostream & operator<<(ostream & stream, Allele &a);
  friend istream & operator>>(istream & stream, Allele &a);

};


///A molecular sequence

typedef std::vector<char> SequenceType ;



///Sequence Allele Class

/**  SeqAllele implements a class that has a variable length sequence
that represents an alleles sequence

This class uses the vector container from the STL for storing the
sites that comprise an sequence based allele */

/// the sequence of a single allele
class SeqAllele : public Allele { 
  SequenceType dnaseq ;
public:

  SeqAllele (int sl=1);                
  ~SeqAllele () ;

  /// The probability that a particular base will change is 1/(sequence length).
  void mutate ();

  ///returns length of sequence
  size_t SeqLen () ;                       
 

  inline void SetState (SequenceType seq)
    {
      dnaseq=seq;
    }
  inline int GetSeqSize()
  {
    return dnaseq.size();
  }
  inline SequenceType GetState ()
    {
      return dnaseq;
    }

  ///set a site to a state
  ///defaults to first site
 
  void SetSite (char st, int sn = 0);

  /** RandomSeq expects the proportions of the bases in a random
  sequence.  They must add to 1.0 (defaults to 0.25,0.25,0.25,0.25
  It then produces a random sequence with the given base composition
  */

  void RandomSeq (double a=0.25, double c=0.25, double t=0.25, double g=0.25);

  /// returns the site at a particular position
  char GetSite (int sn = 0); 

  ///compares two alleleles returns true or false
  bool Compare(SeqAllele obj);		

  void WriteState(ostream &stream);
  void Write(ostream & stream);
  void Scan(istream & stream);

}; // end seqallele

/// this overloaded "==" operator calls the compare method for SeqAllele
bool operator==(SeqAllele allele1, SeqAllele allele2);  


#endif /*ALLELEOBJ*/

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; minor-mode: font-lock ***
;;; End: ***
*/





