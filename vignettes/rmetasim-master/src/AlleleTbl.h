/* This file is part of Metasim
   This file is the decalration of the Allele lookup table
*/

#ifndef ALLELE_LOOKTBL_H
#define ALLELE_LOOKTBL_H

/* includes */
#include <FastAllele.h>
#include <FastSeqAllele.h>

/**
This class implements a landscape-wide collection og sequence allele lookup tables
 */
class AlleleLookTbl: public BaseObj {
protected:
  ///List of pointers to allele tables
  vector <AlleleTbl *> Atbl;
public:
  AlleleLookTbl();
  ~AlleleLookTbl();

  ///Add a new allele tbl to the lookup
  void push_back(AlleleTbl * atp);
  ///clear all tables from the lookup
  void clear();
  ///number of loci
  inline size_t size()
    {
      return Atbl.size();
    }

  ///overload [] so that it returns a reference to the AlleleTbl
  inline AlleleTbl * operator[](int ind)
    {
      return Atbl[ind];
    }
  void DummyFreq(int ps=1000000);
  void ZeroFreq();

  friend ostream &operator<<(ostream & stream, AlleleLookTbl &a);
  friend istream &operator>>(istream & stream, AlleleLookTbl &a);

}; //end AlleleLookTbl

//extern AlleleLookTbl Atbls;

#endif /*ALLELE_LOOKTBL_H*/

/*
;;; Local Variables:        ***
;;; mode: C++               ***
;;; minor-mode:  font-lock  ***
;;; End:                    ***
*/
