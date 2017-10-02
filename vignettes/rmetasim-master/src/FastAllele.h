/* This file is part of Metasim
   This file is the decalration of the Allele lookup table
*/

#ifndef ALLELE_TBL_H
#define ALLELE_TBL_H

/* includes */
#include <metasim.h>
#include <BaseObj.h>
#include <AlleleObj.h>
#include <RandLib.h>
/**
This class implements a general allele lookup table
 */

class AlleleTbl: public BaseObj {
protected:
  double rate;
  int ploidy;
  ///transmission of alleles:  0=biparental, 1=maternal, 2=paternal
  int trans;
  std::vector<int> UNUSED;


public:

  AlleleTbl();
  virtual ~AlleleTbl();

  virtual void WriteAlleleState(int a, ostream &stream)=0;
  virtual void Write(ostream &stream)=0;
  virtual vector<int> getAindices()=0;
  virtual void Scan(istream &stream)=0;
  virtual void dummyfreq(int ps=1000000)=0;
  virtual void zerofreq()=0;
  /// sets mutation rate for this set of alleles
  inline void setMutationRate(double r=0)
    {
      assert(r>=0 && r<=1);
      rate = r;
    }
  /// returns mutation rate for locus
  inline double getMutationRate() 
    {
      return rate;
    }		 

  /// mutates allele based upon the value of rate (calls the mutation
  ///function of theallele) returns the index of the new allele.  This
  ///function should be overridden by descendants

  virtual int mutator (int anum, int t);

  ///Set the ploidy or number of alleles at this locus
  inline void setPloidy(int pl)
    {
      ploidy=pl;
    }
  /// returns the ploidy or number of alleles in this locus
  inline int getPloidy()
    {
      return ploidy;
    }

  /// Sets the transmission mode of the locus
  inline void setTrans(int tr)
    {
      trans=tr;
    }
  /// Gets the transmission mode of the locus
  int getTrans()
    {
      return trans;
    }

  Allele getAllele(int i);

  //virtual function stubs
  virtual int getRandAlleleIndex()=0;		 
  virtual void getAlleleRef(int i, Allele* ptr)=0;
  virtual int addAlleleAndIndexRef(Allele* ptr, int i)=0;
  virtual void AddAlleleFreq (int i)=0;
  virtual void GCAlleles()=0;
  virtual void KillAlleleCopy(int i, int t)=0;
  virtual void CalcProps()=0;
  virtual void clear()=0;
  virtual void setSeqLen(size_t sl)=0;

};  // end AlleleTbl


/**

   This class implements the Infinite allele model.  Each mutation
   produces an allele that is new to the system

 */
class InfAlleleTbl: public AlleleTbl {
protected:
  std::map<int, Allele, less <int> > A;
  //  std::vector<Allele> DeadAlleles;
  int maxstate;

public:
  InfAlleleTbl();
  ~InfAlleleTbl();
  /// returns allele at index i if i is in range.	
  inline Allele getAllele_InfAlleleTbl(int i)
  {
    map<int, Allele, less<int> >::iterator tmpiter;
    if (A.size()>0)
      {
	tmpiter = A.find(i);
	if (tmpiter!=A.end())
	  {
	    return (*tmpiter).second;
	  }
	else
	  {
#ifdef DEBUG
	    cerr << "allele index: "<<i<<" not found in Allele.h::getAllele"<<endl;
#endif
	    assert(tmpiter!=A.end());
	    return (*tmpiter).second; //will never reach this statement.  Present to avoid compile complaints
	  }
      }
    else
      {
#ifdef DEBUG
	    cerr << "allele table empty in Allele.h::getAllele"<<endl;
#endif
	    assert(A.size()>0);

	    return (*tmpiter).second; //will never reach this statement.  Present to avoid compile complaints
      }
  }		 

  inline Allele getAllele(int i)
  {
    return getAllele_InfAlleleTbl(i);
  }
/// returns an allele state pointed to by i
  inline int getAlleleState(int i)
    {
      Allele ta;
      ta=getAllele(i);
      return ta.GetState();
    }
  

  void getAlleleRef(int i, Allele* ptr)
  {
    *ptr = getAllele(i);
  }

  int addAlleleAndIndexRef(Allele* ptr, int i)
  {
    return addAlleleAndIndex(*ptr,i);
  }

  void dummyfreq(int ps=1000000);
  void zerofreq();

  ///Increases the allele counter for the particular index.
  inline void AddAlleleFreq(int i)
    {
      map<int, Allele, less<int> >::iterator tmpiter;
      if (A.size()>0)
	{
	  tmpiter = A.find(i);
	  if (tmpiter!=A.end())
	    {
	      (*tmpiter).second.SetFreq((*tmpiter).second.GetFreq()+1);
//	      cerr << "incrementing the allele pointed by indx: " <<i<<" allele: "<<(*tmpiter).second;
	    }
	  else
	    {
#ifdef DEBUG
	      cerr << "allele index: "<<i<<" not found in Allele.h::getAllele"<<endl;
#endif
	      assert(tmpiter!=A.end());
	    }
	}
      else
	{
#ifdef DEBUG
	  cerr << "allele table empty in Allele.h::getAllele"<<endl;
#endif
	  assert(A.size()>0);
	}
    }		 



  /// returns an allele at random
  Allele getRandAllele();		 
  /// returns an index to an allele at random
  int getRandAlleleIndex();		 
  ///returns the total number of alleles in this table * their respective freq in the entire pop system
  int AlleleTotalCnt();
  /// adds an allele to table.  returns the index of the added allele (or found existing allele)
  int addAllele(Allele  na, int t=0);	 
  /// adds an allele to table.  also specifies the index to use (used when reading landscape files)
  /// returns the index of the added allele (or found existing allele)
  int addAlleleAndIndex(Allele  na, int ai);	 
  ///takes an allele state and adds to the freq of the appropriate allele
  ///returns the allele index
  int addAlleleState(int is, int t=0);	 

  /// returns the number of alleles for this table
  inline int getAlleleNum() 
    {
      return A.size();
    }		 

  inline int getFreq(int i)
    {
      return A[i].GetFreq();
    }
  inline double getProp(int i)

    {
      return A[i].GetProp();
    }

  inline void setFreq(int i,int fr)
    {
      A[i].SetFreq(fr);
    }

  void clear();

  int mutator(int anum, int t);  

  inline void setSeqLen(size_t sl)
  {
    sl++;  ///doesn't do anything for this class, defined as virtual in parent.
  }


  void SetMaxState();
  void GCAlleles();
  void KillAlleleCopy(int i, int t);

  ///take the frequencies stored in the allele table and calculate the
  ///proportions of alleles
  void CalcProps();

  void WriteAlleleState(int a, ostream &stream);
  void Write(ostream & stream);
  vector<int>  getAindices();
  void Scan(istream & stream);


  friend ostream &operator<<(ostream & stream, InfAlleleTbl &a);
  friend istream &operator>>(istream & stream, InfAlleleTbl &a);


}; //end InfAlleleTbl


/**

   This class implements the stepwise mutation model. The state of the allele determines the state of any mutant.  This is a strict model so that an allele can only mutate to it's neighbors.  If the state is zero, then the allele cannot get any smaller.


 */
class StepAlleleTbl: public InfAlleleTbl {

public:
  StepAlleleTbl();
  ~StepAlleleTbl();

  inline Allele getAllele(int i)
  {
    return getAllele_InfAlleleTbl(i);
  }


  int mutator(int anum, int t);

}; //end StepAlleleTbl



/**

This is a GLOBAL variable.  It stores a series of lookup tables that
specify characteristics of alleles at each locus

 */

#endif /*ALLELE_TBL_H*/

/*
;;; Local Variables:        ***
;;; mode: C++               ***
;;; minor-mode:  font-lock  ***
;;; End:                    ***
*/
