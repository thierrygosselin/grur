/* This file is part of Metasim
   This file is the decalration of the Allele lookup table
*/

#ifndef SEQALLELE_TBL_H
#define SEQALLELE_TBL_H

/* includes */
#include <FastAllele.h>


/**
This class implements a sequence allele lookup table
 */
class SeqAlleleTbl: public AlleleTbl {
protected:
  std::map<int, SeqAllele, less <int> > A;
  int seqlen;
public:
  SeqAlleleTbl();
  ~SeqAlleleTbl();

  /// returns allele at index i if i is in range.	
  inline SeqAllele getAllele(int i)
  {
    map<int, SeqAllele, less<int> >::iterator tmpiter;
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
	    cerr << "allele index: "<<i<<" not found in SeqAllele.h::getAllele"<<endl;
#endif
	    assert(tmpiter!=A.end());
	    return (*tmpiter).second; //will never reach this statement.  Present to avoid compile complaints
	  }
      }
    else
      {
#ifdef DEBUG
	    cerr << "allele table empty in SeqAllele.h::getAllele"<<endl;
#endif
	    assert(A.size()>0);

	    return (*tmpiter).second; //will never reach this statement.  Present to avoid compile complaints
      }
  }		 

/// returns an allele state pointed to by i
  inline SequenceType getAlleleState(int i)
    {
      SeqAllele ta(1);
      ta=getAllele(i);
      return ta.GetState();
    }

   void getAlleleRef(int i, Allele* ptr)
  {
    *(dynamic_cast<SeqAllele *>(ptr)) = getAllele(i);
  }

  int addAlleleAndIndexRef(Allele* ptr, int i)
  {
    return addAlleleAndIndex(*(dynamic_cast<SeqAllele *>(ptr)),i);
  }

  ///set up the frequency of alleles based upon an assumed popsize and their proportions
  void dummyfreq(int ps=1000000);

  ///Set the freqnuency of all alleles to zero.
  void zerofreq();

  ///Increases the allele counter for the particular index.
  inline void AddAlleleFreq(int i)
    {
      map<int, SeqAllele, less<int> >::iterator tmpiter;
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
  SeqAllele getRandAllele();		 
  /// returns an index to an allele at random
  int getRandAlleleIndex();		 
  ///returns the total number of alleles in this table * their respective freq in the entire pop system
  int AlleleTotalCnt();
  /// adds an allele to table.  returns the index of the added allele (or found existing allele)
  int addAllele(SeqAllele  na, int t=0);	 
  /// adds an allele to table.  also specifies the index to use (used when reading landscape files)
  /// returns the index of the added allele (or found existing allele)
  int addAlleleAndIndex(SeqAllele  na, int ai);	 
  ///takes an allele state and adds to the freq of the appropriate allele
  ///returns the allele index
  int addAlleleState(SequenceType seq, int t=0);	 

  /// returns the number of alleles for this table
  inline int getAlleleNum() 
    {
      return A.size();
    }		 

  inline int getFreq(int i)
    {
      int tmp = A[i].GetFreq();
      return tmp;
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
    seqlen = sl;
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


  friend ostream &operator<<(ostream & stream, SeqAlleleTbl &a);
  friend istream &operator>>(istream & stream, SeqAlleleTbl &a);


}; //end SeqAlleleTbl


#endif /*SEQALLELE_TBL_H*/

/*
;;; Local Variables:        ***
;;; mode: C++               ***
;;; minor-mode:  font-lock  ***
;;; End:                    ***
*/
