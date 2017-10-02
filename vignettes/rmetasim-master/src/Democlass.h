/*

$Modified: astrand $

Copyright (C) 1999 Allan E. Strand

This file is part of Metasim
*/

#ifndef DEMOCLASS_H
#define DEMOCLASS_H

/*includes
*/

#include <metasim.h>
#include <PackedIndividual.h>

using namespace std;


class DemoClass {

  ///the class 
  int cl;
  ///the next new index 
  int maxind;
  ///the structure containing the actual individuals
  map<int, PackedIndividual, less<int> > I ; 

  ///an iterator that points to the next individual to pull from the data structure
  map<int, PackedIndividual, less<int> >::iterator nextind ; 

int indcnt ;

  ///a stack containing indices that have been allocated, but have then been deleted
  vector <int> UNUSED;

public:
  DemoClass ();
  ~DemoClass ();
  
  inline void SetClass(int c)
    {
      cl = c;
    }

  ///add an individual to the data structure.  Returns the index to the individual.   

int AddIndividual (PackedIndividual & ind);

  ///return an individual pointed to by 'ind' If ind not found, still
  ///returns an individual, but a nonsense individual with a class = -1 ;
inline  PackedIndividual GetIndividual(int ind=0)
    {
      PackedIndividual TI;
      map<int,PackedIndividual,less <int> >::iterator iiter;
      
      TI.SetClass(-1);
      
      iiter = I.find(ind);
	if (iiter!=I.end())
	  {
	    return (*iiter).second ;
	  }
	else
	  {
	    return TI;
	  }
    }

int GetRandomIndex();


inline PackedIndividual GetRandomInd()
    {
      PackedIndividual tmpI;
      tmpI.SetClass(-1);
      int indx= GetRandomIndex();
      if (indx>=0)  // fix bug the JDR found (5/31/2014)
	{
	  return GetIndividual(indx);
	}
      else
	{
	return tmpI;
	}
    }  

  /*
inline PackedIndividual GetRandomInd()
  {
      PackedIndividual tmpI;
      map<int,PackedIndividual,less <int> >::iterator biter,fiter;
      map<int,PackedIndividual,less<int> > res;
      insert_iterator < map<int,PackedIndividual,less<int> > >  resi(res) ;
      tmpI.SetClass(-1);
      biter=I.begin();
      fiter=I.end();
      
      __gnu_cxx::random_sample_n(biter,fiter,resi,1);
      return  tmpI;
    }  
  */
inline void ResetIndividuals()
    {
      nextind = I.begin(); 
      indcnt = 0;
    }


  ///repeated execution of this function gets successive individuals
  ///from the data structure.  Returns an individual with class = -1 if
  ///no more persent in the strucutre


inline  PackedIndividual  GetCurrentIndividual()
{
  PackedIndividual TI;
  TI.SetClass(-1);

  if (nextind!=I.end())
    {
      TI = (*nextind).second;
    }
  return TI;
}


inline  int  GetCurrentIndex()
{
  int i = -1;

  if (nextind!=I.end())
    {
      i = (*nextind).first;
    }
  return i;
}


  ///Clears the class 
void ClearClass (int t, AlleleLookTbl &Atbls);


  ///Removes the individual indicated by 'ind' from the structure.  returns 0 if ind was not found, else 1
inline int RemoveInd (int ind, int t, AlleleLookTbl &Atbls)
{
  int er;
  map<int,PackedIndividual,less <int> >::iterator iiter;
  er = 0;
  iiter = I.find(ind);
  if (iiter!=I.end())
    {
      if (ind==maxind)
	{
	  maxind--;
	}
      else
	{
	  UNUSED.push_back(ind);
	}
      (*iiter).second.Death(t,Atbls);
      I.erase(iiter);
      er = 1;
    }
  return er;
}


inline int NextIndividual()
{
  
  if (nextind!=I.end())
    {
      nextind++;
      indcnt++;
      if (nextind!=I.end())
	{
	  return 0;
	}
      else 
	{
	  return 1;
	}
    }
  else
    {
#ifdef DEBUG
      cerr << "trying advance the nextind iterator even though it is already at the end of the list Democlass.cc:NextIndividual()";
#endif
      return 1;
    }
}

inline int RemoveCurrentInd(int t, AlleleLookTbl &Atbls)
{
  map<int,PackedIndividual,less <int> >::iterator tmpiter;

  if (nextind!=I.end())
    {
      tmpiter = nextind;
      (*tmpiter).second.Death(t,Atbls);

      //      cerr << "In Democlass.cc::RemoveCurrentInd()  deleting:  "<<(*tmpiter).second<<" indcnt = " << indcnt <<endl;

      UNUSED.push_back((*tmpiter).first); //push index onto list of availible indices
      if (I.size()<1)
	{
#ifdef DEBUG
	  cerr << "trying to erase an individual when none exist" << endl;
#endif
	  assert(I.size()>0);
	}
      else
	{
	  I.erase(tmpiter);
	}
      return 1;
    }
  return 0;
}

inline size_t size()
{
  return I.size();
}

void RemoveRandomInd (int t, AlleleLookTbl &Atbls);

inline void ChangeCurrentInd(int t)
{
  (*nextind).second.Change(t);
}
///
inline void ChangeInd(int ind, int t)
{
  if (I.find(ind)!=I.end())
    {
      I[ind].Change(t);
    }
  else
    {
#ifdef DEBUG
      cerr  <<"couldnt find individual to change "<<endl;
#endif
      assert (I.find(ind)!=I.end());
    }
}

inline int GetChangedCurrentInd()
{
  return (*nextind).second.GetChanged();
}

inline void SetCurrentLastRep(int lr)
{
  if (nextind!=I.end())
    {
      (*nextind).second.SetLastRep(lr);
    }
  else
    {
#ifdef DEBUG
      cerr << "past end of individual list"<<endl;
#endif
      assert(nextind!=I.end());
    }
}  
///Return the last time click the current individual reproduced
inline int GetCurrentLastRep()
{
  if (nextind!=I.end())
    {
      return (*nextind).second.GetLastRep();
    }
  else
    {
      return -1;
    }
}
  
inline void SetCurrentNumOff(int no)
{
  if (nextind!=I.end())
    {
      (*nextind).second.SetNumOff(no);
    }
  else
    {
#ifdef DEBUG
      cerr << "past end of individual list"<<endl;
#endif
      assert(nextind!=I.end());
    }
}  
///Return the number of offspring last reproductive event
inline int GetCurrentNumOff()
{
  if (nextind!=I.end())
    {
      return (*nextind).second.GetNumOff();
    }
  else
    {
      return -1;
    }
}  
  

void CompressClass(double frac=0.5);

double GenLength(int t);


friend ostream &operator<<(ostream & stream, DemoClass & DC);

}; // end DemoClass

#endif /*DemoClass */

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; minor-mode: font-lock ***
;;; End: ***
*/
