/**

$Modified: astrand $

Copyright (C) 1999 Allan E. Strand

This file is part of Metasim

The implementation of the demo class object

*/


/*includes
*/

#include <Democlass.h>

DemoClass::DemoClass ()
{
  maxind = 0;
  indcnt = 0;
}
DemoClass::~DemoClass ()
{
}

///add an individual to the data structure.  Returns the index to the individual.
///this method also adds the alleles to the allele table for each locus.  This maintains the allele
///frequency tables.
int DemoClass::AddIndividual (PackedIndividual & PkInd)
{
  int i;
  i = -1;
  if (UNUSED.empty())
    {
      I[maxind]=PkInd;
      i=maxind;
      maxind++;
    } //end     "if UNUSED is empty"
  else
    {
      i = UNUSED.back();
      UNUSED.pop_back();
      I[i]=PkInd;
    }
  return i;
}

void DemoClass::ClearClass (int t,AlleleLookTbl &Atbls)
{
  map<int,PackedIndividual,less <int> >::iterator iiter;
  iiter=I.begin();
  while (iiter!=I.end())
    {
      (*iiter).second.Death(t,Atbls);
      iiter++;
    }
  I.clear();
  UNUSED.clear();
  maxind=0;
}
/*
int DemoClass::GetRandomIndex ()
{
  int c=0;
  //  int f=0;
  int indx;
  map<int,PackedIndividual,less <int> >::iterator iiter;

  if (I.size()>0)
    {
      iiter=I.begin();
      indx = RandLibObj.unirange(I.size());
      c=0;
      while ((c!=indx)&&(iiter!=I.end()))
	  {
	    iiter++;
	    c++;
	  }	

      return (*iiter).first;
    }
  else
    {
      return -1;
    }
}  
*/
int DemoClass::GetRandomIndex ()
{
  int f=0;
  int indx;

  while (f==0) // keep trying until an individual is found
    {

      //      cerr <<"maxind "<< maxind << endl;
      //      indx = RandLibObj.unirange(maxind); //fix the bug JDR found (5/31/2014)
      indx = RandLibObj.unirange((maxind-1));

      //      cerr << "RandomIndex "<<indx<<endl;

      if (I.find(indx)!=I.end())
	{
	  f=1;
	}
      else
	{
	  f=0;
	}
    }
  return indx;
}  

void DemoClass::RemoveRandomInd (int t,AlleleLookTbl &Atbls)
{
  while (!RemoveInd(RandLibObj.unirange(maxind),t,Atbls)) // keep trying until an individual is erased
    {
    }
}  

void DemoClass::CompressClass (double frac)
{
  vector <PackedIndividual> tvec;
  size_t i,sz;
  if (I.size()>0&&(I.size()<maxind*frac))
    {
      tvec.reserve(I.size());
      ResetIndividuals();
      i=0;
      for (nextind=I.begin();nextind!=I.end();nextind++)
	{
	  tvec.push_back((*nextind).second);
	}
      I.clear();
      UNUSED.clear();
      maxind=0;
      sz=tvec.size();
      for (i=0;i<sz;i++)
	{
	  I[i]=tvec[i];
	  maxind++;
	}
      ResetIndividuals();
    }
}

double DemoClass::GenLength (int t)
{
  double genoff, totoff;
  //  int indx;


  int lr, no;

  if (I.size()>0)
    {
      genoff=0.0;
      totoff=0.0;
      ResetIndividuals();
      do
	{
	  //	  indx = GetCurrentIndex();
	  lr = GetCurrentLastRep();
	  no =  GetCurrentNumOff();
	  genoff += ((t - lr) * no);
	  totoff += no;
	}
      while (!NextIndividual());
      if (totoff==0) 
	{
	  return 0;
	}
      else
	{
	  return double(genoff/totoff);
	}
    }
else
  {
    return 0.0;
  }
}

ostream &operator<<(ostream & stream, DemoClass & DC)
{
  size_t i;
  PackedIndividual tmpI;
  DC.ResetIndividuals();
  for (i=0;i<DC.I.size();i++)
    {
      stream << (*DC.nextind).second ;
      DC.NextIndividual() ;
    }
  return stream;
}

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; minor-mode: font-lock ***
;;; End: ***
*/
