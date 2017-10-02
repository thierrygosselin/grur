/* This file is part of Metasim
   This file is the implementation of the Allele lookup table
*/

/* includes */
#include <FastAllele.h>

/**
This class implements a general allele lookup table
 */



AlleleTbl::AlleleTbl()
{
  setMutationRate(0);
  UNUSED.reserve(500);
  setClassType(ALLELETBL);
}

AlleleTbl::~AlleleTbl()
{
#ifdef DEBUG
  cerr << "executing AlleleTbl destructor"<<endl;
#endif
}


///t is the current clock-tick
int AlleleTbl::mutator(int /*anum*/, int /*t*/)
{
#ifdef DEBUG
  cerr << "AlleleTbl mutator called.  This function should be overridden"<<endl;
  assert(1==0);
#endif
  return -1;
}



/**

   Begin implementation of Infinite allele table

 */

///t is the current clock-tick
InfAlleleTbl::InfAlleleTbl()
{
  clear();
  setMutationRate(0);
  UNUSED.reserve(500);
  setClassType(INFALLELETBL);
  maxstate=0;
}

InfAlleleTbl::~InfAlleleTbl()
{

  clear();
  setMutationRate(0);
  UNUSED.resize(0);

#ifdef DEBUG
  cerr << "executing InfAlleleTbl destructor"<<endl;
#endif

  maxstate=0;
}

Allele InfAlleleTbl::getRandAllele()
{
  return getAllele(getRandAlleleIndex());
}

/// returns an allele index at random
int InfAlleleTbl::getRandAlleleIndex()
{
  int sz,i,tofind;
  sz = A.size();
  double *p = new double[sz];
  int *lookup = new int[sz];

  map<int, Allele, less<int> >::iterator tmpiter;

  assert (sz>0);

  i=0;
  for (tmpiter=A.begin();tmpiter!=A.end();tmpiter++)
    {
      p[i]= (*tmpiter).second.GetProp();
      lookup[i]=(*tmpiter).first;
      i++;
    }

//  tmpiter = A.end();

  do 
    { 
      tofind = lookup[RandLibObj.multinomial(p,sz)];
      tmpiter = A.find(tofind);
    }
  while(tmpiter==A.end());

  delete[] lookup;
  delete[] p;

  return (*tmpiter).first;  
}		 

int InfAlleleTbl::AlleleTotalCnt()
{
  int tot=0;
  map<int, Allele, less<int> >::iterator tmpiter;
  for (tmpiter=A.begin();tmpiter!=A.end();tmpiter++)
    {
      tot= tot + (*tmpiter).second.GetFreq();
    }
  return tot;
}

void InfAlleleTbl::dummyfreq(int ps)
{
  map<int, Allele, less<int> >::iterator tmpiter;
  for (tmpiter=A.begin();tmpiter!=A.end();tmpiter++)
    {
      (*tmpiter).second.SetFreq(int(ceil(ps * (*tmpiter).second.GetProp())));
    }
}
void InfAlleleTbl::zerofreq()
{
  map<int, Allele, less<int> >::iterator tmpiter;
  for (tmpiter=A.begin();tmpiter!=A.end();tmpiter++)
    {
      (*tmpiter).second.SetFreq(0);
    }
}

void InfAlleleTbl::CalcProps()
{
  int i,sz;
  double np,tot;
  map<int, Allele, less<int> >::iterator tmpiter;
  sz=A.size();
  tot=0;
  for (tmpiter=A.begin();tmpiter!=A.end();tmpiter++)
    {
      if ((*tmpiter).second.GetFreq()>0)
	{
	  np = double((*tmpiter).second.GetFreq())/double(AlleleTotalCnt());
	}
      else
	{
	  np = 0.0 ;
	}

      (*tmpiter).second.SetProp(np);     //update proportions of alleles
      tot+=np;
    }
  if (tot>1) //proportions are total to more than one.  Can happen from rounding error.  If small, worth scaling
             // entire vector
    {

      assert((1 - (1/tot))< 0.05);
      for (i=0;i<sz;i++)
	{
	  (*tmpiter).second.SetProp((*tmpiter).second.GetProp()/tot) ;
	}      
    }
}

/// adds an allele to table
int InfAlleleTbl::addAllele(Allele  na, int t)
{
  int f;
  int anum;
  Allele a;
  map<int, Allele, less<int> >::iterator tmpiter;

  anum=-1;
  f=0;
  tmpiter=A.begin();
  while (((f==0)&&(tmpiter!=A.end()))&&(A.size()>0))
    {
      if ((*tmpiter).second.GetState()==na.GetState())
	{
	  f=1; //found = true
	  (*tmpiter).second.SetFreq((*tmpiter).second.GetFreq()+1); // always keep count of alleles
	  anum=(*tmpiter).first;
	}
      tmpiter++;
    }
  if (f==0)
    {
      if (na.GetState()>maxstate)
	{
	  maxstate = na.GetState();
	}
      a.SetState(na.GetState());
      a.SetBirth(t);
      if (na.GetFreq()>0)
	{
	  a.SetFreq(na.GetFreq());
	}
      else
	{
	  a.SetFreq(1);
	}
      a.SetProp(double(a.GetFreq())/double(AlleleTotalCnt()));

      //find an unused allele index
      //first check and see if there are allele indices that are unused.  If so, use one and remove it
      //from the unused pile

      if (UNUSED.size()>0)
	{
	  anum = UNUSED.back();
	  UNUSED.pop_back();
	}
      else //if UNUSED is empty
	{
	  anum=0;
	  while (A.find(anum)!=A.end())
	    {
	      anum++;
	    }
	}
      A[anum]=a;
    }
  assert(anum>-1);

  return anum;
 }		 


int InfAlleleTbl::addAlleleAndIndex(Allele  na, int ai)
{
  int f;
  int anum;
  Allele a;
  map<int, Allele, less<int> >::iterator tmpiter;

  anum=-1;
  f=0;
  tmpiter=A.begin();
  while (((f==0)&&(tmpiter!=A.end()))&&(A.size()>0))
    {
      if ((*tmpiter).second.GetState()==na.GetState())
	{
	  if (ai==(*tmpiter).first)
	    {
#ifdef DEBUG
	      cerr << "Allele index: "<<ai<< " already present in table " <<endl;
	      assert(0==1);
#endif
	    }
	}
      tmpiter++;
    }
  if (f==0)
    {
      if (na.GetState()>maxstate)
	{
	  maxstate = na.GetState();
	}
      a.SetState(na.GetState());
      a.SetBirth(na.GetBirth());
      a.SetFreq(na.GetFreq());
      a.SetProp(na.GetProp());

      A[ai]=a;
      anum = ai;
    }
  assert(anum>-1);

  return anum;
 }		 




///takes a state and adds an allele
int InfAlleleTbl::addAlleleState(int is, int t)
{
  Allele a(is,0,0);
  return addAllele(a,t);
}

///Go through the allele table and remove those alleles with freq=0
void InfAlleleTbl::GCAlleles()
{
  map<int, Allele, less<int> >::iterator tmpiter;
  for (tmpiter = A.begin();tmpiter!=A.end();tmpiter++)
    {
      if ((*tmpiter).second.GetFreq()<=0)
	{
	  UNUSED.push_back((*tmpiter).first); //keep the index to resuse later
	  A.erase(tmpiter); // erase the allele
	}
    }
}


void InfAlleleTbl::KillAlleleCopy(int i, int /*t*/)
{
  map<int, Allele, less<int> >::iterator tmpiter;
  tmpiter = A.find(i);
  if (tmpiter!=A.end())
    {
      if ((*tmpiter).second.GetFreq()>0)
	{
	  (*tmpiter).second.SetFreq((*tmpiter).second.GetFreq() - 1);
	}
      if ((*tmpiter).second.GetFreq()<=0)
	{
	  UNUSED.push_back(i); //keep the index to resuse later
	  A.erase(tmpiter); // erase the allele
	}
    }
  else
    {
#ifdef DEBUG
      cerr << "allele index : "<<i<<" not found in Allele.h::KillAlleleCopy"<<endl;
      cerr << "This is the allele table: " <<endl;
      Write(cerr);
#endif
      assert(tmpiter!=A.end());
    }

}


void InfAlleleTbl::SetMaxState()
{
  map<int, Allele, less<int> >::iterator tmpiter;
  int maxstate=0;

  //chug through all of the alleles and look for the max state
  for (tmpiter=A.begin();tmpiter!=A.end();tmpiter++)
    {
      if (maxstate<(*tmpiter).second.GetState()) 
	{ 
	  maxstate = (*tmpiter).second.GetState() ;
	}
    }
}

void InfAlleleTbl::clear()
{
  ploidy=0;
  trans=0;
  rate=0.0;
  A.clear();
  UNUSED.resize(0);
}

int InfAlleleTbl::mutator(int anum, int t)
{
  map<int, Allele, less<int> >::iterator tmpiter;

  assert(anum>=0);
  if (RandLibObj.uniform()<rate)    //a mutation has occurred
    {
  
      Allele na;
      int newanum;


      SetMaxState();
      na.SetState(maxstate+1);  //sets state to a previously unused value.
      na.SetBirth(t);
      na.SetFreq(1);

      newanum = addAllele(na,t);
      return newanum;
    }
  else
    {
  
      tmpiter=A.find(anum);
      if (tmpiter!=A.end())
	{
	  ///update the allele freq 
	  ///without having to lookup allele by value
	  (*tmpiter).second.SetFreq((*tmpiter).second.GetFreq()+1);  
	  return anum;
	}
      else
	{
#ifdef DEBUG
	  cerr <<"Allele number "<<anum<<" not found in allele table: " <<endl<<*this<<endl;
#endif
	  assert(tmpiter!=A.end());
	  return -1;
	}
    }
}

void InfAlleleTbl::WriteAlleleState(int a, ostream &stream)
{
  Allele al;
  al = getAllele(a);
  al.WriteState(stream);
}

void InfAlleleTbl::Write(ostream &stream)
{
  int sz;
  int ai;
  map<int, Allele, less<int> >::iterator tmpiter;
  Allele al;

  //don't clean table while debugging
#ifndef DEBUG
  CalcProps();
  GCAlleles();
#endif

  sz=A.size();

  stream << sz << endl;
  stream << rate << endl;
  stream << ploidy << endl;
  stream << trans << endl;

  for (tmpiter=A.begin();tmpiter!=A.end();tmpiter++)
    {
	  ai = (*tmpiter).first;
	  al = (*tmpiter).second;
	  stream << ai << "  " << al;
    }
  stream << endl;
}		 

vector<int>  InfAlleleTbl::getAindices()
{
  //  int sz;

  vector<int> aindices;
  map<int, Allele, less<int> >::iterator tmpiter;

  //sz=A.size();

  for (tmpiter=A.begin();tmpiter!=A.end();tmpiter++)
    {
	  aindices.push_back((*tmpiter).first);
    }
  return aindices;
}		 

void InfAlleleTbl::Scan(istream &stream)
{
  Allele newa;
  int i;
  int ai;
  //  int tmp;
  int numa;

  double tprop;
  tprop = 0;

  clear();

  stream >> numa;
  stream >> rate;
  stream >> ploidy;
  stream >> trans;
  for (i=0;i<numa;i++)
    {
      stream >> ai >> newa;
      addAlleleAndIndex(newa,ai);
      
      tprop = newa.GetProp() + tprop;
      if (newa.GetState()>maxstate)
	{
	  maxstate=newa.GetState();
	}
    }
  if (tprop != 1.0) //primitive error checking
    {
      //      cerr << "Proportions of alleles at locus do not total to 1! Instead, they total to: " << tprop << endl;
    }
}		 

ostream &operator<<(ostream &stream, InfAlleleTbl &a)
{
  a.Write(stream);
  return stream;
}		 

istream &operator>>(istream &stream, InfAlleleTbl &a)
{
  a.Scan(stream);
  return stream;
}		 


/**

begin implementation of stepwise allele model

*/


StepAlleleTbl::StepAlleleTbl()
{
  clear();
  setMutationRate(0);
  UNUSED.reserve(500);
  maxstate=0;
  setClassType(STEPALLELETBL);
}
StepAlleleTbl::~StepAlleleTbl()
{

  clear();
  setMutationRate(0);
  UNUSED.resize(0);

#ifdef DEBUG
  cerr << "executing StepAlleleTbl destructor"<<endl;
#endif

  maxstate=0;
}

int StepAlleleTbl::mutator(int anum, int t)
{
  map<int, Allele, less<int> >::iterator tmpiter;
  Allele na;
  int newanum;

  assert(anum>=0);
  tmpiter=A.find(anum);
  if (tmpiter!=A.end())
    {
      if (RandLibObj.uniform()<rate)    //a mutation has occurred
	{
	  if ((*tmpiter).second.GetState()>0)
	    {
	      if (RandLibObj.uniform()>0.5)
		{
		  na.SetState((*tmpiter).second.GetState()+1);  
		}
	      else
		{
		  na.SetState((*tmpiter).second.GetState()-1);  
		}
	    }
	  else if ((*tmpiter).second.GetState()==0)
	    {
	      na.SetState((*tmpiter).second.GetState()+1);  
	    }
	  na.SetBirth(t);
	  na.SetFreq(1);
	  newanum = addAllele(na,t);
	  return newanum;
	}
      else
	{
	  (*tmpiter).second.SetFreq((*tmpiter).second.GetFreq()+1);  //update the allele freq without having to lookup allele by value
	  return anum;
	}
    }
  else
    {
#ifdef DEBUG
      cerr <<"Allele number "<<anum<<" not found in allele table: " <<endl<<*this<<endl;
#endif
      assert(tmpiter!=A.end());
      return -1;
    }
}

ostream &operator<<(ostream &stream, StepAlleleTbl &a)
{
  a.Write(stream);
  return stream;
}		 

istream &operator>>(istream &stream, StepAlleleTbl &a)
{
  a.Scan(stream);
  return stream;
}		 

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; minor-mode:  font-lock  ***
;;; End:  ***
*/












































