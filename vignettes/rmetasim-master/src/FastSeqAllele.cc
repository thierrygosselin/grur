/* This file is part of Metasim
   This file is the implementation of the SeqAllele lookup table
*/

/* includes */
#include <FastSeqAllele.h>


/*

   Begin implementation of Sequence allele table

 */


///t is the current clock-tick
SeqAlleleTbl::SeqAlleleTbl()
{
  clear();
  setMutationRate(0);
  UNUSED.reserve(500);
  setClassType(SEQALLELETBL);
}
SeqAlleleTbl::~SeqAlleleTbl()
{
  clear();
  setMutationRate(0);
#ifdef DEBUG
  cerr << "executing SeqAlleleTbl destructor"<<endl;
#endif
  UNUSED.resize(0);
}


SeqAllele SeqAlleleTbl::getRandAllele()
{
  return getAllele(getRandAlleleIndex());
}

/// returns an allele index at random
int SeqAlleleTbl::getRandAlleleIndex()
{
  int sz,i,tofind;
  sz = A.size();
  double *p = new double[sz];
  int *lookup = new int[sz];

  map<int, SeqAllele, less<int> >::iterator tmpiter;

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

int SeqAlleleTbl::AlleleTotalCnt()
{
  int tot=0;
  map<int, SeqAllele, less<int> >::iterator tmpiter;
  for (tmpiter=A.begin();tmpiter!=A.end();tmpiter++)
    {
      tot= tot + (*tmpiter).second.GetFreq();
    }
  return tot;
}
void SeqAlleleTbl::dummyfreq(int ps)
{
  map<int, SeqAllele, less<int> >::iterator tmpiter;
  for (tmpiter=A.begin();tmpiter!=A.end();tmpiter++)
    {
      (*tmpiter).second.SetFreq(int(ceil(ps * (*tmpiter).second.GetProp())));
    }
}

void SeqAlleleTbl::zerofreq()
{
  map<int, SeqAllele, less<int> >::iterator tmpiter;
  for (tmpiter=A.begin();tmpiter!=A.end();tmpiter++)
    {
      (*tmpiter).second.SetFreq(0);
    }
}

void SeqAlleleTbl::CalcProps()
{
  int i,sz;
  double np,tot;
  map<int, SeqAllele, less<int> >::iterator tmpiter;
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
int SeqAlleleTbl::addAllele(SeqAllele  na, int t)
{
  int f;
  int anum;
  SeqAllele a;
  map<int, SeqAllele, less<int> >::iterator tmpiter;

  anum=-1;
  f=0;
  tmpiter=A.begin();
  while (((f==0)&&(tmpiter!=A.end()))&&(A.size()>0))
    {
      if ((*tmpiter).second==na)
	{
	  f=1; //found = true
	  (*tmpiter).second.SetFreq((*tmpiter).second.GetFreq()+1); // always keep count of alleles
	  anum=(*tmpiter).first;
	}
      tmpiter++;
    }
  if (f==0)
    {
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


int SeqAlleleTbl::addAlleleAndIndex(SeqAllele  na, int ai)
{
  int f;
  int anum;
  SeqAllele a;
  map<int, SeqAllele, less<int> >::iterator tmpiter;

  anum=-1;
  f=0;
  tmpiter=A.begin();
  while (((f==0)&&(tmpiter!=A.end()))&&(A.size()>0))
    {
      if ((*tmpiter).second==na)
	{
	  if (ai==(*tmpiter).first)
	    {
#ifdef DEBUG
	      cerr << "Allele index: "<<ai<< " already present in table " <<endl;
#endif
	      assert(0==1);
	    }
	}
      tmpiter++;
    }
  if (f==0)
    {
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
int SeqAlleleTbl::addAlleleState(SequenceType seq, int t)
{
  SeqAllele a(seq.size());
  a.SetState(seq);
  return addAllele(a,t);
}


///Go through the allele table and remove those alleles with freq=0
void SeqAlleleTbl::GCAlleles()
{
  map<int, SeqAllele, less<int> >::iterator tmpiter;
  for (tmpiter = A.begin();tmpiter!=A.end();tmpiter++)
    {
      if ((*tmpiter).second.GetFreq()<=0)
	{
	  UNUSED.push_back((*tmpiter).first); //keep the index to resuse later
	  A.erase(tmpiter); // erase the allele
	}
    }
}


void SeqAlleleTbl::KillAlleleCopy(int i, int /*t*/)
{
  map<int, SeqAllele, less<int> >::iterator tmpiter;
  tmpiter = A.find(i);
  if (tmpiter!=A.end())
    {
      (*tmpiter).second.SetFreq((*tmpiter).second.GetFreq() - 1);
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


void SeqAlleleTbl::clear()
{
  ploidy=0;
  trans=0;
  rate=0.0;
  A.clear();
  UNUSED.resize(0);
}

int SeqAlleleTbl::mutator(int anum, int t)
{
  map<int, SeqAllele, less<int> >::iterator tmpiter;
  int newanum;
  assert(anum>=0);

  if (RandLibObj.uniform()<rate)    //a mutation has occurred
    {
      SeqAllele na;
      na = getAllele(anum);
      
      na.mutate();  //sets state to a previously unused value.
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
	  (*tmpiter).second.SetFreq((*tmpiter).second.GetFreq()+1);  //update the allele freq without having to lookup allele by value
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
void SeqAlleleTbl::WriteAlleleState(int a, ostream &stream)
{
  SeqAllele al;
  al = getAllele(a);
  al.WriteState(stream);
}

void SeqAlleleTbl::Write(ostream &stream)
{
  int sz;
  int i;
  int ai;
  map<int, SeqAllele, less<int> >::iterator tmpiter;
  SeqAllele al;

  //  CalcProps();
  //  GCAlleles();


  sz=A.size();
  stream << sz << endl;
  stream << seqlen << endl;
  stream << rate << endl;
  stream << ploidy << endl;
  stream << trans << endl;
  
  i=0;

  for (tmpiter=A.begin();tmpiter!=A.end();tmpiter++)
    {
      ai = (*tmpiter).first;
      al = (*tmpiter).second;
      i++;
      //      if ((*tmpiter).second.GetFreq()>0)
      //	{
      stream << ai << "  " << al;
      //	}
      //      else
      //	{
      //	  cerr <<"Not written: Allele index ai: "<<ai <<" locus: "<<al<<endl;
      //	}
    }
  stream << endl;

}		 

vector<int>  SeqAlleleTbl::getAindices()
{


  vector<int> aindices;
  map<int, SeqAllele, less<int> >::iterator tmpiter;
  //sz=A.size();
  for (tmpiter=A.begin();tmpiter!=A.end();tmpiter++)
    {
	  aindices.push_back((*tmpiter).first);
    }
  return aindices;
}		 


void SeqAlleleTbl::Scan(istream &stream)
{
  int i;
  int ai;
  //  int tmp;
  int numa;

  double tprop;
  tprop = 0;

  clear();

  stream >> numa;
  stream >> seqlen;

  SeqAllele newa(seqlen);

  newa.SeqLen();

  stream >> rate;
  stream >> ploidy;
  stream >> trans;
  for (i=0;i<numa;i++)
    {
      stream >> ai;
      newa.Scan(stream);
      addAlleleAndIndex(newa,ai);
      tprop = newa.GetProp() + tprop;
    }
  if (tprop != 1.0) //primitive error checking
    {
#ifdef DEBUG
      //      cerr << "Proportions of alleles at locus do not total to 1! Instead, they total to: " << tprop << endl;
#endif
    }
}		 

ostream &operator<<(ostream &stream, SeqAlleleTbl &a)
{
  a.Write(stream);
  return stream;
}		 

istream &operator>>(istream &stream, SeqAlleleTbl &a)
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
