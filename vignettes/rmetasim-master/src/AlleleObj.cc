/*

$Modified: astrand $

Copyright (C) 1999 Allan E. Strand

This file is part of Metasim

This is the implementation of the allele object type.

*/

/*includes
*/
#include<RandLib.h>
#include <AlleleObj.h>

Allele::Allele (double p, int s, int b, int /*d*/)
{
  SetState(s);
  SetBirth(b);
  SetFreq(0);
  SetProp(p);
}

Allele::~Allele ()
{
}


void Allele::WriteState(ostream & stream)
{
  stream  <<state;
}
void Allele::Write(ostream & stream)
{
  stream  << prop <<" "<< birth <<" "<<state  << endl;
}

void Allele::Scan(istream & stream)
{
  char strn[100];
  stream >> prop >> birth >> state ;
  assert (prop >=0 && prop <= 1);
  stream.getline(&strn[0],100);
  freq = 0;
}

ostream & operator<<(ostream & stream, Allele &a)
{

  a.Write(stream);
  return stream;
}

istream & operator>>(istream & stream, Allele &a)
{
  a.Scan(stream);
  return stream;
}


SeqAllele::SeqAllele(int sl)
{
  //  dnaseq.clear();
  //  dnaseq.reserve(sl);
  dnaseq.resize(sl);
}

SeqAllele::~SeqAllele()
{
  
}


size_t SeqAllele::SeqLen()
{
  return dnaseq.size();
}

void SeqAllele::SetSite(char st, int sn)//st is the state, sn is the position
{
  dnaseq[sn]=st;
}

void SeqAllele::RandomSeq(double a, double c, double t, double g)
{
  int numa, numc, numt, numg, i;
  double bc;
  size_t sl;
  sl = SeqLen();
  i = 0;

  if (a+c+t <= 1.0)
    {
      g = 1.0 - (a+c+t) ; // if the sum of the other bases are less than one, set g = 1 - sum of others
      numa = 0; numc = 0; numt = 0; numg = 0;
      while (i<sl)
	{
	  bc = RandLibObj.uniform();
	  //	  if (a!=0.0 && (double(numa)/double(sl))<=a) 
	  if (bc<a)
	    {
	      dnaseq[i]='A';
	      numa++;
	      i++;
	    }
	  //	  if (c!=0.0 && (double(numc)/double(sl))<=c) 
	  if ((a<=bc)&(bc<(a+c))) 
	    {
	      dnaseq[i]='C';
	      numc++;
	      i++;
	    }
	    //	  if (t != 0.0 && (double(numt)/double(sl))<=t) 
	  if (((a+c)<=bc)&(bc<t)) 
	    {
	      dnaseq[i]='T';
	      numt++;
	      i++;
	    }
	    //	  if (g != 0 && (double(numg)/double(sl))<=g) 
	  if ((a+c+t)<=bc) 
	    {
	      dnaseq[i]='G';
	      numg++;
	      i++;
	    }
	}
    }
  else
    {
      // throw an exception 
    }
  random_shuffle(dnaseq.begin(), dnaseq.end(),randWrapper);
}

char SeqAllele::GetSite(int sn)
{
  return dnaseq[sn];
}

void SeqAllele::mutate()
{
  double uni;
  char s1, s2, s3;
  size_t sl = SeqLen(), b;
  if (state=='A')
    {
      s1 = 'G';
      s2 = 'C';
      s3 = 'T';
    }
  else if (state=='G')
    {
      s1 = 'A';
      s2 = 'C';
      s3 = 'T';
    }
  else if (state=='C')
    {
      s1 = 'G';
      s2 = 'A';
      s3 = 'T';
    }
  else
    {
      s1 = 'G';
      s2 = 'C';
      s3 = 'A';
    }

  uni = RandLibObj.uniform() ; 
  b = RandLibObj.unirange(sl - 1);

  //  Rprintf("sl is %d, ",sl);
  //   Rprintf("b is %d ",b);
  
  if (uni < 0.33333)
    {
      dnaseq[b] = s1;
    }
  else if (uni < 0.666667)
    {
      dnaseq[b] = s2;
    }
  else 
    {
      dnaseq[b] = s3;
    }
}

bool SeqAllele::Compare(SeqAllele obj) {
 
   size_t i;
   if (SeqLen() == obj.SeqLen()) 
      {
	for( i=0; i<SeqLen(); i++)
	  {
	    if (GetSite(i) != obj.GetSite(i))
	      {
		return false;
	      }
	  }
	return true;
      }
   else return false;
   }

void SeqAllele::WriteState(ostream & stream)
{
  char sob;
  for (size_t i=0; i<SeqLen();++i)
    {
      sob = GetSite(i);
      stream << sob;
    }
}

void SeqAllele::Write(ostream & stream)
{
  char sob;

  stream << prop <<"  "<< birth <<" ";
  for (size_t i=0; i<SeqLen();++i)
    {
      sob = GetSite(i);
      stream << sob;
    }
  stream << endl;
}

void SeqAllele::Scan(istream & stream)
{
  size_t i;
  char tmp;
  string alist = "agtcAGTC"; 
  stream >> prop >> birth ;
  assert (prop >=0 && prop <= 1);
  freq = 0;
  for (i=0;i<SeqLen();i++)
    {
      stream >> tmp;
      if (alist.find(tmp)==string::npos)
	{
#ifdef DEBUG
	  cerr << "Problem with DNA sequence.  A base other than AGTCagtc read (could signify premature sequence end"<<endl;
	  cerr << "Problem found at base: "<<i <<endl;
	  assert(0==1);
#endif
	}

      dnaseq[i]=tmp;

    }
}

bool operator==(SeqAllele allele1, SeqAllele allele2)
{
  return allele1.Compare(allele2);
}

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; minor-mode: font-lock ***
;;; End: ***
*/
