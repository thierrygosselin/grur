/* This file is part of Metasim
   This file is the implementation of the landscape global Allele lookup table
*/

/* includes */
#include <AlleleTbl.h>


/*

   Begin implementation of allele table

 */

AlleleLookTbl::AlleleLookTbl()
{
}
AlleleLookTbl::~AlleleLookTbl()
{
  clear();
}

void AlleleLookTbl::push_back(AlleleTbl * atp)
{
  Atbl.push_back(atp);
}

void AlleleLookTbl::clear()
{
  int asz,i;
#ifdef DEBUG
  cerr << "deleting Atbl[i]" <<endl;
  cerr << "ntbls"<<Atbl.size()<<endl;
#endif
  asz = size();
  if (asz>0)
    {
      for (i=0;i<asz;i++)
	{
#ifdef DEBUG
	  cerr << "Cleaning  Atbl[i] i= "<<i <<endl;
#endif
	  delete Atbl[i];
	}
      Atbl.resize(0);
    }
}


void AlleleLookTbl::DummyFreq(int ps)
{
  size_t i;
  for (i=0;i<size();i++)
    Atbl[i]->dummyfreq(ps);
}

void AlleleLookTbl::ZeroFreq()
{
  size_t i;
  for (i=0;i<Atbl.size();i++)
    Atbl[i]->zerofreq();
}

ostream &operator<<(ostream & stream, AlleleLookTbl &a)
{
  int i ;
  for (i=0;i<int(a.size());i++)
    {
      a[i]->Write(stream);
    }
  return stream;  
}


istream &operator>>(istream & stream, AlleleLookTbl &a)
{
  AlleleTbl * AT;
  AT =NULL;
  int nloc,i,loctype;

  stream >> nloc;
  ///		    Atbls.reserve(l.nloc);
  for (i=0;i<nloc;i++)
    {
      stream >> loctype;
      if (loctype==0)
	{
	  AT = new InfAlleleTbl;
	}
      else if (loctype==1)
	{
	  AT = new StepAlleleTbl;
	}
      else if (loctype==2)
	{
	  AT = new SeqAlleleTbl;
	}
      else
	{
#ifdef DEBUG
	  cerr << "Cannot identify the locus type in input stream"<<endl;
#endif
	  assert (0==1);
	}
      AT->clear();
      AT->Scan(stream);
      a.push_back(AT);
    }
  return stream;
}

//AlleleLookTbl Atbls;

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; minor-mode:  font-lock  ***
;;; End:  ***
*/




