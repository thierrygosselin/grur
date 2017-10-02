/*

$Modified: astrand $

Copyright (C) 1999 Allan E. Strand

This file is part of Metasim

This is the declaration of the base object type.

*/


/*includes
*/

#include <SiteObj.h>
#include <RandLib.h>

///SiteObj constructor.  sets the state as it goes

SiteObj::SiteObj(char newstate)
{
      state = newstate;
}

SiteObj::~SiteObj()
{
}

char SiteObj::GetState()
{
  return state;
}

void SiteObj::SetState(char newstate)
{
  state=newstate;
}

void SiteObj::Mutate()
{
  double uni;
  char s1, s2, s3;

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
  if (uni < 0.33333)
    {
      state = s1;
    }
  else if (uni < 0.666667)
    {
      state = s2;
    }
  else 
    {
      state = s3;
    }
}

int SiteObj::operator==(SiteObj site)
{
  if (site.state == state) return 1;
  else return 0;
}

int SiteObj::operator!=(SiteObj site)
{
  if (site.state != state) return 1;
  else return 0;
}


// inserter for sites
ostream &operator<<(ostream &stream, SiteObj &site)
{
  stream << site.state;
  return stream;
}

// extractor for sites
istream &operator>>(istream & stream, SiteObj &site)
{
  char tmp;
  string alist = "agtcAGTC";
  stream >> tmp;
  if (alist.find(tmp)==string::npos)
    {
#ifdef DEBUG
      cerr << "Problem with DNA sequence.  A base other than AGTCagtc read (could signify premature sequence end"<<endl;
#endif
      assert(0==1);
    }
  site.state=tmp;
  return stream;
}


/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; End: ***
*/
