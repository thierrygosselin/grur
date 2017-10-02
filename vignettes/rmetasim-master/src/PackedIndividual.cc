/*

$Modified: astrand $

Copyright (C) 1999 Allan E. Strand

This file is part of Metasim

*/


/*includes
*/

#include <PackedIndividual.h>
#include <RandLib.h>

/**
Class PackedIndividual methods
 */

PackedIndividual::PackedIndividual (int c, int sx, int g, int nl)
{
  cl=c;
  sex=sx;
  gen=g;
  assert(nl <= MAXLOCI);
  ///  SetLoci();
  Change(g-1);
  SetLastRep(0);
  SetNumOff(0);
  id=0;mid=0;pid=0;
}

PackedIndividual::~PackedIndividual ()
{
}

/// from a uniform distribution
int PackedIndividual::RandomizeClass(int numclass)
{
  return RandLibObj.unirange(numclass);
}

void  PackedIndividual::SetClass(int numclass)
{
  cl = numclass;
}

int  PackedIndividual::GetClass()
{
  int c;
  c = int(cl);
  return c;
}

void PackedIndividual::SetSex(int newsex)
{
  sex = newsex;
}

int PackedIndividual::GetSex()
{
  return sex;
}

void PackedIndividual::SetGen(int newgen)
{
  gen = newgen;
}

int PackedIndividual::GetGen()
{
  return gen;
}


void PackedIndividual::SetLoci(AlleleLookTbl &Atbls)
{
  int i;
  for (i=0;i<MAXLOCI;i++)
    {
      if (i<int(Atbls.size()))
	{
	  PL[i]=Atbls[i]->getPloidy();
	}
      else
	{
	  PL[i]=0;
	}
    }
  nloc=int(Atbls.size());
}

void PackedIndividual::resetLoci(AlleleLookTbl &Atbls)
{
  int i,j;
  SetLoci(Atbls);
  for (i=0;i<MAXLOCI;i++)
    {
      for (j=0;j<MAXPLOIDY;j++)
	{
	  G[ ((i * MAXPLOIDY) + j) ]   = -1;
	}
    }
}


PackedIndividual PackedIndividual::MakeGamete(AlleleLookTbl &Atbls)
{
  int i;
  int lsize;
  //  int a=0,b=0,c=0;
  PackedIndividual pi;
  pi.resetLoci(Atbls);
  lsize = nloc;
  for (i=0;i<lsize;i++)
    {
      if (Atbls[i]->getTrans()==0)  //biparental inheritance
	{
	  assert(Atbls[i]->getPloidy()==2);
	  pi.G[((i * MAXPLOIDY))] = G [((i * MAXPLOIDY) + RandLibObj.unirange(1))];
	}
      else if (Atbls[i]->getTrans()==1 && GetSex()<2) //maternal inheritance.  female parent
	{
	  //	  b++;
	  pi.G[((i * MAXPLOIDY))] = G [ i * MAXPLOIDY ];
	}
      else if (Atbls[i]->getTrans()==2 && GetSex()>1) //Paternal inheritance. male parent
	{
	  //	  c++;
	  pi.G[((i * MAXPLOIDY))] = G [ i * MAXPLOIDY ]; 
	}
      else
	{
#ifdef DEBUG
	  cerr << "Fell through all of the inheritance types in MakeGamete " << endl;
#endif
	  assert(1==0);
	}
      if (!(pi.G[((i * MAXPLOIDY))]>-1))
	{
	  //	  cerr << "a "<<a<< " b "<<b<<" c "<<c<<endl;
	    assert(pi.G[((i * MAXPLOIDY))]>-1);
	}
    }
  return pi;
}

int PackedIndividual::GetRandAlleleIndex(int l)
{
  int index;
  double ru;

  ru=RandLibObj.uniform();
  if (ru==1) {ru=0.999999999999;} //if ru=1 then index will equal ploidy below, instead of ploidy-1
  index = (int)floor(ru * (PL[l]));

  assert(G [((l * MAXPLOIDY) + index)]>=0);

  return G [((l * MAXPLOIDY) + index)];
}


int PackedIndividual::IsGenotypeSet()
{
  int j,s,i;
  s=1;
  for (i=0;i<nloc;i++)
    {
      for (j=0;j<PL[i];j++)
	{
	  if(G[ ((i * MAXPLOIDY) + j) ]<0)
	    {
	      s=0;
	    }
	}
    }

  return s;

}

void PackedIndividual::SetRandGenotype(AlleleLookTbl &Atbls)
{
  int j,i;


  for (i=0;i<nloc;i++)
    {
      for (j=0;j<PL[i];j++)
	{
	  G[ ((i * MAXPLOIDY) + j) ]   =  Atbls[i]->getRandAlleleIndex();
	}
    }
}

PackedIndividual  PackedIndividual::repro_sex(PackedIndividual & SO1, PackedIndividual & SO2, int /*t*/, AlleleLookTbl &Atbls)
{
  int i;
  int k, l;
  PackedIndividual pi(SO1), ti0, ti1;
  pi.resetLoci(Atbls);

  ///  assert(SO2.IsGenotypeSet());

  ti0 = SO1.MakeGamete(Atbls);
  ti1 = SO2.MakeGamete(Atbls);

  for (i=0;i<nloc;i++)
    {
      k = ti0.GetAllele(i,0);
      assert(k>=0); 
      pi.G[((i * MAXPLOIDY) + 0)] = k;
      l = ti1.GetAllele(i,0);
      assert(l>=0); 
      pi.G[((i * MAXPLOIDY) + 1)] = l;

      ///swap alleles so that diploid heterozygotes are sorted
      if (PL[i]==2)
	{ 
	  if ((pi.G[((i * MAXPLOIDY) + 1)] >= 0) && (pi.G[((i * MAXPLOIDY) + 0)] > pi.G[((i * MAXPLOIDY) + 1)] ))
	    {
	      pi.swap_allele(i);
	    }
	}
    }
  return pi;
}



PackedIndividual  PackedIndividual::repro_asex(PackedIndividual & SO, int /*t*/)
{
  //  int i,j;

  //  int tmpi;

  //  int k;

  PackedIndividual pi(SO);
  //  pi.resetLoci(Atbls);

  return pi;
}

void PackedIndividual::Growth(AlleleLookTbl &Atbls)
{
  int j,i;
  for (i=0;i<nloc;i++)
    {
      for (j=0;j<PL[i];j++)
	{
	  Atbls[i]->AddAlleleFreq(G[ ((i * MAXPLOIDY) + j) ]);
	}
    }
}

void PackedIndividual::Birth(int t, AlleleLookTbl &Atbls)
{
  int i;
  //  cerr<< "in tmpI.Birth"<<endl;
  for (i=0;i<nloc;i++)
    {
      if (PL[i]==1)
	{
	  if (t>=0)
	    {
	      //	      cerr<<"running mutator hap"<<endl;
	      G[ ((i * MAXPLOIDY) + 0) ] = Atbls[i]->mutator(G[ ((i * MAXPLOIDY) + 0) ],t);
	      //cerr<<"mutator hap run"<<endl;
	    }
	  else
	    {
	      Atbls[i]->AddAlleleFreq(G[ ((i * MAXPLOIDY) + 0) ]);
	    }
	}
      if (PL[i]==2)
	{
	  if (t>=0)
	    {
	      //cerr<<"running mutator dip"<<endl;
	      G[ ((i * MAXPLOIDY) + 0) ] = Atbls[i]->mutator(G[ ((i * MAXPLOIDY) + 0) ],t);
	      //cerr<<"running mutator dip 1/2"<<endl;
	      G[ ((i * MAXPLOIDY) + 1) ] = Atbls[i]->mutator(G[ ((i * MAXPLOIDY) + 1) ],t);
	      //cerr<<"mutator dip run"<<endl;
	    }
	  else
	    {
	      Atbls[i]->AddAlleleFreq(G[ ((i * MAXPLOIDY) + 0) ]);
	      Atbls[i]->AddAlleleFreq(G[ ((i * MAXPLOIDY) + 1) ]);
	    }
	}
    }
  //cerr<< "leaving tmpI.Birth"<<endl;
}

void PackedIndividual::Death(int t, AlleleLookTbl &Atbls)
{
  size_t i, sz;
  sz = Atbls.size();
  for (i=0;i<sz;i++)
    {
      if (Atbls[i]->getPloidy()==1)
	{
	  Atbls[i]->KillAlleleCopy(G[ ((i * MAXPLOIDY) + 0) ],t);
	}
      if (Atbls[i]->getPloidy()==2)
	{
	  Atbls[i]->KillAlleleCopy(G[ ((i * MAXPLOIDY) + 0) ],t);
	  Atbls[i]->KillAlleleCopy(G[ ((i * MAXPLOIDY) + 1) ],t);
	}
    }
}


ostream & operator<<(ostream & stream, PackedIndividual &ind)
{
  int j,i,is;

  stream << ind.GetClass() << " " << ind.GetSex() << " " << ind.GetGen() << "  "<<ind.GetID()<<" "<<ind.GetMID()<<" "<<ind.GetPID()<<" ";

  for (j=0;j<ind.nloc;j++)
    {
      for (i=0;i<ind.PL[j]; i++)
	{
	  is=ind.GetAllele(j,i);
	  stream << is << " " ;
	}
      stream << "   " ;
    }
  stream << endl;
  return stream;
}

istream & operator>>(istream & stream, PackedIndividual &ind)
{
  int i, is, j;

  stream >> ind.cl >> ind.sex >> ind.gen >> ind.id >> ind.mid >> ind.pid;

  for (j=0;j<ind.nloc;j++)
    {
      for (i=0;i<ind.PL[j]; i++)
	{
	  stream >> is ;
	  ind.SetAllele(j,i,is);
	}
    }

  return stream;
}


/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; minor-mode: font-lock ***
;;; End: ***
*/
