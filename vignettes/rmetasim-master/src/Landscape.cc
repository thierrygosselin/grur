/*

$Modified: astrand $

Copyright (C) 1999-2006 Allan E. Strand

This file is part of Metasim
*/

/*includes
*/

#include <Landscape.h>
#include <sstream>
#include <unistd.h>

using namespace std;


///LocalMat class

void LocalMat::SetSize(size_t sz)
{
#ifdef DEBUG
  cerr << "SetSize Slocal" <<endl;
#endif

  Slocal.SetSize(sz);

#ifdef DEBUG
  cerr << "SetSize Rlocal" <<endl;
#endif

  Rlocal.SetSize(sz);

#ifdef DEBUG
  cerr << "SetSize Mlocal" <<endl;
#endif

  Mlocal.SetSize(sz);
}



ostream &operator<<(ostream & stream, LocalMat &lm)
{
  stream << lm.Slocal << endl << lm.Rlocal << endl << lm.Mlocal << endl;
  return stream;
}
istream &operator>>(istream & stream, LocalMat &lm)
{
  stream >> lm.Slocal ;
  stream >> lm.Rlocal ;
  stream >> lm.Mlocal ;
  return stream;
}




///end LocalMat class


///Landscape class

///Constructor
Landscape::Landscape (int /*h*/, int /*stg*/, int /*loc*/, int /*ep*/, int /*nd*/, int /*gn*/)
{
  ndemo=1;
#ifdef DEBUG
  cerr << "Constructing landscape object" <<endl;
#endif
  ///  init(h, stg, loc, ep, nd, gn);
}

///Destructor

Landscape::~Landscape()
{
#ifdef DEBUG
  cerr << "Destructing landscape object" <<endl;
#endif
  /*
  for (int i=0;i < s*nhab; i++)
    {
      I[i].ClearClass(0);
    }


#ifdef DEBUG
	  cerr << "Landscape destructor: Deleting vectors and mats i= "<<i <<endl;
#endif
	  S.resize(0);
	  R.resize(0);
	  M.resize(0);
	  for (i=0;i<nep;i++)
	    {
	      demoProbVec[i].resize(0);
	    }
	  evec.resize(0);
	  kvec.resize(0);
	  LM.resize(0);
	  demoProbVec.resize(0);
  */
#ifdef DEBUG
  cerr << "Landscape destructor exiting" <<endl;
#endif

}

/**
   set the appropriate matrices
*/


void Landscape::setepochs(int ep)
{
  int i;
  
  nep=ep;
  epochs.resize(nep);
  epochprobs.resize(nep);
  S.resize(nep);
  R.resize(nep);
  M.resize(nep);
  evec.resize(nep);
  kvec.resize(nep);
  demoProbVec.resize(nep);
  for (i=0;i<nep;i++)
    {
      S[i].SetSize(s*nhab)      ;
      R[i].SetSize(s*nhab)      ;
      M[i].SetSize(s*nhab)      ;
      if (int(demoProbVec[i].size())!=ndemo)
	{
#ifdef DEBUG
	  cerr << "sizing i of nep demoProbVec[i] to ndemo.  i="<<i<<" ndemo="<<ndemo <<endl;
#endif
	  demoProbVec[i].resize(ndemo);
	}
#ifdef DEBUG
  cerr << "done sizing the demoProbVecs" <<endl;
#endif
      evec[i].resize(nhab);
      kvec[i].resize(nhab);
    }
}


void Landscape::setndemo(int nd)
{
  int i;
  ndemo=nd;
#ifdef DEBUG
  cerr << "REserving space for LM" <<endl;
#endif
  LM.resize(ndemo);
#ifdef DEBUG
  cerr << "done REserving space for LM" <<endl;
#endif

#ifdef DEBUG
  cerr << "Setting sizes of LM atrices to: "<<s <<endl;
#endif
  for (i=0;i<ndemo;i++)
    {
      LM[i].SetSize(s);
    }
///KKM 6.2.05..................................................................
#ifdef DEBUG
  cerr << "REserving space for LMK" <<endl;
#endif
  LMK.resize(ndemo);
#ifdef DEBUG
  cerr << "done REserving space for LMK" <<endl;
#endif

#ifdef DEBUG
  cerr << "Setting sizes of LMK matrices to: "<<s <<endl;
#endif
  for (i=0;i<ndemo;i++)
    {
      LMK[i].SetSize(s);
    }
///............................................................................
#ifdef DEBUG
  cerr << "Resizing the demoProbVecs to ndemo="<<ndemo <<endl;
#endif
  for (i=0;i<nep;i++)
    {
      if (int(demoProbVec[i].size())!=ndemo)
	{
	  demoProbVec[i].resize(ndemo);
	}
    }
#ifdef DEBUG
  cerr << "About to exit setndemo" <<endl;
#endif

}


void Landscape::sethabs(int h) 
{
  nhab=h;
}

void Landscape::setstages(int stg) { s=stg; }
void Landscape::setxdim(int x) { xdim=x; }
void Landscape::setydim(int y) { ydim=y; }
void Landscape::setepochprob(int ce, double prob) {epochprobs[ce]=prob;}
void Landscape::setepochstart(int ce, int strt) {epochs[ce]=strt;}

void Landscape::init(int h, int stg, int /*loc*/, int ep, int nd, int gn)
{
#ifdef DEBUG
  cerr << "waiting for return to cont";
  cerr << endl;
  cerr << "Running: sethabs(h)" <<endl;
#endif
  sethabs(h);
#ifdef DEBUG
  cerr << "Running: setstages(stg)" <<endl;
#endif
  setstages(stg);
#ifdef DEBUG
  cerr << "Running: setxdim()" <<endl;
#endif
  setxdim();
#ifdef DEBUG
  cerr << "Running: setydim()" <<endl;
#endif
  setydim();
#ifdef DEBUG
  cerr << "Running: setepochs(ep)" <<endl;
#endif
  setepochs(ep);
#ifdef DEBUG
  cerr << "Running: setndemo(nd)" <<endl;
#endif
  setndemo(nd);
#ifdef DEBUG
  cerr << "Running: setgens(gn)" <<endl;
#endif
  setgens(gn);
#ifdef DEBUG
  cerr << "Running: setMaxLandSize()" <<endl;
#endif
  setMaxLandSize();
#ifdef DEBUG
  cerr << "Running: setself()" <<endl;
#endif
  setself();

#ifdef DEBUG
  cerr << "Running: unsetRandEpoch()" <<endl;
#endif

  unsetRandEpoch();

#ifdef DEBUG
  cerr << "done running: unsetRandEpoch()" <<endl;
#endif
  t=0;
  e=0;

  title = "";

#ifdef DEBUG
  cerr << "Running:   I.resize(s * nhab)" <<endl;
#endif
  I.resize(s * nhab);

  multiple_paternity = 1 ;

  ///HABDELTA is precision when inplementing carrying capacity via Carry()
  habdelta = HABDELTA;

  setnextID(1);

#ifdef DEBUG
  cerr << "end of init" <<endl;
#endif

}

void Landscape::setS(TransMat a, int ep)
{
  int q;
  if (ep<0) 
    {
      for (q=0;q<nep;q++)
	{
	  S[q].SetElement(0,0, 0.4); S[q].SetElement(1,0, 0.0); S[q].SetElement(2,0, 0.0); S[q].SetElement(3,0, 0.0);
	  S[q].SetElement(0,1, 0.3); S[q].SetElement(1,1, 0.6); S[q].SetElement(2,1, 0.0); S[q].SetElement(3,1, 0.2);
	  S[q].SetElement(0,2, 0.0); S[q].SetElement(1,2, 0.0); S[q].SetElement(2,2, 0.4); S[q].SetElement(3,2, 0.0);
	  S[q].SetElement(0,3, 0.0); S[q].SetElement(1,3, 0.1); S[q].SetElement(2,3, 0.3); S[q].SetElement(3,3, 0.6);
	}
    }
  else
    {
      S[ep].SetMat(a);
    }
}

void Landscape::setR(TransMat a, int ep)
{
  int q;
  if (ep<0) 
    {
      for (q=0;q<nep;q++)
	{
      R[q].SetElement(0,0, 0.0); R[q].SetElement(1,0, 5.0); R[q].SetElement(2,0, 0.0); R[q].SetElement(3,0, 0.0);
      R[q].SetElement(0,1, 0.0); R[q].SetElement(1,1, 0.0); R[q].SetElement(2,1, 0.0); R[q].SetElement(3,1, 0.0);
      R[q].SetElement(0,2, 0.0); R[q].SetElement(1,2, 0.0); R[q].SetElement(2,2, 0.0); R[q].SetElement(3,2, 3.0);
      R[q].SetElement(0,3, 0.0); R[q].SetElement(1,3, 0.0); R[q].SetElement(2,3, 0.0); R[q].SetElement(3,3, 0.0);
	}
    }
  else
    {
      R[ep].SetMat(a);
    }
}

void Landscape::setM(TransMat a, int ep)
{
  int q;
  if (ep<0) 
    {
      for (q=0;q<nep;q++)
	{
      M[q].SetElement(0,0, 0.0); M[q].SetElement(1,0, 5.0); M[q].SetElement(2,0, 0.0); M[q].SetElement(3,0, 0.0);
      M[q].SetElement(0,1, 0.0); M[q].SetElement(1,1, 0.0); M[q].SetElement(2,1, 0.0); M[q].SetElement(3,1, 0.0);
      M[q].SetElement(0,2, 0.0); M[q].SetElement(1,2, 0.0); M[q].SetElement(2,2, 0.0); M[q].SetElement(3,2, 3.0);
      M[q].SetElement(0,3, 0.0); M[q].SetElement(1,3, 0.0); M[q].SetElement(2,3, 0.0); M[q].SetElement(3,3, 0.0);
	}
    }
  else
    {
      M[ep].SetMat(a);
    }
}


void Landscape::setextinct(int ep, double *ev)
{
  int i;
    for (i=0;i<nhab;i++)
    {
      evec[ep][i]=ev[i];
    }
}

void Landscape::setk(int ep, int *cv)
{
  int i;
    for (i=0;i<nhab;i++)
    {
      kvec[ep][i]=cv[i];
    }
}

void Landscape::getextinct(int ep, double *ev)
{
  int i;
    for (i=0;i<nhab;i++)
    {
      ev[i]=evec[ep][i];
    }
}
 
void Landscape::getk(int ep, int *cv)
{
  int i;
    for (i=0;i<nhab;i++)
    {
      cv[i]=kvec[ep][i];
    }
}


void Landscape::setldemovector(int ep, double *dv)
{
  int i;
  for (i=0;i<ndemo;i++)
    {
      demoProbVec[ep][i]=dv[i];
    }
}

void Landscape::getldemovector(int ep, double *dv)
{
  int i;
    for (i=0;i<ndemo;i++)
    {
      dv[i]=demoProbVec[ep][i];
    }
}

void Landscape::zeroextinct()
{
  int q,i;
  for (q=0;q<nep;q++)
    for (i=0;i<nhab;i++)
    {
      evec[q][i]=0;
    }
}
 
void Landscape::zerok()
{
  int q,i;
  for (q=0;q<nep;q++)
    for (i=0;i<nhab;i++)
    {
      kvec[q][i]=0;
    }
}

void Landscape::ChooseEpoch()
{
  int i;
  if (randepoch)
    {
      RandomlyChooseEpoch();
    }
  else
    {
      for (i=0;i<nep;i++)
	{
	  if (epochs[i]<=t) 
	    {
	      e=i;
	    }
	}
    }
}

void Landscape::RandomlyChooseEpoch()
 {
   int i;
   if (randepoch>0)
    {
      double *p = new double[nep];

      for (i=0;i<nep;i++)
	{
	  p[i]=epochprobs[i];
	}
      e = RandLibObj.multinomial(p,nep);

      delete[] p;
    }
}

void Landscape::SequentiallyConstructDemoMatrix()
{

  int fr,to;
  int i,rm;
  int newto,newfr;

  rm=0;
  for (i=0;i<nhab;i++)
    {
      if (rm>=ndemo)
	{
	  rm=0;
	}
      for (fr=0;fr<s;fr++)
	{
	  for (to=0;to<s;to++)
	    {
	      newto = (s*i)+to ;
	      newfr = (s*i)+fr ;
	      S[e].SetElement(newfr,newto,LM[rm].GetSlocalVal(fr,to));
	      R[e].SetElement(newfr,newto,LM[rm].GetRlocalVal(fr,to));
	      M[e].SetElement(newfr,newto,LM[rm].GetMlocalVal(fr,to));
	    }
	}
      rm++;
    }
}
///KKM 6.6.05...............................................................
void Landscape::SequentialDensityDependentDemoMatrix()
{

  int fr,to;
  int i,rm;
  int newto,newfr;
  double ValZero,ValK,NewVal;

  rm=0;
  for (i=0;i<nhab;i++)
    {
      if (rm>=ndemo)
	{
	  rm=0;
	}
      for (fr=0;fr<s;fr++)
	{
	  for (to=0;to<s;to++)
	    {
	      newto = (s*i)+to ;
	      newfr = (s*i)+fr ;
	      ValZero = LM[rm].GetSlocalVal(fr,to);
	      ValK = LMK[rm].GetSlocalVal(fr,to);
	      NewVal = (ValK-ValZero)*(double (PopSize(i))/double(kvec[e][i]))+ValZero;
	      S[e].SetElement(newfr,newto,NewVal);
	      ValZero = LM[rm].GetRlocalVal(fr,to);
	      ValK = LMK[rm].GetRlocalVal(fr,to);
	      NewVal = (ValK-ValZero)*(double (PopSize(i))/double(kvec[e][i]))+ValZero;
	      R[e].SetElement(newfr,newto,NewVal);
	      ValZero = LM[rm].GetMlocalVal(fr,to);
	      ValK = LMK[rm].GetMlocalVal(fr,to);
	      NewVal = (ValK-ValZero)*(double (PopSize(i))/double(kvec[e][i]))+ValZero;
	      M[e].SetElement(newfr,newto,NewVal);
	    }
	}
      rm++;
    }
}
///........................................................................
void Landscape::RandomlyConstructDemoMatrix()
{
  double *p = new double[ndemo];
  int fr,to;
  int i,rm;
  int newto,newfr;

  //set the probs of the multinomial distribution to pass to the rng
  for (i=0;i<ndemo;i++)
    {
      p[i]=demoProbVec[e][i];
    }
  
  for (i=0;i<nhab;i++)
    {
      rm = RandLibObj.multinomial(p,ndemo);
      for (fr=0;fr<s;fr++)
	{
	  for (to=0;to<s;to++)
	    {
	      newto = (s*i)+to ;
	      newfr = (s*i)+fr ;
	      S[e].SetElement(newfr,newto,LM[rm].GetSlocalVal(fr,to));
	      R[e].SetElement(newfr,newto,LM[rm].GetRlocalVal(fr,to));
	      M[e].SetElement(newfr,newto,LM[rm].GetMlocalVal(fr,to));
	    }
	}
    }

  delete[] p;
}

///KKM 6.7.05..................................................................
void Landscape::RandomDensityDependentDemoMatrix()
{
  double *p = new double[ndemo];
  int fr,to;
  int i,rm;
  int newto,newfr;
  double ValZero,ValK,NewVal;

  //set the probs of the multinomial distribution to pass to the rng
  for (i=0;i<ndemo;i++)
    {
      p[i]=demoProbVec[e][i];
    }
  
  for (i=0;i<nhab;i++)
    {
      rm = RandLibObj.multinomial(p,ndemo);
      for (fr=0;fr<s;fr++)
	{
	  for (to=0;to<s;to++)
	    {
	      newto = (s*i)+to ;
	      newfr = (s*i)+fr ;
	      ValZero = LM[rm].GetSlocalVal(fr,to);
	      ValK = LMK[rm].GetSlocalVal(fr,to);
	      NewVal = (ValK-ValZero)*(double (PopSize(i))/double(kvec[e][i]))+ValZero;
	      S[e].SetElement(newfr,newto,NewVal);
	      ValZero = LM[rm].GetRlocalVal(fr,to);
	      ValK = LMK[rm].GetRlocalVal(fr,to);
	      NewVal = (ValK-ValZero)*(double (PopSize(i))/double(kvec[e][i]))+ValZero;
	      R[e].SetElement(newfr,newto,NewVal);
	      ValZero = LM[rm].GetMlocalVal(fr,to);
	      ValK = LMK[rm].GetMlocalVal(fr,to);
	      NewVal = (ValK-ValZero)*(double (PopSize(i))/double(kvec[e][i]))+ValZero;
	      M[e].SetElement(newfr,newto,NewVal);
	    }
	}
    }

  delete[] p;
}
///............................................................................

void Landscape::popsizeset(std::vector<int> &ps)
{
  int i,j,psz;
  int totpop ;
  PackedIndividual Ind;
  DemoClass DC;

  psz=ps.size();

  totpop = 0;

  nextID = 1;

  for (i=0; i<psz; i++)
    {
      totpop = totpop + ps[i];
      I.push_back(DC);
    }

  for (i=0; i<psz; i++)
    {
      I[i].SetClass(i);
      for (j=0; j<ps[i]; j++)
	{
	  Ind.SetClass(i);
	  Ind.SetLoci(Atbls);
	  Ind.SetRandGenotype(Atbls);
	  Ind.Change(-1);
	  Ind.SetLastRep(-1);
	  Ind.SetID(nextID);
	  nextID=nextID+1;
	  Ind.SetMID(0);
	  Ind.SetPID(0);
	  Ind.SetNumOff(0);
	  Ind.Birth(-1,Atbls);
	  I[i].AddIndividual(Ind);
	}
    }
}


int Landscape::Habitat(int stage)
{
  int retv;
  double st;
  double ns;

  st = double(stage)*1.0;
  ns = double(s)*1.0;
  retv = int(floor(st/ns));
  return retv;
}


int Landscape::PopSize(int p)
{
  int i,tot;
  int sz = nhab*s;

  tot=0;

  if (p!=-1)
    {
      for (i=0;i<sz;i++)
	{
	  if (Habitat(i)==p)
	    {
	      tot = tot + I[i].size();
	    }
	}
    }
  else
    {
      for (i=0;i<sz;i++)
	{
	  tot = tot + I[i].size();
	}
    }
  return tot;
}



//void Landscape::Survive(double eigenratio)
void Landscape::Survive()
{
  PackedIndividual ind,tmpind;
  vector < int > deadindices, changeindices;
  vector < int >::iterator inditer;
  int indx;
  int rs;
  size_t i, j, isz, sz;


  deadindices.reserve(1000);
  changeindices.reserve(1000);
  sz = nhab * s;
  for (i=0;i<sz;i++)
    {
      S[e].SetFromState(i); //choose a column in the survival/migration matrix
      //S[e].SetRandomToStateVec(eigenratio);
      S[e].SetRandomToStateVec(1);
      I[i].ResetIndividuals();
      isz = I[i].size();
      I[i].ResetIndividuals();
      for (j=0;j<isz;j++)
	{
	  ind = I[i].GetCurrentIndividual();
	  indx = I[i].GetCurrentIndex();
	  if ((indx<0)||(ind.GetClass()<0))
	    {
#ifdef DEBUG
	      cerr << " run off the the end of the individual map for class " << i<<endl;
#endif
	      assert(ind.GetClass()>=0);
	    }
	  if (ind.GetChanged()<t)
	    {
	      rs = S[e].RandomState();
	      if (rs<0)//ind dies
		{
		  deadindices.push_back(indx);
		}
	      else if (rs!=int(i))
		{
		  ind.Change(t);
		  ind.SetClass(rs);
		  ind.Growth(Atbls);
		  I[rs].AddIndividual(ind);
		  changeindices.push_back(indx);
		}
	      else
		{
		  I[i].ChangeInd(indx,t);
		}
	    }
	  else
	    {
	    }
	  if (I[i].NextIndividual()) //advance the individual pointer
	    {
	      break;
	    }
	}

      for (inditer=deadindices.begin();inditer!=deadindices.end();inditer++)
	{
	  I[i].RemoveInd(*inditer,t,Atbls);
	}
      for (inditer=changeindices.begin();inditer!=changeindices.end();inditer++)
	{
	  I[i].RemoveInd(*inditer,t,Atbls);
	}

      deadindices.clear();
      changeindices.clear();
      RandLibObj.FreeDiscreteLookup();
    }
}
/** 

    This method returns a pointer to an array of double "p"  that is always
   the length of th enumber of demographic classes.  It is filled with
   type double 

*/
//int Landscape::CalculateMaleGameteClassVector(PackedIndividual pi)
int Landscape::CalculateMaleGameteClassVector(int k)
{

  int i,sz;
  sz = nhab * s;

  double *n = new double[sz];
  double *p = new double[sz];
  double clsz = 0.0, tmp, tmpct, tot=0;  

  M[e].SetToState(k);

  ///  cerr <<"M[e].GetToState "<< M[e].GetToState()<<endl;
  ///multiply the individual class sizes times the m's
  for (i=0;i<sz;i++)
    { 
      M[e].SetFromState(i);
///      cerr <<"M[e].GetFromState "<< M[e].GetFromState()<<endl;

      tmp =  M[e].Value();
      tmpct = double(I[i].size());
      n[i] = ( tmp * tmpct );
      clsz += n[i];

///      cerr << "tmp, tmpct, clsz :"<<tmp<<","<<tmpct<<","<<clsz<<endl;
    }
  ///calculate the relative probabilities of getting a gamete from a class

  if (clsz>0)
    {
      for (i=0;i<sz;i++)
	{
	  p[i] = (n[i] / clsz );
	  tot += p[i];
	}
      if (tot>1)
	{
	  if (tot > 1.1) //something is very wacky and the program should terminate
	    {
#ifdef DEBUG
	      cerr << "the probabilities of choosing a class total to more than 1: total = "<<tot<<endl;
#endif
	      assert (tot<=1);
	    }
	  else
	    {
	      for (i=0;i<sz;i++)
		{
		  p[i] = p[i]/tot;
		}
	    }
	}
      RandLibObj.SetDiscreteLookup(p,sz);

      delete[] p;
      delete[] n;

      return 1;
    }
  else 
    {
      delete[] p;
      delete[] n;

      RandLibObj.FreeDiscreteLookup();

      return 0;
      //      cerr << "no individuals in any of the donating classes in CalculateMaleGameteClassVector"<<endl;
    }
}

PackedIndividual Landscape::FindMate(PackedIndividual /*pi*/)
{
  PackedIndividual tmpI;

  int mc;

  mc = RandLibObj.PickMultinomial();
  //  tmpI.SetClass(-1);

  //  cerr << "mc "<<mc<<endl;
  //  int sz = I[mc].size();
  //  cerr << "length I[mc] "<< sz <<endl;

  tmpI = I[mc].GetRandomInd();

  //  tmpI = I[mc].GetIndividual(0);
  /*
      
  assert(tmpI.GetClass()>=0);
  assert(tmpI.IsGenotypeSet());

  if (tmpI.GetClass()<0)
    {
      cerr << "no mate found in FindMate, in Landscape.cc " << endl;
      assert(0==1);
    }
  */
  return tmpI;
}

void Landscape::testfindmate(PackedIndividual pi)
{
  PackedIndividual mate;
  size_t i;
#ifdef DEBUG
  cerr << "target ind: " <<pi<<endl;
#endif
  for (i=1;i<100;i++)
    {
      mate = FindMate(pi);
#ifdef DEBUG
      cerr <<"potential mate # "<<i<<" :  "<< mate<<endl;
#endif
    }
}

/**

    This would be the method to override if you were to add a feedback
    between genotype and offspring production/dispersal
    characteristics.  By override, I mean make a class that inherits
    everything from class Landscape_statistics (Landscape would work
    too).  All you need to do is modify Reproduce(), everything else
    would be inherited from the parent class.  Please don't modify
    Landscape or LAndscape_statistics, because there are software that
    depend upon it.

    There are some verbose notes below on how to convert for selection.

 */
//void Landscape::Reproduce(double eigenratio)
void Landscape::Reproduce()
{
  PackedIndividual tmpI, mate, searchI;
  vector < double > pvec;
  //  int err;
  //  int indx;
  int q,noff ;
  size_t j, k, l, sz, lsz ;
  sz = nhab * s;
  CompressInd();
  
  //cerr<<"Entering L.Reproduce"<<endl;

  for (k=0;k<sz;k++)
    {
      //      cerr << "Trying to reproduce from class: "<<k<<endl;
      if (R[e].AnyFrom(k)) ///find out if offspring can be produced by this class
	{
	  //cerr<<"L.Reproduce k:"<< k <<endl;
  
	  R[e].SetFromState(k);
	  I[k].CompressClass(0.5);///save space (may help speed )
	  I[k].ResetIndividuals();///set an internal pointer to I[k] first ind in list
	  
	  lsz=I[k].size();
	  if (CalculateMaleGameteClassVector(k))
	    {  
	      for (l=0;l<lsz;l++)
		{
		  searchI = I[k].GetCurrentIndividual();
		  if (searchI.GetClass()<0)
		    {
#ifdef DEBUG
		      cerr << "no individual returned from deomgraphic class"<<endl;
#endif
		      assert(searchI.GetClass()==0);
		    }
		  //		  indx = I[k].GetCurrentIndex();

		  //iterate through mothers
		  
                  for (j=0;j<sz;j++)
		    {
		      R[e].SetToState(j);
		      ///pick a number of offspring from a Poisson dist with mean=R[tostate,fromstate]
		      //noff = R[e].PoissonOffspring(eigenratio);
		      noff = R[e].PoissonOffspring(1);
		      
		      if (RandLibObj.CheckDiscreteLookup()==0)
			{
			  noff=0;
			}
		      
		      if (noff>0)
			{
			  I[k].SetCurrentLastRep(t);
			  I[k].SetCurrentNumOff(noff);
			  
			  /*
			    choosing mate.  At this point the effects of genotype upon the mates
			    ability to produce pollen could be inserted.
			  */
			  if (!multiple_paternity)///all offspring from one father
			    {
			      if (RandLibObj.uniform()<self)
				{
				  mate = searchI;
				}
			      else
				{
				    mate = FindMate(searchI);
				}
			    }
			  for (q=0;q<noff;q++)
			    {
			      if (multiple_paternity)///each offspring the product of mixed mating
				{
				  if (RandLibObj.uniform()<self)
				    {
				      mate = searchI;
				    }
				  else
				    {
				      mate = FindMate(searchI);
				    }
				}
			      
			      ///"do the deed" between searchI and mate. tmpI is the baby
			      if (mate.GetClass()>-1)//if no mate, no offspring
				{
				  tmpI = searchI.repro_sex(searchI,mate,t,Atbls);
				  
				  tmpI.SetClass(j);
				  
				  ///this could/should be made user selectable
				  tmpI.SetSex(0);///would require some more code modifications, but might be worth it
				  
				  tmpI.SetGen(t);
				  
				  tmpI.SetMID(0);
				  tmpI.SetPID(0);
				  tmpI.SetID(0);
				  
				  tmpI.SetMID(searchI.GetID());
				  tmpI.SetPID(mate.GetID());
				  tmpI.SetID(nextID);
				  if (nextID>MAXIDS) 
				    {
				      nextID=1;
				    }
				  else
				    {
				      nextID=nextID+1;
				    }
				  
				  tmpI.Change(-1);
				  //cerr<<"running birth"<<endl;
				  
				  tmpI.Birth(t,Atbls);
				  
				  //cerr<<"birth run"<<endl;
				  //				  err = 0;
				  if (I[j].AddIndividual(tmpI)<0)
				    {
#ifdef DEBUG
				      cerr << "adding an individual failed" << endl;
#endif
				    }
				  //cerr<<"added an individual.  ID"<<nextID<<endl;
				}
			    } //q
			  
			}//end if noff>0
		    }//j
 RandLibObj.FreeDiscreteLookup();
 I[k].NextIndividual();
		}//l
	}
    } //if R[e].AnyFrom
      
}//k 

}//end of function Reproduce


void Landscape::Extirpate()
{
  int h,etrue;
  size_t cl,sz;
  double rn;

  std::vector<int> p;

  p.resize(nhab);

  sz = nhab * s;

  etrue=0;
  for (h=0;h<nhab;h++)
    {
      rn = RandLibObj.uniform();
      if (evec[e][h]>rn)
	{
	  p[h]=1;
	  etrue=1;
	}
    }
  if (etrue)
    {
      for (cl=0;cl<sz;cl++)
	{
	  if (p[Habitat(cl)])
	    {
	      I[cl].ClearClass(t,Atbls);
	    }
	}
    }
}

void Landscape::CompressInd()
{
  int j,sz;
  sz=I.size();
  for (j=0; j<sz; j++)
    {
      I[j].CompressClass(1);
    }
}

void Landscape::CarryState(size_t maxsz, int i)
{
  int numdel,k;
  //KKM 8.1.05..................................................................
  size_t max;
  
  //setting the hard ceiling on abundance at 110% of carrying capacity when density
  //dependence is active.  That way, long-term average abundance is closer to K.
  if (densdepdemo == 1) max = int(floor(0.5 + maxsz*1.1));
  else max = maxsz;
  //............................................................................

  if (max<I[i].size())
    {
      numdel = (I[i].size()-max);
      for (k=0;k<numdel;k++)
	{
	  I[i].RemoveRandomInd(t,Atbls);
	}
    }
}

void Landscape::HabCarry(int k)
{
  int h;
  size_t j, sz;

  std::vector <double> prop;

  sz = nhab * s;
  prop.resize(nhab);

  for (h=0;h<nhab;h++)
    {
      if (k<0)
	{
	  prop[h] = double(kvec[e][h])/double(PopSize(h));
	}
      else
	{
	  prop[h] = double(k)/double(PopSize(h));
	}
      if (prop[h]>1) {prop[h]=1.0;}
    }
  
  for (j=0;j<sz;j++)
    {
      h = Habitat(j);
      CarryState(size_t((prop[h])*I[j].size()),j);
    }
}

void Landscape::LambdaAdjust(int bypop)
{
  int i, j, k, l, bigto, bigfrom;
  double pred_l, sim_l, adjrate;
  TransMat diag, Spopmat, Rpopmat;
  if (bypop)
    {
      if (bypop==1) 
	{
	  diag.SetSize(s);
	  Spopmat.SetSize(s);
	  Rpopmat.SetSize(s);
	  for (i=0;i<nhab;i++)
	    {
	      for (j=0;j<s;j++)
		for (k=0;k<s;k++)
		{
		  bigto = (i*s)+k;
		  bigfrom = (i*s)+j;
		  Spopmat.SetElement(k,j,S[e].GetElement(bigto,bigfrom));
		  Rpopmat.SetElement(k,j,S[e].GetElement(bigto,bigfrom));
		}
	      pred_l = (Spopmat+Rpopmat).Lambda();
	      sim_l = (Spopmat*(Rpopmat+diag)).Lambda();
	      adjrate = pred_l/sim_l;
	      for (l=(i*s);l<((i*s)+s);l++)
		{
		  CarryState((int)(I[l].size()*adjrate+0.5),l);
		}
	    }
	}
      else
	{
	  diag.SetSize(s*nhab);
	  pred_l = (S[e]+R[e]).Lambda();
	  sim_l = (S[e]*(R[e]+diag)).Lambda();
	  adjrate = pred_l/sim_l;
	  
	  for (i=0;i<(s*nhab);i++)
	    {
	      CarryState((int)(I[l].size()*adjrate+0.5),l);
	    }
	}
    } //if bypop is ne 0;
}
  
void Landscape::LandCarry()
{
  size_t j, sz;
  double pr;

  sz = nhab * s;
  pr =  double(maxlandsz)/double(PopSize(-1));
  for (j=0;j<sz;j++)
    {
      CarryState(size_t(pr*I[j].size()),j);
    }
}

void Landscape::GCAlleles()
{
  int i;
  for (i=0;i<getloci();i++)
    {
      Atbls[i]->CalcProps();
      Atbls[i]->GCAlleles();
    }  
}

ostream &Landscape::WriteLoci(ostream &stream)
  {
    stream << Atbls <<endl;
    return stream;
  }



ostream &operator<<(ostream &stream, Landscape &l)

{
  int ie, id, it, i;
  size_t j, sz;
  
  sz = l.nhab * l.s;
  
  stream << "nhab      " <<" "<< l.nhab << endl;
  stream << "stages    " <<" "<< l.s << endl;
  stream << "ndemo     " << " " << l.ndemo <<endl;
  stream << "rdemo     " << " " << l.rdemo <<endl;
  stream << "nepochs   " <<" "<< l.nep<< endl;
  stream << "cepoch    " <<" "<< l.e<< endl;
  stream << "repoch    " <<" "<< l.randepoch<< endl;
  stream << "ngen      " <<" "<< l.ngen<< endl;
  stream << "cgen      " <<" "<< l.t<< endl;
  stream << "self      " <<" "<< l.self<< endl;
  stream << "nextID    " <<" "<< l.nextID<< endl;
  stream << "multp     " <<" "<< l.multiple_paternity<< endl;
  stream << "maxlandsz " <<" "<< l.maxlandsz<< endl;
  stream << "densdep   " <<" "<< l.densdepdemo<<endl;

  stream << endl;

  stream << "epochvec  " <<" "<< endl;
  for (ie=0;ie<l.nep;ie++)
    {
      stream << l.epochs[ie] << " " ;
    }

  stream << endl;

  stream << "epochprobs" <<" "<< endl;
  for (ie=0;ie<l.nep;ie++)
    {
      stream << l.epochprobs[ie] << " " ;
    }
  stream << endl;

  stream<< "matrices  " << " " << endl;
  for (id=0;id<l.ndemo;id++)
    {
      stream << l.LM[id] << endl;
    }
  for (ie=0;ie<l.nep;ie++)
    {
      stream << l.S[ie] << l.R[ie] << l.M[ie] << endl ;

      for (id=0;id<l.ndemo;id++)
	{
	  stream << l.demoProbVec[ie][id] << " " ;
	}
      stream << endl;
      for (it=0;it<l.nhab ;it++)
	{
	  stream << l.evec[ie][it] << " " ;
	}
      stream << endl ;
      for (it=0;it<l.nhab ;it++)
	{
	  stream << l.kvec[ie][it] << " " ;
	}
      stream << endl << endl;

    }

  stream << "loci        "<<l.getloci()<<endl;

  for (i=0;i<l.getloci();i++)
    {
      l.Atbls[i]->CalcProps();
      l.Atbls[i]->GCAlleles();

      if (l.Atbls[i]->getClassType()==INFALLELETBL)
	{
	  stream << 0 << endl;
	}
      else if (l.Atbls[i]->getClassType()==STEPALLELETBL)
	{
	  stream << 1 << endl;
	}
      else if (l.Atbls[i]->getClassType()==SEQALLELETBL)
	{
	  stream << 2 << endl;
	}
      else
	{
#ifdef DEBUG
	  cerr << "Don't understand the ID of this locus " <<endl;
#endif
	  assert(0==1);
	}
      l.Atbls[i]->Write(stream);
      stream << endl;
    }
  
  stream <<     "individ   " << " " << l.PopSize(-1) << endl;
  for (j=0;j<sz;j++)
    {
	  stream << l.I[j];
    }
  return stream;
}


istream &operator>>(istream & stream, Landscape &l)
{
  //  size_t sz;
  int i, j, ni, h, tmpi, nl ;
  int indflag=0;
  int ktot,maxk;
  std::vector <int> p;
  double totprob;


  PackedIndividual ind;
  
  Allele tmpallele;
  
  string c,tmps,alist,nlist,wlist ;
  int epochprobs = 0;
  char tc;
  char cp[256];

  alist = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
  nlist = "e-+0123456789.";
  wlist = " \n\t\f";

  nl=1;



  l.init();

  while (!stream.eof())
    {
      tc = stream.peek();
      if (short(tc)>=0)
	{
	  if (wlist.find(tc)==string::npos && nl) //no white space
	    {
	      nl=0;
	      if (!(tc=='#')) {
		      
		c.erase(c.begin(),c.end());
		for (i=0;i<TOKENLEN;i++)
		  {
		    stream.get(tc);
		    if (tc!=' ')
		      {
			tmps=tc;
			c.append(tmps);
		      }
		    if (stream.eof()) {break;}
		  }
		      
		if      (c=="nhab")
		  {
		    stream >> l.nhab;
		  }
		else  if (c=="stages")
		  {
		    stream >> l.s;
		  }
		else if      (c=="repoch")
		  {
		    stream >> l.randepoch;
		    epochprobs = 1;
		  }
		else  if (c=="nepochs")
		  {
		    stream >> l.nep;
		    l.setepochs(l.nep);
			  
		  }
		else  if (c=="cepoch")
		  {
		    stream >> l.e;
		  }
		else if      (c=="ndemo")
		  {
		    stream >> l.ndemo;
		    if (l.ndemo==0)
		      {
#ifdef DEBUG
			cerr << "chose not to use the habitat-level demography specification"<<endl;
#endif
		      }
#ifdef DEBUG
		    cerr << "Running l.setndemo in the landcape inserter "<<l.s <<endl;
#endif
		    l.setndemo(l.ndemo);
		  }
		else if      (c=="rdemo")
		  {
		    stream >> l.rdemo;
		    if (l.rdemo)
		      {
#ifdef DEBUG
			cerr << "chose to randomly assign habitat-level demography specification"<<endl;
#endif
		      }
		  }
		else if (c=="ngen")
		  {
		    stream >> l.ngen;
		  }
		else if (c=="cgen")
		  {
		    stream >> l.t;
		  }
		else if (c=="self")
		  {
		    stream >> l.self;
		  }
		else if (c=="multp")
		  {
		    stream >> l.multiple_paternity;
		  }
		else if (c=="nextID")
		  {
		    stream >> l.nextID;
		  }
		else if (c=="maxlandsz")
		  {
		    stream >> l.maxlandsz;
		  }
		else if (c=="densdep")
		  {
		    stream >> l.densdepdemo;
		  }
		else if (c=="epochvec")
		  {
		    for (i=0;i<l.nep;i++)
		      {
			stream >> l.epochs[i] ;
		      }
		  }
		else if (c=="epochprobs")
		  {
		    totprob=0;
		    for (i=0;i<l.nep;i++)
		      {
			stream >> l.epochprobs[i] ;
			totprob += l.epochprobs[i];
		      }
		    epochprobs=0;
		    if (totprob!=1)
		      {
#ifdef DEBUG
			cerr << "epoch probabilities do not sum to one" <<endl;;
#endif
			assert(totprob==1);
		      }
		  }
		else if (c=="matrices")
		  {
		    maxk = 0;
			  
		    for (i=0;i<l.ndemo;i++)
		      {
			stream >> l.LM[i];
		      }
			  
		    for (i=0;i<l.nep;i++)
		      {
			ktot = 0;

			stream >> l.S[i] ;
			stream >> l.R[i] ;
			stream >> l.M[i] ;

			for (j=0;j<l.ndemo;j++)
			  {
			    stream >> l.demoProbVec[i][j];
			  }
			for (j=0;j<l.nhab;j++)
			  {
			    stream >> l.evec[i][j] ;
			  }
			for (j=0;j<l.nhab;j++)
			  {
			    stream >> l.kvec[i][j] ;
			    ktot = l.kvec[i][j] + ktot;
			  }
			if (ktot>maxk) 
			  {
			    maxk=ktot;
			  }
		      }
		  }
		      
		else if (c=="loci")
		  {
		    stream >> l.Atbls;
		    l.setloci();
		  }
		      
		else if (c=="individ")
		  {
		    if (!indflag)
		      {
			stream >> ni ;
			assert(ni > 0);
			l.I.resize( l.nhab * l.s );
			for (i=0;i<ni;i++)
			  {
			    l.SetUpInd(ind);
			    stream >> ind;
			    ind.Birth(-1,l.Atbls);
			    l.I[ind.GetClass()].AddIndividual(ind);
			  }
			indflag=1;
		      }
		    else
		      {
#ifdef DEBUG
			cerr << "already defined popsize vectors, can't define inds";
#endif
			assert(1==0);
		      }
		  }
		else if (c=="popinit")
		  {
		    if (!indflag)
		      {
			stream >> h ;
			for (i=0;i<h;i++)
			  {
			    stream >> tmpi;
			    p.push_back(tmpi);
			  }
			l.popsizeset(p);
			indflag=1;
			l.GCAlleles();
		      }
		    else
		      {
#ifdef DEBUG
			cerr << "already defined individual vectors";
#endif
			assert(0);
		      }
		  }
		else
		  {
#ifdef DEBUG
		    cerr << "unrecognized token `" << c << "' in Landscape inserter"<<endl ;
#endif
		  }
	      }
	      else
		{
		  stream.getline(cp,256,'\n');
		  nl = 1;
		}
	    }
	  else 
	    {
	      stream.get(tc);
	      if (tc=='\n')
		{
		  nl=1;
		}
	    }
	}
      else
	{
	  break;
	}
    }
  if (epochprobs)
    {
#ifdef DEBUG
      cerr << "asked for random epochs, but did not include probabilities" <<endl; ;
#endif
      assert(1==2);
    }
  l.GCAlleles();
  return stream;
}



void Landscape::Advance()
{
  t++; //Increment generation
  ChooseEpoch();
  ConstructDemoMatrix();
}


///begin implementation of Landscape_statistics


Landscape_statistics::Landscape_statistics (int /*h*/, int /*stg*/, int /*loc*/, int /*ep*/, int /*gn*/)
{
  ///  init(h,stg, loc, ep, gn);
}

Landscape_statistics::~Landscape_statistics ()
{
	  
#ifdef DEBUG
  cerr << "Landscape_statistics destructor starting" <<endl;
#endif

#ifdef DEBUG
  cerr << "Landscape_statistics destructor exiting" <<endl;
#endif

}
/*****
void Landscape_statistics::Statistics(ostream & streamout)
{

  int i,sgz  ;
  sgz = s * nhab;

  streamout << t <<"  " <<e <<" " ;
  for (i=0;i<sgz;i++)
    {
      //      streamout.form("%6i ",I[i].size());
      streamout << setw(6) << I[i].size();
    }
  streamout <<"  "<< PopSize(-1) <<"  "<<endl;;
}
*****/

double Landscape_statistics::GenLength()
{
  size_t i,sz;
  double gensiz=0.0,tmp;
  double totparents=0.0;
  sz=s*nhab;
  for (i=0;i<sz;i++)
    {
      if (R[e].AnyFrom(i))
	{
	  gensiz = gensiz + I[i].GenLength(t)*I[i].size();
	  totparents= totparents + I[i].size();
	}
    }
  tmp = gensiz/totparents;
  return tmp;
}


vector <int>  Landscape_statistics::Rmat(int numind)
{
  //numind=0 (default), sample all ind

  int i,j, k, ss, ps, p;
  size_t sz;
  vector <PackedIndividual> IVec;
  vector <int> diptbl;
  vector <int> retval;


  //find the number of occupied habitats
  IVec.reserve(numind);
  diptbl.reserve(nloc);
  retval.reserve(2000);

  for (k=0;k<nloc;k++)
    {
      diptbl.push_back(k);
    }

  ps=0;
  j=0;

  for (i=0;i<nhab;i++)
    {
      if (PopSize(i)>0)
	{
	  j++;
	}
    }

  //  streamout << "pop class individual locus aindex allele"<<endl;

  for (i=0;i<nhab;i++)
    {
      ps = PopSize(i);
      if (ps>0)
	{
	  IVec.resize(0);
	  for (j=i*s;j<((i*s)+s);j++)
	    {
	      I[j].ResetIndividuals();
	      for (sz=0;sz<I[j].size();sz++)
		{
		  IVec.push_back(I[j].GetCurrentIndividual());
		  I[j].NextIndividual();
		}
	    }
	  random_shuffle(IVec.begin(),IVec.end(),randWrapper);


	  if ((numind>0)&&(ps>numind))
	    {
	      ss = numind;
	    }
	  else
	    {
	      ss = ps;
	    }
	  
	  for (j=0;j<ss;j++)
	    {

	      for (sz=0;sz<diptbl.size();sz++)
		{
		  for (p=0;p<Atbls[sz]->getPloidy();p++)
		    {
		      retval.push_back(i); //population
		      retval.push_back(IVec[j].GetClass());//class
		      retval.push_back(j); //individual
		      retval.push_back(sz);//locus
		      retval.push_back(p);//allele index
		      retval.push_back(IVec[j].GetAllele(diptbl[sz],p)) ;//allele
		    }
		}
	    }
	}
    }
  return retval;
}




/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; minor-mode: font-lock ***
;;; End: ***
*/
