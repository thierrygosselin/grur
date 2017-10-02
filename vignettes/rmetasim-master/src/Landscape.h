/*

$Modified: astrand $

Copyright (C) 1999 Allan E. Strand

This file is part of Metasim
*/

#ifndef LANDSCAPE_H
#define LANDSCAPE_H



/*includes
*/

#include <metasim.h>
#include <utilities.h>
#include <PackedIndividual.h>
#include <Democlass.h>
#include <TransMat.h>
#include <RandLib.h>
#include <string>
#include <vector>
#include <list>
#include <iostream> 
//#include <fstream> 


using namespace std;




///Local Matrix information Class
/**

   This class implements a set of matrices that model demography within populations.
   There are three matrices all of type TransMat.

   1) a survival matrix
   2) a reproduction matrix
   3) a male gamete weight matrix (usually if a category can contribute gametes, it is 
      given a 1, but not necessarily) 

 */
class LocalMat {

  ///size of each of the matrices
  size_t sz;
  ///list of survival matrices local to a particular habitat.  
  TransMat Slocal;
  ///list of repro matrices local to a particular habitat
  TransMat Rlocal;
  ///list of male function matrices local to a particular habitat
  TransMat Mlocal;

public:

  LocalMat () {}
  ~LocalMat() {}
  
inline  double GetSlocalVal(size_t from, size_t to)
    {
      assert(from<Slocal.Size());
      assert(to<Slocal.Size());
      Slocal.SetFromState(from);
      Slocal.SetToState(to);
      return Slocal.Value();
    }

inline  double GetRlocalVal(size_t from, size_t to)
    {
      assert(from<Rlocal.Size());
      assert(to<Rlocal.Size());
      Rlocal.SetFromState(from);
      Rlocal.SetToState(to);
      return Rlocal.Value();
    }

inline  double GetMlocalVal(size_t from, size_t to)
    {
      assert(from<Mlocal.Size());
      assert(to<Mlocal.Size());
      Mlocal.SetFromState(from);
      Mlocal.SetToState(to);
      return Mlocal.Value();
    }



inline  void SetSlocalVal(size_t from, size_t to, double val)
    {
      assert(from<Slocal.Size());
      assert(to<Slocal.Size());
      assert(val>=0);
      Slocal.SetElement(from,to,val);
    }

inline  void SetRlocalVal(size_t from, size_t to, double val)
    {
      assert(from<Rlocal.Size());
      assert(to<Rlocal.Size());
      assert(val>=0);
      Rlocal.SetElement(from,to,val);
    }

inline  void SetMlocalVal(size_t from, size_t to, double val)
    {
      assert(from<Mlocal.Size());
      assert(to<Mlocal.Size());
      assert(val>=0);
      Mlocal.SetElement(from,to,val);
    }


  void  SetSize(size_t sz);

  friend ostream &operator<<(ostream & stream, LocalMat &lm);

  friend istream &operator>>(istream & stream, LocalMat &lm);



}; //end localmat



///Landscape Class
  /**
This is the declaration of the Landscape class.  The landscape class
implements the simulation, basically.

o It is built upon the concept of a 2d array of suitable sites.

o The landscape moves through time one year (generation) at a time.  

o every x generations a new epoch can take over, where the migration
  matrices, demography, everything can change.

o the migration matrix provides a value, m_ij which is the
  probablility that an individual will migrate from population i to
  population j for each pair of populations

*/


class Landscape {
protected:
  /// Individuals data structure.
  /**

     This data structure (I) holds all of the individuals in the
     landscape.

     Right now, it's implemented as a vector of lists.  Each of the
     lists refers to the individuals within a particular habitat.
     Each generation individuals will be added and deleted from the
     lists in each habitat, but the number of habitats will not
     change.

     These habitats are determined solely by the state of the
     individual.  So, a landscape with 3 habitats (p) and 3
     demographic stages (s) is a 9x9 matrix.  All individuals with
     states in columns 1-3 are in habitat 1, 4-6 in habitat 2, and 7-9
     in habitat 3.  Individuals don't know how to convert state into
     habitat, so this has to be a Landscape-level function.

  */
  /*
    
    Variable declarations

   */

  ///Title of landscape
  string title;

  /// Individuals
  vector< DemoClass > I;

  ///Lookup table of alleles at all loci
  AlleleLookTbl Atbls;
  /// number of potentially suitable habitats in the landscape.
  int nhab;

  /// number of demographic stages
  int s   ;

  /// number of genetic loci
  int nloc;

  /// number of populations in X-dimension
  int xdim;

  /// number of populations in Y-dimension
  int ydim;
  
  /// selfing rate
  double self;

  /// the number of different types of migration and other demographic conditions (number of epochs)
  int nep ;

  ///number of different within-habitat demographies to choose from.
  ///This variable determines the length of the vector LM defined in
  ///this class.  It also defines the length of each element in the
  ///outer vector of the variable "demoProbVec"

  int ndemo;

  /// a flag to tell whether the demographies are chosen at random from demoProbVec
  int rdemo;

  ///KKM 5.20.05..............................................................
  /// a flag to tell whether demography is density dependent
  int densdepdemo;
  ///........................................................................

  /// a counter of the current epoch
  int e   ;

  /// when true, choose epochs at random at every advance of the clock.  Ignores epoch vector.

  int randepoch;

  /// total number of generations
  int ngen;

  /// the current generation number 
  int t   ;

  ///the percentage of deviation allowed when imposing carrying capacity
  double habdelta;

  ///counter for the next individual.  Usually the total number of births in a landscape, though can roll over at MAXIDS
  int nextID;

  ///The maximum number of individuals allowed in the landscape
  int maxlandsz;
  ///is multiple paternity allowed? 0=single father, 1=every child has a randomly selected father 
  int multiple_paternity ;
  /**

    the next several lines define vectors of "vital matrices" each of
    these vectors is nep long.  Therefore each"epoch" can have a
    drastically different migration model, number of sites allowed to
    be occupied, demography, extinction, etc..

   */

  /// a vector of epoch begin dates
  std::vector<int> epochs;
  /// a vector of probabilities of seeing the conditions of a particular epoch in any given year
  std::vector<double> epochprobs;

  ///Survival matrix
  std::vector<TransMat> S;
  ///Reproduction matrix 
  std::vector<TransMat> R;
  ///Probability of an individual contributing a gamete to another individual in the same or different class
  std::vector<TransMat> M;

  ///A vector of Local Matrix Types.  These represent the demography within populations and can be used to
  ///set the diagonal elements of S, R, and M.
  ///The length of this vector is the same as the number of different within-pop demographies possible

  std::vector<LocalMat> LM;

  ///KKM 5.20.05...............................................................
  ///A vector of Local Matrix Types.  These represent the demography within populations
  ///near carrying capacity when density dependence is used and can be used to
  ///along with the LM matrices to interpolate the diagonal elements of S, R, and M.
  ///The length of this vector is the same as the number of different within-pop demographies possible

  std::vector<LocalMat> LMK;
  ///..........................................................................

  ///This structure is a vector of length nep.  Each element is is a
  ///vector of probabilities of observing a particular local
  ///demography of type LocalMat. The length of each of these
  ///sub-vectors is the number of different possible local demographies

  std::vector<std::vector <double> > demoProbVec ;

  ///"non-demographic" annual extinction probabilities for each site
  std::vector<std::vector<double> > evec;
  /// carrying capacities for each site
  std::vector<std::vector<int> > kvec;

  /**
     Table of allele frequencies rows are alleles and cols are habitats

   */  

public:

  ///Constructor
  Landscape (int h=1, int stg=2, int loc=1, int ep=1, int nd=1, int gn=2);
  ///Destructor
  virtual ~Landscape();

  /*
    set the initial parameters
  */

  inline void setself(double slf=0) {self=slf;}
  inline void setmultp(int mp=1) {multiple_paternity=mp;}
  inline void setranddemo(int rd=1) {rdemo=rd;}
  ///KKM 5.20.05...........................................................
  inline void setdensdep(int dd=0) {densdepdemo=dd;}
  ///......................................................................
  inline void setgens(int gn=2) {ngen=gn;}
  inline void setCgen(int cg) {t=cg;}
  inline void setCepoch(int ce) {e=ce;}

  inline void setMaxLandSize(int mls=300000) {maxlandsz=mls;}
  inline void assignRandEpoch(int re=1) {randepoch=re;}

  inline void setRandEpoch() {randepoch=1;}
  inline void unsetRandEpoch() {randepoch=0;}
  inline void reserveclasses() {I.resize(s*nhab);}

         void setepochs(int ep=1);
         void setndemo(int nd=1);
         void sethabs(int h=1);
         void setstages(int stg=2);
  void setxdim(int x=0);
  void setydim(int y=0);
         
  inline void setnextID(int newid=1) {nextID=newid;}
  inline void setloci() {nloc=Atbls.size();}
         void setepochprob(int ce, double prob);
         void setepochstart(int ce, int strt);

  /*
    get values of the parameters
  */


  inline int gethabs() {return nhab;}
  inline int getstages() {return s;}
    inline int getxdim() {return xdim;}
    inline int getydim() {return ydim;}
  inline double getself() {return self;}
  inline int getmultp() {return multiple_paternity;}
  inline int getloci() {setloci(); return nloc;}
  inline int getepochs() {return nep;}
  inline int getCgen() {return t;}
  inline int getCepoch() {return e;}
  inline double getepochprob(int ce) {return epochprobs[ce];}
  inline int getepochstart(int ce) {return epochs[ce];}
  inline int getgens() {return ngen;}
  inline int getrandepoch() {return randepoch;}
  inline int getranddemo() {return rdemo;}
  ///KKM 5.20.05.............................................................
  inline int getdensdep() {return densdepdemo;}
  ///........................................................................
  inline int getndemo() {return ndemo;}
  inline int getnextID() {return nextID;}
  inline int getMaxLandSize() {return maxlandsz;}

  void init(int h=1, int stg=2, int loc=0, int ep=1, int nd=1, int gn=1);

  /**
    set the appropriate matrices
   */
  void setS(TransMat a, int ep=-1);
  void setR(TransMat a, int ep=-1);
  void setM(TransMat a, int ep=-1);

    inline void setSmatElement(int /*ep*/, int t, int f, double val)  {S[e].SetElement(f,t,val);}
    inline void setRmatElement(int /*ep*/, size_t t, size_t f, double val)  {R[e].SetElement(f,t,val);}
    inline void setMmatElement(int /*ep*/, size_t t, size_t f, double val)  {M[e].SetElement(f,t,val);}

    inline double getSmatElement(int /*ep*/, size_t t, size_t f)  {return S[e].GetElement(f,t);}
    inline double getRmatElement(int /*ep*/, size_t t, size_t f)  {return R[e].GetElement(f,t);}
    inline double getMmatElement(int /*ep*/, size_t t, size_t f)  {return M[e].GetElement(f,t);}


  inline void setLSmatElement(int d, size_t t, size_t f, double val) {LM[d].SetSlocalVal(f,t,val);}
  inline void setLRmatElement(int d, size_t t, size_t f, double val) {LM[d].SetRlocalVal(f,t,val);}
  inline void setLMmatElement(int d, size_t t, size_t f, double val) {LM[d].SetMlocalVal(f,t,val);}

  inline double getLSmatElement(int d, size_t t, size_t f)  {return LM[d].GetSlocalVal(f,t);}
  inline double getLRmatElement(int d, size_t t, size_t f)  {return LM[d].GetRlocalVal(f,t);}
  inline double getLMmatElement(int d, size_t t, size_t f)  {return LM[d].GetMlocalVal(f,t);}

  ///KKM 5.20.05..............................................................
  inline void setLSKmatElement(int d, size_t t, size_t f, double val) {LMK[d].SetSlocalVal(f,t,val);}
  inline void setLRKmatElement(int d, size_t t, size_t f, double val) {LMK[d].SetRlocalVal(f,t,val);}
  inline void setLMKmatElement(int d, size_t t, size_t f, double val) {LMK[d].SetMlocalVal(f,t,val);}

  inline double getLSKmatElement(int d, size_t t, size_t f)  {return LMK[d].GetSlocalVal(f,t);}
  inline double getLRKmatElement(int d, size_t t, size_t f)  {return LMK[d].GetRlocalVal(f,t);}
  inline double getLMKmatElement(int d, size_t t, size_t f)  {return LMK[d].GetMlocalVal(f,t);}
  ///KKM.......................................................................

inline void ConstructDemoMatrix()
  {
    if (ndemo){
	   if (rdemo){
         ///KKM 6.7.05........................................................
         if (densdepdemo){
           RandomDensityDependentDemoMatrix();
         }
         else {
         ///..................................................................
            RandomlyConstructDemoMatrix();
         }
      }
	   else{
         ///KKM 6.7.05........................................................
         if (densdepdemo){
           SequentialDensityDependentDemoMatrix();
         }
         else {
         ///..................................................................
	        SequentiallyConstructDemoMatrix();
         }
      }
    }
  }
  /**

     This function selects sub-matrices from the vector of local
  matrices.  The sub-matrices are chosen in order, and are reused
  until all habitats have a matrix assigned to them.  The method takes
  these matrices and inserts them on the diagonal of the S, R, and M
  matrices.  This function DOES NOT ALTER THE OFF-DIAGONAL
  SUBMATRICES.  These are defined in the current epochs S,R, and M
  matrices.

     The current epoch determines which  S,R,M, and demoProbVec to choose.

   */
void SequentiallyConstructDemoMatrix();

///KKM....................................................................
  /**

  Works the same as SequentiallyConstructDemoMatrix(), except values
  entered into the sub-matrices are interpolated from the values at
  zero population density and carrying capacity.

   */
void SequentialDensityDependentDemoMatrix();
///...........................................................................

  /**

     This function randomly selects sub-matrices from the vector of
     local matrices.  The sub-matrices are chosen from a multinomial
     distribution given by "demoProbVec".  The method takes these
     matrices and inserts them on the diagonal of the S, R, and M
     matrices.  This function DOES NOT ALTER THE OFF-DIAGONAL
     SUBMATRICES.  These are defined in the current epochs S,R, and M matrices.

     The current epoch determines which  S,R,M, and demoProbVec to choose.

   */
void RandomlyConstructDemoMatrix();

/*/KKM 6.7.05.................................................................

  Works the same as SequentiallyConstructDemoMatrix(), except values
  entered into the sub-matrices are interpolated from the values at
  zero population density and carrying capacity.

   */
void RandomDensityDependentDemoMatrix();
///............................................................................
  /*
    set the population characteristic vectors
   */
  ///Set up the extinction vector(s).  'ev' is a vector of type double that correpond
  ///to extinction rates for each of the 0..nhab habitats. 'ep' is the epoch
void setextinct(int ep, double *ev);
  ///Get an extinction vector(s).  for each of the 0..nhab habitats. 'ep' is the epoch
void getextinct(int ep, double *ev);
  ///Set up the carrying capacity vector(s) 'cv' is a vector of type
  //int that correspond to extinction rates for each of the 0..nhab
  ///habitats. 'ep' is the epoch
void setk(int ep, int *cv);
  ///Get carry capacity vector(s).  0..nhab habitats. 'ep' is the epoch
void getk(int ep, int *cv);

  ///Set up the local demography vector for a particular epoch
void setldemovector(int ep, double *dv);
  ///Get the local demography vector for a particular epoch
void getldemovector(int ep, double *dv);

void zeroextinct();
void zerok();


  ///This resets the pointer to members of the individual class 'cl' to the start of the list
  inline void resetStage(int cl)
  {
    I[cl].ResetIndividuals(); 
  }
  ///This function advances the pointer to inds in the demographic class 'cl'.
  ///it returns 1 if there is another ind to grab, otherwise it returns 0

  inline int advanceStagePtr(int cl)
  {
    return I[cl].NextIndividual();
  }
  ///This function tells the size of the demographic class pointed to by 'cl'
  inline int StageSize(int cl=0)
  {
    return I[cl].size();
  }

  inline PackedIndividual getNextInd(int cl)
  {
    return I[cl].GetCurrentIndividual();
  }

  inline void SetUpInd(PackedIndividual &ind)
  {
    ind.resetLoci(Atbls);
  }
  inline int addIndividual(PackedIndividual ind, int t)
  {
    ind.Birth(t,Atbls);
    return (I[ind.GetClass()].AddIndividual(ind));
  }

  ///choose an epoch using some selection criteron
void ChooseEpoch();
  ///randomly chooses an epoch, random epoch selection is in action
void RandomlyChooseEpoch();
  ///Initialize the Populations
void popsizeset(std::vector<int> & ps);

  /*
    functions that allow outside processes to interact with the allele lookup tbl
   */

  inline void Atbl_push_back(AlleleTbl * atp)
  {
    Atbls.push_back(atp);
  }

  inline int LocusGetClassType(int l)
  {
    return Atbls[l]->getClassType();
  }

  inline int LocusGetPloidy(int l)
  {
    return Atbls[l]->getPloidy();
  }

  inline int LocusGetTrans(int l)
  {
    return Atbls[l]->getTrans();
  }

  inline double LocusGetMutRate(int l)
  {
    return Atbls[l]->getMutationRate();
  }

  inline vector<int> LocusGetAindices(int l)
  {
    return Atbls[l]->getAindices();
  }

  inline void LocusGetAlleleRef(int l, int andx, Allele* ptr)
  {
    if (Atbls[l]->getClassType()==SEQALLELETBL)
      {
	Atbls[l]->getAlleleRef(andx,(dynamic_cast<SeqAllele *>(ptr)));
      }
    else if (Atbls[l]->getClassType()==INFALLELETBL)
      {
	Atbls[l]->getAlleleRef(andx,ptr);
      }
    else if (Atbls[l]->getClassType()==STEPALLELETBL)
      {
	Atbls[l]->getAlleleRef(andx,ptr);
      }
    else 
      {
#ifdef DEBUG
	cerr <<"dont know what type of locus this is"<<endl;
#endif
	assert(1==0);
      }
  }

  /**

     This function takes a deomgraphic stage on a s*p x s*p matrix and returns the habitat that that stage belongs to


   */
int Habitat(int stage);


  /** 
      This function reports the population size of the habitat given
      by the integer argument.  If the argument is -1, then the total
      size of all habitats are reported
 */

int PopSize(int p=-1);




  /**

     Survive

     Survive might be better named because this is the method that
     implements survival, growth, and migration.
     
   */



void Survive();

/** 

    This method sets the discrete lookup table in the global
    RandLibObj.  The table is set with probabilities that individual
    demographic classes will contribute a mate to a particular pi
    passed in.  

    ***Important: this function allocates a lookup table.  it must be
    ***freed at some point after the function is invoked by issuing the
    ***command: RandLibObj.FreeDiscreteLookup();

*/
//int CalculateMaleGameteClassVector(PackedIndividual pi);
int CalculateMaleGameteClassVector(int k);

  /**
This function goes through each of the classes and performs internal
cleanups.  May increase speed, but there will be tradeoffs
   */
  void CompressInd();

  /**

     FindMate:

     Find a mate from the landscape for the passed individual and return it.

     The mate is identified by first multiplying the population size
     in each class by the probability that an individual in that class
     will become part of the pool of individuals that could
     potentially mate with the focal individual A mate is then chosen
     at random from the collection of potential gamete donators.

   */

PackedIndividual FindMate(PackedIndividual pi);

void testfindmate(PackedIndividual pi);

  /**
     Reproduce:

     o traipses through all of the individuals 
     o Decides on the number of offspring produced by every individual.
     o Matches individuals to their mates (if necessary)
     
   */

void Reproduce();


  /** 

      extirpate will go through extinction vector pick random numbers
      from a binomial distribution and set the census size at those
      sites to zero

  */
void Extirpate(); 



  /**
     
     Adjust the size of a particular state down to a particular size
     
   */
void CarryState(size_t maxsz, int i);

  /**

     set habitat sizes down to carrying capacity

   */
void HabCarry(int k = -1);

  /**

Adjust for the difference between single matrix model lambda and the
lambda in the simulation.  Adjustment occurs by killing individuals at random every generation.  bypop==1 means use the eigenvalue for each subpop to adjust.  bypop!=1 means treat entire landscape as a single transition matrix.

   */
void LambdaAdjust(int bypop = 1);
  /**

     set overall landscape size down to "maxlandsz" carrying capacity

   */
void LandCarry();


  /**

     shuffles the individuals vector randomly.  Must be used between
     generations if the order in which individuals are chose is
     important.  Neglecting it could simulate inbreeding in the case of sexual reproduction


   */
void Randomize();
  

void Advance();

  /**

     This function goes through and checks the number of copies of each allele present in the entire landscape.  If an allele is missing from the table, an error should be generated.  If there are no individuals with an allele, eliminate that allele from the table

   */
void GCAlleles();  

  ostream &WriteLoci(ostream &stream);

  friend ostream &operator<<(ostream & stream, Landscape &l);

  friend istream &operator>>(istream & stream, Landscape &l);


}; // end Landscape



class Landscape_statistics: public Landscape {

public:
Landscape_statistics (int h=1, int stg=2, int loc=1, int ep=1, int gn=2);
~Landscape_statistics ();

  /** Statistics 

      Reports various staistics about the state of the landscape.
      Sends the report to standard out.  Would be nice to somehow
      determine the nature of the output produced by Statistics with a
      run-time choice...

*/
  //void Statistics(ostream &streamout = cout);

  //Demographic Stats
  /**

     Return generation length by averaging the ages of reproductive individuals weighted by the number of offspring they produce

  */
double GenLength();

  /**
     returns the landscape as a vector of int to be used in some r popgen analyses.  Called by rmetasim
  */
  vector <int> Rmat(int numind=0);

};

#endif


/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; minor-mode: font-lock ***
;;; End: ***
*/
