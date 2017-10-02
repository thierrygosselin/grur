/*

$Modified: astrand $

Copyright (C) 1999 Allan E. Strand

This file is part of Metasim

*/

#ifndef TRANSMAT_H
#define TRANSMAT_H
/*
includes
*/

#include <metasim.h>
#include <RandLib.h>
#include <iostream>

using namespace std;
/// Generic transition matrix
/**  

   TransMat implements a generic transition matrix for simulating
   markov chains it also maintains pointers to the current from and to
   states.  The pointers are also updated by every method that takes
   them as a parameter.  

*/

class TransMat {
private:

  size_t size;
  int f;
  int t;
  std::vector< std::vector<float> > tm;

public:
  TransMat ( size_t s=1 ) ;
  ~TransMat () ;

  /// Sets a matrix cell value (returns a zero if successful)
  inline int SetElement(size_t lf, size_t lt, double val) 
    { 
      //size_t sz;
      f=lf;
      t=lt;
      //sz = tm[t].size();
      return (tm[t][f] = val) > 0; 
    }
  
  
  ///returns value of element at r,c  
  inline double GetElement(size_t lf, size_t lt) 
    { 
      f=lf;
      t=lt;
      return tm[t][f]; 
    } 

  ///Returns the size of the matrix (assumes a square matrix)
  inline size_t Size () 
    { 
	return size; 
    }
  ///Sets the size of the matrix: Resizes the matrix and sets "size"
  void SetSize(size_t sz = 1);
  ///sets the  from  state
  inline void SetFromState(size_t fs=0) 
    { 
      f = fs; 
    }

  ///sets the  to  state (obviously) 
  inline void SetToState(int ts=0) 
    {
        t = ts; 
    }
  ///gets the  from  state
  inline size_t GetFromState()
    {
      return f;
    }
  ///yadda yadda...to state
  inline size_t GetToState()
    {
      return t;
    }

  ///returns the Value at the current from and to coordinates
  inline float Value() 
    {
      return tm[t][f];
    }
  ///Sets an entire TransMat of size s from a 2d array pointed to by a
  void SetMat(TransMat a);

  ///

  /**
     Takes the current from state and makes a vector of probs that the state winds up in any of the 
to states.
   */
void SetRandomToStateVec(double eigenratio=1.0);

  /**
     Takes the current from state and makes a vector of probs that the state winds up in any of the 
to states.
   */
void SetRandomFromStateVec();

  ///Returns the state of an indiviudal in the next generation (-1 means dead)
  /** 
      
      This method treats the from column of the matrix as a
      multinomial prob dist and chooses the appropriate to value from
      the distribution

 */
  size_t RandomState();


  ///Returns the number of offspring (haploid or diploid) in the next generation
  /** 

      This method takes the value at the current from and two
      coordinates and returns the number of offspring produced from
      choosing a random variate from a poisson distribution

 */
  size_t PoissonOffspring(double eigenratio=1.0);

  /**

     returns 0 if there are no outputs from a "from" class. Else return 1

     Useful for ignoring classes that do not reproduce
     
   */
  int AnyFrom(size_t fs);

  void Diag();
  //calculates the leading eigenvalue of the TransMatrix
  double Lambda();

  ///overloaded operators
  TransMat operator+(TransMat TM);
  TransMat operator*(TransMat TM);

  ///Inserter
  friend ostream &operator<<(ostream &stream, TransMat & TM);
  ///extractor)
  friend istream &operator>>(istream &stream, TransMat & TM);

};//end TransMat


class DemoVec {
private:
  std::vector < double > v ;
public:
  DemoVec (int h=1) {v.resize(h);}
  ~DemoVec () {};

  int size() {return v.size();}
  void resize(int h=1) {v.resize(h);}

  double Val(int h=0) const
    {
      return v[h];
    }
  void Set(double d,int h=0)
    {
      v[h]=d;
    }

  friend ostream &operator<<(ostream & stream, DemoVec &d);
  friend istream &operator>>(istream & stream, DemoVec &d);
};



#endif /* TRANSMAT */

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; End: ***
*/
