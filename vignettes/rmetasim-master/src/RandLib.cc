/*

$Modified: astrand $

Copyright (C) 1999 Allan E. Strand

This file is part of Metasim
*/

/**
   
*/


/*includes
*/

#include <RandLib.h>
#include <R.h>
#include <Rmath.h>

/*
 */

RandLib::RandLib () 
{
  init();  
}


void RandLib::init()
{
    GetRNGstate();
}


RandLib::~RandLib () 
{
    PutRNGstate();
}


void RandLib::FreeDiscreteLookup()
{
}

int RandLib::CheckDiscreteLookup()
{
  return 1;
}

void RandLib::SetDiscreteLookup(double *p, int ncat) //returns a randomly chosen
                              //outcome from a multinomial distribution given in array
                              //p, of length n
     {
       double tot=0.0;
       int i;

       lp.resize(ncat);

       for (i=0;i<(ncat-1);i++)
	 {
	   tot = tot + p[i];
	   lp[i]=p[i];
	 }
       ///do different things if the passed vector totals to less than, equal to, or greater than 1
       if ((tot>=0) && (tot < 1)) //the total of the passed vector on [0,1)
	 {
	   lp[ncat-1]=1-tot;
	 }
       else if (tot==1)
	 {
	   lp[ncat-1]=0;
	 }
       else if (tot > 1)
	 {
	   if (tot > 1.5) ///the total of the vector is really big.  something is wrong.
	     {
#ifdef DEBUG
	       cerr <<"In Randlib.cc, the total of a vector passed to multinomial is much greater than 1:  "<<tot<<endl;
#endif
	       assert(tot<1);
	     }
	   else //scale all the probs to sum to one
	     {
	       for (i=0;i<(ncat);i++)
		 {
		   lp[i]= p[i]/tot;
		 }
	       lp[ncat-1]=0;
	     }
	 }
     }

int RandLib::PickMultinomial()
{
  int i;
  int *resvec = new int[lp.size()];
  double *pvec = new double[lp.size()];
  for (i=0;i<int(lp.size());i++)
    {
      pvec[i]=lp[i];
    }
  rmultinom(1,pvec,int(lp.size()),resvec);
  i=0;
  while (resvec[i]<1)
    {
      i++;
    }
  delete [] pvec;
  delete [] resvec;
  return i;
}

int RandLib::multinomial(double *p, int ncat)
{
  int rv;

  SetDiscreteLookup(p,ncat);
  rv = PickMultinomial();
  return rv;
}


///this function should return numbers from the range inclusive  of the endpoints
int RandLib::unirange(int maxval)
{
  int rv;
  double uni;
  uni = runif(0.0,maxval);
  //  cerr << "uni "<<uni<<endl;
  //  rv=int(round(runif(0.0,maxval)));
  //rv=int(fround(uni,0));
  
  //return rv;
  return (int)(uni+0.5);
}


double RandLib::uniform()
{
  return runif(0.0,1.0);
}

int RandLib::poisson(double mu)
{
  int rv;
  //  rv=int(round(rpois(mu)));
  rv=int(rpois(mu));
  return rv;
}

void RandLib::negexp_xy(double ix, double iy, double mu, double &newx, double &newy)
{
  double dir = uniform()* 2 * M_PI;
  double dist = rexp(mu);

  newx = (sin(dir)*dist) + ix;
  newy = (cos(dir)*dist) + iy;
}

double RandLib::negexp(double dist, double mu)
{
  return dexp(dist,mu,0);
}


void RandLib::SetSeed(long int /*sd*/)
{
  //  gsl_rng_set(sd);
}


/**
Declaration  of a global RandLib 'RandLibObj';
 */

RandLib RandLibObj;


/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; minor-mode: font-lock ***
;;; End: ***
*/



