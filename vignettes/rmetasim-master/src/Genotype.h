/* This file is part of Metasim
   This file is teh declaration of the locus object type.
*/

#ifndef GENOTYPE_H
#define GENOTYPE_H

/* includes */
#include <FastSeqAllele.h>

typedef std::vector<int> Genotype ;


/**
   This class maintains a list of loci (Alleletbls)
 */
class GenotypePool: public BaseObj {
protected:
  std::vector<AlleleTbl *> Loci;


public:

  GenotypePool(int nloc=0);
  virtual ~GenotypePool();
  ///return a pointer to an AlleleTbl object.  The object returned is cast to the 
  ///Type given by its getClassType()
    inline AlleleTbl * getAlleleTbl(size_t /*i*/=0)
  {
    return NULL;
  }

  ///
    inline size_t push_back(AlleleTbl * /*AT*/)
  {
    return 0;
  }

  inline size_t size()
  {
    return Loci.size();
  }

  inline void resize(size_t ns=0)
  {
    Loci.resize(ns);
  }

};  // end AlleleTbl

#endif /*GENOTYPE*/

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; minor-mode:  font-lock  ***
;;; End:  ***
*/























