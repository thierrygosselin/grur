/*

$Modified: astrand $

Copyright (C) 1999 Allan E. Strand

This file is part of Metasim

This is the declaration of the base object type.

*/

#ifndef BASEOBJ_H
#define BASEOBJ_H

/*includes
*/

#include <metasim.h>

class BaseObj {
  int classtype; // a index of the object types
public:
  BaseObj () ;
  virtual ~BaseObj ();
  void setClassType(int ct=BASEOBJ)
  {
    classtype = ct;
  }
  inline int getClassType()
  {
    return classtype;
  }
  
  
}; // end baseobj

#endif /*BASEOBJ_H*/
/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; minor-mode: font-lock ***
;;; End: ***
*/
