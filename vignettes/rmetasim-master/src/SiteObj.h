/*

$Modified: astrand $

Copyright (C) 1999 Allan E. Strand

This file is part of Metasim

This is the declaration of the site object type.


In this file the Basic Site object is defined (inherits from BaseObject)
Also a DNASiteObj is defined which inherits from the SiteObj
*/

#ifndef SITEOBJ_H
#define SITEOBJ_H

/*
includes
*/

#include <metasim.h>

class SiteObj {
  char state; 
public:
  SiteObj (char newstate = ' ');   // default state is a blank
  ~SiteObj ();
  
  char GetState ();                 // returns the state 

  void SetState (char newstate);    // sets state

  void Mutate();      // changes the site according to some function
  
  //overload the equality binary operators
  int operator==(SiteObj site);
  int operator!=(SiteObj site);
  
  friend ostream &operator<<(ostream &stream, SiteObj & site);
  friend istream &operator>>(istream &stream, SiteObj & site);

}; // end SiteObj

#endif /*SITEOBJ*/

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; End: ***
*/
