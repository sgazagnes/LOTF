/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#pragma once
#ifndef TRACK_OBJECT_H
#define TRACK_OBJECT_H

//#include <cstdlib>
#include <set>
#include <vector>

#include "utilfunctions.h"

struct TrackObject {
 public:
  // Constructors
 TrackObject();
 
 /* explicit TrackObject(std::set<int> const &sttComp, std::vector<int>
    const &MVD_Comp);*/
 
 // Destructor
  virtual ~TrackObject();

  std::set<int>* m_sttComponent;
  std::vector<int> m_MVD_Component;
  
  //_________ Protected functions and members _______
  //protected:

  //_________ Private functions and members _______
  private:
  TrackObject(TrackObject const &ot);
  TrackObject& operator=(TrackObject const &ot);
  bool operator==(TrackObject const &ot) const;
  bool operator<(TrackObject const &ot) const;
  bool operator>(TrackObject const &ot) const;
};// End interface Trackobjects

///// _____ Interface MC-Track object def. ________
struct MCTrackObject {
 public:
  // Constructors
  MCTrackObject()
  : m_pointSTTCoordList(std::vector < point3D >()),
    m_pointMVDCoordList(std::vector < point3D >()),
    m_STT_Component(std::vector<int>()),
    m_MVD_Component(std::vector<int>()),
    m_matched(false)
  {};
  //Dtor
  virtual ~MCTrackObject(){};

  // member functions and operators
  inline bool operator< (MCTrackObject const &ot) const;
  inline bool operator> (MCTrackObject const &ot) const;
  inline friend bool lessThanLength(MCTrackObject const *lt, MCTrackObject const *rt);
  inline friend bool greaterThanLength(MCTrackObject const *lt, MCTrackObject const *rt);

  /* Dumy Remove me later */
  void print() const {};
  /* End of dummy functions ... */
  
  // member parameters.
  std::vector < point3D > m_pointSTTCoordList;
  std::vector < point3D > m_pointMVDCoordList;

  std::vector<int> m_STT_Component;
  std::vector<int> m_MVD_Component;

  bool m_matched;// If have been matched before
  //_________ Protected functions and members _______
  //protected:

  //_________ Private functions and members _______
  private:
  MCTrackObject(MCTrackObject const &ot);
  MCTrackObject& operator=(MCTrackObject const &ot);
  bool operator==(MCTrackObject const &ot) const;
};// End Interface MC-Track
//=================================================================
//____________ Implementations ________
bool MCTrackObject::operator> (MCTrackObject const &ot) const
{
  return ( (this->m_STT_Component).size() > (ot.m_STT_Component).size() );
}
bool MCTrackObject::operator< (MCTrackObject const &ot) const
{
  return ( (this->m_STT_Component).size() < (ot.m_STT_Component).size() );
}
bool lessThanLength(MCTrackObject const *lt, MCTrackObject const *rt)
{
  return ( (lt->m_STT_Component).size() < (rt->m_STT_Component).size() );
}
bool greaterThanLength(MCTrackObject const *lt, MCTrackObject const *rt)
{
  return ( (lt->m_STT_Component).size() > (rt->m_STT_Component).size() );
}
#endif
