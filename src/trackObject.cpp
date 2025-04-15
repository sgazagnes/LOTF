/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#include "trackObject.h"

TrackObject::TrackObject()
  // : m_sttComponent(new std::set<int>()),
  //   m_MVD_Component(std::vector<int>())
{}

TrackObject::~TrackObject()
{
  m_sttComponent->clear();
  delete (m_sttComponent);
}
