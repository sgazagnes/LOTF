/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#include "mvdMapCreator.h"

MVDMapCreator::MVDMapCreator()
  : m_geoFileName("")
{}

MVDMapCreator::MVDMapCreator(std::string const &geoFile)
  : m_geoFileName(geoFile)
{}

MVDMapCreator::~MVDMapCreator()
{}

