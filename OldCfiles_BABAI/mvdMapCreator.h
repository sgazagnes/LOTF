/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#pragma once
#ifndef MVD_MAP_CREATOR_H
#define MVD_MAP_CREATOR_H

struct MVDMapCreator {
 public:
  /**
   * Constructor
   */
  MVDMapCreator();
  
  /**
   * Constructor.
   *@param geoFile Path to file containing the geometry.
   */
  explicit MVDMapCreator(std::string const &geoFile);
  
  /**
   * Destructor.
   */
  virtual ~MVDMapCreator();
  
  // Parameter members
  std::string m_geoFileName;
  
  // ___________ Protected ____
 protected:
  //___________ Private _______
 private:
  MVDMapCreator(MVDMapCreator const&ot);
  MVDMapCreator& operator=(MVDMapCreator const& ot);
};// END class definition
#endif
