/**************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#pragma once
#ifndef PATH_AUXILIARY_FUNCTIONS_H
#define PATH_AUXILIARY_FUNCTIONS_H

// Start id for virtual nodes
#define START_VIRTUAL_ID 6000

/**
 * Data structure to hold the fit parameters after fitting a circle to
 * determine the curvature of a track candidate.
 */
struct CurvatureParameters{
  // public:
  // Constructors
CurvatureParameters()
: m_a(0.0), m_b(0.0), m_ra(0.0), m_r(0.0), m_E(0.0)
  {};
  
  explicit CurvatureParameters(double a, double b, double ra, double r, double e)
    : m_a(a), m_b(b), m_ra(ra), m_r(r), m_E(e)
  {};
  // D-tor
  virtual ~CurvatureParameters(){};
  
  // Copy constructor
CurvatureParameters(CurvatureParameters const &ot)
: m_a(ot.m_a), m_b(ot.m_b), m_ra(ot.m_ra), m_r(ot.m_r), m_E(ot.m_E)
  {};
  // Assignment operator
  CurvatureParameters& operator=(CurvatureParameters const &ot)
  {
    if(this != &ot) {
      this->m_a = ot.m_a;
      this->m_b = ot.m_b;
      this->m_ra = ot.m_ra;
      this->m_r = ot.m_r;
      this->m_E = ot.m_E;
    }
    return (*this);
  };
  // data Members
  double m_a;// x coord of circle center
  double m_b;// y_coord of circle center
  double m_ra;//(radius)
  double m_r;// (1/Radius)
  double m_E;// Energy value(minimized)
  // protected:
private:
  // Operators to avoid mistakes
  bool operator==(CurvatureParameters const &ot) const;
  bool operator>(CurvatureParameters  const &ot) const;
  bool operator<(CurvatureParameters  const &ot) const;
};
//------------------- END class CurvatureParameters -------------
#endif