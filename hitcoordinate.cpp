/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/

#include "hitcoordinate.h"

/**
 * Default Constructor
 */
HitCoordinate::HitCoordinate()
  : x(0.0), y(0.0), z(0.0),
    mx(0.0), my(0.0), mz(0.0),
    r(0.0), theta(0.0), thetaDeg(0.0),
    mr(0.0), mtheta(0.0), mthetaDeg(0.0),
    isochrone(0.0),
    type(UNKNOWN),
    m_detID(-1),
    m_EvtNum(0),
    m_trackID(-1),
    m_timeStamp(0.00)
{}

/**
 * Constructor
 */
HitCoordinate::HitCoordinate(float const val)
  : x(val), y(val), z(val),
    mx(val), my(val), mz(val),
    r(val), theta(val), thetaDeg(val),
    mr(val), mtheta(val), mthetaDeg(val),
    isochrone(val),
    type(UNKNOWN),
    m_detID(-1),
    m_EvtNum(0),
    m_trackID(val),
    m_timeStamp(0.00)
{}

/**
 * Destructor.
 */
HitCoordinate::~HitCoordinate()
{}

/**
 * Copy constructor.
 */
HitCoordinate::HitCoordinate( HitCoordinate const& ot)
  : x(ot.x), y(ot.y), z(ot.z),
    mx(ot.mx), my(ot.my), mz(ot.mz),
    r(ot.r), theta(ot.theta), thetaDeg(ot.thetaDeg),
    mr(ot.mr), mtheta(ot.mtheta), mthetaDeg(ot.mthetaDeg),
    isochrone(ot.isochrone),
    type(ot.type), m_detID(ot.m_detID),
    m_EvtNum(ot.m_EvtNum),
    m_trackID(ot.m_trackID),
    m_timeStamp(ot.m_timeStamp)
{}

/**
 * Assignment operator.
 */
HitCoordinate& HitCoordinate::operator=(HitCoordinate const& ot)
{
  // check for self-assignment
  if( this != &ot ) {
    // Copy (deep)
    this->x         = ot.x;
    this->y         = ot.y;
    this->z         = ot.z;
    this->mx        = ot.mx;
    this->my        = ot.my;
    this->mz        = ot.mz;
    this->r         = ot.r;
    this->theta     = ot.theta;
    this->thetaDeg  = ot.thetaDeg;
    this->mr        = ot.mr;
    this->mtheta    = ot.mtheta;
    this->mthetaDeg = ot.mthetaDeg;
    this->isochrone = ot.isochrone;
    this->type      = ot.type;
    this->m_detID   = ot.m_detID;
    this->m_EvtNum  = ot.m_EvtNum;
    this->m_trackID = ot.m_trackID;
    this->m_timeStamp = ot.m_timeStamp;
  }
  return (*this);
}
