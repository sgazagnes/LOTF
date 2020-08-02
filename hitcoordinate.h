/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#pragma once
#ifndef HIT_COORDINATE_CLASS_H
#define HIT_COORDINATE_CLASS_H

struct HitCoordinate {
public:
  /**
   * Type of the stored coordinates.
   */
  typedef enum Hit_Detector {
    UNKNOWN  = 0,
    STT_TYPE = 1,
    MVD_TYPE = 2,
  } Hit_DetectorType;
  
  /**
   * Default Constructor
   */
  HitCoordinate();
  
  /**
   * Constructor
   */
  explicit HitCoordinate(float const val);
  
  /**
   * Destructor.
   */
  virtual ~HitCoordinate();

  /**
   * Copy constructor.
   */
  HitCoordinate( HitCoordinate const& ot);
  
  /**
   * Assignment operator.
   */
  HitCoordinate& operator=(HitCoordinate const& ot);
  
  // Coordinates of the current hit
  float x;
  float y;
  float z;
  // MC coordinates of the current hit
  float  mx;
  float  my;
  float  mz;
  // Polar Coordinates
  float r;
  double theta;
  double thetaDeg;
  // MC Polar Coordinates
  float  mr;
  double mtheta;
  double mthetaDeg;
  // If Stt the isochrone
  double isochrone;
  Hit_DetectorType type;
  // Detector id for now just STT, for others we should avoid
  int m_detID;
  int m_EvtNum;// Event number
  // MC trackID. If = -1, it does not belong to any MC track
  int m_trackID;
  float m_timeStamp;// Time based
};
#endif
