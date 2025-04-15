/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#include <iostream>
#include <limits>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include "gridNode.h"

GridNode::GridNode()
  : m_Orig_detID(-1),
    m_detID(-1),
    m_active(false),
    m_x(0.00),
    m_xDet(0.00),
    m_y(0.00),
    m_yDet(0.00),
    m_z(0.00),
    m_zDet(0.00),
    m_r(0.00),
    m_thetaRad(0.00),
    m_thetaDeg(0.00),
    m_WireDirection(TVector3(0.0, 0.0, 0.0)),
    m_halfLength(0.00),
    m_Slope(0.00),
    m_Sector(0),
    m_Layer(0),
    m_LayerLimit(false),
    m_SectorLimit(0),
    m_neighborsOri(std::vector<int>()),
    m_neighbors(std::vector<int>()),
    m_type(UNKNOWN),
    m_weight(1),
    m_parent(-1),
    m_cm(std::vector<int>())
{}

GridNode::GridNode( int det_id, bool const act, float const xin, float const yin, float const zin,
		    TVector3 const& Wdirect, float const halfLength, size_t sector, size_t layer,
		    bool const LayerLimit, short const SectorLimit,
                    std::vector<int> const & nghIDs, DetectorType t, short const w)
  : m_Orig_detID(det_id),
    m_detID(det_id),
    m_active(act),
    m_x(xin),
    m_xDet(xin),
    m_y(yin),
    m_yDet(yin),
    m_z(zin),
    m_zDet(-100.0),
    m_r(0.00),
    m_thetaRad(0.00),
    m_thetaDeg(0.00),
    m_WireDirection(Wdirect),
    m_halfLength(halfLength),
    m_Slope(0.0),
    m_Sector(sector),
    m_Layer(layer),
    m_LayerLimit(LayerLimit),
    m_SectorLimit(SectorLimit),
    m_neighborsOri(nghIDs),
    m_neighbors(nghIDs),
    m_type(t),
    m_weight(w),
    m_parent(-1),
    m_cm(std::vector<int>())
{
  ComputeSlope();
}

GridNode::GridNode (GridNode const &ot)
{
  this->m_Orig_detID = ot.m_Orig_detID;
  this->m_detID = ot.m_detID;
  this->m_active = ot.m_active;
  this->m_x = ot.m_x;
  this->m_xDet = ot.m_xDet;
  this->m_y = ot.m_y;
  this->m_yDet = ot.m_yDet;
  this->m_z = ot.m_z;
  this->m_zDet = ot.m_zDet;
  this->m_r = ot.m_r;
  this->m_thetaRad=ot.m_thetaRad;
  this->m_thetaDeg=ot.m_thetaDeg;
  this->m_WireDirection = TVector3(ot.m_WireDirection);
  this->m_halfLength = ot.m_halfLength;
  this->m_Slope = ot.m_Slope;
  this->m_Sector = ot.m_Sector;
  this->m_Layer = ot.m_Layer;
  this->m_LayerLimit = ot.m_LayerLimit;
  this->m_SectorLimit = ot.m_SectorLimit;
  this->m_neighborsOri = ot.m_neighborsOri;
  this->m_neighbors = ot.m_neighbors;
  this->m_type = ot.m_type;
  this->m_weight = ot.m_weight;
  this->m_parent = ot.m_parent;
  this->m_cm = ot.m_cm;
}

GridNode::~GridNode()
{}

GridNode& GridNode::operator=(GridNode const &ot)
{
  if( this != &ot ) {
    // Copy (deep)
    this->m_Orig_detID = ot. m_Orig_detID;
    this->m_detID      = ot.m_detID;
    this->m_active     = ot.m_active;
    this->m_x          = ot.m_x;
    this->m_xDet       = ot.m_xDet;
    this->m_y          = ot.m_y;
    this->m_yDet       = ot.m_yDet;
    this->m_z          = ot.m_z;
    this->m_zDet       = ot.m_zDet;
    this->m_r          = ot.m_r;
    this->m_thetaRad   = ot.m_thetaRad;
    this->m_thetaDeg   = ot.m_thetaDeg;
    this->m_WireDirection = TVector3(ot.m_WireDirection);
    this->m_halfLength = ot.m_halfLength;
    this->m_Slope      = ot.m_Slope;
    this->m_Sector     = ot.m_Sector;
    this->m_Layer      = ot.m_Layer;
    this->m_LayerLimit = ot.m_LayerLimit;
    this->m_SectorLimit = ot.m_SectorLimit;
    this->m_neighborsOri  = std::vector<int>(ot.m_neighborsOri);
    this->m_neighbors  = std::vector<int>(ot.m_neighbors);
    this->m_type       = ot.m_type;
    this->m_weight     = ot.m_weight;
    this->m_parent      = ot.m_parent;
    this->m_cm         = ot.m_cm;
  }
  return (*this);
}

void GridNode::resetNode()
{}


/* Operator "<". */
bool GridNode::operator<(GridNode const &ot) const
{
  if ( this->m_x < ot.m_x ) {
    return true;
  }
  if( !(this->m_x > ot.m_x) && !(this->m_x < ot.m_x)) { // equal
    return (this->m_y < ot.m_y);
  }
  return false;
}
/* Operator ">". */
bool GridNode::operator>(GridNode const &ot) const
{
  if ( this->m_x > ot.m_x ) {
    return true;
  }
  if( !(this->m_x < ot.m_x) && !(this->m_x > ot.m_x)) { // equal
    return (this->m_y > ot.m_y);
  }
  return false;
}

bool GridNode::IsNeighboring(int tubeID)
{
  std::vector<int>::iterator it;
  it = std::find (m_neighbors.begin(), m_neighbors.end(), tubeID);
  return ( it != m_neighbors.end() );
}

void GridNode::ComputeSlope() {
  //  m_Slope is already set to 0.00
  if( m_WireDirection.Phi() > 0.0 ) {// +3 degrees
    m_Slope = tan( (3.0 * M_PI)/180.00 );
  }
  else if( m_WireDirection.Phi() < 0.0 ) {// -3 degrees
    m_Slope = tan( (-3.0 * M_PI)/180.00 );
  }
  else {// Otherwise the slope is zero with respect to the z-axis
    m_Slope = 1.00;
  }
}

float GridNode::GetDistance2D(GridNode const &ot) const
{
  return sqrt( ( (this->m_x - ot.m_x) * (this->m_x - ot.m_x) ) +
	       ( (this->m_y - ot.m_y) * (this->m_y - ot.m_y) ) );
}


//_____________________________________________________________________________
// Debug functions
void GridNode::PrintNode() const
{
  std::string nodeType = "";
  if( m_type == STT_TYPE_PARA ) {
    nodeType = "Axial";
  }
  else if( m_type == STT_TYPE_SKEW ) {
    nodeType = "Skewed";
  }
  else {
    nodeType = "Virtual";
    std::cout << "Neighbour list: ";
    for(size_t i = 0; i < m_neighborsOri.size(); ++i) {
      std::cout << m_neighborsOri[i] << ", ";
     }
  }
  std::cout << " ID = " << m_detID
	    << " type = " <<  nodeType
	    << " layer = " << m_Layer
	    << " x = "     << this->m_x
	    << " y = "     << this->m_y
	    <<'\n';
}
