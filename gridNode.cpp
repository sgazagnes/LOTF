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
    m_z_Det(0.00),
    m_WireDirection(TVector3(0.0, 0.0, 0.0)),
    m_halfLength(0.00),
    m_Slope(0.00),
    m_Sector(0),
    m_Layer(0),
    m_LayerLimit(false),
    m_SectorLimit(0),
    m_neighbors(std::vector<int>()),
    m_type(UNKNOWN),
    m_weight(1),
    m_length(0),
    m_lengthFW(0),
    m_lengthBW(0),
    m_forwardVisited(false),
    m_backwardVisited(false),
    m_area(0),
    m_visited(false),
    m_maxPathVisited(false),
    m_MVDAssigned(false),
    m_orientations(std::vector< NodeOrientation >()),
    m_OrientBackWard(std::vector< NodeOrientation >()),
    m_maxOrientVal(0),
    m_minOrientVal(0),
    m_orintVisited(false),
    m_maxOrientIndex(0),
    m_minOrientIndex(std::numeric_limits<size_t>::max()),
    m_maxBackWardOrientIndex(0),
    m_maxBackwardOrientVal(0),
    m_zestiVisited(false),
    m_times_visited(0),
    m_fitValue(0.00),
    m_mahalanobisDist(0.00),
    visited(0),
    parent(-1),
    m_cm(std::vector<int>())

{}

GridNode::GridNode( int det_id, bool const act, float const xin, float const yin, float const zin,
		    TVector3 const& Wdirect, float const halfLength, size_t sector, size_t layer,
		    bool const LayerLimit, int const SectorLimit,
                    std::vector<int> const & nghIDs, DetectorType t, short int const w)
  : m_Orig_detID(det_id),
    m_detID(det_id),
    m_active(act),
    m_x(xin),
    m_xDet(xin),
    m_y(yin),
    m_yDet(yin),
    m_z(zin),
    m_z_Det(-100.0),
    m_WireDirection(Wdirect),
    m_halfLength(halfLength),
    m_Slope(0.0),
    m_Sector(sector),
    m_Layer(layer),
    m_LayerLimit(LayerLimit),
    m_SectorLimit(SectorLimit),
    m_neighbors(nghIDs),
    m_type(t),
    m_weight(w),
    m_length(0),
    m_lengthFW(0),
    m_lengthBW(0),
    m_forwardVisited(false),
    m_backwardVisited(false),
    m_area(0),
    m_visited(false),
    m_maxPathVisited(false),
    m_MVDAssigned(false),
    m_orientations(std::vector< NodeOrientation >()),
    m_OrientBackWard(std::vector< NodeOrientation >()),
    m_maxOrientVal(0),
    m_minOrientVal(0),
    m_orintVisited(false),
    m_maxOrientIndex(0),
    m_minOrientIndex(std::numeric_limits<size_t>::max()),
    m_maxBackWardOrientIndex(0),
    m_maxBackwardOrientVal(0),
    m_zestiVisited(false),
    m_times_visited(0),
    m_fitValue(0.00),
    m_mahalanobisDist(0.00),
    visited(0),
    parent(-1),
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
  this->m_z_Det = ot.m_z_Det;
  this->m_WireDirection = TVector3(ot.m_WireDirection);
  this->m_halfLength = ot.m_halfLength;
  this->m_Slope = ot.m_Slope;
  this->m_Sector = ot.m_Sector;
  this->m_Layer = ot.m_Layer;
  this->m_LayerLimit = ot.m_LayerLimit;
  this->m_SectorLimit = ot.m_SectorLimit;
  this->m_neighbors = ot.m_neighbors;
  this->m_type = ot.m_type;
  this->m_weight = ot.m_weight;
  this->m_length = ot.m_length;
  this->m_lengthFW = ot.m_lengthFW;
  this->m_lengthBW = ot.m_lengthBW;
  this->m_forwardVisited = ot.m_forwardVisited;
  this->m_backwardVisited = ot.m_backwardVisited;
  this->m_area = ot.m_area;
  this->m_visited = ot.m_visited;
  this->m_maxPathVisited = ot.m_maxPathVisited;
  this->m_MVDAssigned = ot.m_MVDAssigned;
  this->m_orientations = ot.m_orientations;
  this->m_OrientBackWard = ot.m_OrientBackWard;
  this->m_maxOrientVal = ot.m_maxOrientVal;
  this->m_minOrientVal = ot.m_minOrientVal;
  this->m_orintVisited = ot.m_orintVisited;
  this->m_maxOrientIndex = ot.m_maxOrientIndex;
  this->m_minOrientIndex = ot.m_minOrientIndex;
  this->m_maxBackWardOrientIndex = ot.m_maxBackWardOrientIndex;
  this->m_maxBackwardOrientVal = ot.m_maxBackwardOrientVal;
  this->m_zestiVisited = ot.m_zestiVisited;
  this->m_times_visited = ot.m_times_visited;
  this->m_fitValue = ot.m_fitValue;
  this->m_mahalanobisDist = ot.m_mahalanobisDist;
  this->visited = ot.visited;
  this->parent = ot.parent;
  this->m_cm = ot.m_cm;
}

GridNode::~GridNode()
{}

GridNode& GridNode::operator=(GridNode const &ot)
{
  if( this != &ot ) {
    // Copy (deep)
    this->m_Orig_detID = ot. m_Orig_detID;
    this->m_detID     = ot.m_detID;
    this->m_active    = ot.m_active;
    this->m_x         = ot.m_x;
    this->m_xDet         = ot.m_xDet;
    this->m_y         = ot.m_y;
    this->m_yDet         = ot.m_yDet;
    this->m_z         = ot.m_z;
    this->m_z_Det     = ot.m_z_Det;
    this->m_WireDirection = TVector3(ot.m_WireDirection);
    this->m_halfLength = ot.m_halfLength;
    this->m_Slope      = ot.m_Slope;
    this->m_Sector     = ot.m_Sector;
    this->m_Layer      = ot.m_Layer;
    this->m_LayerLimit = ot.m_LayerLimit;
    this->m_SectorLimit = ot.m_SectorLimit;
    this->m_neighbors = std::vector<int>(ot.m_neighbors);
    this->m_type      = ot.m_type;
    this->m_weight    = ot.m_weight;
    this->m_length    = ot.m_length;
    this->m_lengthFW  = ot.m_lengthFW;
    this->m_lengthBW  = ot.m_lengthBW;
    this->m_forwardVisited  = ot.m_forwardVisited;
    this->m_backwardVisited = ot.m_backwardVisited;
    this->m_area      = ot.m_area;
    this->m_visited   = ot.m_visited;
    this->m_maxPathVisited = ot.m_maxPathVisited;
    this->m_MVDAssigned = ot.m_MVDAssigned;
    this->m_orientations = std::vector< NodeOrientation >(ot.m_orientations);
    this->m_OrientBackWard = std::vector< NodeOrientation >(ot.m_OrientBackWard);
    this->m_maxOrientVal = ot.m_maxOrientVal;
    this->m_minOrientVal = ot.m_minOrientVal;
    this->m_orintVisited = ot.m_orintVisited;
    this->m_maxOrientIndex = ot.m_maxOrientIndex;
    this->m_minOrientIndex = ot.m_minOrientIndex;
    this->m_maxBackWardOrientIndex = ot.m_maxBackWardOrientIndex;
    this->m_maxBackwardOrientVal = ot.m_maxBackwardOrientVal;
    this->m_zestiVisited = ot.m_zestiVisited;
    this->m_times_visited = ot.m_times_visited;
    this->m_fitValue = ot.m_fitValue;
    this->m_mahalanobisDist = ot.m_mahalanobisDist;
    this->visited = ot.visited;
    this->parent = ot.parent;
    this->m_cm = ot.m_cm;
  }
  return (*this);
}

void GridNode::resetNode()
{}

void GridNode::resetNodeOrientation()
{
  m_orientations.clear();
  m_OrientBackWard.clear();
  this->m_maxOrientVal   = 0;
  this->m_minOrientVal   = 0;
  this->m_maxOrientIndex = 0;
  this->m_minOrientIndex = 0;
  this->m_maxBackWardOrientIndex = 0;
  this->m_maxBackwardOrientVal = 0;
}

void GridNode::initNodeOrientation(size_t const numSlices)
{
  m_orientations.clear();
  m_OrientBackWard.clear();
  
  m_orientations   = std::vector< NodeOrientation >(numSlices);
  m_OrientBackWard = std::vector< NodeOrientation >(numSlices);

  for(size_t i = 0; i < m_orientations.size(); ++i) {
    m_orientations[i].m_angle  = 0;
    m_orientations[i].m_radius = 0;
    m_orientations[i].m_act    = false;
    (m_orientations[i].m_memberIds) = std::vector<int>();
    // For the backward(in case of two ways scanning)
    m_OrientBackWard[i].m_angle  = 0;
    m_OrientBackWard[i].m_radius = 0;
    m_OrientBackWard[i].m_act    = false;
    (m_OrientBackWard[i].m_memberIds) = std::vector<int>();
  }
  this->m_maxOrientVal = 0;
  this->m_minOrientVal = 0;
  this->m_maxBackwardOrientVal = 0;
}
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
  // for(size_t i = 0; i < m_neighbors.size(); ++i) {
  //   if( m_neighbors[i] == tubeID ) {
  //     return true;
  //   }
  // }
  // return false;
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

int GridNode::GetOrientationIndex(float val, float tol, bool OuterToInner) const
{
  int output_index = -1;
  /*
   * Note:: This can go wrong if number of orientation larger than
   * "max int"; but we do not expect that. It can practically never
   * happen.
   */
  int numberOfOrient = (OuterToInner)?(static_cast<int>(m_orientations.size())) :
    (static_cast<int>(m_OrientBackWard.size()));
  
  // Select the list to search
  std::vector< NodeOrientation > const &OrList = ((OuterToInner)?(this->m_orientations):(this->m_OrientBackWard));
  
  // If orientation spaces have been computed.
  if( numberOfOrient != 0 ) {
    int index = 0;
    while( index < numberOfOrient ) {
      if( fabs(OrList[index].m_angle - val) < tol) {
	// Found the orientation; assign and exit
	output_index = index;
	// Force to exit while loop
	index = (numberOfOrient + 1);
      }
      else {
        index++;
      }
    }
  }
  return output_index;
}

NodeOrientation& GridNode::GetOrientation(int index, bool OuterToInner)
{
  // Doe een index check
  return (OuterToInner ? m_orientations[index] : m_OrientBackWard[index]);
}

NodeOrientation const& GridNode::GetOrientation(int index, bool OuterToInner) const
{
  // Doe een index check
  return (OuterToInner ? m_orientations[index] : m_OrientBackWard[index]);
}

void GridNode::InsertOrientAtt(NodeOrientation const& att, bool OuterToInner)
{
  if( OuterToInner ){
    m_orientations.push_back(att);
  }
  else{
    m_OrientBackWard.push_back(att);
  }
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
    for(size_t i = 0; i < m_neighbors.size(); ++i) {
      std::cout << m_neighbors[i] << ", ";
    }
  }
  std::cout << " ID = " << m_detID
	    << " type = " <<  nodeType
	    << " num orient = " << m_orientations.size()
	    << " max orient = " << m_maxOrientVal
	    << " index = " << m_maxOrientIndex
	    << " layer = " << m_Layer
	    << " x = "     << this->m_x
	    << " y = "     << this->m_y
	    <<'\n';
}
