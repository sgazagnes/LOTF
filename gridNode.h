/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#pragma once
#ifndef GRID_NODE_CLASS_H
#define GRID_NODE_CLASS_H

#include <vector>
#include <set>
#include <utility>
#include "TVector3.h"

//_______________________________________
struct NodeOrientation {
  // public:

  // Constructor
NodeOrientation()
: m_angle(0.0),
    m_radius(0.0),
    m_act(false),
    m_weight(0),
    m_val(1),
    m_memberIds(std::vector<int>())
  {};

  // Copy constructor
NodeOrientation(NodeOrientation const &ot)
: m_angle(ot.m_angle),
    m_radius(ot.m_radius),
    m_act(ot.m_act),
    m_weight(ot.m_weight),
    m_val(ot.m_val),
    m_memberIds( std::vector<int>(ot.m_memberIds) )
  {};
  
  // Deep copy assignment operator
  inline NodeOrientation& operator=(NodeOrientation const &ot);

  // Destructor  
  virtual ~NodeOrientation(){};

  // Compare functions
  inline friend bool GrThanOrient  (NodeOrientation const &lf, NodeOrientation const &rt);
  inline friend bool LessThanOrient(NodeOrientation const &lf, NodeOrientation const &rt);
  
  //___________ Members
  // Radians (the angle of the current orientation)
  float m_angle;
  // The radius of the constructed set of points participating in the
  // current orientation
  float m_radius;
  // If active for the current angle
  bool m_act;
  // Weight of the orientation(number of members)
  unsigned int m_weight;
  // Value for this orientation. Positioning order in layer-based
  // scheme.
  unsigned int m_val;
  // The list of node id's in the current orientation
  std::vector<int> m_memberIds;
  //protected:
  
private:
  bool operator>(NodeOrientation const &ot);
  bool operator<(NodeOrientation const &ot);
};// END OF NodeOrientation
/////////////// NodeOrientation ///////////////////
NodeOrientation& NodeOrientation::operator=(NodeOrientation const &ot)
{
  if( this != &ot ) {
    this->m_angle     = ot.m_angle;
    this->m_radius    = ot.m_radius;
    this->m_act       = ot.m_act;
    this->m_weight    = ot.m_weight;
    this->m_val       = ot.m_val;
    // Maybe a for loop, but for now, trust copy const. of vector.
    this->m_memberIds = std::vector<int>(ot.m_memberIds);
  }
  return (*this);
}

bool GrThanOrient(NodeOrientation const &lf, NodeOrientation const &rt)
{
  return (lf.m_angle > rt.m_angle);
}

bool LessThanOrient(NodeOrientation const &lf, NodeOrientation const &rt)
{
  return (lf.m_angle < rt.m_angle);
}
//************************************************************
//_________________ Grid node interface definition ___________
struct GridNode {
  // public:
  /**
   * Type of the stored coordinates.
   */
  typedef enum Detector_Type {
    UNKNOWN             = 0,
    STT_TYPE_PARA       = 1,
    STT_TYPE_SKEW       = 2,
    STT_TYPE_SPLIT_SKEW = 3,
    VIRTUAL_NODE        = 4,
    VIRTUAL_NODE_SPLIT  = 5,
    MVD_TYPE_PIXEL      = 6,
    MVD_TYPE_STRIP      = 7
  } DetectorType;
  /*
   * SK_V_SK: Skewd <-> Virtual <-> Skewed
   * AX_V_SK: Axial <-> Virtual <-> Skewed
  */
  typedef enum Virtual_Type {
    SK_V_SK     = 0,
    AX_V_SK     = 1
  } VirtualTubeType;
  
  // Constructors
  GridNode();// Default

  explicit GridNode( int det_id, bool const act,
		     float const xin, float const yin, float const zin,
		     TVector3 const& Wdirect, float const halfLength,
		     size_t sector, size_t layer,
		     bool const LayerLimit, int const SectorLimit,
		     std::vector<int> const &nghIDs, DetectorType type,
		     short int const w);
  // Copy
  GridNode (GridNode const &ot);

  // Dtor
  virtual ~GridNode();

  // Functions
  void resetNode();
  void resetNodeOrientation();
  void initNodeOrientation(size_t const numSlices = 0);

  // Operators
  GridNode& operator=(GridNode const &ot);
  bool operator<(GridNode const &ot) const;
  bool operator>(GridNode const &ot) const;

  // Functions to sort based on detector ID
  inline friend bool LessThanID(GridNode const &lf, GridNode const &rt);
  inline friend bool GreaterThanID(GridNode const &lf, GridNode const &rt);
  // Orientation based
  inline friend bool LessThanOrientLength(GridNode const &lf, GridNode const &rt);
  inline friend bool GreaterThanOrientLength(GridNode const &lf, GridNode const &rt);
  //////////// Sorting on layer numbers.
  inline friend bool GrThanLayer(GridNode const &lf, GridNode const &rt);
  inline friend bool LessThanLayer(GridNode const &lf, GridNode const &rt);
  //////////////// Sector based
  inline friend bool GrThanSec(GridNode const &lf, GridNode const &rt);
  inline friend bool LessThanSec(GridNode const &lf, GridNode const &rt);
  inline friend bool GrThanLayerSec(GridNode const &lf, GridNode const &rt);

  // Neighborings and distance in 2D(Z dimension is not defined
  // properly)
  inline std::vector<int> const& GetNeighbors() const;
  inline std::vector<int>& GetNeighbors();

  bool IsNeighboring(int tubeID);
  
  float GetDistance2D(GridNode const &ot) const;
  
  // Orientation related functions
  inline std::vector<NodeOrientation> const& GetOrientations(bool OuterToInner = true) const;
  inline std::vector<NodeOrientation>& GetOrientations(bool OuterToInner = true);
  /**
   * Return the index of the orientation with the given value "val".
   * Return -1 if not found.
   *@param val Angle of the orientation (rad).
   *@param OuterToInner True, search in outer to inner list. False,
   * search in inner to outer orientations.
   *@return -1 if not found, index of orientation otherwise.
   */
  int GetOrientationIndex(float val, float tol, bool OuterToInner = true) const;

  /**
   * Get a reference to the orientation attribute with a give index.
   *@param index The index of the orientation attribute.
   *@param OuterToInner The list direction
   *@return Reference the the Attribute
   */
  NodeOrientation& GetOrientation(int index, bool OuterToInner = true);
  NodeOrientation const& GetOrientation(int index, bool OuterToInner = true) const;
  
  void InsertOrientAtt(NodeOrientation const& att, bool OuterToInner = true);
  
  // Type functions
  inline bool IsSTTSplitSkewedNode() const;
  inline bool IsVirtualSplitNode() const;
  
  // Tube properties.
  inline float GetHalfLength() const;
  void ComputeSlope();
   
  // Debug functions
  void PrintNode() const;
  
  //______________________ Members
  int   m_Orig_detID;// The detector ID in the original map(geometry)
  int   m_detID; // unique ID of the detector(ID after adding virtuals)
  bool  m_active;// Has fired
  float m_x;    // X ccordinate
  float m_xDet;// X ccordinate(determined)
  float m_y;   // Y coordinate
  float m_yDet;// Y coordinate(determined)
  float m_z;   // Z coordinate (tube half length)
  float m_z_Det; // Z coordinate (determined)
  TVector3 m_WireDirection;// Direction of the internal wire
  float  m_halfLength;// Tube half length
  float  m_Slope;// Slope of the tube in mounting frame.
  size_t m_Sector;//Sector (note this is too large we may use unsigned short int)
  size_t m_Layer;//Layer (note this is too large we may use unsigned short int)
  bool   m_LayerLimit;// If Layer limit
  int    m_SectorLimit;// If sector limit
  // List of detector id's of the direct neigbors.
  std::vector<int> m_neighbors;
  DetectorType m_type;// Node type
  short int m_weight;//
  size_t m_length;// Total Path length
  size_t m_lengthFW;// Forward path length
  size_t m_lengthBW;// Backward path length
  bool   m_forwardVisited;// Visited in forward direction
  bool   m_backwardVisited;// Visited in backward direction
  size_t m_area;// Size of Area (opening)
  bool   m_visited;// If has been added to area
  bool   m_maxPathVisited;// Used for maximum path for a given direction.
  bool   m_MVDAssigned;// If mvd node and if it has been assigned to a track
  /* 
   * Set of orientations in which this node is participating.
   */
  std::vector< NodeOrientation > m_orientations;// OuterToInner
  std::vector< NodeOrientation > m_OrientBackWard;// InnerToOuter
  size_t m_maxOrientVal;// Maximum responce value
  size_t m_minOrientVal;// Minimum responce value
  bool   m_orintVisited;// If node has been visited during current Att. Space computation
  size_t m_maxOrientIndex;// Index of the maximum responce
  size_t m_minOrientIndex;// Index of the minimum responce
  unsigned int m_maxBackWardOrientIndex;
  unsigned int m_maxBackwardOrientVal;
  
  // If visited for z-estimation
  bool m_zestiVisited;
  // How often a node is visited
  size_t m_times_visited;
  // Holding the goodness of the fit for the current node. How good
  // this node fits in a given tracklet.
  double m_fitValue;
  float m_mahalanobisDist;
  int visited;
  int   parent;
  std::vector<int> m_cm;


  //_________****************___________ 
  // protected:
  
 private:
  bool operator==(GridNode const &ot) const;
  bool operator!=(GridNode const &ot) const;
};// END of class definition

//______________ Inline functions
inline float GridNode::GetHalfLength() const
{
  return this->m_halfLength;
}

inline std::vector<int> const& GridNode::GetNeighbors() const
{
  return this->m_neighbors;
}

inline std::vector<int> &GridNode::GetNeighbors()
{
  return this->m_neighbors;
}

inline std::vector<NodeOrientation> const& GridNode::GetOrientations(bool OuterToInner) const
{
  return (OuterToInner ? (this->m_orientations) : (this->m_OrientBackWard));
}

inline std::vector<NodeOrientation>& GridNode::GetOrientations(bool OuterToInner)
{
  return (OuterToInner ? (this->m_orientations) : (this->m_OrientBackWard));
}

bool GridNode::IsSTTSplitSkewedNode() const
{
  return (this->m_type == GridNode::STT_TYPE_SPLIT_SKEW);
}

bool GridNode::IsVirtualSplitNode() const
{
  return (this->m_type == GridNode::VIRTUAL_NODE_SPLIT);
}

// Functions to sort based on detector ID
bool LessThanID(GridNode const &left, GridNode const &right)
{
  return (left.m_detID < right.m_detID);
}

bool GreaterThanID(GridNode const &left, GridNode const &right)
{
  return (left.m_detID > right.m_detID);
}

bool LessThanOrientLength(GridNode const &left, GridNode const &right)
{
  return (left.m_maxOrientVal < right.m_maxOrientVal);
}

bool GreaterThanOrientLength(GridNode const &left, GridNode const &right)
{
  return (left.m_maxOrientVal > right.m_maxOrientVal);
}

bool GrThanLayer(GridNode const &lf, GridNode const &rt)
{
  return (lf.m_Layer > rt.m_Layer );
}

bool LessThanLayer(GridNode const &lf, GridNode const &rt)
{
  return (lf.m_Layer < rt.m_Layer );
}

bool GrThanSec(GridNode const &lf, GridNode const &rt)
{
  return (lf.m_Sector > rt.m_Sector);
}

bool LessThanSec(GridNode const &lf, GridNode const &rt)
{
  return (lf.m_Sector < rt.m_Sector);
}

bool GrThanLayerSec(GridNode const &lf, GridNode const &rt)
{
  // if(lf.m_Layer > rt.m_Layer){
  //  return true;
  // }
  // if(lf.m_Layer == rt.m_Layer) {
  //  return (lf.m_Sector > rt.m_Sector);
  // }
  // return false;
  return (lf.m_Layer == rt.m_Layer) ? (lf.m_Sector > rt.m_Sector) : (lf.m_Layer > rt.m_Layer);
}
#endif//END of interface
