/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#pragma once
#ifndef PATH_CANDIDATE_H
#define PATH_CANDIDATE_H

#include <vector>
#include <set>

#include "auxiliaryfunctions.h"
#include "utilfunctions.h"
#include "gridNode.h"
#include "CoordGrid.h"

//_________________________  PathCandidate ________________
struct PathCandidate{
 public:
  // Constructors
  PathCandidate();
  
 // Destructor
  virtual ~PathCandidate();

  // Copy
  PathCandidate(PathCandidate const &ot);

  // Assignment
  PathCandidate& operator=(PathCandidate const &ot);

  // Member fuctions
  void updateHeadAndTailNodes();
  void insertNewNode(CoordGrid &gr, std::vector< GridNode > &Ingrid, GridNode *node,  std::vector<int>::iterator it);
  void insertNewNodeFinal(CoordGrid &gr, std::vector< GridNode > &Ingrid, GridNode *node,  std::vector<int>::iterator it);
  bool compareTwoPathsLength(PathCandidate *a, PathCandidate *b) ;
  void determineSkewedXY(CoordGrid &gr, std::vector< GridNode > &Ingrid, GridNode *node, size_t vecindex);
  void fEstimateSkewedXY(CoordGrid &gr, std::vector< GridNode > &Ingrid, GridNode *node, size_t vecindex);
  void addToAnchor(CoordGrid &gr, std::vector< GridNode > &Ingrid,GridNode *node, size_t vecindex, int prevLayer);
  void correctPrevAnchor(CoordGrid &gr, std::vector< GridNode > &Ingrid,GridNode *node, size_t vecindex);
  inline bool isInCandidate(int nodeId) const;

  // Comparison functions
  inline friend bool GreaterThanAngle (PathCandidate const &lf, PathCandidate const &rt);
  inline friend bool LessThanAngle (PathCandidate const &lf, PathCandidate const &rt);
  inline friend bool GreaterThanWeight(PathCandidate const &lf, PathCandidate const &rt);
  inline friend bool LessThanWeight(PathCandidate const &lf, PathCandidate const &rt);
  inline friend bool GreaterThanLength(PathCandidate const &lf, PathCandidate const &rt);
  inline friend bool LessThanLength(PathCandidate const &lf, PathCandidate const &rt);
  inline friend bool GreaterThanMaxLayer(PathCandidate const &lf, PathCandidate const &rt);
  inline friend bool GreaterThanMaxLayerP(PathCandidate const *lf, PathCandidate const *rt);
  inline friend bool GreaterThanMinLayer(PathCandidate const &lf, PathCandidate const &rt);
  inline friend bool LessThanMaxLayer(PathCandidate const &lf, PathCandidate const &rt);
  inline friend bool LessThanMinLayer(PathCandidate const &lf, PathCandidate const &rt);
  // bool compareTwoPathsLength(PathCandidate lf, PathCandidate rt);
  bool compareTwoPathsLength (PathCandidate const *lf, PathCandidate const *rt);
  // Member variables and containers
  unsigned int m_id;
  unsigned int m_level;
  unsigned int m_weight;// Number of member nodes
  float  m_orientation;// The orientation in radians

  bool  m_OuterToInner;// Direction of the path
  unsigned int m_length;// Length of the path (>= 1)
  bool m_isValid;
  bool m_isShort;// If it has less length(weight) than min
  bool m_isMerged;
  bool m_hasCurvature;

  CurvatureParameters m_CurV_par;// Estimate of curvature

  // Ids of the head and tail nodes. Non virtual
  int    m_maxlayerNodeId;//id of node with max layer
  int    m_minlayerNodeId;//id of node with min layer
  size_t m_maxLayernodeLayer;// Layer number of max
  size_t m_minlayerNodeLayer;// Layer number of min


  int    m_maxLayerNodeId;//id of node with max layer
  int    m_minLayerNodeId;//id of node with min layer
  int    m_maxLayer; //Layer number of max//id of node with max layer
  int    m_minLayer;// Layer number of min//id of node with min layer
  int    m_nPerLayer;
  int    m_lastNodeVirtual;//id of node with min layer
  // int    m_lastNodeVirtualId;//id of node with min layer

  int    m_finished;
  bool   m_isOnSectorLimit;
  bool   m_startOnSectorLimit;

  int    m_headNode;
  int    m_tailNode;
  int    m_lastNodeAdded;
  int    m_seenVirtual;//id of node with min layer
  int    m_lastVirtual;
  int    m_numVirtual;
  
  std::vector<int>       *m_memberList;// List of members in a vector(delete me)
  std::set<int>          *m_memberIdSet;// The set of member ids
  std::vector<int> 	 m_prevVirtuals;

  std::vector<unsigned int>	 m_toMergeHead;
  std::vector<unsigned int> 	 m_toMergeTail;
  std::vector<unsigned int> 	 m_sectors;
  std::vector<int> 	         m_layers;
  
  std::vector<int> m_headNeigh;
  std::vector<int> m_tailNeigh;

  std::vector<unsigned int> 	 m_listSkewed;
  std::vector<double> m_x;// List of x
  std::vector<double> m_y;// List of y
  std::vector<double> m_z;// List of z
  std::vector<double> m_r;// List of r
  std::vector<double> m_theta;// List of theta
  std::vector<GridNode> m_anchors;// List of theta

  /* Shape and spread parameters */
  
  std::vector< float > m_meanVector;
  std::vector< std::vector<float> > m_covMatrix;
  std::vector< std::vector<float> > m_inveCovMatrix;
  float m_eigevValue1;
  float m_eigevValue2;
  std::vector<float> m_eigenVector1;
  std::vector<float> m_eigenVector2;
  // The list of parents and childeren
  std::vector<unsigned int> m_parents;
  std::vector<unsigned int> m_childeren;


  //_________ Protected functions and members _______
  //protected:

  //_________ Private functions and members _______
  private:
  bool operator==(PathCandidate const &ot) const;
  bool operator<(PathCandidate const &ot) const;
  bool operator>(PathCandidate const &ot) const;
};
//________________________ END PathCandidate ______________
// Inline implementations
bool GreaterThanAngle (PathCandidate const &lf, PathCandidate const &rt)
{
  return (lf.m_orientation > rt.m_orientation);
}

bool LessThanAngle (PathCandidate const &lf, PathCandidate const &rt)
{
  return (lf.m_orientation < rt.m_orientation);
}

bool GreaterThanWeight(PathCandidate const &lf, PathCandidate const &rt)
{
  return (lf.m_weight > rt.m_weight);
}

bool LessThanWeight(PathCandidate const &lf, PathCandidate const &rt)
{
  return (lf.m_weight < rt.m_weight);
}

bool GreaterThanLength(PathCandidate const &lf, PathCandidate const &rt)
{
  return (lf.m_length > rt.m_length);
}

bool LessThanLength(PathCandidate const &lf, PathCandidate const &rt)
{
  return (lf.m_length > rt.m_length);
}

bool GreaterThanMaxLayer(PathCandidate const &lf, PathCandidate const &rt)
{
  return (lf.m_maxLayernodeLayer > rt.m_maxLayernodeLayer);
}
bool GreaterThanMaxLayerP(PathCandidate const *lf, PathCandidate const *rt)
{
  return ( GreaterThanMaxLayer( (*lf), (*rt) ) );
}
bool LessThanMaxLayer(PathCandidate const &lf, PathCandidate const &rt)
{
  return (lf.m_maxLayernodeLayer < rt.m_maxLayernodeLayer);
}

bool GreaterThanMinLayer(PathCandidate const &lf, PathCandidate const &rt)
{
  return (lf.m_minlayerNodeLayer > rt.m_minlayerNodeLayer);
}

bool LessThanMinLayer(PathCandidate const &lf, PathCandidate const &rt)
{
  return (lf.m_minlayerNodeLayer < rt.m_minlayerNodeLayer);
}
bool PathCandidate::isInCandidate(int nodeId) const
{
  return( (m_memberIdSet->find(nodeId)) != m_memberIdSet->end() );
}

bool compareTwoPathsLength (PathCandidate const *lf, PathCandidate const *rt)
{ 
  
    return (lf->m_length > rt->m_length); 
} 

#endif// END header
