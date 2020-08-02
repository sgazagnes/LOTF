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
  void insertNewNode(GridNode *node);

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
  
  // Member variables and containers
  unsigned int m_id;
  unsigned int m_level;
  unsigned int m_weight;// Number of member nodes
  float  m_orientation;// The orientation in radians
  std::vector<int> *m_memberList;// List of members in a vector(delete me)
  std::set<int>    *m_memberIdSet;// The set of member ids
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
  int    m_maxLayer;//id of node with max layer
  int    m_minLayer;//id of node with min layer

  int    m_lastNodeVirtual;//id of node with min layer
  int    m_lastNodeVirtualId;//id of node with min layer

  bool   m_finished;
  bool   m_isOnSectorLimit;
  bool   m_startOnSectorLimit;

  int    m_lastNodeVisited;
  int    m_firstNodeVisited;

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
#endif// END header
