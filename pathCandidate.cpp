/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#include <iterator>
#include "pathCandidate.h"

// Constructors
PathCandidate::PathCandidate()
  : m_id(0),
    m_level(0),
    m_weight(0),
    m_orientation(0.00),
    m_memberList(new std::vector<int>()),
    m_memberIdSet(new std::set<int>()),
    m_OuterToInner(true),
    m_length(1),
    m_isValid(true),
    m_isShort(false),
    m_isMerged(false),
    m_hasCurvature(false),
    m_CurV_par(CurvatureParameters()),
    m_maxlayerNodeId(-1),
    m_minlayerNodeId(-1),
    m_maxLayernodeLayer(100),
    m_minlayerNodeLayer(100),
    m_meanVector(std::vector<float>(2, 0.0)),
    m_covMatrix(std::vector< std::vector<float> >(2, std::vector<float>(2, 0.0))),
    m_inveCovMatrix(std::vector< std::vector<float> >(2, std::vector<float>(2, 0.0))),
    m_eigevValue1(0.0),
    m_eigevValue2(0.0),
    m_eigenVector1(std::vector<float>(2,0.0)),
    m_eigenVector2(std::vector<float>(2,0.0)),
    m_parents(std::vector<unsigned int>()),
    m_childeren(std::vector<unsigned int>()),
    m_maxLayerNodeId(-1),
    m_minLayerNodeId(100000),
    m_maxLayer(-1),
    m_minLayer(100000),
    m_finished(0),
    m_isOnSectorLimit(false),
    m_headNode(-1),
    m_tailNode(-1),
    m_lastNodeVirtual(false),
    m_lastNodeVirtualId(-1),
    m_headNeigh(std::vector<int>()),
    m_tailNeigh(std::vector<int>()),
    m_toMergeHead(std::vector<unsigned int>()),
    m_toMergeTail(std::vector<unsigned int>())

{} 
// Destructor
PathCandidate::~PathCandidate()
{
  delete m_memberList;
  delete m_memberIdSet;
}

// Copy constructor
PathCandidate::PathCandidate(PathCandidate const &ot)
  : m_id(ot.m_id),
    m_level(ot.m_level),
    m_weight(ot.m_weight),
    m_orientation(ot.m_orientation),
    m_memberList (new std::vector<int>( *(ot.m_memberList) ) ),
    m_memberIdSet(new std::set<int>   ( *(ot.m_memberIdSet) ) ),
    m_OuterToInner(ot.m_OuterToInner),
    m_length(ot.m_length),
    m_isValid(ot.m_isValid),
    m_isShort(ot.m_isShort),
    m_isMerged(ot.m_isMerged),
    m_hasCurvature(ot.m_hasCurvature),
    m_CurV_par(ot.m_CurV_par),
    m_maxlayerNodeId(ot.m_maxlayerNodeId),
    m_minlayerNodeId(ot.m_minlayerNodeId),
    m_maxLayernodeLayer(ot.m_maxLayernodeLayer),
    m_minlayerNodeLayer(ot.m_minlayerNodeLayer),
    m_maxLayerNodeId(ot.m_maxLayerNodeId),
    m_minLayerNodeId(ot.m_minLayerNodeId),
    m_maxLayer(ot.m_maxLayer),
    m_minLayer(ot.m_minLayer),
    m_finished(ot.m_finished),
    m_isOnSectorLimit(ot.m_isOnSectorLimit),
    m_headNode(ot.m_headNode),
    m_tailNode(ot.m_tailNode),
    m_lastNodeVirtual(ot.m_lastNodeVirtual),
    m_lastNodeVirtualId(ot.m_lastNodeVirtualId),
    m_headNeigh(ot.m_headNeigh),
    m_tailNeigh(ot.m_tailNeigh),
    m_toMergeHead(ot.m_toMergeHead),
    m_toMergeTail(ot.m_toMergeTail),
    m_parents(ot.m_parents),
    m_childeren(ot.m_childeren),
    m_eigevValue1(ot.m_eigevValue1),
    m_eigevValue2(ot.m_eigevValue2)
{
  this->m_meanVector[0] = ot.m_meanVector[0];
  this->m_meanVector[1] = ot.m_meanVector[1];

  this->m_covMatrix [0][0] = ot.m_covMatrix [0][0];
  this->m_covMatrix [0][1] = ot.m_covMatrix [0][1];
  this->m_covMatrix [1][0] = ot.m_covMatrix [1][0];
  this->m_covMatrix [1][1] = ot.m_covMatrix [1][1];

  this->m_inveCovMatrix [0][0] = ot.m_inveCovMatrix [0][0];
  this->m_inveCovMatrix [0][1] = ot.m_inveCovMatrix [0][1];
  this->m_inveCovMatrix [1][0] = ot.m_inveCovMatrix [1][0];
  this->m_inveCovMatrix [1][1] = ot.m_inveCovMatrix [1][1];

  this->m_eigenVector1 = ot.m_eigenVector1;
  this->m_eigenVector2 = ot.m_eigenVector2;
}

// Deep copy
PathCandidate& PathCandidate::operator=(PathCandidate const &ot)
{
  if( this != &ot ) {
    this->m_id  = ot.m_id;
    this->m_level = ot.m_level;
    this->m_weight = ot.m_weight;
    this->m_orientation = ot.m_orientation;
    
    // Copy member vector id list (May be deleted)
    (this->m_memberList)->clear();// Remove all
    std::copy((ot.m_memberList)->begin(), (ot.m_memberList)->end(),
     	      std::back_inserter( *(this->m_memberList) ) );
    
    // Clear all members(create empty set)
    (this->m_memberIdSet)->clear();
    // Copy members from the other candidate.
    (this->m_memberIdSet)->insert ( (ot.m_memberIdSet)->begin(), (ot.m_memberIdSet)->end());
    // Copy direction
    this->m_OuterToInner = ot.m_OuterToInner;
    // Path length
    this->m_length = ot.m_length;
    // Curvature measure
    this->m_isValid  = ot.m_isValid;
    this->m_isShort  = ot.m_isShort;
    this->m_isMerged = ot.m_isMerged;
    this->m_hasCurvature = ot.m_hasCurvature;
    this->m_CurV_par = ot.m_CurV_par;
    this->m_maxlayerNodeId = ot.m_maxlayerNodeId;
    this->m_minlayerNodeId = ot.m_minlayerNodeId;
    this->m_maxLayernodeLayer = ot.m_maxLayernodeLayer;
    this->m_minlayerNodeLayer = ot.m_minlayerNodeLayer;

    /*new*/
    this->m_maxLayerNodeId = ot.m_maxLayerNodeId;
    this->m_minLayerNodeId = ot.m_minLayerNodeId;
    this->m_maxLayer = ot.m_maxLayer;
    this->m_minLayer = ot.m_minLayer;
    this->m_lastNodeVirtual = ot.m_lastNodeVirtual;
    this->m_lastNodeVirtualId = ot.m_lastNodeVirtualId;
    this->m_headNeigh = ot.m_headNeigh;
    this->m_tailNeigh = ot.m_tailNeigh;
    this->m_toMergeHead = ot.m_toMergeHead;
    this->m_toMergeTail = ot.m_toMergeTail;
    this->m_finished = ot.m_finished;
    this->m_isOnSectorLimit = ot.m_isOnSectorLimit;
    this->m_headNode = ot.m_headNode ;
    this->m_tailNode = ot.m_tailNode;

    /**/
    this->m_meanVector[0] = ot.m_meanVector[0];
    this->m_meanVector[1] = ot.m_meanVector[1];
    
    this->m_covMatrix [0][0] = ot.m_covMatrix [0][0];
    this->m_covMatrix [0][1] = ot.m_covMatrix [0][1];
    this->m_covMatrix [1][0] = ot.m_covMatrix [1][0];
    this->m_covMatrix [1][1] = ot.m_covMatrix [1][1];
    
    this->m_inveCovMatrix [0][0] = ot.m_inveCovMatrix [0][0];
    this->m_inveCovMatrix [0][1] = ot.m_inveCovMatrix [0][1];
    this->m_inveCovMatrix [1][0] = ot.m_inveCovMatrix [1][0];
    this->m_inveCovMatrix [1][1] = ot.m_inveCovMatrix [1][1];

    this->m_parents = ot.m_parents;
    this->m_childeren = ot.m_childeren;

    this->m_eigenVector1 = ot.m_eigenVector1;
    this->m_eigenVector2 = ot.m_eigenVector2;
  }
  return (*this);
}

void PathCandidate::updateHeadAndTailNodes()
{
  std::set<int>::iterator it = m_memberIdSet->begin();
  m_minlayerNodeId = *it;
  
  it = m_memberIdSet->end();
  it--;
  m_maxlayerNodeId = *it;
  // Skip virtuals
  while( (m_maxlayerNodeId >= START_VIRTUAL_ID) &&
         (it != m_memberIdSet->begin()) ) {
    --it;
    m_maxlayerNodeId = *it;
    //m_maxlayerNodeId = *(--it);
  }
  
}

void PathCandidate::insertNewNode(GridNode *node, int direction)
{
  int id = node->m_detID;
  int layer = node->m_Layer;

  

  (node->m_cm).push_back(m_id);

  if(node->m_type != GridNode::VIRTUAL_NODE){
    m_length++;
    m_maxLayerNodeId = MAX(id, m_maxlayerNodeId);
    m_minLayerNodeId = MIN(id, m_minlayerNodeId);
    m_maxLayer = MAX(layer, m_maxLayer);
    m_minLayer = MIN(layer, m_minLayer);
    m_lastNodeVirtual = false;
    m_memberIdSet->insert(id);
    direction == 1? m_memberList->push_back(id): ( void )m_memberList->insert(m_memberList->begin(),id);

  } else {
    m_lastNodeVirtual = true;
    m_lastNodeVirtualId = id;
  }
	  
  m_headNode = id;
}

