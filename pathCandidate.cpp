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
    m_seenVirtual(false),
    m_lastVirtual(-1),
    m_headNeigh(std::vector<int>()),
    m_tailNeigh(std::vector<int>()),
    m_toMergeHead(std::vector<unsigned int>()),
    m_toMergeTail(std::vector<unsigned int>()),
    m_listSkewed(std::vector<unsigned int>()),
    m_x(std::vector<double>()),
    m_y(std::vector<double>()),
    m_z(std::vector<double>()),
    m_r(std::vector<double>()),
    m_theta(std::vector<double>())

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
    //m_lastNodeVirtual(ot.m_lastNodeVirtual),
    // m_lastNodeVirtualId(ot.m_lastNodeVirtualId),
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
    //   this->m_lastNodeVirtual = ot.m_lastNodeVirtual;
    //  this->m_lastNodeVirtualId = ot.m_lastNodeVirtualId;
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

void PathCandidate::insertNewNode(CoordGrid &gr, GridNode *node,  std::vector<int>::iterator it)
{
  int id = node->m_detID;
  int layer = node->m_Layer;
  std::vector< GridNode > &Ingrid  = gr.m_grid;

  int vecindex = std::distance( m_memberList->begin(), it );
  printf("ADDING %d to track %d \n", id, m_id);

  m_memberIdSet->insert(id);
  m_memberList->insert(it,id);


  (node->m_cm).push_back(m_id);

  if(node->m_type != GridNode::VIRTUAL_NODE){
    
    m_length++;
    m_maxLayerNodeId = MAX(id, m_maxlayerNodeId);
    m_minLayerNodeId = MIN(id, m_minlayerNodeId);
    m_maxLayer = MAX(layer, m_maxLayer);
    m_minLayer = MIN(layer, m_minLayer);
    //m_lastNodeVirtual = false;
    //m_memberIdSet->insert(id);
    //  m_memberList->insert(it,id);
    
    if(node->m_type == GridNode::STT_TYPE_SKEW  && m_seenVirtual )
       m_listSkewed.push_back(id);

    //printf("%d, %d, %d \n", vecindex, m_x.size(), m_memberList->size());

    //  if(node->m_type != GridNode::STT_TYPE_SKEW){
      m_x.insert(m_x.begin() + vecindex,node->m_x);
      m_y.insert(m_y.begin() + vecindex,node->m_y);
      m_z.insert(m_z.begin() + vecindex,node->m_z);
      m_r.insert(m_r.begin() + vecindex,node->m_r);
      m_theta.insert(m_theta.begin() + vecindex,node->m_thetaDeg);

      /* } else {
      m_x.insert(m_x.begin() + vecindex, -1);
      m_y.insert(m_y.begin() + vecindex, -1);
      m_z.insert(m_z.begin() + vecindex, -1);
      m_r.insert(m_r.begin() + vecindex, node->m_r);
      m_theta.insert(m_theta.begin() + vecindex, node->m_thetaDeg);

      if(m_listSkewed.size() >= 3){
	int tt = m_listSkewed.size()-3;
	auto it = find(m_memberList->begin(), 
		       m_memberList->end(), m_listSkewed[tt]);   
	int index = distance(m_memberList->begin(), it);
	int idx = gr.Find(m_listSkewed[tt]);
	GridNode &skewedToproc = Ingrid[idx];
	m_x[index] = skewedToproc.m_xDet;
	m_y[index] = skewedToproc.m_yDet;
	m_z[index] = skewedToproc.m_z_Det;
      }
      }*/
    
  } else {

    if(m_seenVirtual == true && m_listSkewed.size() > 0){
      // Correcting skewed
      int virtIdx = gr.Find(m_lastVirtual);
      GridNode &lastVirtNode = Ingrid[virtIdx];
     
      float x_diff = node->m_x - lastVirtNode.m_x;
      float y_diff = node->m_y - lastVirtNode.m_y;
      float z_diff = node->m_z - lastVirtNode.m_z;

      x_diff /= static_cast<float>(m_listSkewed.size()+1);
      y_diff /= static_cast<float>(m_listSkewed.size()+1);
      z_diff /= static_cast<float>(m_listSkewed.size()+1);

      float xInc = lastVirtNode.m_x + x_diff;
      float yInc = lastVirtNode.m_y + y_diff;
      float zInc = lastVirtNode.m_z + z_diff;

      /* Correct xy-coordinates of the skewed nodes */
      for(size_t m = 0; m < m_listSkewed.size(); ++m) {
	int idx = gr.Find(m_listSkewed[m]);
	GridNode &skewedToproc = Ingrid[idx];
	if(layer != lastVirtNode.m_Layer){
	  skewedToproc.m_xDet = xInc;
	  skewedToproc.m_yDet = yInc;
	}
	skewedToproc.m_z_Det = zInc;
	auto it = find(m_memberList->begin(), 
		       m_memberList->end(), m_listSkewed[m]); 
  
	int index = distance(m_memberList->begin(), it); 
    
	m_x[index] = skewedToproc.m_xDet;
	m_y[index] = skewedToproc.m_yDet;
	m_z[index] = skewedToproc.m_z_Det;
	
	std::pair<float, float> r_Theta;
	float theta_deg = Cartesian_To_Polar(skewedToproc.m_xDet, skewedToproc.m_yDet, r_Theta);
	m_r[index] = r_Theta.first;
	m_theta[index] = theta_deg;
	
	printf("New node %d with %lf and %lf, (r %lf, theta %lf) \n",m_listSkewed[m], skewedToproc.m_xDet ,  skewedToproc.m_yDet,m_r[index], theta_deg ); 
	xInc += x_diff;
	yInc += y_diff;

      }
      
      m_listSkewed.clear();
    }
    // printf("%d, %d, %d, %lf, %lf \n", vecindex, m_x.size(), m_memberList->size(),node->m_xDet,node->m_yDet);
    m_x.insert(m_x.begin() + vecindex, node->m_xDet);
    m_y.insert(m_y.begin() + vecindex, node->m_yDet);
    m_z.insert(m_z.begin() + vecindex, node->m_z_Det);

    /*  std::pair<float, float> r_Theta;
    float theta_deg = Cartesian_To_Polar(node->m_xDet, node->m_yDet, r_Theta);
    node->m_r = r_Theta.first;
    node->m_thetaDeg = theta_deg;*/
	
    m_r.insert(m_r.begin() + vecindex,node->m_r);
    m_theta.insert(m_theta.begin() + vecindex,node->m_thetaDeg);
    
    //  printf("Added virtual with r %lf, and theta %lf \n", node->m_r,node->m_thetaDeg);
    m_seenVirtual = true;
    m_lastVirtual = node->m_detID;
  }
	  
  m_headNode = id;
  // printf("END ADDING \n");
}

