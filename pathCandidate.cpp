/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#include <iterator>
#include "pathCandidate.h"
#include "auxfunctions.h"
#include "logc.h"
// Constructors
PathCandidate::PathCandidate()
  : m_id(0),
    m_level(0),
    m_weight(0),
    m_scattAngle(-1),
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
    m_nPerLayer(0),
    m_finished(0),
    m_isOnSectorLimit(false),
    m_headNode(-1),
    m_tailNode(-1),
    m_numVirtual(0),
    m_seenVirtual(false),
    m_lastVirtual(-1),
    m_lastNodeAdded(-1),
    m_headNeigh(std::vector<int>()),
    m_tailNeigh(std::vector<int>()),
    m_toMergeHead(std::vector<unsigned int>()),
    m_toMergeTail(std::vector<unsigned int>()),
    m_listSkewed(std::vector<unsigned int>()),
    m_sectors(std::vector<unsigned int>()),
    m_layers(std::vector< int>()),
    m_anchors(std::vector< GridNode>()),
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
    m_scattAngle(ot.m_scattAngle),

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
    m_numVirtual(ot.m_numVirtual),

    //m_lastNodeVirtual(ot.m_lastNodeVirtual),
    // m_lastNodeVirtualId(ot.m_lastNodeVirtualId),
    m_headNeigh(ot.m_headNeigh),
    m_tailNeigh(ot.m_tailNeigh),
    m_toMergeHead(ot.m_toMergeHead),
    m_toMergeTail(ot.m_toMergeTail),
    m_sectors(ot.m_sectors),
    m_layers(ot.m_layers),

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
    this->m_scattAngle = ot.m_scattAngle;

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
    this->m_numVirtual = ot.m_numVirtual;

    this->m_sectors = ot.m_sectors;
    this->m_layers = ot.m_layers;

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



void PathCandidate::fEstimateSkewedXY(CoordGrid &gr, std::vector< GridNode > &Ingrid, GridNode *node, size_t vecindex){

  // To redetermine the xy (and z) of each node in between two virtual layers
 
 int dir = vecindex == m_memberList->size()-1 ? -1 : 1;
 int off = vecindex == m_memberList->size()-1 ? m_memberList->size()-1: 0;
 int offAnc = vecindex == m_memberList->size()-1 ? m_anchors.size()-1: 0;

 // Determine the new anchor based on the virtual nodes found (max 3, otherwise it is messy)
 // error("size of prev %d", m_prevVirtuals.size());
  size_t nVirt = MIN(m_prevVirtuals.size(), 3);
  GridNode virt;
  if(nVirt > 0){
    virt =  Ingrid[m_prevVirtuals[0]-1];
    for(size_t i =1; i< nVirt; i++){
      virt.m_xDet += Ingrid[m_prevVirtuals[i]-1].m_xDet;
      virt.m_yDet += Ingrid[m_prevVirtuals[i]-1].m_yDet;
      virt.m_zDet += Ingrid[m_prevVirtuals[i]-1].m_zDet;
    }
    virt.m_xDet  /= nVirt;
    virt.m_yDet  /= nVirt;
    virt.m_zDet  /= nVirt;
    //   dbgconnect("The virt anchor nodes is %f, %f, %f",  virt.m_xDet, virt.m_yDet, virt.m_z_Det) ;

    
    //  error("size of prev %d", m_prevVirtuals.size());

    // Determine the last node that xy was previously re-determined 
    int  curLayer = -1, nNodes=0; 
    GridNode secondAnc;
    for(size_t i = m_prevVirtuals.size()+1, j = 0; i < m_memberList->size(); i++){
      GridNode mynode = Ingrid[m_memberList->at(off+dir*i)-1];
      //dbgtrkerror("%d, %d, %d", j, mynode.m_detID, mynode.m_weight);
      //  continue;
      if(curLayer == -1){
	// dbgconnect("First node found is %d (i %d, j %d)",mynode.m_detID, i ,j);
	secondAnc = mynode;
	curLayer = mynode.m_Layer;
	nNodes++;
	//   break;
      } else if(mynode.m_Layer == curLayer){
	//dbgconnect("An other node found is %d",mynode.m_detID);
	secondAnc.m_xDet += mynode.m_xDet;
	secondAnc.m_yDet += mynode.m_yDet;
	secondAnc.m_zDet += mynode.m_zDet;
	nNodes++;
      } else if( curLayer != -1 && (mynode.m_Layer - curLayer != 0)){
	break;
      }
    }
    // error("size of prev %d", m_prevVirtuals.size());

    secondAnc.m_xDet /= nNodes;
    secondAnc.m_yDet /= nNodes;
    secondAnc.m_zDet /= nNodes;

    PointsLineIntersectFinal( *node, secondAnc.m_xDet, virt.m_xDet, secondAnc.m_yDet, virt.m_yDet);
    //   error("The new node pos is %f, %f, %f",  node->m_xDet, node->m_yDet, node->m_z_Det) ;

  }
}



void PathCandidate::determineSkewedXY(CoordGrid &gr, std::vector< GridNode > &Ingrid, GridNode *node, size_t vecindex){

  // To redetermine the xy (and z) of each node in between two virtual layers
 
 int dir = vecindex == m_memberList->size()-1 ? -1 : 1;
 int off = vecindex == m_memberList->size()-1 ? m_memberList->size()-1: 0;
 int offAnc = vecindex == m_memberList->size()-1 ? m_anchors.size()-1: 0;
	
 // Because virtual between Para and Ske are not very well determined, we skip them
 if(m_anchors[offAnc].m_type  != GridNode::STT_TYPE_PARA){
   // Determine the new anchor based on the virtual nodes found (max 3, otherwise it is messy)
   GridNode virt =  Ingrid[m_prevVirtuals[0]-1];
   size_t nVirt = MIN(m_prevVirtuals.size(), 3);
   for(size_t i =1; i< nVirt; i++){
     virt.m_xDet += Ingrid[m_prevVirtuals[i]-1].m_xDet;
     virt.m_yDet += Ingrid[m_prevVirtuals[i]-1].m_yDet;
     virt.m_zDet += Ingrid[m_prevVirtuals[i]-1].m_zDet;
   }
   virt.m_xDet  /= nVirt;
   virt.m_yDet  /= nVirt;
   virt.m_zDet  /= nVirt;
   //  info("The virt anchor nodes is %f, %f, %f",  virt.m_xDet, virt.m_yDet, virt.m_zDet) ;

   // Determine the last node that xy was previously re-determined 
   int lastAnc =0, incAnc = 0, curLayer = -1, nNodes = 0, nToCorrect = 0;; 
   GridNode secondAnc;
   for(size_t i = m_prevVirtuals.size()+1, j = 0; i < m_memberList->size(); i++){
     GridNode mynode = Ingrid[m_memberList->at(off+dir*i)-1];
     //dbgtrkerror("%d, %d, %d", j, mynode.m_detID, mynode.m_weight);

     if(mynode.m_type != GridNode::VIRTUAL_NODE && !mynode.m_weight)
       nToCorrect++;
     //  continue;
     if(mynode.m_weight == 1 && curLayer == -1){
       lastAnc = i;
       incAnc = j;
       // dbgconnect("First node found is %d (i %d, j %d)",mynode.m_detID, i ,j);
       secondAnc = mynode;
       nNodes++;
       curLayer = mynode.m_Layer;
       //   break;
     } else if(mynode.m_weight == 1 && mynode.m_Layer == curLayer){
       //dbgconnect("An other node found is %d",mynode.m_detID);
       secondAnc.m_xDet += mynode.m_xDet;
       secondAnc.m_yDet += mynode.m_yDet;
       secondAnc.m_zDet += mynode.m_zDet;
       nNodes++;
     } else if( curLayer != -1 && (mynode.m_Layer - curLayer != 0)){
       break;
     }
     j++;
   }

   secondAnc.m_xDet /= nNodes;
   secondAnc.m_yDet /= nNodes;
   secondAnc.m_zDet /= nNodes;
   // info("The previous anchor nodes %d is %f, %f, %f", secondAnc.m_detID,  secondAnc.m_xDet, secondAnc.m_yDet, secondAnc.m_zDet) ;
   // info("nToCorrect %d", nToCorrect);
   if(labs(secondAnc.m_Layer - virt.m_Layer < 2) || nToCorrect > 6){
     //  error("%d %d %d",secondAnc.m_Layer,  virt.m_Layer, labs(secondAnc.m_Layer - virt.m_Layer));
     for(size_t i =0; i< 1; i++){
       GridNode &curVirt =  Ingrid[m_prevVirtuals[i]-1];
       curVirt.m_weight = 1;
       auto it = find(m_memberList->begin(), m_memberList->end(), curVirt.m_detID);
       int index = it - m_memberList->begin();
	  //dbgconnect("i %d, ECorrecting %d, %d, %d", i, inter.m_detID, prev == NULL ? 0:prev->m_Layer, inter.m_Layer);
       
       m_x[index] = virt.m_xDet;
       m_y[index] = virt.m_yDet;
       m_z[index] = virt.m_zDet;
     }

     /* for(size_t i =0; i< 1; i++){
       GridNode &curVirt =  Ingrid[gr.Find(m_prevVirtuals[i])];
       //   error("Virtual %d weights %d", curVirt.m_detID,curVirt.m_weight);

       }*/

      for(int i = 0, j = 0; j < incAnc; i++){
	GridNode &inter =  Ingrid[m_memberList->at(off+dir*lastAnc-dir*(i+1))-1];

	//  dbgconnect("i %d, j%d, Correcting %d, %d, %d", i,j, inter.m_detID, prev == NULL ? 0:prev->m_Layer, inter.m_Layer);
	// if(inter.m_type == GridNode::VIRTUAL_NODE)
	//	continue;

	if(inter.m_type != GridNode::VIRTUAL_NODE){
	  auto it = find(m_memberList->begin(), m_memberList->end(), inter.m_detID);
	  int index = it - m_memberList->begin();
	  //dbgconnect("i %d, ECorrecting %d, %d, %d", i, inter.m_detID, prev == NULL ? 0:prev->m_Layer, inter.m_Layer);

	  m_x[index] = inter.m_x;
	  m_y[index] = inter.m_y;
	  m_z[index] = 0.;
	  //	dbgconnect("i %d, Corrected %d, %f, %f", i, inter.m_detID,m_x[index], m_y[index]);

	  inter.m_weight = 0;
	  //nNodeperL++;

	}
	  
	j++;
	//prev = &inter;
      }   
     
     return;
   }
   //dbgconnect("The previous anchor nodes %d is %f, %f, %f", secondAnc.m_detID,  secondAnc.m_xDet, secondAnc.m_yDet, secondAnc.m_z_Det) ;
  std::vector<int> layers;

   // Let's correct all nodes found between the first and last node anchors
    GridNode *prev = NULL;
    int nNodeperL = 0;
    for(int i = 0, j = 0; j < incAnc; i++){
      GridNode &inter =  Ingrid[m_memberList->at(off+dir*lastAnc-dir*(i+1))-1];

      //info("i %d, j%d, Correcting %d, %d, %d", i,j, inter.m_detID, prev == NULL ? 0:prev->m_Layer, inter.m_Layer);
      // if(inter.m_type == GridNode::VIRTUAL_NODE)
      //	continue;

      if(inter.m_type != GridNode::VIRTUAL_NODE){
	if(prev == NULL || prev->m_Layer != inter.m_Layer || nNodeperL < 3){
	  if(prev != NULL && prev->m_Layer != inter.m_Layer)
	    nNodeperL = 0;
	  PointsLineIntersectFinal( inter, secondAnc.m_xDet, virt.m_xDet, secondAnc.m_yDet, virt.m_yDet);
	  auto it = find(m_memberList->begin(), m_memberList->end(), inter.m_detID);
	  int index = it - m_memberList->begin();
	  //dbgconnect("i %d, ECorrecting %d, %d, %d", i, inter.m_detID, prev == NULL ? 0:prev->m_Layer, inter.m_Layer);
	  //  info("i %d, Bef %d, %f, %f", i, inter.m_detID,m_x[index], m_y[index]);

	  m_x[index] = inter.m_xDet;
	  m_y[index] = inter.m_yDet;
	  m_z[index] = inter.m_zDet;
	  // info("i %d, Corrected %d, %f, %f", i, inter.m_detID,m_x[index], m_y[index]);

	  inter.m_weight = 1;
	  nNodeperL++;

	}
      }
      j++;
      prev = &inter;
    }
    //info("HERE");
 } else {

   GridNode virt =  Ingrid[m_prevVirtuals[0]-1];
   size_t nVirt = MIN(m_prevVirtuals.size(), 3);
   for(size_t i =1; i< nVirt; i++){
     virt.m_xDet += Ingrid[m_prevVirtuals[i]-1].m_xDet;
     virt.m_yDet += Ingrid[m_prevVirtuals[i]-1].m_yDet;
     virt.m_zDet += Ingrid[m_prevVirtuals[i]-1].m_zDet;
   }
   virt.m_xDet  /= nVirt;
   virt.m_yDet  /= nVirt;
   virt.m_zDet  /= nVirt;
   for(size_t i =0; i< m_prevVirtuals.size(); i++){
     GridNode &curVirt =  Ingrid[m_prevVirtuals[i]-1];
     curVirt.m_weight = 1;
     //auto it = find(m_memberList->begin(), m_memberList->end(), curVirt.m_detID);
     // int index = it - m_memberList->begin();
     //dbgconnect("i %d, ECorrecting %d, %d, %d", i, inter.m_detID, prev == NULL ? 0:prev->m_Layer, inter.m_Layer);
     
     // m_x[index] = virt.m_xDet;
     // m_y[index] = virt.m_yDet;
     //m_z[index] = virt.m_z_Det;
   }
 }
}

void PathCandidate::addToAnchor(CoordGrid &gr, std::vector< GridNode > &Ingrid, GridNode *node, size_t vecindex, int prevLayer){
  
  int layer = node->m_Layer;
  int dir = vecindex == m_memberList->size()-1 ? -1 : 1;
  int off = vecindex == m_memberList->size()-1 ? m_anchors.size()-1: 0;
  // dbgconnect("node %d with  %f %f %f", node->m_detID, node->m_xDet, node->m_yDet, node->m_z_Det);
  //If we have a new node on an other layer, let's make a new anchor
  if(prevLayer != layer){
    if(node->m_type == GridNode::STT_TYPE_SKEW && m_length > 2){
      PointsLineIntersectLive( *node, m_anchors[off+3*dir].m_xDet, m_anchors[off].m_xDet, //*node
			       m_anchors[off+3*dir].m_yDet, m_anchors[off].m_yDet); //Output
      // dbgconnect("New anchor correct, node %d with  %f %f %f", node->m_detID, node->m_xDet, node->m_yDet, node->m_z_Det);

    }
    m_anchors.push_back(*node);
    m_nPerLayer = 0; // To be corrected for tail and head node
  }
  else {
    // The idea here is to group nodes two by two on the same layer to reduce some scattering in the position    
    if(!(m_nPerLayer%2)){
      m_anchors.push_back(*node);
      //dbgconnect("New anchor on same layer is %f %f %f", node->m_xDet, node->m_yDet, node->m_z_Det);
    } else {
      m_anchors[off].m_xDet =  (m_anchors[off].m_xDet + node->m_xDet)/2.;
      m_anchors[off].m_yDet =  (m_anchors[off].m_yDet + node->m_yDet)/2.;
      m_anchors[off].m_zDet =  (m_anchors[off].m_zDet + node->m_zDet)/2.;
      //dbgconnect("We corrected previous anchor as %f %f %f", m_anchors[off].m_xDet, m_anchors[off].m_yDet, m_anchors[off].m_z_Det);
    }
  }
}

void PathCandidate::correctPrevAnchor(CoordGrid &gr, std::vector< GridNode > &Ingrid,GridNode *node, size_t vecindex){

  int dir = vecindex == m_memberList->size()-1 ? -1 : 1;
  int off = vecindex == m_memberList->size()-1 ? m_anchors.size()-1: 0;

  if(m_anchors[off].m_type  != GridNode::STT_TYPE_PARA){
    
    //   dbgconnect("We have %d virtual nodes", m_prevVirtuals.size());
    GridNode virt =  Ingrid[m_prevVirtuals[0]-1];
    size_t nVirt = MIN(m_prevVirtuals.size(), 3);
    for(size_t i =1; i< nVirt; i++){
      virt.m_xDet += Ingrid[m_prevVirtuals[i]-1].m_xDet;
      virt.m_yDet += Ingrid[m_prevVirtuals[i]-1].m_yDet;
      virt.m_zDet +=Ingrid[m_prevVirtuals[i]-1].m_zDet;
    }
    virt.m_xDet  /= nVirt;
    virt.m_yDet  /= nVirt;
    virt.m_zDet /= nVirt;
    
    //  dbgconnect("The virtual anchor is %f, %f, %f",     virt.m_xDet, virt.m_yDet, virt.m_z_Det) ;

    // Find last correct anchor
    int lastAnc = -1;
    for(size_t i = 1; i < m_anchors.size(); i++){
      if(m_anchors[off+dir*i].m_weight == 1){
	lastAnc = i;
	break;
      }
    }
     
    GridNode &secondAnc = m_anchors[off+dir*lastAnc];   
    // dbgconnect("The previous anchor is %f, %f, %f",     secondAnc.m_xDet, secondAnc.m_yDet, secondAnc.m_z_Det) ;

    if(labs(secondAnc.m_Layer - virt.m_Layer < 2)){
      //   dbgconnect("Both these anchors are almost on the same layer, let's be careful and not correct anything");
      virt.m_weight = 1;
      m_anchors.push_back(virt);
      return;
    }
    
    // Correct all anchors in between
    GridNode *prev = NULL;
    std::vector<int> layers;
    for(int i = 0; i < lastAnc; i++){
      GridNode &inter = m_anchors[off+dir*lastAnc-dir*(i+1)];
      //  dbgconnect("Inter is %f, %f, %f, layer %d",     inter.m_xDet, inter.m_yDet, inter.m_z_Det, inter.m_Layer) ;
      if(prev == NULL  ||std::find(layers.begin(), layers.end(), inter.m_Layer)==layers.end()){
	PointsLineIntersectFinal( inter, secondAnc.m_xDet, virt.m_xDet, secondAnc.m_yDet, virt.m_yDet);
	inter.m_weight = 1;
	//	dbgconnect("Corrected Inter is %f, %f, %f",     inter.m_xDet, inter.m_yDet, inter.m_z_Det, inter.m_Layer) ;
	layers.push_back(inter.m_Layer);
      }
      prev = &inter;
    }

    virt.m_weight = 1;
    m_anchors.push_back(virt);
  }

}

void PathCandidate::insertNewNode(CoordGrid &gr, std::vector< GridNode > &Ingrid, GridNode *node,  std::vector<int>::iterator it)
{
  
  int id = node->m_detID;

  int layer = node->m_Layer, prevLayer;
  //  float xadd = node->m_xDet, yadd = node->m_yDet, zadd = node->m_z_Det;
  size_t vecindex = (size_t) std::distance( m_memberList->begin(), it );
  node->m_weight = node->m_type == GridNode::STT_TYPE_PARA? 1: 0;

  //info("Inserting node %d at %lu", id, vecindex);

  m_memberIdSet->insert(id);
  m_memberList->insert(it,id);
  (node->m_cm).push_back(m_id);

  if(vecindex != m_memberList->size() - 1 && vecindex != 0){
    //dbgconnect("We inserted this node in the middle, so let's not do anything else");
    m_x.insert(m_x.begin() + vecindex, node->m_xDet);
    m_y.insert(m_y.begin() + vecindex, node->m_yDet);
    m_z.insert(m_z.begin() + vecindex, node->m_zDet);     
    m_r.insert(m_r.begin() + vecindex, node->m_r);    
    m_theta.insert(m_theta.begin() + vecindex,node->m_thetaDeg);
    return;
  }
  
  if (vecindex == 0){
    m_tailNode = id;
    prevLayer = m_length > 2? m_layers[0]: -1;
  }
  else if(vecindex == m_memberList->size() - 1){
    m_headNode = id;
    prevLayer = m_length > 2? m_layers.back(): -1;
  }
   

  
  if(node->m_type != GridNode::VIRTUAL_NODE){
    if(layer == prevLayer)
      m_nPerLayer++;

    // Correcting the intersection points on the skewed tubes
    if(m_lastNodeAdded >= gr.firstVirtIdx){ // We might have an issue when we stop on a skewed layer
      //info("Correcting");
      //    info("Correcting");

      //dbgconnect("We need to correct the previous anchors");
      correctPrevAnchor(gr, Ingrid, node, vecindex);
      //  dbgconnect("We need to correct the previous nodes");

      determineSkewedXY(gr, Ingrid, node, vecindex);
      m_listSkewed.clear();
      m_prevVirtuals.clear();
    }
    // dbgconnect("Adding a new anchor");

    // Adding a new anchor, temporary for the skewed nodes
    addToAnchor(gr, Ingrid, node, vecindex, prevLayer);
    // dbgconnect("Adding a new anchor");

    if(vecindex == 0)
      m_layers.insert(m_layers.begin(),layer);
    else
      m_layers.push_back(layer);

    m_sectors.push_back(node->m_Sector);

    m_length++;
	      //dbgconnect("Adding a new anchor");

    // THis is never used, we need to clean the PathCandidate structure
    m_maxLayerNodeId  = MAX(id, m_maxlayerNodeId); 
    m_minLayerNodeId  = MIN(id, m_minlayerNodeId);
    m_maxLayer 	= MAX(layer, m_maxLayer);
    m_minLayer 	= MIN(layer, m_minLayer);

    if(node->m_type == GridNode::STT_TYPE_SKEW &&m_seenVirtual ) // pushing for later correction
      m_listSkewed.push_back(id);    


    
  } else {

    m_seenVirtual = true;
    m_numVirtual++;
    m_lastVirtual = node->m_detID;
    m_prevVirtuals.push_back(id);
  }
    /* if(m_seenVirtual == true && m_listSkewed.size() > 0){
      dbgconnect("We have %d virtuals and %d skewed to correct\n", m_prevVirtuals.size(), m_listSkewed.size());
      determineSkewedXY(gr, Ingrid, node);
      m_listSkewed.clear();
      m_prevVirtuals.clear();
      }*/
    // printf("%d, %d, %d, %lf, %lf \n", vecindex, m_x.size(), m_memberList->size(),node->m_xDet,node->m_yDet);

  double newtheta = node->m_thetaDeg;
  if(m_length > 2){
    double prevtheta = m_theta.back();
    if(prevtheta > 150 && newtheta  < -150)
      newtheta += 360.;
    else if (prevtheta < -150 && newtheta > 150)
      newtheta -= 360.;
  }
  m_x.insert(m_x.begin() + vecindex, node->m_xDet);
  m_y.insert(m_y.begin() + vecindex, node->m_yDet);
  m_z.insert(m_z.begin() + vecindex, node->m_zDet);     
  m_r.insert(m_r.begin() + vecindex, node->m_r);    
  m_theta.insert(m_theta.begin() + vecindex,newtheta);
    
  m_lastNodeAdded = id;
  return;
}



void PathCandidate::insertNewNodeFinal(CoordGrid &gr, std::vector< GridNode > &Ingrid, GridNode *node,  std::vector<int>::iterator it)
{
  
  int id = node->m_detID;
  //error("Inserting node %d ", id);

  int layer = node->m_Layer, prevLayer;
  //  float xadd = node->m_xDet, yadd = node->m_yDet, zadd = node->m_z_Det;
  size_t vecindex = (size_t) std::distance( m_memberList->begin(), it );
  node->m_weight = node->m_type == GridNode::STT_TYPE_PARA? 1: 0;

  // info("Inserting node %d at %lu", id, vecindex);

  m_memberIdSet->insert(id);
  m_memberList->insert(it,id);
  (node->m_cm).push_back(m_id);

  if(vecindex != m_memberList->size() - 1 && vecindex != 0){
    //dbgconnect("We inserted this node in the middle, so let's not do anything else");
    m_x.insert(m_x.begin() + vecindex, node->m_xDet);
    m_y.insert(m_y.begin() + vecindex, node->m_yDet);
    m_z.insert(m_z.begin() + vecindex, node->m_zDet);     
    m_r.insert(m_r.begin() + vecindex, node->m_r);    
    m_theta.insert(m_theta.begin() + vecindex,node->m_thetaDeg);
    return;
  }
  
  if (vecindex == 0){
    m_tailNode = id;
    prevLayer = m_length > 2? m_layers[0]: -1;
  }
  else if(vecindex == m_memberList->size() - 1){
    m_headNode = id;
    prevLayer = m_length > 2? m_layers.back(): -1;
  }
   

  
  if(node->m_type != GridNode::VIRTUAL_NODE){
    if(layer == prevLayer)
      m_nPerLayer++;

    //  error("FAIL HERE");
    if (node->m_type == GridNode::STT_TYPE_SKEW){
      //  error("FAIL HERE");
      fEstimateSkewedXY(gr, Ingrid, node, vecindex);
    }
    //   error("FAIL HERE");
    //    error("The new node pos is %f, %f, %f",  node->m_xDet, node->m_yDet, node->m_z_Det) ;

    // Correcting the intersection points on the skewed tubes
    if(m_lastNodeAdded >= gr.firstVirtIdx){ // We might have an issue when we stop on a skewed layer
      //dbgconnect("We need to correct the previous anchors");
      // error("%d", m_lastNodeAdded);
      correctPrevAnchor(gr, Ingrid, node, vecindex);
      //  dbgconnect("We need to correct the previous nodes");

      determineSkewedXY(gr, Ingrid, node, vecindex);
      m_listSkewed.clear();
      m_prevVirtuals.clear();
    }



    // dbgconnect("Adding a new anchor");

    // Adding a new anchor, temporary for the skewed nodes
    addToAnchor(gr, Ingrid, node, vecindex, prevLayer);
    // dbgconnect("Adding a new anchor");
    //    error("The new node pos is %f, %f, %f",  node->m_xDet, node->m_yDet, node->m_z_Det) ;


    if(vecindex == 0)
      m_layers.insert(m_layers.begin(),layer);
    else
      m_layers.push_back(layer);

    m_sectors.push_back(node->m_Sector);

    m_length++;
	      //dbgconnect("Adding a new anchor");

    // THis is never used, we need to clean the PathCandidate structure
    m_maxLayerNodeId  = MAX(id, m_maxlayerNodeId); 
    m_minLayerNodeId  = MIN(id, m_minlayerNodeId);
    m_maxLayer 	= MAX(layer, m_maxLayer);
    m_minLayer 	= MIN(layer, m_minLayer);

    if(node->m_type == GridNode::STT_TYPE_SKEW &&m_seenVirtual ) // pushing for later correction
      m_listSkewed.push_back(id);    


    
  } else {
    // if(m_lastNodeAdded < START_VIRTUAL_ID)
    //    m_prevVirtuals.clear();
  
    m_seenVirtual = true;
    m_numVirtual++;
    m_lastVirtual = node->m_detID;
    m_prevVirtuals.push_back(id);
  }
    /* if(m_seenVirtual == true && m_listSkewed.size() > 0){
      dbgconnect("We have %d virtuals and %d skewed to correct\n", m_prevVirtuals.size(), m_listSkewed.size());
      determineSkewedXY(gr, Ingrid, node);
      m_listSkewed.clear();
      m_prevVirtuals.clear();
      }*/
    // printf("%d, %d, %d, %lf, %lf \n", vecindex, m_x.size(), m_memberList->size(),node->m_xDet,node->m_yDet);

  double newtheta = node->m_thetaDeg;
  if(m_length > 2){
    double prevtheta = m_theta.back();
    if(prevtheta > 150 && newtheta  < -150)
      newtheta += 360.;
    else if (prevtheta < -150 && newtheta > 150)
      newtheta -= 360.;
  }
  //  error("The new node pos is %f, %f, %f",  node->m_xDet, node->m_yDet, node->m_z_Det) ;

  m_x.insert(m_x.begin() + vecindex, node->m_xDet);
  m_y.insert(m_y.begin() + vecindex, node->m_yDet);
  m_z.insert(m_z.begin() + vecindex, node->m_zDet);     
  m_r.insert(m_r.begin() + vecindex, node->m_r);    
  m_theta.insert(m_theta.begin() + vecindex,newtheta);
    
  m_lastNodeAdded = id;
  //info("Inserted node %d at %lu", id, vecindex);

  return;
}



/*
void PathCandidate::insertNewNode(CoordGrid &gr, std::vector< GridNode > &Ingrid, GridNode *node,  std::vector<int>::iterator it)
{
  
  int id = node->m_detID;
  int layer = node->m_Layer;
  //  std::vector< GridNode > &Ingrid  = gr.m_grid;
  float xadd = node->m_xDet, yadd = node->m_yDet, zadd = node->m_z_Det;
  int vecindex = std::distance( m_memberList->begin(), it );
  //  printf("ADDING %d to track %d, layer %d \n", id, m_id, layer);
  //dbgconnect("Addind %d", neighNode->m_detID);

  m_memberIdSet->insert(id);
  m_memberList->insert(it,id);
  m_layers.insert(m_layers.begin()+vecindex,layer);

  if (std::find(m_sectors.begin(), m_sectors.end(),node->m_Sector)==m_sectors.end())
    m_sectors.push_back(node->m_Sector);

  (node->m_cm).push_back(m_id);

  if(node->m_type != GridNode::VIRTUAL_NODE){
    
    m_length++;
    m_maxLayerNodeId = MAX(id, m_maxlayerNodeId);
    m_minLayerNodeId = MIN(id, m_minlayerNodeId);
    m_maxLayer = MAX(layer, m_maxLayer);
    m_minLayer = MIN(layer, m_minLayer);
    
    if(node->m_type == GridNode::STT_TYPE_SKEW &&m_seenVirtual ){
      m_listSkewed.push_back(id);
      //Let's try to find the closesrt previous tube on previous layer
      /* std::vector<GridNode *> 	 virt;
      std::vector<GridNode *> 	 prevNodes;     

      if(vecindex != 0){ //Adding in head direction
	for (int i = m_memberList->size()-2;  i>=0; i--){
	  int prevId = m_memberList->at(i);
	  int prevIdx = gr.Find(prevId);
	  GridNode &prevNode = Ingrid[prevIdx];
	  //  printf("Id %d, layer %u \n", prevId, prevNode.m_Layer);

	  if (prevNode.m_type == GridNode::VIRTUAL_NODE) {
	    virt.push_back(&prevNode);
	    //    printf("Found this virtual %d \n", prevId);
	  }
	  else if ((prevNode.m_type == GridNode::STT_TYPE_SKEW || prevNode.m_type == GridNode::STT_TYPE_PARA) && prevNode.m_Layer != layer){
	    if(labs(prevNode.m_Layer - (size_t) layer) < 2){
	      //   printf("Found this previous node %d\n", prevId);
	      prevNodes.push_back(&prevNode);
	    } else{
	      //  printf("Diff in layer is too high %d to %d\n", layer, prevNode.m_Layer);
	      break;
	    }
	  }
	}
	if(prevNodes.size() > 0 && virt.size()>0){
	  // printf("Let's find a previous node\n");
	  float mindist = 100;
	  int goodid = -1;
	  GridNode *goodNode;// = prevNodes[0];
	  for(size_t i = 0; i < prevNodes.size(); i++){
	    float curdist = sqrt(pow(node->m_x-prevNodes[i]->m_xDet,2) + pow(node->m_y-prevNodes[i]->m_yDet,2));
	    //	    printf("Node %d curd dist %f\n", prevNodes[i]->m_detID, curdist);

	    if(curdist<mindist){
	      goodid = prevNodes[i]->m_detID;
	      goodNode = prevNodes[i];
	      //	      printf("Found %d\n", goodid);
	    }
	    // goodNode.m_x += Ingrid[gr.Find(m_prevVirtuals[p])].m_x;
	    //  //  goodNode.m_y += Ingrid[gr.Find(m_prevVirtuals[p])].m_y;
	    //  .m_z += Ingrid[gr.Find(m_prevVirtuals[p])].m_z;
	  }
	  //	  anchorPrev.m_x /= (float) m_prevVirtuals.size();
	  //  anchorPrev.m_y /= (float)m_prevVirtuals.size();
	  //  anchorPrev.m_z /= (float) m_prevVirtuals.size();
	
	if(goodid > -1){
	    for(size_t i = 0; i < virt.size(); i++){
	      //   printf("%d, %d\n", virt[i]->m_neighbors[0], virt[i]->m_neighbors[1]);
	      if(virt[i]->IsNeighboring(id) && virt[i]->IsNeighboring(goodid)){
		//	printf("We use virt %d to modifiy this node %d", virt[i]->m_detID, id);
		PointsLineIntersect( *node, goodNode->m_xDet, virt[i]->m_xDet,
					  goodNode->m_yDet, virt[i]->m_yDet); //Output 

		xadd = node->m_xDet;
		yadd = node->m_yDet;
		zadd = node->m_z_Det;
		break;
	      }
	    }
	  }
	  }*/
      /*
      	  for(size_t i = 0; i < prevNodes.size(); i++){
	    float curdist = sqrt(pow(node->m_xDet-prevNodes[i]->m_xDet,2) + pow(node->m_yDet-prevNodes[i]->m_yDet,2));
	    //	    printf("Node %d curd dist %f\n", prevNodes[i]->m_detID, curdist);

	    if(curdist<mindist){
	      goodid = prevNodes[i]->m_detID;
	      goodNode = prevNodes[i];
	      //	      printf("Found %d\n", goodid);
	      }
	    //   //  goodNode.m_x += prevNodes[i]->m_x;
	    //   goodNode.m_y += prevNodes[i]->m_y;
	    //    goodNode.m_z += prevNodes[i]->m_z;
	  }
	  //  goodNode.m_x /= (float) prevNodes.size();
	  //  goodNode.m_y /= (float) prevNodes.size();
	  // goodNode.m_z /= (float) prevNodes.size();

	  GridNode virtNode = *virt[0];
	  int nvirt = 1;
	  // if(goodid > -1){
	  for(size_t i = 0; i < virt.size(); i++){
	    //   printf("%d, %d\n", virt[i]->m_neighbors[0], virt[i]->m_neighbors[1]);
	    if(virt[i]->IsNeighboring(id)){// && virt[i]->IsNeighboring(goodid)){
	      //	printf("We use virt %d to modifiy this node %d", virt[i]->m_detID, id);
	      virtNode.m_x += virt[i]->m_x;
	      virtNode.m_y += virt[i]->m_y;
	      virtNode.m_z += virt[i]->m_z;
	      nvirt++;
	    }
	  }
	  virtNode.m_x /= (float) nvirt;
	  virtNode.m_y /= (float) nvirt;
	  virtNode.m_z /= (float) nvirt;
	  //	PointsLineIntersect( *node, goodNode->m_xDet, virt[i]->m_xDet,
	  //				  goodNode->m_yDet, virt[i]->m_yDet); //Output 
	  //  PointsLineIntersect( *node, goodNode.m_xDet, virtNode.m_xDet,
	  //      goodNode.m_yDet, virtNode.m_yDet); //Output
	   IntersectionPointCoord(*node,goodNode);
	     xadd = node->m_xDet;
	     yadd = node->m_yDet;
	  zadd = node->m_z_Det;
	   //  printf("We found %f, %f, %f\n", xadd, yadd, zadd);
	  //	break;
	  //	}
      } else {
	printf("EEROR, somethign to implement asap \n");
      }
      	         
    }

    
    m_x.insert(m_x.begin() + vecindex,xadd);
    m_y.insert(m_y.begin() + vecindex,yadd);
    m_z.insert(m_z.begin() + vecindex,zadd);
    m_r.insert(m_r.begin() + vecindex,node->m_r);


    double newtheta = node->m_thetaDeg;

    if(m_length > 2){ //Accounting for the 2pi change 
      double prevtheta = m_theta.back();
      if(prevtheta > 150 && newtheta < -150)
	newtheta += 360.;
      else if (prevtheta < -150 && newtheta > 150)
	newtheta -= 360.;
    }
    m_theta.insert(m_theta.begin() + vecindex, newtheta);
    
  } else {

    if(m_seenVirtual == true && m_listSkewed.size() > 0){
      // Correcting skewed
      /* int virtIdx = gr.Find(m_lastVirtual);
      GridNode &lastVirtNode = Ingrid[virtIdx];
     
      float x_diff = node->m_x - lastVirtNode.m_x;
      float y_diff = node->m_y - lastVirtNode.m_y;
      //  float z_diff = node->m_z - lastVirtNode.m_z;

      x_diff /= static_cast<float>(m_listSkewed.size()+1);
      y_diff /= static_cast<float>(m_listSkewed.size()+1);
      //  z_diff /= static_cast<float>(m_listSkewed.size()+1);

      float xInc = lastVirtNode.m_x + x_diff;
      float yInc = lastVirtNode.m_y + y_diff;
      //   float zInc = lastVirtNode.m_z + z_diff;

      /// Correct xyz-coordinates of the skewed nodes 
      for(size_t m = 0; m < m_listSkewed.size(); ++m) {
	int idx = gr.Find(m_listSkewed[m]);
	GridNode &skewedToproc = Ingrid[idx];
	if(layer != lastVirtNode.m_Layer){
	  skewedToproc.m_xDet = xInc;
	  skewedToproc.m_yDet = yInc;
	}
	//	skewedToproc.m_z_Det = zInc;
	auto it = find(m_memberList->begin(), m_memberList->end(), m_listSkewed[m]); 
  
	int index = distance(m_memberList->begin(), it); 
    
	m_x[index] = skewedToproc.m_xDet;
	m_y[index] = skewedToproc.m_yDet;
	//	m_z[index] = skewedToproc.m_z_Det;
	
	std::pair<float, float> r_Theta;
	float theta_deg = Cartesian_To_Polar(skewedToproc.m_xDet, skewedToproc.m_yDet, r_Theta);
	m_r[index] = r_Theta.first;

	double prevtheta = m_theta.back();
	double newtheta = theta_deg;
	if(m_length > 2){
	  //	  printf("%lf, %lf\n",  newtheta,m_theta[0] );
	  if(prevtheta > 150 && newtheta  < -150)
	    newtheta += 360.;
	  else if (prevtheta < -150 && newtheta > 150)
	    newtheta -= 360.;
	  //  printf("%lf, %lf\n",  newtheta,m_theta[0] );
	}
	m_theta[index] = newtheta;
	
	//	printf("New node %d with %lf and %lf, (r %lf, theta %lf) \n",m_listSkewed[m], skewedToproc.m_xDet ,  skewedToproc.m_yDet,m_r[index], theta_deg ); 
	xInc += x_diff;
	yInc += y_diff;
	//	zInc += z_diff;
      }
      

      //printf("We have %d virtuals and %d skewed to correct\n", m_prevVirtuals.size(), m_listSkewed.size());
      GridNode &anchorPrev = Ingrid[gr.Find(m_prevVirtuals[0])];

      if(m_prevVirtuals.size() > 1){
	for(size_t p =1; p < m_prevVirtuals.size();p++){
	  anchorPrev.m_x += Ingrid[gr.Find(m_prevVirtuals[p])].m_x;
	  anchorPrev.m_y += Ingrid[gr.Find(m_prevVirtuals[p])].m_y;
	  anchorPrev.m_z += Ingrid[gr.Find(m_prevVirtuals[p])].m_z;
	}
	anchorPrev.m_x /= (float) m_prevVirtuals.size();
	anchorPrev.m_y /= (float)m_prevVirtuals.size();
	anchorPrev.m_z /= (float) m_prevVirtuals.size();
      }

      if(anchorPrev.m_Layer == layer){
	//	printf("Both virtuals are on same layers, better find something else\n");
      } else{
	//	printf("Find nodes per layer\n");
	std::vector<GridNode *>    layerN1;     
	std::vector<GridNode *>    layerN2;
	int curLayer;
	for(size_t p =0; p < m_listSkewed.size();p++){
	  if(p == 0)
	    curLayer = Ingrid[gr.Find(m_listSkewed[p])].m_Layer;
	  if(Ingrid[gr.Find(m_listSkewed[p])].m_Layer == curLayer){
	    //   printf("Push %d on firt layer nodes\n", m_listSkewed[p]);
	    layerN1.push_back(&Ingrid[gr.Find(m_listSkewed[p])]);
	  } else if (labs(Ingrid[gr.Find(m_listSkewed[p])].m_Layer - curLayer) <2){
	    //    printf("Push %d on second layer nodes\n", m_listSkewed[p]);
	    layerN2.push_back(&Ingrid[gr.Find(m_listSkewed[p])]);
	  } else {
	    //  printf("We have a too large layer ???? \n");
	  }
	}

	GridNode anchorInter;
	if(layerN1.size() > 0){
	  //	  printf("Correcting first layer\n");
	  anchorInter =  *layerN1[0];
	  for (size_t p = 1; p < MIN(layerN1.size(),3); p++){
	    anchorInter.m_x += layerN1[p]->m_x;
	    anchorInter.m_y += layerN1[p]->m_y;
	    anchorInter.m_z += layerN1[p]->m_z;
	  }
	  anchorInter.m_x /= (float) MIN(layerN1.size(),3);//interNodesL1.size();
	  anchorInter.m_y /= (float) MIN(layerN1.size(),3);// interNodesL1.size();
	  anchorInter.m_z /= (float) MIN(layerN1.size(),3);//;interNodesL1.size();
      
	  //  printf("Finding intersection with %f, %f, %f\n", anchorInter.m_x, anchorInter.m_y, anchorInter.m_z);

	  PointsLineIntersect( anchorInter, anchorPrev.m_x, node->m_xDet,
					  anchorPrev.m_y, node->m_yDet);
	
	  //  printf("Intersect point is %f, %f, %f\n",anchorInter.m_xDet, anchorInter.m_yDet, anchorInter.m_z_Det);
	  for (size_t p = 0; p <layerN1.size(); p++){
	    auto it = find(m_memberList->begin(), m_memberList->end(), layerN1[p]->m_detID); 
	    int index = distance(m_memberList->begin(), it); 
	    
	    layerN1[p]->m_xDet = m_x[index] = anchorInter.m_xDet;
	    layerN1[p]->m_yDet = m_y[index] = anchorInter.m_yDet;
	    layerN1[p]->m_z_Det = m_z[index] = anchorInter.m_z_Det;
	  }
 
	}

	if(layerN2.size() > 0){
	  //  printf("Correcting second layer\n");
	  anchorInter =  *layerN2[0];
	  for (size_t p = 1; p < MIN(layerN2.size(),3); p++){
	    anchorInter.m_x += layerN2[p]->m_x;
	    anchorInter.m_y += layerN2[p]->m_y;
	    anchorInter.m_z += layerN2[p]->m_z;
	  }
	  anchorInter.m_x /= (float) MIN(layerN2.size(),3);//interNodesL1.size();
	  anchorInter.m_y /= (float) MIN(layerN2.size(),3);// interNodesL1.size();
	  anchorInter.m_z /= (float) MIN(layerN2.size(),3);//;interNodesL1.size();
      
	  //  printf("Finding intersection with %f, %f, %f\n", anchorInter.m_x, anchorInter.m_y, anchorInter.m_z);

	  PointsLineIntersect( anchorInter, anchorPrev.m_x, node->m_xDet,
			       anchorPrev.m_y, node->m_yDet);
	
	  // printf("Intersect point is %f, %f, %f\n",anchorInter.m_xDet, anchorInter.m_yDet, anchorInter.m_z_Det);
	  for (size_t p = 0; p <layerN2.size(); p++){
	    auto it = find(m_memberList->begin(), m_memberList->end(), layerN2[p]->m_detID); 
	    int index = distance(m_memberList->begin(), it); 
	    
	    layerN2[p]->m_xDet = m_x[index] = anchorInter.m_xDet;
	    layerN2[p]->m_yDet = m_y[index] = anchorInter.m_yDet;
	    layerN2[p]->m_z_Det = m_z[index] = anchorInter.m_z_Det;
	  }
 
	}
      }
      m_listSkewed.clear();
      m_prevVirtuals.clear();
    }
    // printf("%d, %d, %d, %lf, %lf \n", vecindex, m_x.size(), m_memberList->size(),node->m_xDet,node->m_yDet);
    m_x.insert(m_x.begin() + vecindex, node->m_xDet);
    m_y.insert(m_y.begin() + vecindex, node->m_yDet);
    m_z.insert(m_z.begin() + vecindex, node->m_z_Det);     
    m_r.insert(m_r.begin() + vecindex, node->m_r);
    
   
    double newtheta = node->m_thetaDeg;
    if(m_length > 2){
      double prevtheta = m_theta.back();

      //  printf("%lf, %lf\n",  newtheta,m_theta[0] );
      if(prevtheta > 150 && newtheta  < -150)
	newtheta += 360.;
      else if (prevtheta < -150 && newtheta > 150)
	newtheta -= 360.;
      //    printf("%lf, %lf\n",  newtheta,m_theta[0] );
    }
    m_theta.insert(m_theta.begin() + vecindex,newtheta);
    
    m_seenVirtual = true;
    m_lastVirtual = node->m_detID;
    m_prevVirtuals.push_back(id);
  }
	  
  m_headNode = id;
  //    printf("ADDING %d to track %d, layer %d \n", id, m_id, layer);

}

*/
