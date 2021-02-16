
#include <iostream>
#include <set>
#include <sstream>
#include <algorithm>
#include  <cmath>
#include <fstream>
#include <vector>
#include <gsl/gsl_poly.h>

// Root headers
#include "TFile.h"
#include "TNtuple.h"
#include "TStopwatch.h"
#include "TH1.h"

#include "reconstruction.h"
#include "gridNode.h"
#include "logc.h"
#include "simon_functions.h"
#include "path_queue.h"

void findEasyTracks (CoordGrid &gr, std::vector< GridNode > &Ingrid, std::vector < PathCandidate* > &tracklets, std::vector<int> &activeId, char *visited, int *candidateId){

  //std::vector< GridNode > &Ingrid = gr.m_grid;  

  for(unsigned int n = 0; n < activeId.size(); ++n) {
      
    std::vector<int> sameLayer;
    std::vector<int> nextLayer;
    std::vector<int> prevLayer;
    std::vector<int> nextVirt;
    std::vector<int> prevNodes;
    std::vector<int> *v;
      
    int dir 	       	= 0;
    int curId 		= activeId[n];
    int curIdx 		= gr.Find(curId);
    GridNode *currentNode 	= &Ingrid[curIdx];

    /* ++++++++++++++++++++++++++++++++++++++++++++++++++++ */   
    /* BEGIN WITH NODES WITH ONE NEIGHBOR OR IN LAYER LIMIT */     
    /* ++++++++++++++++++++++++++++++++++++++++++++++++++++ */
      
    if(visited[curId] == 0 &&
       ((currentNode->m_LayerLimit == 1 && currentNode->m_neighbors.size() <= 2)
	||  currentNode->m_neighbors.size() == 1)){ 
	
      int n_neighbors = currentNode->m_neighbors.size();
      int curLayer 	= currentNode->m_Layer;

      dbgconnect("Starting node %d has %d neighbors", curId, n_neighbors);

      int    		n_connected = 0;
      int 		neighId;
      int 		neighIdx;
      GridNode 	*neighNode;
      bool   		cond = true;

      /* Find new neighbors */
      cond = sortNeighbors(gr, currentNode, prevLayer, sameLayer, nextLayer, nextVirt, visited,  &dir);
	
      if(cond == false){ // if no simple neighbors, stop 
	resetLists(visited, prevLayer, sameLayer, nextLayer);
	continue;
      }

      // Create a new candidate path
      PathCandidate *cand 	= new PathCandidate();// Create a new candidate
      cand->m_id 		= (*candidateId)++;// Set id
      cand->m_tailNode 	= curId;
      visited[curId] 		= 4;
      cand->insertNewNode(gr,currentNode, cand->m_memberList->end());
      prevNodes.push_back(curId);


      /* ++++                    ++++ */
      /* Start putting nodes together */
      /* ++++                    ++++ */
	
      while(cond){
       	  
	if(n_neighbors == 1){ /* 1 NEIGHBOR */

	  if (dir == UP){
	    v = &nextLayer; // Going up from node %d to node %d", curId, v->at(0));
	  } else if (dir == DOWN){
	    v = &prevLayer; //Going down from node %d to node %d", curId, v->at(0));
	  } else {
	    v = &sameLayer; //Going same level from node %d to node %d", curId, v->at(0));
	  }
	     
	    
	  if(nextVirt.size() > 0){ // If there were some virtual nodes to add	      
	    for (size_t i = 0; i < nextVirt.size(); i++){
	      neighId   = nextVirt[i];
	      neighIdx  = gr.Find(neighId);
	      neighNode = &Ingrid[neighIdx];
	      cand->insertNewNode(gr,neighNode,cand->m_memberList->end());
	      visited[neighId] = 4;
	      n_connected++;
	      // removeIdFromNeigh(neighNode, &prevNodes, curId);
	    }
	    prevNodes.insert(prevNodes.end(),  nextVirt.begin(),  nextVirt.end());
	    nextVirt.clear();
	  } // end if virtual nodes
	    
	  neighId    = v->at(0); neighIdx   = gr.Find(neighId); neighNode  = &Ingrid[neighIdx];
	  cand->insertNewNode(gr,neighNode,cand->m_memberList->end());
	  visited[neighId] = 4;
	  n_connected++;
	    
	  removeIdFromNeigh(neighNode, &prevNodes, curId);

	  curId       = neighId;  curIdx      = neighIdx;  currentNode = neighNode;
	  curLayer    = currentNode->m_Layer;
	    
	  nextLayer.clear();  sameLayer.clear();   prevLayer.clear();    prevNodes.clear();
	  prevNodes.push_back(curId);
	  dir = 0;

	  cond = sortNeighbors(gr, currentNode, prevLayer, sameLayer, nextLayer, nextVirt, visited,  &dir);
	  n_neighbors = sameLayer.size() + prevLayer.size() + nextLayer.size();
	  // info("%d nodes were connected,  %d found for the next step \n", n_connected, n_neighbors);
	  n_connected = 0;	    
	} // end if 1 Neighbor



	else if (sameLayer.size() > 0 ){   /* Same layer neighbors */

	  v = dir == UP ? &nextLayer: &prevLayer; // find direction	    

	  if( !areAdjacent(gr, v) ){ // If next layer neighbors are not adjacent
	    cond = false;	      
	  } else { 

	    int candId 	 = sameLayer[0];
	    int candIdx        = gr.Find(candId);
	    GridNode *candNode = &Ingrid[candIdx];
	    curLayer           = candNode->m_Layer;
	    sameLayer.clear();

	    //info("Still on the same layer, investigating node %d", candId);
	    
	    for(size_t i = 0; i < candNode->m_neighbors.size(); i++){//CHECK if neighbors on next layer are adjacent 
	      neighId = candNode->m_neighbors[i];	      
	      neighIdx  = gr.Find(neighId);
	      neighNode = &Ingrid[neighIdx];
		
	      if (neighId == curId || neighNode->m_type == GridNode::VIRTUAL_NODE) continue;
			      
	      int haveNeigh = 0;
		
	      for (size_t j = 0; j < v->size(); j++){		  
		if(v->at(j) == neighId || (neighNode->IsNeighboring(v->at(j))))
		  haveNeigh = 1;		  
	      }
		
	      if(haveNeigh == 0){ // If not adjecent, we quite
		cond = false;
		break;
	      }
		
	    } // for neighbors of candidates
		
	    if (cond == true){ // All seem to belong to the same track, we can continue
		
	      //	info("All neighbors look good, let's insert this one !");

	      cand->insertNewNode(gr,candNode,cand->m_memberList->end());
	      visited[candId] = 4;
	      n_connected++;

	      removeIdFromNeigh(candNode, &prevNodes, curId);
	      
	      prevNodes.push_back(candId);
	      curId       = candId;
	      curIdx      = candIdx;
	      currentNode = candNode;
	      cond = sortNeighbors(gr, currentNode, prevLayer, sameLayer, nextLayer, nextVirt, visited,  &dir);

	      n_neighbors = sameLayer.size()+prevLayer.size()+nextLayer.size();
	      //	info("%d nodes were connected,  %d found for the next step \n\n", n_connected, n_neighbors);
	      n_connected = 0;		
	    } // IF COND TRUE
	      
	    else {

	      //	info("Neighbors not connected... \n\n", n_connected, n_neighbors);
	      if(visited[candId] != 4)
		visited[candId] = 0;
		   
	      cand->m_headNeigh.push_back(candId); // Push for future fitting
		
	    } // ELSE (COND IS FALSE)
	  } // ELSE (neighbors are adjacent
	} // ELSE IF (same layer size > 0)



	else if (n_neighbors > 1) {  	  /* Let's have a look to the other neighbors */

	  if (dir == UP){ // we are going up 	      
	    v = &nextLayer;
	  }  else if (dir == DOWN){ // we are going down
	    v = &prevLayer;
	  } else // we have no idea where we are going
	    dbgconnect("WHAT IS THE DIRECTION NOW?");
	           
	  if(areAdjacent(gr, v)){ // All neighbors are adjacent ??
	      
	    //nfo("Adding %d nodes to the CM", v->size());
	      
	    if(nextVirt.size() > 0){ // taking care of the virtual nodes
	      for (size_t i = 0; i < nextVirt.size(); i++){
		neighId   = nextVirt[i];
		neighIdx  = gr.Find(neighId);
		neighNode = &Ingrid[neighIdx];
		cand->insertNewNode(gr, neighNode, cand->m_memberList->end());
		visited[neighId] = 4;
		n_connected++;
		//	removeIdFromNeigh(neighNode, &prevNodes, curId);		  
	      }
		
	      prevNodes.insert(prevNodes.end(),  nextVirt.begin(),  nextVirt.end());
	      nextVirt.clear();
	    }
	      
	    std::vector<int> lookneigh;

	    for (size_t i = 0; i < v->size(); i++){ // add the next nodes
	      neighId   = v->at(i);
	      neighIdx  = gr.Find(neighId);
	      neighNode = &Ingrid[neighIdx];
	      cand-> insertNewNode(gr,neighNode, cand->m_memberList->end());
	      visited[neighId] = 4;
	      n_connected++;
	      lookneigh.push_back(neighId);
	      removeIdFromNeigh(neighNode, &prevNodes, curId);		  
	      removeIdFromNeigh(neighNode, v, curId);		  
	    }
	      
	    nextLayer.clear(); sameLayer.clear();    prevLayer.clear();    prevNodes.clear();
	      
	    dir = 0;
	     
	    curId       = lookneigh[0];
	    curIdx      = gr.Find(curId);
	    currentNode = &Ingrid[curIdx];
	    curLayer    = currentNode->m_Layer;
	      
	    n_neighbors = 0;
	      
	    for(size_t i = 0; i < lookneigh.size(); i++){ // Looking for new neighbors
	      int id         = lookneigh[i];
	      int idx        = gr.Find(id);
	      GridNode *node = &Ingrid[idx];
	      prevNodes.push_back(id);
	      cond = sortNeighbors(gr, node, prevLayer, sameLayer, nextLayer, nextVirt, visited,  &dir);
	    }

	    n_neighbors = sameLayer.size()+prevLayer.size()+nextLayer.size();
	    //     info("%d nodes were connected,  %d found for the next step \n", n_connected, n_neighbors);
	    n_connected = 0;
	     
	  } // IF ADJACENT

	  else {
	    //  info("Some of these nodes are node adjacent", v->size());
	    cond = false;
	  }
	    
	      
	} // END IF neighbors size > 0

	else { // We have no more neighbors

	  int firstId         = cand->m_tailNode;
	  int firstIdx        = gr.Find(firstId);
	  GridNode *firstNode = &Ingrid[firstIdx];

	  int lastId         = cand->m_headNode;
	  int lastIdx        = gr.Find(lastId);
	  GridNode *lastNode = &Ingrid[lastIdx];
	    
	  if(n_neighbors == 0) {
	      
	    dbgconnect("No more neighbors in sight, are we finished here ?");
	      
	    if((cand->m_minLayer == 0 && cand->m_maxLayer > 21) || (firstNode->m_LayerLimit == 1 && lastNode->m_LayerLimit == 1)){		 
	      dbgconnect("track goes through all layers or makes a loop, likily finished");		 
	      cand->m_finished = 3;		 
	    } else if(labs(currentNode->m_SectorLimit) > 0 || cand->m_isOnSectorLimit){
	      dbgconnect("Track is on sector limit, might have a connection somewhere else");
	      cand->m_finished = 2;
	      cand->m_isOnSectorLimit= true;		 
	    } else {		 
	      dbgconnect("Candidate has no more neighbors, but doesn't seem finished");
	      cand->m_finished = 2;		 
	    }
	      
	  } // end if n_neighbors == 0

	    
	  cond = false;
	    
	}

	  	  
	if(cond == false){ // This track is finished, but let's push neighbors as we need to fit in the next step

	  for (size_t i = 0; i < sameLayer.size(); i++)
	    cand->m_headNeigh.push_back(sameLayer[i]);
	      
	  for (size_t i = 0; i < nextLayer.size(); i++)
	    cand->m_headNeigh.push_back(nextLayer[i]);

	  for (size_t i = 0; i < prevLayer.size(); i++)
	    cand->m_headNeigh.push_back(prevLayer[i]);

	  cand->m_headNeigh.insert((cand->m_headNeigh).begin(),  (nextVirt).begin(),  (nextVirt).end());	
	  resetLists(visited, prevLayer, sameLayer, nextLayer);
	}

      }// WHILE COND

	
      if(cand->m_length > 2){
	dbgconnect("Pushing cm %d with length %d, tail node %d, head node %d, min layer %d, max layer %d, IsOnSectorLimit %d, finished ? %d. ", cand->m_id, cand->m_length, cand->m_tailNode, cand->m_headNode,cand->m_minLayer, cand->m_maxLayer, cand->m_isOnSectorLimit, cand->m_finished);
	tracklets.push_back(cand);
      } else {
	dbgconnect("Not a good cm %d",cand->m_headNode);
	for(size_t i = 0; i < (cand->m_memberList)->size(); i++){   
	  visited[(cand->m_memberList)->at(i)] = 0;
	}
	  
	delete cand;
	(*candidateId)--;
      } // Else Candidate large enough
	
    } // for node with a single neighbor
       	     
  } // For active nodes

}





/* IF ASSIGNING REMAINING NODES // TO CHEC LATER

   for(unsigned int n = 0; n < nactiveReal; ++n) {
	
      int curId      	= activeId[n];
      std::vector< int > toAdd;
      std::vector< int > possiCand;
      std::vector< double > distToCand;
      std::vector< int > nodeId;
      std::vector< int > curr;
      std::vector< int > next;

	
      if(visited[curId] != 0) continue;
	
      info("Remaining node %d", curId);
      int curIdx = gr.Find(curId);
      GridNode *currentNode = &Ingrid[curIdx];
	
      int n_neighbors = currentNode->m_neighbors.size();
      if(n_neighbors == 0){
	debug("LOST CHILD");
	visited[curId] = 4;
	continue;
      }
	
      curr.insert(curr.end(),  currentNode->m_neighbors.begin(),   currentNode->m_neighbors.end());
      toAdd.push_back(curId);

      bool cond = true;
	
      while(cond){
	
	std::sort (curr.begin(), curr.end());
	    
	for (int i = 0; i < curr.size(); i++){
	  int neighId = curr[i];
	  int neighIdx = gr.Find(neighId);
	  GridNode *neighNode = &Ingrid[neighIdx];
	  
	  // debug("Neighbor %d", neighId);
	  
	  if (visited[neighId] == 4){
	    double dist = distanceBetweenTube(*currentNode, *neighNode);
	    info("Neigh %d, One possible CM is %d, tube distance %lf",neighId, neighNode->m_cm[0], dist);
	    possiCand.push_back(neighNode->m_cm[0]);
	    distToCand.push_back(dist);
	    nodeId.push_back(neighId);
	  }
	  else if(!(std::find(toAdd.begin(), toAdd.end(),neighId)!= toAdd.end())){	      
	    toAdd.push_back(neighId);
	    next.push_back(neighId);
	    info("Node %d added for next loop", neighId);
	  }
	    
	  //std::string dummy;
	  //std::cout << "Enter to continue..." << std::endl;
	  //std::getline(std::cin, dummy);
	}

	if(possiCand.size() == 1){
	  debug("We found a path candidate for this node");
	  int potCm = possiCand[0];
	  const auto p = std::find_if(tracklets.begin(), tracklets.end(), [potCm](const PathCandidate *obj){ return obj->m_id == potCm; } );
	  std::vector<int>::iterator it = std::find(((*p)->m_memberList)->begin(), ((*p)->m_memberList)->end(), nodeId[0]);

	  (*p)->insertNewNode(gr,currentNode, it+1);
	  visited[curId]= 4;
	  cond = false;
	}

	    
	else if (possiCand.size() > 1) {
	  debug("We found several path candidates for this node");
	  std::vector<double>::iterator result = std::min_element(distToCand.begin(),distToCand.end());
	  int minEltIdx = std::distance(distToCand.begin(), result);
	  // int minElement = *std::min_element(v.begin(), v.end());
	  debug("Node %d is closer with dist %lf and belongs to cm %d", nodeId[minEltIdx], distToCand[minEltIdx], possiCand[minEltIdx]);
	  int potCm = possiCand[minEltIdx];
	  const auto p = std::find_if(tracklets.begin(), tracklets.end(), [potCm](const PathCandidate *obj){ return obj->m_id == potCm; } );

	  std::vector<int>::iterator it = std::find(((*p)->m_memberList)->begin(), ((*p)->m_memberList)->end(), nodeId[minEltIdx]);
	  int EltIdx = std::distance(((*p)->m_memberList)->begin(), it);


	  (*p)->insertNewNode(gr,currentNode, it);
	  visited[curId]= 4;
	  cond = false;
	}

	    
	else if (next.size() > 0){
	  debug("WELL let's contiinue");
	  curr.clear();
	  curr.insert(curr.end(),  next.begin(),   next.end());
	  next.clear();
	}

	else {
	  debug("Seems that no neighbors belongs to track?, let's create a new one");
	  PathCandidate *cand 	= new PathCandidate();// Create a new candidate
	  cand->m_id 		= candidateId++;// Set id
	  cand->m_tailNode 	        = curId;
	  visited[curId] 		= 4;

	  for (int i = 0; i < toAdd.size(); i++){
	    int nId = toAdd[i];
	    int nIdx = gr.Find(nId);
	    GridNode *nNode = &Ingrid[nIdx];
	    cand->insertNewNode(gr,nNode, cand->m_memberList->end());
	    visited[nId] = 4;
	  }

	  info("Pushing cm %d: \n               length is %d, \n     tail node %d \n     head node %d \n       Min layer %d, \n              Max layer %d, \n               IsOnSectorLimit %d, \n               Finished ? %d. ", cand->m_id, cand->m_length, cand->m_tailNode, cand->m_headNode,cand->m_minLayer, cand->m_maxLayer, cand->m_isOnSectorLimit, cand->m_finished);

	  tracklets.push_back(cand);
	  cond = false;
	      
	}	      
	     	    	    
      } // WHILE COND
    } // for nodes
      */


//____


void CompZCoordinates(CoordGrid &hitMap, PathCandidate *trk)
{
 
 
  dbgtrkz("Re-determining the z coordinates for the track", trk->m_id);

  // MVD list of available elements in the graph.
  // std::vector< GridNode > &MVDPoints = hitMap.m_MVD_grid;
  //
  // STT lists of elements available in the graph.
  std::vector< GridNode > &Ingrid = hitMap.m_grid;
  // Local variables
  std::set<int>::iterator stt_It;
  // Track (Connected components) )loop

  // Current track object  std::vector<double> z =  std::vector<double>( trk->m_z );

  std::vector<int> const *vect = trk->m_memberList;

  // Fnd the most outer virtual node or the max Z
  // size_t maxLayer = std::numeric_limits<size_t>::min();
  // float maxVirtZVal = std::numeric_limits<float>::min();
  //float  lastVirtualZCoord = 0.00;

  /*dbgtrkz("Current coordinates are:");
  for( size_t i = 0; i < vect->size(); i++){
    int nodeID = vect->at(i);
    size_t node_index = hitMap.Find(nodeID);
    GridNode const &node = Ingrid[node_index];
     dbgtrkz("%d, x %lf, y %lf, z %lf", vect->at(i),node.m_xDet, node.m_yDet,node.m_z_Det);

     }*/
  std::vector<double> &x =  trk->m_x;
  std::vector<double> &y =  trk->m_y;
  std::vector<double> &z =  trk->m_z;
  
  std::vector< GridNode* > prevNodes;
  std::vector< GridNode* > interNodesL1;
  std::vector< GridNode* > interNodesL2;
  std::vector< GridNode* > nextNodes;

  std::vector< GridNode* > prevVirt;
  std::vector< GridNode* > sameLayer;
  std::vector< GridNode* > nextVirt;

  std::vector< double >  xPts;
  std::vector< double >  yPts;
  std::vector< double >  zPts;

  int cntBef = 0, cntAft = 0;
  int prevLayer, nextLayer, curLayer =0;
  GridNode anchorPrev, anchorNext;
  int visitVirt = 0;
  int nextPara = 0, prevPara =0;
  int dir  = 0, loop = 0;
  dbgtrkz("Finding virtuals for interpolation");
  int numElts = 2;
  for(size_t i = 0; i < vect->size(); i++){
    int nodeID = vect->at(i);
    size_t node_index = hitMap.Find(nodeID);
    GridNode  &node = Ingrid[node_index];

    /*    if(i > 0 && node.m_Layer > (*prevNodes.back()).m_Layer)
      dir = 1;
    else if (i > 0 && node.m_Layer < (*prevNodes.back()).m_Layer)
      dir = -1;
    else
    dir = 0;*/
    
    if(node.m_type == GridNode::VIRTUAL_NODE && !visitVirt){ // Firt visit of virtual

      prevVirt.push_back(&node);
      prevLayer = (prevNodes[prevNodes.size()-1])->m_Layer;
      prevPara = (prevNodes[prevNodes.size()-1])->m_type != GridNode::STT_TYPE_SKEW? 1: 0;

      dbgtrkz("Previous node %d before virtual had layer %d and is a parallel ? %d", (prevNodes[prevNodes.size()-1])->m_detID, prevLayer, prevPara);

      while(true){ // find all virtuals on same layer
	nodeID = vect->at(i+1);
	node_index = hitMap.Find(nodeID);
	if(Ingrid[node_index].m_type == GridNode::VIRTUAL_NODE){
	  dbgtrkz("Adding an other virt %d from first layer", nodeID);
	  prevVirt.push_back(&Ingrid[node_index]);
	  i++;
	} else
	  break;
      }
      
      visitVirt = 1;
    }

    else if( visitVirt == 1 && node.m_type != GridNode::VIRTUAL_NODE){

      dbgtrkz("Next node %d on layer %d", nodeID, node.m_Layer);

      if(labs(node.m_Layer - prevLayer) == 1 ){
	dbgtrkz("Pushing on first layer");
	interNodesL1.push_back(&node);
      }
      else if (labs(node.m_Layer - prevLayer) == 2){
	dbgtrkz("Pushing on second layer");
	interNodesL2.push_back(&node);
      }
      else {
	//	if(labs(node.m_Layer - prevLayer) > 3)
	//	  anchorNext.m_detID = -1;
	//	else{
	dbgtrkz("We are too far in layers, missing virtual, creating one");
	GridNode &prevNode = *interNodesL2.back();
	GridNode Dummy_coord;
	// Find intersection point.
	  
	IntersectionPointSkeSke(hitMap, prevNode, node, Dummy_coord);
	dbgtrkz("Created virtual at %f, %f, %f", Dummy_coord.m_x,  Dummy_coord.m_y, Dummy_coord.m_z);
	nextVirt.push_back(&Dummy_coord);
	nextNodes.push_back(&node);
	nextLayer = node.m_Layer;
	nextPara =  node.m_type == GridNode::STT_TYPE_SKEW ? 0 : 1;
	//	}
	visitVirt =2;

      }
    }

    else if(visitVirt && node.m_type == GridNode::VIRTUAL_NODE){
      nextVirt.push_back(&node);
      while(true){ // find al virtuals on same layer
	nodeID = vect->at(i+1);
	node_index = hitMap.Find(nodeID);
	if(Ingrid[node_index].m_type == GridNode::VIRTUAL_NODE){
	  nextVirt.push_back(&Ingrid[node_index]);
	  dbgtrkz("Adding an other virt %d from second layer", nodeID);
	  i++;
	} else {
	  nextLayer = Ingrid[node_index].m_Layer;
	  nextPara =  Ingrid[node_index].m_type == GridNode::STT_TYPE_SKEW ? 0 : 1;
	  nextNodes.push_back(&Ingrid[node_index]);
	  i++;
	  break;
	}
      }
      visitVirt = 2;
    } else{
      prevNodes.push_back(&node);
    }

    if(visitVirt == 2 || (i == vect->size()-1 && visitVirt ==1)){

      if(visitVirt == 2){
	dbgtrkz("Reached second layer of virtuals, we have our sandwhich from layer %d to layer %d", prevLayer, nextLayer);
	
	if(prevPara){
	  cntBef = 0;
	  for(int j = prevNodes.size()-1,  layer = prevLayer;  layer == prevLayer; j--){
	    sameLayer.push_back(prevNodes[j]);	    	    
	    layer = (prevNodes[j-1])->m_Layer;
	    cntBef ++;
	  }
	  dbgtrkz("Previous anchor needs to be a parallel tube, taking average from %d tubes", cntBef);

	  anchorPrev = *(sameLayer[0]);

	  if(sameLayer.size() > 1){
	    for(size_t p =1; p < sameLayer.size();p++){
	      anchorPrev.m_x += sameLayer[p]->m_x;
	      anchorPrev.m_y += sameLayer[p]->m_y;
	      anchorPrev.m_z += sameLayer[p]->m_z;
	    }
	    anchorPrev.m_x /= (float) sameLayer.size();
	    anchorPrev.m_y /= (float) sameLayer.size();
	    anchorPrev.m_z /= (float) sameLayer.size();
	  }
	}
      
	else if (anchorNext.m_detID != -1){
	  dbgtrkz("Previous anchor will be a virtual, we can take the previous anchor");
	  anchorPrev = anchorNext;
	  xPts.push_back((double) anchorPrev.m_x);
	  yPts.push_back((double) anchorPrev.m_y);
	  zPts.push_back((double) anchorPrev.m_z);
	} else {
	  dbgtrkz("track start on skewed...");
	  anchorPrev= *prevVirt[0];
	  printf("%f\t",prevVirt[0]->m_z);

	  for(size_t p =1; p < prevVirt.size();p++){
	    anchorPrev.m_x += prevVirt[p]->m_x;
	    anchorPrev.m_y += prevVirt[p]->m_y;
	    anchorPrev.m_z += prevVirt[p]->m_z;
	    printf("%f\t",prevVirt[p]->m_z);
	  }
	  printf("\n");
	  anchorPrev.m_x /= (float) prevVirt.size();
	  anchorPrev.m_y /= (float) prevVirt.size();
	  anchorPrev.m_z /= (float) prevVirt.size();
	  xPts.push_back((double) anchorPrev.m_x);
	  yPts.push_back((double) anchorPrev.m_y);
	  zPts.push_back((double) anchorPrev.m_z);
	}

	dbgtrkz("First anchor has coord x %f, y %f, z %f", anchorPrev.m_x, anchorPrev.m_y, anchorPrev.m_z);

	if(nextPara){
	  dbgtrkz("Next anchor needs to be a parallel tube");
	  cntAft = 0;
	  while(true){
	    nodeID = vect->at(i+1);
	    node_index = hitMap.Find(nodeID);
	    if( Ingrid[node_index].m_Layer == nextLayer && Ingrid[node_index].m_type != GridNode::VIRTUAL_NODE){
	      nextNodes.push_back(&Ingrid[node_index]);
	      dbgtrkz("Adding %d, layer %d", Ingrid[node_index].m_detID,Ingrid[node_index].m_Layer  );
	      cntAft ++;
	      i++;
	    } else{
	      dbgtrkz("Found node %d on layer %d, too far for now", nodeID ,Ingrid[node_index].m_Layer);
	      break;
	    }
	  }

	  anchorNext= *nextNodes[0];

	  dbgtrkz("Taking average of %lu tubes", nextNodes.size());

	  if(nextNodes.size() > 1){
	    for(size_t p =1; p < nextNodes.size();p++){
	      anchorNext.m_x += nextNodes[p]->m_x;
	      anchorNext.m_y += nextNodes[p]->m_y;
	      anchorNext.m_z += nextNodes[p]->m_z;
	    }
	    anchorNext.m_x /= (float) nextNodes.size();
	    anchorNext.m_y /= (float) nextNodes.size();
	    anchorNext.m_z /= (float) nextNodes.size();
	  }
	}

	else {
	  dbgtrkz("Next anchor needs to be a vitual point, we have %d", nextVirt.size());
	  anchorNext= *nextVirt[0];
	  for(size_t p =1; p < nextVirt.size();p++){
	    anchorNext.m_x += nextVirt[p]->m_x;
	    anchorNext.m_y += nextVirt[p]->m_y;
	    anchorNext.m_z += nextVirt[p]->m_z;
	  }
	  anchorNext.m_x /= (float) nextVirt.size();
	  anchorNext.m_y /= (float) nextVirt.size();
	  anchorNext.m_z /= (float) nextVirt.size();
	
	}

	dbgtrkz("Second anchor location is x %f, y %f, z %f", anchorNext.m_x, anchorNext.m_y, anchorNext.m_z);


	
	dbgtrkz("Using found anchors to find z intersection");

      } else if(anchorNext.m_detID != -1) {
	dbgtrkz("We reached the end of the vector size, but we can use the previous anchors");
	xPts.push_back((double) anchorNext.m_x);
	yPts.push_back((double) anchorNext.m_y);
	zPts.push_back((double) anchorNext.m_z);
      } else {
	dbgtrkz("Not enough points in virtuals for interpolations, going out");
	break;
      }

      if(prevLayer == nextLayer){
	dbgtrkz("We came back to the same layer of virtual, let's continue");
      } else if (	interNodesL1.size() > 3 || 	interNodesL2.size() > 3){
	dbgtrkz("Too many nodes in between, let's pass");
      }else {
	dbgtrkz("First layer %d with %lu nodes",interNodesL1[0]->m_Layer, interNodesL1.size());

	if(interNodesL1.size() == 0)
	  error("No node found on first layer");
	/*	else if(interNodesL1.size() >2){
	  dbgtrkz("We have more than 2 nodes on a same layer, warning");
	  float dist1 = sqrt(pow(anchorPrev.m_x - interNodesL1[0]->m_x,2) + pow(anchorPrev.m_y - interNodesL1[0]->m_y,2));
	  float dist2 = sqrt(pow(anchorPrev.m_x - interNodesL1.back()->m_x,2) + pow(anchorPrev.m_y - interNodesL1.back()->m_y,2));
	  if(dist2 < dist1)
	    std::reverse(interNodesL1.begin(), interNodesL1.end());
	    }*/
      
	GridNode anchorInter =  *interNodesL1[0];
	for (size_t p = 1; p < MIN(interNodesL1.size(),numElts); p++){
	  anchorInter.m_x += interNodesL1[p]->m_x;
	  anchorInter.m_y += interNodesL1[p]->m_y;
	  anchorInter.m_z += interNodesL1[p]->m_z;
	}
	anchorInter.m_x /= (float) MIN(interNodesL1.size(),numElts);//interNodesL1.size();
	anchorInter.m_y /= (float)MIN(interNodesL1.size(),numElts);// interNodesL1.size();
	anchorInter.m_z /= (float) MIN(interNodesL1.size(),numElts);//;interNodesL1.size();
      
	dbgtrkz("Finding intersection with %f, %f, %f", anchorInter.m_x, anchorInter.m_y, anchorInter.m_z);
	float xInt, yInt, zInt;
	LineLineIntersect( anchorPrev, anchorNext, anchorInter, xInt, yInt, zInt);
	dbgtrkz("Intersect point is %f, %f, %f", xInt, yInt, zInt);

	/*	for (size_t p = 0; p < interNodesL1.size(); p++){
	  interNodesL1[p]->m_xDet = xInt;
	  interNodesL1[p]->m_yDet = yInt;
	  interNodesL1[p]->m_z_Det = zInt;
	  }*/

	xPts.push_back((double) xInt);
	yPts.push_back((double) yInt);
	zPts.push_back((double) zInt);

	dbgtrkz("Second layer with %lu nodes", interNodesL2.size());

	if(interNodesL2.size() == 0)
	  error("No node found on second layer");
	/*	else if(interNodesL2.size() >3){
	  dbgtrkz("We have more than 3 nodes on a same layer, interpolation will be biased");
	   float dist1 = sqrt(pow(anchorNext.m_x - interNodesL2[0]->m_x,2) + pow(anchorNext.m_y - interNodesL2[0]->m_y,2));
	  float dist2 = sqrt(pow(anchorNext.m_x - interNodesL2.back()->m_x,2) + pow(anchorNext.m_y - interNodesL2.back()->m_y,2));
	  if(dist2 < dist1)
	  std::reverse(interNodesL2.begin(), interNodesL2.end());
	  }*/
      
	anchorInter =  *interNodesL2[0];

      
	for (size_t p = 1; p < MIN(interNodesL2.size(),numElts); p++){
	  anchorInter.m_x += interNodesL2[p]->m_x;
	  anchorInter.m_y += interNodesL2[p]->m_y;
	  anchorInter.m_z += interNodesL2[p]->m_z;
	}
	anchorInter.m_x /= (float) MIN(interNodesL2.size(),numElts);//interNodesL2.size();
	anchorInter.m_y /= (float) MIN(interNodesL2.size(),numElts);//interNodesL2.size();
	anchorInter.m_z /= (float) MIN(interNodesL2.size(),numElts);//interNodesL2.size();
	dbgtrkz("Finding intersection with %f, %f, %f", anchorInter.m_x, anchorInter.m_y, anchorInter.m_z);

	LineLineIntersect( anchorPrev, anchorNext, anchorInter, xInt, yInt, zInt);
	dbgtrkz("Intersect point is %f, %f, %f", xInt, yInt, zInt);

	/*	for (size_t p = 0; p < interNodesL2.size(); p++){
	  interNodesL2[p]->m_xDet = xInt;
	  interNodesL2[p]->m_yDet = yInt;
	  interNodesL2[p]->m_z_Det = zInt;
	  }*/

	xPts.push_back((double) xInt);
	yPts.push_back((double) yInt);
	zPts.push_back((double) zInt);
      }


      dbgtrkz("Prev Lay %d, next Laye %d, dir %d", prevLayer, nextLayer, dir);

      if(prevLayer < nextLayer)
	dir = 1;
      else if(nextLayer < prevLayer)
	dir = -1;
      else
	dir *=-1;

      dbgtrkz("Prev Lay %d, next Laye %d, dir %d", prevLayer, nextLayer, dir);

      prevNodes.clear();
      if(nextLayer == prevLayer){	
	prevNodes.insert(prevNodes.end(), interNodesL1.begin(), interNodesL1.end());
      }else {
	prevNodes.insert(prevNodes.end(), interNodesL2.begin(), interNodesL2.end());
      }
      if(!nextPara){
	dbgtrkz("Continuing with virtuals, dir %d", dir);
	prevLayer = dir == 1? nextLayer-1: nextLayer+1;//(prevNodes[prevNodes.size()-1])->m_Layer;
	prevPara = nextPara;//prevNodes[prevNodes.size()-1]).m_type != GridNode::STT_TYPE_SKEW ? 1:0;
	interNodesL1.clear();
	interNodesL2.clear();
	interNodesL1.insert(interNodesL1.end(), nextNodes.begin(), nextNodes.end());
	nextNodes.clear();
	nextVirt.clear();
	visitVirt = 1;
      } else {
	dbgtrkz("We reached the next layer of parallel tubes, resetting");
	interNodesL1.clear();
	interNodesL2.clear();
	interNodesL1.insert(interNodesL1.end(), nextNodes.begin(), nextNodes.end());
	nextNodes.clear();
	nextVirt.clear();
	visitVirt = 0;
	anchorPrev = anchorNext =  GridNode();
      }
      dbgtrkz("End for, prevLayer is %d, nextLayer is %d", prevLayer, nextLayer);
    }
   
  }

  if(xPts.size()>1){
    
    dbgtrkz("Coordinates found for fitting");
    for(size_t i = 0; i < zPts.size(); i++)
      dbgtrkz("x %lf, y %lf, z %lf", xPts[i],  yPts[i],  zPts[i]);


    std::vector<double> p;
      if(xPts.size() > 10){
      zPts.erase(zPts.begin());
      yPts.erase(yPts.begin());
      zPts.pop_back();
      yPts.pop_back();
      }

    p.push_back(0.);
    for (size_t i = 0; i <  zPts.size()-1; i++){ 
      double newval = p[i] + sqrt(pow(zPts[i+1]-zPts[i],2.) + pow(yPts[i+1]-yPts[i],2.));
      p.push_back(newval);
    }

    double *z_coef = polyFit(p, zPts, 1);
    double *y_coef = polyFit(p, yPts, 1);

    dbgtrkz("Fitted track is %f + p%f", z_coef[0], z_coef[1]);
    for(size_t i = 0; i < vect->size(); i++){
      //    int nodeID = vect->at(i);
      //    size_t node_index = hitMap.Find(nodeID);
      //    GridNode  &node = Ingrid[node_index];
      double pval = (y[i] - y_coef[0])/y_coef[1];
      z[i] = z_coef[0] + z_coef[1]*pval;
      dbgtrkz("%f, %f", y[i], z[i]);
      //node.m_z_Det = z_coef[0] + z_coef[1]*pval;
    }
    
  } else
    dbgtrkz("No enough good coord found");

  dbgtrkz("Finished cur track \n\n\n");
  
}

void TrackZ_CoordinatesDistNorm(CoordGrid &hitMap, std::vector<TrackObject*>* trks)
{
  if( (trks == 0) || (trks->size() == 0) ) {
    dbgtrkz("Input track list is empty. Terminating program.");
    //exit(EXIT_FAILURE);
    return;
  }
  dbgtrkz("Determining the z coordinates for points in tracks (Dist Norm), num of imput tracks %lu",trks->size());
	   
  
  // MVD list of available elements in the graph.
  // std::vector< GridNode > &MVDPoints = hitMap.m_MVD_grid;
  //
  // STT lists of elements available in the graph.
  std::vector< GridNode > &Ingrid = hitMap.m_grid;
  // Local variables
  std::set<int>::iterator stt_It;
  // Track (Connected components) )loop
  for(size_t tr = 0; tr < trks->size(); ++tr) {
    // Current track object
    TrackObject *track = trks->at(tr);
    // Set of the ID's of the participating STT nodes. Orderd
    // increasing (set property).
    std::set<int> const *sttTubes = track->m_sttComponent;
    if(!sttTubes) {
      continue;
    }
    // Fnd the most outer virtual node or the max Z
    size_t maxLayer = std::numeric_limits<size_t>::min();
    float maxVirtZVal = std::numeric_limits<float>::min();
    float  lastVirtualZCoord = 0.00;
    
    for(stt_It = sttTubes->begin(); stt_It != sttTubes->end();++stt_It){
      int nodeID = *stt_It;
      size_t node_index = hitMap.Find(nodeID);
      GridNode const &node = Ingrid[node_index];
      dbgtrkz("Node %d, z coord %f, layer %d", nodeID, node.m_z_Det,node.m_Layer);
      if( (node.m_type == GridNode::VIRTUAL_NODE || node.m_type ==GridNode::STT_TYPE_SKEW) &&
       	  (node.m_Layer > maxLayer) ) {
       	maxLayer = node.m_Layer;
       	lastVirtualZCoord = node.m_z_Det;
	dbgtrkz("We take this one (max layer bigger");
      }// END IF virtual
      // Determine the largest Z value (Virtual node)
      if( (node.m_type == GridNode::VIRTUAL_NODE|| node.m_type ==GridNode::STT_TYPE_SKEW) &&
       	  (node.m_z_Det > maxVirtZVal) ) {
	maxVirtZVal = node.m_z_Det;
	dbgtrkz("We take this one (max z coord bigger");

      }// END IF virtual

    }// END FOR (node loop)
    // Found the most outer virtual node
    /* Delta_z is equal to lastVirtualZCoord. Because we believe that
       the tracks start at (0,0,0)*/

    dbgtrkz("Last z value (dZ) %f, max layer %d", lastVirtualZCoord, maxLayer);
    // We need the correct order of tubes by their layer numbers. The
    // ids of virtual tubes does not correspond with their places in
    // the graph.
    if(true){
      dbgtrkz("Determine Z. Reordering nodes by layer.");

      PathQueue tubeIndexQueue;
      for(stt_It = sttTubes->begin(); stt_It != sttTubes->end();++stt_It){
	int nodeID = *stt_It;
	size_t node_index = hitMap.Find(nodeID);
	GridNode &tube = Ingrid[node_index];
	if( !tube.m_zestiVisited ) {
	  tube.m_zestiVisited = true;
	  tubeIndexQueue.inQueue(node_index);
	  // Fetch list of neigbors
	  std::vector<int> &neighboList = tube.m_neighbors;
	  for( size_t n = 0; n < neighboList.size(); ++n) {
	    int nid = neighboList[n];
	    size_t neigh_idx  = hitMap.Find(nid);
	    GridNode &Neigh_tube = Ingrid[neigh_idx];	  
	    if( (sttTubes->find(nid) != sttTubes->end() ) &&// Node in the current component
		(!Neigh_tube.m_zestiVisited) ) {// Not visited to avoid visiting multiple times
	      Neigh_tube.m_zestiVisited = true;
	      tubeIndexQueue.inQueue(neigh_idx);
	    }
	  }// END FOR neighborlist
	}// END IF not visited
      }// FOR STT_ITER, Node loop
      // Here, all participating nodes in the current component are
      // reordered by their layers.

      // We need to determine the total distance of the component with
      // respect to (0,0,0);
      std::vector<int> const &QnodeList = tubeIndexQueue.m_queueCont;
      float totaldist = 0.00;// Length of the component in XY-plane
      point3D prevPoint(0.00, 0.00, 0.00);// Temporary structure

      for(size_t n = 0; n < QnodeList.size(); ++n) {
	size_t Node_Index = QnodeList[n];
	GridNode const &node = Ingrid[Node_Index];
	// Doe we really need to modify virtuals as well?
	//if( node.m_type != GridNode::VIRTUAL_NODE ) {
	totaldist += sqrt( (node.m_xDet - prevPoint.m_x) * (node.m_xDet - prevPoint.m_x) +
			   (node.m_yDet - prevPoint.m_y) * (node.m_yDet - prevPoint.m_y) );
	//}
	// Update prevPoint
	prevPoint.m_x = node.m_xDet;
	prevPoint.m_y = node.m_yDet;
      }
      dbgtrkz("Total dist is %f",totaldist);
      // Process Queue
      point3D prevNodeCoord(0.00, 0.00, 0.00);
      float Node_distance = 0;

      while( !tubeIndexQueue.isEmpty() ) {
	size_t nodeIndex = tubeIndexQueue.popFront();
	GridNode &node = Ingrid[nodeIndex];

	// We do not need to modify virtuals
	Node_distance += sqrt( (node.m_xDet - prevNodeCoord.m_x) * (node.m_xDet - prevNodeCoord.m_x) +
			       (node.m_yDet - prevNodeCoord.m_y) * (node.m_yDet - prevNodeCoord.m_y) );
      
	if( node.m_type != GridNode::VIRTUAL_NODE) {
	  // node.m_z_Det = (Node_distance / totaldist) * lastVirtualZCoord;
	  node.m_z_Det = (Node_distance / totaldist) * maxVirtZVal;
	  dbgtrkz("New z coord dist of %d is %f", node.m_detID, node.m_z_Det);

	}
	prevNodeCoord.m_x = node.m_xDet;
	prevNodeCoord.m_y = node.m_yDet;
      }
    }
    else{
      dbgtrkz("Determine Z normalized by prev Dist.");

      /*
       * We need to modify the z coordinates of the non virtual nodes
       * nomalized by the distance from the center relatieve to the
       * previous node.*/
      float step_size = lastVirtualZCoord / static_cast<float>(sttTubes->size());

      float current_Zvalue = 0;
      for(stt_It = sttTubes->begin(); stt_It != sttTubes->end();++stt_It){
	int nodeID = *stt_It;
	size_t node_index = hitMap.Find(nodeID);
	GridNode &tube = Ingrid[node_index];
	// We do not need to modify virtuals
	if( tube.m_type != GridNode::VIRTUAL_NODE) {
	  current_Zvalue += step_size;
	  tube.m_z_Det = current_Zvalue;
	}
	else{
	  current_Zvalue = tube.m_z_Det + step_size;
	}
      }//END FOR STT_IT
    }
  }//END tracks loop
}
