
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


void findEasyTracks (CoordGrid &gr, std::vector < PathCandidate* > &tracklets, std::vector<int> &activeId, char *visited, int *candidateId){

  std::vector< GridNode > &Ingrid = gr.m_grid;  

  for(unsigned int n = 0; n < activeId.size(); ++n) {
      
    std::vector<int> sameLayer;
    std::vector<int> nextLayer;
    std::vector<int> prevLayer;
    std::vector<int> nextVirt;
    std::vector<int> prevNodes;
    std::vector<int> *v;
      
    int dir 			= 0;
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

      info("Starting node %d has %d neighbors", curId, n_neighbors);

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
	    for (int i = 0; i < nextVirt.size(); i++){
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
	    
	    for(int i = 0; i < candNode->m_neighbors.size(); i++){ //CHECK if neighbors on next layer are adjacent 
	      neighId = candNode->m_neighbors[i];	      
	      neighIdx  = gr.Find(neighId);
	      neighNode = &Ingrid[neighIdx];
		
	      if (neighId == curId || neighNode->m_type == GridNode::VIRTUAL_NODE) continue;
			      
	      int haveNeigh = 0;
		
	      for (int j = 0; j < v->size(); j++){		  
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
	    error("WHAT IS THE DIRECTION NOW?");
	           
	  if(areAdjacent(gr, v)){ // All neighbors are adjacent ??
	      
	    //nfo("Adding %d nodes to the CM", v->size());
	      
	    if(nextVirt.size() > 0){ // taking care of the virtual nodes
	      for (int i = 0; i < nextVirt.size(); i++){
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

	    for (int i = 0; i < v->size(); i++){ // add the next nodes
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
	      
	    for(int i = 0; i < lookneigh.size(); i++){ // Looking for new neighbors
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
	      
	    info("No more neighbors in sight, are we finished here ?");
	      
	    if((cand->m_minLayer == 0 && cand->m_maxLayer > 21) || (firstNode->m_LayerLimit == 1 && lastNode->m_LayerLimit == 1)){		 
	      info("track goes through all layers or makes a loop, likily finished");		 
	      cand->m_finished = 3;		 
	    } else if(abs(currentNode->m_SectorLimit) > 0 || cand->m_isOnSectorLimit){
	      info("Track is on sector limit, might have a connection somewhere else");
	      cand->m_finished = 2;
	      cand->m_isOnSectorLimit= true;		 
	    } else {		 
	      info("Candidate has no more neighbors, but doesn't seem finished");
	      cand->m_finished = 2;		 
	    }
	      
	  } // end if n_neighbors == 0

	    
	  cond = false;
	    
	}

	  	  
	if(cond == false){ // This track is finished, but let's push neighbors as we need to fit in the next step

	  for (int i = 0; i < sameLayer.size(); i++)
	    cand->m_headNeigh.push_back(sameLayer[i]);
	      
	  for (int i = 0; i < nextLayer.size(); i++)
	    cand->m_headNeigh.push_back(nextLayer[i]);

	  for (int i = 0; i < prevLayer.size(); i++)
	    cand->m_headNeigh.push_back(prevLayer[i]);

	  cand->m_headNeigh.insert((cand->m_headNeigh).begin(),  (nextVirt).begin(),  (nextVirt).end());	
	  resetLists(visited, prevLayer, sameLayer, nextLayer);
	}

      }// WHILE COND

	
      if(cand->m_length > 2){
	info("Pushing cm %d: \n \t length is %d, \n \t tail node %d \n \t head node %d \n \t Min layer %d, \n \t Max layer %d, \n \tIsOnSectorLimit %d, \n \t Finished ? %d. ", cand->m_id, cand->m_length, cand->m_tailNode, cand->m_headNode,cand->m_minLayer, cand->m_maxLayer, cand->m_isOnSectorLimit, cand->m_finished);
	tracklets.push_back(cand);
      } else {
	info("Not a good cm %d",cand->m_headNode);
	for(int i = 0; i < (cand->m_memberList)->size(); i++){   
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

void fitZCoordinates(CoordGrid &hitMap, PathCandidate *trk)
{
 
 
  info("Determining the z coordinates for the track", trk->m_id);

  // MVD list of available elements in the graph.
  // std::vector< GridNode > &MVDPoints = hitMap.m_MVD_grid;
  //
  // STT lists of elements available in the graph.
  std::vector< GridNode > &Ingrid = hitMap.m_grid;
  // Local variables
  std::set<int>::iterator stt_It;
  // Track (Connected components) )loop

  // Current track object  std::vector<double> z =  std::vector<double>( trk->m_z );

  std::vector<double> &x =  trk->m_x;
  std::vector<double> &y =  trk->m_y;
  std::vector<double> &z =  trk->m_z;
  std::vector<int> const *vect = trk->m_memberList;

  // Fnd the most outer virtual node or the max Z
  // size_t maxLayer = std::numeric_limits<size_t>::min();
  // float maxVirtZVal = std::numeric_limits<float>::min();
  //float  lastVirtualZCoord = 0.00;
  for( int i = 0; i < z.size(); i++){
    printf("%d, %lf \t", vect->at(i), z[i]);
  }
  printf("\n");

  std::vector<GridNode*>  virt;
  std::vector<GridNode*>  prevLayer;
  std::vector<GridNode*>  nextLayer;
  std::vector<GridNode*>  interLayer1;
  std::vector<GridNode*>  interLayer2;
  int detect = 0;
  int idxFirstVirt = 0, idxFirstSkewed;
  double xvirtfirst = 0, yvirtfirst = 0, zvirtfirst = 0;
  double xvirtsecond = 0, yvirtsecond = 0, zvirtsecond = 0;
  double lastx=0, lasty =0;
  GridNode Dummy_coord;

  for(int i = 0; i < vect->size(); i++){
    int nodeID = vect->at(i);
    size_t node_index = hitMap.Find(nodeID);
    GridNode  *node = &Ingrid[node_index];

    if( detect < 2 && node->m_type == GridNode::VIRTUAL_NODE) {
      virt.push_back(node);
      int curId = node->m_neighbors[0];
      size_t curIdx = hitMap.Find(curId);
      GridNode  *curNode = &Ingrid[curIdx];
      lastx = curNode->m_x;
      lasty = curNode->m_y;
      if(curNode->m_type == GridNode::STT_TYPE_SKEW){
	curId = node->m_neighbors[1];
	curIdx = hitMap.Find(curId);
	curNode = &Ingrid[curIdx];
	lastx = curNode->m_x;
	lasty = curNode->m_y;
      }
      if(detect == 0){
	idxFirstVirt = i;
	//	prevLayer.push_back(prevnode);
	//	nextLayer.push_back(nextnode);
	detect = 1;
      } 
	
    }	
    else if(detect == 1){
      
      debug("Number of virtual nodes found %d", virt.size());


      for(int j = 0; j< virt.size(); j++){
	xvirtfirst += virt[j]->m_x;
	yvirtfirst += virt[j]->m_y;
	zvirtfirst += virt[j]->m_z_Det;
      }
      xvirtfirst /= virt.size();
      yvirtfirst /= virt.size();
      zvirtfirst /= virt.size();
      //xvirtfirst = (x[i] + lastx)/2.0;
      //yvirtfirst = (y[i] + lasty)/2.0;

      debug("average coord of virt layer %lf %lf", xvirtfirst, yvirtfirst);
      for(int j = idxFirstVirt; j < idxFirstVirt + virt.size(); j++){
	x[j] = xvirtfirst;
	y[j] = yvirtfirst;
	z[j] = zvirtfirst;
      }
      
      virt.clear();
      // interLayer1.insert(insertLayer1.end(),  nextLayer->begin(),  nextLayer->end());
      //nextLayer.clear();
      interLayer1.push_back(node);
      idxFirstSkewed = i;
      debug("Adding %d on interlayer 1",node->m_detID);
      detect = 2;
      
    } else if(detect == 2 &&  node->m_type != GridNode::VIRTUAL_NODE){
      if (node->m_Layer == interLayer1[0]->m_Layer && !(std::find(std::begin(interLayer1), std::end(interLayer1), node) != std::end(interLayer1))) {
	interLayer1.push_back(node);
	debug("Adding %d on interlayer 1",node->m_detID);
      } else if (labs(node->m_Layer - interLayer1[0]->m_Layer)<2) {
	interLayer2.push_back(node);
	debug("Adding %d on interlayer 2",node->m_detID);
      } else {
	debug("We are missing a virtual for %d",node->m_detID);
	GridNode tubeA = *(interLayer2.back());
	GridNode tubeB = *(node);
	info("%d %d", tubeA.m_detID , tubeB.m_detID); 

	IntersectionPoint(hitMap, tubeA, tubeB, Dummy_coord);
	Dummy_coord.m_detID      = 100000;
	info("Dummy_coord %d",Dummy_coord.m_detID ); 
	std::pair<float, float> r_Theta;
	float theta_deg = Cartesian_To_Polar(Dummy_coord.m_xDet, Dummy_coord.m_yDet, r_Theta);
	Dummy_coord.m_r = r_Theta.first;
	Dummy_coord.m_thetaDeg = theta_deg;
	xvirtsecond = Dummy_coord.m_x;
	yvirtsecond = Dummy_coord.m_y;
	zvirtsecond = Dummy_coord.m_z;
	idxFirstVirt = -1;
	detect=3;
	i--;
	virt.push_back(&Dummy_coord);
	//	interLayer2.push_back(node);// to remove
      }
    } else if(detect == 2  && node->m_type == GridNode::VIRTUAL_NODE){
      idxFirstVirt = i;
      virt.push_back(node);
      detect = 3;
      lastx = x[i-1];
      lasty = y[i-1];
    } else if(detect == 3 && node->m_type == GridNode::VIRTUAL_NODE){
      virt.push_back(node);
    } else if(detect == 3 && node->m_type != GridNode::VIRTUAL_NODE){

      
      debug("We have the intermediate layer. 1 %d, 2 %d, virtual further  %d", interLayer1.size(), interLayer2.size(), virt.size());

      if(idxFirstVirt != -1){
	for(int j = 0; j< virt.size(); j++){
	  xvirtsecond += virt[j]->m_x;
	  yvirtsecond += virt[j]->m_y;
	  zvirtsecond += virt[j]->m_z_Det;
	}
	xvirtsecond /= virt.size();
	yvirtsecond /= virt.size();
	zvirtsecond /= virt.size();
	//xvirtsecond = (x[i] + lastx)/2.0;
	//yvirtsecond = (y[i] + lasty)/2.0;

	for(int j = idxFirstVirt; j < idxFirstVirt + virt.size(); j++){
	  x[j] = xvirtsecond;
	  y[j] = yvirtsecond;
	  z[j] = zvirtsecond;
	}
      }
      virt.clear();


      debug("Virt point 1 %lf %lf %lf \n virt point 2 %lf, %lf, %lf", xvirtfirst, yvirtfirst, zvirtfirst,xvirtsecond, yvirtsecond, zvirtsecond);
	
	if(interLayer1.size() <= 3 ){

	  double layerX = interLayer1[0]->m_x, layerY= interLayer1[0]->m_y, layerZ =interLayer1[0]->m_z ;
	  TVector3 dir = interLayer1[0]->m_WireDirection; 
	  double hl = interLayer1[0]->m_halfLength;
	  if(interLayer1.size() > 1){
	    for(int j = 1; j < interLayer1.size(); j++){
	      layerX += interLayer1[j]->m_x;
	      layerY += interLayer1[j]->m_y;
	      layerZ += interLayer1[j]->m_z;
	    }
	    layerX /= interLayer1.size();
	    layerY /= interLayer1.size();
	    layerZ /= interLayer1.size();
	  }



	  // debug(" Start x1 %lf, end x1 %lf, start y1 %lf, end y 1 %lf", layerX - hl*dir[0], layerX + hl*dir[0],  layerY - hl*dir[1],  layerY + hl*dir[1]);
	  double coef = IntersectionXY(layerX - hl*dir[0], layerX + hl*dir[0],  layerY - hl*dir[1],  layerY + hl*dir[1], xvirtfirst,xvirtsecond ,yvirtfirst, yvirtsecond);
	  double newX =  (1.0-coef) * (layerX - hl*dir[0]) + coef * (layerX + hl*dir[0]);
	  double newY =  (1.0-coef) * (layerY - hl*dir[1]) + coef * (layerY + hl*dir[1]);
	  double newZ =  (1.0-coef) * (layerZ - hl*dir[2]) + coef * (layerZ + hl*dir[2]);
      
	  info("Coef %lf, Intersection point on layer 1 is %lf, %lf, %lf", coef, newX, newY, newZ);

	  double diffX = newX - layerX;
	  double diffY = newY - layerY;
	  double diffZ = newZ - layerZ;
      
	  for(int j = idxFirstSkewed, k =0; j < idxFirstSkewed + interLayer1.size(); j++, k++){
	    x[j] = interLayer1[k]->m_x + diffX;
	    y[j] = interLayer1[k]->m_y + diffY;
	    z[j] = newZ;
	    interLayer1[k]->m_z_Det = newZ;
	    debug("Node %d has new coord %lf, %lf %lf", vect->at(j), x[j], y[j], z[j]);
	  }
      

	  if(interLayer2.size()> 0){
      
	  layerX = interLayer2[0]->m_x;
	  layerY = interLayer2[0]->m_y;
	  layerZ = interLayer2[0]->m_z;
	  dir    = interLayer2[0]->m_WireDirection; 
	  hl     = interLayer2[0]->m_halfLength;
	  if(interLayer2.size()> 1){
	    for(int j = 1; j < interLayer2.size(); j++){
	      layerX += interLayer2[j]->m_x;
	      layerY += interLayer2[j]->m_y;
	      layerZ += interLayer2[j]->m_z;
	    }
	    layerX /= interLayer2.size();
	    layerY /= interLayer2.size();
	    layerZ /= interLayer2.size();
	  }
	  coef = IntersectionXY(layerX - hl*dir[0], layerX + hl*dir[0],  layerY - hl*dir[1],  layerY + hl*dir[1], xvirtfirst,xvirtsecond ,yvirtfirst, yvirtsecond);
	  newX =  (1.0-coef) * (layerX - hl*dir[0]) + coef * (layerX + hl*dir[0]);
	  newY =  (1.0-coef) * (layerY - hl*dir[1]) + coef * (layerY + hl*dir[1]);
	  newZ =  (1.0-coef) * (layerZ - hl*dir[2]) + coef * (layerZ + hl*dir[2]);
      
	  info("Coef %lf, Intersection point on layer 2 is %lf, %lf, %lf", coef, newX, newY, newZ);

	  diffX = newX - layerX;
	  diffY = newY - layerY;
	  diffZ = newZ - layerZ;
      
	  for(int j = idxFirstSkewed+interLayer1.size(), k =0; j < idxFirstSkewed +interLayer1.size()+ interLayer2.size(); j++, k++){
	    x[j] = interLayer2[k]->m_x + diffX;
	    y[j] = interLayer2[k]->m_y + diffY;
	    z[j] = newZ;
	    interLayer2[k]->m_z_Det = newZ;

	    debug("Node %d has new coord %lf, %lf %lf", vect->at(j), x[j], y[j], z[j]);
		  
	  }
	  }
	  
      } else {
	info("Ambiguous, too many nodes ?");
	for(int j = idxFirstSkewed, k =0; j < idxFirstSkewed + interLayer1.size(); j++, k++){
	  x[j] = interLayer1[k]->m_x;
	  y[j] = interLayer1[k]->m_y;
	  z[j] = -1;
	  debug("Node %d has new coord %lf, %lf %lf", vect->at(j), x[j], y[j], z[j]);
	}
      
      }
      
      interLayer1.clear();
      interLayer2.clear();
      xvirtfirst = xvirtsecond;
      yvirtfirst = yvirtsecond;
      xvirtsecond = yvirtsecond = 0;

      if(node->m_type != GridNode::STT_TYPE_SKEW){
	info("Need to interpolate now\n\n");
	break;
      }
      interLayer1.push_back(node);
      debug("Adding %d on interlayer 1",node->m_detID);
      idxFirstSkewed = i;
      detect = 2;
    }
      
    
    // END IF virtual
  // Determine the largest Z value (Virtual node)
  }

  /* interpolation */
  /*  point3D prevPoint(0.00, 0.00, 0.00);// Temporary structure

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
    //std::cout << " Total dist is " << totaldist << std::endl;
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
    }
    prevNodeCoord.m_x = node.m_xDet;
    prevNodeCoord.m_y = node.m_yDet;
    }*/

  /*
   * We need to modify the z coordinates of the non virtual nodes
   * nomalized by the distance from the center relatieve to the
   * previous node.*/
  /*float step_size = lastVirtualZCoord / static_cast<float>(sttTubes->size());

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
  /*  


      for(stt_It = sttTubes->begin(); stt_It != sttTubes->end();++stt_It){
      int nodeID = *stt_It;
      size_t node_index = hitMap.Find(nodeID);
      GridNode const &node = Ingrid[node_index];
      if( (node.m_type == GridNode::VIRTUAL_NODE) &&
      (node.m_Layer > maxLayer) ) {
      maxLayer = node.m_Layer;
      lastVirtualZCoord = node.m_z_Det;
      }// END IF virtual
      // Determine the largest Z value (Virtual node)
      if( (node.m_type == GridNode::VIRTUAL_NODE) &&
      (node.m_z_Det > maxVirtZVal) ) {
      maxVirtZVal = node.m_z_Det;
      }// END IF virtual

      }// END FOR (node loop)
  */
  // Found the most outer virtual node
  /* Delta_z is equal to lastVirtualZCoord. Because we believe that
     the tracks start at (0,0,0)*/

  // We need the correct order of tubes by their layer numbers. The
  // ids of virtual tubes does not correspond with their places in
  // the graph.
  /*

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
   
    }//END FOR STT_IT*/
}
//____
