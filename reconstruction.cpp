
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
	      removeIdFromNeigh(neighNode, &prevNodes, curId);
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



	else if (sameLayer.size() > 0){   /* Same layer neighbors */

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
		removeIdFromNeigh(neighNode, &prevNodes, curId);		  
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
