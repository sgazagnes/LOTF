#include <iostream>
#include <set>
#include <sstream>
#include <algorithm>
#include  <cmath>
#include <fstream>
#include <vector>
#include <gsl/gsl_poly.h>
#include <numeric>
// Root headers
#include "TFile.h"
#include "TNtuple.h"
#include "TStopwatch.h"
#include "TH1.h"

//#include "reconstruction.h"
#include "gridNode.h"
#include "logc.h"
#include "simon_functions.h"
#include "path_queue.h"
#include "phconnect.h"


/* removeIdFromNeigh */

void removeIdFromNeigh(GridNode *neighNode, std::vector<int> *prevNodes, int curId){
  for(size_t i = 0; i < prevNodes->size(); i++){
    if(prevNodes->at(i) != neighNode->m_detID){ 		
      (neighNode->m_neighbors).erase(std::remove((neighNode->m_neighbors).begin(), (neighNode->m_neighbors).end(),prevNodes->at(i)), (neighNode->m_neighbors).end());
      if(i == 0)
	neighNode->parent = curId;
    }
  }
}


/* areAdjacent */


bool areAdjacent(CoordGrid &gr,  std::vector< GridNode > &Ingrid, std::vector<int> *v){
  size_t adjacent = 0;
  //std::vector< GridNode > &Ingrid = gr.m_grid;

  for (size_t i = 0; i < v->size(); i++){
    int neighId   = v->at(i);
    int neighIdx  = gr.Find(neighId);
    GridNode &neighNode = Ingrid[neighIdx];
    //prevNodes.push_back(neighId);
    if(neighNode.m_cm.size() >0) break;
    for (size_t j = i+1; j < v->size (); j++){
      // info("Are %d and %d connected?", neighId, v->at(j));
      if(neighId == v->at(j)) error("areAdjacent: should not be in vector");
      else if(neighNode.IsNeighboring(v->at(j)))
	adjacent++;
      // else
    }
  }
  
  if(adjacent >= v->size() -1)
    return true;
  else
    return false;
}

/* sortNeighbors */

bool sortNeighbors(CoordGrid &gr, GridNode *currentNode,  PathCandidate &cand, std::vector<int> &prev, std::vector<int> &same, std::vector<int> &next, std::vector<int> &virt, char *visited, int *dir){

  int curDir = *dir;
  std::vector< GridNode > &Ingrid  = gr.m_grid;
  size_t curLayer 	= currentNode->m_Layer;
  int curId   		= currentNode->m_detID;
  bool cond 		= true;

  for(size_t i = 0; i < currentNode->m_neighbors.size(); i++){
    int neighId 	= currentNode->m_neighbors[i];
    GridNode &neighNode = Ingrid[gr.Find(neighId)];
    
    if(cand.isInCandidate(neighId))
      continue;
    // dbgconnect("Node %d has one neig %d", curId, neighId);
    
    if(neighNode.m_type == GridNode::VIRTUAL_NODE){
      virt.push_back(neighId);
      continue;
    }	 

    if(neighNode.m_Layer > curLayer){
      //   dbgconnect("Node %d has one neigh up %d", curId, neighId);
      next.push_back(neighId);
      curDir |= UP;
    }	
    else if( neighNode.m_Layer < curLayer) {
      //dbgconnect("Node %d has one neigh down %d", curId, neighId);
      prev.push_back(neighId);
      curDir |= DOWN;
    }

    else {
      // dbgconnect("Node %d has one neigh on the same %d", curId, neighId);
      same.push_back(neighId);
      curDir |= SAME;
    }
  } 
   

  sort( next.begin(), next.end() );
  next.erase( unique( next.begin(), next.end() ), next.end() );
  sort( prev.begin(), prev.end() );
  prev.erase( unique( prev.begin(), prev.end() ), prev.end() );
  sort( same.begin(), same.end() );
  same.erase( unique( same.begin(), same.end() ), same.end() );
  	
  if( curDir > 6){
    //  dbgconnect("Neighbors of %d, a complex case to be solved later", curId);
    cond = false;
  }
  
  *dir = curDir;
  return cond;
}

void addNodesToCand (CoordGrid &gr, std::vector< GridNode > &Ingrid,  PathCandidate &cand, char *visited, std::vector<int> &v){

  for (size_t i = 0; i < v.size(); i++){
    int neighId   = v[i];
    GridNode *neighNode = &Ingrid[gr.Find(neighId)];
    cand.insertNewNode(gr, Ingrid, neighNode, cand.m_memberList->end());
    visited[neighId] = 1;
  }
}

void findEasyTracks (CoordGrid &gr, std::vector< GridNode > &Ingrid, std::vector < PathCandidate* > &tracklets,     std::vector<pair<int, unsigned short>> idToProcess, char *visited, int *candidateId){

  //std::vector< GridNode > &Ingrid = gr.m_grid;  

  for(unsigned int n = 0; n < idToProcess.size(); ++n) {

    //BETER REPLACE THIS WITH SETS
    std::vector<int> sameLayer; // For nodes on same layer
    std::vector<int> nextLayer; // For nodes on next layer
    std::vector<int> prevLayer; // For nodes on previous layer
    std::vector<int> nextVirt; // For virtual nodes
    std::vector<int> prevNodes; // For nodes processed during the previous step
    std::vector<int> *v; // To point to the correct vector depending where we are going
      
    int curDir 	       	= 0; // Current direction
    int nextDir 	= 0; // Next direction (when finding neighbors)
    int curId 		= idToProcess[n].first; // Id to process 
    GridNode *curNode   = &Ingrid[gr.Find(curId)]; // Current node 

    /* ++++++++++++++++++++++++++++++++++++++++++++++++++++ */   
    /* BEGIN WITH NODES WITH ONE NEIGHBOR OR IN LAYER LIMIT */     
    /* ++++++++++++++++++++++++++++++++++++++++++++++++++++ */

    // If we have to many neighbors, it might be hard to start from this point
    // Needs to be relaxed because some event are not processed because of this.
    
    if(curNode->m_neighbors.size() > 5) 
      break;

    // If node is not visited, and has either 1 neighbor, or 2 neighbors and is in the layer limit.
    if(!visited[curId] && (curNode->m_neighbors.size() == 1 || (curNode->m_LayerLimit == 1))){  // && curNode->m_neighbors.size() == 2
	
      int       n_neighbors = curNode->m_neighbors.size();
      size_t    curLayer    = curNode->m_Layer;
      int    	n_connected = 0;
      int 	neighId;
      GridNode *neighNode;
      bool      cond        = true; // To check whether we keep adding node or we stop

      PathCandidate *cand 	= new PathCandidate();// Create a new tracklet candidate

      // Find the next neighbors of this node
      cond = sortNeighbors(gr, curNode, *cand, prevLayer, sameLayer, nextLayer, nextVirt, visited,  &nextDir);

      // if no simple neighbors to add, we stop already 
      if(cond == false){ 
	delete(cand);
	continue;
      }

      // Create a new candidate path
      cand->m_id 		= (*candidateId)++;// tracklet id
      cand->m_tailNode 		= curId;
      visited[curId] 		= 1;
      cand->insertNewNode(gr, Ingrid, curNode, cand->m_memberList->end());
      prevNodes.push_back(curId); // Add to previously processed nodes

      //Force first encounter
      /*   if (nextDir & UP){           // we are going up 	      
	v = &nextLayer;
	n_neighbors -= sameLayer.size();
	sameLayer.clear();
      }  else if (nextDir & DOWN){ // we are going down
	v = &prevLayer;
	n_neighbors -= sameLayer.size();
	sameLayer.clear();
      } else                       // we are going on the same layer
      v = &sameLayer;*/

      
      int startIt = -1;
      // Stqrt the loop
      while(cond){
	startIt++;
       	//dbgconnect("With my buddy %d, we have %d neighbors, nextDir is %d", curId, n_neighbors, nextDir);

	//Choosing the next direction based on previous search of neaighbors    
	if (nextDir & UP){           // we are going up 	      
	  v = &nextLayer;
	  /* if(startIt < 2){
	    n_neighbors -= sameLayer.size();
	    sameLayer.clear();
	    }*/
	}  else if (nextDir & DOWN){ // we are going down
	  v = &prevLayer;
	  /* if(startIt < 2){
	    n_neighbors -= sameLayer.size();
	    sameLayer.clear();
	    }*/
	} else                       // we are going on the same layer
	  v = &sameLayer;

	
	if(n_neighbors == 1){ // Easy case, only one neighbor

	  if(visited[v->at(0)])
	    cond = false;
	  else{
	    // If there are some virtual nodes to add	  
	    if(nextVirt.size() > 0){ 	      
	      addNodesToCand(gr, Ingrid, *cand, visited, nextVirt); 
	      n_connected += nextVirt.size();	    
	      prevNodes.insert(prevNodes.end(),  nextVirt.begin(),  nextVirt.end());
	      nextVirt.clear();
	    }

	    //Adding the next neighbor
	    neighId    = v->at(0);
	    neighNode  = &Ingrid[gr.Find(neighId)];
	    cand->insertNewNode(gr, Ingrid, neighNode,cand->m_memberList->end());
	    visited[neighId] = 1;
	    n_connected++;

	    //Removing the previous nodes from the neighbor list of this node (to make things faster)
	    removeIdFromNeigh(neighNode, &prevNodes, curId);

	    //Setting the current node to the one added
	    curId       = neighId;
	    curNode = neighNode;
	    curLayer    = curNode->m_Layer;

	    //Cleaning lists, updating list of previous nodes, setting next direction
	    nextLayer.clear();  sameLayer.clear();   prevLayer.clear();    prevNodes.clear();
	    prevNodes.push_back(curId);
	    curDir = nextDir;
	    nextDir = 0;

	    //Finding next neighbors
	    cond = sortNeighbors(gr, curNode, *cand, prevLayer, sameLayer, nextLayer, nextVirt, visited, &nextDir);
	    n_neighbors = sameLayer.size() + prevLayer.size() + nextLayer.size();
	    //dbgconnect("%d nodes connected, %d found for next step (cond %d)\n", n_connected, n_neighbors, cond);
	    n_connected = 0;
	  }
	} // end if 1 Neighbor case



	else if (sameLayer.size() == 1 ){ // Same layer neighbor to handle... (NEED TO CHECK IS CAN BE MORE THAN 1)

	  // dbgconnect("Cur id%d, Check the same layer with the next node %d", curId,sameLayer[0], visited[sameLayer[0]]);
	  
	  if(visited[sameLayer[0]])  // If the same layer node has already been visited, we stop
	    cond = false;	      
	  
	  else { 

	    int candId 	 	= sameLayer[0];
	    GridNode *candNode  = &Ingrid[gr.Find(candId)];
	    curLayer            = candNode->m_Layer;

	    //CHECK if neighbors of the current candidate are all neighbors to one of the node on the next layer
	    for(size_t i = 0; i < candNode->m_neighbors.size(); i++){
	      neighId = candNode->m_neighbors[i];	      
	      neighNode = &Ingrid[gr.Find(neighId)];
	      if (neighId == curId || neighNode->m_type == GridNode::VIRTUAL_NODE)
		continue;
	      
	      int haveNeigh = 0; // Check variable to know if we found a neighbor		
	      for (size_t j = 0; j < v->size(); j++){
		//dbgconnect("Test adjacency between neigh ID %d and v at j %d", neighId, v->at(j));

		// If the node in the next layer is a neighbor, we can stop
		if(v->at(j) == neighId || (neighNode->IsNeighboring(v->at(j)))){
		  haveNeigh = 1; 
		  break;
		}

		// if the node in the next layer is not a neighbor,
		// that might be because of a peculier geometrical configuration
		// We can test this by checking if all this nod eneighbors are either visited or
		// in the next layer vector of nodes
		else{
		  GridNode &mynode = Ingrid[gr.Find(v->at(j))];
		  bool addthisnode = true;
		  for(size_t k = 0; k < mynode.m_neighbors.size(); k++){
		    int oth = mynode.m_neighbors[k];

		    if(oth!= candId && !cand->isInCandidate(oth) && std::find(v->begin(), v->end(), oth) == v->end() && std::find( candNode->m_neighbors.begin(),  candNode->m_neighbors.end(), oth) ==  candNode->m_neighbors.end()  ){
		      //dbgconnect("This node %d has other neighbor %d not visited", v->at(j), oth);
		      addthisnode = false;
		      break; //The node did not fullfilled the condition, this is a hard case, we should break here
		    }
		  }
		  if(addthisnode){
		    // The node fulfilled the condition, let's add it to the list and keep looking for more
		    if(nextVirt.size() > 0){ // Add the virtual if it exists
		      for(size_t i = 0; i < nextVirt.size();){
			GridNode &virt = Ingrid[gr.Find(nextVirt[i])];
			if(virt.IsNeighboring(v->at(j))){
			  cand->insertNewNode(gr, Ingrid, &virt,cand->m_memberList->end());
			  visited[virt.m_detID] = 1;
			  n_connected++;
			  nextVirt.erase(nextVirt.begin() + j);
			}
		      }
		    }
		    //	    dbgconnect("We are adding the node %d", v->at(j));

		    cand->insertNewNode(gr, Ingrid, &mynode,cand->m_memberList->end());
		    visited[mynode.m_detID] = 1;
		    n_connected++;
		    removeIdFromNeigh(&mynode, &prevNodes, curId);
		    prevNodes.push_back(mynode.m_detID);
		    v->erase(v->begin()+j);
		    j--;
		    haveNeigh = 1;
		    break;
		  } 
		}
	      }

	      // If we did not find that all nodes of the candidate were connected
	      // to other neighbors, we should stop here
	      if(haveNeigh == 0){ 
		cond = false;
	       	break;
	      }	else if ( v->size() == 0)// Just to solve a case where the next list is empty
		break;
	    } // End FOR neighbors of candidates

	    // If we did not find a good reason to stop, then we add the candidate	    
	    if (cond == true){ 		
	      // dbgconnect("All neighbors of candidate look good, let's insert it!");

	      cand->insertNewNode(gr, Ingrid, candNode,cand->m_memberList->end());
	      visited[candId] = 1;
	      n_connected++;

	      removeIdFromNeigh(candNode, &prevNodes, curId);
	      
	      prevNodes.push_back(candId);
	      curId       = candId;
	      curNode     = candNode;
	      curDir 	  = SAME;
	      nextDir 	  = 0;
	      sameLayer.clear();

	      cond        = sortNeighbors(gr, curNode, *cand, prevLayer, sameLayer, nextLayer,
					  nextVirt, visited, &nextDir);

	      n_neighbors = sameLayer.size()+prevLayer.size()+nextLayer.size();
	      // dbgconnect("%d nodes connected,  %d found for next step (cond %d) \n", n_connected, n_neighbors, cond);
	      n_connected = 0;		
	    } 

	    // If we had to stop
	    else {
	      // dbgconnect("Neighbors not connected... \n", n_connected, n_neighbors);
	    } 
	  } // ELSE the same layer not was visited already
	} // END of IF SAME LAYER == 1



	// IF we have less than 5 neighbors on the next/prev layer (5 is empirical limit set by me...)
	
	else if (n_neighbors < 5) {  	  
	  // dbgconnect("We have %d neighbors in the direction %d", n_neighbors, nextDir);

	   // The idea here is not check if all the nodes in the next direction are neighbors to each other
	  // We first test whether there are not, because then we can check whether one of this node is
	  // an obvious choice or not

	       
	  if(!areAdjacent(gr, Ingrid, v)){
	    // dbgconnect("Nodes are not adjacent but we can check wheter one of this node only has a single neighbor");
	    std::vector <char> nother;

	    // I need to comment properly this otherwise nobody will understand what is happening
	    for(size_t i = 0; i < v->size(); i++){
	      GridNode &mynode = Ingrid[gr.Find(v->at(i))];
	      char n = 0;
	      for(size_t j = 0; j < mynode.m_neighbors.size(); j++){
		GridNode &thenode = Ingrid[gr.Find(mynode.m_neighbors[j])];
		if(thenode.m_type == GridNode::VIRTUAL_NODE)
		  continue;
		if(thenode.m_Layer == curLayer && !(cand->isInCandidate(thenode.m_detID))){
		  n = 1;
		  //	  dbgconnect("This node %d has a neighbord %d that does not belong to the track",
		  //    v->at(i), thenode.m_detID);
		}
	      }
	      nother.push_back(n);
	    }
	    
	    char sum = std::accumulate(nother.begin(), nother.end(), decltype(nother)::value_type(0));
	    if(sum && sum != (char) v->size()){
	      //dbgconnect("We found at least one node with only previous nodes as neighbors");
	      std::vector <int> toAdd;
	      for(size_t i =0; i < nother.size(); ){
		if(nother[i] == 0){
		  i++;
		} else {
		  v->erase(v->begin()+i);
		  nother.erase(nother.begin()+i);
		}
	      }

	      // Delete virtuals that are not neighbors to these nodes
	      for(size_t i = 0; i < nextVirt.size();){
		GridNode &virt = Ingrid[gr.Find(nextVirt[i])];
		int neighId1 = virt.m_neighbors[0];
		int neighId2 = virt.m_neighbors[1];
		if(std::find(v->begin(), v->end(), neighId2) != v->end() ||
		  std::find(v->begin(), v->end(), neighId1) != v->end()){
		  i++;
		}
		else
		  nextVirt.erase(nextVirt.begin()+i);
	      }
	    }
	    else {
	      //  dbgconnect("Too many neighbors and too complicated");
	      cond = false;
	    }
	  } else {
	    //dbgconnect("All Adjacent, adding %lu nodes to the CM (and %lu virtuals)", v->size(), nextVirt.size());
	  }
	    // if(areAdjacent(gr, Ingrid, v)){ // All neighbors are adjacent ??
	  for(size_t i = 0; i < v->size(); i++){
	    if(visited[v->at(i)])
	      cond = false;
	  }     
	  if(cond){


	    
	    if(nextVirt.size() > 0){ // taking care of the virtual nodes
	      if(nextVirt.size()> 3){
		//Cleaning virtual nodes
		//	error("Changing virt lsit");
		int sizeT = prevNodes.size();
		prevNodes.insert(prevNodes.end(),  nextVirt.begin(),  nextVirt.end());

		//	for(size_t f = 0; f < nextVirt.size(); f++)
		  //	  error("Virt %d", nextVirt[f]);
		nextVirt.clear();
		int goodV[2] = {-1,-1};
		for(size_t f = 0; f < sizeT; f++){
		  float mindist[2] = {1000,1000};
		  GridNode &prevC = Ingrid[gr.Find(prevNodes[f])];
		  for(size_t g = 0; g < v->size(); g++){
		    GridNode &nextC = Ingrid[gr.Find(v->at(g))];
		    float curdis = sqrt(pow(prevC.m_x-nextC.m_x,2) + pow(prevC.m_y-nextC.m_y,2));
		    if(curdis < mindist[0]){
		      goodV[1]= goodV[0];
		      mindist[1]= mindist[0];
		      for(size_t h = 0; h< nextC.m_neighbors.size(); h++){
			GridNode &neighC = Ingrid[gr.Find(nextC.m_neighbors[nextC.m_neighbors.size()- 1 -h])];
			if(neighC.m_type != GridNode::VIRTUAL_NODE)
			  continue;
			if(prevC.IsNeighboring(neighC.m_detID))
			  goodV[0] = neighC.m_detID;
		      }
		      mindist[0]=curdis;
		    } else if( curdis <mindist[1]){
		      for(size_t h = 0; h< nextC.m_neighbors.size(); h++){
			GridNode &neighC = Ingrid[gr.Find(nextC.m_neighbors[nextC.m_neighbors.size()- 1 -h])];
			if(neighC.m_type !=  GridNode::VIRTUAL_NODE)
			  continue;
			if(prevC.IsNeighboring(neighC.m_detID))
			  goodV[1] = neighC.m_detID;
		      }
		      //   goodV[1] = v->at(g);
		      mindist[1]=curdis;
		    }
		  }
		  if (goodV[0] != -1)
		    nextVirt.push_back(goodV[0]);
		  //  if (goodV[1] != -1)
		  //  nextVirt.push_back(goodV[1]);
		}
		//	for(size_t f = 0; f < nextVirt.size(); f++)
		//	  error("New Virt %d", nextVirt[f]);
	      }
		  
		  
	      addNodesToCand(gr, Ingrid, *cand, visited, nextVirt); 
	      n_connected += nextVirt.size();	    	  		
	      nextVirt.clear();
	    }
	      
	    std::vector<int> lookneigh(*v);

	    for (size_t i = 0; i < v->size(); i++){ // add the next nodes
	      neighId   = v->at(i);
	      neighNode = &Ingrid[gr.Find(neighId)];
	      cand-> insertNewNode(gr, Ingrid, neighNode, cand->m_memberList->end());
	      visited[neighId] = 1;
	      n_connected++;
	      removeIdFromNeigh(neighNode, &prevNodes, curId);		  
	    }
	      
	    nextLayer.clear(); sameLayer.clear();    prevLayer.clear();    prevNodes.clear();
	      
	    curDir = nextDir;
	    nextDir = 0;
	    
	    curId       = lookneigh[0];
	    curNode     = &Ingrid[gr.Find(curId)];
	    curLayer    = curNode->m_Layer;
	      
	    n_neighbors = 0;
	      
	    for(size_t i = 0; i < lookneigh.size(); i++){ // Looking for new neighbors
	      int id         = lookneigh[i];
      	      GridNode *node = &Ingrid[ gr.Find(id)];
	      prevNodes.push_back(id);
	      cond = sortNeighbors(gr, node, *cand, prevLayer, sameLayer, nextLayer, nextVirt, visited,  &nextDir);
	    }

	    n_neighbors = sameLayer.size()+prevLayer.size()+nextLayer.size();
	    // dbgconnect("%d nodes connected, %d found for next step (cond %d)\n", n_connected, n_neighbors, cond);
	    n_connected = 0;
	     
	  } 	      
	} // END IF neighbors size > 0

	// We have no more neighbors or too many
	if (n_neighbors == 0 || n_neighbors >= 5 || !cond){	    
	  if(n_neighbors == 0) {
	    int firstId         = cand->m_tailNode;
	    GridNode &firstNode = Ingrid[ gr.Find(firstId)];
	    int lastId          = cand->m_headNode;
	    GridNode &lastNode  = Ingrid[gr.Find(lastId)];
	    dbgconnect("No more neighbors in sight, checking if could be finished already ?");
	      
	    if((firstNode.m_LayerLimit == 1 && lastNode.m_LayerLimit == 1)){		 
	      dbgconnect("track goes through all layers or makes a loop, likily finished");		 
	      cand->m_finished = FINISHED;		 
	    } else if(labs(curNode->m_SectorLimit) > 0 || cand->m_isOnSectorLimit){
	      dbgconnect("Track is on sector limit, might have a connection somewhere else");
	      cand->m_finished = ONGOING;
	      cand->m_isOnSectorLimit= true;		 
	    } else {		 
	      dbgconnect("Candidate has no more neighbors, but doesn't seem finished (First %d, last %d)", firstId, lastId);
	      cand->m_finished = ONGOING;		 
	    }	      
	  } // end if n_neighbors == 0
	    
	  cond = false;
	    
	}

	if(nextDir == 5){
	  //dbgconnect("We could go either up or down in next round, let's check the previous direction");
	  if(curDir&UP){
	    //dbgconnect("Let's go UP");
	    nextDir = 4;
	    n_neighbors = nextLayer.size();
	    for(size_t i = 0; i < nextVirt.size();){
	      GridNode &virt = Ingrid[gr.Find(nextVirt[i])];
	      int neighId1 = virt.m_neighbors[0];
	      int neighId2 = virt.m_neighbors[1];
	      if((cand->isInCandidate(neighId1) && std::find(nextLayer.begin(), nextLayer.end(), neighId2) != nextLayer.end()) || (cand->isInCandidate(neighId2) && std::find(nextLayer.begin(), nextLayer.end(), neighId1) != nextLayer.end())){
		i++;
	      }
	      else
		nextVirt.erase(nextVirt.begin()+i);
	    }


	    // CLEAN VIRTUALs
	  } else if(curDir&DOWN){
	    //dbgconnect("Let's go DOWN");
	    nextDir = 1;
	    n_neighbors = prevLayer.size();
	    for(size_t i = 0; i < nextVirt.size();){
	      GridNode &virt = Ingrid[gr.Find(nextVirt[i])];
	      int neighId1 = virt.m_neighbors[0];
	      int neighId2 = virt.m_neighbors[1];
	      if((cand->isInCandidate(neighId1) && std::find(prevLayer.begin(), prevLayer.end(), neighId2) != prevLayer.end()) || (cand->isInCandidate(neighId2) && std::find(prevLayer.begin(), prevLayer.end(), neighId1) != prevLayer.end())){
		i++;
	      }
	      else
		nextVirt.erase(nextVirt.begin()+i);
	    }
	  } else {
	    //dbgconnect("Can't decide, stop");
	    cond = false;
	  }
	}
	  	  
	if(cond == false && n_neighbors != 0){ // This track is finished, but let's push neighbors as we need to fit in the next step

	  // dbgconnect("ADDING lot of nodes");
	  for (size_t i = 0; i < sameLayer.size(); i++)
	    cand->m_headNeigh.push_back(sameLayer[i]);
	      
	  for (size_t i = 0; i < nextLayer.size(); i++)
	    cand->m_headNeigh.push_back(nextLayer[i]);

	  for (size_t i = 0; i < prevLayer.size(); i++)
	    cand->m_headNeigh.push_back(prevLayer[i]);

	  cand->m_headNeigh.insert((cand->m_headNeigh).begin(),  (nextVirt).begin(),  (nextVirt).end());	
	  // resetLists(visited, prevLayer, sameLayer, nextLayer);
	  //for(size_t i = 0; i < cand->m_headNeigh.size(); i++)
	  // dbgconnect("Node %d",cand->m_headNeigh[i]);

	  cand->m_finished = ONGOING;		 
	  break;
	}
       		
      }// END OF COND WHILE LOOP (MEANS WE STOP THE CONNECT PHASE FOR THIS TRACK

	
      if(cand->m_length > 2){
	dbgconnect("Pushing cm %d with length %d, tail node %d, head node %d, min layer %d, max layer %d, IsOnSectorLimit %d, status  %d. ", cand->m_id, cand->m_length, cand->m_tailNode, cand->m_headNode,cand->m_minLayer, cand->m_maxLayer, cand->m_isOnSectorLimit, cand->m_finished);
	tracklets.push_back(cand);
      } else {
	error("Not a good cm %d, has length < 2",cand->m_headNode);
	for(size_t i = 0; i < (cand->m_memberList)->size(); i++){   
	  visited[(cand->m_memberList)->at(i)] = 0;
	  GridNode &toDel = Ingrid[gr.Find((cand->m_memberList)->at(i))];
	  toDel.m_cm.clear();

	}
	  
	delete cand;
	(*candidateId)--;
      } // Else Candidate large enough
	
    } // for node with a single neighbor
       	     
  } // For active nodes

}

