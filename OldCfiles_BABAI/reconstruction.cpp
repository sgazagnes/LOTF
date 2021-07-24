
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
  int numElts = 2, firstid = -1, lastid = 0;
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
      if(firstid == -1) firstid = i;
      visitVirt = 1;
    }

    else if( visitVirt == 1 && node.m_type != GridNode::VIRTUAL_NODE){

      dbgtrkz("Next node %d on layer %d", nodeID, node.m_Layer);

      if(labs(node.m_Layer - prevLayer) == 1 ){
	dbgtrkz("Pushing on first layer");
	interNodesL1.push_back(&node);
	lastid = i;
      }
      else if (labs(node.m_Layer - prevLayer) == 2){
	dbgtrkz("Pushing on second layer");
	interNodesL2.push_back(&node);
	lastid = i;
      }
      else {
	//	if(labs(node.m_Layer - prevLayer) > 3)
	//	  anchorNext.m_detID = -1;
	//	else{
	dbgtrkz("We are too far in layers, missing virtual, creating one");
	if(interNodesL2.size() == 0){
	  error("Error");
	  break;
	}
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
	  dbgtrkz("Next node after virtual is %d and is para %d?", nodeID, nextPara);

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
	else{
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
      }



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
      if(xPts.size() > 10){
	zPts.erase(zPts.begin());
	yPts.erase(yPts.begin());
	zPts.pop_back();
      yPts.pop_back();
      }
       if(dir == -1){
      xPts.push_back(0.);
      yPts.push_back(0.);
      zPts.push_back(0.);
    } else {
      xPts.insert(xPts.begin(), 0.);
      yPts.insert(yPts.begin(), 0.);
      zPts.insert(zPts.begin(), 0.);
      }
    for(size_t i = 0; i < zPts.size(); i++)
      dbgtrkz("x %lf, y %lf, z %lf", xPts[i],  yPts[i],  zPts[i]);


    std::vector<double> p;
    std::vector<double> distxy;

    double newval;
    double distxx = 0;
    p.push_back(0.);
    distxy.push_back(0.);
    int start =1;

    for (size_t i = 0; i <  zPts.size()-1; i++){ 
       newval = p[i] + sqrt(pow(zPts[i+1]-zPts[i],2.) + pow(yPts[i+1]-yPts[i],2.));
      p.push_back(newval);
      //if(i >=start){
	//	distxx +=  sqrt(pow(xPts[i+1]-xPts[i],2.) + pow(yPts[i+1]-yPts[i],2.));
      //	distxy.push_back(newval);
      //  }
      //  dbgtrkz("%lf", newval);
    }

    if(dir == -1){
      distxx += sqrt(x.back()*x.back()+y.back()*y.back());
      for(int i =vect->size()-1; i >=lastid; i--)
	distxx +=  sqrt(pow(x[i-1]-x[i],2.) + pow(y[i-1]-y[i],2.));
    }else{
      distxx += sqrt(x[0]*x[0]+y[0]*y[0]);
      for(size_t i =0; i <lastid; i++)
	distxx +=  sqrt(pow(x[i+1]-x[i],2.) + pow(y[i+1]-y[i],2.));    
    }
    double *z_coef = polyFit(p, zPts, 1);
    double *y_coef = polyFit(p, yPts, 1);

    dbgtrkz("Fitted z track is %f + p%f", z_coef[0], z_coef[1]);
    dbgtrkz("Fitted y track is %f + p%f", y_coef[0], y_coef[1]);
    dbgtrkz("p: %lf, %lf", p[start], p[p.size()-1]);

    double pfirst = (y[0] - y_coef[0])/y_coef[1];
    double Node_distance =0;
    double z0 = dir == -1? zPts.back(): zPts[0];//z_coef[0]+z_coef[1]*p[start];
    double zf = dir == -1? z_coef[0]:z_coef[0]+z_coef[1]*p[p.size()-1];
    dbgtrkz("z: %lf, %lf", z0, zf);
    double y0 =y_coef[0]+y_coef[1]*p[start];
    double yf =y_coef[0]+y_coef[1]*p[p.size()-1];
    dbgtrkz("y: %lf, %lf", y0, yf);

    double totxydist =  distxx;//sqrt( pow(xPts.back()-xPts[start],2) + pow(yPts.back()-yPts[start],2));// * (node.m_xDet - prevNodeCoord.m_x) +
    //      (node.m_yDet - prevNodeCoord.m_y) * (node.m_yDet - prevNodeCoord.m_y) );
    //	Node_distance += sqrt( (node.m_xDet - prevNodeCoord.m_x) * (node.m_xDet - prevNodeCoord.m_x) +
    //			       (node.m_yDet - prevNodeCoord.m_y) * (node.m_yDet - prevNodeCoord.m_y) );
    dbgtrkz("Total fist %lf, slope %lf", totxydist, (zf-z0)/totxydist);
    bool change =0;

    /* for(int i = firstid; i >= 0; i--){
      Node_distance += sqrt(pow(x[i]-x[i+1],2)+pow(y[i]-y[i+1],2));//-xPts[start]-yPts[start]
      z[i] = z0 - (zf-z0)*Node_distance/totxydist;
      dbgtrkz("Dist cur %lf", Node_distance);
      dbgtrkz("New %f, %f, %f", x[i], y[i], z[i]);

      }*/

    /* Node_distance = 0;
    for(size_t i = 0; i <vect->size(); i++){
      if(i > firstid) Node_distance += sqrt(pow(x[i]-x[i-1],2)+pow(y[i]-y[i-1],2));//-xPts[start]-yPts[start]
      z[i] = z0 + (zf-z0)*Node_distance/totxydist;
      dbgtrkz("Dist cur %lf", Node_distance);
      dbgtrkz("New %f, %f, %f", x[i], y[i], z[i]);
      }*/
    if(dir == -1){
      for(int i = vect->size()-1; i >=0; i--){
	if(i < vect->size()-1) Node_distance += sqrt(pow(x[i]-x[i+1],2)+pow(y[i]-y[i+1],2));//-xPts[start]-yPts[start]
	else Node_distance += sqrt(pow(x[i],2)+pow(y[i],2));
	
	z[i] = z0 + (zf-z0)*Node_distance/totxydist;
	dbgtrkz("Dist cur %lf", Node_distance);
	dbgtrkz("New %f, %f, %f", x[i], y[i], z[i]);
      }
    }else{
      for(size_t i = 0; i <vect->size(); i++){
	if(i >= 0) Node_distance += sqrt(pow(x[i]-x[i-1],2)+pow(y[i]-y[i-1],2));//-xPts[start]-yPts[start]
	else Node_distance += sqrt(pow(x[i],2)+pow(y[i],2));
	
	z[i] = z0 + (zf-z0)*Node_distance/totxydist;
	dbgtrkz("Dist cur %lf", Node_distance);
	dbgtrkz("New %f, %f, %f", x[i], y[i], z[i]);
      }
  }
    //  Node_distance = 0;
    //z_coef[0] + z_coef[1]*(100.*i/vect->size());
      // dbgtrkz("%f",pval);
      //node.m_z_Det = z_coef[0] + z_coef[1]*pval;
    /*for(size_t i = 0; i < vect->size(); i++){
      //    int nodeID = vect->at(i);
      //    size_t node_index = hitMap.Find(nodeID);
      //    GridNode  &node = Ingrid[node_index];
      //dbgtrkz("%f, %f, %f", x[i], y[i], z[i]);
      if(i == firstid)
	change = 1;
      //  double pval = (y[i] - y_coef[0])/y_coef[1];
      //      if(i != 0)
      //	Node_distance += sqrt(pow(x[i],2)+pow(y[i],2));x
      //   else
      Node_distance = sqrt(pow(x[i]-xPts[start],2)+pow(y[i]-yPts[start],2));//-xPts[start]-yPts[start]
    dbgtrkz("Dist cur %lf", Node_distance);

    if(change) z[i] = z0 + (zf-z0)*Node_distance/totxydist;
    else z[i] = z0 - (zf-z0)*Node_distance/totxydist;//z_coef[0] + z_coef[1]*(100.*i/vect->size());
      // dbgtrkz("%f",pval);
      //node.m_z_Det = z_coef[0] + z_coef[1]*pval;
            dbgtrkz("New %f, %f, %f", x[i], y[i], z[i]);

	    }*/
    
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
