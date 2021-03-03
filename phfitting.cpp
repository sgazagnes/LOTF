
#include <iostream>
#include <set>
#include <sstream>
#include <algorithm>
#include  <cmath>
#include <fstream>
#include <vector>
#include <gsl/gsl_poly.h>
#include <cstdlib>
#include <stdlib.h> 
// Root headers
#include "TFile.h"
#include "TNtuple.h"
#include "TStopwatch.h"
#include "TH1.h"

#include "simon_functions.h"
#include "logc.h"
#include "path_queue.h"
#include "phfitting.h"



double nodeDistanceToLinearFit(double xdet, double ydet, double *x_coef, double *y_coef){
  
  //double xdet = (double) node->m_r/  sqrt( 2*pow(40,2));//node->m_xDet;
  //double ydet = (double)  (node->m_thetaDeg+180.) /360.;

  double newx_coef[2] = {x_coef[0] - xdet, x_coef[1]};
  double newy_coef[2] = {y_coef[0] - ydet, y_coef[1]};

  double vectortanx[2] = {x_coef[1], 0.};
  double vectortany[2] = {y_coef[1], 0.};
	      
  double d[3] = {newx_coef[0]*vectortanx[0] + newy_coef[0]*vectortany[0],
		 newx_coef[0]*vectortanx[1] + newx_coef[1]*vectortanx[0] +
		 newy_coef[0]*vectortany[1] + newy_coef[1]*vectortany[0],
		 newx_coef[1]*vectortanx[1] + newy_coef[1]*vectortany[1]  };

  double x0[2];
  int nroot = gsl_poly_solve_quadratic(d[2], d[1], d[0], x0, x0+1);

  double xIntersect, yIntersect, currDist;
      
  if(nroot ==1){
    xIntersect = gsl_poly_eval(x_coef, 2, x0[0]);
    yIntersect = gsl_poly_eval(y_coef, 2, x0[0]);
    currDist = sqrt(pow(xIntersect - xdet,2) + pow(yIntersect - ydet,2));
  }

		
  for (int j = 0; j < nroot; j++){
    double newx = gsl_poly_eval(x_coef, 2, x0[j]);
    double newy = gsl_poly_eval(y_coef, 2, x0[j]);
    double newdist = sqrt(pow(newx - xdet,2) + pow(newy - ydet,2));
    if(j == 0){
      xIntersect = newx;
      yIntersect = newy;
      currDist = newdist;
    } else if( currDist > newdist) {
      xIntersect = newx;
      yIntersect = newy;
      currDist = newdist;
    }
  }

  return currDist;
}


double nodeDistanceToQuadFit(double xdet, double ydet, double *x_coef, double *y_coef){

  double newx_coef[3] = {x_coef[0] - xdet, x_coef[1], x_coef[2]};
  double newy_coef[3] = {y_coef[0] - ydet, y_coef[1], y_coef[2]};

  double vectortanx[3] = {x_coef[1], 2*x_coef[2], 0.};
  double vectortany[3] = {y_coef[1], 2*y_coef[2], 0.};
	      
  double d[4] = {newx_coef[0]*vectortanx[0] + newy_coef[0]*vectortany[0],
		 newx_coef[0]*vectortanx[1] + newx_coef[1]*vectortanx[0] +
		 newy_coef[0]*vectortany[1] + newy_coef[1]*vectortany[0],
		 newx_coef[0]*vectortanx[2] + newx_coef[1]*vectortanx[1] + newx_coef[2]*vectortanx[0] +
		 newy_coef[0]*vectortany[2] + newy_coef[1]*vectortany[1] + newy_coef[2]*vectortany[0],
		 newx_coef[2]*vectortanx[1] + newy_coef[2]*vectortany[1]};

  double x0[3];
  int nroot = gsl_poly_solve_cubic(d[2]/d[3], d[1]/d[3], d[0]/d[3], x0, x0+1, x0+2);

  //debug("Polynomial coeff : %lf + %lf x + %lf x^2 + %lf x^3", d[0], d[1], d[2], d[3]);
  //	debug("Real roots %d : x1 = %lf x2 =  %lf x3 =  %lf \n\n", nroot, x0[0], x0[1], x0[2]);

  double xIntersect, yIntersect, currDist = -1;
  if(nroot ==1){
    xIntersect = gsl_poly_eval(x_coef, 3, x0[0]);
    yIntersect = gsl_poly_eval(y_coef, 3, x0[0]);
    currDist = sqrt(pow(xIntersect - xdet,2) + pow(yIntersect - ydet,2));
  }

		
  for (int j = 0; j < nroot; j++){
    double newx = gsl_poly_eval(x_coef, 3, x0[j]);
    double newy = gsl_poly_eval(y_coef, 3, x0[j]);
    double newdist = sqrt(pow(newx - xdet,2) + pow(newy - ydet,2));
    if(j == 0){
      xIntersect = newx;
      yIntersect = newy;
      currDist = newdist;
    } else if( currDist > newdist) {
      xIntersect = newx;
      yIntersect = newy;
      currDist = newdist;
    }
  }
  return currDist;
}

   
/* fitNextId */

int fitNextId(CoordGrid &gr, std::vector< GridNode > &Ingrid, PathCandidate &cand, std::vector<int> &next, int k){
  
  // std::vector< GridNode > &Ingrid = gr.m_grid;
  int goodId     = -1;
  int method, degree, nElts;
  std::vector<double> x     =  std::vector<double>( cand.m_x );
  std::vector<double> y     =  std::vector<double>( cand.m_y ); 
  std::vector<GridNode> anchors = std::vector<GridNode>( cand.m_anchors );
  std::vector<int> layers = std::vector<int>( cand.m_layers );
  
  if(k == 0){
    std::reverse(x.begin(),x.end());
    std::reverse(y.begin(),y.end());
    std::reverse(anchors.begin(),anchors.end());
    std::reverse(layers.begin(),layers.end());
  }

   std::vector<int> plausible;
   std::vector<int> uncertain;
   std::vector<int> unlikely;
   std::vector<int> *tocheck;


   //To faster the search, if it's clear that we are going straight
   bool testLayer = false;
   int cum  = layers[layers.size()-1] - layers[layers.size()-3];
   int dirLayer = 0;
   if(cum == 2){
     testLayer = true;
     dirLayer =1;
   } else if(cum == -2){
     testLayer = true;
     dirLayer = -1;
   }
  
  for (int i = 0; i <next.size(); i++){
    int curId = next[i];
    GridNode *node = &Ingrid[gr.Find(curId)];
       // Replacing the node by the correct Intersection point based on its virtual node if it exists
    if(node->m_type == GridNode::VIRTUAL_NODE){
      GridNode *neigh = &Ingrid[gr.Find(node->m_neighbors[0])];
      if(neigh->m_cm.size()>0 && std::find(neigh->m_cm.begin(), neigh->m_cm.end(), cand.m_id) != neigh->m_cm.end()){
	neigh = &Ingrid[gr.Find(node->m_neighbors[1])];
      }
      if(neigh->m_cm.size() > 0 && neigh->m_cm[0] == cand.m_id){
	continue;	
      }

      if(anchors[anchors.size()-1].m_type != GridNode::STT_TYPE_PARA){
	dbgfit("Replacing node %d with %d, and recorrecting the impact coordinates", node->m_detID, neigh->m_detID);
	PointsLineIntersectLive( *neigh, anchors[anchors.size()-1].m_xDet, node->m_xDet,
				 anchors[anchors.size()-1].m_yDet, node->m_yDet); //Output
      }
      node = neigh;
      if(cand.isInCandidate(node->m_detID))
	continue;
    }


    if(testLayer){
      if(node->m_Layer - layers[layers.size()-1] == dirLayer){
	plausible.push_back(node->m_detID);
      } else if(node->m_Layer - layers[layers.size()-1] == 0) {
	uncertain.push_back(node->m_detID);
      } else
	unlikely.push_back(node->m_detID);
    } else
      plausible.push_back(node->m_detID);


    sort( plausible.begin(), plausible.end() );
    plausible.erase( unique( plausible.begin(), plausible.end() ), plausible.end() );
    sort( uncertain.begin(), uncertain.end() );
    uncertain.erase( unique( uncertain.begin(), uncertain.end() ), uncertain.end() );
  }
  if(plausible.size() >= 1){
    tocheck = &plausible;
  } else {
    tocheck = &uncertain;
  }


  // Checking Track angle in polar coord 

  // debug("Points %lf, %lf \t %lf, %lf \t %lf, %lf",r[0], theta[0], r[r.size()/2], theta[theta.size()/2],r[r.size()-1], theta[theta.size()-1]);
  
  // double curv = returnCurvature(x[0], x[x.size()/2], x[x.size()-1], y[0], y[y.size()/2], y[y.size()-1]);

  // float angle = returnAngle(r[0], r[r.size()/2], r[r.size()-1],  theta[0], theta[theta.size()/2], theta[theta.size()-1]);
  //debug("Angle of track is %f", angle);
  
  //
  //if(fabs(angle) < 170){
  //   dbgfit("Quadratic Fit");
  //  method = 1;
    // degree = 2;
    // nElts  = x.size() -1;
    // } else {
  // dbgfit("Linear Fit");

  double minDist = std::numeric_limits<double>::max();

  if(tocheck->size() > 0){
    
    method = 0;
    degree = 1;
    nElts = 3;
    int firstElt = MAX(0, (int) anchors.size()-nElts);
    //}
 
    std::vector<double> p;
    std::vector<double> xfit;
    std::vector<double> yfit;

    
    p.push_back(0.);
    xfit.push_back(anchors[ firstElt].m_xDet);
    yfit.push_back(anchors[ firstElt].m_yDet);

    for (int i = firstElt, inc=0; i <  anchors.size()-1; i++, inc++){ 
      double newval = p[inc] + sqrt(pow(anchors[i+1].m_xDet-anchors[i].m_xDet,2.)
				    + pow(anchors[i+1].m_yDet-anchors[i].m_yDet,2.));
      p.push_back(newval);
      xfit.push_back(anchors[i+1].m_xDet);
      yfit.push_back(anchors[i+1].m_yDet);
    }

    double *x_coef = polyFit(p, xfit, degree);
    double *y_coef = polyFit(p, yfit, degree);   
  
    for (size_t i = 0; i < tocheck->size(); i++){
      int curId = tocheck->at(i);
      int curIdx = gr.Find(curId);
      GridNode *node = &Ingrid[curIdx];
      double xdet = (double) node->m_xDet;
      double ydet = (double) node->m_yDet;  
   
      double newp = p[p.size()-1] + sqrt(pow(xdet-anchors[anchors.size()-1].m_xDet,2.)
					 + pow(ydet-anchors[anchors.size()-1].m_yDet,2.));
      double xest =  method == 0? x_coef[0]+x_coef[1]*newp: x_coef[0]+x_coef[1]*newp+x_coef[2]*newp*newp;
      double yest = method == 0? y_coef[0]+y_coef[1]*newp: (y_coef[0]+y_coef[1]*newp+y_coef[2]*newp*newp);
 
      //dbgfit("Estimated new coord %lf, %lf", xest,yest);
      // dbgfit("Node coord %lf, %lf",xdet, ydet);  

      double currDist =sqrt(pow(xest -xdet, 2) + pow(yest -ydet,2));
      dbgfit("Node %d is at %lf", curId, currDist);
		
      if(minDist > currDist && currDist < 6){ // Empircal threshold...
	minDist = currDist;		  
	goodId = curId;	
      }        
    }
  }
  if(goodId != -1) dbgfit("The good id is %d, and is at distance %lf", goodId,minDist);
  else dbgfit("No good ID found");
  
  return goodId;
}


   
void fittingPhase(CoordGrid &gr, std::vector< GridNode > &Ingrid, std::vector < PathCandidate* > &tracklets, std::vector<pair<int, unsigned short>> idToProcess, char *visited, int **toMergeWith){


  for(unsigned int l = 0; l < tracklets.size(); l++){ // Go for each tracklet
	
    PathCandidate &curCand = *(tracklets[l]);
    dbgfit("Track %d, status %d, length %d", curCand.m_id, curCand.m_finished,curCand.m_length);
	
    if (curCand.m_finished == 3  ) continue;
    //|| curCand.m_length < 5
    // HEAD and TAIL nodes
    GridNode &firstNode = Ingrid[gr.Find(curCand.m_tailNode)];
    GridNode &lastNode  = Ingrid[gr.Find(curCand.m_headNode)];

    dbgfit("Tail node %d (num neigh %d),  head node %d (num neigh %d)",
	   firstNode.m_detID, curCand.m_tailNeigh.size(), lastNode.m_detID, curCand.m_headNeigh.size());

    //We might have to fit in the tail or the head direction, check both
     for(int k = 0; k < 2; k++){
       // k == 0 Tail (first node added), k == 1 Head (last node added, most often)
       //	int prevId = k == 1? curCand.m_headNode: curCand.m_tailNode;
	GridNode *prevNode = k == 1? &lastNode: &firstNode;

	//curent list of neighbors from previous phase
	std::vector<int> *curNeigh = k == 1? &(curCand.m_headNeigh): &(curCand.m_tailNeigh);
	// Whether the tracklets needs to be merged already
	std::vector<unsigned int> *curMerge = k == 1? &(curCand.m_toMergeHead):&(curCand.m_toMergeTail);
	std::vector<int> next;//A vector for the next nodes in the loop

	// If there is no neighbors, or if it already needs to be merged, we pass (for now)
	if(curNeigh->size() == 0){
	  //  dbgfit("k = %d, there is no neighbors in the list, is it consistent with the layer of the node ? %d",
	  //	 k, prevNode->m_Layer);
	  continue;
	}

	if(curMerge->size() > 0){
	  dbgfit("This tracklets has already %d merging partner(s) in that direction ",curMerge->size() );
	  continue;
	}

	
	if (curCand.m_finished == 1 && curNeigh->size() == 0){
	  continue;
	  dbgfit("We previously assumed that this track could be continued, we should check for second order neighbors in the list of %lu remaining ones", idToProcess.size());

	  for(unsigned int n = 0; n < idToProcess.size(); ++n) {
	    int testID 		= idToProcess[n].first;
	    GridNode &testNode  = Ingrid[gr.Find(testID)];
	    float currDist      = sqrt((prevNode->m_x - testNode.m_x) * (prevNode->m_x - testNode.m_x) +
				 (prevNode->m_y - testNode.m_y) * (prevNode->m_y - testNode.m_y));

	    if(currDist<5.){
	      curNeigh->push_back(testID);
	      dbgfit("Adding node %d to potential neighbors, into k %d", testID, k);	     
	    }
	  }

	  dbgfit("Let's now look into the %lu tracklets we found previously", tracklets.size());

	  for(unsigned int n = 0; n < tracklets.size(); ++n) {
	    PathCandidate &testCand = *(tracklets[n]);
	
	    if (testCand.m_finished == 3 || n == l) continue;

	    GridNode &tailNode = Ingrid[gr.Find(testCand.m_tailNode)];
	    GridNode &headNode = Ingrid[gr.Find(testCand.m_headNode)];
	    GridNode Dummy;
	    
	    double currDistTail = IntersectionPointSkeSke(gr, *prevNode, tailNode, Dummy);
	    double currDistHead = IntersectionPointSkeSke(gr, *prevNode, headNode, Dummy);
	    //dbgfit("%d %f, %d %f", tailId, currDistTail, headId, currDistHead);
	    if(currDistTail < 5. || currDistHead < 5.){
	      if(currDistTail <= currDistHead){
		curNeigh->push_back(tailNode.m_detID);
		dbgfit("Adding tail node %d to potential neighbors", tailNode.m_detID);
	      } else {
		curNeigh->push_back(headNode.m_detID);
		dbgfit("Adding head node %d to potential neighbors of head", headNode.m_detID);
	      }
	    }
	  }
	  
	} // End test if cand is status is unclear


	
	if(curNeigh->size() == 0){
	  dbgfit("Still no good candidate has been found");
	  continue;
	}
	
	k == 1? dbgfit("HEAD : Starting fitting next neighbors"):
	  dbgfit("TAIL : Starting fitting next neighbors");

	
	next.insert(next.end(),  curNeigh->begin(),  curNeigh->end());

	bool cond   = next.size() > 0? true: false;

	std::vector<int> virt;
	std::vector<int> *trk = curCand.m_memberList;

	//	int potCm = -1;
	int id = k == 1? trk->at(trk->size() - 2) : 1;
	GridNode &node = Ingrid[gr.Find(id)];		
	   
	// Starting the big loop

	while (cond){

	  std::vector<int> idToRemove;	
	  for(size_t i = 0; i  < next.size(); i++){ 
	    GridNode &node = Ingrid[gr.Find(next[i])];
	    if(node.m_type == GridNode::VIRTUAL_NODE){
	      // Remove second order neighbors from virtual
	      idToRemove.push_back(node.m_neighbors[0]);
	      idToRemove.push_back(node.m_neighbors[1]);
	    }
	  }
	  for(size_t i =0; i < idToRemove.size(); i++)
	    next.erase(std::remove(next.begin(), next.end(), idToRemove[i]), next.end());
	  

	  GridNode *goodNode;	    	    
	  int goodId = fitNextId(gr, Ingrid, curCand, next, k);
	    
	  if (goodId == -1) {
	    dbgfit("No good candidates have been found, stop");
	    dbgfit("Current cm %d: length is %d,  tail node %d  head node %d  Min layer %d, Max layer %d. ", curCand.m_id, curCand.m_length, curCand.m_tailNode, curCand.m_headNode, curCand.m_minLayer, curCand.m_maxLayer);
	    if((curCand.m_minLayer == 0 && curCand.m_maxLayer > 21) || (firstNode.m_LayerLimit == 1 && lastNode.m_LayerLimit == 1)){		 
	      dbgfit("track goes through all layers or makes a loop, likily finished");		 
	      curCand.m_finished = 3;		 
	    } else {
	      dbgfit("Are we missing something?");		 
	      curCand.m_finished = 1;
	    }
	    cond = false;
	    break;
	  }

	  goodNode = &Ingrid[gr.Find(goodId)];

	  //Check that we did not forget a virtual node before
	    
	  if(goodNode->m_type != GridNode::VIRTUAL_NODE){
	    for(size_t i = 0; i < goodNode->m_neighbors.size(); i++){
	      int neighId = goodNode->m_neighbors[i];
	      GridNode *comNode = &Ingrid[gr.Find(neighId)];
	      if(prevNode->IsNeighboring(neighId) && comNode->m_type == GridNode::VIRTUAL_NODE){
		curCand.insertNewNode(gr, Ingrid, &Ingrid[gr.Find(neighId)], k == 0? curCand.m_memberList->begin(): curCand.m_memberList->end());
		visited[neighId] = 1;
	      }
	    }
	  }



	  // Check if the node found belongs to an other track	     
	  if(visited[goodId] == 1){ // To remove Layler limit cond?
	    // if(goodNode->m_cm.size() > 1)
	    dbgfit("This node %d belongs to %d tracks, testing all", goodId, goodNode->m_cm.size());

	    for (size_t i = 0; i< goodNode->m_cm.size(); i++){	      
	      int potCCtoMerge = goodNode->m_cm[i];
	      if(potCCtoMerge == curCand.m_id)
		continue;
	      
	      const auto p = std::find_if(tracklets.begin(), tracklets.end(),
					  [potCCtoMerge](const PathCandidate *obj){ return obj->m_id == potCCtoMerge; } );

	      PathCandidate &neighCand = *(*p); // The CC that the node belongs to	   	   

	      //Find where the node is in the list of the other CC
	      std::vector<int>::iterator it = std::find((neighCand.m_memberList)->begin(),
							(neighCand.m_memberList)->end(), goodId);
	    
	      int index              = std::distance((neighCand.m_memberList)->begin(), it);
	      int nodeIdNeighCand    = neighCand.m_memberList->at(index);


	      dbgfit("Testing the possibility to merge with the node %d of CM %d, tail Node %d, and head Node %d",
		     nodeIdNeighCand, potCCtoMerge, neighCand.m_tailNode,neighCand.m_headNode);

	      // Check if we were already planning a merging
	      if (std::find(curCand.m_toMergeHead.begin(), curCand.m_toMergeHead.end(), potCCtoMerge) != curCand.m_toMergeHead.end() || std::find(curCand.m_toMergeTail.begin(), curCand.m_toMergeTail.end(), potCCtoMerge) != curCand.m_toMergeTail.end() ){
		dbgfit("We were already planning to merge with this track, then let's stop here");
		cond = false;
		break;
	      }  

	      //Check if the node is somewhere in the middle, so it is unlikely that both tracks go together
	      if(nodeIdNeighCand != neighCand.m_headNode && nodeIdNeighCand != neighCand.m_tailNode){
		dbgfit("The index is neither the tail or the head, we should continue");
		// Maybe we should do some tests here as well
		// Can we add this id to say that we should never merge with it?
	      }

	      else{

		int dirNeigh       = 0;
		int caseMerge      = 0;
		int offAnc         = 0;
		int nextNeighId    = -1;
		
		if(nodeIdNeighCand == neighCand.m_headNode){
		  dbgfit("We found a match with the head");
		  //dirNeigh    = -1;           // The position of the next node is at -1
		  offAnc      = neighCand.m_anchors.size()-2;
		  caseMerge   = k == 1? 3: 1; // head to head or tail to head
		}

		else if (nodeIdNeighCand == neighCand.m_tailNode) {
		  dbgfit("We found a match with the tail");
		  //  dirNeigh    = +1;            // The position of the next node is at +1
		  offAnc      = 1;
		  caseMerge   = k == 1? 2: 0;  // head to tail or tail to tail
		}

		dbgfit("We shall test if both tracklets direction are consistent");
		GridNode &prevAnc = k == 1? curCand.m_anchors[curCand.m_anchors.size()-2]: curCand.m_anchors[1];
		//
		//float angle_r = returnAngle(prevNode->m_r, goodNode->m_r, neighCand.m_r[nextNeigh], (prevNode->m_thetaDeg+180.)/360., (goodNode->m_thetaDeg+180.)/360., (neighCand.m_theta[nextNeigh]+180.)/360.);
		dbgfit("PrevAncho %d, %f, %f", prevAnc.m_detID,prevAnc.m_xDet, prevAnc.m_yDet);
		dbgfit("Next neigh Anc %d, %f, %f", neighCand.m_anchors[offAnc].m_detID,neighCand.m_anchors[offAnc].m_xDet, neighCand.m_anchors[offAnc].m_yDet);

		float angle_xy = returnAngle(prevAnc.m_xDet, goodNode->m_xDet, neighCand.m_anchors[offAnc].m_xDet,
					     prevAnc.m_yDet, goodNode->m_yDet, neighCand.m_anchors[offAnc].m_yDet);
		  
		// dbgfit("Angle r with track %f", angle_r);
		dbgfit("Angle xy with track %f", angle_xy);
		  
		if(fabs(angle_xy) > 105){
		  dbgfit("The direction looks consistent between two tracklets, let's merge them");
		  k == 1? curCand.m_toMergeHead.push_back(potCCtoMerge):
		    curCand.m_toMergeTail.push_back(potCCtoMerge);
		  caseMerge % 2 ? neighCand.m_toMergeHead.push_back(curCand.m_id):
		    neighCand.m_toMergeTail.push_back(curCand.m_id);		  		    
		  toMergeWith[curCand.m_id][potCCtoMerge] = caseMerge;
		  curCand.m_finished = 2;
		  //Something to implement for later, we should check that the tracks are consistent
		  cond = false;
		  break;
		} else {
		  dbgfit("The direction is not consistent. We should decide what to do");
		}
	      } // END OF ELSE WE FOUND A MATCH WITH THE HEAD OR TAIL	      
	    } // END OF FOR BOUCLE IN THE CM OF THE NODE
	    if(!cond)
	      break;
	  } // End of check for node belonging to other track

	  //If we did not break before, we should then add this node to the track (head or tail depends on k)
	  curCand.insertNewNode(gr, Ingrid, goodNode, k == 0?
				curCand.m_memberList->begin(): curCand.m_memberList->end());
	  visited[goodId] = 1;
	  next.clear();
     	  
	  //Finding next neighbors (We shoudl add something to stop if all neighbors are assigned?)
	  for(size_t i = 0; i < goodNode->m_neighbors.size(); i++){	      
	    int neighId = goodNode->m_neighbors[i];
	    if(curCand.isInCandidate(neighId)) continue;	      
	    int neighIdx = gr.Find(neighId);
	    GridNode *neighNode = &Ingrid[neighIdx];
	    dbgfit("Pushing this node %d to the list", neighId );
	    next.push_back(neighId);		
	  }

	  // If the node had a parent from an other track but we decided not to match
	  // with the track, then we should add this parent	  
	  if(goodNode->parent != -1){
	    next.push_back(goodNode->parent);
	    dbgfit("Adding the parent node %d", goodNode->parent );
	  }

	  // If we found some neighbors, then we can continue;
	  if(next.size() > 0){	      
	    //  prevId   = goodId;
	    prevNode = goodNode;
	    virt.clear(); // Clearing the list of virtuals that we did not use?
	    if(goodNode->m_Layer == 0){ // This is to be removed
	      dbgfit("We found neighbors but we reached the end ?, should check this");
	      curCand.m_finished = 1;
	      cond = false;
	    }			 
	  }
	  
	  else {	      
	    dbgfit("we have no more neighbors");
	    cond = false;	     
	  }
	} // END OF BIG LOOP
     } // END OF FOR k = 0 or 1

     dbgfit("Finished with current tracklet");
     dbgfit("Current cm %d: length is %d,  tail node %d  head node %d  Min layer %d, Max layer %d. ", curCand.m_id, curCand.m_length, curCand.m_tailNode, curCand.m_headNode, curCand.m_minLayer, curCand.m_maxLayer);
     if((firstNode.m_LayerLimit == 1 && lastNode.m_LayerLimit == 1)){		 
       dbgfit("track goes through all layers or makes a loop, likily finished");		 
       curCand.m_finished = 3;		 
     } else {
       dbgfit("Are we missing something?");		 
       curCand.m_finished = 1;
     }

     dbgfit("Moving to next track\n");
  } // END OF ALL TRACKLETS
}
   
