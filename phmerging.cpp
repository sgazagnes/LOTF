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

#include "gridNode.h"
#include "logc.h"
#include "simon_functions.h"
#include "path_queue.h"
#include "phmerging.h"


void addTracklets (CoordGrid &gr, std::vector< GridNode > &Ingrid, PathCandidate *newCand, PathCandidate &mergeCand,  int curdir, int mergedir){


  if(mergedir == 1){
    for(int i =  (mergeCand.m_memberList)->size()-1; i >=0; i--){
      int curid = (mergeCand.m_memberList)->at(i);
      //   debug("adding id %d", curid);
      if(newCand->isInCandidate(curid)) continue;
      int curidx = gr.Find(curid);
      GridNode* node = &Ingrid[curidx];
      newCand->insertNewNode(gr, Ingrid, node, curdir == 1? newCand->m_memberList->end(): newCand->m_memberList->begin());
    }
    
  } else {

    for(size_t i =  0; i < mergeCand.m_memberList->size(); i++){
      int curid = mergeCand.m_memberList->at(i);
      // debug("adding id %d", curid);

      if(newCand->isInCandidate(curid)) continue;
      int curidx = gr.Find(curid);
      GridNode* node = &Ingrid[curidx];
      newCand->insertNewNode(gr, Ingrid, node, curdir == 1? newCand->m_memberList->end(): newCand->m_memberList->begin());
    }
  }
  
  mergeCand.m_isMerged = 1;
  mergeCand.m_isValid = 0;

}


void mergeTracks (CoordGrid &gr, std::vector< GridNode > &Ingrid, std::vector < PathCandidate* > &tracklets, int *candidateId ){

  int curCandId = *candidateId;
  for(unsigned int l = 0; l < tracklets.size(); l++){
     
    PathCandidate &curCand = *(tracklets[l]);
    dbgmerge("Current tracklet %d, isValid %d", curCand.m_id, curCand.m_isValid );
      
    if (curCand.m_isValid == 0 && (curCand.m_finished ==  3 || (curCand.m_toMergeHead.size() == 0 && curCand.m_toMergeTail.size() ==0     && curCand.m_isMerged == 0))){
      dbgmerge("This CM is either finished or has no further close neighbors");
      //check that it has not merging candidates
      if (tracklets[l]->m_toMergeHead.size() != 0  || tracklets[l]->m_toMergeTail.size() != 0)
	error("This looks finished but has to be merged ???? ");

      PathCandidate *newCand 	= new PathCandidate();// Create a new candidate
	       
      newCand->m_id 		= curCandId++;// Set id
      newCand->m_tailNode       = curCand.m_tailNode;
	  
      for(size_t i = 0; i < (curCand.m_memberList)->size(); i++){
	int curid = (curCand.m_memberList)->at(i);
	int curidx = gr.Find(curid);
	GridNode* node = &Ingrid[curidx];
	newCand->insertNewNodeFinal(gr, Ingrid, node,newCand->m_memberList->end());
      }

      dbgmerge("Pushing new merged cm %d:  length is %d, tail node %d, head node %d, min layer %d, max layer %d,    IsOnSectorLimit %d, finished ? %d. ", newCand->m_id, newCand->m_length, newCand->m_tailNode, newCand->m_headNode,newCand->m_minLayer, newCand->m_maxLayer, newCand->m_isOnSectorLimit, newCand->m_finished);
      
      curCand.m_isValid = 0;
      curCand.m_isMerged = 1;
      newCand->m_isValid = 1;
      newCand->m_finished = 3;
	


      tracklets.push_back(newCand);
      
      //   curCand.m_isValid = 1;	
    } else if (curCand.m_isMerged == 1){
      dbgmerge("This CM has already been merged");
      curCand.m_isValid = 0;
    }

    else {
 
      int sizeMergeHead =  curCand.m_toMergeHead.size();
      int sizeMergeTail =  curCand.m_toMergeTail.size();

      if(sizeMergeHead == 0 && sizeMergeTail == 0){
	dbgmerge("Cand does not look like it is done, but we don't have merging to do");
	continue;
      }
      
      for(size_t i = 0; i < tracklets[l]->m_toMergeHead.size(); i++){

	dbgmerge("Merge in head with %d", tracklets[l]->m_toMergeHead[i]);
      }

      for(size_t i = 0; i < tracklets[l]->m_toMergeTail.size(); i++){

	dbgmerge("Merge in tail with %d", tracklets[l]->m_toMergeTail[i]);
      }

      // continue;
      PathCandidate *newCand 	= new PathCandidate();// Create a new candidate
	       
      newCand->m_id 		= curCandId++;// Set id
      newCand->m_tailNode       = curCand.m_tailNode;
	  
      for(size_t i = 0; i < (curCand.m_memberList)->size(); i++){
	int curid = (curCand.m_memberList)->at(i);
	int curidx = gr.Find(curid);
	GridNode* node = &Ingrid[curidx];
	newCand->insertNewNodeFinal(gr, Ingrid, node,newCand->m_memberList->end());
      }

      int curCandId =  curCand.m_id;

      if(sizeMergeHead != 0){

	int idToMerge = curCand.m_toMergeHead[0];
	bool cond = true;
	  
	while(cond){
	    
	  dbgmerge("Tracklets %d needs to be merged in head with %d", curCand.m_id, idToMerge);
	  const auto p = std::find_if(tracklets.begin(), tracklets.end(),
				      [idToMerge](const PathCandidate *obj){ return obj->m_id == idToMerge;});
	    
	  PathCandidate &mergeCand = *(*p);
	  
	  if(std::find(mergeCand.m_toMergeHead.begin(), mergeCand.m_toMergeHead.end(), curCandId) != mergeCand.m_toMergeHead.end()) {
	      
	    dbgmerge("Merging head with head");
	    addTracklets(gr, Ingrid, newCand, mergeCand, 1, 1);	    
	    if(mergeCand.m_toMergeTail.size() > 0){
	      curCandId = mergeCand.m_id;
	      idToMerge = mergeCand.m_toMergeTail[0];
	      dbgmerge("ONE MORE TO PUSH %d", idToMerge);
	      if(idToMerge == curCand.m_id)
		cond = false;
	    } else
	      cond = false;
	      
	  } else {
	      
	    dbgmerge("Merging head with tail");	    
	    addTracklets(gr, Ingrid, newCand, mergeCand, 1, 0);
	    
	    if(mergeCand.m_toMergeHead.size()){
	      curCandId = mergeCand.m_id;
	      idToMerge = mergeCand.m_toMergeHead[0];
	      if(idToMerge == curCand.m_id)
		cond = false;
	      dbgmerge("ONE MORE TO PUSH %d", idToMerge);
	    } else
	      cond = false;

	  }
	  curCand.m_isValid = 0;
	  curCand.m_isMerged = 1;
	} // while continu
      } // HEAD MERGING

      if(sizeMergeTail != 0){

	int idToMerge = curCand.m_toMergeTail[0];
	bool cond = true;
	  
	while(cond){

	  dbgmerge("Tracklets %d needs to be merged in tail with %d", curCand.m_id, idToMerge);

	  const auto p = std::find_if(tracklets.begin(), tracklets.end(),
				      [idToMerge](const PathCandidate *obj){ return obj->m_id == idToMerge; } );
	    
	  PathCandidate &mergeCand = *(*p);
	  
	  if(std::find(mergeCand.m_toMergeHead.begin(), mergeCand.m_toMergeHead.end(), curCandId) != mergeCand.m_toMergeHead.end()) {
	      
	    dbgmerge("Merging tail with head");	    
	    addTracklets(gr, Ingrid, newCand, mergeCand, 0, 1);
	    
	    if(mergeCand.m_toMergeTail.size() > 0 ){
	      curCandId = mergeCand.m_id;
	      idToMerge = mergeCand.m_toMergeTail[0];
	      if(idToMerge == curCand.m_id)
		cond = false;
	      dbgmerge("ONE MORE TO PUSH %d (not implemented yet)", idToMerge);
	    } else
	      cond = false;
	      
	  } else {
	      
	    dbgmerge("Merging head with tail");	    
	    addTracklets(gr, Ingrid, newCand, mergeCand, 0, 0);
	    
	    if(mergeCand.m_toMergeHead.size() > 0){
	      curCandId = mergeCand.m_id;
	      idToMerge = mergeCand.m_toMergeHead[0];
	      if(idToMerge == curCand.m_id)
		cond = false;
	      dbgmerge("ONE MORE TO PUSH %d (not implemented yet)", idToMerge);
	    } else
	      cond = false;

	  }
	  curCand.m_isValid = 0;
	  curCand.m_isMerged = 1;
	} // while continu

      } // Tail MERGING
	  
      //	else
      //	  error("ISSUE, no one to merge with ?");
	
      curCand.m_isValid = 0;
      curCand.m_isMerged = 1;
      newCand->m_isValid = 1;
      newCand->m_finished = 3;
	
      dbgmerge("Pushing new merged cm %d:  length is %d, tail node %d, head node %d, min layer %d, max layer %d,    IsOnSectorLimit %d, finished ? %d. ", newCand->m_id, newCand->m_length, newCand->m_tailNode, newCand->m_headNode,newCand->m_minLayer, newCand->m_maxLayer, newCand->m_isOnSectorLimit, newCand->m_finished);

      tracklets.push_back(newCand);
    } // ELSE NOT FINISHED

  }
  *candidateId = curCandId;
}// FOR TRACKLETS
// }
