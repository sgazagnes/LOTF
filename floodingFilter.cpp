
#include <iostream>
#include <set>
#include <sstream>
#include <algorithm>
#include  <cmath>
#include <fstream>
#include <vector>

// Root headers
#include "TFile.h"
#include "TNtuple.h"
#include "TStopwatch.h"
#include "TH1.h"

// Local headers
#include "simon_functions.h"
#include "reconstruction.h"
// DEBUG AND STORE definitions

#define EXCLUDE_STTSKEWED_PLOT   0
#define EXCLUDE_VIRTUALS_IN_PLOT 0

#define WRITE_CONNECTED_COMPONENTS 1
#define EVALUATE_ERROR 0
#define INCLUDE_MVD_INOUTPUT_TRACK 0
#define PRINT_TRACKS_DEBUG_INFO 1
#define PRINT_DEBUG_INFO_COMP_MATCH 0
#define PRINT_ERROR_EXTRA_DEBUG_INFO 1
#define WRITE_CONNECTED_COMPONENTS_JSON 0
#define WRITE_CM_ASCII_FILE 1


//______________________ BEGIN MCTrackPoints _____________________________
std::vector< std::vector < MCTrackObject* >* >*
MCTrackPoints( std::vector < std::vector<HitCoordinate*>* > const &evtData)
{
  info("Extracting MC tracks for %d events", evtData.size());
  //	    << " events.\n";
  // Output Parameter
  std::vector< std::vector < MCTrackObject* >* >* outVar =
    new std::vector< std::vector < MCTrackObject* >* >();

  int numTracks = 0;
  std::vector< int > idtracks;

   for(size_t e = 0; e < evtData.size(); ++e) {
    std::vector<HitCoordinate*> const *Current_Event = evtData[e];
    // numTracks = 0;
    // Find out how many MC tracks are available
    for(size_t i = 0; i < Current_Event->size(); ++i) {
      HitCoordinate const *currentHit = Current_Event->at(i);
      if( currentHit->m_trackID != HIT_EXCLUSION && !(std::find(idtracks.begin(), idtracks.end(), currentHit->m_trackID) != idtracks.end())) {
	printf("%d \n", currentHit->m_trackID);
	//	numTracks = currentHit->m_trackID;
	idtracks.push_back(currentHit->m_trackID);
      }
    }// END current event loop
    numTracks = idtracks.size();
    debug("Event %d contains %d tracks", e,(numTracks));
    
   // Now, we know the number of available tracks for the current
    // event; We can allocate memory.
    std::vector < MCTrackObject* >* evtTracks = new std::vector < MCTrackObject* >();
    for(int j = 0; j < numTracks; ++j) {
      MCTrackObject *trk = new MCTrackObject(); 
      evtTracks->push_back(trk);
    }
    // Fill the list

    int curTrack = 0;
    for(size_t k = 0; k < Current_Event->size(); ++k) {
      HitCoordinate const *currentHit = Current_Event->at(k);
      if(currentHit->m_trackID != HIT_EXCLUSION) {
	std::vector<int>::iterator it = std::find(idtracks.begin(), idtracks.end(), currentHit->m_trackID);
	int index = std::distance(idtracks.begin(), it);
	//	printf("%d, %d \n", currentHit->m_trackID,index);
	int trackPos = std::distance(idtracks.begin(), it); //currentHit->m_trackID;
	point3D spacePoint;
	spacePoint.m_x = currentHit->mx;
	spacePoint.m_y = currentHit->my;
	spacePoint.m_z = currentHit->mz;
	if(currentHit->type == HitCoordinate::STT_TYPE) {
	  ((evtTracks->at(trackPos))->m_pointSTTCoordList).push_back(spacePoint);
	  ((evtTracks->at(trackPos))->m_STT_Component).push_back(currentHit->m_detID);
	}
	else if(currentHit->type == HitCoordinate::MVD_TYPE) {
	  ((evtTracks->at(trackPos))->m_pointMVDCoordList).push_back(spacePoint);
	  ((evtTracks->at(trackPos))->m_MVD_Component).push_back(currentHit->m_detID);
	}
      }// END if not HIT_EXCLUSION
    }// END
    outVar->push_back(evtTracks);
  }// END events loop
  return outVar;
  /*  int realnumTracks = 0;

  std::vector<int> IDtracks ;


  for(size_t e = 0; e < evtData.size(); ++e) {
    std::vector<HitCoordinate*> const *Current_Event = evtData[e];
    // numTracks = 0;
    // Find out how many MC tracks are available
    for(size_t i = 0; i < Current_Event->size(); ++i) {
      HitCoordinate const *currentHit = Current_Event->at(i);
      if (std::find(IDtracks.begin(), IDtracks.end(), currentHit->m_trackID) == IDtracks.end())
	IDtracks.push_back(currentHit->m_trackID);

      if( currentHit->m_trackID > numTracks &&
	  currentHit->m_trackID != HIT_EXCLUSION) {
	numTracks = currentHit->m_trackID;
      }

    }// END current event loop

#if (PRINT_TRACKS_DEBUG_INFO > 0)
    std::cout << "\t<-I-> Event " << e
	      << " Contains " << (IDtracks.size())
	      << " Tracks.\n";
#endif
    // Now, we know the number of available tracks for the current
    // event; We can allocate memory.
    std::vector < MCTrackObject* >* evtTracks = new std::vector < MCTrackObject* >();
    for(int j = 0; j <= numTracks; ++j) {
      MCTrackObject *trk = new MCTrackObject(); 
      evtTracks->push_back(trk);
    }
    // Fill the list
    for(size_t k = 0; k < Current_Event->size(); ++k) {
      HitCoordinate const *currentHit = Current_Event->at(k);
      if(currentHit->m_trackID != HIT_EXCLUSION) {
	int trackPos = currentHit->m_trackID;
	point3D spacePoint;
	spacePoint.m_x = currentHit->mx;
	spacePoint.m_y = currentHit->my;
	spacePoint.m_z = currentHit->mz;
	if(currentHit->type == HitCoordinate::STT_TYPE) {
	  ((evtTracks->at(trackPos))->m_pointSTTCoordList).push_back(spacePoint);
	  ((evtTracks->at(trackPos))->m_STT_Component).push_back(currentHit->m_detID);
	}
	else if(currentHit->type == HitCoordinate::MVD_TYPE) {
	  ((evtTracks->at(trackPos))->m_pointMVDCoordList).push_back(spacePoint);
	  ((evtTracks->at(trackPos))->m_MVD_Component).push_back(currentHit->m_detID);
	}
      }// END if not HIT_EXCLUSION
    }// END
    outVar->push_back(evtTracks);
  }// END events loop
  return outVar;*/
}
//______________________ END MCTrackPoints _______________________________







void floodingFilter(std::string const &OutFileName,int firstEvt, int lastEvt)
{
  TStopwatch timer;

  // Structure to hold the detector data (grid)
  std::vector < GridNode > detNodes;
  // File and structure to hold the output coordinates.
  TFile Out_Put_File(OutFileName.c_str(),"RECREATE","Outputfile Created by performFilter", 9);

  // Collected coordinates for tracks.
  TNtuple coord ("CoordCollected" , "Collected Coordinates in x y plane", "x:y:z:x_Det:y_Det:z_Det");

  // Ntuple to hold Error values for all events available in the
  // current events set. The value is evaluated per image.
  std::string errorParameter = "Error_underMerge:Error_overMerge:TotalError";
  errorParameter += ":Error_underMergeNorm:Error_overMergeNorm:TotalErrorNorm";
  // Create Ntuple to hold parameters.
  TNtuple ErrorNtuple("ErrorEstimate","Segmentation error values", errorParameter.c_str());

  // Second error type. Per track error value. Based on curvature
  // data.
  std::string PerTrkErrPars = "misMatched:BestMatchMCLength:CurrentTrackLength";
  PerTrkErrPars += ":MCMinCurrentLength:CurrentMinMCLength";
  PerTrkErrPars += ":UnderMergeError:OverMergeError:MC_a:MC_b:MC_r:MC_E:tr_a:tr_b:tr_r:tr_E";
  TNtuple ErrorNtuplePerTrack("PerTrackError","Per track values of error", PerTrkErrPars.c_str());
  
  // NTuple to hold the coordinates of all connected components.
  std::string ConnCompPar = "EvtNum:CompNum:tubeId:x:y:z:r:thetaDeg:x_Det:y_Det:z_Det";
  TNtuple ConnectedCoord ("ConnectedComponents", "Connected component Coordinates",
			  ConnCompPar.c_str());
  // Hold number of components per event
  TNtuple ComponentPerEvt ("ComponentPerEvt", "Component per event","evtNum:numComponents");


  /* Read all data directly from sim, digi and parameter files
     (OLD) */
  std::vector < std::vector<HitCoordinate*>* >* Hit_coords = 
    CollectSttMvdPoints(detNodes, Out_Put_File, firstEvt, lastEvt);

  std::vector< std::vector < MCTrackObject* >* > *MC_Tracks = MCTrackPoints(*Hit_coords);
  
  // Write event info to output  
  WriteEventPlotsToFile( (*Hit_coords), Out_Put_File);
  
  // Create an empty grid object
  CoordGrid gr;
  
  // Init Grid for STT detector nodes (fill the map).
  gr.Initialize(detNodes);

  TNtuple Layers("LayerLimits","Layer Limits.","x:y:det_z:z");
  TNtuple Sections("SectionsLimits","Section Limits.","x:y:det_z:z");
  // Isolate Sector and Layer Limits
  isolateSectorAndLayerLimits(gr, Sections, Layers);
  Sections.Write();
  Layers.Write();

  TNtuple* OrigGrid = GridToNtuple(detNodes, "OrigGridCoord");
  OrigGrid->SetMarkerStyle(8);
  OrigGrid->SetMarkerSize(0.2);
  OrigGrid->SetMarkerColor(kBlack);
  OrigGrid->Write();

  info("Fix neighbouring before extension.");
  fixNeighboring(gr);
  std::vector < GridNode > VNodes;
  ///Compute_Add_VirtualNodes_Neigbor2(gr, VNodes);// Old method
  Add_VirtualNodes(gr, VNodes);

  // Compute_Virtual_InterSector_Nodes(gr, 6,VNodes);
   TNtuple* virtualTubes = GridToNtuple(VNodes, "VirtualNodes");
  virtualTubes->SetMarkerStyle(8);
  virtualTubes->SetMarkerSize(0.2);
  virtualTubes->SetMarkerColor(kMagenta);
  virtualTubes->Write();
      /*
   * Extend the grid with new virtual points and fix the missing
   * neighboring relations
   */
  // Extend the grid by virtual nodes between the layer.


  info("Extending the grid by %d virtual nodes between the layers.", VNodes.size());
  gr.ExtendedGrid(VNodes);
  // Fix_InterSector_Nodes(gr, 6);

  info("Fixing neighbouring after extension.");
  fixNeighboring(gr);
  TNtuple* extendedGrid = GridToNtuple(gr.m_grid, "ExtendedGrid");
  extendedGrid->SetMarkerStyle(7);//8
  extendedGrid->SetMarkerSize(0.3);
  extendedGrid->SetMarkerColor(17);//41
  extendedGrid->Write();

  info("Total number of tubes after extension = %d", gr.GetNumNodes());
  // Delete allocated memory
  delete (OrigGrid);
  delete virtualTubes;
  delete extendedGrid;
  unsigned int totalnumEvt = Hit_coords->size();
  // Start the timer.

  if(true){

  std::vector< std::set<int>* >* connectedComp = 0;
  connectedComp = new std::vector< std::set<int>* >();
  
  timer.Start();

  /* Fill the grid with the current hit points and process. Handles
     each event separately.*/
  // Event loop

  for(size_t k = 0; k < Hit_coords->size(); ++k) {
    // Data for the current event
    info("Processing event: %d", k);
    std::vector<HitCoordinate*> const *dd = 0;
    dd = Hit_coords->at(k);
    if(dd) {
      gr.FillGrid(*dd);
    }

    
    // std::vector< GridNode > oriIngrid = gr.m_grid;  
    std::vector< GridNode > &Ingrid = gr.m_grid;  
    std::vector< int > activeId;
    std::vector< int > remainingActiveId;


    /* Pushing all active detectors into queue */
    
    int nactiveAll = 0, nactiveReal = 0;
    for(unsigned int n = 0; n < Ingrid.size(); ++n) {
      GridNode &curNode = Ingrid[n];
      if(curNode.m_active){
	int NodeId = Ingrid[n].m_detID;
	if(curNode.m_type != GridNode::VIRTUAL_NODE  ) {
	  activeId.push_back(NodeId);
	  nactiveReal++;
	  int n_neighbors = curNode.m_neighbors.size();    /* Keep only active neighbors  */
	  for  ( int i = 0; i < curNode.m_neighbors.size(); i++){
	    int neigh_ID = curNode.m_neighbors[i];

	    int neigh_index = gr.Find(neigh_ID);
	    if(!Ingrid[neigh_index].m_active){

	      (curNode.m_neighbors).erase((curNode.m_neighbors).begin()+i);
	      i--;
	    }
	  }
	}	
	nactiveAll++;
      }      
    }
    
    info("Found %d active detectors (%d with virtuals)", nactiveReal, nactiveAll);

    bool firsttrack =true;
    bool fitting = true;
    bool merging = true;

    std::vector < PathCandidate* > tracklets;
    int candidateId 	= 0;
    char *visited 	= (char *) calloc(15000, sizeof(char));

    info("First step, let's find the obvious tracks");

    if(firsttrack) findEasyTracks (gr, tracklets, activeId, visited, &candidateId);

    info("After first step, we found %d tracklets", candidateId);

    std::sort(tracklets.begin(), tracklets.end(), compareTwoPathsLength); 
    for(unsigned int l = 0; l < tracklets.size(); l++){
      PathCandidate &curCand = *(tracklets[l]);
      info("Cm %d has length %d", curCand.m_id, curCand.m_length);
    }

    //   char *visitedTracks = (char *) calloc(tracklets.size(), sizeof(char));


    
    /* ++++++++++++++++++++++++++++++++++++++++++++++++++++ */
      
    /*                Let'sdo some fitting !!!!             */
      
    /* ++++++++++++++++++++++++++++++++++++++++++++++++++++ */

    
    if(fitting){

      int **sayYes =  (int **) calloc(tracklets.size(), sizeof(int*));
      for (int i =0; i < tracklets.size(); i++)
	sayYes[i] = (int *) calloc(tracklets.size(), sizeof(int));

      
      for(unsigned int l = 0; l < tracklets.size(); l++){ // Go for each tracklet
	
	PathCandidate &curCand = *(tracklets[l]);
	info("\n\n\n");
	info("NEW Track %d, status %d, length %d", curCand.m_id, curCand.m_finished,curCand.m_length);
	
	if (curCand.m_finished >= 2 || curCand.m_length < 5 ) continue;
	
	int first =  curCand.m_tailNode;
	int firstIdx = gr.Find(first);
	GridNode &firstNode = Ingrid[firstIdx];

	int last =  curCand.m_headNode;
	int lastIdx = gr.Find(last);
	GridNode &lastNode = Ingrid[lastIdx];
	
	info("Tail node %d (num neigh %d),  head node %d (num neigh %d)", first, curCand.m_tailNeigh.size(), last, curCand.m_headNeigh.size());

	for(int k = 0; k < 2; k++){ // First Tail, then Head
	  int prevId = k == 1? curCand.m_headNode: curCand.m_tailNode;
	  //int curLayer = k == 1? lastNode.m_Layer: firstNode.m_Layer;
	  // int layerCurDiff;
	  GridNode *prevNode = k == 1? &lastNode: &firstNode;
	
	  std::vector<int> *curNeigh = k == 1? &(curCand.m_headNeigh): &(curCand.m_tailNeigh);
	  std::vector<unsigned int> *curMerge = k == 1? &(curCand.m_toMergeHead):&(curCand.m_toMergeTail);
	  std::vector<int> next;

	  if(curNeigh->size() == 0 || curMerge->size() > 0)
	    continue;

	  k == 1? info("HEAD : Fitting next neighbors "): info("TAIL : Fitting next neighbors");
	      
	  for(int i = 0; i  < curNeigh->size(); i++){ 
	    int id = curNeigh->at(i);
	    int idx = gr.Find(id);
	    GridNode &node = Ingrid[idx];
	    for (int j = 0; j < node.m_neighbors.size(); j++) {
	      if(node.m_type == GridNode::VIRTUAL_NODE){ // Remove second order neighbors from virtual
		int neigh1 = node.m_neighbors[0];
		int neigh2 = node.m_neighbors[1];
		curNeigh->erase(std::remove(curNeigh->begin(), curNeigh->end(), neigh1), curNeigh->end());
		curNeigh->erase(std::remove(curNeigh->begin(), curNeigh->end(), neigh2), curNeigh->end());
	      }
	    }
	  }
	  
	  next.insert(next.end(),  curNeigh->begin(),  curNeigh->end());

	  bool cond = next.size() > 0? true: false;
	  bool adding = true;

	  std::vector<int> virt;
	  std::vector<int> *trk = curCand.m_memberList;

	  int potCm = -1;
	  int id = k == 1? trk->at(trk->size() - 2) : 1;
	  int idx = gr.Find(id);
	  GridNode &node = Ingrid[idx];		
	   

	  while (cond){
	    
	    GridNode *goodNode;	    
	    
	    int goodId = fitNextId(gr, curCand, next, k);
	    
	    if (goodId == -1) {
	      info("No good candidates have been found, stop");
	      info("Current cm %d: \n \t length is %d, \n \t tail node %d \n \t head node %d \n \t Min layer %d, \n \t Max layer %d. ", curCand.m_id, curCand.m_length, curCand.m_tailNode, curCand.m_headNode, curCand.m_minLayer, curCand.m_maxLayer);
	      if((curCand.m_minLayer == 0 && curCand.m_maxLayer > 21) || (firstNode.m_LayerLimit == 1 && lastNode.m_LayerLimit == 1)){		 
		info("track goes through all layers or makes a loop, likily finished");		 
		curCand.m_finished = 3;		 
	      } else {
		curCand.m_finished = 2;
		cond = false;
		break;
	      }
	    }


	      
	    int goodIdx = gr.Find(goodId);
	    goodNode = &Ingrid[goodIdx];

	      
	    	    //Check that we did not forget a virtual node before
	    
	    if(goodNode->m_type != GridNode::VIRTUAL_NODE){
	      for(int i = 0; i < goodNode->m_neighbors.size(); i++){
		int neighId = goodNode->m_neighbors[i];
		int idx = gr.Find(neighId);
		GridNode *comNode = &Ingrid[idx];
		if(prevNode->IsNeighboring(neighId) && comNode->m_type == GridNode::VIRTUAL_NODE){
		  int idx = gr.Find(neighId);
		  curCand.insertNewNode(gr, &Ingrid[idx], k == 0? curCand.m_memberList->begin(): curCand.m_memberList->end());
		  visited[neighId] = 4;

		}
	      }
	    }
	    
	    // Automatically add next neighbord if we just found a virtual node
	    
	    if(goodNode->m_type == GridNode::VIRTUAL_NODE){
	      curCand.insertNewNode(gr,goodNode, k == 0? curCand.m_memberList->begin(): curCand.m_memberList->end());
	      visited[goodId] = 4;
	      int neighId = curCand.isInCandidate(goodNode->m_neighbors[0]) ? goodNode->m_neighbors[1]:goodNode->m_neighbors[0];
	      int neighIdx = gr.Find(neighId);
	      goodNode = &Ingrid[neighIdx];
	      info("We added a virtual node, also adding the next neighbord %d", neighId);
	      goodId = neighId;
	    }
	    
	    
	    if(visited[goodId] == 4 && !goodNode->m_LayerLimit ){

	      potCm = goodNode->m_cm[0];
	      const auto p = std::find_if(tracklets.begin(), tracklets.end(),
					  [potCm](const PathCandidate *obj){ return obj->m_id ==potCm; } );

	      PathCandidate &neighCand = *(*p);

	      bool willMerge = neighCand.m_toMergeHead.size() == 0 && neighCand.m_toMergeTail.size() == 0?false:true;
	      /*     for(int i = 0; i<tracklets.size(); i++){
		if(sayYes[potCm][i] != 0)
		  willMerge = true;
		  }*/

	      if(!willMerge){

		std::vector<int>::iterator it = std::find((neighCand.m_memberList)->begin(), (neighCand.m_memberList)->end(), goodId);
		int index = std::distance((neighCand.m_memberList)->begin(), it);
		int id = neighCand.m_memberList->at(index);
		int nextNeigh;
		int dir = 0;

		info("This node %d already belongs to a CM (%d), tail Node %d, and head Node %d", id, goodNode->m_cm[0],neighCand.m_tailNode,neighCand.m_headNode);


		if(id != neighCand.m_headNode && id != neighCand.m_tailNode){
		  info("The index is neither the tail or the head, we should continue");
		}else{
		  if(id == neighCand.m_headNode){
		    info("We are merging in the tail direction");
		    nextNeigh = index-1;
		    dir = k == 1? 3: 1; // head to head or tail to head
		    k == 1? curCand.m_toMergeHead.push_back(potCm):  curCand.m_toMergeTail.push_back(potCm);
		    neighCand.m_toMergeHead.push_back(curCand.m_id);
		    sayYes[curCand.m_id][potCm] = dir;
		  } else if ( id == neighCand.m_tailNode ) {
		    info("We are merging in the head direction");
		    nextNeigh= index+1;
		    dir = k == 1? 2: 0; // head to tail or tail to tail
		    k == 1? curCand.m_toMergeHead.push_back(potCm):  curCand.m_toMergeTail.push_back(potCm);
		    neighCand.m_toMergeTail.push_back(curCand.m_id);
		    sayYes[curCand.m_id][potCm] = dir;
		  }
		  
		  float angle_r = returnAngle(prevNode->m_r, goodNode->m_r, neighCand.m_r[nextNeigh], (prevNode->m_thetaDeg+180.)/360., (goodNode->m_thetaDeg+180.)/360., (neighCand.m_theta[nextNeigh]+180.)/360.);
		  
		  float angle_xy = returnAngle(prevNode->m_x, goodNode->m_x, neighCand.m_x[nextNeigh], prevNode->m_y, goodNode->m_y, neighCand.m_y[nextNeigh]);
		  
		  info("Angle r track %f", angle_r);
		  info("Angle xy track %f", angle_xy);
		  
		  break;
		}
	      } else 
		info("we did not merge because we were in the middle of an other track");
	    } // Visited good id 4
		/*else {
		  debug("we must find which is the direction");
		  int currr =  neighCand.m_r[index-1] - goodNode->m_r;
		  int prevr = goodNode->m_r - prevNode->m_r;
		  if((currr >= 0 && prevr >= 0) || (currr < 0 && prevr < 0)){
		    debug("Middle node, We are checking in the tail direction");
		    nextNeigh= index-1;
		    dir = k == 1? 3: 1; // head to head or tail to head
		  } else {
		    debug("Middle node, We are checking in the head direction");
		    nextNeigh= index+1;
		    dir = k == 1? 2: 0; // head to head or tail to head
		  }
		  }*/



	       
		 /*  if(potCm == goodNode->m_cm[0]){
		  
		debug("Two nodes in a row belong to the same cm %d, let's stop and merge them later", potCm);
		cond = false;
		const auto p = std::find_if(tracklets.begin(), tracklets.end(),
					    [potCm](const PathCandidate *obj){ return obj->m_id == potCm; } );
		

		curCand.m_toMergeHead.push_back(potCm);

		if(prevId == (*p)->m_headNode){
		  (*p)->m_toMergeHead.push_back(curCand.m_id);
		  debug("Merging head to head");
		}
		  
		else if(prevId ==  (*p)->m_tailNode){
		  (*p)->m_toMergeTail.push_back(curCand.m_id);
		  debug("Merging head to tail");
		}

		else {
		  std::vector<int>::iterator it = std::find(((*p)->m_memberList)->begin(), ((*p)->m_memberList)->end(), prevId);
		  int index = std::distance(((*p)->m_memberList)->begin(), it);
		  if(index >= (tracklets[potCm]->m_memberList)->size()/2){
		    debug("Merging head to head (to check)");
		    ( (*p)->m_toMergeHead).push_back(curCand.m_id) ;
		  } else{
		    debug("Merging head to tail (to check)");
		    ( (*p)->m_toMergeTail).push_back(curCand.m_id);
		  }
		}
		break; 
		//visitedTracks[curCand.m_id] = 1;
		// visitedTracks[potCm] = 1;
	      } else {
		potCm = goodNode->m_cm[0];
		}	*/	
	   
	    //else
	    // potCm = -1;

	    // goodIdx = gr.Find(goodId);
	    curCand.insertNewNode(gr, goodNode, k == 0? curCand.m_memberList->begin(): curCand.m_memberList->end());
	    visited[goodId] = 4;


	    //   curLayer = goodNode->m_Layer;

	    next.clear();

	    for(int i = 0; i < goodNode->m_neighbors.size(); i++){
	      
	      int neighId = goodNode->m_neighbors[i];
	      //  debug("Node %d has one neig %d, visited %d", goodId, neighId, visited[neighId]);

	      if(curCand.isInCandidate(neighId)) continue;

	      

	      int neighIdx = gr.Find(neighId);
	      GridNode *neighNode = &Ingrid[neighIdx];

	      next.push_back(neighId);
		
	    }// for neighbors

	    if(goodNode->parent != -1){
	      next.push_back(goodNode->parent);
	      debug("Also adding this node %d", goodNode->parent );
	    }


	    if(next.size() > 0){	      
	      prevId = goodId;
	      prevNode = goodNode;
	      virt.clear();
	      if(goodNode->m_Layer == 0){
		error("Reached end ?");
		curCand.m_finished = 2;
		cond = false;
	      }			 
	    }

	    else {
	      
	      info("we have no more neighbors");
	      cond = false;
	      
	      /*if(potCm != -1){
		
		debug("Last node added belonged to a cm %d, let's see if we should merge them", potCm);
		
		curCand.m_toMergeHead.push_back(potCm);

		const auto p = std::find_if(tracklets.begin(), tracklets.end(), [potCm](const PathCandidate *obj){ return obj->m_id == potCm; } );
		

		if(goodId == (*p)->m_headNode){
		  ((*p)->m_toMergeHead).push_back(curCand.m_id);
		  debug("Merging head to head");
		}
		else if(goodId == (*p)->m_tailNode){
		  ((*p)->m_toMergeTail).push_back(curCand.m_id);
		  debug("Merging head to tail");
		}else{
		  std::vector<int>::iterator it = std::find(((*p)->m_memberList)->begin(), ((*p)->m_memberList)->end(), goodId);
		  int index = std::distance(((*p)->m_memberList)->begin(), it);
		  if(index >= ((*p)->m_memberList)->size()/2){
		    debug("Merging head to head (to check)");
		    ((*p)->m_toMergeHead).push_back(curCand.m_id) ;
		  } else{
		    debug("Merging head to tail (to check)");
		    ((*p)->m_toMergeTail).push_back(curCand.m_id) ;
		  }
		}
		//	visitedTracks[curCand.m_id] = 1;
		//	visitedTracks[potCm] = 1;
	      } else {
		info("No potential neighbors ?");
		curCand.m_finished = 2;
		}*/
		     
	    } // next,size() == 0
			  	      
	    //   cond = false;
	  } // WHILE COND
	} // for k 2
      } // for k trackets
         
    } // IF TRUE

    



    if(merging){
      // Merging potential candidates
      info("Starting Merging phase\n");
      for(unsigned int l = 0; l < tracklets.size(); l++){
      
	PathCandidate &curCand = *(tracklets[l]);
	debug("\n\n");
	debug("Current tracklet %d", curCand.m_id);
      
	if (curCand.m_finished ==  3){
	  info("This CM is either finished or has no further close neighbors");
	  //check that it has not merging candidates
	  if (tracklets[l]->m_toMergeHead.size() != 0  || tracklets[l]->m_toMergeTail.size() != 0)
	    debug("This looks finished but has to be merged ???? ");
	  else
	    debug("No potential merging candidates, we're done");
	  curCand.m_isValid = 1;	
	} else if (curCand.m_finished ==  2){
	  info("This CM needs to be checked");
	  curCand.m_isValid = 1;
	} else if (curCand.m_isMerged == 1){
	  info("This CM has already been merged");
	  curCand.m_isValid = 0;
	}

	else {
 
	  int sizeMergeHead =  curCand.m_toMergeHead.size();
	  int sizeMergeTail =  curCand.m_toMergeTail.size();

	  if(sizeMergeHead == 0 && sizeMergeTail == 0)
	    error("NEW CAND FOR NOTHIIIING");
	

	  PathCandidate *newCand 	= new PathCandidate();// Create a new candidate
	  newCand->m_id 		= candidateId++;// Set id
	  newCand->m_tailNode           = curCand.m_tailNode;
	  
	  for(int i = 0; i < (curCand.m_memberList)->size(); i++){
	    int curid = (curCand.m_memberList)->at(i);
	    int curidx = gr.Find(curid);
	    GridNode* node = &Ingrid[curidx];
	    newCand->insertNewNode(gr,node,newCand->m_memberList->end());
	  }

	  int curCandId =  curCand.m_id;

	  if(sizeMergeHead != 0){

	    int idToMerge = curCand.m_toMergeHead[0];
	    bool continu = true;
	  
	    while(continu){
	    
	      debug("Tracklets %d needs to be merged in head with %d", curCand.m_id, idToMerge);
	      const auto p = std::find_if(tracklets.begin(), tracklets.end(),
					  [idToMerge](const PathCandidate *obj){ return obj->m_id == idToMerge; } );
	    
	      PathCandidate &mergeCand = *(*p);
	  
	      if(std::find(mergeCand.m_toMergeHead.begin(), mergeCand.m_toMergeHead.end(), curCandId) != mergeCand.m_toMergeHead.end()) {
	      
		debug("Merging head with head");
		addTracklets(gr, newCand, mergeCand, 1, 1);
	    
		if(mergeCand.m_toMergeTail.size() > 0){
		  curCandId = mergeCand.m_id;
		  idToMerge = mergeCand.m_toMergeTail[0];
		  debug("ONE MORE TO PUSH %d (not implemented yet)", idToMerge);
		} else
		  continu = false;

	      
	      } else {
	      
		debug("Merging head with tail");	    
		addTracklets(gr, newCand, mergeCand, 1, 0);
	    
		if(mergeCand.m_toMergeHead.size()){
		  curCandId = mergeCand.m_id;

		  idToMerge = mergeCand.m_toMergeHead[0];
		  debug("ONE MORE TO PUSH %d (not implemented yet)", idToMerge);
		} else
		  continu = false;

	      }
	    } // while continu
	  } // HEAD MERGING

	  if(sizeMergeTail != 0){

	    int idToMerge = curCand.m_toMergeTail[0];
	    bool continu = true;
	  
	    while(continu){

	      debug("Tracklets %d needs to be merged in tail with %d", curCand.m_id, idToMerge);

	      const auto p = std::find_if(tracklets.begin(), tracklets.end(),
					  [idToMerge](const PathCandidate *obj){ return obj->m_id == idToMerge; } );
	    
	      PathCandidate &mergeCand = *(*p);
	  
	      if(std::find(mergeCand.m_toMergeHead.begin(), mergeCand.m_toMergeHead.end(), curCandId) != mergeCand.m_toMergeHead.end()) {
	      
		debug("Merging tail with head");	    
		addTracklets(gr, newCand, mergeCand, 0, 1);
	    

		if(mergeCand.m_toMergeTail.size() > 0 ){
		  curCandId = mergeCand.m_id;
		  idToMerge = mergeCand.m_toMergeTail[0];
		  debug("ONE MORE TO PUSH %d (not implemented yet)", idToMerge);
		} else
		  continu = false;
	      
	      } else {
	      
		debug("Merging head with tail");	    
		addTracklets(gr, newCand, mergeCand, 0, 0);
	    
		if(mergeCand.m_toMergeHead.size() > 0){
		  curCandId = mergeCand.m_id;
		  idToMerge = mergeCand.m_toMergeHead[0];
		  debug("ONE MORE TO PUSH %d (not implemented yet)", idToMerge);
		} else
		  continu = false;

	      }
	    } // while continu
	  } // Tail MERGING
	  
	  //	else
	  //	  error("ISSUE, no one to merge with ?");
	
	  curCand.m_isValid = 0;
	  curCand.m_isMerged = 1;
	  newCand->m_isValid = 1;
	  newCand->m_finished = 2;
	
	  info("Pushing new merged cm %d: \n               length is %d, \n     tail node %d \n     head node %d \n       Min layer %d, \n              Max layer %d, \n               IsOnSectorLimit %d, \n               Finished ? %d. ", newCand->m_id, newCand->m_length, newCand->m_tailNode, newCand->m_headNode,newCand->m_minLayer, newCand->m_maxLayer, newCand->m_isOnSectorLimit, newCand->m_finished);

	  tracklets.push_back(newCand);
	} // ELSE NOT FINISHED

      } // FOR TRACKLETS
      // }
    }

    //fitZCoordinates(gr, tracklets);

    for(unsigned int l = 0; l < tracklets.size(); l++){
      PathCandidate *curCand = tracklets[l];
      std::set<int> const *trk = curCand->m_memberIdSet;
      std::vector<int> const *vect = curCand->m_memberList;

       if(curCand->m_isValid) {
	 // fitZCoordinates(gr, curCand);
	 std::set<int> *comp = new std::set<int>((*trk));	    
	 connectedComp->push_back(comp);	 
      }
    }
    
    int NumConnComp = connectedComp->size();
    info("Number of connected components: %d", NumConnComp);    
    std::vector<TrackObject*>* MVDMergedTraks = MergeConnectedComponentsWithMVD(gr, connectedComp);
    TrackZ_CoordinatesDistNorm(gr, MVDMergedTraks);

    //__________________ Determind eth Z-coordinate values.
    
    ComponentPerEvt.Fill(k, NumConnComp);


    int cm =0;
    for(unsigned int l = 0; l < tracklets.size(); l++){
      PathCandidate &curCand = *(tracklets[l]);
      std::set<int> const *trk = curCand.m_memberIdSet;
      std::vector<int> const *vect = curCand.m_memberList;
      std::vector<double> &x =  curCand.m_x;
      std::vector<double> &y =  curCand.m_y;
      std::vector<double> &z =  curCand.m_z; 

      std::vector<double> &r =  curCand.m_r;
      std::vector<double> &theta =  curCand.m_theta; 

      if(curCand.m_isValid) {
	for( int i = 0;  i < vect->size(); ++i) {
	  int detID = vect->at(i);// Id of the current detector
	  //	printf("CM %d, id %d \n", cm, detID);
	  int d_Index = gr.Find(detID);// Index in the grid
	  GridNode  &node = Ingrid[d_Index];
	  //  debug("%lf, %lf, %lf", x[i], y[i], z[i]);
	  //if(node.m_type!= GridNode::STT_TYPE_SKEW)
	  // .. ConnectedCoord.Fill(k, cm, detID, node.m_x, node.m_y, node.m_z , node.m_r,node.m_thetaDeg,
	  //			x[i], y[i], 0);
	//else
	    ConnectedCoord.Fill(k, cm, node.m_detID, node.m_x, node.m_y, node.m_z, node.m_r,node.m_thetaDeg,
				x[i], y[i], node.m_z_Det);
	  /*	for( int i = 0;  i < vect->size(); ++i) {
		int detID = vect->at(i);// Id of the current detector
		//	printf("CM %d, id %d \n", cm, detID);
		int d_Index = gr.Find(detID);// Index in the grid
		GridNode  &node = Ingrid[d_Index];	
		if(curCand.m_x[i] == -1)
		ConnectedCoord.Fill(k, cm, node.m_detID, node.m_x, node.m_y, node.m_z, node.m_r,node.m_thetaDeg,
		node.m_xDet, node.m_yDet, node.m_z_Det);
		else
		ConnectedCoord.Fill(k, cm, node.m_detID, node.m_x, node.m_y, node.m_z, node.m_r,node.m_thetaDeg,
		curCand.m_x[i], curCand.m_y[i], node.m_z_Det);
		}*/
	}
	cm++;
		  
      }
    }
    // Store the data for each constructed component
    /* for(size_t cm = 0 ; cm < connectedComp->size(); ++cm) {
       std::set<int> const* idset = connectedComp->at(cm);
      if(!idset){
	continue;
      }
      std::set<int>::iterator it;
      for( it = idset->begin(); it != idset->end(); ++it) {
     	int detID = *it;// Id of the current detector
    	int d_Index = gr.Find(detID);// Index in the grid
     	GridNode  &node = Ingrid[d_Index];	
	// k = event number, cm = component number,....
	
	//	ConnectedCoord.Fill(k, cm, node.m_detID, node.m_x, node.m_y, node.m_z,
	//		    node.m_xDet, node.m_yDet, node.m_z_Det);
	//printf("%d, %d, %d, %f, %f, %f, %f, %f, %f \n", k, cm, node.m_detID, node.m_x, node.m_y, node.m_z,
	//		    node.m_xDet, node.m_yDet, node.m_z_Det);


      }
    // Print Info
    }*/

     // ------------------------------------------------------------------------
    #if (WRITE_CM_ASCII_FILE > 0)  
    // Write to an ASCII file
    //////////// Comments to write to file
    std::string header = "cm,x,y,z,mx,my,mz\n";
    //////// End of comments.
    std::ofstream OutTxtFile;
    OutTxtFile.open ("ConnectedComponents.csv");
  
    if (OutTxtFile.is_open()) {
      OutTxtFile << header;
      for(unsigned int l = 0; l < tracklets.size(); l++){
	PathCandidate &curCand = *(tracklets[l]);
	std::vector<int>  *trk = curCand.m_memberList;
	for(size_t cm = 0 ; cm < trk->size(); ++cm) {
	  int detID = trk->at(cm);// Id of the current detector
	  //	printf("CM %d, id %d \n", cm, detID);
	  int d_Index = gr.Find(detID);// Index in the grid
	  GridNode  &node = Ingrid[d_Index];
	  OutTxtFile << l << "," << node.m_x <<"," << node.m_y <<"," << node.m_z <<"," <<	  node.m_xDet<<"," << node.m_yDet<<"," << node.m_z_Det<< '\n';
	  // k = event number, cm = component number,....
	  //	ConnectedCoord.Fill(k, cm, node.m_detID, node.m_x, node.m_y, node.m_z,
	  //			    node.m_xDet, node.m_yDet, node.m_z_Det);
	  //printf("%d, %d, %d, %f, %f, %f, %f, %f, %f \n", k, cm, node.m_detID, node.m_x, node.m_y, node.m_z,
	  //		    node.m_xDet, node.m_yDet, node.m_z_Det);


	}
	// Print Info
      }
      // Writ
    }

  
  
    //OutTxtFile.close();
    #endif

    CollectGridToTree(gr, coord);

    if(connectedComp != 0) {
      for(size_t c = 0; c < connectedComp->size(); ++c) {
        delete connectedComp->at(c);
      }
      delete connectedComp;
    }

  }
  timer.Stop();
  ComponentPerEvt.Write();
  ConnectedCoord.Write();


  
  // Write coordinates ntuple
  coord.Write();
  Out_Put_File.Close();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  std::cout <<"=======================================\n"
            << "Macro finished succesfully.\n"
            << "Real time " << (rtime/totalnumEvt)
	    << " (s/Event), CPU time " << (ctime/totalnumEvt)
            << " (s/Event).\n"
            << '\n';
  }
  // Extract MC track values
  
}
