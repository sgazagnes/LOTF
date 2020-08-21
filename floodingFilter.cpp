
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
  std::string ConnCompPar = "EvtNum:CompNum:tubeId:x:y:z:x_Det:y_Det:z_Det";
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
  Compute_Add_VirtualNodes_Neigbor(gr, VNodes);

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
  Fix_InterSector_Nodes(gr, 6);

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



    std::vector < PathCandidate* > tracklets;
    int candidateId 	= 0;
    char *visited 	= (char *) calloc(5000, sizeof(int));

    info("First step, let's find the obvious tracks");

    
    for(unsigned int n = 0; n < nactiveReal; ++n) {
      
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

	cond = sortNeighbors(gr, currentNode, prevLayer, sameLayer, nextLayer, nextVirt, visited,  &dir);
	
	if(cond == false){
	  // info("Cond is false, we go to next node");
	  resetLists(visited, prevLayer, sameLayer, nextLayer);
	  continue;
	}
	
	PathCandidate *cand 	= new PathCandidate();// Create a new candidate
	cand->m_id 		= candidateId++;// Set id
	cand->m_tailNode 	= curId;
	visited[curId] 		= 4;

	cand->insertNewNode(gr,currentNode, cand->m_memberList->end());
	prevNodes.push_back(curId);


	while(cond){

       	  
	  if(n_neighbors == 1){ /* 1 NEIGHBOR */

	    if (dir == UP){
	      v = &nextLayer;
	      //  info("Going up from node %d to node %d", curId, v->at(0));
	    } else if (dir == DOWN){
	      v = &prevLayer;
	      //  info("Going down from node %d to node %d", curId, v->at(0));
	    } else {
	      v = &sameLayer;
	      //  info("Going same level from node %d to node %d", curId, v->at(0));
	    }
	     
	    
	    if(nextVirt.size() > 0){
	      // debug("Handling the virtuals");
	      
	      for (int i = 0; i < nextVirt.size(); i++){
		neighId   = nextVirt[i];
		neighIdx  = gr.Find(neighId);
		neighNode = &Ingrid[neighIdx];
		cand->insertNewNode(gr,neighNode,cand->m_memberList->end());
		visited[neighId] = 4;
		n_connected++;

		removeIdFromNeigh(neighNode, &prevNodes, curId);
	      }
	      
	      // prevNodes.clear();
	      prevNodes.insert(prevNodes.end(),  nextVirt.begin(),  nextVirt.end());
	      nextVirt.clear();
	    }
	      
	    

	    neighId    = v->at(0);
	    neighIdx   = gr.Find(neighId);
	    neighNode  = &Ingrid[neighIdx];
	    cand->insertNewNode(gr,neighNode,cand->m_memberList->end());
	    visited[neighId] = 4;
	    n_connected++;
	    
	    removeIdFromNeigh(neighNode, &prevNodes, curId);

	    curId       = neighId;
	    curIdx      = neighIdx;
	    currentNode = neighNode;
	    curLayer    = currentNode->m_Layer;
	    
	    nextLayer.clear();  sameLayer.clear();   prevLayer.clear();
	    prevNodes.clear();
	    prevNodes.push_back(curId);
	    dir = 0;

	    cond = sortNeighbors(gr, currentNode, prevLayer, sameLayer, nextLayer, nextVirt, visited,  &dir);

	    n_neighbors = sameLayer.size() + prevLayer.size() + nextLayer.size();
	    info("%d nodes were connected,  %d found for the next step \n", n_connected, n_neighbors);
	    n_connected = 0;

	    
	  }


	  /* 1 Same layer neighbor */

	  else if (sameLayer.size() > 0){

	    v = dir == UP ? &nextLayer: &prevLayer;
	    

	    if( !areAdjacent(gr, v) ){
	      
	      //info("Not all further nodes are adjacent, we must stop");
	      cond = false;
	      
	    } else { 


	      int candId 	 = sameLayer[0];
	      int candIdx        = gr.Find(candId);
	      GridNode *candNode = &Ingrid[candIdx];
	      curLayer           = candNode->m_Layer;
	      sameLayer.clear();

	      info("Still on the same layer, investigating node %d", candId);

	      //  bool toadd = true;
	      // int cursize = v->size();

	    
	      for(int i = 0; i < candNode->m_neighbors.size(); i++){
		neighId = candNode->m_neighbors[i];	      
		neighIdx  = gr.Find(neighId);
		neighNode = &Ingrid[neighIdx];
		
		if (neighId == curId || neighNode->m_type == GridNode::VIRTUAL_NODE) continue;
			      
		int haveNeigh = 0;
		
		for (int j = 0; j < v->size(); j++){
		  
		  //	  debug("Connection between %d and %d ?", neighId, v->at(j));

		  if(v->at(j) == neighId || (neighNode->IsNeighboring(v->at(j))))
		    haveNeigh = 1;
		  // else 
		  // debug("%d and %d are not (obviously) connected", v->at(j), neighId);
		  
		}

		
		if(haveNeigh == 0){
		  //  debug("Not neighbor found for %d", neighId);
		  cond = false;
		  break;
		}
		
	      } // for in neighbors
		
	      if (cond == true){
		
		info("All neighbors look good, let's insert this one !");

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
		info("%d nodes were connected,  %d found for the next step \n\n", n_connected, n_neighbors);
		n_connected = 0;
		
	      }
	      
	      else {

		info("Neighbors not connected... \n\n", n_connected, n_neighbors);
		if(visited[candId] != 4)
		   visited[candId] = 0;
		   
		cand->m_headNeigh.push_back(candId);
		
	      }
	    }
	  }



	  /* MORE NEIGHBOOOORS */


	  else if (n_neighbors > 1) {

	    if (dir == UP){	      
	      //   info("Next nodes are up and we have %d of them", nextLayer.size());
	      v = &nextLayer;
	    }  else if (dir == DOWN){
	      //   info("NExt nodes are down and we have %d of them", prevLayer.size());
	      v = &prevLayer;
	    } else
	      error("WHAT IS THE DIRECTION NOW?");
	    
	    
	    if(areAdjacent(gr, v)){
	      
	      info("Adding %d nodes to the CM", v->size());

	      if(nextVirt.size() > 0){
		//	debug("Handling the virtual");
		for (int i = 0; i < nextVirt.size(); i++){
		  neighId   = nextVirt[i];
		  neighIdx  = gr.Find(neighId);
		  neighNode = &Ingrid[neighIdx];
		  cand->insertNewNode(gr, neighNode, cand->m_memberList->end());
		  visited[neighId] = 4;
		  n_connected++;

		  removeIdFromNeigh(neighNode, &prevNodes, curId);		  
		}
		
		//	prevNodes.clear();
		prevNodes.insert(prevNodes.end(),  nextVirt.begin(),  nextVirt.end());
		nextVirt.clear();

	      }
	      
	      std::vector<int> lookneigh;

	      for (int i = 0; i < v->size(); i++){
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
	      
	      nextLayer.clear();
	      sameLayer.clear();
	      prevLayer.clear();
	      prevNodes.clear();
	      
	      dir = 0;
	     
	      curId       = lookneigh[0];
	      curIdx      = gr.Find(curId);
	      currentNode = &Ingrid[curIdx];
	      curLayer    = currentNode->m_Layer;
	      // info("New current node is %d, looking for neighbors of %d nodes", curId, lookneigh.size());
	      
	      n_neighbors = 0;
	      
	      for(int i = 0; i < lookneigh.size(); i++){
		int id         = lookneigh[i];
		int idx        = gr.Find(id);
		GridNode *node = &Ingrid[idx];
		prevNodes.push_back(id);

		cond = sortNeighbors(gr, node, prevLayer, sameLayer, nextLayer, nextVirt, visited,  &dir);
	      }

	      n_neighbors = sameLayer.size()+prevLayer.size()+nextLayer.size();
	      info("%d nodes were connected,  %d found for the next step \n", n_connected, n_neighbors);
	      n_connected = 0;
	     

	    } // IF ADJACENT

	    else {
	      info("Some of these nodes are node adjacent", v->size());
	      cond = false;
	    }
	    
	      
	  } else {


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
	    } 
	    cond = false;
	    
	  }

	  	  
	  if(cond == false){
	    //    info("We are going out of the boucle");

	    for (int i = 0; i < sameLayer.size(); i++)
	      cand->m_headNeigh.push_back(sameLayer[i]);
	      
	    for (int i = 0; i < nextLayer.size(); i++)
	      cand->m_headNeigh.push_back(nextLayer[i]);

	    for (int i = 0; i < prevLayer.size(); i++)
	      cand->m_headNeigh.push_back(prevLayer[i]);
	    
	    resetLists(visited, prevLayer, sameLayer, nextLayer);
	    //  continue;
	  }

	}// WHILE COND
	
	  info("Pushing cm %d: \n               length is %d, \n     tail node %d \n     head node %d \n       Min layer %d, \n              Max layer %d, \n               IsOnSectorLimit %d, \n               Finished ? %d. ", cand->m_id, cand->m_length, cand->m_tailNode, cand->m_headNode,cand->m_minLayer, cand->m_maxLayer, cand->m_isOnSectorLimit, cand->m_finished);

	tracklets.push_back(cand);

	} // for Node
       	     
    }








    

    info("Refinement");

    std::sort(tracklets.begin(), tracklets.end(), compareTwoPathsLength); 
    for(unsigned int l = 0; l < tracklets.size(); l++){
      PathCandidate &curCand = *(tracklets[l]);
      info("Length %d", curCand.m_length);
      
      //  computePathCurvature(gr,curCand);
      //info("%lf", (curCand.m_CurV_par).m_r);
    }

    char *visitedTracks = (char *) calloc(tracklets.size(), sizeof(char));

    if(true){
      
      for(unsigned int l = 0; l < tracklets.size(); l++){
	
	PathCandidate &curCand = *(tracklets[l]);
	debug("Cur track %d Is finished ? %d, length %d", curCand.m_id, curCand.m_finished,curCand.m_length  );
	
	if (curCand.m_finished >= 2 || curCand.m_length < 6 ) continue;
	
	int first =  curCand.m_tailNode;
	int firstIdx = gr.Find(first);
	GridNode &firstNode = Ingrid[firstIdx];

	int last =  curCand.m_headNode;
	int lastIdx = gr.Find(last);
	GridNode &lastNode = Ingrid[lastIdx];
	
	debug("This cm has tail node %d  and head node %d, size head neigh %d tail %d", first, last,curCand.m_headNeigh.size() , curCand.m_tailNeigh.size());


	for(int k = 0; k < 2; k++){
	  std::vector<int> next;
	  int prevId, curLayer, layerCurDiff;
	  GridNode *prevNode;
	  
	  if(k == 1 && curCand.m_headNeigh.size() > 0 && curCand.m_toMergeHead.size() == 0 ){
	  
	    info("HEAD : Finding potential continuation with previous neighbors in head direction ");

	    for(int i = 0; i  < curCand.m_headNeigh.size(); i++)
	      debug("%d", curCand.m_headNeigh[i]);

	    next.insert(next.end(),  (curCand.m_headNeigh).begin(),  (curCand.m_headNeigh).end());
	    prevId = curCand.m_headNode;
	    prevNode = &lastNode;
	    curLayer = lastNode.m_Layer;
	    
	  } else if (k == 0 && curCand.m_tailNeigh.size() > 0 && curCand.m_toMergeTail.size() == 0 ){
	  
	    info("TAIL : Finding potential continuation with previous neighbors in tail direction ");

	    for(int i = 0; i  < curCand.m_tailNeigh.size(); i++)
	      debug("%d", curCand.m_tailNeigh[i]);

	    next.insert(next.end(),  (curCand.m_tailNeigh).begin(),  (curCand.m_tailNeigh).end());
	    prevId = curCand.m_tailNode;
	    prevNode = &firstNode;
	    curLayer = firstNode.m_Layer;

	  } else
	    continue;
	    
	  bool cond = next.size() > 0? true: false;
	  int potCm = -1;
	    
	  std::vector<int>  *trk = curCand.m_memberList;

	  std::vector<double> x;
	  std::vector<double> y;
	  std::vector<int> virt;

	  if(k == 1){
	    for (int i = trk->size() - MIN(trk->size(), 10) ; i < trk->size() ; i++){
	      int id = trk->at(i);
	      int idx = gr.Find(id);
	      GridNode &node = Ingrid[idx];		
	      x.push_back(node.m_xDet);
	      y.push_back(node.m_yDet);
	      if(i == trk->size() -2)
		layerCurDiff = curLayer - node.m_Layer;
	    }
	  } else {
	    for (int i = MIN(trk->size(), 10) ; i >= 0 ; i--){
	      int id = trk->at(i);
	      int idx = gr.Find(id);
	      GridNode &node = Ingrid[idx];		
	      x.push_back(node.m_xDet);
	      y.push_back(node.m_yDet);
	      if(i == trk->size() -1)
		layerCurDiff = curLayer - node.m_Layer;
	    }
	  }

	  // finding previous direction //
	  while (cond){
	    
	    GridNode *goodNode;
	    int goodId = fitNextId(gr, x, y, next, curLayer, layerCurDiff, 1);
	    
	    if (goodId == -1) {
	      info("No good candidates have been found, stop");
	      curCand.m_finished = 2;
	      cond = false;
	      break;
	    }
	    
	      
	    int goodIdx = gr.Find(goodId);
	    goodNode = &Ingrid[goodIdx];

	    if(visited[goodId] == 4){
	      info("This node already belongs to a CM (%d)", goodNode->m_cm[0]);

	      if(potCm == goodNode->m_cm[0]){
		  
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
	      }		
	    } // Visited good id 4
	    else
	      potCm = -1;
	      
	  
	    next.clear();
	    
	    for(int i = 0; i < goodNode->m_neighbors.size(); i++){
	      
	      int neighId = goodNode->m_neighbors[i];
	      if(curCand.isInCandidate(neighId)) continue;

	      debug("Node %d has one neig %d, visited %d", goodId, neighId, visited[neighId]);

	      int neighIdx = gr.Find(neighId);
	      GridNode *neighNode = &Ingrid[neighIdx];

	      if(neighNode->m_type == GridNode::VIRTUAL_NODE){
		virt.push_back(neighId);
		continue;
	      }
		  
	      if(visited[neighId] == 4){
		debug("Node %d has already been connected, let see", neighId);
		// potTracks.push_back(neighNode->m_cm[0]);
	      }

	      next.push_back(neighId);
	      debug("Adding Node %d ", neighId);
		  
		
	    }// for neighbors

	    if(goodNode->parent != -1){
	      next.push_back(goodNode->parent);
	      debug("Also adding this node %d", goodNode->parent );
	    }

	    if(next.size() > 0) {

	      if(virt.size() > 0){
		for(int i = 0; i < virt.size(); i++){
		  if(goodNode->IsNeighboring(virt[i]) && prevNode->IsNeighboring(virt[i])){
		    
		    debug("Node %d is common neighbor, add before", virt[i]);
		    int virtIdx = gr.Find(virt[i]);
		    GridNode *virtNode = &Ingrid[virtIdx];
		    curCand.insertNewNode(gr,virtNode, k == 0? curCand.m_memberList->begin(): curCand.m_memberList->end());
		  }
		}
	      }

		
	      debug("Let's add this new node and continue");
	      
	      curCand.insertNewNode(gr,goodNode, k == 0? curCand.m_memberList->begin(): curCand.m_memberList->end());
	      visited[goodId] = 4;
	      prevId = goodId;
	      layerCurDiff = goodNode->m_Layer - curLayer;
	      curLayer = goodNode->m_Layer;
	      prevNode = goodNode;
	      x.erase(x.begin());
	      y.erase(y.begin());
	      x.push_back(goodNode->m_xDet);
	      y.push_back(goodNode->m_yDet);
	      virt.clear();
	    }

	    else {
	      
	      info("we have no more neighbors");
	      cond = false;
	      
	      if(potCm != -1){
		
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
	      }
		     
	    } // next,size() == 0
		
	  
	      
	    //   cond = false;
	  } // WHILE COND
	} // for k 2
      } // for k trackets
      
     // if true
	




    // Assigning remaining nodes

      /*   for(unsigned int n = 0; n < nactiveReal; ++n) {
	
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

    } // IF TRUE

    
    if(false){

    
      // Merging potential candidates

    for(unsigned int l = 0; l < tracklets.size(); l++){
      
      PathCandidate &curCand = *(tracklets[l]);
      debug("Current tracklet %d", curCand.m_id);
      
      if (curCand.m_finished ==  3){
	info("This CM is either finished or has no further close neighbors");
	//check that it has not merging candidates
	if (tracklets[l]->m_toMergeHead.size() != 0  ||tracklets[l]->m_toMergeTail.size() != 0)
	  debug("This looks finished but has to be merged ???? ");
	else
	  debug("No potential merging candidates, we're done");
	curCand.m_isValid = 1;
	
      }
      else if (curCand.m_finished ==  2){
	info("This CM needs to be checked");
	curCand.m_isValid = 1;
      }

      else if (curCand.m_isMerged == 1){
	info("This CM has already been merged");
	curCand.m_isValid = 0;
      }

      else {
 
	int sizeMergeHead =  curCand.m_toMergeHead.size();
	int sizeMergeTail =  curCand.m_toMergeTail.size();

	if(sizeMergeHead == 0 && sizeMergeTail == 0)
	  error("NEW CAND FOR NOTHIIIING");
	
	//int idToMergeHead = curCand.m_toMergeHead;
	//int idToMergeTail = curCand.m_toMergeTail;

	PathCandidate *newCand 	= new PathCandidate();// Create a new candidate
	newCand->m_id 		= candidateId++;// Set id
	newCand->m_tailNode 	= curCand.m_tailNode;
	//	cand->m_memberList      
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
	    
	      // (mergeCand.m_toMergeHead).erase(std::remove((mergeCand.m_toMergeHead).begin(), (mergeCand.m_toMergeHead).end(), curCandId), (mergeCand.m_toMergeHead).end());

	      if(mergeCand.m_toMergeTail.size() > 0){
		curCandId = mergeCand.m_id;
		idToMerge = mergeCand.m_toMergeTail[0];
		debug("ONE MORE TO PUSH %d", idToMerge);
	      } else
		continu = false;

	      
	    } else {
	      
	      debug("Merging head with tail");
	    
	      addTracklets(gr, newCand, mergeCand, 1, 0);

	    
	      //(mergeCand.m_toMergeTail).erase(std::remove((mergeCand.m_toMergeTail).begin(), (mergeCand.m_toMergeTail).end(),curCand.m_id), (mergeCand.m_toMergeTail).end());
	    
	      if(mergeCand.m_toMergeHead.size()){
		curCandId = mergeCand.m_id;

		idToMerge = mergeCand.m_toMergeHead[0];
		debug("ONE MORE TO PUSH %d", idToMerge);
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
	    
	      // (mergeCand.m_toMergeHead).erase(std::remove((mergeCand.m_toMergeHead).begin(), (mergeCand.m_toMergeHead).end(), curCandId), (mergeCand.m_toMergeHead).end());

	      if(mergeCand.m_toMergeTail.size() > 0 ){
		curCandId = mergeCand.m_id;
		idToMerge = mergeCand.m_toMergeTail[0];
		debug("ONE MORE TO PUSH %d", idToMerge);
	      } else
		continu = false;

	      
	    } else {
	      
	      debug("Merging head with tail");
	    
	      addTracklets(gr, newCand, mergeCand, 0, 0);

	    
	      //(mergeCand.m_toMergeTail).erase(std::remove((mergeCand.m_toMergeTail).begin(), (mergeCand.m_toMergeTail).end(),curCand.m_id), (mergeCand.m_toMergeTail).end());
	    
	      if(mergeCand.m_toMergeHead.size() > 0){
		curCandId = mergeCand.m_id;

		idToMerge = mergeCand.m_toMergeHead[0];
		debug("ONE MORE TO PUSH %d", idToMerge);
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
	newCand->m_finished = 3;
	
	info("Pushing new merged cm %d: \n               length is %d, \n     tail node %d \n     head node %d \n       Min layer %d, \n              Max layer %d, \n               IsOnSectorLimit %d, \n               Finished ? %d. ", newCand->m_id, newCand->m_length, newCand->m_tailNode, newCand->m_headNode,newCand->m_minLayer, newCand->m_maxLayer, newCand->m_isOnSectorLimit, newCand->m_finished);

	tracklets.push_back(newCand);
      } // ELSE NOT FINISHED

    } // FOR TRACKLETS
    // }
    }


  
    for(unsigned int l = 0; l < tracklets.size(); l++){
      PathCandidate &curCand = *(tracklets[l]);
      std::set<int> const *trk = curCand.m_memberIdSet;
      
       if(curCand.m_isValid) {
	 std::set<int> *comp = new std::set<int>((*trk));	    
	 connectedComp->push_back(comp);    
      }
    }
    
    int NumConnComp = connectedComp->size();
    info("Number of connected components: %d", NumConnComp);    
    std::vector<TrackObject*>* MVDMergedTraks = MergeConnectedComponentsWithMVD(gr, connectedComp);

    //__________________ Determind eth Z-coordinate values.
   TrackZ_CoordinatesDistNorm(gr, MVDMergedTraks);
    
    ComponentPerEvt.Fill(k, NumConnComp);
    // Store the data for each constructed component
    for(size_t cm = 0 ; cm < connectedComp->size(); ++cm) {
      std::set<int> const* idset = connectedComp->at(cm);
      if(!idset){
	continue;
      }
      std::set<int>::iterator it;
      for( it = idset->begin(); it != idset->end(); ++it) {
     	int detID = *it;// Id of the current detector
	//	printf("CM %d, id %d \n", cm, detID);
     	int d_Index = gr.Find(detID);// Index in the grid
     	GridNode  &node = Ingrid[d_Index];	
	// k = event number, cm = component number,....
	ConnectedCoord.Fill(k, cm, node.m_detID, node.m_x, node.m_y, node.m_z,
			    node.m_xDet, node.m_yDet, node.m_z_Det);
	//printf("%d, %d, %d, %f, %f, %f, %f, %f, %f \n", k, cm, node.m_detID, node.m_x, node.m_y, node.m_z,
	//		    node.m_xDet, node.m_yDet, node.m_z_Det);


      }
    // Print Info
    }

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
