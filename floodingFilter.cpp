
/*************************************
 * Author: S. Gazagnes               *
 * Version:                          *
 * License:                          *
 *************************************/

#include <iostream>
#include <set>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include<unistd.h>

// Root headers
#include "TFile.h"
#include "TNtuple.h"
#include "TStopwatch.h"
#include "TH1.h"

// Local headers
#include "simon_functions.h"
#include "CollectSttMvdPoints.h"
#include "hitcoordinate.h"
#include "utilfunctions.h"
#include "pathopen.h"
#include "reconstruction.h"
#include "error.h"
#include "logc.h"
#include "phconnect.h"
#include "phfitting.h"
#include "phmerging.h"
#include "phzinterp.h"

// DEBUG AND STORE definitions
#define EVALUATE_ERROR 1
#define READ_GRID_FROM_FILE 0
#define DO_RECONSTRUCTION 1
#define DO_CONNECT 1
#define DO_FITTING 1
#define DO_MERGING 1
#define DO_ZRECONS 1
#define WRITE_CONNECTED_COMPONENTS 1
#define INCLUDE_MVD_INOUTPUT_TRACK 0
#define WRITE_CONNECTED_COMPONENTS_JSON 0
#define WRITE_CM_ASCII_FILE 1


void floodingFilter(std::string const &OutFileName,int firstEvt, int lastEvt)
{
  
  TStopwatch timer;

  //Setting verbosity level, put 1 for the debug you want
  //error/time/info/collect/grid/connect/fit/merge/trkz/trkerror
  bool v[10] = {1,1,1,0,0,0,0,0,0,0};//{0,0,0,0,0,0,0,0,1,1}{1,1,1,0,0,0,0,0,0,0};
  set_verbosity(v);
  // Structure to hold the detector data (grid)
  std::vector < GridNode > detNodes;
  // File and structure to hold the output coordinates.
  TFile Out_Put_File(OutFileName.c_str(),"RECREATE","Outputfile Created by performFilter", 9);

  // Collected coordinates for tracks.
  TNtuple coord ("CoordCollected" , "Collected Coordinates in x y plane", "x:y:z:x_Det:y_Det:z_Det");

  // Ntuple to hold Error values for all events available in the
  // current events set. The value is evaluated per image.
  std::string errorParameter = "nMC:nComp:Error_underMerge:Error_overMerge:TotalError";
  errorParameter += ":Error_underMergeNorm:Error_overMergeNorm:TotalErrorNorm";
  // Create Ntuple to hold parameters.
  TNtuple ErrorNtuple("ErrorEstimate","Segmentation error values", errorParameter.c_str());

  // Second error type. Per track error value. Based on curvature data.
  std::string PerTrkErrPars = "complex:misMatched:bestIdx";
   PerTrkErrPars += ":MCLength";
  PerTrkErrPars += ":Jacardsingle:Jacardaverage";
  PerTrkErrPars += ":UnderMergeError:OverMergeError:disX:disY:disZ";
  TNtuple ErrorNtuplePerTrack("PerTrackError","Per track values of error", PerTrkErrPars.c_str());

  std::string CurvTrak = "MC_a:MC_b:MC_r:MC_E:tr_a:tr_b:tr_r:tr_E";

  std::string DisTrak = "disx:disy:disz";

  TNtuple ErrorNtupleDisPerTrack("DisPerTrackError","Per track values of displacement",DisTrak.c_str());

  // PerTrkErrPars += ":MCMinCurrentLength:CurrentMinMCLength";
  TNtuple CurvNtuplePerTrack("PerTrackCurv","Per track values of circle fit", CurvTrak.c_str());
  // NTuple to hold the coordinates of all connected components.
  std::string ConnCompPar = "EvtNum:CompNum:bestIdx:tubeId:x:y:z:r:thetaDeg:x_Det:y_Det:z_Det";
  TNtuple ConnectedCoord ("ConnectedComponents", "Connected component Coordinates", ConnCompPar.c_str());
   std::string AnchorCCPar = "EvtNum:CompNum:x_Det:y_Det:z_Det";
  TNtuple AnchorCCCoord ("ConnectedComponentsAnchors", "Connected component anchors Coordinates", AnchorCCPar.c_str());
  // Hold number of components per event
  TNtuple ComponentPerEvt ("ComponentPerEvt", "Component per event","evtNum:numComponents");

  //geo 2 1572086365 // geo 1 1583944737 // Muon_z0 1611761510 // Muon_z30 1611844116 // Muon_z120 1611844771 /// 1000 1613554871 // 20000: 1614788215

  /* Read all data directly from sim, digi and parameter files */
  std::vector < std::vector<HitCoordinate*>* >* Hit_coords = 
    CollectSttMvdPoints(detNodes, "../rootfiles/evtcomplete20000",Out_Put_File, 1614788215, firstEvt, lastEvt);
  //evtmuonz120
  std::vector< std::vector < MCTrackObject* >* > *MC_Tracks = MCTrackPoints(*Hit_coords);
  
  // Write event info to output  
  WriteEventPlotsToFile( (*Hit_coords), Out_Put_File);
  
  // Create an empty grid object
  CoordGrid gr;
  
  // Init Grid for STT detector nodes (fill the map).
  gr.Initialize(detNodes);
  gr.CorrectLayerLimit();
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

  dbggrid("Fix neighbouring before extension.");
  
  fixNeighboring(gr);
  std::vector < GridNode > VNodesLayer;
  std::vector < GridNode > VNodesSector;

  Add_VirtualNodes(gr, VNodesLayer, VNodesSector);

  // Compute_Virtual_InterSector_Nodes(gr, 6,VNodes);
  TNtuple* virtualTubesLayer = GridToNtuple(VNodesLayer, "VirtualNodesLayer");
  TNtuple* virtualTubesSector = GridToNtuple(VNodesSector, "VirtualNodesSector");

  virtualTubesLayer->SetMarkerStyle(8);
  virtualTubesLayer->SetMarkerSize(0.2);
  virtualTubesLayer->SetMarkerColor(kMagenta);
  virtualTubesLayer->Write();

  virtualTubesSector->SetMarkerStyle(8);
  virtualTubesSector->SetMarkerSize(0.2);
  virtualTubesSector->SetMarkerColor(kCyan);
  virtualTubesSector->Write();
  
  /* Extend the grid with new virtual points and fix the missing neighboring relations   */

  dbggrid("Extending the grid by %d virtual nodes between the layers.", VNodesLayer.size());
  dbggrid("Extending the grid by %d virtual nodes between the sectors.", VNodesSector.size());

  gr.ExtendedGrid(VNodesLayer);
  gr.ExtendedGrid(VNodesSector);

  dbggrid("Fixing neighbouring after extension.");
  fixNeighboring(gr);
  
  TNtuple* extendedGrid = GridToNtuple(gr.m_grid, "ExtendedGrid");
  extendedGrid->SetMarkerStyle(7);//8
  extendedGrid->SetMarkerSize(0.3);
  extendedGrid->SetMarkerColor(17);//41
  extendedGrid->Write();

  dbggrid("Total number of tubes after extension = %d", gr.GetNumNodes());
  // Delete allocated memory
  delete OrigGrid;
  delete virtualTubesLayer;
  delete virtualTubesSector;
  delete extendedGrid;

  
  unsigned int totalnumEvt = Hit_coords->size();
  info("There will be %u event(s) to process", totalnumEvt);
  // Start the timer.

  //  CoordGrid grCopy;

#if(DO_RECONSTRUCTION == 1)

  //  std::vector< std::set<int>* >* connectedComp = 0;
  //  connectedComp = new std::vector< std::set<int>* >();

  timer.Start();
  auto t1 = std::chrono::high_resolution_clock::now();

  /* Fill the grid with the current hit points and process. Handles   each event separately.*/
  // Event loop
  int nEvproc = 0;

  for(size_t k = 0; k < Hit_coords->size(); ++k) {
    //   if(k == 28 || k == 90 || k == 42||k == 55||k==97) continue;

    if(MC_Tracks->at(k)->size() == 0){
      info("This event did not contain anay tracks");
      continue;
    }
    
    // info("MC_Tracks->at(k)->size() is %d", MC_Tracks->at(k)->size());
  
    // Data for the current event
    info("Processing event: %d", k);
    std::vector<HitCoordinate*> const *dd = 0;
    dd = Hit_coords->at(k);
    if(dd) 
      gr.FillGrid(*dd);
    
    //   CoordGrid grCopy (gr);


 
    /* Pushing all active detectors into queue */
    timing("Fill grid phase ended. Time %lf s",
	   std::chrono::duration<double>( std::chrono::high_resolution_clock::now() - t1 ).count());

    std::vector< GridNode > &Ingridori = gr.m_grid;
    std::vector< GridNode > Ingrid(Ingridori);  
    std::vector< std::set<int>* >* connectedComp = 0;
    connectedComp = new std::vector< std::set<int>* >();
    std::vector< int > activeId, remaining;
    std::vector<pair<int, unsigned short>> idToProcess;

    int nactiveAll = 0, nactiveReal = 0;

    for(size_t n = 0; n < Ingrid.size(); ++n) {
      GridNode &curNode = Ingrid[n];
      if(curNode.m_active){
	int NodeId = Ingrid[n].m_detID;
	//info("ACTIVE %d", NodeId);
	if(curNode.m_type != GridNode::VIRTUAL_NODE  ) {
	  activeId.push_back(NodeId);
	  nactiveReal++;
	  /* Keep only active neighbors  */
	  for  ( int i = 0; i < curNode.m_neighbors.size();){
	    int neigh_ID = curNode.m_neighbors[i];
	    int neigh_index = gr.Find(neigh_ID);
	    if(Ingrid[neigh_index].m_type == GridNode::VIRTUAL_NODE){
	      int nID = Ingrid[neigh_index].m_neighbors[0] == NodeId? Ingrid[neigh_index].m_neighbors[1]:
		Ingrid[neigh_index].m_neighbors[0];
	      GridNode const &First  = Ingrid[gr.Find(nID)];
	      if(First.m_active)
		Ingrid[neigh_index].m_active = true;
	      else{
		(curNode.m_neighbors).erase((curNode.m_neighbors).begin()+i);
		continue;
	      }
		
	    } else if(!Ingrid[neigh_index].m_active){
	      (curNode.m_neighbors).erase((curNode.m_neighbors).begin()+i);
	      continue;
	    }
	    i++;
	  }
	  //pair<int, unsigned short> p1 = {NodeId, (unsigned short) curNode.m_neighbors.size()};
	  idToProcess.push_back(make_pair(NodeId, (unsigned short) curNode.m_neighbors.size()));
	}	
	nactiveAll++;
      }      
    }

    if(nactiveReal <= 10){
      info("This event did not contain any tracks, continue");
      gr.ResetGrid();
      continue;
    } else {
      nEvproc++;
    }
    

    //   remaining = activeId;
    dbgconnect("Found %d active detectors (%d with virtuals)", nactiveReal, nactiveAll);

    sort(idToProcess.begin(), idToProcess.end(), sortbysec);

    //Find coomplex sectors (OLD not used)
    // std::vector< int > sectorToCheck;
    //    complexSectors(gr, activeId, &sectorToCheck);
    //  info("Found %d complex sectors", sectorToCheck.size());

    //Find complex tracks
    std::vector< int > idComplex;  

    #if(EVALUATE_ERROR == 1)
    //info("Finding complex tracks");
    complexTracks(gr, MC_Tracks->at(k), &idComplex);
    #endif


    
    
    std::vector < PathCandidate* > tracklets;
    int candidateId 	= 0;
    char *visited 	= (char *) calloc(Ingrid[Ingrid.size()-1].m_detID+1, sizeof(char));

    dbgconnect("First step is finding the obvious tracklets");
    
#if(DO_CONNECT == 1)
    findEasyTracks (gr, Ingrid, tracklets, idToProcess, visited, &candidateId);
    timing("Connect phase ended. Time %lf s",  std::chrono::duration<double>( std::chrono::high_resolution_clock::now() - t1 ).count());
#endif

    if(tracklets.size() > 0)
      dbgconnect("We found %d tracklets", candidateId);
    // else{
    //   error("FIRST PHASE DID NOT WORK, NEED TO ABORT");
    //   break;
    //  }

    size_t sizeBef = idToProcess.size();
    for(unsigned int n = 0; n < idToProcess.size();) {
      int curId  = idToProcess[n].first;
      if(visited[curId] == 1){
	idToProcess.erase(idToProcess.begin() + n);
	continue;
      }
      n++;
    }
    
    dbgconnect("We matched %lu detectors (Remaining %lu)", sizeBef-idToProcess.size(),idToProcess.size() );


    std::sort(tracklets.begin(), tracklets.end(), compareTwoPathsLength);

    
    /* ++++++++++++++++++++++++++++++++++++++++++++++++++++ */
      
    /*               Let's do some fitting !!!!             */
      
    /* ++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#if(DO_FITTING)
    dbgfit("Starting fitting phase");

    int **toMergeWith =  (int **) calloc(tracklets.size(), sizeof(int*));
    for (size_t i =0; i < tracklets.size(); i++)
      toMergeWith[i] = (int *) calloc(tracklets.size(), sizeof(int));
    
    fittingPhase(gr, Ingrid, tracklets, idToProcess, visited, toMergeWith);


    for(size_t n = 0; n < idToProcess.size(); ++n) {
      if(visited[idToProcess[n].first] == 1){
	idToProcess.erase(idToProcess.begin() + n);
	n--;
      }
    }
    
    dbgfit("After fitting, we have %lu remaining active nodes\n",  idToProcess.size());
    timing("Fitting phase ended. Time %lf s",  std::chrono::duration<double>( std::chrono::high_resolution_clock::now() - t1 ).count());


    // Trying to connect easy remaining nodes

    
    sort( idToProcess.begin(), idToProcess.end() );
	
    for(size_t n = 0; n < idToProcess.size(); ++n) {
      std::vector<int> sameLayer;
      std::vector<int> otherLayer;

      int nodeId = idToProcess[n].first;
      if(visited[nodeId])
	continue;
      GridNode &myNode = Ingrid[gr.Find(nodeId)];
      //  info("Checking remaining node %d, visited %d", nodeId, visited[nodeId]);
      //visited[nodeId] = 3;
      
      bool allvisited = true;
      for(size_t m = 0; m < myNode.m_neighbors.size(); m++){
	int neighId = myNode.m_neighbors[m];
	GridNode &myNeigh = Ingrid[gr.Find(neighId)];
	if(myNeigh.m_type == GridNode::VIRTUAL_NODE)
	  continue;
	//	info("Has neighbor %d, which is connect to %d CC, and is layer %d", neighId, myNeigh.m_cm.size(), myNeigh.m_Layer);
	if(!visited[neighId]){
	  allvisited = false;
	  // visited[neighId] = 3;
	}else if(visited[neighId] == 1){
	  if(myNeigh.m_Layer == myNode.m_Layer){
	    //	    info("Same Layer, put in list (connect with CC %d)", myNeigh.m_cm[0] );
	    if(std::find(sameLayer.begin(), sameLayer.end(), neighId) == sameLayer.end())
	      sameLayer.push_back(neighId);
	  } else{
	    //	    info("Other Layer, put in list (connect with CC %d)", myNeigh.m_cm[0]);
	    if(std::find(otherLayer.begin(), otherLayer.end(), neighId) == otherLayer.end())
	      otherLayer.push_back(neighId);
	  }
	}
      }


      
      if(sameLayer.size() > 0){

	//	info("We have %d nodes on same layer, and all visited, we connect to %d", sameLayer.size(),sameLayer[0] );
	
	GridNode &myNeigh = Ingrid[gr.Find(sameLayer[0])];
	//	info("We have %d nodes on same layer, and all visited, we connect to %d, cm is %d", sameLayer.size(),myNeigh.m_detID, myNeigh.m_cm[0] );

	int potCCtoMerge = myNeigh.m_cm[0];

	//	info("We connect %d to CC %d",nodeId,potCCtoMerge   );

	const auto p = std::find_if(tracklets.begin(), tracklets.end(),
				    [potCCtoMerge](const PathCandidate *obj){ return obj->m_id == potCCtoMerge; } );

	PathCandidate &neighCand = *(*p); // The CC that the node belongs to	   	   

	//Find where the node is in the list of the other CC

	std::vector<int>::iterator it = std::find((neighCand.m_memberList)->begin(),
						  (neighCand.m_memberList)->end(), myNeigh.m_detID);
	//	info("We connect %d to CC %d",nodeId,potCCtoMerge   );

	neighCand.insertNewNode(gr, Ingrid, &myNode, it);
	visited[nodeId] = 1;
	//	info("Connected");
      } else if(allvisited && otherLayer.size() >  0){
	//	info("We have %d nodes on other layers, and all visited", otherLayer.size());
	std::vector<int> CC;

	for(size_t m = 0; m < otherLayer.size(); m++){
	  GridNode &myNeigh = Ingrid[gr.Find(otherLayer[m])];
	  if(std::find(CC.begin(), CC.end(), myNeigh.m_cm[0]) == CC.end()){
	    CC.push_back(myNeigh.m_cm[0]);
	  }
	}
	if(CC.size() == 0){
	  // info("All neighbors are connected to single CC %d", CC[0]);
	  int potCCtoMerge = CC[0];
	  //	  dbgconnect("We connect %d to CC %d",nodeId,potCCtoMerge   );

	  const auto p = std::find_if(tracklets.begin(), tracklets.end(),
				      [potCCtoMerge](const PathCandidate *obj){ return obj->m_id == potCCtoMerge;});

	  PathCandidate &neighCand = *(*p); // The CC that the node belongs to	   	   

	  //Find where the node is in the list of the other CC
	  std::vector<int>::iterator it = std::find((neighCand.m_memberList)->begin(),
						    (neighCand.m_memberList)->end(), otherLayer[0]);
	  neighCand.insertNewNode(gr, Ingrid, &myNode, it);
	  visited[nodeId] = 1;
	} 
      }
    }


    for(size_t n = 0; n < idToProcess.size(); ++n) {
      if(visited[idToProcess[n].first] == 1){
	idToProcess.erase(idToProcess.begin() + n);
	n--;
      }
    }

    
    for(size_t n = 0; n < idToProcess.size(); ++n) {
      std::vector<int> sameLayer;
      std::vector<int> otherLayer;
      std::vector<int> connected; 
      std::vector<GridNode*> inQueue; 

      int nodeId = idToProcess[n].first;
      if(visited[nodeId])
	continue;
      GridNode &myNode = Ingrid[gr.Find(nodeId)];
      //  info("Checking remaining node %d, visited %d", nodeId, visited[nodeId]);
      connected.push_back(nodeId);
      visited[nodeId] = 3;
      for(size_t m = 0; m < myNode.m_neighbors.size(); m++){
	int neighId = myNode.m_neighbors[m];
	GridNode &myNeigh = Ingrid[gr.Find(neighId)];
	if(myNeigh.m_type == GridNode::VIRTUAL_NODE)
	  continue;
	//	info("Has neighbor %d, which is connect to %d CC, and is layer %d", neighId, myNeigh.m_cm.size(), myNeigh.m_Layer);
	if(!visited[neighId]){
	  connected.push_back(neighId);
	  inQueue.push_back(&myNeigh);
	  visited[neighId] = 3;
	}
      }
      while(inQueue.size()>0){
	GridNode *cur = inQueue.back();
	//	info("Checking %d", cur->m_detID);
		
	inQueue.pop_back();
	for(size_t m = 0; m < cur->m_neighbors.size(); m++){
	  int neighId = cur->m_neighbors[m];
	  // info("Neigh %d", neighId);

	  GridNode &myNeigh = Ingrid[gr.Find(neighId)];
	    if(myNeigh.m_type == GridNode::VIRTUAL_NODE)
	      continue;

	    if(!visited[neighId]){
	      connected.push_back(neighId);
	      inQueue.push_back(&myNeigh);
	      visited[neighId] = 3;
	      //   info("Adding neighbor %d", neighId, myNeigh.m_cm.size());
	    } else if(visited[neighId] == 1){
	      //info("Neigh %d is connected to CC  %d",neighId, myNeigh.m_cm[0]);
	    }
	  }
	}
	if(connected.size() > 5){
	  //  info("we found %d nodes connected, set a new cand", connected.size());
	  sort( connected.begin(),  connected.end() );
	  PathCandidate *cand 	= new PathCandidate();// Create a new tracklet candidate
	  cand->m_id 		= (candidateId)++;// tracklet id
	  int prevLayer = -1;
	  std::vector<int> virt; 

	  for(size_t m = 0; m < connected.size(); m++){
	    //  info("Contain node %d", connected[m]);
	    int curId = connected[m];
	    GridNode *addNode = &Ingrid[gr.Find(curId)];
	    if(virt.size() > 0 && addNode->m_Layer != prevLayer){
	      addNodesToCand (gr, Ingrid, *cand, visited, virt);
	    }
	    
	    visited[curId] = 1;
	    cand->insertNewNode(gr, Ingrid, addNode, cand->m_memberList->end());
	    for(size_t p = 0; p < addNode->m_neighbors.size(); p++){
	      int neighId =  addNode->m_neighbors[p];
	      GridNode &myNeigh = Ingrid[gr.Find(neighId)];
	      if(myNeigh.m_type == GridNode::VIRTUAL_NODE)
		virt.push_back(neighId);
	      
	    }
	    prevLayer = addNode->m_Layer;
	  }
	  int firstId         = cand->m_tailNode;
	  GridNode &firstNode = Ingrid[ gr.Find(firstId)];
	  int lastId          = cand->m_headNode;
	  GridNode &lastNode  = Ingrid[gr.Find(lastId)];
	      
	  if((firstNode.m_LayerLimit == 1 && lastNode.m_LayerLimit == 1)){		 
	    // dbgconnect("track goes through all layers or makes a loop, likily finished");		 
	    cand->m_finished = FINISHED;		 
	  }  else {		 
	    //  dbgconnect("Candidate has no more neighbors, but doesn't seem finished");
	    cand->m_finished = ONGOING;	   
	  }
	  //dbgconnect("Pushing cm %d with length %d, tail node %d, head node %d, first layer %d, last layer %d, IsOnSectorLimit %d, status  %d. ", cand->m_id, cand->m_length, cand->m_tailNode, cand->m_headNode,cand->m_layers[0], cand->m_layers[cand->m_layers.size()-1], cand->m_isOnSectorLimit, cand->m_finished);
	  tracklets.push_back(cand);
	} /*else {
	  info("we found %d nodes connected, not enough for a new cand, let's just connect to closest CC", connected.size());

	  for(size_t m = 0; m < connected.size(); m++){
	    info("node %d", connected[m]);
	    float mindist = 10000;
	    int goodcc = -1;
	    int goodNeigh = -1;
	    int curId = connected[m];
	    GridNode &addNode = Ingrid[gr.Find(curId)];
	    for(size_t p = 0; p < addNode.m_neighbors.size(); p++){
	      int neighId =  addNode.m_neighbors[p];
	      GridNode &myNeigh = Ingrid[gr.Find(neighId)];
	      if(myNeigh.m_type == GridNode::VIRTUAL_NODE || visited[neighId] != 1)
		continue;
	      double currDist = distanceBetweenTube(addNode, myNeigh);
	      if(currDist < mindist){
		goodcc = myNeigh.m_cm[0];
		goodNeigh = neighId;
		mindist=currDist;
	      }
	    }

	    info("Best match is %d with CC %d", goodNeigh, goodcc);
	    const auto p = std::find_if(tracklets.begin(), tracklets.end(),
					[goodcc](const PathCandidate *obj){ return obj->m_id==goodcc;});

	    PathCandidate &neighCand = *(*p); // The CC that the node belongs to	   	   

	    //Find where the node is in the list of the other CC
	    std::vector<int>::iterator it = std::find((neighCand.m_memberList)->begin(),
						      (neighCand.m_memberList)->end(), goodNeigh);
	    neighCand.insertNewNode(gr, Ingrid, &addNode, it);
	    visited[nodeId] = 1;
	  }
	}
	}*/
    }


    dbgfit("Ended the connecttion of remaining nodes");

    std::sort(tracklets.begin(), tracklets.end(), compareTwoPathsLength);

    using nodeDist = std::pair<int, float>;



    
    for(unsigned int l = 0; l < tracklets.size(); l++){
     
      PathCandidate &curCand = *(tracklets[l]);
      if(!curCand.m_isValid || curCand.m_finished > 2)
	continue;

      GridNode &firstNode = Ingrid[gr.Find(curCand.m_tailNode)];
      GridNode &lastNode  = Ingrid[gr.Find(curCand.m_headNode)];
      dbgmerge("NEW Tracklet %d is unfinished, firstNode %d, lastNode %d", curCand.m_id, firstNode.m_detID, lastNode.m_detID);

      if(!(firstNode.m_LayerLimit == 1 || curCand.m_toMergeTail.size()> 0 || firstNode.m_neighbors.size() > 0)) { 
	dbgmerge("Let's look into tail direction, with the %lu tracklets we found previously", tracklets.size());
	std::vector<nodeDist> toCheck;
	for(unsigned int n = 0; n < tracklets.size(); ++n) {
	  PathCandidate &testCand = *(tracklets[n]);
	  if (testCand.m_finished == 3 || n == l) continue;

	  GridNode checkNode;
	  if(labs(firstNode.m_detID - testCand.m_tailNode) < labs(firstNode.m_detID - testCand.m_headNode)){
	    checkNode = Ingrid[gr.Find(testCand.m_tailNode)];
	  }else{
	    checkNode = Ingrid[gr.Find(testCand.m_headNode)];
	  }

	  if(labs(checkNode.m_Sector - firstNode.m_Sector) > 1)
	    continue;
	  dbgmerge("Testing with node %d from tracklet %d",checkNode.m_detID, testCand.m_id);

	  double currDist = sqrt(pow(firstNode.m_x- checkNode.m_x,2) +pow(firstNode.m_y- checkNode.m_y,2)) ;
	  dbgmerge("Distance %lf", currDist);
	  if(currDist < 10.){
	    toCheck.push_back(make_pair(checkNode.m_detID,currDist));
	   } else
	    dbgmerge("Too far, no possible merging");
	  
	}

	dbgmerge("Checking now \n");

	if(toCheck.size() > 0){
	  // info("Found nodes to Cgeck");

	  std::sort(toCheck.begin(), toCheck.end(), [](nodeDist const &a, nodeDist const &b) { 
						      return a.second < b.second;	  });
	  
	  for(size_t i = 0; i < MIN(toCheck.size(),5); i++){
	    GridNode &check = Ingrid[gr.Find(toCheck[i].first)];
	    int potCCtoMerge = check.m_cm[0];
	      
	    const auto p = std::find_if(tracklets.begin(), tracklets.end(),
					[potCCtoMerge](const PathCandidate *obj){ return obj->m_id == potCCtoMerge; } );

	    PathCandidate &testCand = *(*p); // The CC that the node belongs to
	    dbgmerge("Check CC %d", potCCtoMerge);

	    int caseMerge      = 0;
	    int offAnc         = 0;
	    if(check.m_detID == testCand.m_tailNode && testCand.m_toMergeTail.size() == 0){
	      dbgmerge("Possibility in tail, testing angle");
	      offAnc      = 0;
	      caseMerge   = 0;
	      GridNode &prevAnc = curCand.m_anchors[1];
	      dbgmerge("PrevAncho %d, %f, %f", prevAnc.m_detID,prevAnc.m_xDet, prevAnc.m_yDet);
	      dbgmerge("Next neigh Anc %d, %f, %f", testCand.m_anchors[offAnc].m_detID,testCand.m_anchors[offAnc].m_xDet, testCand.m_anchors[offAnc].m_yDet);

	      float angle_xy = returnAngle(prevAnc.m_xDet, firstNode.m_xDet, testCand.m_anchors[offAnc].m_xDet,
	      		   prevAnc.m_yDet, firstNode.m_yDet, testCand.m_anchors[offAnc].m_yDet);
	      dbgmerge("Angle xy with track %f", angle_xy);


	      float dot = (firstNode.m_xDet - prevAnc.m_xDet)*(testCand.m_anchors[offAnc+1].m_xDet - testCand.m_anchors[offAnc].m_xDet) +(firstNode.m_yDet - prevAnc.m_yDet)*(testCand.m_anchors[offAnc+1].m_yDet - testCand.m_anchors[offAnc].m_yDet);    //# dot product between [x1, y1] and [x2, y2]
	      float det = (firstNode.m_xDet - prevAnc.m_xDet)*(testCand.m_anchors[offAnc+1].m_yDet - testCand.m_anchors[offAnc].m_yDet) -(firstNode.m_yDet - prevAnc.m_yDet)*(testCand.m_anchors[offAnc+1].m_xDet - testCand.m_anchors[offAnc].m_xDet); //x1*y2 - y1*x2    //  # determinant
	      float angle = atan2(det, dot)* 180 / 3.14;
	      dbgmerge("Vector Angle xy with track %f", angle);

	      // # atan2(y, x) or atan2(sin, cos)
	      if(fabs(angle) < 60){
	      	dbgmerge("We should merge %d and %d \n", curCand.m_id, testCand.m_id);
		curCand.m_toMergeTail.push_back(potCCtoMerge);
		testCand.m_toMergeTail.push_back(curCand.m_id);		  		    
		//	toMergeWith[curCand.m_id][potCCtoMerge] = caseMerge;
		curCand.m_finished = 2;
		testCand.m_finished = 2;

	      	break;
	      }
	    } else if(check.m_detID == testCand.m_headNode && testCand.m_toMergeHead.size() == 0){
	      dbgmerge("Possibility in head, testing angle");
	      offAnc      = testCand.m_anchors.size()-1;
	      caseMerge   = 1;
	      GridNode &prevAnc = curCand.m_anchors[1];
	      dbgmerge("PrevAncho %d, %f, %f", prevAnc.m_detID,prevAnc.m_xDet, prevAnc.m_yDet);
	      dbgmerge("Next neigh Anc %d, %f, %f", testCand.m_anchors[offAnc].m_detID,testCand.m_anchors[offAnc].m_xDet, testCand.m_anchors[offAnc].m_yDet);

	      float angle_xy = returnAngle(prevAnc.m_xDet, firstNode.m_xDet, testCand.m_anchors[offAnc].m_xDet,
					   prevAnc.m_yDet, firstNode.m_yDet, testCand.m_anchors[offAnc].m_yDet);
	      dbgmerge("Angle xy with track %f", angle_xy);

	      float dot = (firstNode.m_xDet - prevAnc.m_xDet)*(testCand.m_anchors[offAnc-1].m_xDet - testCand.m_anchors[offAnc].m_xDet) +(firstNode.m_yDet - prevAnc.m_yDet)*(testCand.m_anchors[offAnc-1].m_yDet - testCand.m_anchors[offAnc].m_yDet);    //# dot product between [x1, y1] and [x2, y2]
	      float det = (firstNode.m_xDet - prevAnc.m_xDet)*(testCand.m_anchors[offAnc-1].m_yDet - testCand.m_anchors[offAnc].m_yDet) -(firstNode.m_yDet - prevAnc.m_yDet)*(testCand.m_anchors[offAnc-1].m_xDet - testCand.m_anchors[offAnc].m_xDet); //x1*y2 - y1*x2    //  # determinant
	      float angle = atan2(det, dot)* 180 / 3.14;
	      if(fabs(angle) < 60){
		dbgmerge("We should merge %d and %d \n", curCand.m_id, testCand.m_id);
		curCand.m_toMergeTail.push_back(potCCtoMerge);
		testCand.m_toMergeHead.push_back(curCand.m_id);		  		    
		//	toMergeWith[curCand.m_id][potCCtoMerge] = caseMerge;
		curCand.m_finished = 2;
		testCand.m_finished = 2;

		break;
	      }
	    }
	  }
	} else {
	  dbgmerge("No good candidate found");
	}
      }

      

      
      if(!(lastNode.m_LayerLimit == 1 || curCand.m_toMergeHead.size()> 0 || lastNode.m_neighbors.size() > 0)) {
	dbgmerge("Let's look into head direction, with the %lu tracklets we found previously", tracklets.size());
	std::vector<nodeDist> toCheck;

	for(unsigned int n = 0; n < tracklets.size(); ++n) {
	  PathCandidate &testCand = *(tracklets[n]);
	  if (testCand.m_finished == 3 || n == l) continue;
	
	  GridNode checkNode;
	  if(labs(lastNode.m_detID - testCand.m_tailNode) < labs(lastNode.m_detID - testCand.m_headNode)){
	    checkNode = Ingrid[gr.Find(testCand.m_tailNode)];
	  }else{
	    checkNode = Ingrid[gr.Find(testCand.m_headNode)];
	  }
	  
	  dbgmerge("Testing with node %d from tracklet %d",checkNode.m_detID, testCand.m_id);

	  double currDist = sqrt(pow(lastNode.m_x- checkNode.m_x,2) +pow(lastNode.m_y- checkNode.m_y,2)) ;
	  if(labs(checkNode.m_Sector - lastNode.m_Sector) > 1 && labs(checkNode.m_Sector - lastNode.m_Sector) != 5 )
	    continue;
	  dbgmerge("Distance %lf", currDist);

	  if(currDist < (double) 10.){
	    toCheck.push_back(make_pair(checkNode.m_detID,currDist));
	  } else
	    dbgmerge("Too far, no possible merging");
	  
	}
	if(toCheck.size() > 0){
	  //  info("Found nodes to Cgeck");

	  std::sort(toCheck.begin(), toCheck.end(), [](nodeDist const &a, nodeDist const &b) { 
						      return a.second < b.second;	  });
	  
	  for(size_t i = 0; i < MIN(toCheck.size(),5); i++){
	    GridNode &check = Ingrid[gr.Find(toCheck[i].first)];
	    int potCCtoMerge = check.m_cm[0];
	      
	    const auto p = std::find_if(tracklets.begin(), tracklets.end(),
					[potCCtoMerge](const PathCandidate *obj){ return obj->m_id == potCCtoMerge; } );
	    dbgmerge("Check CC %d", potCCtoMerge);

	    PathCandidate &testCand = *(*p); // The CC that the node belongs to
	    int caseMerge      = 0;
	    int offAnc         = 0;
	    if(check.m_detID == testCand.m_tailNode && testCand.m_toMergeTail.size() == 0){
	      dbgmerge("Possibility in tail, testing angle");
	      offAnc      = 0;
	      caseMerge   = 2;
	      GridNode &prevAnc = curCand.m_anchors[curCand.m_anchors.size()-2];
	      dbgmerge("PrevAncho %d, %f, %f", prevAnc.m_detID,prevAnc.m_xDet, prevAnc.m_yDet);
	      dbgmerge("Next neigh Anc %d, %f, %f", testCand.m_anchors[offAnc].m_detID,testCand.m_anchors[offAnc].m_xDet, testCand.m_anchors[offAnc].m_yDet);

	      float angle_xy = returnAngle(prevAnc.m_xDet, lastNode.m_xDet, testCand.m_anchors[offAnc].m_xDet,
					   prevAnc.m_yDet, lastNode.m_yDet, testCand.m_anchors[offAnc].m_yDet);
	      dbgmerge("Angle xy with track %f", angle_xy);

	      float dot = (lastNode.m_xDet - prevAnc.m_xDet)*(testCand.m_anchors[offAnc+1].m_xDet - testCand.m_anchors[offAnc].m_xDet) +(lastNode.m_yDet - prevAnc.m_yDet)*(testCand.m_anchors[offAnc+1].m_yDet - testCand.m_anchors[offAnc].m_yDet);    //# dot product between [x1, y1] and [x2, y2]
	      float det = (lastNode.m_xDet - prevAnc.m_xDet)*(testCand.m_anchors[offAnc+1].m_yDet - testCand.m_anchors[offAnc].m_yDet) -(lastNode.m_yDet - prevAnc.m_yDet)*(testCand.m_anchors[offAnc+1].m_xDet - testCand.m_anchors[offAnc].m_xDet); //x1*y2 - y1*x2    //  # determinant
	      float angle = atan2(det, dot)* 180 / 3.14;
	      dbgmerge("Vector Angle xy with track %f", angle);
	      if(fabs(angle) < 60){
		dbgmerge("We should merge %d and %d \n", curCand.m_id, testCand.m_id);
		curCand.m_toMergeHead.push_back(potCCtoMerge);
		testCand.m_toMergeTail.push_back(curCand.m_id);		  		    
		//toMergeWith[curCand.m_id][potCCtoMerge] = caseMerge;
		curCand.m_finished = 2;
		testCand.m_finished = 2;
		break;
	      }
	    } else if(check.m_detID == testCand.m_headNode && testCand.m_toMergeHead.size() == 0){
	      dbgmerge("Possibility in head, testing angle");
	      offAnc      = testCand.m_anchors.size()-1;
	      caseMerge   = 1;
	      GridNode &prevAnc = curCand.m_anchors[curCand.m_anchors.size()-2];
	      dbgmerge("PrevAncho %d, %f, %f", prevAnc.m_detID,prevAnc.m_xDet, prevAnc.m_yDet);
	      dbgmerge("Next neigh Anc %d, %f, %f", testCand.m_anchors[offAnc].m_detID,testCand.m_anchors[offAnc].m_xDet, testCand.m_anchors[offAnc].m_yDet);

	      float angle_xy = returnAngle(prevAnc.m_xDet, lastNode.m_xDet, testCand.m_anchors[offAnc].m_xDet,
					   prevAnc.m_yDet, lastNode.m_yDet, testCand.m_anchors[offAnc].m_yDet);
	      dbgmerge("Angle xy with track %f", angle_xy);

	      float dot = (lastNode.m_xDet - prevAnc.m_xDet)*(testCand.m_anchors[offAnc-1].m_xDet - testCand.m_anchors[offAnc].m_xDet) +(lastNode.m_yDet - prevAnc.m_yDet)*(testCand.m_anchors[offAnc-1].m_yDet - testCand.m_anchors[offAnc].m_yDet);    //# dot product between [x1, y1] and [x2, y2]
	      float det = (lastNode.m_xDet - prevAnc.m_xDet)*(testCand.m_anchors[offAnc-1].m_yDet - testCand.m_anchors[offAnc].m_yDet) -(lastNode.m_yDet - prevAnc.m_yDet)*(testCand.m_anchors[offAnc-1].m_xDet - testCand.m_anchors[offAnc].m_xDet); //x1*y2 - y1*x2    //  # determinant
	      float angle = atan2(det, dot)* 180 / 3.14;
	      dbgmerge("Vector Angle xy with track %f", angle);
	      if(fabs(angle) < 60){
		dbgmerge("We should merge %d and %d \n", curCand.m_id, testCand.m_id);
		curCand.m_toMergeHead.push_back(potCCtoMerge);
		testCand.m_toMergeHead.push_back(curCand.m_id);		  		    
		//toMergeWith[curCand.m_id][potCCtoMerge] = caseMerge;
		curCand.m_finished = 2;
		testCand.m_finished = 2;

		break;
	      }
	    }
	  }
	} else {
	  dbgmerge("No good candidate found");
	}

      }

      //if(curCand.m_finished < 2){
      //	dbgmerge("We found no good candidate to merge with");
      //	curCand.m_finished = 3;
     // }
	  
     // dbgmerge("\n\n\n");
    }
   
 #endif // IF TRUE

    



#if(DO_MERGING)
      // Merging potential candidates
      dbgmerge("Starting Merging phase\n");
      mergeTracks (gr, Ingrid, tracklets, &candidateId);
      // timing("Merging phase ended. Real time %f s, CPU time %f s.", timer.RealTime(),timer.CpuTime());
      timing("Merging phase ended. Time %lf s",  std::chrono::duration<double>( std::chrono::high_resolution_clock::now() - t1 ).count());

#endif

#if (DO_ZRECONS)  
      //CompZCoordinates(gr, curCand);
      ZCoordinates(gr, Ingrid, tracklets);
#endif
      
    //fitZCoordinates(gr, tracklets);

      //    if(tracklets.size() > 0){
      for(unsigned int l = 0; l < tracklets.size(); l++){
	PathCandidate &curCand = *(tracklets[l]);
	if(!curCand.m_isValid){
	  tracklets.erase(tracklets.begin() + l);
	  l--;
	  continue;
	}
      }

      info("Number of valid tracklets: %d", tracklets.size());    

      std::vector < MCTrackObject* > *mcTracksCurrentEvent = MC_Tracks->at(k);
      // SOrt MC tracks increasing length
      std::sort(mcTracksCurrentEvent->begin(), mcTracksCurrentEvent->end(), greaterThanLength);
      std::vector<int> matchedId = BestCompIdToMCTracks( mcTracksCurrentEvent, &tracklets);
      
      for(unsigned int l = 0; l < tracklets.size(); l++){
	PathCandidate *curCand = tracklets[l];
	std::set<int> const *trk = curCand->m_memberIdSet;
	  std::vector<int> const *vect = curCand->m_memberList;

	  if(curCand->m_isValid) {

	    std::set<int> *comp = new std::set<int>((*trk));	    
	    connectedComp->push_back(comp);
	    // info("Best match id for track %d is %d",l, matchedId[l]);
	  }
	}
    
	int NumConnComp = connectedComp->size();
	//info("Number of connected components: %d", NumConnComp);    
	std::vector<TrackObject*>* MVDMergedTraks = MergeConnectedComponentsWithMVD(gr, connectedComp);
	// TrackZ_CoordinatesDistNorm(gr, MVDMergedTraks);
	//	timer.Stop();
  
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
	    for( size_t i = 0;  i < vect->size(); ++i) {
	      int detID = vect->at(i);// Id of the current detector
	      //	printf("CM %d, id %d \n", cm, detID);
	      int d_Index = gr.Find(detID);// Index in the grid
	      GridNode  &node = Ingrid[d_Index];
	      //  debug("%lf, %lf, %lf", x[i], y[i], z[i]);
	      //if(node.m_type!= GridNode::STT_TYPE_SKEW)
	      // .. ConnectedCoord.Fill(k, cm, detID, node.m_x, node.m_y, node.m_z , node.m_r,node.m_thetaDeg,
	      //			x[i], y[i], 0);
	      //else
	      ConnectedCoord.Fill(k, cm, matchedId[l], node.m_detID, node.m_x, node.m_y, node.m_z,
				  node.m_r,node.m_thetaDeg,  x[i],y[i],z[i]);

	    }
	    for( size_t i = 0;  i < curCand.m_anchors.size(); ++i) {
	      //if(curCand.m_anchors[i].m_weight == 1)
		AnchorCCCoord.Fill(k, cm,  curCand.m_anchors[i].m_xDet,curCand.m_anchors[i].m_yDet,curCand.m_anchors[i].m_z_Det);
	    }
	    cm++;
		  
	  }
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


	#if(EVALUATE_ERROR)
    
	dbgtrkerror("Error for all sectors");

	// Determine the segmentation error. Based on total area.
	//	if(connectedComp->size() > 0){
	MCMatchingError *MC_match_error =  MatchMCTracksWithConnectedComponents(MC_Tracks->at(k), connectedComp);
	  
	  
	ErrorNtuple.Fill(MC_match_error->NumberOfMCTracks,
			 MC_match_error->NumberOfTracklets,
			 MC_match_error->Error_underMerge,
			 MC_match_error->Error_overMerge,
			 MC_match_error->TotalError,
			 MC_match_error->Error_underMergeNorm,
			 MC_match_error->Error_overMergeNorm,
			 MC_match_error->TotalErrorNorm);
	delete MC_match_error;
	  
	  //	} 

	dbgtrkerror("Error per track");

	//std::vector < MCTrackObject* > *mcTracksCurrentEvent = MC_Tracks->at(k);
	// SOrt MC tracks increasing length
	//	std::sort(mcTracksCurrentEvent->begin(), mcTracksCurrentEvent->end(), greaterThanLength);
	std::vector< MCMatchingErrorStruct* > *match_error2 =
	  MatchPerTrackWithMCTracks(gr, mcTracksCurrentEvent, &tracklets, idComplex, matchedId);
	if(match_error2 != 0) {
	  for(size_t f = 0; f < match_error2->size(); ++f) {
	    MCMatchingErrorStruct const *erObj = match_error2->at(f);
	    ErrorNtuplePerTrack.Fill(static_cast<float>(erObj->Complex),
				     static_cast<float>(erObj->isNotmatched),
				     static_cast<float>(erObj->matchIndex),
				     static_cast<float>(erObj->BestMatchMCLength),
				     // static_cast<float>(erObj->CurrentTrackLength),
				     //  static_cast<float>(erObj->MCMinCurrentLength),
				     //   static_cast<float>(erObj->CurrentMinMCLength),
				     erObj->Jacardsingle,
				     erObj->Jacardaverage,
				     erObj->Error_underMerge,
				     erObj->Error_overMerge,
				     static_cast<float>(erObj->disX),
				     static_cast<float>(erObj->disY),
				     static_cast<float>(erObj->disZ));
	    CurvNtuplePerTrack.Fill( erObj->MC_a, erObj->MC_b, erObj->MC_r,
				     erObj->MC_E, erObj->tr_a, erObj->tr_b,
				     erObj->tr_r, erObj->tr_E);
	    if(!erObj->isNotmatched){
	      for(size_t p = 0; p < erObj->alldisx.size(); p++){
		ErrorNtupleDisPerTrack.Fill(erObj->alldisx[p],erObj->alldisy[p],erObj->alldisz[p]);
	      }
	    }
	  }//
	  // Clean memory for now. Maybe better to put all lists in a main
	  // list and fill the ntuple later.(HINT FIXME later)
	  for(size_t r = 0; r < match_error2->size(); ++r) {
	    delete match_error2->at(r);
	  }
	  delete match_error2;
	}
	//}

	#endif

      if(connectedComp != 0) {
	for(size_t c = 0; c < connectedComp->size(); ++c) {
	  delete connectedComp->at(c);
	}
	delete connectedComp;
      }

      CollectGridToTree(gr, coord);
      
      gr.ResetGrid();
  }
  #endif
  
  timer.Stop();


  ComponentPerEvt.Write();
  ConnectedCoord.Write();
  AnchorCCCoord.Write();


  
  // Write coordinates ntuple
  coord.Write();

  
  //Write error estimations.
  ErrorNtuple.Write();
  ErrorNtuplePerTrack.Write();
  ErrorNtupleDisPerTrack.Write();

  CurvNtuplePerTrack.Write();

  Out_Put_File.Close();

  for(size_t i = 0; i < Hit_coords->size(); ++i) {
    std::vector<HitCoordinate*>* dd = (*Hit_coords)[i];
    for(size_t j = 0; j < dd->size(); ++j) {
      delete dd->at(j);
    }
    delete (*Hit_coords)[i];
  }
  Hit_coords->clear();
 
  for(size_t j = 0; j < MC_Tracks->size(); ++j) {
    std::vector < MCTrackObject* >* MCtracksDel = MC_Tracks->at(j);
    for(size_t k = 0; k < MCtracksDel->size(); ++k) {
      delete MCtracksDel->at(k);
    }
    delete (MC_Tracks->at(j));
  }
  delete MC_Tracks;

  info("Number of events processed is %d", nEvproc);
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  timing("Macro finished successfully. Time %lf s",  std::chrono::duration<double>( std::chrono::high_resolution_clock::now() - t1 ).count());
  timing("Real time %f (s/Event), CPU time %f (s/Event).", (rtime), (ctime));

  timing("Real time %f (s/Event), CPU time %f (s/Event).", (rtime/nEvproc), (ctime/nEvproc));
}
  // Extract MC track values
  

