
/*************************************
 * Author: S. Gazagnes               *
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
#include "auxfunctions.h"
#include "gridNode.h"
#include "CoordGrid.h"
#include "CollectSttMvdPoints.h"
#include "hitcoordinate.h"
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
#define WRITE_CM_ASCII_FILE 0
#define WRITE_TIME_ASCII_FILE 1
#define WRITE_LIST_RECO_ID 0

void floodingFilter(std::string const &OutFileName,int firstEvt, int lastEvt)
{
  
  /* Files and Tuples definitions */
  
  // File and structure to hold the output coordinates.
  TFile OutputFile(OutFileName.c_str(), "RECREATE", "Outputfile created by floodingFilter macro", 9);

  // Collected coordinates for tracks.
  TNtuple GridCoordNtuple("CoordCollected" , "Collected Coordinates in x y plane", "x:y:z:x_Det:y_Det:z_Det");

  // Ntuple to hold Global Error values 
  std::string GlobalErrorPars = "evtID:nMC:nReco:Error_underMerge:Error_overMerge:TotalError";
  GlobalErrorPars += ":Error_underMergeNorm:Error_overMergeNorm:TotalErrorNorm";
  TNtuple GlobalErrorNtuple("GlobalEvtError","Segmentation error values", GlobalErrorPars.c_str());

  // NTuple to hold error per MC track
  std::string PerMCTrkErrPars = "evtID:Complex:Mismatched:MatchIdx:MCLength:RecoLength:InterLength:UnionLength";
  PerMCTrkErrPars += ":F1score:MeanDiffX:MeanDiffY:MeanDiffZ";
  TNtuple ErrorPerMCTrackNtuple("PerMCTrackError", "Per track values of error", PerMCTrkErrPars.c_str());
  
  //NTuple to hold curvature/momentum parameters per pair of MC and Reco track
  std::string CurvTrkPars = "evtID:MC_px:MC_py:MC_pz:MC_r1:MC_r2";
  CurvTrkPars+= ":tr_rIsoRand:tr_rIsoRand16:tr_rIsoRand84:tr_rAnc:tr_rPts:tr_scattAngle";
  TNtuple CurvPerTrackNtuple("PerTrackCurv", "Per track values of circle fit and MC momentum", CurvTrkPars.c_str());

  // NTuple to hold panda error metric per reco track
  std::string PerRecoTrkErrPars = "evtID:Rank:Clone";
  TNtuple ErrorPerRecoTrackNtuple("PerRecoTrackError", "Per track values of error", PerRecoTrkErrPars.c_str());


  //NTuple to hold coordinates differences between MC and Reco
  std::string DiffTrkPars = "disx:disy:disz";
  TNtuple CoordDiffPerTrackNtuple("CoordDiffPerTrackError", "Per track values of coordinates errors",
				 DiffTrkPars.c_str());

  // NTuple to hold the coordinates of all reco'd tracks.
  std::string ConnCompPar = "EvtNum:CompNum:bestIdx:tubeId:x:y:z:r:thetaDeg:x_Det:y_Det:z_Det";
  TNtuple ConnectedCoord ("ConnectedComponents", "Connected component Coordinates", ConnCompPar.c_str());

    // NTuple to hold the coordinates of the anchors of the reco'd tracks.
  std::string AnchorCCPar = "EvtNum:CompNum:x_Det:y_Det:z_Det";   
  TNtuple AnchorCCCoord ("ConnectedComponentsAnchors", "Connected component anchors Coordinates",
			 AnchorCCPar.c_str());
  
  // NTuple to Hold number of components per event
  TNtuple ComponentPerEvt ("ComponentPerEvt", "Component per event","evtNum:numComponents");

  // Ntuple to hold MC tubes IDS and track number
  TNtuple MCDetID("MCDetID", "MC track detectors per event","evtNum:MCtrkID:detID");
  TNtuple RecoDetID("RecoDetID", "Reco'd track detectors per event","evtNum:RecotrkID:detID");

  
  // Setting verbosity level, put 1 for the debug you want
  // {error,time,info,collect,grid,connect,fit,merge,trkz,trkerror}
  // examples: {0,0,0,0,0,0,0,0,1,1}{1,1,1,0,0,0,0,0,0,0};
  bool v[10] = {0,0,0,0,0,0,0,0,0,0};
  set_verbosity(v);
  

  // Reading the parameters of given simulations
  char *SimName = "../rootfiles/evtcomplete20000Beam15";   
  // geo 2 1572086365
  // geo 1 1583944737
  // Muon_z0 1611761510
  // Muon_z30 1611844116
  // Muon_z120 1611844771
  // 20000: 1614788215
  // 20000B15: 1618615353
  // Muon B1 1619780508
  // Beam 3 1619749644

  // Read all data directly from sim, digi and parameter files 
  std::vector < GridNode > DetectorNodes;
  std::vector < std::vector<HitCoordinate*>* >* HitsData = 
    CollectSttMvdPoints(DetectorNodes, SimName, OutputFile, 1618615353, firstEvt, lastEvt);
  
  std::vector< std::vector < MCTrackObject* >* > *MCTracks = MCTrackPoints(*HitsData);
  
  // Write event info to output  
  WriteEventPlotsToFile( (*HitsData), OutputFile);
  
  // Create an empty grid object and init it for STT
  CoordGrid GridStruct;
  GridStruct.Initialize(DetectorNodes);
  GridStruct.CorrectLayerLimit();
  TNtuple Layers("LayerLimits","Layer Limits.","x:y:det_z:z");
  TNtuple Sections("SectionsLimits","Section Limits.","x:y:det_z:z");
  
  // Isolate Sector and Layer Limits
  GridStruct.isolateSectorAndLayerLimits(Sections, Layers);
  Sections.Write();
  Layers.Write();

  // Write Original grid setup to root file
  TNtuple* OrigGrid = GridToNtuple(DetectorNodes, "OrigGridCoord");
  OrigGrid->SetMarkerStyle(7);
  OrigGrid->SetMarkerSize(0.3);
  OrigGrid->SetMarkerColor(17);
  OrigGrid->Write();

  // Fix some neighbors before extending to Virtual Node  
  GridStruct.fixNeighboring();

  // Extend grid with Virtual Nodes
  std::vector < GridNode > VNodesLayer;
  std::vector < GridNode > VNodesSector;
  GridStruct.AddVirtualNodes(VNodesLayer, VNodesSector);

  // Old extension for intersector Nodes: Compute_Virtual_InterSector_Nodes(gr, 6,VNodes);
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

  dbggrid("Extending the grid by %d virtual nodes between the layers.", VNodesLayer.size());
  dbggrid("Extending the grid by %d virtual nodes between the sectors.", VNodesSector.size());

  GridStruct.ExtendedGrid(VNodesLayer);
  GridStruct.ExtendedGrid(VNodesSector);

  // Fix some neighbors after extending to Virtual Node  
  GridStruct.fixNeighboring();

  // Write extended grid
  TNtuple* extendedGrid = GridToNtuple(GridStruct.m_grid, "ExtendedGrid");
  extendedGrid->SetMarkerStyle(7);//8
  extendedGrid->SetMarkerSize(0.3);
  extendedGrid->SetMarkerColor(17);//41
  extendedGrid->Write();

  dbggrid("Total number of tubes after extension = %d", GridStruct.GetNumNodes());
  
  // Clean up
  delete OrigGrid;
  delete virtualTubesLayer;
  delete virtualTubesSector;
  delete extendedGrid;
  VNodesLayer.clear();
  VNodesSector.clear();

  // Finding total number of events collected
  unsigned int totalnumEvt = HitsData->size();
  info("We have %u event(s) to process", totalnumEvt);
  
  // Start the timer.
  auto TimeAll        = std::chrono::high_resolution_clock::now();
  double RecTotalTime = 0;


  /*  +++++++++++++++++++++++++++++++++++++++  */
  /*  ++++++++ Start reconstruction +++++++++  */
  /*  +++++++++++++++++++++++++++++++++++++++  */

  TStopwatch timer;

  #if(DO_RECONSTRUCTION == 1)
  timer.Start();
  int lastevtID = -1;
  int nEvtProc  = 0;
  int nTrkProc  = 0;
  int nTrk5hitsProc  = 0;

  int nHitsProc = 0;

  // Open file to write timings
  #if (WRITE_TIME_ASCII_FILE > 0)  
  std::string header = "evt,trks,hits,filltime,contime,fittime,mergtime,ztime,tottime\n";
  std::ofstream TimeTxtFile;
  TimeTxtFile.open ("Timeperevent.csv");
  if (TimeTxtFile.is_open()) 
    TimeTxtFile << header;
  #endif

  #if (WRITE_LIST_RECO_ID > 0)
  std::string MCheader = "evt,trkid,detid\n";
  std::ofstream MCOutTxtFile;
  MCOutTxtFile.open ("MCTrackIDList.csv");
  if (MCOutTxtFile.is_open()) 
    MCOutTxtFile << MCheader;

  std::string Recoheader = "evt,trkid,detid\n";
  std::ofstream RecOutTxtFile;
  RecOutTxtFile.open ("RecoTrackIDList.csv");
  if (RecOutTxtFile.is_open()) 
    RecOutTxtFile << Recoheader; 
  #endif
  
  // Checking that the event has enough number of hits to process (5 minimum)
  for(size_t evt = 0; evt < HitsData->size(); ++evt) {
    int totHits = 0;
    std::set<int> nIndTubes;
    
    for(size_t i = 0; i < MCTracks->at(evt)->size(); ++i) {
      MCTrackObject const *MCtrack = MCTracks->at(evt)->at(i);
      std::set<int> MCSttComp((MCtrack->m_STT_Component).begin(), (MCtrack->m_STT_Component).end());
      nIndTubes.insert(MCSttComp.begin(), MCSttComp.end());
    }
    
    totHits = nIndTubes.size();

    if(totHits == 0){ // 
      info("This event did not contain any tracks");
      continue;
    } else if( totHits < 5){
      info("This event contains only %d hits, skip", totHits);
      continue;
    } else {
      info("Processing event: %d", evt);
      nHitsProc += totHits;
      nTrkProc  += MCTracks->at(evt)->size();
      for(size_t i = 0; i < MCTracks->at(evt)->size(); ++i) {
	MCTrackObject const *MCtrack = MCTracks->at(evt)->at(i);
	std::set<int> MCSttComp((MCtrack->m_STT_Component).begin(), (MCtrack->m_STT_Component).end());
    
	if(MCtrack->m_STT_Component.size()>5){ // Only considering MC tracks with >5 hits
	  nTrk5hitsProc++;
	}
      }
      nIndTubes.clear();
      nEvtProc++;
      lastevtID = evt;
    }

    
    //Starting Evt clock
    auto TimeEvt = std::chrono::high_resolution_clock::now();   
    double FillingTime;
    double ConnectTime;
    double FittingTime;
    double MergingTime;
    double zRecTime;
    double EvtTotalTime;
    
    int nSTTHits  = (int) totHits; 
    int nMCTracks = (int) MCTracks->at(evt)->size();

    // Filling the grid
    if(HitsData->at(evt)) 
      GridStruct.FillGrid(*(HitsData->at(evt)));
        
    auto TimeNow = std::chrono::high_resolution_clock::now();
    FillingTime  = std::chrono::duration<double>(TimeNow - TimeEvt).count();
    timing("Time elapsed to fill the grid: %.6lf s", FillingTime);

    std::vector< GridNode > &DetGridCopy = GridStruct.m_grid; // Copying the grid to avoid re-init for each event
    std::vector< GridNode > DetGrid(DetGridCopy);              // Starting from a fresh already init grid :)
  
    std::vector< std::set<int>* >* RecoTrackListDetID = 0; // Vector of reconstruction tracks
    RecoTrackListDetID = new std::vector< std::set<int>* >(); 
    
    int nTubesFired = 0;
    
    std::vector<pair<int, unsigned short>> ListActiveTubes;
    std::vector<pair<int, unsigned short>> ListPrioTubes;
    std::vector<pair<int, unsigned short>> ListOtherTubes;
    std::set<int>::iterator it = GridStruct.m_STT_idx.begin();

     while (it != GridStruct.m_STT_idx.end()) {
       int Idx = *it;
       GridNode &CurNode = DetGrid[Idx-1];
       it++;

       for  ( size_t i = 0; i < CurNode.m_neighbors.size() ; ){
	 int NeighTubeIdx = CurNode.m_neighbors[i];
	 
	 if(DetGrid[NeighTubeIdx-1].m_type == GridNode::VIRTUAL_NODE){
	   int nID = DetGrid[NeighTubeIdx-1].m_neighbors[0] == Idx?
	     DetGrid[NeighTubeIdx-1].m_neighbors[1] : DetGrid[NeighTubeIdx-1].m_neighbors[0];
	   GridNode const &First  = DetGrid[nID-1];
	   
	   if(First.m_active)
	     DetGrid[NeighTubeIdx-1].m_active = true;
	   else{
	     (CurNode.m_neighbors).erase((CurNode.m_neighbors).begin()+i);
	     continue;
	   }
	    
	 } else if(!DetGrid[NeighTubeIdx-1].m_active){
	   (CurNode.m_neighbors).erase((CurNode.m_neighbors).begin()+i);
	   continue;
	 }
	 
	 i++;
       } // End For CurNode neighbors
	  
       if(CurNode.m_LayerLimit == 1 )
	 ListPrioTubes.push_back(make_pair(Idx, (unsigned short) CurNode.m_neighbors.size()));
       else if (CurNode.m_neighbors.size() == 1)
	 ListOtherTubes.push_back(make_pair(Idx, (unsigned short) CurNode.m_neighbors.size()));
       else
	 ListActiveTubes.push_back(make_pair(Idx, (unsigned short) CurNode.m_neighbors.size()));
       
       nTubesFired++;         
     } // End while GridStruct

    
    dbgconnect("Found %d active detectors (%d with virtuals)",ListActiveTubes.size(), nTubesFired);
    sort(ListPrioTubes.begin(), ListPrioTubes.end(), sortbysec);

    ListActiveTubes.insert( ListActiveTubes.begin(), ListOtherTubes.begin(), ListOtherTubes.end() );
    ListActiveTubes.insert( ListActiveTubes.begin(), ListPrioTubes.begin(), ListPrioTubes.end() );
    ListPrioTubes.clear();
    ListOtherTubes.clear();

    std::vector< int > ListComplexTracks;  
    #if(EVALUATE_ERROR == 1)
    ComplexTracks(GridStruct, MCTracks->at(evt), &ListComplexTracks);
    #endif
    
    
    std::vector < PathCandidate* > RecoTracks;
    int candidateId  = 0;
    char *visited    = (char *) calloc(DetGrid.size()+1, sizeof(char));


    /* ++++++++++++++++++++++++++++++++++++++++++++++++++++ */
      
    /*        First Reconstruction Phase: Connect           */
      
    /* ++++++++++++++++++++++++++++++++++++++++++++++++++++ */

    #if(DO_CONNECT == 1)
    info("First phase: finding obvious tracklets");
    findEasyTracks (GridStruct, DetGrid, RecoTracks, ListActiveTubes, visited, &candidateId);
    TimeNow     =  std::chrono::high_resolution_clock::now();
    ConnectTime = std::chrono::duration<double>( TimeNow - TimeEvt ).count() - FillingTime;
    timing("Time elapsed in connect phase: %.6lf s",  ConnectTime);
    
    if(RecoTracks.size() > 0)
      dbgconnect("We found %d RecoTracks", candidateId);

    size_t OldnTubes = ListActiveTubes.size();
    
    for(unsigned int n = 0; n < ListActiveTubes.size();) {
      
      if(visited[ListActiveTubes[n].first] == 1){
	ListActiveTubes.erase(ListActiveTubes.begin() + n);
	continue;
      }
      
      n++;
    }
    
    dbgconnect("We matched %lu detectors (remaining %lu)",
	       OldnTubes-ListActiveTubes.size(), ListActiveTubes.size());

    std::sort(RecoTracks.begin(), RecoTracks.end(), compareTwoPathsLength);    
    #endif

    
    /* ++++++++++++++++++++++++++++++++++++++++++++++++++++ */
      
    /*        First Reconstruction Phase: Fitting           */
      
    /* ++++++++++++++++++++++++++++++++++++++++++++++++++++ */

    #if(DO_FITTING)
    // Array to store the IDs of tracks to be merged
    int **TrackIDToMerge =  (int **) calloc(RecoTracks.size(), sizeof(int*));
    for (size_t i = 0; i < RecoTracks.size(); i++)
      TrackIDToMerge[i] = (int *) calloc(RecoTracks.size(), sizeof(int));
    
    info("Second phase: fitting to find complex tracks");
    fittingPhase(GridStruct, DetGrid, RecoTracks, ListActiveTubes, visited, TrackIDToMerge);

    OldnTubes = ListActiveTubes.size();
    for(size_t n = 0; n < ListActiveTubes.size(); ++n) {
      
      if(visited[ListActiveTubes[n].first] == 1){
	ListActiveTubes.erase(ListActiveTubes.begin() + n);
	n--;
      }
      
    }
    
    dbgfit("We matched %lu detectors (remaining %lu)",
	       OldnTubes-ListActiveTubes.size(), ListActiveTubes.size());

    TimeNow     =  std::chrono::high_resolution_clock::now();
    FittingTime = std::chrono::duration<double>( TimeNow - TimeEvt ).count() - ConnectTime - FillingTime;
    timing("Time elapsed in fitting phase: %.6lf s",  FittingTime);



    
    // Trying to connect remaining nodes if possible FIXME because I'm a mess
    sort( ListActiveTubes.begin(), ListActiveTubes.end() );	
    for(size_t n = 0; n < ListActiveTubes.size(); ++n) {
      std::vector<int> SameLayerTubes;
      std::vector<int> OtherLayerTubes;
      int TubeID = ListActiveTubes[n].first;
      
      if(visited[TubeID])
	continue;
      
      GridNode &ThisTube     = DetGrid[TubeID-1];
      bool     AreAllVisited = true;
      
      for(size_t m = 0; m < ThisTube.m_neighbors.size(); m++){
	int NeighTubeID     = ThisTube.m_neighbors[m];
	GridNode &NeighTube = DetGrid[NeighTubeID-1];
	
	if(NeighTube.m_type == GridNode::VIRTUAL_NODE)
	  continue;

	if(!visited[NeighTubeID]){
	  AreAllVisited = false;
	} else if(visited[NeighTubeID] == 1){
	  if(NeighTube.m_Layer == ThisTube.m_Layer){
	    if(std::find(SameLayerTubes.begin(), SameLayerTubes.end(), NeighTubeID) == SameLayerTubes.end())
	      SameLayerTubes.push_back(NeighTubeID);
	  } else{
	    if(std::find(OtherLayerTubes.begin(), OtherLayerTubes.end(), NeighTubeID) == OtherLayerTubes.end())
	      OtherLayerTubes.push_back(NeighTubeID);
	  }
	}
	
      } // End For tube neighbors

      
      if(SameLayerTubes.size() > 0){
	GridNode &NeighTube = DetGrid[SameLayerTubes[0]-1];
	int TrackToMerge    = NeighTube.m_cm[0];

	const auto p = std::find_if(RecoTracks.begin(), RecoTracks.end(),
				    [TrackToMerge](const PathCandidate *obj){ return obj->m_id == TrackToMerge; } );

	PathCandidate &neighCand = *(*p); // The CC that the node belongs to	   	   

	//Find where the node is in the list of the other CC
	std::vector<int>::iterator it = std::find((neighCand.m_memberList)->begin(),
						  (neighCand.m_memberList)->end(),
						  NeighTube.m_detID);

	neighCand.insertNewNode(GridStruct, DetGrid, &ThisTube, it);
	visited[TubeID] = 1;
      } else if(AreAllVisited && OtherLayerTubes.size() >  0){
	std::vector<int> PotTrackList;
	
	for(size_t m = 0; m < OtherLayerTubes.size(); m++){
	  GridNode &NeighTube = DetGrid[OtherLayerTubes[m]-1];
	  
	  if(std::find(PotTrackList.begin(), PotTrackList.end(), NeighTube.m_cm[0]) == PotTrackList.end()){
	    PotTrackList.push_back(NeighTube.m_cm[0]);
	  }
	  
	}
	
	if(PotTrackList.size() == 0){
	  int TrackToMerge = PotTrackList[0];
	  const auto p = std::find_if(RecoTracks.begin(), RecoTracks.end(),
				      [TrackToMerge](const PathCandidate *obj){ return obj->m_id == TrackToMerge;});

	  PathCandidate &neighCand = *(*p); // The PotTrackList that the node belongs to	   	   
	  //Find where the node is in the list of the other PotTrackList
	  std::vector<int>::iterator it = std::find((neighCand.m_memberList)->begin(),
						    (neighCand.m_memberList)->end(), OtherLayerTubes[0]);
	  neighCand.insertNewNode(GridStruct, DetGrid, &ThisTube, it);
	  visited[TubeID] = 1;
	}
	
      } // End IF ELSE IF same Layer size
      
    } // End For List Active Tubes


    for(size_t n = 0; n < ListActiveTubes.size(); ++n) {
      
      if(visited[ListActiveTubes[n].first] == 1){
	ListActiveTubes.erase(ListActiveTubes.begin() + n);
	n--;
      }
      
    }

    
    // Connecting to the tracklets under certain conditions (FIXME ALSO I'm extremely brutal)
    for(size_t n = 0; n < ListActiveTubes.size(); ++n) {
      std::vector<int> SameLayerTubes;
      std::vector<int> OtherLayerTubes;
      std::vector<int> TubesConnect; 
      std::vector<GridNode*> inQueue; 
      int TubeID = ListActiveTubes[n].first;
      
      if(visited[TubeID])
	continue;
      
      GridNode &ThisTube = DetGrid[TubeID-1];
      TubesConnect.push_back(TubeID);
      visited[TubeID] = 3;
      
      for(size_t m = 0; m < ThisTube.m_neighbors.size(); m++){
	int NeighTubeID     = ThisTube.m_neighbors[m];
	GridNode &NeighTube = DetGrid[NeighTubeID-1];
	
	if(NeighTube.m_type == GridNode::VIRTUAL_NODE)
	  continue;

	if(!visited[NeighTubeID]){
	  TubesConnect.push_back(NeighTubeID);
	  inQueue.push_back(&NeighTube);
	  visited[NeighTubeID] = 3;
	}
      }
      
      while(inQueue.size()>0){
	GridNode *cur = inQueue.back();		
	inQueue.pop_back();
	
	for(size_t m = 0; m < cur->m_neighbors.size(); m++){
	  int NeighTubeID = cur->m_neighbors[m];

	  GridNode &NeighTube = DetGrid[NeighTubeID-1];
	  if(NeighTube.m_type == GridNode::VIRTUAL_NODE)
	    continue;

	  if(!visited[NeighTubeID]){
	    TubesConnect.push_back(NeighTubeID);
	    inQueue.push_back(&NeighTube);
	    visited[NeighTubeID] = 3;
	  }/* else if(visited[NeighTubeID] == 1){
	   //info("Neigh %d is TubesConnect to PotTrackList  %d",NeighTubeID, NeighTube.m_cm[0]);
	   }*/
	}
      }
      
      if(TubesConnect.size() > 5){
	sort( TubesConnect.begin(),  TubesConnect.end() );
	PathCandidate *cand 	= new PathCandidate();// Create a new tracklet candidate
	cand->m_id 		= (candidateId)++;// tracklet id
	int prevLayer = -1;
	std::vector<int> virt; 

	for(size_t m = 0; m < TubesConnect.size(); m++){
	  int curId         = TubesConnect[m];
	  GridNode *addNode = &DetGrid[curId-1];
	  
	  if(virt.size() > 0 && addNode->m_Layer != prevLayer){
	    addNodesToCand (GridStruct, DetGrid, *cand, visited, virt);
	  }
	    
	  visited[curId] = 1;
	  cand->insertNewNode(GridStruct, DetGrid, addNode, cand->m_memberList->end());
	  for(size_t p = 0; p < addNode->m_neighbors.size(); p++){
	    int NeighTubeID     =  addNode->m_neighbors[p];
	    GridNode &NeighTube = DetGrid[NeighTubeID-1];
	    
	    if(NeighTube.m_type == GridNode::VIRTUAL_NODE)
	      virt.push_back(NeighTubeID);
	      
	  }
	  prevLayer = addNode->m_Layer;
	}
	
	int firstId         = cand->m_tailNode;
	int lastId          = cand->m_headNode;
	GridNode &firstNode = DetGrid[firstId-1];
	GridNode &lastNode  = DetGrid[lastId-1];
	      
	if((firstNode.m_LayerLimit == 1 && lastNode.m_LayerLimit == 1)){		 
	  cand->m_finished = FINISHED;		 
	}  else {		 
	  cand->m_finished = ONGOING;	   
	}
	  RecoTracks.push_back(cand);
	}

      /// FIX THIS PART (CODE still works well without it)
      /*else {

	  info("we found %d nodes TubesConnect, not enough for a new cand, let's just connect to closest PotTrackList", TubesConnect.size());

	  for(size_t m = 0; m < TubesConnect.size(); m++){
	    info("node %d", TubesConnect[m]);
	    float mindist = 10000;
	    int goodcc = -1;
	    int goodNeigh = -1;
	    int curId = TubesConnect[m];
	    GridNode &addNode = DetGrid[gr.Find(curId)];
	    for(size_t p = 0; p < addNode.m_neighbors.size(); p++){
	      int NeighTubeID =  addNode.m_neighbors[p];
	      GridNode &NeighTube = DetGrid[gr.Find(NeighTubeID)];
	      if(NeighTube.m_type == GridNode::VIRTUAL_NODE || visited[NeighTubeID] != 1)
		continue;
	      double currDist = distanceBetweenTube(addNode, NeighTube);
	      if(currDist < mindist){
		goodcc = NeighTube.m_cm[0];
		goodNeigh = NeighTubeID;
		mindist=currDist;
	      }
	    }

	    info("Best match is %d with PotTrackList %d", goodNeigh, goodcc);
	    const auto p = std::find_if(RecoTracks.begin(), RecoTracks.end(),
					[goodcc](const PathCandidate *obj){ return obj->m_id==goodcc;});

	    PathCandidate &neighCand = *(*p); // The PotTrackList that the node belongs to	   	   

	    //Find where the node is in the list of the other PotTrackList
	    std::vector<int>::iterator it = std::find((neighCand.m_memberList)->begin(),
						      (neighCand.m_memberList)->end(), goodNeigh);
	    neighCand.insertNewNode(gr, DetGrid, &addNode, it);
	    visited[TubeID] = 1;
	  }
	}
	}*/
	//	SameLayerTubes.clear();
	//	OtherLayerTubes.clear();
//	TubesConnect.clear();
	  
    }


    dbgfit("Ended the connecttion of remaining nodes");
    std::sort(RecoTracks.begin(), RecoTracks.end(), compareTwoPathsLength);
    using nodeDist = std::pair<int, float>;
    
    for(unsigned int l = 0; l < RecoTracks.size(); l++){
     
      PathCandidate &curCand = *(RecoTracks[l]);
      if(curCand.m_finished > 2) //!curCand.m_isValid || 
	continue;

      GridNode &firstNode = DetGrid[curCand.m_tailNode-1];
      GridNode &lastNode  = DetGrid[curCand.m_headNode-1];
      dbgmerge("NEW Tracklet %d is unfinished, firstNode %d, lastNode %d",
	       curCand.m_id, firstNode.m_detID, lastNode.m_detID);

      if(!(firstNode.m_LayerLimit == 1 || curCand.m_toMergeTail.size()> 0 || firstNode.m_neighbors.size() > 0)) { 
	dbgmerge("Let's look into tail direction, with the %lu tracklets found previously", RecoTracks.size());
	std::vector<nodeDist> toCheck;
	for(unsigned int n = 0; n < RecoTracks.size(); ++n) {
	  PathCandidate &testCand = *(RecoTracks[n]);
	  if (testCand.m_finished == 3 || n == l) continue;

	  GridNode checkNode;
	  if(labs((int) firstNode.m_detID - testCand.m_tailNode)
	     < labs((int) firstNode.m_detID - testCand.m_headNode)){
	    checkNode = DetGrid[testCand.m_tailNode-1];
	  }else{
	    checkNode = DetGrid[testCand.m_headNode-1];
	  }

	  if(labs((int) checkNode.m_Sector - firstNode.m_Sector) > 1)
	    continue;
	  dbgmerge("Testing with node %d from tracklet %d",checkNode.m_detID, testCand.m_id);

	  double currDist = sqrt(pow(firstNode.m_x- checkNode.m_x,2) +pow(firstNode.m_y- checkNode.m_y,2)) ;
	  dbgmerge("Distance %lf", currDist);
	  if(currDist < 6.){ // Empirical Criterion
	    toCheck.push_back(make_pair(checkNode.m_detID,currDist));
	   } else
	    dbgmerge("Too far, no possible merging");
	  
	}

	dbgmerge("Checking now \n");
	if(toCheck.size() > 0){
	  std::sort(toCheck.begin(), toCheck.end(), [](nodeDist const &a, nodeDist const &b) { 
						      return a.second < b.second;	  });
	  
	  for(size_t i = 0; i < MIN(toCheck.size(),5); i++){
	    GridNode &check = DetGrid[toCheck[i].first-1];
	    int TrackToMerge = check.m_cm[0];
	      
	    const auto p = std::find_if(RecoTracks.begin(), RecoTracks.end(),
					[TrackToMerge](const PathCandidate *obj){ return obj->m_id == TrackToMerge; } );

	    PathCandidate &testCand = *(*p); // The PotTrackList that the node belongs to
	    dbgmerge("Check PotTrackList %d", TrackToMerge);

	    int caseMerge      = 0;
	    int offAnc         = 0;
	    
	    if(check.m_detID == testCand.m_tailNode && testCand.m_toMergeTail.size() == 0){
	      dbgmerge("Possibility in tail, testing angle");
	      offAnc      = 0;
	      caseMerge   = 0;
	      GridNode &prevAnc = curCand.m_anchors[1];
	      dbgmerge("PrevAncho %d, %f, %f", prevAnc.m_detID,prevAnc.m_xDet, prevAnc.m_yDet);
	      dbgmerge("Next neigh Anc %d, %f, %f", testCand.m_anchors[offAnc].m_detID,
		       testCand.m_anchors[offAnc].m_xDet, testCand.m_anchors[offAnc].m_yDet);

	      float angle_xy = returnAngle(prevAnc.m_xDet, firstNode.m_xDet, testCand.m_anchors[offAnc].m_xDet,
	      		   prevAnc.m_yDet, firstNode.m_yDet, testCand.m_anchors[offAnc].m_yDet);
	      dbgmerge("Angle xy with track %f", angle_xy);

	      //# dot product between [x1, y1] and [x2, y2]
	      float dot = (firstNode.m_xDet - prevAnc.m_xDet)*
		(testCand.m_anchors[offAnc+1].m_xDet - testCand.m_anchors[offAnc].m_xDet) +
		(firstNode.m_yDet - prevAnc.m_yDet)*
		(testCand.m_anchors[offAnc+1].m_yDet - testCand.m_anchors[offAnc].m_yDet);
	       //x1*y2 - y1*x2    //  # determinant
	      float det = (firstNode.m_xDet - prevAnc.m_xDet)*
		(testCand.m_anchors[offAnc+1].m_yDet - testCand.m_anchors[offAnc].m_yDet) -
		(firstNode.m_yDet - prevAnc.m_yDet)*
		(testCand.m_anchors[offAnc+1].m_xDet - testCand.m_anchors[offAnc].m_xDet);
	      float angle = atan2(det, dot)* 180 / 3.14;
	      dbgmerge("Vector Angle xy with track %f", angle);

	      if(fabs(angle) < 70){
	      	dbgmerge("We should merge %d and %d \n", curCand.m_id, testCand.m_id);
		curCand.m_toMergeTail.push_back(TrackToMerge);
		testCand.m_toMergeTail.push_back(curCand.m_id);		  		    
		//	TrackIDToMerge[curCand.m_id][TrackToMerge] = caseMerge;
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
	      dbgmerge("Next neigh Anc %d, %f, %f", testCand.m_anchors[offAnc].m_detID,
		       testCand.m_anchors[offAnc].m_xDet, testCand.m_anchors[offAnc].m_yDet);

	      float angle_xy = returnAngle(prevAnc.m_xDet, firstNode.m_xDet, testCand.m_anchors[offAnc].m_xDet,
					   prevAnc.m_yDet, firstNode.m_yDet, testCand.m_anchors[offAnc].m_yDet);
	      dbgmerge("Angle xy with track %f", angle_xy);

	      float dot = (firstNode.m_xDet - prevAnc.m_xDet)*
		(testCand.m_anchors[offAnc-1].m_xDet - testCand.m_anchors[offAnc].m_xDet) +
		(firstNode.m_yDet - prevAnc.m_yDet)*
		(testCand.m_anchors[offAnc-1].m_yDet - testCand.m_anchors[offAnc].m_yDet);  
	      float det = (firstNode.m_xDet - prevAnc.m_xDet)*
		(testCand.m_anchors[offAnc-1].m_yDet - testCand.m_anchors[offAnc].m_yDet) -
		(firstNode.m_yDet - prevAnc.m_yDet)*
		(testCand.m_anchors[offAnc-1].m_xDet - testCand.m_anchors[offAnc].m_xDet); 
	      float angle = atan2(det, dot)* 180 / 3.14;
	      
	      if(fabs(angle) < 70){
		dbgmerge("We should merge %d and %d \n", curCand.m_id, testCand.m_id);
		curCand.m_toMergeTail.push_back(TrackToMerge);
		testCand.m_toMergeHead.push_back(curCand.m_id);		  		    
		//	TrackIDToMerge[curCand.m_id][TrackToMerge] = caseMerge;
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
	dbgmerge("Let's look into head direction, with the %lu RecoTracks we found previously", RecoTracks.size());
	std::vector<nodeDist> toCheck;

	for(unsigned int n = 0; n < RecoTracks.size(); ++n) {
	  PathCandidate &testCand = *(RecoTracks[n]);
	  if (testCand.m_finished == 3 || n == l) continue;
	
	  GridNode checkNode;
	  if(labs((int) lastNode.m_detID - testCand.m_tailNode) <
	     labs((int) lastNode.m_detID - testCand.m_headNode)){
	    checkNode = DetGrid[testCand.m_tailNode-1];
	  }else{
	    checkNode = DetGrid[testCand.m_headNode-1];
	  }
	  
	  dbgmerge("Testing with node %d from tracklet %d",checkNode.m_detID, testCand.m_id);

	  double currDist = sqrt(pow(lastNode.m_x- checkNode.m_x,2) + pow(lastNode.m_y- checkNode.m_y,2));
	  if(labs((int) checkNode.m_Sector - lastNode.m_Sector) > 1 &&
	     labs((int) checkNode.m_Sector - lastNode.m_Sector) != 5 )
	    continue;
	  dbgmerge("Distance %lf", currDist);

	  if(currDist < (double) 6.)
	    toCheck.push_back(make_pair(checkNode.m_detID,currDist));
	  else
	    dbgmerge("Too far, no possible merging");
	  
	}
	
	if(toCheck.size() > 0){
	  std::sort(toCheck.begin(), toCheck.end(), [](nodeDist const &a, nodeDist const &b) { 
						      return a.second < b.second;	  });
	  
	  for(size_t i = 0; i < MIN(toCheck.size(),5); i++){
	    GridNode &check = DetGrid[GridStruct.Find(toCheck[i].first)];
	    int TrackToMerge = check.m_cm[0];
	      
	    const auto p = std::find_if(RecoTracks.begin(), RecoTracks.end(),
					[TrackToMerge](const PathCandidate *obj){ return obj->m_id == TrackToMerge; } );
	    dbgmerge("Check PotTrackList %d", TrackToMerge);

	    PathCandidate &testCand = *(*p); // The PotTrackList that the node belongs to
	    int caseMerge      = 0;
	    int offAnc         = 0;
	    
	    if(check.m_detID == testCand.m_tailNode && testCand.m_toMergeTail.size() == 0){
	      dbgmerge("Possibility in tail, testing angle");
	      offAnc      = 0;
	      caseMerge   = 2;
	      GridNode &prevAnc = curCand.m_anchors[curCand.m_anchors.size()-2];
	      dbgmerge("PrevAncho %d, %f, %f", prevAnc.m_detID,prevAnc.m_xDet, prevAnc.m_yDet);
	      dbgmerge("Next neigh Anc %d, %f, %f", testCand.m_anchors[offAnc].m_detID,
		       testCand.m_anchors[offAnc].m_xDet, testCand.m_anchors[offAnc].m_yDet);

	      float angle_xy = returnAngle(prevAnc.m_xDet, lastNode.m_xDet, testCand.m_anchors[offAnc].m_xDet,
					   prevAnc.m_yDet, lastNode.m_yDet, testCand.m_anchors[offAnc].m_yDet);
	      dbgmerge("Angle xy with track %f", angle_xy);

	      float dot = (lastNode.m_xDet - prevAnc.m_xDet)
		*(testCand.m_anchors[offAnc+1].m_xDet - testCand.m_anchors[offAnc].m_xDet) +
		(lastNode.m_yDet - prevAnc.m_yDet)*
		(testCand.m_anchors[offAnc+1].m_yDet - testCand.m_anchors[offAnc].m_yDet);   
	      float det = (lastNode.m_xDet - prevAnc.m_xDet)*
		(testCand.m_anchors[offAnc+1].m_yDet - testCand.m_anchors[offAnc].m_yDet) -
		(lastNode.m_yDet - prevAnc.m_yDet)*
		(testCand.m_anchors[offAnc+1].m_xDet - testCand.m_anchors[offAnc].m_xDet);
	      float angle = atan2(det, dot)* 180 / 3.14;
	      dbgmerge("Vector Angle xy with track %f", angle);
	      
	      if(fabs(angle) < 70){
		dbgmerge("We should merge %d and %d \n", curCand.m_id, testCand.m_id);
		curCand.m_toMergeHead.push_back(TrackToMerge);
		testCand.m_toMergeTail.push_back(curCand.m_id);		  		    
		//TrackIDToMerge[curCand.m_id][TrackToMerge] = caseMerge;
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
	      dbgmerge("Next neigh Anc %d, %f, %f", testCand.m_anchors[offAnc].m_detID,
		       testCand.m_anchors[offAnc].m_xDet, testCand.m_anchors[offAnc].m_yDet);

	      float angle_xy = returnAngle(prevAnc.m_xDet, lastNode.m_xDet, testCand.m_anchors[offAnc].m_xDet,
					   prevAnc.m_yDet, lastNode.m_yDet, testCand.m_anchors[offAnc].m_yDet);
	      dbgmerge("Angle xy with track %f", angle_xy);

	      float dot = (lastNode.m_xDet - prevAnc.m_xDet)*
		(testCand.m_anchors[offAnc-1].m_xDet - testCand.m_anchors[offAnc].m_xDet) +
		(lastNode.m_yDet - prevAnc.m_yDet)*
		(testCand.m_anchors[offAnc-1].m_yDet - testCand.m_anchors[offAnc].m_yDet); 
	      float det = (lastNode.m_xDet - prevAnc.m_xDet)*
		(testCand.m_anchors[offAnc-1].m_yDet - testCand.m_anchors[offAnc].m_yDet) -
		(lastNode.m_yDet - prevAnc.m_yDet)*
		(testCand.m_anchors[offAnc-1].m_xDet - testCand.m_anchors[offAnc].m_xDet);
	      float angle = atan2(det, dot)* 180 / 3.14;
	      dbgmerge("Vector Angle xy with track %f", angle);
	      
	      if(fabs(angle) < 70){
		dbgmerge("We should merge %d and %d \n", curCand.m_id, testCand.m_id);
		curCand.m_toMergeHead.push_back(TrackToMerge);
		testCand.m_toMergeHead.push_back(curCand.m_id);		  		    
		//TrackIDToMerge[curCand.m_id][TrackToMerge] = caseMerge;
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
    }
   
    #endif // IF FITTING TRUE


    
    /* ++++++++++++++++++++++++++++++++++++++++++++++++++++ */
      
    /*        Third Reconstruction Phase: Merging           */
      
    /* ++++++++++++++++++++++++++++++++++++++++++++++++++++ */



    #if(DO_MERGING)
    info("Third phase: merging tracks");

    for(unsigned int l = 0; l < RecoTracks.size(); l++){
      PathCandidate &curCand = *(RecoTracks[l]);
      curCand.m_isValid =false;
    }
    
    dbgmerge("Starting Merging phase\n");
    mergeTracks (GridStruct, DetGrid, RecoTracks, &candidateId);
    TimeNow     =  std::chrono::high_resolution_clock::now();
    MergingTime = std::chrono::duration<double>(TimeNow - TimeEvt ).count() -
      FittingTime - ConnectTime - FillingTime;

    for(unsigned int l = 0; l < RecoTracks.size(); l++){
      PathCandidate &curCand = *(RecoTracks[l]);
      if(!curCand.m_isValid){
	RecoTracks.erase(RecoTracks.begin() + l);
	l--;
	continue;
      }
    }

    for(unsigned int l = 0; l < RecoTracks.size(); l++){
      PathCandidate &curCand = *(RecoTracks[l]);
      std::vector<int>  curTrk( *(curCand.m_memberList));
      std::sort(curTrk.begin(), curTrk.end());

      std::vector<int>::iterator it;
      it = remove_if(curTrk.begin(), curTrk.end(), bind2nd(greater<int>(), GridStruct.firstVirtIdx-1));
      curTrk.erase(it,curTrk.end());

      for(unsigned int m = l+1; m < RecoTracks.size(); m++){
	PathCandidate &testCand = *(RecoTracks[m]);
	std::vector<int>  testTrk (*(testCand.m_memberList));
	std::sort(testTrk.begin(), testTrk.end());
	it = remove_if(testTrk.begin(), testTrk.end(), bind2nd(greater<int>(), GridStruct.firstVirtIdx-1));
	testTrk.erase(it,testTrk.end());

	std::vector<int> IntersectionList( (curTrk.size() + testTrk.size()), 0 );
	std::vector<int>::iterator it2;
	//	  std::vector<int> IntersectionList( (Cur_Comp_list.size() + MCSttComp.size()), 0 );
	it2 = std::set_intersection(curTrk.begin(), curTrk.end(),
				    testTrk.begin(), testTrk.end(),
				    IntersectionList.begin());
	IntersectionList.resize(it2 - IntersectionList.begin());
	if((float) IntersectionList.size() > 0.66* MIN(curTrk.size(), testTrk.size())){
	  //error("Track %d and track %d, length %ld and %ld, intersection %ld", l,m,curTrk.size(), testTrk.size(),		  IntersectionList.size());
	  if(curTrk.size() < testTrk.size()){
	    RecoTracks.erase(RecoTracks.begin() + l);
	    l--;
	    break;
	  } else {
	    RecoTracks.erase(RecoTracks.begin() + m);
	    m--;
	  }
	}
      }
    }

    dbgmerge("Number of valid RecoTracks: %d", RecoTracks.size());  
    timing("Time elapsed in merging phase: %.6lf s.", MergingTime);
#endif
   

#if (DO_ZRECONS)  
    info("Fourth phase: z reconstruction");
   ZCoordinates(GridStruct, DetGrid, RecoTracks);
   TimeNow  = std::chrono::high_resolution_clock::now();
   zRecTime = std::chrono::duration<double>(TimeNow - TimeEvt).count() -
     MergingTime - FittingTime - ConnectTime - FillingTime;
   timing("Time elapsed in z-reconstruction phase: %.6lf s.", zRecTime);
#endif
      

      
   //    if(RecoTracks.size() > 0){
   //for (size_t i =0; i < RecoTracks.size(); i++)
   //	free(TrackIDToMerge[i]);
   //free(TrackIDToMerge);
  

   std::vector < MCTrackObject* > *MCTracksEvt = MCTracks->at(evt);
   std::sort(MCTracksEvt->begin(), MCTracksEvt->end(), greaterThanLength);


   #if (WRITE_LIST_RECO_ID > 0)
   if (MCOutTxtFile.is_open()) {
     for(unsigned int l = 0, cm =0; l < MCTracksEvt->size(); l++){
       MCTrackObject const *CurMCtrack = MCTracksEvt->at(l);
       std::set<int> MCSttComp((CurMCtrack->m_STT_Component).begin(), (CurMCtrack->m_STT_Component).end());
       if(MCSttComp.size()>5){
	 for(auto it = MCSttComp.begin(); it != MCSttComp.end(); ++it){
	   int id = *it;
	   MCOutTxtFile << evt+firstEvt << "," << cm << "," << id << '\n';
	   MCDetID.Fill(evt+firstEvt,
			  cm,
			  id);	  
	 }
	 cm++;
       }
     }
   }
   #endif
   
   for(unsigned int l = 0; l < RecoTracks.size(); l++){
     PathCandidate *Cand          = RecoTracks[l];
     std::set<int> const *trk     = Cand->m_memberIdSet;
     std::vector<int> const *vect = Cand->m_memberList;

     if(Cand->m_isValid) {
       std::set<int> *comp = new std::set<int>((*trk));	    
       RecoTrackListDetID->push_back(comp);
     }
	
   }
    
   int NumConnComp = RecoTrackListDetID->size();
   info("Number of tracks reconstructed: %d", NumConnComp);
	      
   ComponentPerEvt.Fill(evt, NumConnComp);

   TimeNow      =std::chrono::high_resolution_clock::now() ;
   EvtTotalTime = std::chrono::duration<double>(TimeNow - TimeEvt).count();
   RecTotalTime += EvtTotalTime;
   timing("Time elapsed for event %d (%d hits, %d tracks): %.6lf s",
	  evt+firstEvt, totHits, nMCTracks, EvtTotalTime);

   std::vector<int> BestRecoMCMatchIDs = MatchBestRecoToMC( GridStruct, MCTracksEvt, &RecoTracks);

	
   for(unsigned int l = 0, cm = 0; l < RecoTracks.size(); l++){
     PathCandidate &curCand       = *(RecoTracks[l]);
     std::set<int>  trkset        = *curCand.m_memberIdSet;
     std::vector<int> const *vect = curCand.m_memberList;
     std::vector<double> &x       = curCand.m_x;
     std::vector<double> &y       = curCand.m_y;
     std::vector<double> &z       = curCand.m_z; 
     std::vector<double> &r       = curCand.m_r;
     std::vector<double> &theta   = curCand.m_theta; 

     if(curCand.m_isValid) {
       for( size_t i = 0;  i < vect->size(); ++i) {
	 int detID       = vect->at(i);
	 GridNode  &node = DetGrid[detID-1];
	 ConnectedCoord.Fill(evt,
			     cm,
			     BestRecoMCMatchIDs[l],
			     node.m_detID,
			     node.m_x,
			     node.m_y,
			     node.m_z,
			     node.m_r,
			     node.m_thetaDeg,
			     x[i],
			     y[i],
			     z[i]);

       }
       
       for( size_t i = 0;  i < curCand.m_anchors.size(); ++i) {
	 AnchorCCCoord.Fill(evt,
			    cm,
			    curCand.m_anchors[i].m_xDet,
			    curCand.m_anchors[i].m_yDet,
			    curCand.m_anchors[i].m_zDet);
       }      
       cm++;		  
     }
   }

   // Write for QA
   #if (WRITE_LIST_RECO_ID > 0)  
   if (RecOutTxtFile.is_open()) {
     for(unsigned int l = 0, cm = 0; l < RecoTracks.size(); l++){
       PathCandidate &curCand       = *(RecoTracks[l]);
       if(curCand.m_isValid) {
	 std::set<int>  trkset        = *curCand.m_memberIdSet;
	 auto it = trkset.upper_bound(GridStruct.firstVirtIdx);
	 trkset.erase(it, trkset.end());
	 if(trkset.size()>5){
	   for(auto it = trkset.begin(); it != trkset.end(); ++it){
	     int id = *it;
	     RecOutTxtFile << evt+firstEvt << "," << cm << "," << id << '\n';
	     RecoDetID.Fill(evt+firstEvt,
			    cm,
			    id);	  
	   }
	   cm++;
	 }
       }
     }
   }
   #endif // End IF writing to ascii file

   
   // Write to an ASCII file
   #if (WRITE_CM_ASCII_FILE > 0)  
   std::string header = "trkid,detid,x,y,z,cx,cy,cz\n";
   std::ofstream OutTxtFile;
   OutTxtFile.open ("ReconstructedTracksList.csv");
    if (OutTxtFile.is_open()) {
     OutTxtFile << header;
     for(unsigned int l = 0; l < RecoTracks.size(); l++){
       PathCandidate &curCand = *(RecoTracks[l]);
       if(curCand.m_isValid) {
	 std::vector<int>  *trk = curCand.m_memberList;
	 for(size_t cm = 0 ; cm < trk->size(); ++cm) {
	   int detID = trk->at(cm);
	   GridNode  &node = DetGrid[detID-1];
	   OutTxtFile << l << "," << detID<< "," << node.m_x << "," << node.m_y << "," << node.m_z << ","
		      << node.m_xDet<< "," << node.m_yDet << "," << node.m_zDet << '\n';

	 }
       }
     }
   }
   #endif // End IF writing to ascii file


   #if (WRITE_TIME_ASCII_FILE > 0)  
   if (TimeTxtFile.is_open()) {
     TimeTxtFile << evt+firstEvt << "," << nMCTracks << "," << GridStruct.m_STT_idx.size() << ","
		 << FillingTime << "," << ConnectTime << ","<< FittingTime << ","
		 << MergingTime << "," << zRecTime<< "," <<EvtTotalTime << '\n';	  
   }
   #endif


	
   #if(EVALUATE_ERROR)
	
   EvtErrorStruct *EvtError =  ComputeGlobalEvtErrors(GridStruct, MCTracksEvt, RecoTrackListDetID);
   if(EvtError != 0){
     GlobalErrorNtuple.Fill(static_cast<float>(evt+firstEvt),
			  static_cast<float>(EvtError->nMCTracks),
			  static_cast<float>(EvtError->nRecoTracks),
			  static_cast<float>(EvtError->UnderMerge),
			  static_cast<float>(EvtError->OverMerge),
			  static_cast<float>(EvtError->TotalError),
			  static_cast<float>(EvtError->UnderMergeNorm),
			  static_cast<float>(EvtError->OverMergeNorm),
			  static_cast<float>(EvtError->TotalErrorNorm));
	
     delete EvtError;
   } else {
     GlobalErrorNtuple.Fill(static_cast<float>(evt+firstEvt),
			    static_cast<float>(nMCTracks),
			    static_cast<float>(0),
			    static_cast<float>(1),
			    static_cast<float>(0),
			    static_cast<float>(1),
			    static_cast<float>(1),
			    static_cast<float>(0),
			    static_cast<float>(1));
   }

   std::vector< TrackErrorStruct* > *PandaErrors = PandaErrorMetric(GridStruct, MCTracksEvt, &RecoTracks);

   if(PandaErrors != 0) {
     for(size_t f = 0; f < PandaErrors->size(); ++f) {
       TrackErrorStruct const *erPandaObj = PandaErrors->at(f);	    
       ErrorPerRecoTrackNtuple.Fill(static_cast<float>(evt+firstEvt),
				    static_cast<float>(erPandaObj->tr_rank),
				    static_cast<float>(erPandaObj->tr_isClone));
     }
   

     for(size_t r = 0; r < PandaErrors->size(); ++r) {
       delete PandaErrors->at(r);
     }
   }
   delete PandaErrors;
	
   std::vector< TrackErrorStruct* > *TrackErrors = ComputeErrorPerRecoTrack(GridStruct, MCTracksEvt,
									    &RecoTracks, ListComplexTracks,
									    BestRecoMCMatchIDs);

   if(TrackErrors != 0) {
     for(size_t f = 0; f < TrackErrors->size(); ++f) {
       TrackErrorStruct const *erObj      = TrackErrors->at(f);

       ErrorPerMCTrackNtuple.Fill(static_cast<float>(evt+firstEvt),
				  static_cast<float>(erObj->isComplex),
				  static_cast<float>(erObj->isNotmatched),
				  static_cast<float>(erObj->MatchIndex),
				  static_cast<float>(erObj->MCTrackLength),
				  static_cast<float>(erObj->RecoTrackLength),
				  static_cast<float>(erObj->IntersectionLength),
				  static_cast<float>(erObj->UnionLength),
				  static_cast<float>(erObj->F1score),
				  static_cast<float>(erObj->MeanDiffX),
				  static_cast<float>(erObj->MeanDiffY),
				  static_cast<float>(erObj->MeanDiffZ));
	    
       CurvPerTrackNtuple.Fill(static_cast<float>(evt+firstEvt),
			       static_cast<float>(erObj->MC_px),
			       static_cast<float>(erObj->MC_py),
			       static_cast<float>(erObj->MC_pz),
			       static_cast<float>(erObj->MC_r1),
			       static_cast<float>(erObj->MC_r2),
			       static_cast<float>(erObj->tr_rIsoRand),
			       static_cast<float>(erObj->tr_rIsoRand16),
			       static_cast<float>(erObj->tr_rIsoRand84),
			       static_cast<float>(erObj->tr_rAnc),
			       static_cast<float>(erObj->tr_rPts),
			       static_cast<float>(erObj->tr_scattAngle));
	    
       if(!erObj->isNotmatched)
	 for(size_t p = 0; p < erObj->DiffX.size(); p++){
	   CoordDiffPerTrackNtuple.Fill(erObj->DiffX[p],
					erObj->DiffY[p],
					erObj->DiffZ[p]);
	 } // End For Diff Coord errors	    
	    
     } // End For TrackErrors.size()
	  
     for(size_t r = 0; r < TrackErrors->size(); ++r) {
       delete TrackErrors->at(r);
     }	  
     delete TrackErrors;
	  
   } // End If Writing errors 
	

   #endif // End IF error computation


   // Cleaning

   if(RecoTracks.size()>0){
     for(size_t c = 0; c < RecoTracks.size(); ++c) {
       delete(RecoTracks[c]);
     } 
     RecoTracks.clear();
   }
   
   if(RecoTrackListDetID != 0) {
     for(size_t c = 0; c < RecoTrackListDetID->size(); ++c) {
       delete RecoTrackListDetID->at(c);
     }
     delete RecoTrackListDetID;
   }

   CollectGridToTree(GridStruct, GridCoordNtuple);     
   GridStruct.ResetGrid();
   free(visited);
     
   if(nEvtProc == 15000)

     //if(nEvtProc == 4950) // 4 events batch
     break;
  }
  #endif
  
  timer.Stop();
  double timeTotal = std::chrono::duration<double>( std::chrono::high_resolution_clock::now() - TimeAll ).count();

  ComponentPerEvt.Write();
  ConnectedCoord.Write();
  AnchorCCCoord.Write();
  GridCoordNtuple.Write();

  #if(EVALUATE_ERROR)

  //Write error estimations.
  GlobalErrorNtuple.Write();
  ErrorPerMCTrackNtuple.Write();
  ErrorPerRecoTrackNtuple.Write();
  CoordDiffPerTrackNtuple.Write();
  CurvPerTrackNtuple.Write();
  #endif
  
  OutputFile.Close();
  TFile DetecIDFile("LOTFRecoTrackIDs.root", "RECREATE",
   		     "List of tube IDS per Reco and MC track, file created by floodingFilter macro", 9);
  RecoDetID.Write();
  MCDetID.Write();
  DetecIDFile.Close();
  
  for(size_t i = 0; i < HitsData->size(); ++i) {
    std::vector<HitCoordinate*>* dd = (*HitsData)[i];
    for(size_t j = 0; j < dd->size(); ++j) {
      delete dd->at(j);
    }
    delete (*HitsData)[i];
  }
  HitsData->clear();
 
  for(size_t j = 0; j < MCTracks->size(); ++j) {
    std::vector < MCTrackObject* >* MCtracksDel = MCTracks->at(j);
    for(size_t k = 0; k < MCtracksDel->size(); ++k) {
      delete MCtracksDel->at(k);
    }
    delete (MCTracks->at(j));
  }
  delete MCTracks;

  
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();

  
  printf("Macro finished successfully. \n");
  printf("/************ LOTF Processing SUMMARY ****************/\n");
printf("Number of events processed is %d (last evt id %d)\n", nEvtProc, lastevtID);
printf("Number of tracks processed is %d (with > 5 hits %d)\n", nTrkProc, nTrk5hitsProc);
  printf("Number of hits processed is %d\n", nHitsProc);
  printf("Reconstruction time %.5lf s, (%.0lf hits/s, %.5f s/Event)\n",
	 RecTotalTime, nHitsProc/RecTotalTime, RecTotalTime/nEvtProc);
  printf("Real time including other processes: %.5f s (%.5f s/Event).\n", (rtime), (rtime/nEvtProc));
  printf("CPU time %.5f s (%.5f s/Event).\n", (ctime), (ctime/nEvtProc));
  printf("/************ LOTF Processing ENDED ****************/\n");

}
  // Extract MC track values
  

