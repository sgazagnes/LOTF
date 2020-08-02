
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
#include "auxiliaryfunctions.h"
#include "CollectSttMvdPoints.h"
#include "SttMVDEventDataReader.h"
#include "pathopen.h"
#include "hitcoordinate.h"
#include "floodingFilter.h"
#include "utilfunctions.h"
#include "trackObject.h"
#include "logc.h"
#include "queue.h"
//#include "performFilter.h"
#include "pathCandidate.h"


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

typedef enum {
  DOWN = 1,
  SAME = 2,
  UP = 4
} Direction;

bool sortNeighbors(CoordGrid *gr, GridNode *currentNode, std::vector<int> *prev, std::vector<int> *same, std::vector<int> *next, std::vector<int> *virt, char *visited, int *dir){

  int curDir = *dir;
  std::vector< GridNode > &Ingrid = gr->m_grid;
  int curLayer = currentNode->m_Layer;
  int curId =  currentNode->m_detID;
  bool cond = true;

  
  for(int i = 0; i < currentNode->m_neighbors.size(); i++){
    int neighId = currentNode->m_neighbors[i];
    int neighIdx = gr->Find(neighId);
    GridNode *neighNode = &Ingrid[neighIdx];
    
    if(neighNode->m_type == GridNode::VIRTUAL_NODE){
      virt->push_back(neighId);
      debug("%d is a virtual node, add to list and find next neighbor", neighId);
      neighId    = neighNode->m_neighbors[0] == curId ?  neighNode->m_neighbors[1]: neighNode->m_neighbors[0];
      neighIdx   = gr->Find(neighId);
      neighNode  = &Ingrid[neighIdx];
    }	 

    if(visited[neighId] < 1){
      if(neighNode->m_Layer > curLayer){
	debug("Node %d has one neigh up %d", curId, neighId);
	next->push_back(neighId);
	curDir |= UP;
	visited[neighId] = 2;
      } else if( neighNode->m_Layer < curLayer) {
	debug("Node %d has one neigh down %d", curId, neighId);
	prev->push_back(neighId);
	curDir |= DOWN;
	visited[neighId] = 2;
      } else {
	debug("Node %d has one neigh on the same %d", curId, neighId);
	same->push_back(neighId);
	//	curDir |= SAME;
	visited[neighId] = 2;
      }
    } else if(visited[neighId] == 4){
      debug("Node %d has already been connected, tricky", neighId);
      cond = false;
    } else
      debug("Node %d has already been added, tricky", neighId);
  }

  	
  if( curDir > 4 || same->size() > 1){
    info("Too many possibilities, let skip it for now");
    cond = false;
  }
  *dir = curDir;

  return cond;
  
}

void resetLists(char *visited, std::vector<int> *prev, std::vector<int> *same, std::vector<int> *next){
  for (int i = 0; i < same->size(); i++)
    visited[same->at(i)] =0;
	      
  for (int i = 0; i < next->size(); i++)
    visited[next->at(i)] =0;

  for (int i = 0; i < prev->size(); i++)
    visited[prev->at(i)] =0;
}

void removeIdFromNeigh(GridNode *neighNode, std::vector<int> *prevNodes, int curId){
  for(int i = 0; i < prevNodes->size(); i++){
    if(prevNodes->at(i) != neighNode->m_detID){ 		
      (neighNode->m_neighbors).erase(std::remove((neighNode->m_neighbors).begin(), (neighNode->m_neighbors).end(),prevNodes->at(i)), (neighNode->m_neighbors).end());
      if(i == 0)
	neighNode->parent = curId;
    }
  }
}

bool areAdjacent(CoordGrid *gr, std::vector<int> *v){
  int adjacent = 0;
  std::vector< GridNode > &Ingrid = gr->m_grid;

  for (int i = 0; i < v->size(); i++){
    int neighId   = v->at(i);
    int neighIdx  = gr->Find(neighId);
    GridNode *neighNode = &Ingrid[neighIdx];
    //prevNodes.push_back(neighId);
    
    for (int j = i+1; j < v->size (); j++){
      debug("Are %d and %d connected?", neighId, v->at(j));
      if(neighId == v->at(j)) error("HOUSTON");
      else if(neighNode->IsNeighboring(v->at(j))){
	debug("Yes");
	adjacent++;
      } else
	debug("No");
    }
  }
  
  if(adjacent >= v->size() -1)
    return true;
  else
    return false;
}

//______________________ BEGIN MCTrackPoints _____________________________
std::vector< std::vector < MCTrackObject* >* >*
MCTrackPoints( std::vector < std::vector<HitCoordinate*>* > const &evtData)
{
  info("Extracting MC tracks for %d events", evtData.size());
  //	    << " events.\n";
  // Output Parameter
  std::vector< std::vector < MCTrackObject* >* >* outVar =
    new std::vector< std::vector < MCTrackObject* >* >();

  int numTracks = -1;
   for(size_t e = 0; e < evtData.size(); ++e) {
    std::vector<HitCoordinate*> const *Current_Event = evtData[e];
    // numTracks = 0;
    // Find out how many MC tracks are available
    for(size_t i = 0; i < Current_Event->size(); ++i) {
      HitCoordinate const *currentHit = Current_Event->at(i);
      if( currentHit->m_trackID > numTracks &&
	  currentHit->m_trackID != HIT_EXCLUSION) {
	numTracks = currentHit->m_trackID;
      }
    }// END current event loop
    
    debug("Event %d contains %d tracks", e,(numTracks + 1));
    
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

    std::vector< GridNode > &Ingrid = gr.m_grid;  
    std::vector< int > activeId;
    std::vector< int > remainingActiveId;


    /* Pushing all active detectors into queue */
    
    int nactive = 0, nactive_queue = 0;
    for(unsigned int n = 0; n < Ingrid.size(); ++n) {
      if( Ingrid[n].m_active){
	int NodeId = Ingrid[n].m_detID;
	if(Ingrid[n].m_type != GridNode::VIRTUAL_NODE  ) {
	  activeId.push_back(NodeId);
	  nactive_queue++;
	}
	nactive++;
      }
      
    }

    info("Found the active detectors");

    /* Keep in neighbors only active ones */
    
    for(unsigned int n = 0; n < nactive_queue; ++n) {
      int curid    = activeId[n];
      int curindex = gr.Find(curid);
      GridNode &current_Node = Ingrid[curindex];
      int n_neighbors = current_Node.m_neighbors.size();
      //  debug("%d, has %d neighbors", current_Node.m_detID, n_neighbors);
      for  ( int i = 0; i < current_Node.m_neighbors.size(); i++){
	int neigh_ID = current_Node.m_neighbors[i];
	//	debug("%d curr", neigh_ID);
	int neigh_index = gr.Find(neigh_ID);
	if(!Ingrid[neigh_index].m_active){
	  //  debug("Removing");
	  (current_Node.m_neighbors).erase((current_Node.m_neighbors).begin()+i);
	  i--;
	}
      }
      n_neighbors = current_Node.m_neighbors.size();
      //    debug("%d, has %d neighbors", current_Node.m_detID, n_neighbors);

    }

    info("Removed non active neighbors");

    
    int candidateId = 0;
    std::vector < PathCandidate* > tracklets;
    char *visited = (char *)calloc(50000, sizeof(int));

    info("First step, let's find the obvious tracks");

    
    for(unsigned int n = 0; n < nactive_queue; ++n) {
      
      std::vector<int> sameLayer;
      std::vector<int> nextLayer;
      std::vector<int> prevLayer;
      std::vector<int> nextVirt;
      std::vector<int> prevNodes;
      std::vector<int> *v;

      int curId = activeId[n];
      int curIdx = gr.Find(curId);
      GridNode *currentNode = &Ingrid[curIdx];
      
      int dir = 0;

      if(visited[curId] == 0 && (currentNode->m_LayerLimit == 1 ||  currentNode->m_neighbors.size() == 1)){	
	int n_neighbors = currentNode->m_neighbors.size();
	int curLayer = currentNode->m_Layer;

	info("\n Starting node %d, neigh %d", curId, n_neighbors);

	int neighId;
	int neighIdx;
	GridNode *neighNode;
	bool cond = true;

	cond = sortNeighbors(&gr, currentNode, &prevLayer, &sameLayer, &nextLayer, &nextVirt, visited,  &dir);
	
	if(cond == false){
	  info("Too many possibilities, let skip it for now");
	  resetLists(visited, &prevLayer, &sameLayer, &nextLayer);
	  continue;
	}
	
	PathCandidate *cand = new PathCandidate();// Create a new candidate
	cand->m_id = candidateId++;// Set id
	cand->m_isValid = true;
	cand->insertNewNode(currentNode);
	cand->m_firstNodeVisited = curId;
	if(abs(currentNode->m_SectorLimit) > 0)
	  cand->m_isOnSectorLimit = true;
	visited[curId] = 4;
	prevNodes.push_back(curId);

	int n_connected = 0;

	while(cond){


	  /* 1 NEIGHBOR */
	  
	  if(n_neighbors == 1){
	    info("dir %d", dir);
	    if (dir == UP){
	      v = &nextLayer;
	      info("Going up from node %d to node %d", curId, v->at(0));
	    } else if (dir == DOWN){
	      v = &prevLayer;
	      info("Going down from node %d to node %d", curId, v->at(0));
	    } else {
	      v = &sameLayer;
	      info("Going same level from node %d to node %d", curId, v->at(0));
	    }
	     
	    
	    if(nextVirt.size() > 0){
	      debug("Handling the virtuals");
	      
	      for (int i = 0; i <nextVirt.size(); i++){
		neighId   = nextVirt[i];
		neighIdx  = gr.Find(neighId);
		neighNode = &Ingrid[neighIdx];
		cand->insertNewNode(neighNode);
		visited[neighId] = 4;
		n_connected++;

		removeIdFromNeigh(neighNode, &prevNodes, curId);
	      }
	      
	      prevNodes.clear();
	      prevNodes.insert(prevNodes.end(),  nextVirt.begin(),  nextVirt.end());
	      nextVirt.clear();
	    }
	      
	    

	    neighId    = v->at(0);
	    neighIdx   = gr.Find(neighId);
	    neighNode  = &Ingrid[neighIdx];
	    cand->insertNewNode(neighNode);
	    visited[neighId] = 4;
	    n_connected++;
	    
	    removeIdFromNeigh(neighNode, &prevNodes, curId);

	    curId       = neighId;
	    curIdx      = neighIdx;
	    currentNode = neighNode;
	    curLayer    = currentNode->m_Layer;
	    
	    nextLayer.clear();  sameLayer.clear();   prevLayer.clear();
	    dir = 0;
	    prevNodes.clear();
	    prevNodes.push_back(curId);
	    
	    cond = sortNeighbors(&gr, currentNode, &prevLayer, &sameLayer, &nextLayer, &nextVirt, visited,  &dir);

	    n_neighbors = sameLayer.size() + prevLayer.size() + nextLayer.size();
	    info("%d nodes were connected,  %d found for the next step \n", n_connected, n_neighbors);
	    n_connected = 0;

	    if(cond == false){
	      info("Too many possibilities, let skip it for now");
	      resetLists(visited, &prevLayer, &sameLayer, &nextLayer);
	      continue;
	    }
	    
	  }


	  /* 1 Same layer neighbor */

	  else if (sameLayer.size() > 0){

	    v = dir == UP ? &nextLayer: &prevLayer;

	    int candId = sameLayer[0];
	    int candIdx        = gr.Find(candId);
	    GridNode *candNode = &Ingrid[candIdx];
	    curLayer           = candNode->m_Layer;
	    sameLayer.clear();

	    info("Still on the same layer, investigating node %d", candId);


	    bool toadd = true;
	    int cursize = v->size();
	    
	    for(int i = 0; i < candNode->m_neighbors.size(); i++){
	      neighId = candNode->m_neighbors[i];
	      if (neighId == curId) continue;
	      
	      neighIdx = gr.Find(neighId);
	      neighNode = &Ingrid[neighIdx];
	     
	      if(neighNode->m_type == GridNode::VIRTUAL_NODE){
		nextVirt.push_back(neighId);
		neighId = neighNode->m_neighbors[0] == candId ? neighNode->m_neighbors[1]: neighNode->m_neighbors[0];
		neighIdx = gr.Find(neighId);
		neighNode  = &Ingrid[neighIdx];
	      }
	      
	      for (int j  = 0; j < cursize; j++){
		debug("Connection between %d and %d ?", neighId, v->at(j));

		if(v->at(j) == neighId)
		  debug("They are the same");
		else if(neighNode->IsNeighboring(v->at(j))){
		  debug("They are neighbors",  neighId, v->at(j));
		  if(visited[neighId] < 1){
		    debug("Not visited yet");
		    if(neighNode->m_Layer > curLayer){
		      debug("Pushing %d on the list nextLayer", neighId);
		      nextLayer.push_back(neighId);
		      dir |= UP;
		      visited[neighId] = 2;
		    }  else if( neighNode->m_Layer < curLayer){
		      debug("Pushing %d on the list prevLayer", neighId);
		      prevLayer.push_back(neighId);
		      dir |= DOWN;
		      visited[neighId] = 2;
		    } else{
		      debug("Pushing %d on list sameLayer", neighId);
		      sameLayer.push_back(neighId);
		      visited[neighId] = 2;
		    }
		  } else
		    debug("%d has already been visited",neighId);
		} else {
		  debug("%d and %d are not (obviously) connected", v->at(j), neighId);
		  toadd = false;
		  cond = false;
		  break;
		}
	      }

	      if (cond == false) break;
	    } // for in neighbors

	    if (toadd == true){
	      cand->insertNewNode(candNode);
	      visited[candId] = 4;
	      n_connected++;

	      removeIdFromNeigh(candNode, &prevNodes, curId);
	      
	      prevNodes.push_back(candId);
	      curId       = candId;
	      curIdx      = candIdx;
	      currentNode = candNode;
	      n_neighbors = sameLayer.size()+prevLayer.size()+nextLayer.size();
	      info("%d nodes were connected,  %d found for the next step \n\n", n_connected, n_neighbors);
	      n_connected = 0;
	    }
	    else {

	      info("To many possibilities, we'll see what's next \n\n", n_connected, n_neighbors);

	      visited[candId] = 0;

	      resetLists(visited, &prevLayer, &sameLayer, &nextLayer);

	    }

	  }



	  /* MORE NEIGHBOOOORS */




	  else if (n_neighbors > 1) {

	    info("%d", dir);
	    if (dir == UP){	      
	      info("Next nodes are up and we have %d of them", nextLayer.size());
	      v = &nextLayer;
	    }  else if (dir == DOWN){
	      info("NExt nodes are down and we have %d of them", prevLayer.size());
	      v = &prevLayer;
	    } else
	      info("WHAT IS THE DIRECTION NOW?");
	    
	    
	    if(areAdjacent(&gr, v)){
	      info("Adding %d nodes to the CM", v->size());

	      if(nextVirt.size() > 0){
		debug("Handling the virtual");
		for (int i = 0; i < nextVirt.size(); i++){
		  neighId   = nextVirt[i];
		  neighIdx  = gr.Find(neighId);
		  neighNode = &Ingrid[neighIdx];
		  cand->insertNewNode(neighNode);
		  visited[neighId] = 4;
		  n_connected++;

		  removeIdFromNeigh(neighNode, &prevNodes, curId);		  
		}
		
		prevNodes.clear();
		prevNodes.insert(prevNodes.end(),  nextVirt.begin(),  nextVirt.end());
		nextVirt.clear();

	      }
	      
	      std::vector<int> lookneigh;

	      for (int i = 0; i < v->size(); i++){
		neighId = v->at(i);
		neighIdx = gr.Find(neighId);
		neighNode = &Ingrid[neighIdx];
		cand->insertNewNode(neighNode);
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
	      info("New current node is %d, looking for neighbors of %d nodes", curId, lookneigh.size());
	      
	      n_neighbors = 0;
	      
	      for(int i = 0; i < lookneigh.size(); i++){
		int id         = lookneigh[i];
		int idx        = gr.Find(id);
		GridNode *node = &Ingrid[idx];
		prevNodes.push_back(id);

		cond = sortNeighbors(&gr, node, &prevLayer, &sameLayer, &nextLayer, &nextVirt, visited,  &dir);
	      }

	      n_neighbors = sameLayer.size()+prevLayer.size()+nextLayer.size();
	      info("%d nodes were connected,  %d found for the next step \n", n_connected, n_neighbors);
	      n_connected = 0;
	     

	    } else {
	      info("Some of these nodes are node adjacent", v->size());
	      cond = false;
	    }
	    
	      
	  } else {

	    if(n_neighbors == 0) {
	      info("No more neighbors in sight, are we finished here ? \n\n");
	       if(cand->m_minLayer == 0 && cand->m_maxLayer > 21){
		 info("track goes through all layers, likily finished");
		 cand->m_finished = 1;
	       } else if(abs(currentNode->m_SectorLimit) > 0 || cand->m_isOnSectorLimit){
		 
		 cand->m_isOnSectorLimit= true;
		 info("Track is on sector limit, might have a connection somewhere else");
	       } else {
		 int firstId = cand->m_firstNodeVisited;
		 int firstIdx = gr.Find(firstId);
		 GridNode *firstNode = &Ingrid[firstIdx];
		 if(firstNode->m_LayerLimit == 1 && currentNode->m_LayerLimit == 1) {
		   info("starting and ending on same layer");
		   cand->m_finished = 1;

		 }
	       }
	    } else {
	      info("This track has other neighbors, to correct later");
	    }
	    cond = false;
	    
	  }


	}
	
       	info("Pushing cm %d: \n               length is %d, \n               Min layer %d, \n              Max layer %d, \n               IsOnSectorLimit %d, \n               Finished ? %d. ", cand->m_id, cand->m_length,  cand->m_minLayer, cand->m_maxLayer, cand->m_isOnSectorLimit, cand->m_finished);

	tracklets.push_back(cand);

      }
      // elser
      //	remainingActiveId.push_back(curId);
       	     
    }





    

 
    info(" Find remaining obvious tracklets");

    for(unsigned int n = 0; n < activeId.size(); ++n) {
      int curId = activeId[n];
      int oriId = curId;
      int prevId = curId;
      int curIdx = gr.Find(curId);
      GridNode *currentNode = &Ingrid[curIdx];
      GridNode *firstNode = &Ingrid[curIdx];
      std::vector<int> *nextIds;// List of members in a vector(delete me)

      int n_neighbors = currentNode->m_neighbors.size();

      if(n_neighbors == 2 && visited[curId] == 0) {

	int neighId1 = currentNode->m_neighbors[0];
	int neighIdx1 = gr.Find(neighId1);
	GridNode *neighNode1 = &Ingrid[neighIdx1];

	int neighId2 = currentNode->m_neighbors[1];
	int neighIdx2 = gr.Find(neighId2);
	GridNode *neighNode2 = &Ingrid[neighIdx2];
	
	if(!neighNode1->IsNeighboring(neighId2)){
	  
	  info("%d has two neighbors not connected %d %d, sounds like a good tracklet", curId, neighId1, neighId2);

	
	  PathCandidate *cand = new PathCandidate();// Create a new candIdate
	  cand->m_id = candidateId++;// Set id
	  cand->m_isValid = true;
	  cand->insertNewNode(currentNode);
	  visited[curId] = 4;

	  for (int k = 0; k < 2; k++){
	    std::vector<int> sameLayer;
	    std::vector<int> nextLayer;
	    std::vector<int> prevLayer;
	    std::vector<int> nextVirt;
	    std::vector<int> prevNodes;
	    std::vector<int> *v;

	    int dir = 0;
      
	    int curLayer = firstNode->m_Layer;
	    
	    int neighId =  firstNode->m_neighbors[k];
	    int neighIdx = gr.Find(neighId);
	    GridNode *neighNode = &Ingrid[neighIdx];
	    
	    info("Handling first node %d", neighId);

	    cand->insertNewNode(neighNode);
	    visited[neighId] = 4;
	    
	    (neighNode->m_neighbors).erase(std::remove((neighNode->m_neighbors).begin(), (neighNode->m_neighbors).end(), oriId), (neighNode->m_neighbors).end());

	    //  prevId = curId;
	    curId = neighId;
	    curIdx = neighIdx;
	    currentNode = neighNode;
	    curLayer = currentNode->m_Layer;
	    
	    bool cond = true;
	    
	    cond = sortNeighbors(&gr, currentNode, &prevLayer, &sameLayer, &nextLayer, &nextVirt, visited,  &dir);

	    
	    if(cond == false){
	      if(k == 0) cand->m_firstNodeVisited = curId;
	      else cand->m_lastNodeVisited = curId;	      
	      continue;
	    }

	  

	    prevNodes.push_back(curId);
	    n_neighbors = prevLayer.size() + sameLayer.size() + nextLayer.size();
	    int n_connected = 0;

	    info("Next number of neighbors is %d", n_neighbors);
	    
	    while(cond){


	      /* 1 NEIGHBOR */
	  
	      if(n_neighbors == 1){
		info("dir %d", dir);
		if (dir == UP){
		  v = &nextLayer;
		  info("Going up from node %d to node %d", curId, v->at(0));
		} else if (dir == DOWN){
		  v = &prevLayer;
		  info("Going down from node %d to node %d", curId, v->at(0));
		} else {
		  v = &sameLayer;
		  info("Going same level from node %d to node %d", curId, v->at(0));
		}
	     
	    
		if(nextVirt.size() > 0){
		  debug("Handling the virtuals");
	      
		  for (int i = 0; i <nextVirt.size(); i++){
		    neighId   = nextVirt[i];
		    neighIdx  = gr.Find(neighId);
		    neighNode = &Ingrid[neighIdx];
		    cand->insertNewNode(neighNode);
		    visited[neighId] = 4;
		    n_connected++;

		    removeIdFromNeigh(neighNode, &prevNodes, curId);
		  }
	      
		  prevNodes.clear();
		  prevNodes.insert(prevNodes.end(),  nextVirt.begin(),  nextVirt.end());
		  nextVirt.clear();
		}
	      
	    

		neighId    = v->at(0);
		neighIdx   = gr.Find(neighId);
		neighNode  = &Ingrid[neighIdx];
		cand->insertNewNode(neighNode);
		visited[neighId] = 4;
		n_connected++;
	    
		removeIdFromNeigh(neighNode, &prevNodes, curId);

		curId       = neighId;
		curIdx      = neighIdx;
		currentNode = neighNode;
		curLayer    = currentNode->m_Layer;
	    
		nextLayer.clear();  sameLayer.clear();   prevLayer.clear();
		dir = 0;
		prevNodes.clear();
		prevNodes.push_back(curId);
	    
		cond = sortNeighbors(&gr, currentNode, &prevLayer, &sameLayer, &nextLayer, &nextVirt, visited,  &dir);

		n_neighbors = sameLayer.size() + prevLayer.size() + nextLayer.size();
		info("%d nodes were connected,  %d found for the next step \n", n_connected, n_neighbors);
		n_connected = 0;

		if(cond == false){
		  info("Too many possibilities, let skip it for now");
		  resetLists(visited, &prevLayer, &sameLayer, &nextLayer);
		  if(k == 0) cand->m_firstNodeVisited = curId;
		  else cand->m_lastNodeVisited = curId;	      
		  
		  
		  continue;
		}
	    
	      }


	      /* 1 Same layer neighbor */

	      else if (sameLayer.size() > 0){

		v = dir == UP ? &nextLayer: &prevLayer;

		int candId = sameLayer[0];
		int candIdx        = gr.Find(candId);
		GridNode *candNode = &Ingrid[candIdx];
		curLayer           = candNode->m_Layer;
		sameLayer.clear();

		info("Still on the same layer, investigating node %d", candId);


		bool toadd = true;
		int cursize = v->size();
	    
		for(int i = 0; i < candNode->m_neighbors.size(); i++){
		  neighId = candNode->m_neighbors[i];
		  if (neighId == curId) continue;
	      
		  neighIdx = gr.Find(neighId);
		  neighNode = &Ingrid[neighIdx];
	     
		  if(neighNode->m_type == GridNode::VIRTUAL_NODE){
		    nextVirt.push_back(neighId);
		    neighId = neighNode->m_neighbors[0] == candId ? neighNode->m_neighbors[1]: neighNode->m_neighbors[0];
		    neighIdx = gr.Find(neighId);
		    neighNode  = &Ingrid[neighIdx];
		  }
	      
		  for (int j  = 0; j < cursize; j++){
		    debug("Connection between %d and %d ?", neighId, v->at(j));

		    if(v->at(j) == neighId)
		      debug("They are the same");
		    else if(neighNode->IsNeighboring(v->at(j))){
		      debug("They are neighbors",  neighId, v->at(j));
		      if(visited[neighId] < 1){
			debug("Not visited yet");
			if(neighNode->m_Layer > curLayer){
			  debug("Pushing %d on the list nextLayer", neighId);
			  nextLayer.push_back(neighId);
			  dir |= UP;
			  visited[neighId] = 2;
			}  else if( neighNode->m_Layer < curLayer){
			  debug("Pushing %d on the list prevLayer", neighId);
			  prevLayer.push_back(neighId);
			  dir |= DOWN;
			  visited[neighId] = 2;
			} else{
			  debug("Pushing %d on list sameLayer", neighId);
			  sameLayer.push_back(neighId);
			  visited[neighId] = 2;
			}
		      } else
			debug("%d has already been visited",neighId);
		    } else {
		      debug("%d and %d are not (obviously) connected", v->at(j), neighId);
		      toadd = false;
		      cond = false;
		      break;
		    }
		  }

		  if (cond == false) break;
		} // for in neighbors

		if (toadd == true){
		  cand->insertNewNode(candNode);
		  visited[candId] = 4;
		  n_connected++;

		  removeIdFromNeigh(candNode, &prevNodes, curId);
	      
		  prevNodes.push_back(candId);
		  curId       = candId;
		  curIdx      = candIdx;
		  currentNode = candNode;
		  n_neighbors = sameLayer.size()+prevLayer.size()+nextLayer.size();
		  info("%d nodes were connected,  %d found for the next step \n\n", n_connected, n_neighbors);
		  n_connected = 0;
		}
		else {

		  info("To many possibilities, we'll see what's next \n\n", n_connected, n_neighbors);

		  visited[candId] = 0;

		  resetLists(visited, &prevLayer, &sameLayer, &nextLayer);
		  if(k == 0) cand->m_firstNodeVisited = curId;
		  else cand->m_lastNodeVisited = curId;	      
		  

		}

	      }



	      /* MORE NEIGHBOOOORS */




	      else if (n_neighbors > 1) {

		info("%d", dir);
		if (dir == UP){	      
		  info("Next nodes are up and we have %d of them", nextLayer.size());
		  v = &nextLayer;
		}  else if (dir == DOWN){
		  info("NExt nodes are down and we have %d of them", prevLayer.size());
		  v = &prevLayer;
		} else
		  info("WHAT IS THE DIRECTION NOW?");
	    
	    
		if(areAdjacent(&gr, v)){
		  info("Adding %d nodes to the CM", v->size());

		  if(nextVirt.size() > 0){
		    debug("Handling the virtual");
		    for (int i = 0; i < nextVirt.size(); i++){
		      neighId   = nextVirt[i];
		      neighIdx  = gr.Find(neighId);
		      neighNode = &Ingrid[neighIdx];
		      cand->insertNewNode(neighNode);
		      visited[neighId] = 4;
		      n_connected++;

		      removeIdFromNeigh(neighNode, &prevNodes, curId);		  
		    }
		
		    prevNodes.clear();
		    prevNodes.insert(prevNodes.end(),  nextVirt.begin(),  nextVirt.end());
		    nextVirt.clear();

		  }
	      
		  std::vector<int> lookneigh;

		  for (int i = 0; i < v->size(); i++){
		    neighId = v->at(i);
		    neighIdx = gr.Find(neighId);
		    neighNode = &Ingrid[neighIdx];
		    cand->insertNewNode(neighNode);
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
		  info("New current node is %d, looking for neighbors of %d nodes", curId, lookneigh.size());
	      
		  n_neighbors = 0;
	      
		  for(int i = 0; i < lookneigh.size(); i++){
		    int id         = lookneigh[i];
		    int idx        = gr.Find(id);
		    GridNode *node = &Ingrid[idx];
		    prevNodes.push_back(id);

		    cond = sortNeighbors(&gr, node, &prevLayer, &sameLayer, &nextLayer, &nextVirt, visited,  &dir);
		  }

		  n_neighbors = sameLayer.size()+prevLayer.size()+nextLayer.size();
		  info("%d nodes were connected,  %d found for the next step \n", n_connected, n_neighbors);
		  n_connected = 0;
	     

		} else {
		  info("Some of these nodes are node adjacent", v->size());
		  cond = false;
		}

		if(cond == false){
		  if(k == 0) cand->m_firstNodeVisited = curId;
		  else cand->m_lastNodeVisited = curId;	      
		}
	      
	      } else {

		if(n_neighbors == 0) {
		  info("No more neighbors in sight, are we finished here ? \n\n");
		  if(cand->m_minLayer == 0 && cand->m_maxLayer > 21){
		    info("track goes through all layers, likily finished");
		    cand->m_finished = 1;
		  } else if(abs(currentNode->m_SectorLimit) > 0 || cand->m_isOnSectorLimit){
		 
		    cand->m_isOnSectorLimit= true;
		    info("Track is on sector limit, might have a connection somewhere else");
		  } else {
		    int firstId = cand->m_firstNodeVisited;
		    int firstIdx = gr.Find(firstId);
		    GridNode *firstNode = &Ingrid[firstIdx];
		    if(firstNode->m_LayerLimit == 1 && currentNode->m_LayerLimit == 1) {
		      info("starting and ending on same layer");
		      cand->m_finished = 1;

		    }
		  }
		} else {
		  info("This track has other neighbors, to correct later");
		}
		cond = false;
	    
	      }

	    }
	

	  }
       	     
	

	  tracklets.push_back(cand);
	  info("Pushing cm %d: \n               length is %d, \n               Min layer %d, \n              Max layer %d, \n               IsOnSectorLimit %d, \n               Finished ? %d. ", cand->m_id, cand->m_length,  cand->m_minLayer, cand->m_maxLayer, cand->m_isOnSectorLimit, cand->m_finished);

	}

      }
    }

  
    
    info("Number of connected components: %d", tracklets.size());    
  
    for(unsigned int l = 0; l < tracklets.size(); l++){
      PathCandidate &curCand = *(tracklets[l]);
      std::set<int> const *trk = curCand.m_memberIdSet;
      //  if(curCand.m_isValid) {
      std::set<int> *comp = new std::set<int>((*trk));
	    
      connectedComp->push_back(comp);

      //}
    }
    
    int NumConnComp = connectedComp->size();
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
  // Extract MC track values
  
}
