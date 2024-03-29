/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#include <iostream>
#include <set>
#include <sstream>
#include <algorithm>
#include  <cmath>
#include <fstream>

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
#include "performFilter.h"
#include "utilfunctions.h"
#include "trackObject.h"

// DEBUG AND STORE definitions

#define EXCLUDE_STTSKEWED_PLOT   0
#define EXCLUDE_VIRTUALS_IN_PLOT 0

#define WRITE_CONNECTED_COMPONENTS 1
#define EVALUATE_ERROR 0
#define PRINT_DEBUG_INFO 0
#define INCLUDE_MVD_INOUTPUT_TRACK 1
#define PRINT_TRACKS_DEBUG_INFO 1
#define PRINT_DEBUG_INFO_COMP_MATCH 0
#define PRINT_ERROR_EXTRA_DEBUG_INFO 0
#define WRITE_CONNECTED_COMPONENTS_JSON 1
///_______________________________________________________________////
//______________________ BEGIN MCTrackPoints _____________________________
std::vector< std::vector < MCTrackObject* >* >*
MCTrackPoints( std::vector < std::vector<HitCoordinate*>* > const &evtData)
{
  std::cout << "<INFO> Extracting MC tracks for " << evtData.size()
	    << " events.\n";
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
#if (PRINT_TRACKS_DEBUG_INFO > 0)
    std::cout << "\t<-I-> Event " << e
	      << " Contains " << (numTracks + 1)
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

////////////// Perform filtering ////////
void performFilter( size_t plength, size_t area,
		    size_t MinResponce,
		    float lambda, float tol, size_t gapSizeInput,
		    std::string const &OutFileName,
		    int firstEvt, int lastEvt)
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

  // Extract event data and the detector map from the parameter,
  // simulations and digitization files.
  /*
   * For now only the stt tubes are loaded into "detNodes". This needs
   * to be fixed by including MVD geometry into the same grid as well.
   */
  
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

  //_______________________ Extending the grid with virtual nodes __________
  TNtuple Layers("LayerLimits","Layer Limits.","x:y:det_z:z");
  TNtuple Sections("SectionsLimits","Section Limits.","x:y:det_z:z");
  // Isolate Sector and Layer Limits
  isolateSectorAndLayerLimits(gr, Sections, Layers);
  Sections.Write();
  Layers.Write();

  // Collect points in the original grid
  TNtuple* OrigGrid = GridToNtuple(detNodes, "OrigGridCoord");
  OrigGrid->SetMarkerStyle(8);
  OrigGrid->SetMarkerSize(0.2);
  OrigGrid->SetMarkerColor(kBlack);
  OrigGrid->Write();

  /*
   * Add virtual nodes to the tubes grid. Virtual nodes between layers
   * with different slopes.
   */
   std::cout << "<INFO> Fix neighbouring before extension.\n";
  fixNeighboring(gr);
  std::vector < GridNode > VNodes;
  ///Compute_Add_VirtualNodes_Neigbor2(gr, VNodes);// Old method
  Compute_Add_VirtualNodes_Neigbor(gr, VNodes);
  
  // Collect and write virtual nodes to the output file
  TNtuple* virtualTubes = GridToNtuple(VNodes, "VirtualNodes");
  virtualTubes->SetMarkerStyle(8);
  virtualTubes->SetMarkerSize(0.2);
  virtualTubes->SetMarkerColor(kMagenta);
  virtualTubes->Write();
    
  //_________________ EXtend the grid with virtual nodes _________
  /* 
   * Extend the grid with new virtual points and fix the missing
   * neighboring relations
   */
  // Extend the grid by virtual nodes between the layer.
  std::cout << "<INFO> Extending the grid by " << VNodes.size()
	    << " virtual nodes between the layers.\n";
  gr.ExtendedGrid(VNodes);

  //______________________________________
  /*
   * We need to fix neigboring in case of using neigbor based grid
   * extension. Note: This procedure fixes only the missing symmetric
   * neighborhood issues. This must be done, because the grid is
   * extended by virtuals and the parent nodes do not have the
   * virtuals in their neighbour list.
   */
  std::cout << "<INFO> Fix neighbouring after extension.\n";
  fixNeighboring(gr);

  // Collecte and write
  TNtuple* extendedGrid = GridToNtuple(gr.m_grid, "ExtendedGrid");
  extendedGrid->SetMarkerStyle(7);//8
  extendedGrid->SetMarkerSize(0.3);
  extendedGrid->SetMarkerColor(17);//41
  extendedGrid->Write();

  std::cout << "<INFO> Total number of tubes after extension = "
	    << gr.GetNumNodes() << ".\n";
  
  // Delete allocated memory
  delete (OrigGrid);
  delete virtualTubes;
  delete extendedGrid;
  //___________ Done extending the grid _______________
 
 //___________________________________________________________________
   size_t NumOrientations = 0;
  size_t NumConnComp = 0;
  std::vector<float> OrientationsList;
  
  OrientationsList.push_back(0.0/tol);
  // OrientationsList.push_back(0.2/tol);
  //OrientationsList.push_back(0.4/tol);
  OrientationsList.push_back(0.6/tol);
  //OrientationsList.push_back(0.8/tol);
  //OrientationsList.push_back(1.0/tol);
  OrientationsList.push_back(1.2/tol);
  //OrientationsList.push_back(1.4/tol);
  //OrientationsList.push_back(1.6/tol);
  //OrientationsList.push_back(1.8/tol);
  OrientationsList.push_back(1.9/tol);
  //OrientationsList.push_back(2.0/tol);
  //OrientationsList.push_back(2.2/tol);
  //OrientationsList.push_back(2.4/tol);
  //OrientationsList.push_back(2.6/tol);
  //OrientationsList.push_back(2.8/tol);
  //OrientationsList.push_back(3.0/tol);
  //OrientationsList.push_back(3.1/tol);
  
  //OrientationsList.push_back(0.3/tol);
  // OrientationsList.push_back(0.5/tol);
  // OrientationsList.push_back(0.7/tol);
  //OrientationsList.push_back(0.9/tol);
  // OrientationsList.push_back(1.1/tol);
  // OrientationsList.push_back(1.3/tol);
  //OrientationsList.push_back(1.5/tol);
  // OrientationsList.push_back(1.7/tol);
  // OrientationsList.push_back(1.9/tol);
  //OrientationsList.push_back(2.1/tol);
  // OrientationsList.push_back(2.3/tol);
  // OrientationsList.push_back(2.5/tol);
  //OrientationsList.push_back(2.7/tol);
  // OrientationsList.push_back(2.9/tol);
  
   std::sort(OrientationsList.begin(), OrientationsList.end());

  unsigned int totalnumEvt = Hit_coords->size();
  // Start the timer.
  timer.Start();

  /* Fill the grid with the current hit points and process. Handles
     each event separately.*/
  // Event loop
    for(size_t k = 0; k < Hit_coords->size(); ++k) {
    // Data for the current event
    std::cout << "\n\t<INFO> Processing event: " << k
	      << "\n\n";
    std::vector<HitCoordinate*> const *dd = 0;
    dd = Hit_coords->at(k);
    if(dd) {
      gr.FillGrid(*dd);
    }
    // gr.WriteGrid("eventGridDump.txt");

    //____ Perform Attribute space connected component analysis.
    std::vector< std::set<int>* >* connectedComp = 0;
    
    /* _____ Determine orientations by pre-sensing the grid. */
    //NumOrientations = computeSoftOrienAttributes(gr, lambda, tol, gapSizeInput);
    
    /* ___ Compute orientation based attribute space. Soft */
    //NumOrientations = computeSoftOrienAttributes(gr, lambda, tol, gapSizeInput, &OrientationsList);

    /*
     * Connected components in Attribute Space (Slow but works, needs
     * to be changed to allow overlaps and shared nodes)
     */
    // connectedComp = AttSpaceConnectedComp(gr, MinResponce);
    
    /* ____ Compute orientation attribute layer based */
    // ONLY One layer at a time is processed
    
    // Orientation attr. Alternative method.
    ////// compOrientAttLayerBased_Alt(gr, tol, true);
    
    // Outer -> Inner
  compOrientAttLayerBased_Local(gr, lambda, tol, gapSizeInput, true);
    
    // Inner-> Outer
    //// compOrientAttLayerBased_Local(gr, lambda, tol, gapSizeInput, false);

    // Determine connected components
    connectedComp = AttSpaceConCompLayerBasedSoft(gr, MinResponce);

    //////========== Does not need pre-determined orientations( Work Todo) ====
    // connectedComp = ConCompLayerBasedLocalOrient(gr,  tol);

    // Set the number of connected components (used for plots)
    NumConnComp = connectedComp->size();

    // Print Info
    std::cout << "<INFO> Number of connected components: " << NumConnComp
	      << '\n';

    // If evaluate the error
#if(EVALUATE_ERROR > 0)
    // Determine the segmentation error. Based on total area.
    MCMatchingError *MC_match_error =
      MatchMCTracksWithConnectedComponents(MC_Tracks->at(k), connectedComp);
      
    ErrorNtuple.Fill(MC_match_error->Error_underMerge,
		     MC_match_error->Error_overMerge,
		     MC_match_error->TotalError,
		     MC_match_error->Error_underMergeNorm,
		     MC_match_error->Error_overMergeNorm,
		     MC_match_error->TotalErrorNorm);

    delete MC_match_error;
    // FIXME This prcedure is not complete yet It returns 0 (intentionaly).
    /* Evaluate error per track. Start from MC-Tracks and match to
       components. */
   std::vector < MCTrackObject* > *mcTracksCurrentEvent = MC_Tracks->at(k);
    // SOrt MC tracks increasing length
    std::sort(mcTracksCurrentEvent->begin(), mcTracksCurrentEvent->end(), greaterThanLength);
    std::vector< MCMatchingErrorStruct* > *match_error2 =
      MatchPerTrackWithMCTracks(gr, mcTracksCurrentEvent, connectedComp);
    if(match_error2 != 0) {
      for(size_t f = 0; f < match_error2->size(); ++f) {
        MCMatchingErrorStruct const *erObj = match_error2->at(f);
        ErrorNtuplePerTrack.Fill(static_cast<float>(erObj->isNotmatched),
                                 static_cast<float>(erObj->BestMatchMCLength),
                                 static_cast<float>(erObj->CurrentTrackLength),
                                 static_cast<float>(erObj->MCMinCurrentLength),
                                 static_cast<float>(erObj->CurrentMinMCLength),
                                 erObj->Error_underMerge,
                                 erObj->Error_overMerge,
                                 erObj->MC_a, erObj->MC_b, erObj->MC_r,
                                 erObj->MC_E, erObj->tr_a, erObj->tr_b,
                                 erObj->tr_r, erObj->tr_E);
      }//
      // Clean memory for now. Maybe better to put all lists in a main
      // list and fill the ntuple later.(HINT FIXME later)
      for(size_t r = 0; r < match_error2->size(); ++r) {
        delete match_error2->at(r);
      }
      delete match_error2;
    }
#endif// END (EVALUATE_ERROR > 0)
    //_________________ Merge tracks with MVD components (Tracklets).
    // (This does not work yet.)
    std::vector<TrackObject*>* MVDMergedTraks = MergeConnectedComponentsWithMVD(gr, connectedComp);

    //__________________ Determind eth Z-coordinate values.
    TrackZ_CoordinatesDistNorm(gr, MVDMergedTraks);

    //__________________ Collect component details _________________  
#if(WRITE_CONNECTED_COMPONENTS > 0)
    // First: event number
    // Second: number of components
    ComponentPerEvt.Fill(k, NumConnComp);
    // Collect connected components
    std::vector< GridNode > const &Ingrid = gr.m_grid;
    // Store the data for each constructed component
    for(size_t cm = 0 ; cm < connectedComp->size(); ++cm) {
      std::set<int> const* idset = connectedComp->at(cm);
      if(!idset){
	continue;
      }
      std::set<int>::iterator it;
      for( it = idset->begin(); it != idset->end(); ++it) {
     	int detID = *it;// Id of the current detector
     	int d_Index = gr.Find(detID);// Index in the grid
     	GridNode const &node = Ingrid[d_Index];	
#if(EXCLUDE_STTSKEWED_PLOT > 0)
	if( node.m_type != GridNode::STT_TYPE_SKEW ) {
#elif (EXCLUDE_VIRTUALS_IN_PLOT > 0)
          //if( node.m_type != GridNode::VIRTUAL_NODE) {// This needs to be done(FIXME)
#pragma message ("Excluding skewed or virtuals")
          //#warning "Excluding skewed or virtuals"
#endif
	  // k = event number, cm = component number,....
          ConnectedCoord.Fill(k, cm, node.m_detID, node.m_x, node.m_y, node.m_z,
                              node.m_xDet, node.m_yDet, node.m_z_Det);

	  
#if( (EXCLUDE_STTSKEWED_PLOT > 0) || (EXCLUDE_VIRTUALS_IN_PLOT > 0) )
        }// END IF skewed
#endif
      }// END FOR it = 
      if( idset->size() < MinResponce) {
	std::cout << "<WARNING> index = " << cm << " size = " << idset->size()
		  << " With IDs = ";
	for( it = idset->begin(); it != idset->end(); ++it) {
	  std::cout << *it << ", ";
	}
	std::cout << '\n';
      }
    }	
#endif//(WRITE_CONNECTED_COMPONENTS > 0)

#if(WRITE_CONNECTED_COMPONENTS_JSON > 0)

  // Write to an ASCII file
  //////////// Comments to write to file
  std::string formatt = "{\n \t \"Events\": [\n";
  //////// End of comments.
  std::ofstream OutTxtFilee;
  std::string name = "Connected_comp_met2_evt_" + std::to_string(firstEvt) + ".json";
  OutTxtFilee.open (name);
  
  if (OutTxtFilee.is_open()) {
    OutTxtFilee << formatt;
    int curevt = -1;
    int oldcurevt = curevt;
    int curid = -1;
    int oldcurid = curid;
    // Write
    OutTxtFilee << "\t\t {\n \t\t\t \"ID\": "<<firstEvt<<",\n \t\t\t\"Trajectories\": [\n";

    std::vector< GridNode > const &Ingrid = gr.m_grid;
    // Store the data for each constructed component
    for(size_t cm = 0 ; cm < connectedComp->size(); ++cm) {
      std::set<int> const* idset = connectedComp->at(cm);
      if(!idset){
	continue;
      }
      OutTxtFilee << "\t\t\t\t {\n \t\t\t\t\t \"ID\": "<<cm<<",\n \t\t\t\t\t \"Hits\": \n \t\t\t\t\t\t[\n";
      //  OutTxtFilee << "\t\t\t\t\t\t" << curht->m_detID;
      std::set<int>::iterator it;
      int ii = 0;
      for( it = idset->begin(); it != idset->end(); ++it, ii++) {
     	int detID = *it;// Id of the current detector
     	int d_Index = gr.Find(detID);// Index in the grid
     	GridNode const &node = Ingrid[d_Index];
	if( node.m_type == GridNode::VIRTUAL_NODE) {
	  printf("ALLELUIE \n");
	  continue;
	}// This needs to be done(FIXME)
	if(ii == 0)
	  OutTxtFilee << "\t\t\t\t\t\t" << detID;
	else
	  OutTxtFilee << ",\n \t\t\t\t\t\t" << detID;
      	
      }
      OutTxtFilee << "\n \t\t\t\t\t\t] \n \t\t\t\t},\n";
    }
    OutTxtFilee << "\n \t\t\t]\n \t\t}\n \t]\n }";

    
  }
  
   
 #endif//(WRITE_CONNECTED_COMPONENTS > 0)

    // Collect data for the current event into "coord" to write after
    // modifications (filtering).
    CollectGridToTree(gr, coord);
    
    // Reset grid for the next event.
    gr.ResetGrid();
    
    // Clean connected components
    if(connectedComp != 0) {
      for(size_t c = 0; c < connectedComp->size(); ++c) {
        delete connectedComp->at(c);
      }
      delete connectedComp;
    }
    // Clean mvdMerged structure
    if(MVDMergedTraks != 0) {
      for(size_t i = 0; i < MVDMergedTraks->size(); ++i) {
        delete (MVDMergedTraks->at(i));
      }
      delete MVDMergedTraks;
    }
    // End cleaning
    // GO to next event
  }//End of event loop
  timer.Stop();
  
  // Write connected components for all events
  ComponentPerEvt.Write();
  ConnectedCoord.Write();
  
  // Write coordinates ntuple
  coord.Write();
  
  //Write error estimations.
  ErrorNtuple.Write();
  ErrorNtuplePerTrack.Write();
  
  // Close open file (wel zo netjes)
  Out_Put_File.Close();
  //_________________________________________________________________
  //_______ Clean allocated memory
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
  //______________ Running Time information ____________
  
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  std::cout <<"=======================================\n"
            << "Macro finished succesfully.\n"
            << "Real time " << (rtime/totalnumEvt)
	    << " (s/Event), CPU time " << (ctime/totalnumEvt)
            << " (s/Event).\n"
            << "plength = " << plength
            << ", area = " << area
            << ", Num of orientations = " << NumOrientations
            << ", lambda = " << lambda
            << ", tolerance = " << tol
            << '\n';
  // Extract MC track values
  
  
}
//__________________ MatchMCTracksWithConnectedComponents ___________________
/* Function to evaluate matching error per image. */
MCMatchingError* MatchMCTracksWithConnectedComponents(std::vector< MCTrackObject* >  const *MCTracks,
						      std::vector< std::set<int>* >  const *connectedComp)
{
  if( (MCTracks == 0) || (connectedComp == 0) ) {
    std::cerr << "<ERROR> One of the input parameters for matching is empty.\n";
    return 0;
  }
  // Compute the total MC area (true area) for the current event.
  float TotalArea = 0;
  float Error_underMerge = 0;
  float Error_overMerge  = 0;
  float TotalError = 0;
  // Per error type normalized values
  float Error_underMergeNorm = 0;
  float Error_overMergeNorm  = 0;
  float TotalErrorNorm = 0;
  
  float mc_length = 0;

  // Determine the total area of the all tracks in the current event.
  for(size_t i = 0; i < MCTracks->size(); ++i) {
    MCTrackObject const *MCtrack = MCTracks->at(i);
    TotalArea += (MCtrack->m_STT_Component).size();
  }
  // FORALL elements in "connectedComp"(T_j) find R_k in "MCTracks"
  // such that (T_j intersection R_k) is maximum.
  std::vector<int>::iterator it;
  std::set<int>::iterator compIt;
  // Component loop(T_j)
  for(size_t i = 0; i < connectedComp->size(); ++i) {
    float matchValue = std::numeric_limits<float>::min();
    // Current connected component to analyse
    std::set<int> const *Cur_comp = connectedComp->at(i);
    std::vector<int> Cur_Comp_list;
    for(compIt = Cur_comp->begin(); compIt != Cur_comp->end(); ++compIt) {
      int id = *compIt;
      // Only Real stt tubes, MC does not know anything of the
      // virtuals
      if( (id < START_VIRTUAL_ID ) && (id >= 0) ) {
	Cur_Comp_list.push_back(id);
      }
    }
    std::sort(Cur_Comp_list.begin(), Cur_Comp_list.end());
    // Find which MC track has the largest overlap (R_k)
    for(size_t j = 0; j < MCTracks->size(); ++j) {
      // Current MC track
      MCTrackObject const *MCtrack = MCTracks->at(j);
      // Stt component of the current MC track.
      std::vector<int> MCSttComp(MCtrack->m_STT_Component);
      std::sort(MCSttComp.begin(), MCSttComp.end());
      // Intersection result output
      std::vector<int> IntersectionList( (Cur_Comp_list.size() + MCSttComp.size()), 0 );

      // Determine the sequence intersection. Only real tubes.
      it = std::set_intersection( Cur_Comp_list.begin(), Cur_Comp_list.end(),
				  MCSttComp.begin(), MCSttComp.end(),
				  IntersectionList.begin());

      // Resize. it points to the position of the last common element
      IntersectionList.resize(it - IntersectionList.begin());
      if(static_cast<float>(IntersectionList.size()) > matchValue) {
	matchValue = IntersectionList.size();
	// A_k
	mc_length = static_cast<float>(MCSttComp.size());
      }
    }// MC tracks loop

#if(PRINT_DEBUG_INFO_COMP_MATCH > 0)
    std::cout << "<INFO> MatchValue = " << matchValue
	      << "      <-I-> mc_length = " << mc_length
	      << '\n';
#endif
    
    // We have found ( (T_j intersection R_k) is maximum )
    Error_underMerge +=  ( ( ( mc_length - matchValue ) * matchValue )/ mc_length );
    // Normalized Under segment
    Error_underMergeNorm += ( ( ( mc_length - matchValue ) * matchValue )/ mc_length ) / TotalArea;

    // (alpha_j = Cur_Comp_list.size())
    Error_overMerge  += (static_cast<float>(Cur_Comp_list.size()) - matchValue);

    // Normalized over segment
    Error_overMergeNorm += (static_cast<float>(Cur_Comp_list.size()) - matchValue) / TotalArea;

  }// Components loop
  // compute total error for the current set and normalize by the
  // total area of the actual instance (MC Truth)
  TotalError = (Error_underMerge * Error_underMerge) +
               (Error_overMerge  * Error_overMerge);
  // Normalize by total Area of MC.
  TotalError = sqrt(TotalError)/TotalArea;
  
  // Sub terms are already normalised
  TotalErrorNorm = (Error_underMergeNorm * Error_underMergeNorm) +
                   (Error_overMergeNorm  * Error_overMergeNorm);
  TotalErrorNorm = sqrt(TotalErrorNorm);
  
#if(PRINT_DEBUG_INFO_COMP_MATCH > 0)
  std::cout << "<INFO> "
	    << "Error_underMerge " << Error_underMerge
	    << " Error_overMerge " << Error_overMerge
	    << " TotalError "        << TotalError
	    << '\n';
#endif

  // Prepare output
  MCMatchingError *outError = new MCMatchingError();
  outError->Error_underMerge = Error_underMerge;
  outError->Error_overMerge  = Error_overMerge;
  outError->TotalError       = TotalError;
  // Normalised version
  outError->Error_underMergeNorm = Error_underMergeNorm;
  outError->Error_overMergeNorm  = Error_overMergeNorm;
  outError->TotalErrorNorm       = TotalErrorNorm;
  // Return
  return outError;
}
//__________________ END MatchMCTracksWithConnectedComponents _______________
//___________________________ MatchPerTrackWithMCTracks ____________________
/* Function to evaluate error per track. MC is used as ground truth. */
std::vector< MCMatchingError* >* MatchPerTrackWithMCTracks(CoordGrid const &hitMap,
                                                           std::vector < MCTrackObject* > const *MCTracks,
                                                           std::vector< std::set<int>* > const *connectedComp)
{
  if( (MCTracks == 0) || (connectedComp == 0) ) {
    std::cerr << "<ERROR> One of the input parameters for matching is empty.\n";
    return 0;
  }
  // Create output parameter
  std::vector< MCMatchingError* >* outPutPar = 0;
  outPutPar = new std::vector< MCMatchingError* >();
  if( outPutPar == 0) {
    std::cerr << "<ERROR> Could not allocate memory for output list in MatchPerTrackWithMCTracks.\n"
              << " The program is exiting.(return 0)\n";
    delete outPutPar;
    return 0;
  }
  std::cout << "<INFO> input MC has " << MCTracks->size() << " members and components has "
            << connectedComp->size() << '\n'
            << "\tfirst MC has " << ((MCTracks->at(0))->m_STT_Component).size()
            << " elements and last has " << ((MCTracks->at(MCTracks->size() - 1))->m_STT_Component).size()
            << '\n';
  std::vector<int> matchedTrackletIndices;
  std::vector<int>::iterator FindIndexIt;
  // For all MC-tracklets find the best match in connected componets
  // list.(R_k interset T_j) maximum
  std::vector<int>::iterator it;
  std::set<int>::iterator compIt;
  // MC tracks loop
  for(size_t i = 0; i < MCTracks->size(); ++i) {
    // Current MC track(R_k)
    MCTrackObject const *MCtrack = MCTracks->at(i);
    // Stt component of the current MC track.
    std::set<int> MCSttComp((MCtrack->m_STT_Component).begin(), (MCtrack->m_STT_Component).end());
    // Not empty MC-track (How is this possible??, empty MC-Tracks.)
    if( (MCSttComp.size() > 0) ) {
      // Local variables for comp....
      int    matchTrackIndex = -1;
      float  match_length = 0;
      float  matchValue = std::numeric_limits<float>::min();
      // Tracklets loop (connected components)
      for(size_t j = 0; j < connectedComp->size(); ++j) {
        FindIndexIt = std::find(matchedTrackletIndices.begin(), matchedTrackletIndices.end(), j);
        if( FindIndexIt == matchedTrackletIndices.end() ) {// Was not assigned before
          // Current component (T_j)
          std::set<int> const *Cur_comp = connectedComp->at(j);
          std::set<int> Cur_Comp_list;
          for(compIt = Cur_comp->begin(); compIt != Cur_comp->end(); ++compIt) {
            int id = *compIt;
            // Only Real stt tubes, MC does not know anything of the
            // virtuals
            if( (id < START_VIRTUAL_ID ) && (id >= 0) ) {
              Cur_Comp_list.insert(id);
            }
          }//
          // Set intersection
          std::vector<int> IntersectionList( (Cur_Comp_list.size() + MCSttComp.size()), 0 );
          // Determine the sequence intersection. Only real tubes.
          it = std::set_intersection( Cur_Comp_list.begin(), Cur_Comp_list.end(),
                                      MCSttComp.begin(), MCSttComp.end(),
                                      IntersectionList.begin());
          // Resize. "it" points to the position of the last common element
          IntersectionList.resize(it - IntersectionList.begin());
          /* FIXME hier veranderd om exacte match te vinden */
          // if(static_cast<float>(IntersectionList.size()) > matchValue) {
          if(static_cast<float>(IntersectionList.size()) == MCSttComp.size()) {
            matchValue = static_cast<float>(IntersectionList.size());
            // Area of matched tracklet (a_j)
            match_length = static_cast<float>(Cur_Comp_list.size());
            // Index of tracklet with largest intersection
            matchTrackIndex = j;
          }
        }// End was not assigned before
      }// Components loop (T_j)
      // We found the best matching tracklet for current MC track
      if(matchTrackIndex >= 0) {
        /* Met of zonder teruglegging */
        // Mark as assigned(to avoid multiple assignements)
        ////matchedTrackletIndices.push_back(matchTrackIndex);
        
        std::set<int> const *bestMatchComponent = connectedComp->at(matchTrackIndex);
        std::vector<int> bestMatchlist;
        for(compIt = bestMatchComponent->begin(); compIt != bestMatchComponent->end(); ++compIt) {
          int id = *compIt;
          // Only Real stt tubes, MC does not know anything of the
          // virtuals
          if( (id < START_VIRTUAL_ID ) && (id >= 0) ) {
            bestMatchlist.push_back(id);
          }
        }
        std::sort(bestMatchlist.begin(), bestMatchlist.end());
        // Determine differences
        // All x in MC and not in Tracklet (R_k - T_j)
        std::vector<int> MCdiffCurComp((bestMatchlist.size() + MCSttComp.size()), 0);
        it = std::set_difference( MCSttComp.begin(), MCSttComp.end(),
                                  bestMatchlist.begin(), bestMatchlist.end(),
                                  MCdiffCurComp.begin());
        MCdiffCurComp.resize(it - MCdiffCurComp.begin());
        // All x in Tracklet and not in MC
        std::vector<int> CurCompdiffMC((bestMatchlist.size() + MCSttComp.size()), 0);
        it = std::set_difference( bestMatchlist.begin(), bestMatchlist.end(),
                                  MCSttComp.begin(), MCSttComp.end(),
                                  CurCompdiffMC.begin());
        CurCompdiffMC.resize(it - CurCompdiffMC.begin());
        // Determine curvature parameters for tracklet and MC
        /*______ Curvature for MC points _____*/
        /* HIER was je mee bezig. Het ziet er OK uit. Andere functie
           om MC data te fitten is gemaakt. */
        //MCtrack->print();// FIXME hier Heb je veranderd(FIXME, 24 nov)
        //std::vector<int> MCSttCompList(MCSttComp.begin(), MCSttComp.end());
        CurvatureParameters mcCurvPars;
        /* De functie hieronder niet voor MC gebruiken. Hij gebruikt
         * de in de grid beschikbare coordinaten en niet de MC punten.*/
        //ComputeCurvatureForListOfnodes(hitMap, MCSttCompList, mcCurvPars);

        //int CircleFit( std::vector<point3D> const &points, CurvatureParameters &curvature);
        std::vector<point3D> const &AllSttpoints = MCtrack->m_pointSTTCoordList;
        int blabla = CircleFit(AllSttpoints, mcCurvPars);
        std::cout << "blabla = " << blabla << '\n';

        /*________ Determine curvature for bet match tracklet _____*/
        CurvatureParameters componentCurvPars;
        // FIXME: Maybe better to include virtuals. Better gradual ....
        ComputeCurvatureForListOfnodes(hitMap, bestMatchlist, componentCurvPars);
        // Create error object for current sub-path
        MCMatchingError *erroObject = new MCMatchingError();
        erroObject->BestMatchMCLength  = MCSttComp.size();
        erroObject->CurrentTrackLength = bestMatchlist.size();
        // Set diffs
        erroObject->MCMinCurrentLength = MCdiffCurComp.size();
        erroObject->CurrentMinMCLength = CurCompdiffMC.size();
        // Per track values according to definitions in the book.
        /* {A_k - (T_j ^ R_k)} * (T_j ^ R_k) / A_k */
        erroObject->Error_underMerge = -1;// Not Computed here
        erroObject->Error_overMerge  = -1;// Not Computed here
        // Store curvature parameters. MC
        erroObject->MC_a = mcCurvPars.m_a;
        erroObject->MC_b = mcCurvPars.m_b;
        erroObject->MC_r = mcCurvPars.m_r;// is 1/r
        erroObject->MC_E = mcCurvPars.m_E;
        // Matched tracklet
        erroObject->tr_a = componentCurvPars.m_a;
        erroObject->tr_b = componentCurvPars.m_b;
        erroObject->tr_r = componentCurvPars.m_r;// is 1/r
        erroObject->tr_E = componentCurvPars.m_E;
        erroObject->TotalError       = -1;// Not Computed here
        // Add to output list
        outPutPar->push_back(erroObject);
      }
      else{// not (matchTrackIndex >= 0)
        /* There was not a matching tracklet. Totally mis matched.*/
        MCMatchingError *noMatchError = new MCMatchingError();
        noMatchError->isNotmatched = 1.0;
        noMatchError->BestMatchMCLength  = MCSttComp.size();
        noMatchError->CurrentTrackLength = 0;
        noMatchError->MCMinCurrentLength = MCSttComp.size();
        noMatchError->CurrentMinMCLength = 0;
        // Add to output list
        outPutPar->push_back(noMatchError);
        printIntList(MCSttComp);
      }
      // DEBUG INFO
#if(PRINT_ERROR_EXTRA_DEBUG_INFO > 0)
      std::cout << "<DEBUG> MC track " << i << " best match index = " << matchTrackIndex
                << '\n';
#endif
      //___ End Debug
    }//END if(MCSttComp.size() != 0) 
    //Go to next MC-Track
  }//MC tracks loop R_j
  // Debug info before return
#if(PRINT_ERROR_EXTRA_DEBUG_INFO > 0)
  std::cout << "<INFO> Size of per track error list is " << outPutPar->size()
            << '\n';
#endif
  return outPutPar;
}
//_______________________ END MatchPerTrackWithMCTracks ___________________
// DEBUG functie, mag weg later. Voegt niets toe.
void printIntList(std::vector<int> const &inp)
{
  std::cout << "Print int vector with " << inp.size()
            << " members.\n";
  for(size_t i = 0; i < inp.size(); ++i) {
    std::cout << inp[i] << ", ";
  }
  std::cout << '\n';
}
void printIntList(std::set<int> const &inp)
{
  std::cout << "Print int set with " << inp.size()
            << " members.\n";
  std::set<int>::iterator it;
  for(it = inp.begin(); it != inp.end(); ++it) {
    std::cout << (*it) << ", ";
  }
  std::cout << '\n';
}
// DEBUG, mag weg later
