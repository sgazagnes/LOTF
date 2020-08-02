/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
/**
 * Example macro created by: M. Babai
 * Find and collect STT, MDV hits.
 */
// C and C++ standard headers
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <utility>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <sstream>
#include <algorithm>

// PANDA and RAIR Root
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "FairParRootFileIo.h"
#include "PndEmcGeoPar.h"
#include "FairBaseParSet.h"
#include "FairParAsciiFileIo.h"
#include "PndGeoSttPar.h"
#include "PndSttMapCreator.h"
#include "PndSttTube.h"
#include "PndSttPoint.h"
#include "PndSttHit.h"
#include "PndSdsHit.h"
#include "PndSdsMCPoint.h"
//#include "PndTrkSttAdjacencies.h"
#include "PndSttSkewedHit.h"
#include "PndTrack.h"
#include "PndGeoHandling.h"
#include "FairLink.h"

// Root headers
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TClonesArray.h"
#include "TGeoManager.h"
#include "TCanvas.h"
#include "TStopwatch.h"

// Local headers
#include "CollectSttMvdPoints.h"
#include "hitcoordinate.h"
#include "gridNode.h"
#include "CoordGrid.h"
#include "pathopen.h"
#include "utilfunctions.h"
#include "logc.h"

#define INCLUDE_STT_POINTS  1
#define INCLUDE_MVD_POINTS  0
#define DEBUG_PRINT     0
#define WRITE_TO_ASCII_FILE 0
#define WRITE_TO_JSON_FILE 0
#define READ_NEIGHBORS 0

void CollectGridToTree( CoordGrid const &gr, TNtuple &out)
{
  std::vector< GridNode > const& grid = gr.m_grid;
  
  for(size_t j = 0; j < grid.size(); ++j) {
    GridNode const& curPoint = grid[j];
    if( (curPoint.m_active)
	//&& (curPoint.m_type == GridNode::VIRTUAL_NODE)
	) {// If fired
      out.Fill(curPoint.m_x, curPoint.m_y, curPoint.m_z,
	       curPoint.m_xDet, curPoint.m_yDet, curPoint.m_z_Det);
    }
  }
}

/*
 * Collect all STT (and other) detector information.(space position
 * and neighboring info). This is a helper function to extract
 * detector data from the tube data available inside the simulation
 * output file.
 */
void CollectSttDetecorCoords(TClonesArray const &TubeArray, // In par
                             std::vector < GridNode > &detNodes)// Out par
{
  info("Total number of tubes = %d ",TubeArray.GetEntries());
          
  // Tubes numbering starts with 1 ???
  PndSttTube *tube = NULL;
  size_t sLeft, sRight, nSkewed, nPara;
  nSkewed = nPara = sLeft = sRight = 0;

  //  std::string header = "ID,x,y,z,wire_direction_x,wire_direction_y,wire_direction_z,halflength,type,sector,layer,layer_lim,sector_lim,num_neighbors\n";
   // std::string header = "ID,neighbors\n";
  //////// End of comments.
  // std::ofstream OutTxtFile;
  // OutTxtFile.open ("TubeGeometry2.csv");

  //  OutTxtFile << header;

#if ( READ_NEIGHBORS > 0 )
  int **tube_neigh = (int **) malloc(TubeArray.GetEntries() * sizeof(int*));
  int cur_ind[100];
  char buffer[1024] ;
  char *record,*line;
  int *num_neigh = (int *) malloc(TubeArray.GetEntries() * sizeof(int));
  FILE *fstream = fopen("PANDA_STT_Neighbours/PANDA_STT_tol01.csv","r");
  //FILE *fstream = fopen("neighbours.csv","r");
  if(fstream == NULL)
    {
      printf("\n file opening failed \n");
    }
  int i = 0;
  while((line=fgets(buffer,sizeof(buffer),fstream))!=NULL)
    {
      record = strtok(line,",");
      int count = 0;
      while(record != NULL)
	{
	  //  printf("record : %s",record) ;    //here you can put the record into the array as per your requirement.
	  // printf("%d\n", atoi(record));
	  if(atoi(record) != 0){
	    cur_ind[count++] = atoi(record) ;
	  }
	  record = strtok(NULL,",");
	}
      //	 printf("%d i, count %d \n", i, count);
      tube_neigh[i] =  (int *) malloc((count)*sizeof(int));
      if (i+1 != cur_ind[0])
	printf("TOBE CHANGEEEEDD \n\n\n");
      for (int j = 1; j < count; j++){
	tube_neigh[i][j-1] = cur_ind[j];
	//  printf("%d, cur_ind[j] %d \n",cur_ind[0],  cur_ind[j]);
      }
      num_neigh[i] = count-1;
      ++i ;

      //	 printf("\n\n\n") ;    //here you can put the record into the array as per your requirement.

    }
  //TubeArray.GetEntries()
#endif
  // Note tube count starts at 1.
  for(int itube = 1; itube <=  TubeArray.GetEntries(); itube++) {
    tube = (PndSttTube*) TubeArray.At(itube);
    GridNode node;
    node.m_detID  = tube->GetTubeID();// This may change later
    //  OutTxtFile << node.m_detID <<",";

    
    node.m_Orig_detID  = tube->GetTubeID();// This id remains always the same
    node.m_active = false;

    TVector3 nodePosition = tube->GetPosition();
    node.m_x      = nodePosition.X();
    node.m_y      = nodePosition.Y();
    node.m_z      = nodePosition.Z();
    // OutTxtFile << node.m_x <<"," <<node.m_y <<","<< node.m_z <<",";

    
    node.m_xDet   = node.m_x;
    node.m_yDet   = node.m_y;
    
    node.m_WireDirection = tube->GetWireDirection();
    //  OutTxtFile << node.m_WireDirection.X() <<"," <<node.m_WireDirection.Y() <<","<< node.m_WireDirection.Z() <<",";
    // // OutTxtFile << node.m_WireDirection.rho() <<"," <<node.m_WireDirection.Theta() <<","<< node.m_WireDirection.Phi() <<",";


    node.ComputeSlope();
    node.m_halfLength    = tube->GetHalfLength();
    // OutTxtFile << node.m_halfLength <<",";
    if( tube->IsParallel() ) {
      node.m_type  = GridNode::STT_TYPE_PARA;
      node.m_z_Det = 0.0;
      nPara++;
    }
    else {
     node.m_type = GridNode::STT_TYPE_SKEW;
     node.m_z_Det = node.m_z;
     nSkewed++;
    }
    // OutTxtFile << node.m_type <<",";

    
    node.m_Sector = tube->GetSectorID();
    // OutTxtFile << node.m_Sector <<",";

    node.m_Layer  = tube->GetLayerID();
    //OutTxtFile << node.m_Layer <<",";

    node.m_LayerLimit = tube->IsLayerLimit();
    //OutTxtFile << node.m_LayerLimit <<",";
    
    if( tube->IsSectorLimit() == 1 ) {
      node.m_SectorLimit = 1;
      sRight++;
    }
    else if( tube->IsSectorLimit() == -1 ) {
      node.m_SectorLimit = -1;
      sLeft++;
    }
    else {
      node.m_SectorLimit = 0;
    }
    //OutTxtFile << node.m_SectorLimit <<",";

    // Collect and set neighbor detector id's.
    // Fetch index list of neighors for current detector.
     TArrayI neighbs = tube->GetNeighborings();
    // Neighbor Loop
    // OutTxtFile <<  neighbs.GetSize() <<"\n";
    // printf("%d", itube);
     for(int n = 0; n < neighbs.GetSize(); n++) {
      int neigh_index = neighbs.At(n);
      
      //  printf("%d, %d", neigh_index, neigh_index2);
      PndSttTube *neighb_tube = (PndSttTube*) TubeArray.At(neigh_index);
      int n_Id2;
      /* if(n < num_neigh[itube -1]){
      int neigh_index2 = tube_neigh[itube-1][n];//neighbs.At(n);
      PndSttTube *neighb_tube2 = (PndSttTube*) TubeArray.At(neigh_index2);
      n_Id2 = neighb_tube2->GetTubeID();
      (node.m_neighbors).push_back(n_Id2);

      }else
      n_Id2 = -1;*/
      
      int n_Id = neighb_tube->GetTubeID();
      // printf("%d, %d \n", n_Id, n_Id2);

      // OutTxtFile <<  n_Id <<",";
       (node.m_neighbors).push_back(n_Id);
    }// END neighbor loop

     //OutTxtFile << "\n";*/
    // Add node to the node list
    detNodes.push_back(node);
  }// Tube loop
  info("Number of sector limits = %d, right-limit = %d, left-limit = %d, number of skewed = %d, number of para = %d, total = %d", (sRight + sLeft), sRight, sLeft, nSkewed, nPara, nSkewed+nPara);
  
}

//___________ Main function in this macro __________________
std::vector < std::vector<HitCoordinate*>* >*
CollectSttMvdPoints( std::vector < GridNode >& detNodes, TFile &OutFile,
		     int firstEvt, int lastEvt)
{
   TStopwatch timer;
  timer.Start();

  // --- MC INFO --------------------------------------------------
  // Sim File
  TFile sf("./rootfiles/evtcomplete_sim.root","READ");
  
  // Digi File
  TFile df("./rootfiles/evtcomplete_digi.root","READ");

  // Get file of parameters
  TFile parfile("./rootfiles/evtcomplete_par.root", "READ");

  //______________________________________________
  size_t mvdStrCntMC, mvdPixCntMC, SttCntMC;
  mvdStrCntMC = mvdPixCntMC = SttCntMC = 0;

  // Cartesian coordinates
  TNtuple MCposition  ("MCpos","Cartesian MCpos All","mcx:mcy:mcz:evtNum");
  TNtuple position    ("Pos","Cartesian pos All","x:y:z:evtNum");

  // Polar (Cylindrical)
  TNtuple MCposPolar  ("MCposPolar","Polar MCpos All","mcr:mcthetaRad:mcthetaDeg:z:evtNum");
  TNtuple posPolar    ("PosPolar"  ,"Polar pos All","r:thetaRad:thetaDeg:evtNum");

  // Stt Data only
  TNtuple STTposition   ("STTPos","Cartesian pos STT","x:y:z:isochrone:evtNum");
  TNtuple STTposPolar   ("STTPosPolar","Polar pos STT","r:thetaRad:thetaDeg:z:isochrone:evtNum");

  // MC STT only
  TNtuple STTMCposition   ("STTMCPos","Cartesian MC pos STT" ,"mcx:mcy:mcz:evtNum");
  TNtuple STTMCposPolar   ("STTMCPosPolar","Polar MC pos STT","mcr:mcthetaRad:mcthetaDeg:z:evtNum");

  //______________________________________________

  // Simu tree
  TTree* simuTr = (TTree *) sf.Get("pndsim");

  // Stt MC points
  TClonesArray *SttMCPointAr = new TClonesArray("PndSttPoint");
  simuTr->SetBranchAddress("STTPoint", &SttMCPointAr);

  // MVD MC points
  TClonesArray *MvdMCPointAr = new TClonesArray("PndSdsMCPoint");
  simuTr->SetBranchAddress("MVDPoint", &MvdMCPointAr);
  
  // -------------------------------------------------------------- 
  // Digi tree
  TTree* digiTr = (TTree *) df.Get("pndsim");

  //__________________ STT Hits
  TClonesArray *SttHitArray = new TClonesArray("PndSttHit");
  digiTr->SetBranchAddress("STTHit", &SttHitArray);

  //__________________ MVD hits
  // MVD Pixel hits
  TClonesArray *MvdPixelHitArray = new TClonesArray("PndSdsHit");
  digiTr->SetBranchAddress("MVDHitsPixel", &MvdPixelHitArray);

  // MVD Strip hits
  TClonesArray *MvdStripHitArray = new TClonesArray("PndSdsHit");
  digiTr->SetBranchAddress("MVDHitsStrip", &MvdStripHitArray);

  // load Fair base pars.
  FairBaseParSet *fairbasepar = 0;
  fairbasepar = (FairBaseParSet*) parfile.Get("FairGeoParSet");
  // Load geoManager
  TGeoManager *geoMan = 0;
  geoMan = gGeoManager;
  // if (!gGeoManager) {
  //   fairbasepar = (FairBaseParSet*) parfile.Get("FairBaseParSet");
  //   geoMan = gGeoManager;
  // }
  if(!geoMan) {
    error("Could not find valid GeoManager. Abort now!.");
    exit(1);
  }
  
  fairbasepar->GetName();
  PndGeoHandling geoH(1583944737 , "./rootfiles/evtcomplete_par.root"); //geo 2 1572086365 // geo 1 1583944737
  geoH.FillSensorMap();
  
#if (DEBUG_PRINT > 0)
  fairbasepar->Print();
#endif

  PndGeoSttPar *SttParameters = 0; 
  // Retrieve stt parameters
  SttParameters = (PndGeoSttPar*) parfile.Get("PndGeoSttPar");
  
  //  -------------------------------------------------------------
  PndSttMapCreator *mapper = 0;
  mapper = new PndSttMapCreator(SttParameters);
  
  if(mapper == 0) {
    error("Could not Create the STT mapper");
    exit(1);
  }

  TClonesArray *SttTubeArray = 0;
  SttTubeArray = mapper->FillTubeArray();
  
  if( SttTubeArray == 0) {
    error(" Could not get  STT TubeArray");
    exit(2);
  }
  
  info("Total number of events in digi tree = %d", digiTr->GetEntriesFast());

  // pair to hold the polar coordinates output.
  std::pair<float, float> r_Theta;
  float theta_deg;
  r_Theta.first = r_Theta.second = theta_deg = 0.00;


  // Vector to store hits
  //std::vector < std::vector<HitCoordinate*>* > hitcoords;
  std::vector < std::vector<HitCoordinate*>* >* outHit_coords
    = new std::vector < std::vector<HitCoordinate*>* >();

  //______________  Event loop. ________________________
  int startEvent = 0;
  int LastEvent  = 0;
  if( (firstEvt < 0) && (lastEvt < 0) ) {
    LastEvent = digiTr->GetEntriesFast();
  }
  else {
    startEvent = firstEvt;
    LastEvent = lastEvt;
  }

  for(int e = startEvent; e < LastEvent; ++e)
  {
    simuTr->GetEntry(e);
    digiTr->GetEntry(e);

    // New vector for the current event.
    std::vector<HitCoordinate*>* cuEvtData = new std::vector<HitCoordinate*>();

#if ( INCLUDE_STT_POINTS > 0 )
    // Fetch STT data points.
    for(int i = 0; i < SttHitArray->GetEntriesFast(); ++i)
    {
      HitCoordinate* currentHit = new HitCoordinate();

      // Set the event number
      currentHit->m_EvtNum = e;

      // Fetch the stt hit point
      PndSttHit *pSttHit = (PndSttHit*) SttHitArray->At(i);
      
      //___ *********** HIT/POINT POSITION MC ***********
      int refindex = pSttHit->GetRefIndex();
      PndSttPoint *pSttMCPoint = (PndSttPoint*) SttMCPointAr->At(refindex);

#if (DEBUG_PRINT > 0)
      debug("TrackID = %d", pSttMCPoint->GetTrackID());
      // pSttMCPoint->Print("");
#endif
      // Set the trackID for the current point
      currentHit->m_trackID = pSttMCPoint->GetTrackID();
      
      // MC point coordinates
      TVector3 mcposition;
      pSttMCPoint->Position(mcposition);

#if (DEBUG_PRINT > 0)
      mcposition.Print();
#endif
      // MC data Convert tot Polar coordinates
      theta_deg = Cartesian_To_Polar(mcposition.X(), mcposition.Y(), r_Theta);

      // MC position Polar
      MCposPolar.Fill   (r_Theta.first, r_Theta.second, theta_deg, mcposition.Z(), e);
      // Stt Only
      STTMCposPolar.Fill(r_Theta.first, r_Theta.second, theta_deg, mcposition.Z(), e);
      
      // MC position Cartesian
      MCposition.Fill   (mcposition.X(), mcposition.Y(), mcposition.Z(), e);
      // Stt Only
      STTMCposition.Fill(mcposition.X(), mcposition.Y(), mcposition.Z(), e);

      SttCntMC++;

      // Update current hit
      // Polar coordinates
      currentHit->mr = r_Theta.first;
      currentHit->mtheta = r_Theta.second;
      currentHit->mthetaDeg = theta_deg;
      // Cartesian
      currentHit->mx = mcposition.X();
      currentHit->my = mcposition.Y();
      currentHit->mz = mcposition.Z();

      // ----------------------------------------------------------
#if (DEBUG_PRINT > 0)
      debug("id = %d", pSttHit->GetTubeID());
#endif
      //_______ ********* Not MC Hit coordinates ********
      // Which tube
      PndSttTube *pSttTube = (PndSttTube *) SttTubeArray->At(pSttHit->GetTubeID());
      
      // drift radius
      Double_t isochrone   = pSttHit->GetIsochrone();

      // Coordinates of the center of the tube
      TVector3 tubecenter  = pSttTube->GetPosition();

      /////// DEBUG DEBUG
      // std::cout << " TubeID = " << pSttHit->GetTubeID()
      //           << " number of neigbours = " << (pSttTube->GetNeighborings()).GetSize()
      //           << ", ";
      // tubecenter.Print();
      // std::cout << '\n';
      /////// DEBUG DEBUG

      // Polar coordinates Data.
      theta_deg = Cartesian_To_Polar(tubecenter.X(), tubecenter.Y(), r_Theta);

      // means Skew (NOT Z parallel)
      // if ( pSttTube->GetWireDirection().Z() != 1)
      posPolar.Fill   (r_Theta.first, r_Theta.second, theta_deg, e);
      STTposPolar.Fill(r_Theta.first, r_Theta.second, theta_deg, tubecenter.Z(), isochrone, e);
      
      // Cartesian postion of the STT hit.
      position.Fill   (tubecenter.X(), tubecenter.Y(), tubecenter.Z(), e);
      STTposition.Fill(tubecenter.X(), tubecenter.Y(), tubecenter.Z(), isochrone, e);

      // Update point
      currentHit->r = r_Theta.first;
      currentHit->theta = r_Theta.second;
      currentHit->thetaDeg = theta_deg;

      currentHit->x = tubecenter.X();
      currentHit->y = tubecenter.Y();
      currentHit->z = tubecenter.Z();// Center of the tube

      currentHit->isochrone = isochrone;
      currentHit->type = HitCoordinate::STT_TYPE;
      
      //currentHit->m_detID = pSttTube->GetTubeID();
      currentHit->m_detID = pSttHit->GetTubeID();
      
      // Add to container
      cuEvtData->push_back(currentHit);
    }// END STT hits loop
#endif// end of STT Hit Data inclusion

#if ( INCLUDE_MVD_POINTS > 0 )
    //=========== MVD 3D points
    //+++ Pixel MVD
    for( int i= 0; i< MvdPixelHitArray->GetEntriesFast(); i++)
    {
      HitCoordinate* currentHit = new HitCoordinate();

      // Set the event number
      currentHit->m_EvtNum = e;

      // Get MVD Pixel hit
      PndSdsHit *pMvdPixelHit = (PndSdsHit*) MvdPixelHitArray->At(i);
      // DEBUG INFO
      //std::cout << "Sensor pixel ids " << pMvdPixelHit->GetSensorID() <<'\n' ;

      // ___ ************* Extract MC- Data
      int refindex = pMvdPixelHit->GetRefIndex();
      
      // If equal -1, then noise hits
      if( refindex != -1)
      {
        PndSdsMCPoint *sdsMCPoint = (PndSdsMCPoint*) MvdMCPointAr->At(refindex);
        TVector3 mvdHitPos = sdsMCPoint->GetPosition();
	
	// Set MC trackID
	currentHit->m_trackID = sdsMCPoint->GetTrackID();

#if (DEBUG_PRINT > 0)
        mvdHitPos.Print();
#endif
        // MC data Polar position
        theta_deg = Cartesian_To_Polar(mvdHitPos.X(), mvdHitPos.Y(), r_Theta);
        
        MCposPolar.Fill(r_Theta.first, r_Theta.second, theta_deg, mvdHitPos.Z(), e);

        // MC data Cartesian hit
        MCposition.Fill(mvdHitPos.X(), mvdHitPos.Y(), mvdHitPos.Z(), e);
        mvdPixCntMC++;

        // Update current hit
        // Polar coordinates
        currentHit->mr = r_Theta.first;
        currentHit->mtheta = r_Theta.second;
        currentHit->mthetaDeg = theta_deg;
        // Cartesian
        currentHit->mx = mvdHitPos.X();
        currentHit->my = mvdHitPos.Y();
        currentHit->mz = mvdHitPos.Z();
      }
      else{// Points to exclude
        MCposition.Fill(HIT_EXCLUSION, HIT_EXCLUSION, HIT_EXCLUSION, e);
        MCposPolar.Fill(HIT_EXCLUSION, HIT_EXCLUSION, HIT_EXCLUSION, e);
        (*currentHit) = HitCoordinate(HIT_EXCLUSION);
      }
      //_____ ************ Non MC data MVD Pixel.
      TVector3 tmpPos = pMvdPixelHit->GetPosition();
      
      // Polar coordinates.
      theta_deg = Cartesian_To_Polar(tmpPos.X(), tmpPos.Y(), r_Theta);
      
      posPolar.Fill(r_Theta.first, r_Theta.second, theta_deg, e);
      
      // Cartesian coordinates.
      position.Fill(tmpPos.X(), tmpPos.Y(), tmpPos.Z(), e);
    
      // Update point
      currentHit->r = r_Theta.first;
      currentHit->theta = r_Theta.second;
      currentHit->thetaDeg = theta_deg;

      currentHit->x = tmpPos.X();
      currentHit->y = tmpPos.Y();
      currentHit->z = tmpPos.Z();

      currentHit->type = HitCoordinate::MVD_TYPE;
      currentHit->m_detID = pMvdPixelHit->GetSensorID();
      // Add to container
      cuEvtData->push_back(currentHit);
    }// End MVD Pixel hits.
    
    //+++++++ Strip MVD
    for(int i= 0; i< MvdStripHitArray->GetEntriesFast(); ++i)
    {
      HitCoordinate* currentHit = new HitCoordinate();
      
      // Set the event number
      currentHit->m_EvtNum = e;

      // Get MVD Strip hit
      PndSdsHit *pMvdStripHit = (PndSdsHit*) MvdStripHitArray->At(i);
      // DEBUG INFO
      //std::cout << "Sensor strip ids " << pMvdStripHit->GetSensorID() <<'\n' ;

      //______************* Extract MC data.
      int refindex = pMvdStripHit->GetRefIndex();
      
      if( refindex != -1)// None noise hits
      {
        PndSdsMCPoint *sdsMCPoint = (PndSdsMCPoint*) MvdMCPointAr->At(refindex);
        TVector3 mvdHitPos = sdsMCPoint->GetPosition();

	// Set MC trackID
	currentHit->m_trackID = sdsMCPoint->GetTrackID();
        
#if (DEBUG_PRINT > 0)
        mvdHitPos.Print();
#endif
        // MC data Polar MVD Strip
        theta_deg = Cartesian_To_Polar(mvdHitPos.X(), mvdHitPos.Y(), r_Theta);
      
        MCposPolar.Fill(r_Theta.first, r_Theta.second, theta_deg, mvdHitPos.Z(), e);
        
        // MC data Cartesian coordinates.
        MCposition.Fill(mvdHitPos.X(), mvdHitPos.Y(), mvdHitPos.Z(), e);
        mvdStrCntMC++;

        // Update current hit
        // Polar coordinates
        currentHit->mr = r_Theta.first;
        currentHit->mtheta = r_Theta.second;
        currentHit->mthetaDeg = theta_deg;
        // Cartesian
        currentHit->mx = mvdHitPos.X();
        currentHit->my = mvdHitPos.Y();
        currentHit->mz = mvdHitPos.Z();
      }
      else {// Points to exclude
        MCposition.Fill(HIT_EXCLUSION, HIT_EXCLUSION, HIT_EXCLUSION, e);
        MCposPolar.Fill(HIT_EXCLUSION, HIT_EXCLUSION, HIT_EXCLUSION, e);
        (*currentHit) = HitCoordinate(HIT_EXCLUSION);
      }
      
      //______*************  Not MC data MVD Strip.
      TVector3 tmpPos = pMvdStripHit->GetPosition();
      
      // Polar coordinates.
      theta_deg = Cartesian_To_Polar(tmpPos.X(), tmpPos.Y(), r_Theta);

      posPolar.Fill(r_Theta.first, r_Theta.second, theta_deg, e);
      
      // Cartesian 
      position.Fill(tmpPos.X(), tmpPos.Y(), tmpPos.Z(), e);

      // Update point
      currentHit->r = r_Theta.first;
      currentHit->theta = r_Theta.second;
      currentHit->thetaDeg = theta_deg;

      currentHit->x = tmpPos.X();
      currentHit->y = tmpPos.Y();
      currentHit->z = tmpPos.Z();

      currentHit->type = HitCoordinate::MVD_TYPE;
      currentHit->m_detID = pMvdStripHit->GetSensorID();
      // Add to container
      cuEvtData->push_back(currentHit);
    }// End of MVD Strip hits
#endif// MVD inclusion

    /*
     * We have seen all detectors for the current event, add event
     * data to main container.
     */
    outHit_coords->push_back(cuEvtData);
  } // End of event loop
  //////// Number of hits per detector.
  info("mvdStrip_Cnt = %d, mvdPixel_Cnt = %d, Stt_Cnt = %d", mvdStrCntMC, mvdPixCntMC, SttCntMC);
  
  // Write Output data to out file.
  if( OutFile.IsOpen()) {
    OutFile.cd();
    MCposition.Write();
    MCposPolar.Write();
    
    position.Write();
    posPolar.Write();
    
    STTposition.Write();
    STTposPolar.Write();
    
    STTMCposition.Write();
    STTMCposPolar.Write();
  }

  // ------------------------------------------------------------------------
#if (WRITE_TO_ASCII_FILE > 0)  
  // Write to an ASCII file
  //////////// Comments to write to file
  std::string format = "#######\n#\n# prefix ,m, means Monte Carlo, theta is in rad,";
  format += " thetaDeg is theta converted to deg\n";
  format += "# Values equal to -10000.0 correspond to hits with a bad MC MVD value (noise)";
  format += "\n# those values have to be excluded (not real hits, ghosts or compute errors).\n";
  format += "#\n# EvtNum x y z mx  my mz r theta thetaDeg mr mtheta mthetaDeg";
  format = format + " isochrone type trackID timeStamp detID \n#\n#######\n";
  std::string header = "EvtNum,x,y,z,mx,my,mz,r,theta,thetaDeg,mr,mtheta,mthetaDeg,isochrone,type,trackID,timeStamp,detID\n";
  //////// End of comments.
  std::ofstream OutTxtFile;
  OutTxtFile.open ("HitCoordinates.csv");
  
  if (OutTxtFile.is_open()) {
    OutTxtFile << format;
    OutTxtFile << header;

    // Write
    for(size_t h = 0; h < outHit_coords->size(); ++h) {
      std::vector<HitCoordinate*> const* dd = (*outHit_coords)[h];
      for(size_t j = 0; j < dd->size(); ++j) {
        HitCoordinate const* curht = dd->at(j);
        OutTxtFile << curht->m_EvtNum <<","
                   << curht->x <<","
                   << curht->y <<","
                   << curht->z <<","
                   << curht->mx <<","
                   << curht->my <<","
                   << curht->mz <<","
                   << curht->r <<","
                   << curht->theta <<","
                   << curht->thetaDeg <<","
                   << curht->mr <<","
                   << curht->mtheta <<","
                   << curht->mthetaDeg <<","
                   << curht->isochrone <<","
		   << curht->type <<","
                   << curht->m_trackID <<","
                   << curht->m_timeStamp <<","
		   << curht->m_detID
                   << '\n';
      }
    }
  }
  
  //OutTxtFile.close();
#endif

 #if (WRITE_TO_JSON_FILE > 0)  
  // Write to an ASCII file
  //////////// Comments to write to file
  std::string formatt = "{\n \t \"Events\": [\n";
  //////// End of comments.
  std::ofstream OutTxtFilee;
  OutTxtFilee.open ("Hits.json");
  
  if (OutTxtFilee.is_open()) {
    OutTxtFilee << formatt;
    int curevt = -1;
    int oldcurevt = curevt;
    int curid = -1;
    int oldcurid = curid;
    // Write
    for(size_t h = 0; h < outHit_coords->size(); ++h) {
      std::vector<HitCoordinate*> const* dd = (*outHit_coords)[h];
      for(size_t j = 0; j < dd->size(); ++j) {
	HitCoordinate const* curht = dd->at(j);
	
	curevt = curht->m_EvtNum;
	if (curevt != oldcurevt){
	  if(oldcurevt != -1)
	    OutTxtFilee << "\n\t\t\t\t\t\t] \n \t\t\t\t } \n \t\t\t ]\n \t\t },\n";
	  OutTxtFilee << "\t\t {\n \t\t\t \"ID\": "<<curevt<<",\n \t\t\t\"Trajectories\": [\n";
	}

	curid =  curht->m_trackID ;
	if (curid != oldcurid){
	  if (oldcurid != -1 && curevt==oldcurevt){
	    OutTxtFilee << "\n \t\t\t\t\t\t] \n \t\t\t\t},\n";
	  }
	  OutTxtFilee << "\t\t\t\t {\n \t\t\t\t\t \"ID\": "<<curid<<",\n \t\t\t\t\t \"Hits\": \n \t\t\t\t\t\t[\n";
	  OutTxtFilee << "\t\t\t\t\t\t" << curht->m_detID;

	}
	else
	  OutTxtFilee << ",\n \t\t\t\t\t\t" << curht->m_detID;

	oldcurid = curid;
	oldcurevt = curevt;
      }
      
    }
     OutTxtFilee << "\n \t\t\t\t\t\t] \n \t\t\t\t }\n \t\t\t]\n \t\t}\n \t]\n }";
  }
  
  //OutTxtFile.close();
#endif

  // Collect the map into the output parameter
  CollectSttDetecorCoords( *SttTubeArray, detNodes);

  // Note the numbers differ by one(Fortran legacy??). STT array
  // counts from 1 and detNodes from 0;
  info("SttTubeArray = %d, detnodes = %d", SttTubeArray->GetEntries(), detNodes.size());
  
  timer.Stop();
  //______________ Running Time information ____________
  
   Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  info(" Macro finished succesfully. The time for reading event data and  constructing the detector map is %lf (real) and %lf (cpu)", rtime, ctime);
  // Return
  return outHit_coords;
}
//___________________________ WriteEventPlotsToFile __________________________________
void WriteEventPlotsToFile(std::vector < std::vector<HitCoordinate*>* > const &evtData,
			   TFile &OutFile)
{
  if(evtData.size() == 0) {
    std::cerr << "<ERROR> Empty input\n";
    exit(EXIT_FAILURE);
  }
  if( !OutFile.IsOpen()) {
    std::cerr << "<Warning> OutputFile is not open and writable.\n";
    exit(EXIT_FAILURE);
  }
  if(OutFile.GetDirectory("InputEvents") ) {
    OutFile.cd("InputEvents");
  }
  else {
    OutFile.mkdir("InputEvents");
    OutFile.cd("InputEvents");
  }
  // Create plots and write to outputfile
  for(size_t e = 0; e < evtData.size(); ++e) {
    std::stringstream evN;
    std::string tupName = "Evt";
    size_t evetNumber = e;
    evN << evetNumber;
    tupName = tupName + "_" + evN.str() + "_CoordsTuple";
    std::vector<HitCoordinate*> const *currentEvt = evtData[evetNumber];
    //std::cout << tupName << "\n";
    TNtuple collection (tupName.c_str(), "Collected read data from event.", "EvtNum:x:y:z:mx:my:mz");
    for(size_t h = 0; h < currentEvt->size(); ++h) {
      HitCoordinate const *CurrentHit = currentEvt->at(h);
      collection.Fill(evetNumber, CurrentHit->x, CurrentHit->y, CurrentHit->z,
		      CurrentHit->mx, CurrentHit->my, CurrentHit->mz);
    }// END hit list loop
    collection.SetMarkerColor(2);
    collection.SetMarkerStyle(6);
    collection.Write();
  }//END Events loop
  // Go to top dir
  OutFile.cd();
}
//___________________________ END WriteEventPlotsToFile ______________________________