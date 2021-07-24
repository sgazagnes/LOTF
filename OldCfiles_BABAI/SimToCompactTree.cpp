// C and C++ standard headers


#include <iostream>
#include <vector>
#include <utility>
#include <set>
#include <string>

// PANDA and RAIR Root
#include "FairRunAna.h"
#include "FairRootManager.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"
#include "FairHit.h"
#include "FairEventHeader.h"
#include "FairRootManager.h"
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
#include "PndSttSkewedHit.h"
#include "PndTrack.h"
#include "PndGeoHandling.h"
#include "FairLink.h"

// Root headers
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TArc.h"
#include "TLine.h"
#include "TMath.h"

#include "SimToCompactTree.h"


#define INCLUDE_STT_POINTS  1
#define INCLUDE_MVD_POINTS  1

#define WRITE_TO_ASCII_FILE 0

#define NOT_AVAILABLE  -10000.0
#define HIT_EXCLUSION  -10000.0


struct HitCoordinate {
  /**
   * Type of the stored coordinates.
   */
  typedef enum Hit_Detector {
    UNKNOWN  = 0,
    STT_TYPE = 1,
    MVD_TYPE = 2,
  } Hit_DetectorType;

HitCoordinate()
  : x(0.0), y(0.0), z(0.0),
    mx(0.0), my(0.0), mz(0.0),
    r(0.0), theta(0.0), thetaDeg(0.0),
    mr(0.0), mtheta(0.0), mthetaDeg(0.0),
    m_isochrone(0.0),
    m_type(UNKNOWN),
    m_detID(-1),
    m_EvtNum(0),
    m_trackID(-1),
    m_timeStamp(0.00)
{};

HitCoordinate(float const val)
  : x(val), y(val), z(val),
    mx(val), my(val), mz(val),
    r(val), theta(val), thetaDeg(val),
    mr(val), mtheta(val), mthetaDeg(val),
    m_isochrone(val),
    m_type(UNKNOWN),
    m_detID(-1),
    m_EvtNum(0),
    m_trackID(val),
    m_timeStamp(0.00)
{};

  // Coordinates of the current hit
  float x;
  float y;
  float z;
  // MC coordinates of the current hit
  float  mx;
  float  my;
  float  mz;
  // Polar Coordinates
  float r;
  double theta;
  double thetaDeg;
  // MC Polar Coordinates
  float  mr;
  double mtheta;
  double mthetaDeg;
  // If Stt the isochrone
  double m_isochrone;
  Hit_DetectorType m_type;
  // Detector id for now just STT, for others we should avoid
  int m_detID;
  int m_EvtNum;// Event number
  // MC trackID. If = -1, it does not belong to any MC track
  int m_trackID;
  float m_timeStamp;
};

double Cartesian_To_Polar(float const x, float const y, std::pair<float,float>& polarOut, bool useSign)
{
  float radius = sqrt( x*x + y*y);
  
  float theta = 0;
  
  /* Expressed in radians. Taken into account the sign of the
     coordinates.*/
  if(useSign) {
    theta = atan2(y, x);
  }
  else {// Do not take into account the sign
    // x == 0
    if( !(x > 0.0) && !(x < 0.0) ){
      theta = (M_PI / 2.0);
    }
    else if((y/x) >= 0 ) {
      theta = atan(y/x);
    }
    else {
      theta = M_PI + (atan(y/x));
    }
  }
  
  polarOut.first  = radius;
  polarOut.second = theta;
  
  // Return theta in degrees.
  return ( (theta * 180) / M_PI);// 
}

void TreeToOut(std::string const &SimInFile,// Sim input file
				     std::string const &digiInFile,// Digi input file
				     std::string const &SimParamInfile,// Simparam inputfile
				     std::string const &OutFileName,// Outputfile name
				     std::string const &TreeName, // Output tree name
				     std::string const &NodeTreeName // Tree name to hold node name
				     )
{
  // Sim File
  TFile sf(SimInFile.c_str(),"READ");
  
  // Digi File
  TFile df(digiInFile.c_str(),"READ");
  
  // Get file of parameters
  TFile parfile(SimParamInfile.c_str(), "READ");
  
  TFile outFile(OutFileName.c_str(), "RECREATE","Outputfile containing event parameters", 9);
  
  // Simu tree
  TTree* simuTr = (TTree *) sf.Get("pndsim");
  // simuTr->Print();

  TClonesArray *SttMCPointAr = new TClonesArray("PndSttPoint");
  simuTr->SetBranchAddress("STTPoint", &SttMCPointAr);

  // MVD MC points
  TClonesArray *MvdMCPointAr = new TClonesArray("PndSdsMCPoint");
  simuTr->SetBranchAddress("MVDPoint", &MvdMCPointAr);

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

   if(!geoMan) {
    cerr << "<ERROR> Could not find valid GeoManager. Abort now!."
	 << endl;
    exit(1);
  }
  fairbasepar->GetName();
  PndGeoHandling geoH(1583944737, SimParamInfile.c_str()); //1572086365
  geoH.FillSensorMap();
  
  PndGeoSttPar *SttParameters = 0; 
  // Retrieve stt parameters
  SttParameters = (PndGeoSttPar*) parfile.Get("PndGeoSttPar");

  //  -------------------------------------------------------------
  PndSttMapCreator *mapper = 0;
  mapper = new PndSttMapCreator(SttParameters);
  
  if(mapper == 0) {
    std::cout << "Could not Create the STT mapper\n";
    exit(1);
  }

  TClonesArray *SttTubeArray = 0;
  SttTubeArray = mapper->FillTubeArray();

  
  if( SttTubeArray == 0) {
    std::cout << " Could not get  STT TubeArray\n";
    exit(2);
  } else
    std::cout << "STT TubeArray filled \n";

    

  std::cout << "\n<INFO> Total number of events in digi tree = "
            << digiTr->GetEntriesFast()
            << '\n';

  // Open Output file
  TFile OutputFile(OutFileName.c_str(),"RECREATE","File containing event parameters.", 9);

  // Parameters
  float EvtNum, isochrone,x,y,z,mx,my,mz,r,theta,thetaDeg,mr,mtheta,mthetaDeg;
  int detID,trackID, type; 
  float timeStamp;
  
  TTree *OutTree = new TTree(TreeName.c_str(),"Collected event parameter");
  OutTree->Branch("detID", &detID);
  OutTree->Branch("trackID", &trackID);
  OutTree->Branch("type", &type);
  OutTree->Branch("EvtNum", &EvtNum);
  OutTree->Branch("isochrone", &isochrone);
  OutTree->Branch("x", &x);
  OutTree->Branch("y", &y);  
  OutTree->Branch("z", &z);
  OutTree->Branch("mx", &mx);
  OutTree->Branch("my", &my);
  OutTree->Branch("mz", &mz);
  OutTree->Branch("r", &r);
  OutTree->Branch("theta", &theta);
  OutTree->Branch("thetaDeg", &thetaDeg);
  OutTree->Branch("mr", &mr);
  OutTree->Branch("mtheta", &mtheta);
  OutTree->Branch("mthetaDeg", &mthetaDeg);
  OutTree->Branch("timeStamp",&timeStamp);

  // Temporary container to hold event data.
  std::vector<HitCoordinate*> Temp_container;

  //______________________ Processing events ____________________________
  std::pair<float, float> r_Theta;
  float theta_deg;
  r_Theta.first = r_Theta.second = theta_deg = 0.00;

  std::cout << "<INFO> Processing " << digiTr->GetEntriesFast()
	    << " events.\n";
  // Event loop
  for(int e = 0; e < digiTr->GetEntriesFast(); ++e) {
    //  cout << e << endl;
    simuTr->GetEntry(e);
    digiTr->GetEntry(e);
    if( !(e % 100) ) {
      std::cout << " *";
    }

    // Fetch STT data points.
    for(int i = 0; i < SttHitArray->GetEntriesFast(); ++i) {
      //  cout << i <<" \t" << SttHitArray->GetEntriesFast() << endl;
      HitCoordinate* currentHit = new HitCoordinate();
      // Set the event number
      currentHit->m_EvtNum = e;
      // Fetch the stt hit point
      PndSttHit *pSttHit = (PndSttHit*) SttHitArray->At(i);
      //___ *********** HIT/POINT POSITION MC ***********
      int refindex = pSttHit->GetRefIndex();
      PndSttPoint *pSttMCPoint = (PndSttPoint*) SttMCPointAr->At(refindex);
      // Set the trackID for the current point
      currentHit->m_trackID = pSttMCPoint->GetTrackID();
      // MC point coordinates
      TVector3 mcposition;
      pSttMCPoint->Position(mcposition);
      // MC data Convert tot Polar coordinates
      theta_deg = Cartesian_To_Polar(mcposition.X(), mcposition.Y(), r_Theta);

      // Update current hit
      // Polar coordinates
      currentHit->mr = r_Theta.first;
      currentHit->mtheta = r_Theta.second;
      currentHit->mthetaDeg = theta_deg;
      // Cartesian
      currentHit->mx = mcposition.X();
      currentHit->my = mcposition.Y();
      currentHit->mz = mcposition.Z();
      //_______ ********* Not MC Hit coordinates ********
      // Which tube
      PndSttTube *pSttTube = (PndSttTube *) SttTubeArray->At(pSttHit->GetTubeID());
      //cout<< pSttHit->GetRefIndex()<< endl;

      if(pSttHit->GetTubeID() == 0){
	cout<< "ISSUE WITH event: " << e << "\t, hit " << i << "\t, isochrone is" << pSttHit->GetIsochrone() << "\t, tube ID is" <<  pSttHit->GetTubeID() <<endl;
	continue;
      }
      // drift radius
      Double_t Tmp_isochrone   = pSttHit->GetIsochrone();
      
      // Coordinates of the center of the tube

      // COMMENT THESE LINES TO UNDERSTAND WHAT THE ISSUE IS WITH TUBE 0
      
      TVector3 tubecenter  = pSttTube->GetPosition();
      // Polar coordinates Data.
      //cout<< pSttHit->GetRefIndex() << endl;

      theta_deg = Cartesian_To_Polar(tubecenter.X(), tubecenter.Y(), r_Theta);

      // Update point
      currentHit->r = r_Theta.first;
      currentHit->theta = r_Theta.second;
      currentHit->thetaDeg = theta_deg;

      currentHit->x = tubecenter.X();
      currentHit->y = tubecenter.Y();
      currentHit->z = tubecenter.Z();// Center of the tube

      currentHit->m_isochrone = Tmp_isochrone;
      currentHit->m_type = HitCoordinate::STT_TYPE;
      
      //currentHit->m_detID = pSttTube->GetTubeID();
      currentHit->m_detID = pSttHit->GetTubeID();
      // Add to container

      // END OF COMMENT FOR TUBE 0
      
      Temp_container.push_back(currentHit);
    }// END STT
    ////________________________________ MVD PIXELS __________________________
    //+++ Pixel MVD
    for( int i= 0; i< MvdPixelHitArray->GetEntriesFast(); i++) {
      HitCoordinate* currentHit = new HitCoordinate();
      
      // Set the event number
      currentHit->m_EvtNum = e;
      
      // Get MVD Pixel hit
      PndSdsHit *pMvdPixelHit = (PndSdsHit*) MvdPixelHitArray->At(i);
      // ___ ************* Extract MC- Data
      int refindex = pMvdPixelHit->GetRefIndex();
      
      // If equal -1, then noise hits
      if( refindex != -1) {
        PndSdsMCPoint *sdsMCPoint = (PndSdsMCPoint*) MvdMCPointAr->At(refindex);
        TVector3 mvdHitPos = sdsMCPoint->GetPosition();
	
	// Set MC trackID
	currentHit->m_trackID = sdsMCPoint->GetTrackID();

        // MC data Polar position
        theta_deg = Cartesian_To_Polar(mvdHitPos.X(), mvdHitPos.Y(), r_Theta);
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
        (*currentHit) = HitCoordinate(HIT_EXCLUSION);
      }
      //_____ ************ Non MC data MVD Pixel.
      TVector3 tmpPos = pMvdPixelHit->GetPosition();
      
      // Polar coordinates.
      theta_deg = Cartesian_To_Polar(tmpPos.X(), tmpPos.Y(), r_Theta);

      // Update point
      currentHit->r = r_Theta.first;
      currentHit->theta = r_Theta.second;
      currentHit->thetaDeg = theta_deg;

      currentHit->x = tmpPos.X();
      currentHit->y = tmpPos.Y();
      currentHit->z = tmpPos.Z();

      currentHit->m_type = HitCoordinate::MVD_TYPE;
      currentHit->m_detID = pMvdPixelHit->GetSensorID();
      // Add to container
       Temp_container.push_back(currentHit);
    }// END MVD Pixel
    ////// ++++++++++++++++++++++++++ MVD STRIP +++++++++++++++
    //+++++++ Strip MVD
    for(int i= 0; i< MvdStripHitArray->GetEntriesFast(); ++i) {
      HitCoordinate* currentHit = new HitCoordinate();
      
      // Set the event number
      currentHit->m_EvtNum = e;

      // Get MVD Strip hit
      PndSdsHit *pMvdStripHit = (PndSdsHit*) MvdStripHitArray->At(i);

      //______************* Extract MC data.
      int refindex = pMvdStripHit->GetRefIndex();
      
      if( refindex != -1)// None noise hits
      {
        PndSdsMCPoint *sdsMCPoint = (PndSdsMCPoint*) MvdMCPointAr->At(refindex);
        TVector3 mvdHitPos = sdsMCPoint->GetPosition();

	// Set MC trackID
	currentHit->m_trackID = sdsMCPoint->GetTrackID();

        // MC data Polar MVD Strip
        theta_deg = Cartesian_To_Polar(mvdHitPos.X(), mvdHitPos.Y(), r_Theta);
	
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
        (*currentHit) = HitCoordinate(HIT_EXCLUSION);
      }
      
      //______*************  Not MC data MVD Strip.
      TVector3 tmpPos = pMvdStripHit->GetPosition();
      
      // Polar coordinates.
      theta_deg = Cartesian_To_Polar(tmpPos.X(), tmpPos.Y(), r_Theta);

      // Update point
      currentHit->r = r_Theta.first;
      currentHit->theta = r_Theta.second;
      currentHit->thetaDeg = theta_deg;

      currentHit->x = tmpPos.X();
      currentHit->y = tmpPos.Y();
      currentHit->z = tmpPos.Z();

      currentHit->m_type = HitCoordinate::MVD_TYPE;
      currentHit->m_detID = pMvdStripHit->GetSensorID();
      // Add to container
      Temp_container.push_back(currentHit);
    }// END MVD Strip
  }
  
  std::cout <<'\n';
  std::cout << "<INFO> Total number of hits in container = " << Temp_container.size()
	    << '\n';
  for(size_t i = 0; i < Temp_container.size(); ++i){
    HitCoordinate* ht = Temp_container[i];
    //float type = 0;
    // STT tubes
    if( ht->m_type == HitCoordinate::STT_TYPE ) {
      type = 1;
    }// MVDs
    else if( ht->m_type == HitCoordinate::MVD_TYPE ) {
      type = 2;
    }
    else {// UNKNOWN
      type = 0;
    }
  
    EvtNum = ht->m_EvtNum;
    detID  = ht->m_detID;
    trackID = ht->m_trackID;
    isochrone = ht->m_isochrone;
    x = ht->x;
    y = ht->y;
    z = ht->z;
    mx = ht->mx;
    my = ht->my;
    mz = ht->mz;
    r = ht->r;
    theta = ht->theta;
    thetaDeg = ht->thetaDeg;
    mr = ht->mr;
    mtheta = ht->mtheta;
    mthetaDeg = ht->mthetaDeg;
    timeStamp = ht->m_timeStamp;
    // Fill output tree
    OutTree->Fill();
  }
  // Write Event data tree
  OutTree->Write();

  
  std::cout << '\n';
  //___________________________________ Geometry tree ____________
  // Now we need to collect geometry parameters.
  std::cout << "<INFO> Collecting detector geometry parameters.\n";

  // Variables
  int nodeID;
  int nodeType;
  int SectorLimit;
  unsigned int Sector;
  unsigned int Layer;
  bool LayerLimit;
  double nodeX;
  double nodeY;
  double nodeZ;
  float halfLength;
  float WdirX;
  float WdirY;
  float WdirZ;
  std::vector<int> neigborList;

   TTree *Geometry = new TTree (NodeTreeName.c_str(),"Collected geometry parameters");

  // Bind Parameters
  Geometry->Branch("nodeID", &nodeID);
  Geometry->Branch("nodeType", &nodeType);
  Geometry->Branch("SectorLimit", &SectorLimit);
  Geometry->Branch("Sector", &Sector);
  Geometry->Branch("Layer", &Layer);
  Geometry->Branch("LayerLimit", &LayerLimit);
  Geometry->Branch("nodeX", &nodeX); 
  Geometry->Branch("nodeY", &nodeY);
  Geometry->Branch("nodeZ", &nodeZ); 
  Geometry->Branch("halfLength", &halfLength);
  Geometry->Branch("WdirX", &WdirX);
  Geometry->Branch("WdirY", &WdirY);
  Geometry->Branch("WdirZ", &WdirZ);
  Geometry->Branch("neigborList", &neigborList);


  // Fill tree
  //___________=====================================================================
  // Tubes numbering starts with 1 ???
  PndSttTube *tube = 0;
  size_t sLeft, sRight, nSkewed, nPara;
  nSkewed = nPara = sLeft = sRight = 0;

  // Note tube count starts at 1.
  for(int itube = 1; itube <= SttTubeArray->GetEntries(); itube++) {
    // Current tube
    tube = (PndSttTube*) SttTubeArray->At(itube);

    nodeID        = tube->GetTubeID();// This may change later
    TVector3 nodePosition  = tube->GetPosition();
    nodeX         = nodePosition.X();
    nodeY         = nodePosition.Y();
    nodeZ         = nodePosition.Z();

    WdirX = (tube->GetWireDirection()).X();
    WdirY = (tube->GetWireDirection()).Y();
    WdirZ = (tube->GetWireDirection()).Z();
    
    halfLength    = tube->GetHalfLength();
    
    if( tube->IsParallel() ) {
      nodeType  = 1;
      nPara++;
    }
    else {
      nodeType = 2;
      nSkewed++;
    }
    Sector     = tube->GetSectorID();
    Layer      = tube->GetLayerID();
    LayerLimit = tube->IsLayerLimit();

    if( tube->IsSectorLimit() == 1 ) {
      SectorLimit = 1;
      sRight++;
    }
    else if( tube->IsSectorLimit() == -1 ) {
      SectorLimit = -1;
      sLeft++;
    }
    else {
      SectorLimit = 0;
    }
    // Collect and set neighbor detector id's.
    // Fetch index list of neighors for current detector.
    TArrayI nbr = tube->GetNeighborings();
    // Copy to vector
    for(int i = 0; i < nbr.GetSize(); ++i) {
      int neigh_index = nbr.At(i);
      PndSttTube *neighb_tube = (PndSttTube*) SttTubeArray->At(neigh_index);
      int n_Id = neighb_tube->GetTubeID();
      neigborList.push_back(n_Id);
    }
    
    // Fill tree
    Geometry->Fill();

    // Empty vector for new input.
    neigborList.clear();
  }// Tube loop
  std::cout << "<INFO> Number of sector limits = " << (sRight + sLeft)
	    << ", Right_Limit = " << sRight
	    << "  Left_Limit = "  << sLeft
	    << "\n\t<-I-> number of skewed = " << nSkewed
	    << " number of para.. = " << nPara << " sum = " << (nSkewed + nPara)
	    << "\n";
  //___________=====================================================================
  // Write geometry
  Geometry->Write();

  // Clean memory
  for(size_t i = 0; i < Temp_container.size(); ++i) {
    delete Temp_container[i];
  }
  Temp_container.clear();

  // Close file and Exit
  OutputFile.Close();
  //exit(0);
}
