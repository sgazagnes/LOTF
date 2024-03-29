/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#include <iostream>
#include <cstdlib>

#include "SttMVDEventDataReader.h"
#include "hitcoordinate.h"


#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TVector3.h"

#define PRINT_READER_DEBUG_INFO 0

std::vector <GridNode>* ReadGridGeometry(std::string const &InFileName,// Inputfile name
					 std::string const &GeoTree // Input node geometry
					 )
{
  TFile inFile(InFileName.c_str(),"READ");
  TTree *Geometry = (TTree*) inFile.Get(GeoTree.c_str());
  if(!Geometry){
    std::cerr << "<ERROR> Could not find geometry tree.\n";
    exit(EXIT_FAILURE);    
  }
  std::cout << "<INFO> Loading geometry from " << InFileName
	    << '\n';
  std::vector <GridNode> *outPut = new std::vector <GridNode>();
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
  std::vector<int> *neigborList = 0;


  // Bind Parameters
  Geometry->SetBranchAddress("nodeID", &nodeID);
  Geometry->SetBranchAddress("nodeType", &nodeType);
  Geometry->SetBranchAddress("SectorLimit", &SectorLimit);
  Geometry->SetBranchAddress("Sector", &Sector);
  Geometry->SetBranchAddress("Layer", &Layer);
  Geometry->SetBranchAddress("LayerLimit", &LayerLimit);
  Geometry->SetBranchAddress("nodeX", &nodeX);
  Geometry->SetBranchAddress("nodeY", &nodeY);
  Geometry->SetBranchAddress("nodeZ", &nodeZ);
  Geometry->SetBranchAddress("halfLength", &halfLength);
  Geometry->SetBranchAddress("WdirX", &WdirX);
  Geometry->SetBranchAddress("WdirY", &WdirY);
  Geometry->SetBranchAddress("WdirZ", &WdirZ);
  Geometry->SetBranchAddress("neigborList", &neigborList);

  std::cout << " Geometry->GetEntries() = " << Geometry->GetEntries() << std::endl;

  for(unsigned int i = 0; i < Geometry->GetEntriesFast(); ++i) {
    // Set tree level
    Geometry->GetEntry(i);

    GridNode node;// = detNodes[i];
    node.m_detID = nodeID;
    node.m_Orig_detID = nodeID;
#if( PRINT_READER_DEBUG_INFO > 1)
    std::cout << "Reading node " << i << " id = " << node.m_detID
	      << " Orig = " << node.m_Orig_detID << " ";
#endif
    node.m_x = nodeX;
    node.m_y = nodeY;
    node.m_z = nodeZ;
    
    node.m_halfLength = halfLength;
    (node.m_WireDirection).SetX(WdirX);
    (node.m_WireDirection).SetY(WdirY);
    (node.m_WireDirection).SetZ(WdirZ);

    node.m_Sector = Sector;
    node.m_Layer  = Layer;
    node.m_LayerLimit = LayerLimit;
    node.m_SectorLimit = SectorLimit;
    
    // Add neigbors
#if( PRINT_READER_DEBUG_INFO > 1 )
    std::cout << " numNei = " << neigborList->size() << " ";
#endif
    for(unsigned int n = 0; n < neigborList->size(); ++n) {
      int neib = neigborList->at(n);
      node.m_neighbors.push_back(neib);
    }
    // Set type and set weight accordingly
    if( nodeType == 1 ){// STT parallel (Real)
      node.m_type = GridNode::STT_TYPE_PARA;
      node.m_weight = 1;
    }
    else if( nodeType == 2 ) {// Skewed (real)
      node.m_type = GridNode::STT_TYPE_SKEW;
      node.m_weight = 1;
    }
    else if( nodeType == 6 ) {// Virtual node
      node.m_type = GridNode::VIRTUAL_NODE;
      node.m_weight = 0;
    }
     
    // HIER Ben je bezig FIXME
    node.ComputeSlope();
    // add to container
    outPut->push_back(node);
#if( PRINT_READER_DEBUG_INFO > 1 )
    std::cout << " Added\n";
#endif
  }
  std::cout << "<INFO> Read geometry from file: " << InFileName
	    << ". Added " << outPut->size()
	    << " nodes to grid.\n";
  // Close file
  inFile.Close();
  return outPut;
}

void SttMVDParametersTreeFileRead(std::string const &InFileName,// input file name
				  std::string const &TreeName, // Input tree name
				  int firstEvt, int lastEvt,
				  GroupingOrder grouping,
				  // Output parameter.
				  std::vector < std::vector<HitCoordinate*>* >* DataContainer
				  )
{
  TFile inFile(InFileName.c_str(),"READ");
  TTree *inTree = (TTree*) inFile.Get(TreeName.c_str());
  
  // Clean the container to avoid mixing
  if( !DataContainer) {
    std::cerr << "<ERROR> Output container is not initialized properly.\n";
    exit(EXIT_FAILURE);
  }
  
  if( DataContainer->size() > 0 ){
    std::cerr << "<INFO> Container seems to have " << DataContainer->size()
	      << " Members. We will delete all elements.\n";
    for(size_t i = 0; i < DataContainer->size(); ++i) {
      std::vector < HitCoordinate* > *hitList = DataContainer->at(i);
      for(size_t j = 0; j < hitList->size(); ++j) {
	delete hitList->at(j);
      }
      delete DataContainer->at(i);
    }
    DataContainer->clear();
  }

  //local temporary container befor sorting events.
  std::vector <HitCoordinate*> temporaryHitContainer;
  
  // Variables to read
  float EvtNum, isochrone,x,y,z,mx,my,mz,r,theta,thetaDeg,mr,mtheta,mthetaDeg;
  int detID,trackID, type;
  float timeStamp;

  inTree->SetBranchAddress("detID", &detID);
  inTree->SetBranchAddress("trackID", &trackID);
  inTree->SetBranchAddress("type", &type);
  inTree->SetBranchAddress("EvtNum", &EvtNum);
  inTree->SetBranchAddress("isochrone", &isochrone);
  inTree->SetBranchAddress("x", &x);
  inTree->SetBranchAddress("y", &y);  
  inTree->SetBranchAddress("z", &z);
  inTree->SetBranchAddress("mx", &mx);
  inTree->SetBranchAddress("my", &my);
  inTree->SetBranchAddress("mz", &mz);
  inTree->SetBranchAddress("r", &r);
  inTree->SetBranchAddress("theta", &theta);
  inTree->SetBranchAddress("thetaDeg", &thetaDeg);
  inTree->SetBranchAddress("mr", &mr);
  inTree->SetBranchAddress("mtheta", &mtheta);
  inTree->SetBranchAddress("mthetaDeg", &mthetaDeg);
  inTree->SetBranchAddress("timeStamp", &timeStamp);

  // Number of available events.
  int startEvent, stopEvent;
  
  if( (firstEvt < 0) && (lastEvt < 0) ) {// Process whole tree.
    startEvent = 0;
    stopEvent = inTree->GetEntries();
  }
  else {// Process selected events.
    startEvent = firstEvt;
    stopEvent  = lastEvt;
  }
  //  for(int e = startEvent; e < stopEvent; ++e) {
  for(int e = startEvent; e < inTree->GetEntries(); ++e) {
    // Set tree level
    inTree->GetEntry(e);
    // Create object
    HitCoordinate *current_hit = new HitCoordinate();
    current_hit->x = x;
    current_hit->y = y;
    current_hit->z = z;
    current_hit->mx = mx;
    current_hit->my = my;
    current_hit->mz = mz;
    current_hit->r = r;
    current_hit->theta = theta;
    current_hit->thetaDeg = thetaDeg;
    current_hit->mr = mr;
    current_hit->mtheta = mtheta;
    current_hit->mthetaDeg = mthetaDeg;
    current_hit->isochrone = isochrone;
    if(type == 1) {
      current_hit->type = HitCoordinate::STT_TYPE;
    }
    else if(type == 2) {
      current_hit->type = HitCoordinate::MVD_TYPE;
    }
    else{
      current_hit->type = HitCoordinate::UNKNOWN;
    }
    current_hit->m_detID = detID;
    current_hit->m_EvtNum = EvtNum;
    current_hit->m_trackID = trackID;
    current_hit->m_timeStamp = timeStamp;
    // Add to temporary container.
    temporaryHitContainer.push_back(current_hit);
  }
  // Now we need to group the hits either event based or time based.
  switch(grouping){
  case BY_TIMESTAMP:
    GroupByTimeStamp(startEvent, stopEvent, 0.00, 1.00, temporaryHitContainer, *DataContainer);
    break;
  default:
    GroupeByEvent(startEvent, stopEvent, temporaryHitContainer, *DataContainer);
    break;
  }
  // Clean temporary allocated container
  for(size_t i = 0; i < temporaryHitContainer.size(); ++i) {
    delete temporaryHitContainer[i];
  }
  temporaryHitContainer.clear();

  inFile.Close();
}

void GroupeByEvent( int firstEvt, int lastEvt,
		    std::vector < HitCoordinate* > const &HitContainer,
		    std::vector < std::vector<HitCoordinate*>* > &Output)
{
  std::cout << "<INFO> Grouping events \"event based\".\n"
	    << "      <-I-> input contains " << HitContainer.size()
	    << " hits.\n";
  // Events loop
  for(int e = firstEvt; e < lastEvt; ++e) {
    // Allocate memory for hits of current event
    std::vector< HitCoordinate* > *eventHits = new std::vector< HitCoordinate* >();
    // Sellect all hits belonging to the current event
    for(size_t i = 0; i < HitContainer.size(); ++i){
      HitCoordinate const *currentHit =  HitContainer[i];
      if( currentHit->m_EvtNum == e ){
	eventHits->push_back( new HitCoordinate( (*currentHit) ) );
      }
    }// END HIT loop
    if( eventHits->size() > 0){
      Output.push_back(eventHits);
    }
    else{
      delete eventHits;
    }
  }//END event loop
  std::cout << "<INFO> Output contains " << Output.size() << " events\n";
}

void GroupByTimeStamp( int firstEvt, int lastEvt, float tZero, float offset,
		       std::vector <HitCoordinate*> const &HitContainer,
		       std::vector < std::vector<HitCoordinate*>* > &Output)
{
  std::cout << "<WARNING> Time based is not implemented yet.\n There are "
	    << " Start = " << firstEvt << " last = " << lastEvt
	    << HitContainer.size() << " Input hits and "
	    << Output.size() << " output frames. tzero = "
	    << tZero << " offset = " << offset << '\n';
}
