#include <iostream>
#include <algorithm>

#include "path_queue.h"
#include "pathopen.h"
#include "gridNode.h"
#include "CoordGrid.h"
#include "hitcoordinate.h"
#include "trackObject.h"
#include "utilfunctions.h"
#include "SttMVDEventDataReader.h"
#include "performFilter.h"

#include "TStopwatch.h"

// int main(int argc, char** argv)
int main()
{
  std::cout << "<INFO> Running tracking\n";

  std::string infile_Events = "rootfiles/CompactEventData.root";
  std::string EvenTree = "EventData";
  std::string GeoTree  = "Geometry";
  
  float lambda = 1.0;
  float tol = 0.8;
  unsigned int gapSizeInput = 0;
  unsigned int MinResponce = 6;
  std::vector < std::vector<HitCoordinate*>* > *eventsData =
    new std::vector < std::vector<HitCoordinate*>* >();
  
  // Load geometry
  std::vector <GridNode> *GeoStt = (std::vector < GridNode >*) ReadGridGeometry( infile_Events, GeoTree);
  // Read event parameters.
  SttMVDParametersTreeFileRead(infile_Events, EvenTree, -1, -1, BY_EVENT, eventsData);
  // Create grid structure
  CoordGrid gr;
  // Init grid for the given nodes
  gr.Initialize(*GeoStt);

  // Determine virtual nodes
  std::vector < GridNode > VNodes;
  Compute_Add_VirtualNodes_Neigbor(gr, VNodes);
  //Extend the grid.
  gr.ExtendedGrid(VNodes);
  // Fix the missing symmetric neighbouring relations
  fixNeighboring(gr);

  unsigned int startEvt = 0;
  unsigned int lastEvt  = 1;//eventsData->size();

  // Process events
  TStopwatch timer;
  timer.Start();
  for(unsigned int e = startEvt; e < lastEvt; e++) {
    std::cout << "<INFO> Processing event: " << e << '\n';
    // Current event coordinates
    std::vector<HitCoordinate*> const *ev = eventsData->at(e);
    if(ev) {
      gr.FillGrid(*ev);
    }
    std::vector< std::set<int>* >* connectedComp = 0;

    // Start from most outer layer
    compOrientAttLayerBased_Local(gr, lambda, tol, gapSizeInput, true);
    // Start from most inner layer
    //compOrientAttLayerBased_Local(gr, lambda, tol, gapSizeInput, false);

    // Determine connected components
    connectedComp = AttSpaceConCompLayerBasedSoft(gr, MinResponce);

    std::cout << "<INFO> Found " << connectedComp->size()
	      << " components.\n";
    // Free memory
    for(unsigned int cmp = 0; cmp < connectedComp->size(); ++cmp) {
      delete connectedComp->at(cmp);
    }
    delete connectedComp;

    // Reset grid to process the next event.
    gr.ResetGrid();
    
  }// Event loop
  timer.Stop();
  double rtime = timer.RealTime();
  double ctime = timer.CpuTime();
  double totalNumEvent = static_cast<double>( lastEvt  - startEvt );
  std::cout <<"=======================================\n"
            << "Macro finished succesfully.\n"
            << "\tReal time " << (rtime/totalNumEvent)
	    << " (s/event), CPU time " << (ctime/totalNumEvent)
            << " (s/event)."
            << " lambda = " << lambda
            << ", tolerance = " << tol
            << '\n';
  
  // Clean up memory before returning
  for(unsigned int i = 0; i < eventsData->size(); i++) {
    std::vector<HitCoordinate*>* member = eventsData->at(i);
    for(unsigned int j = 0; j < member->size(); ++j) {
      delete member->at(j);
    }
    delete member;
  }
  delete eventsData;
  return 0;
}