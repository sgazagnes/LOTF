/**************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#include <iostream>
#include <algorithm>
#include <set>
#include <list>
#include <limits>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cassert>
#include <cmath>
//#include <ctgmath>
#include <utility>

// LOCAL headers
#include "pathopen.h"
#include "CoordGrid.h"
#include "path_queue.h"
#include "gridNode.h"
#include "utilfunctions.h"
#include "trackObject.h"
#include "pathCandidate.h"

// ROOT
#include "TPrincipal.h"
#include "TMatrixT.h"
#include "TVectorT.h"
#include "TMatrixTBase.h"

/// Local tollerance. This needs to be defined global???
#define SKEWED_VIRTUAL_ANGLE_TOLLERANCE 1.06 // ~60 degrees

// If we can use c++11 with Root, then use almost equal from utility
// function header
#define Local_ConComp_Epsilon 0.0001

/// DEBUG AND DEFINES
#define INCLUDE_LEFTLEFT_RIGHTRIGHT 0 // Hyper connectivity in orient space.
#define PATH_DEBUG_PRINT 1
#define ORIENT_DEBUG_INFO 1
#define MVD_MERGE_DEBUG_PRINT 1
#define Z_DETERMINATION_DEBUG 0
#define CONNECTED_COMPONENT_DEBUG 1
#define ADD_SHORTQUEUE_TO_TAIL 1
#define REORDER_AND_NORM_NODES_FOR_Z 1
#define LAYERBASED_DEBUG_PRINTS 1
//0: Outer->inner, 1: Inner->outer
#define LAYERBASED_OUTER_TO_INNER 0
// if (val > 0), then do not merge
#define DO_NOT_MERGE_TRACKLETS 0

/**
 * Performs path opening on a give detector grid.
 *@param hitMap The detector (STT) map containing the detector hits.
 *@param length The minimum path length required. The shorter paths
 * are removed from the map.
 */
// FIXME FIXME. The computed path is not correct. We need to fix
// this. The reason is maybe in counting forward and backward
// pixels. We need to memorize the largest path so far per taken step.
void pathopen( CoordGrid &hitMap, size_t pl)
{
  std::cout << "<INFO> Performing compelete path opening.\n";

  // Init Queue
  PathQueue forward;
  PathQueue backward;

  /*
   * Insert all the active detectors in queue. Both forward and
   * backward. We loop once forward (from inner circle to outside) and
   * once backward (from outer circle toward the center.
   */
  // std::vector< GridNode > &Ingrid = hitMap.GetGrid();
  std::vector< GridNode > &Ingrid = hitMap.m_grid;
  
  /*
   * Sort the nodes increasing order. From inside to outside (forward
   * direction). The nodes are sorted based on the detector
   * ID. Detector ID = 1 is the first one in the most inner circle and
   * the IDs increase circular toward the outer circle
   * (counterclockwise).
   */
  std::sort(Ingrid.begin(), Ingrid.end(), LessThanID);
  // Decreasing
  // std::sort(Ingrid.begin(), Ingrid.end(), GrThanLayer);
  // Increasing
  std::sort(Ingrid.begin(), Ingrid.end(), LessThanLayer);

  // Insert all active detectors (pixels) in both queues
  for(size_t i = 0; i < Ingrid.size(); ++i) {
    GridNode const &node = Ingrid[i];
    // Insert if active
    if ( node.m_active ) {
      forward.inQueue(node.m_detID);
      backward.inQueue(node.m_detID);
    }
  }

#if (PATH_DEBUG_PRINT > 1)
  std::cout << "\t<-I-> Total number of elements: Forward = "
            << forward.getNumElement()
            << ", Backward = "
            << backward.getNumElement()
            << '\n';
#endif

  // First walk in the forward direction. Note the ids are increasing
  // and circular.
  for(size_t i = 0; i < forward.getNumElement(); i++) {
    // Detector ID  of the first active element in the forward queue.
    int CurrDetID = forward.GetElement(i);//popFront();
    // Find detector index in grid.
    int d_Index = hitMap.Find(CurrDetID);

#if ( PATH_DEBUG_PRINT > 0 )
    assert(d_Index != -1);
#endif

    // Fetch the detctor ref from the grid to modify
    GridNode &node = Ingrid[d_Index];
    // Fetch the neighbor list for the current detector.
    std::vector<int> const& neighList = node.GetNeighbors();
    //update current node.
    size_t lambda = 0;
    // Neighbor loop
    for(size_t j = 0; j < neighList.size(); j++) {
      // Index of the neighbor in the map (Grid).
      size_t neigh_Index = hitMap.Find(neighList[j]);
      // Neighbor node
      GridNode &neigh_node = Ingrid[neigh_Index];
      if( neigh_node.m_active ) {
        if( neigh_node.m_lengthFW > lambda ) {
          lambda = neigh_node.m_lengthFW;
        }
      }
      // Lambda is maximum
      node.m_forwardVisited = true;
      node.m_lengthFW = node.m_weight + lambda;
    }// End of neighbor loop
  }// End of forward queue handeling

  /*
   * Process Backward direction (reverse the sorted queue of detector
   * id's).
   */
  backward.Reverse();

  for(size_t j = 0; j < backward.getNumElement(); ++j) {
    // Detector ID  of the first active element in the forward queue.
    int CurrDetID = backward.GetElement(j);//popFront();
    // Find detector index in grid.
    int d_Index = hitMap.Find(CurrDetID);

#if ( PATH_DEBUG_PRINT > 0 )
    assert(d_Index != -1);
#endif

    // Fetch the detctor ref from the grid to modify
    GridNode &node = Ingrid[d_Index];
    // Fetch the neighbor list for the current detector.
    std::vector<int> const& neighList = node.GetNeighbors();
    // update current node.
    size_t lambda = 0;
    // Neighbor loop
    for(size_t k = 0; k < neighList.size(); k++) {
      // Index of the neighbor in the map.
      size_t neigh_Index = hitMap.Find(neighList[k]);
      // Neighbor node
      GridNode &neigh_node = Ingrid[neigh_Index];
      if( neigh_node.m_active ) {
        if( neigh_node.m_lengthBW > lambda) {
          lambda = neigh_node.m_lengthBW;
        }
      }
      // Lambda is maximum
      node.m_backwardVisited = true;
      node.m_lengthBW = node.m_weight + lambda;
    }// End of neighbor loop
  }// End of backward direction

  // Update the whole map
  for(size_t i = 0; i < Ingrid.size(); ++i) {
    GridNode &Det_Node = Ingrid[i];
    if(Det_Node.m_active) {
      Det_Node.m_length = (Det_Node.m_lengthFW + Det_Node.m_lengthBW) - 1;
    }
    // Deactivate paths shorter than "pl".
    if( Det_Node.m_length < pl) {
      Det_Node.m_active = false;
    }
  }
}
//________________ End Path Opening ________________________
//________________ Area Opening _____________________________
/**
 * Performs area opening on a give detector grid.
 *@param hitMap The detector (STT) map containing the detector hits.
 *@param area The minimum required area (count of pixels). The smaller
 * areas are removed from the map.
 */
std::vector< std::set<int>* >* areaOpen( CoordGrid &hitMap, size_t area)
{
  std::cout << "<INFO> Performing Area Opening with area = " << area
            << '\n';

  // Init Queue. Start from outer tubes
  PathQueue areaQue;
  
  std::vector< std::set<int>* >* Output = new std::vector<std::set<int>* >();
  
  /*
   * Insert all the active detectors in queue. Process from outside
   * toward the center of the detector.
   */
  std::vector< GridNode > &Ingrid = hitMap.m_grid;
  /*
   * Sort the nodes decreasing order. From outside to inside (backward
   * direction).
   */
  // std::sort(Ingrid.begin(), Ingrid.end(), GreaterThanID);
  // std::sort(Ingrid.begin(), Ingrid.end(), LessThanID);
  // Decreasing
  // std::sort(Ingrid.begin(), Ingrid.end(), GrThanLayer);
  // Increasing
  std::sort(Ingrid.begin(), Ingrid.end(), LessThanLayer);

  for(size_t i = 0; i < Ingrid.size(); ++i) {
    GridNode const &node = Ingrid[i];
    // Insert if active
    if ( node.m_active ) {
      areaQue.inQueue(node.m_detID);
    }
  }// All actief detectors are in queue.
  std::cout << "\t<-I-> Total number of elements: in queue = "
            << areaQue.getNumElement()
            << '\n';

  // Process Queue elements.  
  while ( !areaQue.isEmpty() ) {
    // Detector ID  of the first active element in the queue.
    int CurrDetID = areaQue.popFront();

    // Find detector index in grid.
    int d_Index = hitMap.Find(CurrDetID);

#if ( PATH_DEBUG_PRINT > 0 )
    assert(d_Index != -1);
#endif

    // Fetch the detctor ref from the grid to modify
    GridNode &node = Ingrid[d_Index];
    
    // If not visited already, then process as new starting point
    if( !node.m_visited ) {
      node.m_visited = true;
      // Seed for new area (connected comp).
      std::set<int>* currDetIDSet = new std::set<int>();
      currDetIDSet->insert(node.m_detID);
      // Add to the list of connected compontents (area's)
      Output->push_back(currDetIDSet);
    }// End if( !node.m_visited )

    /* Else: node was visited previously, We need to find out in witch
       list (connected component) it was added before. */
    std::set<int>::iterator it;
    std::vector <size_t> componentList;

    for(size_t k = 0; k < Output->size(); ++k) {
      // Current component
      std::set<int> const *curset = Output->at(k);
      it = curset->find(node.m_detID);
      if( it != curset->end() ) {
        componentList.push_back(k);
      }
    }
    
    // Fetch the neighbor list for the current detector.
    std::vector<int> const& neighList = node.GetNeighbors();

    // Process All neigbors.
    for(size_t k = 0; k < neighList.size(); k++) {
      // Index of the neighbor in the map.
      size_t neigh_Index = hitMap.Find(neighList[k]);
      GridNode &neigh_node = Ingrid[neigh_Index];
      // Check if actief.
      if( neigh_node.m_active ) {
        // Set visited
        neigh_node.m_visited = true;
        // Add to the set of participating points in the current
        // connected component.
        for(size_t c = 0; c < componentList.size(); ++c) {
          size_t indexToAdd = componentList[c];
          std::set<int>* connectedComp = Output->at(indexToAdd);
          connectedComp->insert(neigh_node.m_detID);
        }
      }
    }// End of neigbor processing
  }// End of queueprocessing(WHILE)
#if ( PATH_DEBUG_PRINT > 1 )
  std::cout << "\n\t<DEBUG_INFO> Number of connected components = "
            << Output->size()
            << "  With the following Areas:\n";
  for(size_t l = 0; l < Output->size(); ++l) {
    std::cout << (Output->at(l))->size() << '\t';
  }
  std::cout << '\n';
#endif
  // Update detector map and Output before return, conform the required
  // size "area"
  for(size_t i = 0; i < Output->size(); ++i) {
    // Current connected component Detector IDs
    std::set<int> const* curSet = Output->at(i);
    // Update area in all detectors in the current set.
    std::set<int>::iterator it;
    for( it = curSet->begin(); it != curSet->end(); ++it) {
      int detID = *it;// Id of the current detector
      int d_Index = hitMap.Find(detID);// Index in the grid
      GridNode &node = Ingrid[d_Index];
      // Set the area in which the current detector is participating
      node.m_area = curSet->size();
      // Trun off detectors participating in this connected component
      // if the area smaller than required.
      if( area > curSet->size()) {
        node.m_active = false;
      }
    }
  }
  return Output;
}
//________________ End Area Opening _____________________________

//________________ PathOpenTransform ____________________________
void PathOpenTransform(CoordGrid &hitMap, GridNode &inputnode,
		       float orint, float tol,
		       std::pair<size_t, float>& Leng_Radius,
		       std::vector<int> &ListOfTubes,
		       size_t gapSize)
{
  if(ListOfTubes.size() != 0) {
    ListOfTubes.clear();
  }

  // Dummy to hold the results of slope determination
  std::pair<float,float> Rad_deg;
  
  // Lenght of the longest line (in this case the area)
  size_t currentLength = 1;
  float max_radius = std::numeric_limits<size_t>::min();
  float radius = 0.0;

  // Detector (pixel) queue.
  PathQueue detectorQueue;
  
  // List of all nodes available in the image (detector map)
  std::vector< GridNode > &Ingrid = hitMap.m_grid;
  
  // Insert input node in queue ( first pixel to process)
  detectorQueue.inQueue(inputnode.m_detID);

  // Add to the output list. Length 1 is true for all orientations for
  // all tubes.
  ListOfTubes.push_back(inputnode.m_detID);

  // Process queue
  while (!detectorQueue.isEmpty() ) {
    int CurrDetID = detectorQueue.popFront();
    // Find detector index in grid.
    int d_Index = hitMap.Find(CurrDetID);
    // Fetch the detctor ref from the grid.
    GridNode &currentNode = Ingrid[d_Index];
    // Mark as seen
    currentNode.m_orintVisited = true;

    /*
     * Fetch the List of neigbors (make a local copy). (Slow but for
     * now good enough)
     */
    // std::vector<int> neighList = currentNode.GetNeighbors();
    std::vector<int> neighList(currentNode.GetNeighbors());
    //____________ ******************************************* _______/////
    // FIXME FIXME HIER BEN JE NOGBEZIG (nog niet urgent)
    // Select all active neigbors that have not been visited yet.
    std::set<int> activeNeigIdSet;
    
    for(size_t k = 0; k < neighList.size(); ++k) {
      // Index of neighbor in the map (Grid).
      size_t neigh_Index = hitMap.Find(neighList[k]);
      // Neighbor node
      GridNode &neigh_node = Ingrid[neigh_Index];
      // If active and it was not seen before. Otherwise we might walk
      // backwards.
      if( (neigh_node.m_active) &&
	  (!neigh_node.m_orintVisited) ){
	activeNeigIdSet.insert(neighList[k]);
      }
    }
    // Allow jumps of gap size
    std::set<int>::iterator it;
    if( gapSize > 0 &&
	(activeNeigIdSet.size() == 0) ) {// No active neighbors.
      // It is possible that there are no active neighbors. So the
      // list is empty. We need to copy all available neighbors to the
      // set. So we add all neighbors for the current orientation as
      // the first layer to skip.
      for(size_t i = 0; i < neighList.size(); ++i) {
	int IdxInMap = hitMap.Find(neighList[i]);
	GridNode &Nd = Ingrid[IdxInMap];
	ComputeSlope(currentNode.m_x, currentNode.m_y, Nd.m_x, Nd.m_y, Rad_deg);
	if( fabs(Rad_deg.first - orint) <=  tol ) {
	  activeNeigIdSet.insert(neighList[i]);
	}
      }
      size_t gap = 0;
      while( gap < gapSize) {
	std::set<int> AllNeighbors;
	for(it = activeNeigIdSet.begin(); it != activeNeigIdSet.end(); ++it) {
	  // Node ID
	  int nodeID = *it;
	  // Index of the neighbor in the map (Grid).
	  int NodeIdxInMap= hitMap.Find(nodeID);
	  // Neighbor node
	  GridNode &GrNd = Ingrid[NodeIdxInMap];
	  // If it is active we need to add it to the output list
	  if( GrNd.m_active && (!GrNd.m_orintVisited) ){
	    AllNeighbors.insert(nodeID);
	  }
	  std::vector<int> &NList = GrNd.GetNeighbors();
	  // Add neighbors to the list
	  for(size_t nn = 0; nn < NList.size(); ++nn) {
	    int idx = hitMap.Find(NList[nn]);
	    // Neighbor node
	    GridNode &NeiNd = Ingrid[idx];
	    ComputeSlope(currentNode.m_x, currentNode.m_y, NeiNd.m_x, NeiNd.m_y, Rad_deg);
	    if( fabs(Rad_deg.first - orint) <=  tol ) {
	      AllNeighbors.insert(NList[nn]);
	    }
	  }//FOR nn
	}//FOR n
	activeNeigIdSet = AllNeighbors;
	gap++;
      }//END WHILE
    }// END IF gapsize > 0 The list of all neighbors for a give
     // orientation is collected. We need to compile the list of
     // active nodes (tubes)
    neighList.clear();// Clear local copy
    for(it = activeNeigIdSet.begin(); it != activeNeigIdSet.end(); ++it) {
      int InMapIndex= hitMap.Find((*it));
      GridNode const &Nd = Ingrid[InMapIndex];
      if( Nd.m_active && (!Nd.m_orintVisited) ) {
     	neighList.push_back((*it));
      }
    }
    //____________ ******************************************* _______/////
    // Neigbor loop (only active)
    for(size_t i = 0; i < neighList.size(); ++i) {
      // Index of the neighbor in the map (Grid).
      size_t neigh_Index = hitMap.Find(neighList[i]);
      // Neighbor node
      GridNode &neigh_node = Ingrid[neigh_Index];
      // If active and it was not seen before. Otherwise we might walk
      // backwards.
      if( (neigh_node.m_active) &&
	  (!neigh_node.m_orintVisited) && ( neigh_node.m_type != GridNode::VIRTUAL_NODE)
	  ){
	// Compute angle with the x-axis
	ComputeSlope(currentNode.m_x, currentNode.m_y,
		     neigh_node.m_x, neigh_node.m_y,
		     Rad_deg);
	// If the slope is equal to the given orientation. (FIXME)
	// Attention: this could go wrong if there are more than one
	// neigbours with the same angle+/- tol. Needs to be
	// investigated!!??
	if( fabs(Rad_deg.first - orint) <=  tol ) {
	  // Increase length.
	  currentLength++;
	  // Add to the list of participants.
	  ListOfTubes.push_back(neigh_node.m_detID);
	  // Insert neighbor in queue
	  detectorQueue.inQueue(neigh_node.m_detID);
	  // determine the current radius;
	  radius = ( (inputnode.m_x - neigh_node.m_x) * (inputnode.m_x - neigh_node.m_x) ) + 
	           ( (inputnode.m_y - neigh_node.m_y) * (inputnode.m_y - neigh_node.m_y) );
	  radius = sqrt(radius);
	  //maximum radius
	  if( radius > max_radius) {
	    max_radius = radius;
	  }
	}// If inside tollerance window
      }//IF Neighbor active && Not visited && not virtual
      else if( (neigh_node.m_active) &&
	       (!neigh_node.m_orintVisited) &&
	       ( neigh_node.m_type == GridNode::VIRTUAL_NODE)
	       ){
	// Add to the list of participants.
	ListOfTubes.push_back(neigh_node.m_detID);
	neigh_node.m_orintVisited = true;
      }// If virtual just add the node and mark as seen
    }// END neigbor loop
  }// End while queue
  
  // We need to reset nodes.
  for(size_t j = 0; j < Ingrid.size(); ++j) {
    GridNode &node = Ingrid[j];
    node.m_orintVisited = false;
  }
  // Write to output parameter
  Leng_Radius.first = currentLength;
  Leng_Radius.second = max_radius;
}
//________________ END PathOpenTransform ________________________
//________________ DetermineAngleNthNeighDepth With pre defined number of neigbor layers ________
float DetermineAngleNthNeighDepth(CoordGrid &hitMap, GridNode &InPutNode,
				  size_t depthOflevel, std::vector<int> &NeighborsTubesList,
				  bool usePCA)
{
  if( NeighborsTubesList.size() != 0 ) {
    NeighborsTubesList.clear();
  }
  // List of all tubes in the N- level neigbors. Used for moving PCA
  // computing.
  std::vector<int> NthLevelNeighList;

  // Detector (pixel) queue.
  PathQueue TubeQueue;
  TubeQueue.inQueue(InPutNode.m_detID);

  // Add the currentNode to the list
  NthLevelNeighList.push_back(InPutNode.m_detID);

  // List of all nodes available in the image (detector map)
  std::vector< GridNode > &Ingrid = hitMap.m_grid;
  
  size_t level = 0;
  
  // Compile a list of the nth neigbours of the input node.
  while(level < depthOflevel) {
    // Temporary queue
    PathQueue neighQueue;
    while(!TubeQueue.isEmpty()) {
      // Fetch the first in the row
      int DetID = TubeQueue.popFront();
      // Find detector index in grid.
      int d_Idx = hitMap.Find(DetID);
      // Fetch the tube
      GridNode &current_Node = Ingrid[d_Idx];
      current_Node.m_orintVisited = true;
      // List of neigbors
      std::vector<int> const& neighList = current_Node.GetNeighbors();
      // Neigbor loop
      for(size_t i = 0; i < neighList.size(); ++i) {
	// Index of the neighbor in the map (Grid).
	int neigh_Index = hitMap.Find(neighList[i]);
	// Neighbor node
	GridNode const &neigh_node = Ingrid[neigh_Index];
	// Insert into the queue if active
	if ( (!neigh_node.m_orintVisited) &&
	     neigh_node.m_active) {
	  // Insert neighbor in queue
          neighQueue.inQueue(neigh_node.m_detID);
	  // Insert in the total list of tubes.
	  NthLevelNeighList.push_back(neigh_node.m_detID);
	}// END if active
      }// END Neighbor loop
    }// Queue processing
    TubeQueue = neighQueue;
    level++;
  }// WHILE Depth

  // Reset visited tube fields.
  for(size_t t = 0; t < Ingrid.size(); ++t) {
    (Ingrid[t]).m_orintVisited = false;
  }
  
  // Dummy to hold the results of slope determination
  std::pair<float,float> Rad_deg;
  float LocalMeanAngle = 0;
  
  // TubeQueue Contains the nth-depth neigbor list
  NeighborsTubesList = TubeQueue.GetListOfElements();

  // Number of neigbors.(FIXME this is not correct)
  size_t numberOfTubes = NeighborsTubesList.size();
  if( !usePCA ) {
    for(size_t i = 0; i < numberOfTubes; ++i) {
      size_t L_N_Idx = NeighborsTubesList[i];
      GridNode const &Cur_Node = Ingrid[hitMap.Find(L_N_Idx)];
      ComputeSlope(InPutNode.m_x, InPutNode.m_y,
		   Cur_Node.m_x,  Cur_Node.m_y,
		   Rad_deg);
      // temporaryResults.push_back(Rad_deg.first);
      LocalMeanAngle += Rad_deg.first; 
    }// END neigbor loop
    // Determine mean of all angles. (Dit is niet helemaal correct je
    // moet er wat beters voor verzinnen.
    LocalMeanAngle = (LocalMeanAngle / static_cast<float>(numberOfTubes));
    LocalMeanAngle = (std::isnan(LocalMeanAngle) ? 0 : LocalMeanAngle);
  }// Not PCA
  else{// USE PCA
    // Principal componen (empty options to avoid normalization)
    TPrincipal principale(2,"");
    
    // for(size_t i = 0; i < NeighborsTubesList.size(); ++i) {
    // size_t L_N_Idx = NeighborsTubesList[i];
    for(size_t i = 0; i < NthLevelNeighList.size(); ++i) {
      size_t L_N_Idx = NthLevelNeighList[i];
      GridNode const &this_Node = Ingrid[hitMap.Find(L_N_Idx)];
      double *point = new double[2];
      point[0] = this_Node.m_x;
      point[1] = this_Node.m_y;
      principale.AddRow(point);
    }// END neigbor loop
    
    // Compute PCA
    principale.MakePrincipals();
    TMatrixD const *eiVec = principale.GetEigenVectors();
    
    // std::cout << "GetNcols = " << eiVec->GetNcols() << " GetNrows =
    // 	      " << eiVec->GetNrows() << std::endl;
    
    // eiVec->Print();
    double *data = new double[2];
    
    // First Eigen Vector
    eiVec->ExtractRow(0, 0, data);
    
    // Second Eigen Vector
    // eiVec->ExtractRow(1, 0, data);
    // std::cout << data[0] << " , " << data[1] << std::endl;
    // eiVec->Print();

    if( (!std::isnan(data[0])) && (!std::isnan(data[1])) ) {
      float theta = 0;
      if( !(data[0] > 0) && !(data[0] < 0) ) {// data[0] == 0
	theta = (M_PI / 2.0);
      }
      else {
	theta = atan2(data[1], data[0]);
	if(theta < 0){
	  theta += M_PI;
	}
	if( fabs(theta - M_PI) < LOCAL_PND_TRACKING_EPSILON ) {// Theta == PI
	  theta = 0.0;
	}
      }// ELSE
      LocalMeanAngle = theta;
    }
    // Second Eigen Vector
    // eiVec->ExtractRow(1, 0, data);
    // std::cout << data[0] << " , " << data[1] << std::endl;
    delete data;
    NeighborsTubesList = NthLevelNeighList;
  }// USE PCA
  return LocalMeanAngle;
}
  // TMatrixD const *cov   = principale.GetCovarianceMatrix();
  // TVectorD const *eiVal = principale.GetEigenValues();
  //principale.Print("MSEV");
  //principale.Print("V");
  // for(size_t j = 0; j < temporaryResults.size(); ++j) {
  //   std::cout << temporaryResults[j] / tol << " ";
  // }
  // std::cout << std::endl;
//______________ END DetermineAngleNthNeighDepth With pre defined number of neigbor layers ______
//_________________ Dynamic (SOFT) orientation attribute space _______
size_t computeSoftOrienAttributes(CoordGrid &hitMap, float lambda, float tol, size_t gapSize,
				  std::vector<float> const *OrientationList)
{
  // List of orientation for this snapshot.
  std::vector<float> OrientValList;

  // Dummy vars to hold the temporary results
  std::pair<float,float> Rad_deg;
  std::pair<size_t, float> Leng_Radius;
  float curr_orint = 0.00;

  // Container to store possible orientations for the current
  // snapshot.
  std::set <float> orientations;
  std::set <float> orientationsII;
  
  // Set tollerance for the current grid
  hitMap.m_AttSpaceTolerance = tol;
  // Fetch the list of all detectors.
  std::vector< GridNode > &Ingrid = hitMap.m_grid;

  // Sort the nodes decreasing order. From outside to inside
  // (backward direction).
  // std::sort(Ingrid.begin(), Ingrid.end(), GreaterThanID);
  // std::reverse(Ingrid.begin(), Ingrid.end());
  
  // Decreasing
  std::sort(Ingrid.begin(), Ingrid.end(), GrThanLayer);
  
  // Increasing
  // std::sort(Ingrid.begin(), Ingrid.end(), LessThanLayer);
  
  // std::stable_sort (begin, end, comp)
  // std::sort(Ingrid.begin(), Ingrid.end(), GrThanLayerSec);

  // Determine orientations
  if( !OrientationList) { 
  std::cout << "<INFO> Determining the number of orientations for the current snapshot.\n";
  // Tube loop
  for(size_t i = 0; i < Ingrid.size(); ++i) {
    GridNode &currentNode = Ingrid[i];
    // Active?
    if(currentNode.m_active) {
      std::vector<int> ListOfTubes;
      // List of neighbors
      std::vector<int> const& neighList = currentNode.GetNeighbors();
      // Fetch neigbors and determine the local orientation
      for( size_t j = 0; j < neighList.size(); ++j) {
	size_t neigh_Index = hitMap.Find(neighList[j]);
	// Neighbor node
	GridNode &neigh_node = Ingrid[neigh_Index];
	if( neigh_node.m_active ) {
	  // Compute angle with the x-axis
	  ComputeSlope(currentNode.m_x, currentNode.m_y,
		       neigh_node.m_x, neigh_node.m_y,
		       Rad_deg);
	  curr_orint = Rad_deg.first;
	  PathOpenTransform(hitMap, currentNode, curr_orint, tol, Leng_Radius, ListOfTubes);
	  if( (Leng_Radius.first > 1) ) {
	    // orientations.insert(round(curr_orint/tol));
	    orientations.insert(floor(curr_orint/tol));
	  }
	}// IF neighbor active
      }// Neighbor loop
      // Second method. Use PCA or jump a number of levels. PCA is
      // moving window. Jump method needs to be adapted.
      std::vector<int> LstOfNeighTbs;
      double angle_local = DetermineAngleNthNeighDepth(hitMap, currentNode, 2, LstOfNeighTbs, true);
      orientationsII.insert( floor(angle_local/tol) );
    }// If current node
  }// Tube loop

  // PRINT DEBUG INFO
#if (ORIENT_DEBUG_INFO > 0)
  std::cout << "<DEBUG> Computing attribute space for " << orientations.size()
	    << " orientations."
	    << "\n\t<Second method> = " << orientationsII.size()
	    << '\n';
#endif
  // We have found all available orientations for the current
  // snapshot.
  //std::copy(orientations.begin(), orientations.end(), std::back_inserter(OrientValList));
  std::copy(orientationsII.begin(), orientationsII.end(), std::back_inserter(OrientValList));
  }//if( !OrientationList)
  //____ Use pre defined orientation list
  else {
    std::cout << "<INFO> Processing current snapshot with pre defined orientations.\n";
    OrientValList = std::vector<float> (*OrientationList);
  }
  // Compute orientation space for all nodes.
  std::cout << "<INFO> Computing attribute space for " << OrientValList.size()
	    << " orientations. With gapSize = " << gapSize << '\n';
  
  // PRINT SOME DEBUG INFO
#if (ORIENT_DEBUG_INFO > 0)
  std::cout << "\t List of orientations: ";
  for(size_t r = 0; r < OrientValList.size(); ++r) {
    std::cout << (OrientValList[r] * tol) << ": ";
  }
  std::cout << '\n';
#endif
  ////////// __________________ END DEBUG INFO
  // Init node orientations for all tubes in the grid. __________________
  for(size_t gn = 0; gn < Ingrid.size(); ++gn){
    GridNode &Gr_node = Ingrid[gn];
    Gr_node.initNodeOrientation(OrientValList.size());
  }

  // We need to trim the orientations (ToDo).
  for(size_t k = 0; k < OrientValList.size(); ++k) {
    curr_orint = (OrientValList[k]) * tol;
    // Node loop
    for(size_t i = 0; i < Ingrid.size(); ++i){
      // Fetch the detector ref. from the grid to modify
      GridNode &curr_node = Ingrid[i];
      // If node is active
      if( curr_node.m_active && curr_node.m_type != GridNode::VIRTUAL_NODE
	  ) {
	std::vector<int> ListOfNeighborTubes;
        // Determine the maximum length of the line starting from
        // current node.
	PathOpenTransform(hitMap, curr_node, curr_orint, tol, Leng_Radius, ListOfNeighborTubes, gapSize);

	// Update current tube
	if( ListOfNeighborTubes.size() > (curr_node.m_orientations)[k].m_memberIds.size() ) {
	  (curr_node.m_orientations)[k].m_angle     = curr_orint;
	  (curr_node.m_orientations)[k].m_radius    = Leng_Radius.second;
	  (curr_node.m_orientations)[k].m_memberIds = ListOfNeighborTubes;
	}
	// Update neigbors as well
	//========== Maybe remove
	for(size_t nb = 0; nb < ListOfNeighborTubes.size(); ++nb) {
	  int nbIx = ListOfNeighborTubes[nb];
	  int inGrid_Index = hitMap.Find(nbIx);
	  GridNode &Neigbor_Node = Ingrid[inGrid_Index];	  
	  if( ListOfNeighborTubes.size() > (((Neigbor_Node.m_orientations)[k]).m_memberIds).size()){
	    (Neigbor_Node.m_orientations)[k].m_angle     = curr_orint;
	    (Neigbor_Node.m_orientations)[k].m_radius    = Leng_Radius.second;
	    (Neigbor_Node.m_orientations)[k].m_memberIds = std::vector<int>(ListOfNeighborTubes);
	  }
	}
	//======= Maybe remove
      }// IF active
    }// END Node loop
  }// Angle loop
  
  // All the attribute values are computes and the participating nodes
  // are updated.  Visit all nodes and determine the max and min
  // attribute values.
  size_t indexMin, indexMax;
  indexMin = indexMax = 0 ;
  size_t localMin = std::numeric_limits<size_t>::max();
  size_t localMax = 0;
  
  for(size_t i = 0; i < Ingrid.size(); ++i) {
    // Fetch the detector ref. from the grid to modify
    GridNode &curr_node = Ingrid[i];
    // If node is active
    if(curr_node.m_active) {
      std::vector<NodeOrientation> const &orientLists = curr_node.m_orientations;
      for(size_t j = 0; j < orientLists.size(); ++j) {
	NodeOrientation const &NodeOrient = orientLists[j];
	// Update maximum
	if( (NodeOrient.m_memberIds).size() > localMax) {
	  localMax = (NodeOrient.m_memberIds).size();
	  indexMax = j;
	}
	// Update minimum
	if( (NodeOrient.m_memberIds).size() < localMin) {
	  localMin = (NodeOrient.m_memberIds).size();
	  indexMin = j;
	}
      }// END Orientation loop
      // Update current node
      curr_node.m_minOrientVal   = localMin;
      curr_node.m_minOrientIndex = indexMin;
      curr_node.m_maxOrientVal   = localMax;
      curr_node.m_maxOrientIndex = indexMax;
    }// END If active
     // Reset
    localMin = std::numeric_limits<size_t>::max();
    localMax = 0;
  }// End of grid loop
  // Assign each node to corresponding orientation space. Note: One
  // node might be assigned to one or more different orientations.
  //___________
  // Current value of f(x,a)
  size_t AttrValCurrOrient = 0;
  for(size_t i = 0; i < Ingrid.size(); ++i) {
    GridNode &curr_node = Ingrid[i];
    // If node is active
    if(curr_node.m_active) {
      std::vector< NodeOrientation > &AllOrientationList = curr_node.m_orientations;
      for(size_t j = 0; j < AllOrientationList.size(); ++j) {
	AttrValCurrOrient = (AllOrientationList[j]).m_memberIds.size();

#if (PATH_DEBUG_PRINT > 1)
	std::cout << "\n<DEBUG INFO>:"
		  << " Orient = " << (AllOrientationList[j]).m_angle
		  << " Index = "  << j
		  << " Val = "    << AttrValCurrOrient
		  << " MinOrVal = " << curr_node.m_minOrientVal
		  << " MaxOrVal = " << curr_node.m_maxOrientVal
		  << " lambda = " << lambda;
#endif
	if( (AttrValCurrOrient > (lambda * curr_node.m_minOrientVal)) ||
	    (AttrValCurrOrient == curr_node.m_maxOrientVal)
	    ){
	  // Node is active in the current orientation space
	  (AllOrientationList[j]).m_act = true;
	}
      }// Orientation list loop
    }// If active
  }// Tubes loop
  // Return the total number of different angles that have been
  // selected.
  return orientations.size();
}
//_________________ END Dynamic orientation attribute space __________
///+++++++++++++++++ EDITING HERE  +++++++++++++++++++++++++++++++++++++
//__________________ Process Skewed in XY-Plane _________________________
int determineSkewed_XYPlane( CoordGrid &hitMap, GridNode const &VNode,
                            std::vector<int> &ListOfSkewedNodesIndex,
                            std::vector<int> &ListOfVirtualNodesIndex,
                            bool OuterToInner)
{
  std::cout << "<INFO> Correcting xy-coordinates of skewed nodes for this instance"
	    << "\n\tStrarting from V_Node " << VNode.m_detID
	    << ", OuterToInner = "
	    << (OuterToInner ? "True.\n" : "False.\n");
  // Fetch all graph nodes
  std::vector< GridNode > &Ingrid = hitMap.m_grid;
  // Neighbours of the input node
  std::vector<int> const &ListOfNeighbors = VNode.m_neighbors;
  // Local variables
  int LastVirtualNodeIdx = -10;
  std::pair<float, float> Rad_deg;
  std::vector<int>::iterator FindIt;
  PathQueue SkewedNodeIndexQueue;
  // There are exactly 2 neighbours (different slopes). And there are
  // two types of virtual nodes (Ax_Vi_Sk or Sk_Vi_S).
  assert( ListOfNeighbors.size() == 2);
  int FNode_id    = ListOfNeighbors[0];//First ID
  int FNode_index = hitMap.Find(FNode_id);//First index
  GridNode &First_Neigh = Ingrid[FNode_index];
  int SNode_id    = ListOfNeighbors[1];//Second ID
  int SNode_index = hitMap.Find(SNode_id);//Second Index
  GridNode &Second_Neigh = Ingrid[SNode_index];
  // All neighbours must be active
  assert( First_Neigh.m_active && Second_Neigh.m_active );
  size_t LocalCurrentLayer, LocalNextLayer;
  LocalCurrentLayer = LocalNextLayer = 0;
  //Axial->Virtual->Skewed
  if( First_Neigh.m_type  == GridNode::STT_TYPE_SKEW &&
      Second_Neigh.m_type != GridNode::STT_TYPE_SKEW ) {
    if( OuterToInner && ( First_Neigh.m_Layer < Second_Neigh.m_Layer ) ) {
      SkewedNodeIndexQueue.inQueue(FNode_index);
      LocalCurrentLayer = First_Neigh.m_Layer;
    }
    else if( !OuterToInner && ( First_Neigh.m_Layer > Second_Neigh.m_Layer ) ) {
      SkewedNodeIndexQueue.inQueue(FNode_index);
      LocalCurrentLayer = First_Neigh.m_Layer;
    }
  }
  else if( First_Neigh.m_type != GridNode::STT_TYPE_SKEW &&
	   Second_Neigh.m_type == GridNode::STT_TYPE_SKEW) {
    if( OuterToInner && ( Second_Neigh.m_Layer < First_Neigh.m_Layer) ){
      SkewedNodeIndexQueue.inQueue(SNode_index);
      LocalCurrentLayer = Second_Neigh.m_Layer;
    }
    else if( !OuterToInner && ( Second_Neigh.m_Layer > First_Neigh.m_Layer) ){
      SkewedNodeIndexQueue.inQueue(SNode_index);
      LocalCurrentLayer = Second_Neigh.m_Layer;
    }
  }
  // SKewed->Virtual->SKewed, Both neighbours are skewed
  else {
    if(OuterToInner){//Add smaller layer to the queue
      if( First_Neigh.m_Layer < Second_Neigh.m_Layer) {
	SkewedNodeIndexQueue.inQueue(FNode_index);
        LocalCurrentLayer = First_Neigh.m_Layer;
      }
      else {
	SkewedNodeIndexQueue.inQueue(SNode_index);
        LocalCurrentLayer = Second_Neigh.m_Layer;
      }
    }//END Outer to inner
    else {// Inner to outer
      if( First_Neigh.m_Layer > Second_Neigh.m_Layer) {
	SkewedNodeIndexQueue.inQueue(FNode_index);
        LocalCurrentLayer = First_Neigh.m_Layer;
      }
      else {
	SkewedNodeIndexQueue.inQueue(SNode_index);
        LocalCurrentLayer = Second_Neigh.m_Layer;
      }
    }// Inner to outer
  }
  //assert(LocalCurrentLayer > 0);
  if(OuterToInner){
    if( LocalCurrentLayer > 0) {
      LocalNextLayer = LocalCurrentLayer - 1;
    }
    else{
      LocalNextLayer = 0;
    }
  }
  else{// if( !OuterToInner) {// Inner to outer
    LocalNextLayer = LocalCurrentLayer + 1;
  }
  // Process skewed node index queue
  while( !SkewedNodeIndexQueue.isEmpty() ) {
    int sk_idx = SkewedNodeIndexQueue.popFront();
    GridNode &SK_Node = Ingrid[sk_idx];
    FindIt = std::find(ListOfSkewedNodesIndex.begin(), ListOfSkewedNodesIndex.end(), sk_idx);
    if( FindIt == ListOfSkewedNodesIndex.end() ) {
      ListOfSkewedNodesIndex.push_back(sk_idx);
    }
    // List of neighbours
    std::vector<int> const &Neighbours = SK_Node.m_neighbors;
    for( size_t l = 0; l < Neighbours.size(); ++l) {
      int sknID  = Neighbours[l];
      int sknIdx = hitMap.Find(sknID);
      GridNode &SK_NeighNode = Ingrid[sknIdx];
      // Active and skewed
      if( (SK_NeighNode.m_active) && (SK_NeighNode.m_type == GridNode::STT_TYPE_SKEW) ) {
        // Current or next layer
        if( (SK_NeighNode.m_Layer == LocalCurrentLayer) ||
            (SK_NeighNode.m_Layer == LocalNextLayer) ) {
          FindIt = std::find(ListOfSkewedNodesIndex.begin(), ListOfSkewedNodesIndex.end(), sknIdx);
          // Was not added before
          if( FindIt == ListOfSkewedNodesIndex.end() ) {
            SkewedNodeIndexQueue.inQueue(sknIdx);
            ListOfSkewedNodesIndex.push_back(sknIdx);
          }// End If not added before
        }// Current or next layer
      }// END if active and skewed.
      // If neighbour node is a virtual node. We can start
      // correcting and computing the orientations.
      else if( (SK_NeighNode.m_active) && (SK_NeighNode.m_type == GridNode::VIRTUAL_NODE) ) {
	FindIt = std::find(ListOfVirtualNodesIndex.begin(), ListOfVirtualNodesIndex.end(), sknIdx);
	if( (FindIt == ListOfVirtualNodesIndex.end()) && (sknID != VNode.m_detID) ) {
	  // Outer -> inner layer
	  if( (OuterToInner) && (SK_NeighNode.m_Layer <= SK_Node.m_Layer)) {
	    ListOfVirtualNodesIndex.push_back(sknIdx);
	    LastVirtualNodeIdx = sknIdx;
	  }
	  // Inner -> outer layer
	  if( (!OuterToInner) && (SK_NeighNode.m_Layer >= SK_Node.m_Layer)) {
	    ListOfVirtualNodesIndex.push_back(sknIdx);
	    LastVirtualNodeIdx = sknIdx;
	  }
	}
      }// END if active and virtual
    }//END neighbour loop
  }//END WHILE (SkewedNodeIndexQueue)
  int lastVirtualID = -10;
  if (LastVirtualNodeIdx >= 0) {
    GridNode &LastVirtual_node = Ingrid[LastVirtualNodeIdx];
    lastVirtualID = LastVirtual_node.m_detID;
    /* A list of all skewed nodes between the two virtual nodes is
       created by now and we can proceed with corrections of the
       xy-coordinates. Determine the dx and dy between current node and
       last virtual node.*/
    float x_diff = fabs(VNode.m_x - LastVirtual_node.m_x);
    float y_diff = fabs(VNode.m_y - LastVirtual_node.m_y);
    x_diff /= static_cast<float>(ListOfSkewedNodesIndex.size());
    y_diff /= static_cast<float>(ListOfSkewedNodesIndex.size());
    if( VNode.m_x > LastVirtual_node.m_x ) {
      x_diff *= -1;
    }
    if( VNode.m_y > LastVirtual_node.m_y ) {
      y_diff *= -1;
    }
    ////__________________ DEBUG PRINTS ______________________________________
    std::cout << "<DEBUG>\tListOfSkewedNodesIndex.size() " << ListOfSkewedNodesIndex.size()
              << " ListOfVirtualNodesIndex " << ListOfVirtualNodesIndex.size()
              << " x_diff = " << x_diff << ", y_diff = " << y_diff
              << '\n';
    std::cout << " Starting from " << VNode.m_detID << "("<< VNode.m_Layer << "): ";
    for(size_t v = 0; v < ListOfVirtualNodesIndex.size(); ++v) {
      int indexof = ListOfVirtualNodesIndex[v];
      GridNode const &Vnd = Ingrid[indexof];
      std::cout << Vnd.m_detID << "("<< Vnd.m_Layer << "), ";
    }
    std::cout << " Last V_id = " << LastVirtual_node.m_detID << std::endl;
    ////__________________ DEBUG PRINTS ______________________________________
    float xInc = VNode.m_x + x_diff;
    float yInc = VNode.m_y + y_diff;
    /* Correct xy-coordinates of the skewed nodes */
    for(size_t m = 0; m < ListOfSkewedNodesIndex.size(); ++m) {
      GridNode &skewedToproc = Ingrid[ListOfSkewedNodesIndex[m]];
      skewedToproc.m_xDet = xInc;
      skewedToproc.m_yDet = yInc;
      xInc += x_diff;
      yInc += y_diff;
    }
  }
  /* Skewed nodes are pre-processed. RETURN*/
  return(lastVirtualID);
}
//________________ End Process Skewed in XY-Plane _______________________
void compOrientAttLayerBased_Alt(CoordGrid &hitMap, float tol, bool OuterToInner)
{
  std::cout << "<INFO> Determining Layer Based Orientation spaces(Alt).\n"
	    << "\t tol = " << tol
            << " OuterToInner = " << (OuterToInner?"true\n":"false\n");
  // Fetch all graph nodes
  std::vector< GridNode > &Ingrid = hitMap.m_grid;
  // Store all active nodes in queue
  PathQueue ActiveNodeQueue;
  unsigned int maxLayer = std::numeric_limits<unsigned int>::min();
  unsigned int minLayer = std::numeric_limits<unsigned int>::max();
  // Reset visiting-parameters and push in queue if active.
  for(unsigned int i = 0; i < Ingrid.size(); ++i) {
    Ingrid[i].m_orintVisited    = false;
    Ingrid[i].m_visited         = false;
    Ingrid[i].m_times_visited   = 0;
    Ingrid[i].m_forwardVisited  = false;
    Ingrid[i].m_backwardVisited = false;
    // Add if active.
    if( Ingrid[i].m_active ) {
      int NodeId = Ingrid[i].m_detID;
      ActiveNodeQueue.inQueue(NodeId);
      // Find max layer in the current image
      if( Ingrid[i].m_Layer > maxLayer ) {
	maxLayer = Ingrid[i].m_Layer;
      }
      // Find min layer in the current image
      if( Ingrid[i].m_Layer < minLayer) {
	minLayer = Ingrid[i].m_Layer;
      }
    }
  }
  // All active nodes are in queue; minimum and maximum layer are
  // determined. Compute step size (direction)
  unsigned int currentLayer = (OuterToInner ? maxLayer : minLayer);
  int layerStp = ( OuterToInner ? -1 : 1);
  unsigned int nextLayer = currentLayer + layerStp;
  // Variables for local (re)use
  float orientation = 0.0;
  int   indexOfOrient = -1;
  int   neighborIndexOrinet = -1;
  int   V_Node_orientIdx = -1;
  int   Sk_ND_OrientIdx  = -1;
  std::pair<float, float> RadDeg;
  // Process queue
  while( !ActiveNodeQueue.isEmpty() ) {
    // List of all elements in queue
    std::vector<int> const &InQuElements = ActiveNodeQueue.m_queueCont;
    // Active Members in the current layer (index)
    std::list <int> LayerMemberIndexList;
    // To be removed from the queue (node ID)
    std::vector<int> ToDeleteNodeID;
    // Populate the lists of nodes in current layer.
    for(unsigned int i = 0; i < InQuElements.size(); i++) {
      int nodeID = InQuElements[i];
      // Index in grid
      int nodeIndex = hitMap.Find(nodeID);
      if( Ingrid[nodeIndex].m_Layer == currentLayer ) {
        // Add to list to be processed
        LayerMemberIndexList.push_back(nodeIndex);
        // Mark to be deleted from the queue
        ToDeleteNodeID.push_back(nodeID);
      }
    }// END (Populate current layer node list)
     /* All active nodes in current layer are added to the
        list. Process list.*/
    while(!LayerMemberIndexList.empty()) {
      int curNodeIndex = LayerMemberIndexList.front();// Read front Node
      LayerMemberIndexList.pop_front();//remove node
      GridNode &currentNode = Ingrid[curNodeIndex];// Get node
      std::vector<int> const &ListOfNeighbors = currentNode.m_neighbors;//neighbour list
      /*
       * If Axial or virtual:   Compute the orientations.
       * If Skewed: Remember to add afterwards. We do not add
       * orientation objects to the skewed nodes.
       */
      if(currentNode.m_type != GridNode::STT_TYPE_SKEW ) {
        // Process neighbors in current- or next layer.
        for(size_t nb = 0; nb < ListOfNeighbors.size(); ++nb) {
          int neighbourbId    = ListOfNeighbors[nb];
          int neighbourIndex = hitMap.Find(neighbourbId);
          GridNode &neighborNode = Ingrid[neighbourIndex];
          if(neighborNode.m_active && (neighborNode.m_type != GridNode::STT_TYPE_SKEW) ){
            if( (neighborNode.m_Layer == currentLayer) || (neighborNode.m_Layer == nextLayer) ) {
              // Compute orientation Attribute
              ComputeSlope(currentNode.m_x, currentNode.m_y, neighborNode.m_x, neighborNode.m_y, RadDeg);
              orientation = RadDeg.first;// In Rad
              // If the orientation Att already exists
	      indexOfOrient       = currentNode.GetOrientationIndex(orientation, tol, OuterToInner);
	      neighborIndexOrinet = neighborNode.GetOrientationIndex(orientation, tol, OuterToInner);
              if( indexOfOrient >= 0) {// Found Att
                NodeOrientation &nodeOrinet = currentNode.GetOrientation(indexOfOrient, OuterToInner);
		(nodeOrinet.m_memberIds).push_back(neighborNode.m_detID);
              }//Found orientation current node
              else{// Not found create new
                NodeOrientation ndOrinet;
                ndOrinet.m_angle = orientation;
                (ndOrinet.m_memberIds).push_back(currentNode.m_detID);
                (ndOrinet.m_memberIds).push_back(neighborNode.m_detID);
                (ndOrinet.m_val)++;
                currentNode.InsertOrientAtt(ndOrinet, OuterToInner);
              }
            }// Cuurent or next layer
          }//Active neighbour node
        }// Node loop
      }// END current not skewed
      // END if not skewed
    }// END list of current layer members.
    // Remove processed nodes from the queue
    for(size_t d = 0; d < ToDeleteNodeID.size(); ++d) {
      ActiveNodeQueue.deQueue(ToDeleteNodeID[d]);
    }
    // Go to next layer
    currentLayer += layerStp;
    nextLayer    += layerStp;
  }//End while ! empty queue
}
//___________________ compOrientAttLayerBased_Local ______________________
void compOrientAttLayerBased_Local( CoordGrid &hitMap, float lambda, float tol,
				    size_t gapSize, bool OuterToInner)
{
  // NOTE: We do not use "gap size" Yet.
  std::cout << "<INFO> Determining Layer Based Orientation spaces(Local).\n"
	    << "\t lambda = " << lambda << " tol = " << tol
            << " gapSize = " << gapSize
            << " OuterToInner = " << (OuterToInner?"true\n":"false\n");
  // Fetch all graph nodes
  std::vector< GridNode > &Ingrid = hitMap.m_grid;
  // Store all active nodes in queue and determine min and max layer
  PathQueue ActiveNodeQueue;
  unsigned int max_layer = std::numeric_limits<unsigned int>::min();
  unsigned int min_layer = std::numeric_limits<unsigned int>::max();
  
  // Reset visiting-parameters and push in queue if active.
  for(unsigned int n = 0; n < Ingrid.size(); ++n) {
    Ingrid[n].m_orintVisited    = false;
    Ingrid[n].m_visited         = false;
    Ingrid[n].m_times_visited   = 0;
    Ingrid[n].m_forwardVisited  = false;
    Ingrid[n].m_backwardVisited = false;
    // Add if active.
    if( Ingrid[n].m_active ) {
      int NodeId = Ingrid[n].m_detID;
      ActiveNodeQueue.inQueue(NodeId);
      // Find max layer in the current image
      if( Ingrid[n].m_Layer > max_layer ) {
	max_layer = Ingrid[n].m_Layer;
      }
      // Find min layer in the current image
      if( Ingrid[n].m_Layer < min_layer) {
	min_layer = Ingrid[n].m_Layer;
      }
    }
  }
  // All active nodes are in queue; minimum and maximum layer are
  // determined. Compute step size (direction)
  unsigned int current_layer = (OuterToInner ? max_layer : min_layer);
  int layerStep = ( OuterToInner ? -1 : 1);
  unsigned int nextLayer = current_layer + layerStep;
  // Variables for local (re)use
  float orientation = 0.0;
  int   indexOfOrient = -1;
  int   neighborIndexOrinet = -1;
  int   V_Node_orientIdx = -1;
  int   Sk_ND_OrientIdx  = -1;
  std::pair<float, float> Rad_deg;
  
  // Print debug info
#if(LAYERBASED_DEBUG_PRINTS > 0)
  std::cerr << "<DEBUG> Added " << ActiveNodeQueue.getNumElement()
	    << " in queue. Max layer = " << max_layer << ", Min layer = " << min_layer
	    << ", Current layer = " << current_layer  << ", Stepsize = " << layerStep
	    << '\n';
#endif
  // Process queue
  while( !ActiveNodeQueue.isEmpty() ) {
    // List of all elements in queue
    std::vector<int> const &InQuElements = ActiveNodeQueue.m_queueCont;
    // Active Members in the current layer (indices)
    std::list <int> LayerMemberIndexList;
    // To be removed from the queue (node IDs)
    std::vector<int> ToDeleteNodeID;
    // Populate the lists of nodes in current layer.
    for(unsigned int i = 0; i < InQuElements.size(); i++) {
      int nodeID = InQuElements[i];
      // Index in grid
      int nodeIndex = hitMap.Find(nodeID);
      if( Ingrid[nodeIndex].m_Layer == current_layer ) {
	// Add to list to be processed
	LayerMemberIndexList.push_back(nodeIndex);
	// Mark to be deleted from the queue
	ToDeleteNodeID.push_back(nodeID);
      }
    }// END (Populate current layer node list)
    /* All active nodes in current layer are added to the list.*/
    while(!LayerMemberIndexList.empty()) {
      int cur_node_Index = LayerMemberIndexList.front();// Read front node
      LayerMemberIndexList.pop_front();// remove from list
      GridNode &current_Node = Ingrid[cur_node_Index];// Get node
      std::vector<int> const &ListOfNeighbors = current_Node.m_neighbors;//neighbours
      /*
       * If axial:   Compute the orientation.
       * If virtual: Coming to a virtual from axial: Process normally.
       *    Note: Virtuals may not start a new path.
       * If Skewed:  Correct xy-using virtuals.
       */
      if( current_Node.m_type == GridNode::STT_TYPE_PARA ) {
	// Process neighbors in current- or next layer if not skewed.
	for(size_t nb = 0; nb < ListOfNeighbors.size(); ++nb) {
	  int nb_id    = ListOfNeighbors[nb];
	  int nb_index = hitMap.Find(nb_id);
	  GridNode &neighborNd = Ingrid[nb_index];
	  // If neighbour active & not skewed.
	  if( neighborNd.m_active && (neighborNd.m_type != GridNode::STT_TYPE_SKEW) ) {
	    // Current or next layer(avoid walking back).
	    if( (neighborNd.m_Layer == current_layer) || (neighborNd.m_Layer == nextLayer) ) {
	      // Compute orientation attribute.
	      ComputeSlope(current_Node.m_x, current_Node.m_y, neighborNd.m_x, neighborNd.m_y, Rad_deg);
	      orientation = Rad_deg.first;// In Rad
	      // If the orientation Attribute already exists.
	      indexOfOrient       = current_Node.GetOrientationIndex(orientation, tol, OuterToInner);
	      neighborIndexOrinet = neighborNd.GetOrientationIndex(orientation, tol, OuterToInner);
	      if( indexOfOrient >= 0) {// Orientation exists -> Update
                NodeOrientation &nodeOrinet = current_Node.GetOrientation(indexOfOrient, OuterToInner);
		(nodeOrinet.m_memberIds).push_back(neighborNd.m_detID);
		// Update neighbour node
		if(neighborIndexOrinet >= 0) {
                  NodeOrientation &neighbourNodeOrinet = neighborNd.GetOrientation(neighborIndexOrinet, OuterToInner);
                  (neighbourNodeOrinet.m_memberIds).push_back(current_Node.m_detID);
		  // Increment neighbour node value
		  (neighbourNodeOrinet.m_val)++;
		}
		else{// Create new for the neighbour node
		  NodeOrientation neighbourNodeOrinet(nodeOrinet);
		  (neighbourNodeOrinet.m_val)++;
                  neighborNd.InsertOrientAtt(neighbourNodeOrinet, OuterToInner);
		}
		std::cerr <<" Found, "<< orientation << " Rad., update. " << (nodeOrinet.m_memberIds).size()
			  << " Value = " << nodeOrinet.m_val << "\n";
	      }// END found(Thus update)
	      else {// Not found, create new.
		NodeOrientation ndOrinet;
		// Set Angle
		ndOrinet.m_angle = orientation;
		// Add current node to the list of participants
		(ndOrinet.m_memberIds).push_back(current_Node.m_detID);
		// Add neighbor node
		(ndOrinet.m_memberIds).push_back(neighborNd.m_detID);
                (ndOrinet.m_val)++;
		// Add to the list of current node
                current_Node.InsertOrientAtt(ndOrinet, OuterToInner);
		// Add to the neighbour, the value is one higher
		(ndOrinet.m_val)++;
                neighborNd.InsertOrientAtt(ndOrinet, OuterToInner);
		std::cout <<"Create: ID = " << current_Node.m_detID << " neighbourID = "
			  << neighborNd.m_detID << " with " << orientation
			  <<" rad, value = " << ndOrinet.m_val <<'\n';
	      }// END create new orientation object
	    }// If in current layer or the next one
	  }// If active and not skewed
          // Reset for the next neighbour
          indexOfOrient = -1;
          neighborIndexOrinet = -1;
	}// End For Neighbors
      }// IF Para (Axial)
      /*
       * If node is a virtual node. Skip layers to find the next
       * virtual and correct skewed xy-coordinates. (it is not
       * completely correct yet, FIXME later)
       */
      else if( current_Node.m_type == GridNode::VIRTUAL_NODE ) {
        std::vector<int> SkewedNodesIdx;
        std::vector<int> VirtualNodesIndex;
        // Determine xy-values for skewed nodes between virtuals.
        int lastVNodeId = determineSkewed_XYPlane( hitMap, current_Node, SkewedNodesIdx,
                                                   VirtualNodesIndex, OuterToInner);
	ToDeleteNodeID.push_back(lastVNodeId);
        
        for(size_t vn = 0; vn < VirtualNodesIndex.size(); ++vn) {
          GridNode &VN_Node = Ingrid[VirtualNodesIndex[vn]];
          ComputeSlope( current_Node.m_x, current_Node.m_y, VN_Node.m_x, VN_Node.m_y, Rad_deg);
          orientation = Rad_deg.first;// In Rad
          // Find orient indices.
          indexOfOrient    = current_Node.GetOrientationIndex(orientation, tol, OuterToInner);
          V_Node_orientIdx = VN_Node.GetOrientationIndex(orientation, tol, OuterToInner);
          if( indexOfOrient >= 0 ) {// If exists Current node, Update
            NodeOrientation &nodeOrinet = current_Node.GetOrientation(indexOfOrient, OuterToInner);
            (nodeOrinet.m_memberIds).push_back(VN_Node.m_detID);
            if( V_Node_orientIdx >= 0 ) { // If exists virtual, update
              NodeOrientation &VNodeOrinet = VN_Node.GetOrientation(V_Node_orientIdx, OuterToInner);
              // Increment neighbour node value
              (VNodeOrinet.m_val)++;
              (VNodeOrinet.m_memberIds).push_back(current_Node.m_detID);
            }
            else{// Create new for the virtual node
              NodeOrientation VirtualNodeOrinet(nodeOrinet);
              (VirtualNodeOrinet.m_memberIds).push_back(current_Node.m_detID);
              (VirtualNodeOrinet.m_memberIds).push_back(VN_Node.m_detID);
              (VirtualNodeOrinet.m_val)++;
              VN_Node.InsertOrientAtt(VirtualNodeOrinet, OuterToInner);
            }
          }// END If orientation exists for the current node.
          else{// Not found, create new and add the nodes
            NodeOrientation ndOrinet;
            // Set Angle
            ndOrinet.m_angle = orientation;
            // Add current node to the list of participants
            (ndOrinet.m_memberIds).push_back(current_Node.m_detID);
            // Add next virtual node
            (ndOrinet.m_memberIds).push_back(VN_Node.m_detID);
            (ndOrinet.m_val)++;
            // Add to the list of current node
            current_Node.InsertOrientAtt(ndOrinet, OuterToInner);
            // Add to the neighbour, the value is one higher
            (ndOrinet.m_val)++;
            VN_Node.InsertOrientAtt(ndOrinet, OuterToInner);
          }//END Not found, create
          /* Now we need to handle skewed graph nodes. */
          indexOfOrient    = current_Node.GetOrientationIndex(orientation, tol, OuterToInner);
          V_Node_orientIdx = VN_Node.GetOrientationIndex(orientation, tol, OuterToInner);
          /* Both must be >= 0*/
          assert( (indexOfOrient >= 0) && (V_Node_orientIdx >= 0) );
          // Try to add skewed nodes to the list
          int nofitCnt = 0;
          for(size_t s = 0; s < SkewedNodesIdx.size(); ++s) {
            GridNode &Skewed_Node = Ingrid[SkewedNodesIdx[s]];
            // ComputeSlope( current_Node.m_x, current_Node.m_y,
            //               Skewed_Node.m_xDet, Skewed_Node.m_yDet, Rad_deg);
            ComputeSlope( current_Node.m_x, current_Node.m_y,
                          Skewed_Node.m_x, Skewed_Node.m_y, Rad_deg);
            float SkAbsOrient = fabs(Rad_deg.first - orientation);// In Rad
            if( (SkAbsOrient < SKEWED_VIRTUAL_ANGLE_TOLLERANCE) ) {
              NodeOrientation &Cur_nd_Orinet = current_Node.GetOrientation(indexOfOrient, OuterToInner);
              (Cur_nd_Orinet.m_memberIds).push_back(Skewed_Node.m_detID);
              (Cur_nd_Orinet.m_val)++;
              NodeOrientation &Vir_nd_Orinet = VN_Node.GetOrientation(V_Node_orientIdx, OuterToInner);
              (Vir_nd_Orinet.m_memberIds).push_back(Skewed_Node.m_detID);
              (Vir_nd_Orinet.m_val)++;
              // Add orientation space to the skewed node
              Sk_ND_OrientIdx = Skewed_Node.GetOrientationIndex(orientation, tol, OuterToInner);
              // Exists?, Update
              if(Sk_ND_OrientIdx >= 0){
                NodeOrientation &SKEW_nd_Orinet = Skewed_Node.GetOrientation(Sk_ND_OrientIdx, OuterToInner);
                (SKEW_nd_Orinet.m_memberIds).push_back(current_Node.m_detID);
                (SKEW_nd_Orinet.m_memberIds).push_back(VN_Node.m_detID);
                SKEW_nd_Orinet.m_val += 2;
              }
              else{//Create new
                /* FIXME The Skewed do not have any orientation
                   attributes. Is this correct?????
                // NodeOrientation SKEWNdNewOrinet;
                // SKEWNdNewOrinet.m_angle = orientation;
                /*
                  (SKEWNdNewOrinet.m_memberIds).push_back(current_Node.m_detID);
                  (SKEWNdNewOrinet.m_memberIds).push_back(VN_Node.m_detID);
                  SKEWNdNewOrinet.m_val += 2;
                */
                // Skewed_Node.InsertOrientAtt(SKEWNdNewOrinet, OuterToInner);
              }
              std::cout << "\t<DEBUG INFO> ID = " << Skewed_Node.m_detID << ", fits with Tollerance "
                        << SKEWED_VIRTUAL_ANGLE_TOLLERANCE << " Determined " << SkAbsOrient << '\n';
            }
            else{
              nofitCnt++;
            }
          }
          std::cout << " Totaal aantal nofit = " << nofitCnt << '\n';
          // Reset for next virtual
          indexOfOrient = -1;
          V_Node_orientIdx = -1;
        }// END Virtual node loop
        SkewedNodesIdx.clear();
        VirtualNodesIndex.clear();
      }// END if current node virtual
    }// END WHILE process nodes in current layer.
    /* Remove processed nodes from the queue*/
    for(size_t d = 0; d < ToDeleteNodeID.size(); ++d) {
      ActiveNodeQueue.deQueue(ToDeleteNodeID[d]);
    }
    // Go to next layer
    current_layer += layerStep;
    nextLayer     += layerStep;
  }// END Process queue
  // Reset node visiting. Mark orientation with max and min responce
  // and update weights
  unsigned int indexOfMax, indexOfMin;
  indexOfMax = indexOfMin = 0;
  unsigned int ValueOfMax, ValueOfMin;
  ValueOfMax = 0;
  ValueOfMin = std::numeric_limits <unsigned int>::max();
  // Visit all the nodes.
  for(size_t i = 0; i < Ingrid.size(); ++i) {
    GridNode &G_node = Ingrid[i];
    G_node.m_orintVisited   = false;
    G_node.m_forwardVisited = false;
    if(G_node.m_active) {
      std::vector< NodeOrientation > &OrientationAttList = G_node.GetOrientations(OuterToInner);
      /* Sort orientation attributes in ascending order of angle, Neighbourhood relations. */
      std::stable_sort(OrientationAttList.begin(), OrientationAttList.end(), LessThanOrient);
      for(size_t j = 0; j < OrientationAttList.size(); ++j){
        NodeOrientation &orient = OrientationAttList[j];
        std::vector<int> &oriMembers = orient.m_memberIds;
        // Remove duplicates. Slow but OK for now
        std::set<int> memberset(oriMembers.begin(), oriMembers.end());
        oriMembers.clear();
        oriMembers.insert(oriMembers.end(), memberset.begin(), memberset.end());
        assert( oriMembers.size() == memberset.size());
        orient.m_weight = oriMembers.size();
        // Find max
        if( orient.m_val > ValueOfMax ) {
          ValueOfMax = orient.m_val;
          indexOfMax = j;
        }
        // Find min
        if( orient.m_val < ValueOfMin ) {
          ValueOfMin = orient.m_val;
          indexOfMin = j;
        }
      }//END orient.Att. loop
      // Update G_node(Outer->Inner)
      if(OuterToInner) {
        G_node.m_minOrientVal   = ValueOfMin;
        G_node.m_minOrientIndex = indexOfMin;
        G_node.m_maxOrientVal   = ValueOfMax;
        G_node.m_maxOrientIndex = indexOfMax;
      }
      else {//(Inner -> Outer)
        G_node.m_maxBackwardOrientVal   = ValueOfMax;
        G_node.m_maxBackWardOrientIndex = indexOfMax;
      }
      // Rest variables
      indexOfMax = indexOfMin = 0;
      ValueOfMax = 0;
      ValueOfMin = std::numeric_limits <unsigned int>::max();
    }//END if active
  }//END Grid loop
}
//___________________ END compOrientAttLayerBased_Local __________________
//___________________     AttSpaceConCompLayerBasedSoft _________________
std::vector< std::set<int>* >* AttSpaceConCompLayerBasedSoft( CoordGrid &hitMap, size_t MinResponce)
{
  std::cout << "<INFO> Connected component analysis Per layer orderd (SOFT).\n"
	    << " MinResponce = " << MinResponce << '\n';
  /* Each (sub)path needs to have a unique id. */
  unsigned int candidateId = 0;
  
  // Fetch all graph nodes (Full graph)
  std::vector< GridNode > &Ingrid = hitMap.m_grid;
  
  // Init queues
  PathQueue ForwardQueue;// Outer->Inner
  PathQueue BackwardQueue;//Inner->Outer
  unsigned int max_layer = std::numeric_limits<unsigned int>::min();
  unsigned int min_layer = std::numeric_limits<unsigned int>::max();
  // Store all active nodes in queue
  for(unsigned int n = 0; n < Ingrid.size(); n++) {
    // Reset visiting variables.
    Ingrid[n].m_visited = false;
    Ingrid[n].m_orintVisited = false;
    Ingrid[n].m_times_visited = 0;
    // Add to queue if active.
    if( Ingrid[n].m_active ) {
      int NodeId = Ingrid[n].m_detID;
      ForwardQueue.inQueue(NodeId);
      BackwardQueue.inQueue(NodeId);
      // Find max available layer in the current image
      if( Ingrid[n].m_Layer > max_layer ) {
	max_layer = Ingrid[n].m_Layer;
      }
      // Find min layer in the current image
      if( Ingrid[n].m_Layer < min_layer) {
	min_layer = Ingrid[n].m_Layer;
      }
    }
  }
  // Outer -> Inner
  bool OuterToInner = true;
  unsigned int current_layer = max_layer;
#if(LAYERBASED_OUTER_TO_INNER == 0)
#warning "Outer to inner only."
  // Temporary list of Path candidates
  std::vector < PathCandidate* > temCandidContainerOuterToInner;
  // Process queue (outer -> Inner)
  while( !ForwardQueue.isEmpty() ) {
    // Direct access to elements in forward queue.
    std::vector<int> const &queue_Elements = ForwardQueue.m_queueCont;
    // Active Members in current layer (index)
    std::list <int> LayerActiveMemberIdxLst;
    // To be removed from the queue (node ID)
    std::vector<int> ToRemoveNodeID;
    // Select members of current layer
    for(size_t i = 0; i < queue_Elements.size(); i++) {
      int nodeID    = queue_Elements[i];
      int nodeIndex = hitMap.Find(nodeID);
      // Select only elements in current layer
      if(Ingrid[nodeIndex].m_Layer == current_layer) {
	// Add map index to list, to be processed
	LayerActiveMemberIdxLst.push_back(nodeIndex);
	// Mark to be deleted from the queue, node id's
	ToRemoveNodeID.push_back(nodeID);
      }
    }
    // Process nodes in current layer
    while( !LayerActiveMemberIdxLst.empty() ) {
      int cur_node_Index = LayerActiveMemberIdxLst.front();// Get index
      LayerActiveMemberIdxLst.pop_front();// Remove
      GridNode &CurrentNode = Ingrid[cur_node_Index];// Read node
      // The list of orientAtt (Outer-> Inner)
      std::vector< NodeOrientation > const &orientList = CurrentNode.GetOrientations(OuterToInner);
      /* If visiting for the first time we need to create a candidate
         for each orientation Attribute. Virtuals can not start new
         candidates (paths).*/
      if( (CurrentNode.m_times_visited == 0) && (CurrentNode.m_type != GridNode::VIRTUAL_NODE) ) {
        CurrentNode.m_times_visited++;// Mark as seen
        // If Isolated node (sigletone)
        if( orientList.size() == 0 ){
          PathCandidate *singleton = new PathCandidate();
          singleton->m_id = candidateId;
          candidateId++;
          singleton->m_orientation  = -1.0;
          singleton->m_OuterToInner = OuterToInner;
          // Add current node to the member list
          (singleton->m_memberIdSet)->insert(CurrentNode.m_detID);
          singleton->m_weight = 1;
          singleton->m_length = 1;
          // Add to the list of all candidates
          temCandidContainerOuterToInner.push_back(singleton);
        }
        else {// (First visit) && (not isolated) && (not Virtual)
          // Create a candidate for all orientations.
          for(size_t p = 0; p < orientList.size(); ++p) {
            NodeOrientation  const &ori = orientList[p];
            std::vector<int> const &ori_mem = ori.m_memberIds;
            PathCandidate *cand = new PathCandidate();// Create a new candidate
            cand->m_id = candidateId;// Set id
            candidateId++;// Increment global
            cand->m_orientation  = ori.m_angle;// Set the orientation
            cand->m_OuterToInner = OuterToInner;// Direction
            // Add current node to the member list
            (cand->m_memberIdSet)->insert(CurrentNode.m_detID);
            // Add members of this orientation group
            for(size_t m = 0; m < ori_mem.size(); ++m) {
              int memID  = ori_mem[m];
              int memIdx = hitMap.Find(memID);
              GridNode &member = Ingrid[memIdx];
              member.m_times_visited++;// Set node as visited
              (cand->m_memberIdSet)->insert(memID);
            }
            cand->m_weight += ori_mem.size();
            cand->m_length  = ori.m_val;
            // Add to the list of all candidates
            temCandidContainerOuterToInner.push_back(cand);
          }// END for
        }//End else
      }// End first time and not virtual
      // If node was visited before && not virtual. It has been added
      // to one of the candidates; find.
      else if( (CurrentNode.m_times_visited != 0) && (CurrentNode.m_type != GridNode::VIRTUAL_NODE) ) {
        std::vector<size_t> candidateIndices;
        std::set<int>::iterator FindIdIt;
        // Find all candidates which have the current node in their
        // members list.
        for(size_t cn = 0; cn < temCandidContainerOuterToInner.size(); ++cn){
          PathCandidate *candidate = temCandidContainerOuterToInner[cn];
          std::set<int> *candidateMemIdSet  = candidate->m_memberIdSet;
          FindIdIt = std::find(candidateMemIdSet->begin(), candidateMemIdSet->end(), CurrentNode.m_detID);
          // Found current node in current candidates members list
          if( FindIdIt != candidateMemIdSet->end() ) {
            candidateIndices.push_back(cn);
          }
        }// The list of candidates with current node is compiled.
        /* If skewed, there are no orient. Att. But it has been added
           before, thus it is a member of at least one list; so skip
           for now. */
        // Loop the orientation att. and handle.
        for(size_t i = 0; i < orientList.size(); ++i) {
          NodeOrientation  const &oriAtt = orientList[i];
          std::vector<int> const &oriAtt_mem = oriAtt.m_memberIds;
          bool makeNewCandidate = true;
          // If a candidate track exists with the same orientation
          for(size_t j = 0; j < candidateIndices.size(); ++j) {
            size_t candidate_index = candidateIndices[j];
            PathCandidate *TmpCand = temCandidContainerOuterToInner[candidate_index];
            float orientationDiff = fabs(TmpCand->m_orientation - oriAtt.m_angle);
            // Found candidate with same orientation.
            if( orientationDiff < Local_ConComp_Epsilon ) {
              makeNewCandidate = false; //Do not create new
              /* Add members to the list and mark all as visited */
              for(size_t n = 0; n < oriAtt_mem.size(); ++n) {
                int mId  = oriAtt_mem[n];
                int mIdx = hitMap.Find(mId);
                GridNode &memb = Ingrid[mIdx];
                memb.m_times_visited++;
                (TmpCand->m_memberIdSet)->insert(mId);
                // MAYBE THIS IS not correct But for now!!!!! FIX ME LATER
                TmpCand->m_length = (TmpCand->m_memberIdSet)->size();
                TmpCand->m_orientation = oriAtt.m_angle;
              }
            }// END if found
          }// END candidate index list loop
          // Found no candidate with the same orientation.
          // Create new for current orient. Att.
          if( makeNewCandidate) {
            CurrentNode.m_times_visited++;
            PathCandidate *newCand = new PathCandidate();
            newCand->m_id = candidateId;
            candidateId++;
            newCand->m_orientation  = oriAtt.m_angle;
            newCand->m_OuterToInner = OuterToInner;
            (newCand->m_memberIdSet)->insert(CurrentNode.m_detID);
            newCand->m_weight += oriAtt_mem.size();
            newCand->m_length  = oriAtt.m_val;
            // Add members and mark all as seen
            for(size_t m = 0; m < oriAtt_mem.size(); ++m) {
              int memID  = oriAtt_mem[m];
              int memIdx = hitMap.Find(memID);
              GridNode &member = Ingrid[memIdx];
              member.m_times_visited++;
              (newCand->m_memberIdSet)->insert(memID);
            }
            // Add to the list of all candidates
            temCandidContainerOuterToInner.push_back(newCand);
          }// End Add new candidate
        }// Orientation Attribute loop
      }// prev visited and non virtual
    }// END WHILE (Layer member list not empty)
    // Remove processed tubes from Queue.
    for(size_t d = 0; d < ToRemoveNodeID.size(); ++d) {
      ForwardQueue.deQueue(ToRemoveNodeID[d]);
    }
    // Go to the next layer
    current_layer--;
  }// END WHILE(Outer -> Inner)

#elif (LAYERBASED_OUTER_TO_INNER == 1)
#warning "Inner to outer only"
#endif
  /*
   * Path candidates are ready . Compile the output.
   */
  // Prepair output variable
  std::vector< std::set<int>* >* outPutVar = 0;
  outPutVar = new std::vector< std::set<int>* >();
  if(!outPutVar){
    std::cerr << "<ERROR> Could not allocate memory for connected components\n";
    exit(EXIT_FAILURE);
  }

  // Merge sub-paths If possible (dit kan beter naar perform filtering
  // verplaatst worden)
#if(LAYERBASED_OUTER_TO_INNER == 0)
  /* Mark shared singletons as invalid.*/
  cleanupDuplicateSingletons(temCandidContainerOuterToInner);
  
  /* Update heads and tails (max and min layer nodes)*/
  updateAllCandidateHeadAndTailNodes(temCandidContainerOuterToInner);
  
  /* Compute mean, coVar and invCovar for all paths */
  computeCovMatForAllCandidates(hitMap, temCandidContainerOuterToInner);
  
  /* Set id for the first merged path candidate */
  unsigned int firstOutputPathId = temCandidContainerOuterToInner.size();
  /* _________________    Try to merge ______________________ */
#if(DO_NOT_MERGE_TRACKLETS == 0)
  /* Merge on Mahalanobis distances */
  bool graphBased = true;
  bool useMahalanobis = false;
  size_t maxdepth = 3;
  // Test code
  // graphBased = false;

  // Merge tracklets, select candidates graph-paths based(O(N) methode).
  firstOutputPathId = mergeCandidatesMahalanobisDist( hitMap, temCandidContainerOuterToInner,
                                                      (CAND_MERGE_DISTANCE_SELECTION_TOLERANCE),
                                                      firstOutputPathId, graphBased, useMahalanobis, maxdepth);
  
  firstOutputPathId = mergeCandidatesMahalanobisDist( hitMap, temCandidContainerOuterToInner,
                                                      (CAND_MERGE_DISTANCE_SELECTION_TOLERANCE + 0.2),
                                                      firstOutputPathId, graphBased, useMahalanobis, maxdepth);
  
  firstOutputPathId = mergeCandidatesMahalanobisDist( hitMap, temCandidContainerOuterToInner,
                                                      (CAND_MERGE_DISTANCE_SELECTION_TOLERANCE + 0.5),
                                                      firstOutputPathId, graphBased, useMahalanobis, maxdepth);

  firstOutputPathId = mergeCandidatesMahalanobisDist( hitMap, temCandidContainerOuterToInner,
                                                      (CAND_MERGE_DISTANCE_SELECTION_TOLERANCE + 1.0),
                                                      firstOutputPathId, graphBased, useMahalanobis, maxdepth);

  firstOutputPathId = mergeCandidatesMahalanobisDist( hitMap, temCandidContainerOuterToInner,
                                                      (CAND_MERGE_DISTANCE_SELECTION_TOLERANCE + 1.5),
                                                      firstOutputPathId, graphBased, useMahalanobis, maxdepth);

  // Merge distancebased (O(n^2) methode)
  // graphBased = false;
  firstOutputPathId = mergeCandidatesMahalanobisDist( hitMap, temCandidContainerOuterToInner,
                                                      (CAND_MERGE_DISTANCE_SELECTION_TOLERANCE + 2.0),
                                                      firstOutputPathId, graphBased, useMahalanobis, maxdepth);
  
  firstOutputPathId = mergeCandidatesMahalanobisDist( hitMap, temCandidContainerOuterToInner,
                                                      (CAND_MERGE_DISTANCE_SELECTION_TOLERANCE + 2.5),
                                                      firstOutputPathId, graphBased, useMahalanobis, maxdepth);
  
  firstOutputPathId = mergeCandidatesMahalanobisDist( hitMap, temCandidContainerOuterToInner,
                                                      (CAND_MERGE_DISTANCE_SELECTION_TOLERANCE + 3.0),
                                                      firstOutputPathId, graphBased, useMahalanobis, maxdepth);

  firstOutputPathId = mergeCandidatesMahalanobisDist( hitMap, temCandidContainerOuterToInner,
                                                      (CAND_MERGE_DISTANCE_SELECTION_TOLERANCE + 3.5),
                                                      firstOutputPathId, graphBased, useMahalanobis, maxdepth);
  // Is this really needed? It is still the most effective. Consider
  // to remove after running some tests.
  graphBased = false;// Merge distancebased (O(n^2) methode)
  firstOutputPathId = mergeCandidatesMahalanobisDist( hitMap, temCandidContainerOuterToInner,
                                                      (CAND_MERGE_DISTANCE_SELECTION_TOLERANCE + 4.0),
                                                      firstOutputPathId, graphBased, useMahalanobis, maxdepth);

  // Finished distance based merging, clean overlapping
  mergeMarkOverlappingAsInvalid(temCandidContainerOuterToInner, -1);
  //+++++++++++++++++++++++++++++++++++++
  /* Merge curvature difference based. */
  // FIXME HIER BEN JE BEZIG
  firstOutputPathId = mergeCandidatesCurvatureCompat( hitMap, temCandidContainerOuterToInner,
                                                      0.1, firstOutputPathId);
  //+++++++++++++++++++++++++++++++++++++
  /* Merge on curvature compatibility properties.*/

  // Use Euclidean dist for pre-selection
  // firstOutputPathId = mergeCandidatesCurvaturePriority( hitMap, temCandidContainerOuterToInner,
  //                                                       CANDIDATE_DISTANCE_CURVATURE_MERGE,
  //                                                       firstOutputPathId, false);
  
  /* Use mahalanobis dist for pre-selection*/
  // firstOutputPathId = mergeCandidatesCurvaturePriority( hitMap, temCandidContainerOuterToInner,
  //                                                       CANDIDATE_MAHALANOBIS_DISTANCE_CURVATURE_MERGE,
  //                                                       firstOutputPathId, true);

  // std::sort(temCandidContainerOuterToInner.begin(), temCandidContainerOuterToInner.end(), GreaterThanMaxLayerP);
  
  // Merge overlapping(if one is fully contained in the other.
  //mergeMarkOverlappingAsInvalid(temCandidContainerOuterToInner, -1);
  //mergeMarkOverlappingAsInvalid(temCandidContainerOuterToInner, 2);
  //  mergeMarkOverlappingAsInvalid(temCandidContainerOuterToInner, 1);
  // Remove invalid paths(remove only for drawing, after we are done).
  removeInvalidSubPaths(temCandidContainerOuterToInner);
  //++++++++++++++++++++++++++++ End try to merge +++++++++++++++++++
#endif//(DO_NOT_MERGE_TRACKLETS)
  
  // Voor nu alleen maar van buiten naar binnen.
#elif(LAYERBASED_OUTER_TO_INNER == 1)
  std::cout << "<INFO> Merge Inner -> Outer"
            << "\t Not Implemented yet.\n";
#endif
  
  // Copy all valid paths to output
#if(LAYERBASED_OUTER_TO_INNER == 0)
  
  for(unsigned int l = 0; l < temCandidContainerOuterToInner.size(); l++){
    PathCandidate &curCand = *(temCandidContainerOuterToInner[l]);
    std::set<int> const *trk = curCand.m_memberIdSet;
    if(curCand.m_isValid) {
      std::set<int> *comp = new std::set<int>((*trk));
      outPutVar->push_back(comp);
    }
  }
  
#elif(LAYERBASED_OUTER_TO_INNER == 1)
  std::cout << "<Warning> Not implemented.\n";
#endif

  // Clean allocated temporary memory
#if(LAYERBASED_OUTER_TO_INNER == 0)
  for( size_t n = 0; n < temCandidContainerOuterToInner.size(); ++n) {
    delete temCandidContainerOuterToInner[n];
  }
  temCandidContainerOuterToInner.clear();
#elif(LAYERBASED_OUTER_TO_INNER == 1)
  std::cout << "<Warning> Not implemented.\n";
#endif
  // Not needed but for now (REMOVE later)
  for(size_t i = 0; i < Ingrid.size(); ++i) {
    Ingrid[i].m_times_visited = 0;
  }
  return outPutVar;
}
//___________________ END AttSpaceConCompLayerBasedSoft _________________
//_____________________  ConCompLayerBasedLocalOrient ___________________
/*
 * The attribute space is not pre-determined (SOFT). At each step the
 * possible candidate paths are determined and are given a
 * probability. From each starting point the highest probability might
 * be the winning one (Rough Set idea, fuzzy set).
 *
 * (NOT FINISHED YET), FIXME FIXME. HIER BEN IK MEE BEZIG.
 */
std::vector< std::set<int>* >* ConCompLayerBasedLocalOrient( CoordGrid &hitMap, float tol)
{
  std::cout << "<INFO> Connected component analysis Per layer orderd (local Orientation).\n"
	    << "\t<-I-> tolerance = " << tol
	    << '\n';

  // Fetch all graph nodes
  std::vector< GridNode > &Ingrid = hitMap.m_grid;

  // Sort the nodes decreasing layer order.
  // Decreasing
  std::sort(Ingrid.begin(), Ingrid.end(), GrThanLayer);

  // Local variables
  PathQueue ActiveQueue;// All active nodes
  unsigned int max_layer = std::numeric_limits<unsigned int>::min();
  unsigned int min_layer = std::numeric_limits<unsigned int>::max();
  
  // Temporary Path candidate list
  std::vector < PathCandidate* > temporaryComponents;
  
  // Clollect all active nodes and determine the maximum and minimum
  // layer
  for(size_t i = 0; i < Ingrid.size(); ++i) {
    // Reset visiting variables.
    Ingrid[i].m_orintVisited = false;
    Ingrid[i].m_visited = true;
    // Virtual nodes cannot start a new candidate path (connected
    // component)
    if( Ingrid[i].m_type == GridNode::VIRTUAL_NODE) {
      Ingrid[i].m_times_visited = 1;
    }
    else {
      Ingrid[i].m_times_visited = 0;
    }
    // Put in queue if active.
    if(Ingrid[i].m_active){
      ActiveQueue.inQueue(Ingrid[i].m_detID);
      // Find max available layer in the current image
      if( Ingrid[i].m_Layer > max_layer ) {
	max_layer = Ingrid[i].m_Layer;
      }
      // Find min layer in the current image
      if( Ingrid[i].m_Layer < min_layer) {
	min_layer = Ingrid[i].m_Layer;
      }
    }
  }

  // All active nodes are in queue and minimum and maximum layer
  // numbers are determined.
  unsigned int current_layer  = max_layer;
  unsigned int LookAheadLayer = 0;

  // Dummy local variables
  // std::pair< size_t, float> Leng_Radius;
  std::pair<float,float> Rad_deg;
  float orintentation = 0.00;

  // DEBUG May be deleted later
  PathCandidate *singles = new PathCandidate();

  // Process queue
  while(!ActiveQueue.isEmpty()) {
    std::vector<int> const &QuElements = ActiveQueue.m_queueCont;
    // Active Members in the current layer (index)
    std::list<int>   LayerMemberIndexList;
    // To be removed from the queue (node ID)
    std::vector<int> ToDeleteNodeID;
    // Select all active nodes in current layer
    for(size_t i = 0; i < QuElements.size(); ++i) {
      int nodeID    = QuElements[i];
      int nodeIndex = hitMap.Find(nodeID);
      if(Ingrid[nodeIndex].m_Layer == current_layer) {
	// Add to list to be processed
	LayerMemberIndexList.push_back(nodeIndex);
	// Mark to be deleted from the queue
	ToDeleteNodeID.push_back(nodeID);
      }
    }
    // Set the next layer to look at (Look ahead layer)
    LookAheadLayer = current_layer - 1;// From outer toward inner

    // Process the nodes in the current layer
    while(!LayerMemberIndexList.empty()) {
      // Read front element
      int cur_node_Index = LayerMemberIndexList.front();
      // Remove from list
      LayerMemberIndexList.pop_front();
      // Read node
      GridNode &currentNode = Ingrid[cur_node_Index];
      // Fetch the list of neighbors of the current node
      std::vector<int> &NeighbourList = currentNode.m_neighbors;
      // Not visited before, Seed for new component(s).
      if( (currentNode.m_times_visited == 0) ) {
	// Mark as visited
	currentNode.m_times_visited++;
	// For all active neigbours determine the orientation.
	for(size_t i = 0; i < NeighbourList.size(); ++i) {
	  int nId    = NeighbourList[i];
	  int nIndex = hitMap.Find(nId);
	  GridNode &NeighbourNode = Ingrid[nIndex];
	  // If neighbour node active, then start a new path
	  // candidate.
	  if(NeighbourNode.m_active &&
	     ( (NeighbourNode.m_Layer == current_layer) ||
	       (NeighbourNode.m_Layer == LookAheadLayer) ) ) {
	    // Mark as visited
	    NeighbourNode.m_times_visited++;
	    // Create new candidate
	    PathCandidate *cand = new PathCandidate();
	    // Add current node to the member list
	    (cand->m_memberIdSet)->insert(currentNode.m_detID);
	    cand->m_weight++;
	    // Add neighbour node to member list
	    (cand->m_memberIdSet)->insert(NeighbourNode.m_detID);
	    cand->m_weight++;
	    // Compute angle with the x-axis. The line between current
	    // node and the active neighbour node
	    ComputeSlope(currentNode.m_x,   currentNode.m_y,
			 NeighbourNode.m_x, NeighbourNode.m_y,
			 Rad_deg);
	    orintentation = Rad_deg.first;
	    // Update candidate
	    cand->m_orientation = orintentation;
	    // Add current candidate to the list of possible
	    // components.
	    temporaryComponents.push_back(cand);
	  }// END Neighbour active and (in the same layer or the next one)
	}// END visit all neighbours
      }// END (not visited before)
      else {// Node was visited before. The node was seen before (at
	    // least one time). Find in which candidate(s)
	std::vector<int> candidateIndex;
	std::set<int>::const_iterator FindIT;
	for(unsigned int m = 0; m < temporaryComponents.size(); ++m) {
	  PathCandidate const *cand = temporaryComponents[m];
	  std::set<int> const *MemberSet = cand->m_memberIdSet;
	  FindIT = std::find(MemberSet->begin(), MemberSet->end(), currentNode.m_detID);
	  // If the node was found
	  if( FindIT != MemberSet->end() ) {
	    candidateIndex.push_back(m);
	  }
	}// Found all candidates with the current node in their member
	 // list
	if(candidateIndex.size() == 0) {// Singleton
	  (singles->m_memberIdSet)->insert(currentNode.m_detID);
	}
#if(LAYERBASED_DEBUG_PRINTS > 1)
	std::cout << "Found node in " << candidateIndex.size()
		  << " candidates. ";
	for(unsigned int cnd = 0; cnd < candidateIndex.size(); cnd++) {	
	std:cout << candidateIndex[cnd] << ", ";
	}
	std::cout << '\n';
#endif
	// Loop through the neighbours
	for(unsigned int n = 0; n < NeighbourList.size(); ++n) {
	  int nId    = NeighbourList[n];
	  int nIndex = hitMap.Find(nId);
	  GridNode &Neighb_Node = Ingrid[nIndex];
	  // If neighbour node active, then start a new path
	  // candidate.
	  if(Neighb_Node.m_active &&
	     ( (Neighb_Node.m_Layer == current_layer) ||
	       (Neighb_Node.m_Layer == LookAheadLayer) ) ) {
	    // Mark as visited
	    Neighb_Node.m_times_visited++;

	    // Compute angle with the x-axis. The line between current
	    // node and the active neighbour node
	    ComputeSlope(currentNode.m_x,   currentNode.m_y,
			 Neighb_Node.m_x, Neighb_Node.m_y,
			 Rad_deg);
	    orintentation = Rad_deg.first;
	    // TODO: Compatibility needs to be defined (orientation
	    // based).
	    // Candidates loop (FIXME)(later)
	    for(unsigned int k = 0; k < candidateIndex.size(); ++k) {
	      unsigned int candIndex = candidateIndex[k];
	      PathCandidate *cand = temporaryComponents[candIndex];
	      (cand->m_memberIdSet)->insert(Neighb_Node.m_detID);
	      cand->m_weight++;
	    }
	  }// If active && (current or next layer)
	}// END neighbour loop
      }// End else (node was visited before
    }// END WHILE layer members
    // Remove processed nodes from the queue
    for(size_t d = 0; d < ToDeleteNodeID.size(); ++d) {
      ActiveQueue.deQueue(ToDeleteNodeID[d]);
    }
    // Go to next layer
    current_layer--;
  }// End While ! queue.empty
  /* ************************ */
  // Prepair output variable
  std::vector< std::set<int>* >* outPutVar = 0;
  outPutVar = new std::vector< std::set<int>* >();
  
  if(!outPutVar){
    std::cerr << "<ERROR> Could not allocate memory for connected components\n";
    exit(EXIT_FAILURE);
  }
  // Add singletons
  temporaryComponents.push_back(singles);

#if(LAYERBASED_DEBUG_PRINTS > 0)
  std::cout << "Nodes in singles = " << (singles->m_memberIdSet)->size()
	    << '\n';
#endif
  // Copy Member components to output sets.
  for( size_t n = 0; n < temporaryComponents.size(); ++n) {
    std::set<int> const *trk = (temporaryComponents[n])->m_memberIdSet;
    std::set<int> *comp = new std::set<int>((*trk));
    outPutVar->push_back(comp);
  }
  // Clean memory
  for( size_t n = 0; n < temporaryComponents.size(); ++n) {
    delete temporaryComponents[n];
  }
  temporaryComponents.clear();
  
  // Not needed but for now (REMOVE later)
  for(size_t i = 0; i < Ingrid.size(); ++i) {
    Ingrid[i].m_orintVisited = false;
    Ingrid[i].m_times_visited = 0;
  }
  return outPutVar;
}
//_____________________ END ConCompLayerBasedLocalOrient _______________
//_____________________  AttSpaceConCompLayerBased HARD (Static) _______
/*
 * The attribute space is pre-determined (HARD) and used during the
 * component construction operations.
 * 
 * (NOT FINISHED YET), FIXME FIXME. Deze moet ik nog afmaken.
 */
std::vector< std::set<int>* >* AttSpaceConCompLayerBasedHard( CoordGrid &hitMap, size_t MinResponce)
{
  std::cout << "<INFO> Connected component analysis Per layer orderd (HARD)."
	    << " MinPathOrint responce = " << MinResponce
	    << '\n';
 
  // Get the map nodes STT only at this moment
  std::vector< GridNode > &Ingrid = hitMap.m_grid;
  // Sort the nodes decreasing layer order.
  // Decreasing
  std::sort(Ingrid.begin(), Ingrid.end(), GrThanLayer);

  // Local variables
  PathQueue ActiveQueue;// All active nodes
  PathQueue ShortCompQueue;// Short reponce nodes
  size_t max_layer = std::numeric_limits<size_t>::min();
  size_t min_layer = std::numeric_limits<size_t>::max();

  // Path candidate list
  std::vector < PathCandidate* > temporaryComponents;

  // Clollect all active nodes and determine the maximum and minimum
  // layer
  for(size_t i = 0; i < Ingrid.size(); ++i) {
    // Reset visiting variables.
    Ingrid[i].m_orintVisited = false;
    Ingrid[i].m_visited = true;
    Ingrid[i].m_times_visited = 0;
    // Put in queue if active.
    if(Ingrid[i].m_active){
      ActiveQueue.inQueue(Ingrid[i].m_detID);
      // Find max available layer in the current image
      if( Ingrid[i].m_Layer > max_layer ) {
	max_layer = Ingrid[i].m_Layer;
      }
      // Find min layer in the current image
      if( Ingrid[i].m_Layer < min_layer) {
	min_layer = Ingrid[i].m_Layer;
      }
    }
  }
  // All active nodes are in queue and minimum and maximum layer
  // numbers are determined.
  size_t current_layer  = max_layer;
  size_t LookAheadLayer = 0;
  // Process queue
  while(!ActiveQueue.isEmpty()){
    std::vector<int> const &QuElements = ActiveQueue.m_queueCont;
    std::list<int>   LayerMemberIndexList;
    std::vector<int> ToDeleteNodeID;
    // Select all active nodes in current layer
    for(size_t i = 0; i < QuElements.size(); ++i) {
      int nodeID    = QuElements[i];
      int nodeIndex = hitMap.Find(nodeID);
      if(Ingrid[nodeIndex].m_Layer == current_layer) {
	// Add to list to be processed
	LayerMemberIndexList.push_back(nodeIndex);
	// Mark to be deleted from the queue
	ToDeleteNodeID.push_back(nodeID);
      }
    }
    // Set the next layer to look at (Look ahead layer)
    LookAheadLayer = current_layer - 1;// From outer toward inner
    
    // Process the nodes in the current layer
    while(!LayerMemberIndexList.empty()) {
      int cur_node_Index = LayerMemberIndexList.front();// Read front
      LayerMemberIndexList.pop_front();// Remove from list
      GridNode &currentNode = Ingrid[cur_node_Index];
      if(currentNode.m_times_visited == 0) {// not visited before
	// Mark as visited
	currentNode.m_times_visited++;
	// Short responce node. Delete the node from queue
	if(currentNode.m_maxOrientVal < MinResponce) {
	  // Insert in short queue. Will process later
	  ShortCompQueue.inQueue(currentNode.m_detID);
	}
	else{// Not a short responce node && it was not visited before
	  // Seed for new components
	  // List of neighbors of the current node
	  std::vector<int> &NeighbourList = currentNode.m_neighbors;
	  // List of its orientations
	  std::vector< NodeOrientation > &NodeOrientList = currentNode.m_orientations;
	  // For all orientations with responce larger than Minimum
	  // responce create a component candidate
	  for(size_t n = 0; n < NeighbourList.size(); ++n) {
	    int neighbourIndex = hitMap.Find(NeighbourList[n]);
	    GridNode  &neighbourNode = Ingrid[neighbourIndex];
	    if( neighbourNode.m_active &&
		(neighbourNode.m_Layer == current_layer ||
		 neighbourNode.m_Layer == LookAheadLayer) ) {
	      // Find the orientation space where the neighbour node
	      // is participating
	      for(size_t k = 0; k < NodeOrientList.size(); ++k) {
		NodeOrientation  &orient = NodeOrientList[k];
		std::vector<int> &Orient_members = orient.m_memberIds;
		std::vector<int>::iterator findIt = std::find( Orient_members.begin(),
							       Orient_members.end(),
							       neighbourNode.m_detID);
		if( (findIt != Orient_members.end()) &&
		    (Orient_members.size() >= MinResponce)
		   ) {
		  PathCandidate *cand = new PathCandidate();
		  cand->m_orientation = orient.m_angle;
		  // cand->m_memberList  = new std::vector<int>();
		  cand->m_memberIdSet = new std::set<int>();
		  // (cand->m_memberList)->push_back(currentNode.m_detID);
		  (cand->m_memberIdSet)->insert(currentNode.m_detID);
		  cand->m_weight++;
		  // (cand->m_memberList)->push_back(neighbourNode.m_detID);
		  (cand->m_memberIdSet)->insert(neighbourNode.m_detID);
		  cand->m_weight++;
		  neighbourNode.m_times_visited++;
		  // Add new candidate to the list of possible candidates.
		  temporaryComponents.push_back(cand);
		}
	      }// END FOR orientations list(NodeOrientList)
	    }// END IF Active && (current layer or the next one)
	  }// END FOR NeighbourList
	  // DEBUG Maybe delete later
#if(CONNECTED_COMPONENT_DEBUG > 0)
	  std::cout << "<INFO> Layer Based new seed: node_ID = " << currentNode.m_detID
		    << " Layer number = " << currentNode.m_Layer
		    << '\n';
#endif
	}// End else not short responce && not visited (seed for new
	 // component)
      }// End not visited
      else{// Node was visited before
	//_______________________________ FIXME FIXME (later)
	/*
	 * The current node has been visited before. Find to which
	 * component(s) it was added.
	 */
	std::set<int>::iterator findIterator;
	// The list of all candidates which contain the current node.
	std::vector<size_t> canditateIndexList;
	for(size_t i = 0; i < temporaryComponents.size(); ++i) {
	  PathCandidate *Current_cand  = temporaryComponents[i];
	  std::set<int>  *memberList = Current_cand->m_memberIdSet;
	  findIterator = std::find(memberList->begin(), memberList->end(), currentNode.m_detID);
	  if( findIterator != memberList->end()) {
	    canditateIndexList.push_back(i);
	  }
	}
	// Fetch the list of neighbours and orientations of the
	// current node
	std::vector<int> const &CurrNodeNeighbors = currentNode.m_neighbors;
	std::vector< NodeOrientation > const &CurrNodeOrientList = currentNode.m_orientations;
	std::vector<int>::const_iterator findIt;
	// Loop through the neighbours
	for(unsigned int  i = 0; i < CurrNodeNeighbors.size(); ++i) {
	  int neighbourID = CurrNodeNeighbors[i];
	  int neighbIndex = hitMap.Find(neighbourID);
	  // If node active and (in current layer or in the next layer)
	  if( Ingrid[neighbIndex].m_active &&
	      (
	       (Ingrid[neighbIndex].m_Layer == current_layer ) ||
	       (Ingrid[neighbIndex].m_Layer == LookAheadLayer)
	      )
	      ) {
	    std::vector<unsigned int> orientationIndices;
	    std::vector<unsigned int> Orient_candidateIndices;
	    // Find in which orientations
	    for(unsigned int  j = 0; j < CurrNodeOrientList.size(); ++j){
	      NodeOrientation const &ori = CurrNodeOrientList[j];
	      findIt = std::find( (ori.m_memberIds).begin(), (ori.m_memberIds).end(),
				  Ingrid[neighbIndex].m_detID);
	      if( (findIt != (ori.m_memberIds).end()) &&
		  ( (Ingrid[neighbIndex].m_orientations)[j].m_act )) {
		// Record the orientation index
		orientationIndices.push_back(j);
		// Find if there is already a candidate with the same
		// orientation including the current node.
		for(size_t k = 0; k < canditateIndexList.size(); ++k) {
		  PathCandidate *candidate = temporaryComponents[canditateIndexList[k]];
		  if( !(candidate->m_orientation < ori.m_angle) &&
		      !(candidate->m_orientation > ori.m_angle) ) {
		    Orient_candidateIndices.push_back(k);
		  }
		}// END candidate loop
		// whether there is a candidate with the same
		// orientation
		if( Orient_candidateIndices.size() > 0) {
		  // Add node to all candidate
		  for(unsigned int c = 0; c < Orient_candidateIndices.size(); ++c) {
		    PathCandidate *candidate = temporaryComponents[Orient_candidateIndices[c]];
		    (candidate->m_memberIdSet)->insert(Ingrid[neighbIndex].m_detID);
		    candidate->m_weight++;
		    Ingrid[neighbIndex].m_times_visited++;
		  }
		}
		else {// No candidate found. Look in the neighbourhood
		  // unsigned int nextOrientIndex = ((j < CurrNodeOrientList.size() - 1)? (j + 1) : 0);
		  // unsigned int prevOrientIndex = ((j > 0)? (j - 1) : (CurrNodeOrientList.size() - 1));
		  // float nextOrientAngle = CurrNodeOrientList[nextOrientIndex].m_angle;
		  // float prevOrientAngle = CurrNodeOrientList[prevOrientIndex].m_angle;
		  // for(size_t k = 0; k < canditateIndexList.size(); ++k) {
		  //   PathCandidate *candidate = temporaryComponents[canditateIndexList[k]];
		  //   if( !(candidate->m_orientation < nextOrientAngle) &&
		  // 	!(candidate->m_orientation > nextOrientAngle) ) {
		  //     (candidate->m_memberIdSet)->insert(Ingrid[neighbIndex].m_detID);
		  //     candidate->m_weight++;
		  //     Ingrid[neighbIndex].m_times_visited++;
		  //     candidate->m_orientation = nextOrientAngle;
		  //     std::cout << "Found next\n";
		  //   }
		  //   if( !(candidate->m_orientation < prevOrientAngle) &&
		  // 	!(candidate->m_orientation > prevOrientAngle) ) {
		  //     //(candidate->m_memberIdSet)->insert(Ingrid[neighbIndex].m_detID);
		  //     candidate->m_weight++;
		  //     Ingrid[neighbIndex].m_times_visited++;
		  //     candidate->m_orientation = prevOrientAngle;
		  //     std::cout << "Found prev\n";
		  //   }
		  // }
		  // END candidate loop
		}// ELSE
	      }// END if found in an orientation
	    }// END node orientation loop
	  }// END if active and in the correct layer
	}// END Neighbour loop

	//_______________________________ FIXME FIXME
      }// End Else visited
    }// End While (!LayerMemberIndexList.empty()
    // Remove processed nodes from the queue
    for(size_t d = 0; d < ToDeleteNodeID.size(); ++d) {
      ActiveQueue.deQueue(ToDeleteNodeID[d]);
    }
    // Go to next layer
    current_layer--;
  }//End while ( !ActiveQueue.isEmpty() )

  //==================================
  // Prepair output variable
  std::vector< std::set<int>* >* outPutVar = 0;
  outPutVar = new std::vector< std::set<int>* >();
  
  if(!outPutVar){
    std::cerr << "<ERROR> Could not allocate memory for connected components\n";
    exit(EXIT_FAILURE);
  }
  // Copy Members to sets.
  for( size_t n = 0; n < temporaryComponents.size(); ++n) {
    // std::vector<int> *trk = (temporaryComponents[n])->m_memberList;
    std::set<int> const *trk = (temporaryComponents[n])->m_memberIdSet;
    // std::set<int> *comp = new std::set<int>(trk->begin(), trk->end());
    std::set<int> *comp = new std::set<int>((*trk));
    outPutVar->push_back(comp);
  }
  // Clean memory
  for( size_t n = 0; n < temporaryComponents.size(); ++n) {
    delete temporaryComponents[n];
  }
  temporaryComponents.clear();

  // Not needed but for now(REMOVE)
  for(size_t i = 0; i < Ingrid.size(); ++i) {
    Ingrid[i].m_orintVisited = false;
    Ingrid[i].m_times_visited = 0;
  }
  // return
  return outPutVar;
}
//_____________________ END AttSpaceConCompLayerBased HARD _______________
///+++++++++++++++++ END EDITING HERE  +++++++++++++++++++++++++++++++++
/*_________________ AttSpaceConnectedComp ____________________________*/
// HINT: We can change this code to allow overlapping paths (connected
// components). Use how often visited instead of visited.
std::vector< std::set<int>* >* AttSpaceConnectedComp( CoordGrid &hitMap, size_t MinResponce)
{
  std::cout << "<INFO> Connected component analysis, Min Responce ="
	    << MinResponce
	    << '\n';

  // Tubes queue
  PathQueue ActiveQueue;
  PathQueue ShortCompQueue;
  size_t numOrientations = 0;
  
  // Local temporary variables.
  size_t LeftOrint, RightOrient, LLeftOrint, RRightOrient, maxResponsIdx;
  LeftOrint = LLeftOrint = RightOrient = RRightOrient = maxResponsIdx = 0;
  
  // Init output variable
  std::vector< std::set<int>* >* outPutVar = 0;
  outPutVar = new std::vector< std::set<int>* >();
  if(!outPutVar){
    std::cerr << "<ERROR> Could not allocate memory for connected components\n";
    assert (outPutVar != 0);
  }
  // Get the map nodes STT only at this moment
  std::vector< GridNode > &Ingrid = hitMap.m_grid;
  
  // Sort the nodes decreasing ID order.
  // Decreasing
  std::sort(Ingrid.begin(), Ingrid.end(), GrThanLayer);
  
  // Increasing
  //std::sort(Ingrid.begin(), Ingrid.end(), LessThanLayer);

  /* Insert all active tubes in tube queue. Set all tubes to NOT
   * visited.
   */
  for(size_t i = 0; i < Ingrid.size(); ++i) {
    if( (Ingrid[i].m_active) ) {
      ActiveQueue.inQueue(Ingrid[i].m_detID);
      Ingrid[i].m_orintVisited = false;
    }
  }// All active detectors are in queue.
  
  // Init the number of orientations
  if( ActiveQueue.getNumElement() > 0) {
    int id = ActiveQueue.GetElement(0);
    int Idex = hitMap.Find(id);
    numOrientations = (Ingrid[Idex].m_orientations).size();
    std::cout << "\tNumber of elements in Queue = "
	      << ActiveQueue.getNumElement()
	      << " Number of orientations = " << numOrientations
	      << '\n';
  }
  else {
    std::cerr << "<ERROR> Could not find any active tube.\n";
  }

  // Process queue for all active graph nodes.
  while( !ActiveQueue.isEmpty() ) {
    int DetID = ActiveQueue.popFront();
    int Index = hitMap.Find(DetID);
    GridNode &Cur_Node = Ingrid[Index];
    // Not visited before
    if( !(Cur_Node.m_orintVisited) ) {
      // Short responce nodes?
      if (Cur_Node.m_maxOrientVal >= MinResponce) {
	// Mark as Visited
	(Ingrid[Index]).m_orintVisited = true;

	// Seed for new component
	std::set<int>* TubesIdSet = new std::set<int>();

	// Add Id to the component list
	TubesIdSet->insert(Cur_Node.m_detID);

#if (PATH_DEBUG_PRINT > 0)
	std::cout << "<DEBUG> New seed id = " << Cur_Node.m_detID
		  << " inserted.\n";
#endif
	// Max responce orientation index
	maxResponsIdx = Cur_Node.m_maxOrientIndex;
	// Fetch the list of participants in the current orientation
	// with the current node
	std::vector<NodeOrientation> const &OrientList = Cur_Node.m_orientations;
	std::vector<int> const &membIds  = (OrientList[maxResponsIdx]).m_memberIds;
	/* Add all members to the connected component. If they are
	 * active in the current orientation.
	 */
	for(size_t i = 0; i < membIds.size(); ++i) {
	  int tID = membIds[i];
	  int TubeIndex = hitMap.Find(tID);
	  GridNode &TubeToAdd = Ingrid[TubeIndex];
	  // If active in current orientation
	  if((TubeToAdd.m_orientations)[maxResponsIdx].m_act && (!TubeToAdd.m_orintVisited) ) {
	    TubesIdSet->insert(TubeToAdd.m_detID);
	    // Mark as visited
	    (Ingrid[TubeIndex]).m_orintVisited = true;
	    // if the maximum responce is the current responce value
	    // and remove from the queue
	    if(TubeToAdd.m_maxOrientIndex == maxResponsIdx) {
	      // Fully processed remove from the queue
	      ActiveQueue.deQueue(tID);
	    }
	  }// If active in current orientation
	}// Orient._Memebr loop
	// Add component to the output variable.
	outPutVar->push_back(TubesIdSet);
      }// (Cur_Node.m_maxOrientVal >= MinResponce)
      else {// (Cur_Node.m_maxOrientVal < MinResponce)
	    // Insert in short responce queue and dequeue from active
	ShortCompQueue.inQueue(Cur_Node.m_detID);
	// Mark as Visited
	Cur_Node.m_orintVisited = true;

#if (PATH_DEBUG_PRINT > 0)
	std::cerr<<"<WARNING>: Inserted short component to short queue."
		 << " id = " << Cur_Node.m_detID << '\n';
#endif
      }// END else (insert in short queue)
    }// END If not visited yet
    /* ELSE: It was visited before, Find in which list(s).
     *
     * OR (m_maxOrientVal < MinResponce), was added to the short queue
     */
    // List to keep track of indices of all components containing the
    // ID of the current node(tube)
    std::vector<size_t> ComponentIndexList;
    for(size_t cmpIdx = 0; cmpIdx < outPutVar->size(); ++cmpIdx){
      std::set<int> const* curConnectedlist = outPutVar->at(cmpIdx);
      std::set<int>::iterator FindIter = curConnectedlist->find(Cur_Node.m_detID);
      if( FindIter != curConnectedlist->end()) {
	ComponentIndexList.push_back(cmpIdx);
      }
    }
#if(CONNECTED_COMPONENT_DEBUG > 1)
    std::cout << "ComponentIndexList size = " << ComponentIndexList.size()
	      << '\n';
#endif
    for(size_t cmp = 0; cmp < ComponentIndexList.size(); ++cmp) {
      size_t index = ComponentIndexList[cmp];
      std::set<int>* ActiveConnectedComp = outPutVar->at(index);
      // Determine the values for Left and Left left orientations
      // (allow hyper connectivity in orientation space attr.)
      if(maxResponsIdx == 0) {// 360 deg
	LeftOrint = numOrientations - 1;
	LLeftOrint = LeftOrint - 1;
      }
      else if(maxResponsIdx > 0) {
	LeftOrint = maxResponsIdx - 1;
	LLeftOrint = (LeftOrint > 0)?(LeftOrint - 1):(numOrientations - 1);
      }
      else{// Overbodig
	LeftOrint = 0;
	LLeftOrint = numOrientations - 1;
      }
      // Right and Right Right.
      if(maxResponsIdx == numOrientations){
	RightOrient = 0;
	RRightOrient = RightOrient + 1;
      }
      else if(maxResponsIdx < (numOrientations - 1) ) {
	RightOrient = maxResponsIdx + 1;
	RRightOrient = (RightOrient == (numOrientations - 1)) ? 0 : (RightOrient + 1);
      }
      else{// Ovebodig kan niet gebeuren
	RightOrient = 0;
	RRightOrient = RightOrient + 1;
      }
      //_______________ LEFT orientation neighbors
      // Move to the left
      while( LeftOrint != 0) {
	std::vector<int> TmpLeftList;
	for(std::set<int>::iterator it = ActiveConnectedComp->begin(); it != ActiveConnectedComp->end(); ++it) {
	  int detecotId = *it;
	  int mapIdx    = hitMap.Find(detecotId);
	  GridNode &TMPLeftNode = Ingrid[mapIdx];
	  //Set as visited
	  Ingrid[mapIdx].m_orintVisited = true;
	  
	  // Add tubes to the temporary list
	  std::vector<int> const &partLeftLst = ( (TMPLeftNode.m_orientations)[LeftOrint]).m_memberIds;
#if(INCLUDE_LEFTLEFT_RIGHTRIGHT > 0)
	  std::vector<int> const &partLeftLeftLst = ( (TMPLeftNode.m_orientations)[LLeftOrint]).m_memberIds;
#endif
	  for(size_t k = 0; k < partLeftLst.size(); ++k) {
	    int tubeLeftID = partLeftLst[k];
	    int TubeIndex = hitMap.Find(tubeLeftID);
	    GridNode &TubeToAdd = Ingrid[TubeIndex];
	    // If active in orientation (left)
	    if( (TubeToAdd.m_orientations)[LeftOrint].m_act &&
		(!TubeToAdd.m_orintVisited) ) {
	      TmpLeftList.push_back(tubeLeftID);
	      // Set as visited.
	      (Ingrid[TubeIndex]).m_orintVisited = true;
	      // Remove from the queue
	      if( (Ingrid[TubeIndex]).m_maxOrientIndex == LeftOrint ) {
		ActiveQueue.deQueue(tubeLeftID);
	      }
	    }
	  }// FOR partLeftLst
#if(INCLUDE_LEFTLEFT_RIGHTRIGHT > 0)
	  // LeftLeft list
	  for(size_t n = 0; n < partLeftLeftLst.size(); ++n) {
	    int tubeLeftID = partLeftLeftLst[n];
	    int TubeIndex = hitMap.Find(tubeLeftID);
	    GridNode &TubeToAdd = Ingrid[TubeIndex];
	    // If active in orientation (leftleft)
	    if( (TubeToAdd.m_orientations)[LLeftOrint].m_act &&
		(!TubeToAdd.m_orintVisited) ) {
	      TmpLeftList.push_back(tubeLeftID);
	      // Set as visited.
	      (Ingrid[TubeIndex]).m_orintVisited = true;
	      // Remove from the queue
	      if( (Ingrid[TubeIndex]).m_maxOrientIndex == LLeftOrint ) {
		ActiveQueue.deQueue(tubeLeftID);
	      }
	    }
	  }// FOR partLeftLeftLst
#endif
	}// For it = ...
	
	// Insert IDs to the component set
	for(size_t t = 0; t < TmpLeftList.size(); ++t) {
	  ActiveConnectedComp->insert(TmpLeftList[t]);
	}
	LeftOrint--;
	
	// Update left left
	if( LeftOrint == 0) {// 360 deg
	  LLeftOrint = numOrientations - 1;
	}
	else {
	  LLeftOrint = LeftOrint - 1;
	}
	///////////////////////
      }// While LeftOrint

      // RIGHT orientation neighbors
      // Move to the right
      while( RightOrient < numOrientations) {
	std::vector<int> TmpRightLst;
	for(std::set<int>::iterator it = ActiveConnectedComp->begin(); it != ActiveConnectedComp->end(); ++it) {
	  int detecotId = *it;
	  int mapIdx    = hitMap.Find(detecotId);
	  GridNode &TMPRightNode = Ingrid[mapIdx];
	  Ingrid[mapIdx].m_orintVisited = true;
	  // Add tubes to the temporary list
	  std::vector<int> const &partRightLst = ( (TMPRightNode.m_orientations)[RightOrient]).m_memberIds;
#if(INCLUDE_LEFTLEFT_RIGHTRIGHT > 0)
	  std::vector<int> const &partRightRightLst = ( (TMPRightNode.m_orientations)[RRightOrient]).m_memberIds;
#endif
	  for(size_t k = 0; k < partRightLst.size(); ++k) {
	    int tID = partRightLst[k];
	    int TubeIndex = hitMap.Find(tID);
	    GridNode &TubeToAdd = Ingrid[TubeIndex];
	    // If active in orientation (right)
	    if((TubeToAdd.m_orientations)[RightOrient].m_act && (!TubeToAdd.m_orintVisited)
	      ) {
	      TmpRightLst.push_back(TubeToAdd.m_detID);
	      // Set as visited.
	      (Ingrid[TubeIndex]).m_orintVisited = true;
	      // Remove from the queue
	      if( (Ingrid[TubeIndex]).m_maxOrientIndex == RightOrient) {
		ActiveQueue.deQueue(TubeToAdd.m_detID);
	      }
	    }// If active in current orientation
	  }// Right neigbor list. Added all tubes in current
	   // orientation
#if(INCLUDE_LEFTLEFT_RIGHTRIGHT > 0)
	  // RightRight
	  for(size_t n = 0; n < partRightRightLst.size(); ++n) {
	    int tID = partRightRightLst[n];
	    int TubeIndex = hitMap.Find(tID);
	    GridNode &TubeToAdd = Ingrid[TubeIndex];
	     // If active in orientation (right)
	     if((TubeToAdd.m_orientations)[RRightOrient].m_act &&
	        (!TubeToAdd.m_orintVisited) ) {
	       TmpRightLst.push_back(TubeToAdd.m_detID);
	       // Mark as visited.
	       (Ingrid[TubeIndex]).m_orintVisited = true;
	       // Remove from the queue
	       if( (Ingrid[TubeIndex]).m_maxOrientIndex == RRightOrient) {
	     	ActiveQueue.deQueue(TubeToAdd.m_detID);
	       }
	    }// If active in current orientation
	  }// RightRight
#endif
	}// for (it = ) Visited all tubes in current component
	// Copy to the component ID set
	for(size_t t = 0; t < TmpRightLst.size(); ++t) {
	  ActiveConnectedComp->insert(TmpRightLst[t]);
	}
	RightOrient++;
	// Update right right
	if( RightOrient < (numOrientations - 1) ){
	  RRightOrient = RightOrient + 1;
	}
	else {
	  RRightOrient = 0;
	}
	//////
      }// While RightOrient
    }// END ComponentIndexList loop
  }// While (ActiveQueue)
  //______________________________________________________________________

  /*
   * If short queue not empty then insert ids to the tail of the
   * output list.
   */
  if( !ShortCompQueue.isEmpty() ){
    std::cout << "<INFO> Processing Short queue:\n"
	      << "      <-I-> Number of elements = " << ShortCompQueue.getNumElement()
	      << '\n';
    std::set<int>::iterator it;
    std::vector<int>::const_iterator find_it;
    for(size_t sh = 0; sh < outPutVar->size(); ++sh) {
      // Elements in short queue
      std::vector<int> const &qContainer = ShortCompQueue.m_queueCont;
      // Fetch connected component
      std::set<int> *CurrentComponent = outPutVar->at(sh);
      // Loop, participating tubes in the current component
      for(it = CurrentComponent->begin(); it != CurrentComponent->end(); ++it) {
	int tid = *it;
	// If current node has been added to another component.(double counting)
	find_it = std::find(qContainer.begin(), qContainer.end(), tid);
	// It is also in another connected component. Delete from short queue
	if( find_it != qContainer.end()) {
	  ShortCompQueue.deQueue(tid);
	}
	else {// Fetch node, if it is a direct neighbor of one of the
	      // elements in the current connected component.
	  size_t node_index = hitMap.Find(tid);
	  GridNode const &node = Ingrid[node_index];
	  // List of neigbors
	  std::vector<int> const &ListOf_Neighbors = node.m_neighbors;
	  // Neighbor loop
	  for( size_t nb = 0; nb < ListOf_Neighbors.size(); ++nb) {
	    // If direct neigbor in the short queue
	    find_it = std::find(qContainer.begin(), qContainer.end(), ListOf_Neighbors[nb]);
	    if( find_it != qContainer.end() ) {// Found in the short queue
	      int idToAdd = *find_it;// Fount neighbor
	      CurrentComponent->insert(idToAdd);// Add to the current component
	      ShortCompQueue.deQueue(idToAdd);//remove fromshort queue
	    }//END IF found in neighbor list
	  }//END for nb
	}// END ELSE
      }// Current component loop
    }// All components loop
#if (ADD_SHORTQUEUE_TO_TAIL > 0)
    std::cerr << "<DEBUG> Short id set added to the output.\n"
	      << "        <-I-> It contains still: " << ShortCompQueue.getNumElement()
	      << " orphan elements.\n";
    if(!ShortCompQueue.isEmpty()) {
      std::set<int>* shortSet = new std::set<int>();
      while( !ShortCompQueue.isEmpty()) {
	shortSet->insert(ShortCompQueue.popFront());
      }
      outPutVar->push_back(shortSet);
    }
#endif
  }
  return outPutVar;
}
//_________________ END AttSpaceConnectedComp ________________________
//_____________________ AttSpaceConnectedComp_Union __________________
std::vector< std::set<int>* >* AttSpaceConnectedComp_Union( CoordGrid &hitMap, size_t MinResponce)
{
  std::cout << "<INFO> Connected component analysis(Union Find), MinPathOrint = "
	    << MinResponce
	    << '\n';
  // Init output variable
  std::vector< std::set<int>* >* outPutVar = 0;
  outPutVar = new std::vector< std::set<int>* >();
  if(!outPutVar){
    std::cerr << "<ERROR> Could not allocate memory for connected components\n";
    assert (outPutVar != 0);
  }
  // Get the map nodes
  std::vector< GridNode > &Ingrid = hitMap.m_grid;
  std::cout << "Num Element = " << Ingrid.size() << '\n';
  return outPutVar;
}
//_________________ END AttSpaceConnectedComp_Union __________________
//_____________________ ComputeOrientationsMVD _______________________
void FindCompatibleMVD(CoordGrid &hitMap, GridNode const &tube,
			    float angle, std::vector<int> &indxs)
{
  std::cout << "<INFO> Finding compatible MVD nodes for current components.\n";
  // Fetch the MVD points for the current snapshto. Each grid instance
  // contains one snapshot.
  std::vector< GridNode > &MVDPoints = hitMap.m_MVD_grid;
  float nodex,nodey;
  nodex = tube.m_x;
  nodey = tube.m_y;
  std::pair<float,float> Rad_deg;
  for(size_t i = 0; i < MVDPoints.size(); ++i) {
    GridNode &mvdNode = MVDPoints[i];
    if( !mvdNode.m_MVDAssigned ) {// Not assigned before
      ComputeSlope(nodex, nodey, mvdNode.m_x, mvdNode.m_y, Rad_deg);
      if( (Rad_deg.first - angle) < 0.2) {// Set acceptance for MVD
	indxs.push_back(mvdNode.m_detID);
	mvdNode.m_MVDAssigned = true;
	nodex = mvdNode.m_x;
	nodey = mvdNode.m_y;
      }
    }// If not assigned before
  }// MVD node loop
}
//_____________________ END ComputeOrientationsMVD ____________________
//______________________ MergeConnectedComponentsWithMVD ______________
// FIXME Deze heb je in commentaar gezet
std::vector<TrackObject*>* MergeConnectedComponentsWithMVD(CoordGrid &hitMap,
							   std::vector< std::set<int>* >* SttcompList)
{
  std::cout << "<INFO> Merge STT components with MVD nodes.\n";
  // If the number of track candidates is zero.
  if( !SttcompList || (SttcompList->size() == 0) ){
    std::cerr << "<ERROR> Connected component list is empty.\n";
    return 0;
  }
  
  // Fetch the MVD points for the current snapshto. Each grid instance
  // contains one snapshot.
  // std::vector< GridNode > &MVDPoints = hitMap.m_MVD_grid;
  
  // Sort MVD points in increasing order.
  // std::sort(MVDPoints.begin(), MVDPoints.end());
  // std::reverse(MVDPoints.begin(), MVDPoints.end());
  // std::cout << "\n\nTotal number of MVD = " << MVDPoints.size() << "\n\n";
  
  // Determine orientation space attributes for MVD nodes.
  // MARKED
  // ComputeMVD_OrientationAttr(hitMap);

  // Get the map nodes(STT tubes)
  std::vector< GridNode > &Ingrid = hitMap.m_grid;

  // Allocate memory for output.
  std::vector<TrackObject*>* outPutTracks = new std::vector<TrackObject*>();

  // Connected component loop. Loop and merge
  /*
  std::set<int>::iterator it;
  size_t index_Nearest_mvd = 0;
  float  shortestDist = std::numeric_limits<float>::max();
  */
  for(size_t i = 0; i < SttcompList->size(); ++i) {
    // Is already sorted increasing. Set property
    //std::set<int>* cur_comp = SttcompList->at(i);
    /*
    // Inner most tube id in current component. Note ids are circular
    // assigned.
    it = cur_comp->begin();
    // Fetch the actual tube (grid node)
    int InMapIdex = hitMap.Find(*it);
    GridNode &InnerTube = Ingrid[InMapIdex];

#if(MVD_MERGE_DEBUG_PRINT > 0)
    std::cout << "\t<-I-> STT Component number " << i
	      << " inner tube ID = " << InnerTube.m_detID
	      << " ";
#endif
*/
    /*
      float MaxOrientAngle = InnerTube.m_orientations[InnerTube.m_maxOrientIndex].m_angle;
      std::vector<int> mvdIndices;
      FindCompatibleMVD(hitMap, InnerTube, MaxOrientAngle, mvdIndices);
    */
    // Find MVD node with the minimum distance to the current tube
    // Reset distance temporary var
    /*
    shortestDist = std::numeric_limits<float>::max();
    for( size_t j = 0; j < MVDPoints.size(); ++j) {
      GridNode &mvd = MVDPoints[j];
      if ( InnerTube.GetDistance2D(mvd) < shortestDist ) {
	shortestDist = InnerTube.GetDistance2D(mvd);
	index_Nearest_mvd = j;
      }
    }//END MVD loop
    // We have found the mvd detector with the shortest distance to
    // the most inner STT tube of the current track (connected component)
    GridNode &mvdDetector = MVDPoints[index_Nearest_mvd];
#if(MVD_MERGE_DEBUG_PRINT > 0)
    std::cout << " nearest mvd id = " << mvdDetector.m_detID << std::endl;
#endif
  */
    TrackObject* trk = new TrackObject();
    trk->m_sttComponent = new std::set<int>(*SttcompList->at(i));
    // MARKED
    //trk->m_MVD_Component = ( (mvdDetector.m_orientations)[mvdDetector.m_maxOrientIndex]).m_memberIds;
    outPutTracks->push_back(trk);
  }

  // Print MVD debug info
#if(MVD_MERGE_DEBUG_PRINT > 1)
  std::cout << "<DEBUG> MVD coordinates:\n";
  for(size_t i = 0; i < MVDPoints.size(); ++i) {
    GridNode const &node = MVDPoints[i];
    std::cout << "\t<-D-> ID = " << node.m_detID 
	      << " : x = " << node.m_x << ", "
	      << " y = " << node.m_y << ", "
	      << " z = " << node.m_z
	      << '\n';
  }
#endif
  return outPutTracks;
}
//______________________ END MergeConnectedComponentsWithMVD __________
//__________________________ ComputeMVD_OrientationAttr _______________
void ComputeMVD_OrientationAttr(CoordGrid &hitMap)
{
  std::cout << "<INFO> Determine orientation space attribute values for MVD.\n";
  // Find the list of angles from STT tube set
  std::vector< GridNode > &Ingrid = hitMap.m_grid;
  
  if(Ingrid.size() == 0) {
    std::cerr << "<ERROR> STT tube graph is empty. The program will terminate NOW.\n";
    exit(EXIT_FAILURE);
  }
  std::vector< GridNode > &MVDPoints = hitMap.m_MVD_grid;
  std::vector<float> OrientValList;
  
  size_t index = 0;
  while( (index < Ingrid.size()) &&
	 !(Ingrid[index].m_active)// Not reponding tube
	 ) {
    index++;
  }
  // Found first active tube
  GridNode const &tube = Ingrid[index];
  std::vector< NodeOrientation > const &Node_orinetList = tube.m_orientations;
  
  for(size_t i = 0; i < Node_orinetList.size(); ++i) {
    OrientValList.push_back(Node_orinetList[i].m_angle);
  }
  // Init container for MVD nodes.
  for(size_t j = 0; j < MVDPoints.size(); ++j) {
    MVDPoints[j].initNodeOrientation(OrientValList.size());
  }
  // Set MVD node orientations angles values
  for(size_t orIdx = 0; orIdx < OrientValList.size(); ++orIdx) {
    float curr_orient = OrientValList[orIdx];
    for(size_t gn = 0; gn < MVDPoints.size(); ++gn) {
      GridNode &node = MVDPoints[gn];
      NodeOrientation &ori = node.m_orientations[orIdx];
      ori.m_angle = curr_orient;
      std::vector<int> &memList = ori.m_memberIds;
      memList.push_back(node.m_detID);
    }
  }
  // Dummy to hold the results of slope determination
  
  std::pair<float,float> Rad_deg;
  std::vector<int>::iterator it;
  
  for(size_t i = 0; i < OrientValList.size(); ++i) {
    float Curr_orient = OrientValList[i];
    // Which combination has the current orienation
    for(size_t j = 0; j < MVDPoints.size(); ++j) {
      GridNode &first = MVDPoints[j];
      for(size_t k = 0; k < MVDPoints.size(); ++k) {
	if(k != j) {
	  GridNode &second = MVDPoints[k];
	  ComputeSlope(first.m_x, first.m_y, second.m_x, second.m_y, Rad_deg);
	  if( fabs(Rad_deg.first - Curr_orient) < 0.03 ) {// Set tol for MVD (MVDTOL)
	    // Update both nodes.
	    it = std::find( first.m_orientations[i].m_memberIds.begin(),
			    first.m_orientations[i].m_memberIds.end(), second.m_detID);
	    if( it == first.m_orientations[i].m_memberIds.end()) {
	      first.m_orientations[i].m_memberIds.push_back(second.m_detID);
	    }
	    it = std::find( second.m_orientations[i].m_memberIds.begin(),
			    second.m_orientations[i].m_memberIds.end(), first.m_detID);
	    if( it == second.m_orientations[i].m_memberIds.end()) {
	      second.m_orientations[i].m_memberIds.push_back(first.m_detID);
	    }
	  }// IF fabs
	}// IF k != j
      }// FOR K = 
    }// FOR J = 
  }//FOR orientation list
  // All nodes have been visited for all orientations.
  // Update lists for all MVD detectors.
  for(size_t i = 0; i < MVDPoints.size(); ++i) {
    GridNode &first = MVDPoints[i];
    for(size_t j = 0; j < MVDPoints.size(); ++j) {
      if(i != j) {
	GridNode &second = MVDPoints[j];
	std::vector< NodeOrientation > &orientList = second.m_orientations;
	for(size_t k = 0; k < orientList.size(); ++k) {
	  NodeOrientation &currentOrient = orientList[k];
	  std::vector<int> &memList = currentOrient.m_memberIds;
	  it = std::find(memList.begin(), memList.end(), first.m_detID);
	  if( it != memList.end()){// Element available in the second list
	    for(size_t l = 0; l < memList.size(); ++l) {
	      int id_toAdd = memList[l];
	      it = std::find(first.m_orientations[k].m_memberIds.begin(),
			     first.m_orientations[k].m_memberIds.end(), id_toAdd);
	      if( it == memList.end() ){// Not added before
		first.m_orientations[k].m_memberIds.push_back(id_toAdd);
	      }
	    }// FOR L
	  }// IF found
	}// FOR K
      }// IF
    }// FOR J
  }// FOR I
  // Update node orientation info.
  size_t maxOrientValue = 0;
  size_t maxOrientIndex = 0;
  size_t minOrientValue = 0;
  size_t minOrientIndex = 0;
  
  for(size_t i = 0; i < MVDPoints.size(); ++i) {
    GridNode &mvdNode = MVDPoints[i];
    std::vector< NodeOrientation > &orientList = mvdNode.m_orientations;
    maxOrientValue = 0;
    minOrientValue = std::numeric_limits<size_t>::max();
    for(size_t j = 0; j < orientList.size(); ++j) {
      NodeOrientation &currentOrient = orientList[j];
      if( (currentOrient.m_memberIds.size()) > maxOrientValue ) {//Maximum
	maxOrientValue = (currentOrient.m_memberIds).size();
	maxOrientIndex = j;
      }// MAX
      if( (currentOrient.m_memberIds.size()) < minOrientValue ) {//minimum
	minOrientValue = (currentOrient.m_memberIds).size();
	minOrientIndex = j;
      }//MIN
    }// END orientations loop
    // Update node
    mvdNode.m_maxOrientIndex = maxOrientIndex;
    mvdNode.m_maxOrientVal   = maxOrientValue;
    mvdNode.m_minOrientIndex = minOrientIndex;
    mvdNode.m_minOrientVal   = minOrientValue;
  }// END MVD node loop
  //_____ DEBUG INFO
#if ( MVD_MERGE_DEBUG_PRINT > 1)
  std::cout << "<DEBUG> MVD Orientations.\n";
  for(size_t i = 0; i < MVDPoints.size(); ++i) {
    GridNode &node = MVDPoints[i];
    std::cout << "\t<-D-> ID = " << node.m_detID << '\n';
    std::vector< NodeOrientation > &orientList = node.m_orientations;
    for(size_t j = 0; j < orientList.size(); ++j) {
      NodeOrientation &currentOrient = orientList[j];
      std::cout << currentOrient.m_angle << " rad, ";
      std::vector<int> &memList = currentOrient.m_memberIds;
      for(size_t l = 0; l < memList.size(); ++l) {
	std::cout << memList[l] << " ";
      }
      std::cout << std::endl;
    }
  }
#endif
}
//______________________ END ComputeMVD_OrientationAttr _______________
//_____________________________ TrackZ_CoordinatesDistNorm ____________
/*
 * Determine the track z coordinates normalized by the distance
 * between tubes. In z direction the step size is constant as function
 * of the distance between the points (nodes). Again hier the virtual
 * tubes are not modified.
 */
//_________________________ END TrackZ_CoordinatesDistNorm ____________
//_____________________________ MergeConnectedComponents   ____________
void MergeConnectedComponents(CoordGrid const &hitMap, std::vector< std::set<int>* >* sttTrks)
{
  if( (sttTrks == 0) || (sttTrks->size() == 0) ) {
    std::cerr << "<ERROR> Input track list is empty. Terminating.\n";
    exit(EXIT_FAILURE);
  }
  // STT lists of elements available in the graph.
  std::vector< GridNode > const &Ingrid = hitMap.m_grid;
  std::cout << "Ingrid.size() = " << Ingrid.size() << '\n';
}
//_________________________ END MergeConnectedComponents   ____________
//____________________________ SplitSharedNodes _______________________
void SplitSharedNodes(CoordGrid &hitMap, size_t threshold, std::vector < GridNode > &outNodes)
{
  std::cout << "<INFO> Split shared nodes.\n"
	    << "       With threshold = " << threshold
	    << '\n';
  // Detector graph
  std::vector< GridNode > &SttList = hitMap.m_grid;
  int idToAdd = 80000;
  //
  for(size_t i = 0; i < SttList.size(); ++i) {
    GridNode &Current_Node = SttList[i];
    size_t MaxOrientIndex = Current_Node.m_maxOrientIndex;
    size_t MaxOrientVal   = Current_Node.m_maxOrientVal;
    std::vector< NodeOrientation > &OrientationList = Current_Node.m_orientations;
    if( MaxOrientVal >= threshold) {
      std::vector<size_t> OrientIndices;
      for(size_t j = 0; j < OrientationList.size(); ++j) {
	NodeOrientation &orintation = OrientationList[j];
	if( ( (orintation.m_memberIds).size() >= threshold) &&
	    (orintation.m_act) ){
	  OrientIndices.push_back(j);
	}
      }// END for J
      // Now we need to split nodes.
      for(size_t k = 0; k < OrientIndices.size(); ++k) {
	GridNode splitNode = Current_Node;
	// Modify orientations list
	std::vector< NodeOrientation > &OrntLst = splitNode.m_orientations;
	splitNode.m_Orig_detID = idToAdd;
	splitNode.m_detID = idToAdd;
	splitNode.m_type = GridNode::STT_TYPE_SPLIT_SKEW;
	splitNode.m_maxOrientIndex = OrientIndices[k];
	splitNode.m_maxOrientVal   = (OrntLst[splitNode.m_maxOrientIndex]).m_memberIds.size();
	for(size_t l = 0; l < OrntLst.size(); ++l) {
	  if( l != splitNode.m_maxOrientIndex) {
	    (OrntLst[l]).m_act = false;
	    //(OrntLst[l]).m_memberIds.clear();
	  }// IF L
	}// Split node orient list loop
	// Add node to grid.
	outNodes.push_back(splitNode);
	idToAdd++;
      }// END Indices list loop
    }// If higher than "threshold"
    for(size_t j = 0; j < OrientationList.size(); ++j) {
      if( j != MaxOrientIndex) {
	OrientationList[j].m_act = false;
	//OrientationList[j].m_memberIds.clear();
      }
    }
  }// END for I

  std::cout << "<INFO> Done splitsing. Added "<< (idToAdd - 80000) << " nodes.\n";
}
//________________________ END SplitSharedNodes _______________________
