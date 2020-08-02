/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <set>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cassert>
#include <utility>
#include <limits>
#include <fstream>

// LOCAL
#include "utilfunctions.h"
#include "CoordGrid.h"
#include "gridNode.h"
#include "pathCandidate.h"

// ROOT
#include "TNtuple.h"
#include "TVector3.h"

// Local define
#define UTILITY_DEBUG_PRINT 0
#define UTILITY_WARNING_PRINT 1

/**
 * Polarout <R, theta>: To compute the value, the function takes into
 * account the sign of both arguments in order to determine the
 * quadrant.
 */
double Cartesian_To_Polar(float const x, float const y,
                          std::pair<float,float>& polarOut,
                          bool useSign)
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
//___________________ END of Cartesian_To_Polar ______________________
//____________________ComputeSlope ___________________________________
/**
 * Computes angle of a line between two given points in the space and
 * the x-axis. NOte that the output is in the range of [0, pi).
 *
 *@param x1,y1 define x and y coordinates of the first point.
 *@param x2,y2 define x and y coordinates of the second point.
 *
 *@param Rad_deg Holds the output in degrees and
 * radiants. Rad_deg.first = radiant, Rad_deg.second = degrees.
 */
void ComputeSlope(float const x1, float const y1,
                  float const x2, float const y2,
                  std::pair<float,float>& Rad_deg,
		  float epsilon)
{
  float DeltaX = x2 - x1;
  float DeltaY = y2 - y1;
  float theta;// = 0.00;
  
  ///++++++++++++++++++++++++++++++++++++++++++
  // x == 0
  // if( !(DeltaX > 0.0) && !(DeltaX < 0.0) ){
  //   theta = (M_PI / 2.0);
  // }
  // else if( (DeltaY/ DeltaX) >= 0 ) {
  //   theta = atan(DeltaY/DeltaX);
  // }
  // else {
  //   theta = M_PI + (atan(DeltaY/DeltaX));
  // }
  ///++++++++++++++++++++++++++++++++++++++++++
  // x == 0 (dx == 0)
  // if( !(DeltaX > 0.0) && !(DeltaX < 0.0) ){
  if ( fabs(DeltaX) < LOCAL_PND_TRACKING_EPSILON ) {
    theta = (M_PI / 2.0);
  }
  else {
    theta = atan2(DeltaY, DeltaX);
  }
  
  if(theta < 0){
    theta += M_PI;
  }
  
  float tolerance = (epsilon < 0.0) ? LOCAL_PND_TRACKING_EPSILON : epsilon;
  // Theta == PI, thus set theta to zero (mapping to upper half
  // circle)
  if( fabs(theta - M_PI) < tolerance) {
    theta = 0.0;
  }

  // Fill the output object.
  Rad_deg.first  = theta;
  Rad_deg.second = ( (theta * 180.00) / M_PI);

#if(UTILITY_DEBUG_PRINT > 1)
  float radius = sqrt( (DeltaX * DeltaX) + (DeltaY * DeltaY));
  std::cout << " Radius = " << radius
            << " angle = " << ( (theta * 180) / M_PI)
            << " deg and " << theta << " radiants.\n";
#endif
}
//__________________________ END OF ComputeSlope _______________________
//____________________ FindNodeBetweenLayerNodePairs  ___________________
void FindNodeBetweenLayerNodePairs(CoordGrid const &hitMap,
				   std::vector< TubeLayerPairProperty > &NodePairSet)
{
  std::cout << "<INFO> Finding tube pairs to compute virtual nodes.\n";

  // List of all nodes available in the image (detector map)
  std::vector< GridNode > const &Ingrid = hitMap.m_grid;
  size_t layerDiff = 0;
  
  // Grid node loop to find the possible pairs.
  for(size_t i = 0; i < Ingrid.size(); ++i) {
    GridNode const &current_Node = Ingrid[i];
    // List of neighbors
    std::vector<int> const &neighList = current_Node.GetNeighbors();
    // Neighbors loop
    for( size_t j = 0; j < neighList.size(); ++j) {
      int neighbourID = neighList[j];// Neighbour ID
      int neigh_Index = hitMap.Find(neighbourID);// Index in the map
      // Neighbor node
      GridNode const &neigh_node = Ingrid[neigh_Index];
      /*
       * If the nodes have different slopes and not in the same layer.
       */
      // Not the same slope
      if( (current_Node.m_Slope > neigh_node.m_Slope) ||
	  (current_Node.m_Slope < neigh_node.m_Slope) ) {
	// Not the same layer
	layerDiff = (current_Node.m_Layer > neigh_node.m_Layer) ?
	            (current_Node.m_Layer - neigh_node.m_Layer) :
	            (neigh_node.m_Layer - current_Node.m_Layer);
	// Exactly one layer difference
	if( layerDiff == 1 ) {
	  // Create a pair
	  TubeLayerPairProperty pp;
	  pp.firstNodeID     = current_Node.m_detID;
	  pp.firstNodeIndex  = i;
	  pp.secondNodeID    = neighbourID;
	  pp.secondNodeIndex = neigh_Index;
	  // Add to pairs list
	  NodePairSet.push_back(pp);
	}
      }// If not same slope
    }// Neighbour loop
  }// Grid node loop
  /*
   * Now we need to mark duplicates as invalid combinations. Not so
   * efficient. Since geometry is nog changing, we don't realy
   * mind.
   */
  for(size_t f = 0; f < NodePairSet.size(); ++f) {
    // First pair
    TubeLayerPairProperty &first = NodePairSet[f];
    for(size_t s = (f + 1); s < NodePairSet.size(); ++s) {
      // Second pair
      TubeLayerPairProperty &second = NodePairSet[s];
      // If pairs are equivalent mark the second as invalid
      if( equalTubePairProperty(first, second) ) {
	// Mark the second as invalid
	second.isValid = false;
      }
    }
  }
  size_t cnt = 0;
  for(size_t i = 0; i < NodePairSet.size(); ++i) {
    if( NodePairSet[i].isValid ) {
      cnt++;
    }
  }
  std::cout << "\t<-I-> Found " << NodePairSet.size() << " Pairs.\n"
	    << "\t<-I-> and " << cnt << " valid pairs.\n";
}
//____________________ END FindNodeBetweenLayerNodePairs _____________
//____________________ Compute_Add_VirtualNodesNeigbor_________________
void Compute_Add_VirtualNodes_Neigbor(CoordGrid &hitMap, std::vector < GridNode > &VNodes)
{
  std::cout << "<INFO> Computing grid virtual tubes (Neigbor based).\n";
  // Dummy local counters.
  int NumTubesAdded  = 0;
  // Start node ID for virtual nodes.
  int StartVirtualID = START_VIRTUAL_ID;//6000;

  // List of all nodes available in the detector map.
  std::vector< GridNode > &Ingrid = hitMap.m_grid;
  // Sort all nodes by their layer number
  std::cout << "\t<-I-> Sorting nodes, increasing layer number.\n";
  std::sort(Ingrid.begin(), Ingrid.end(), LessThanLayer);
  // Determine the pairs of nodes between the layers for which we want
  // to compute virtual nodes and discard duplicates.
  std::vector< TubeLayerPairProperty > NodePairSet;
  FindNodeBetweenLayerNodePairs(hitMap, NodePairSet);

  // Loop through the possible combinations and create virtual nodes
  // if the actual tubes virtually intersect.
  for(size_t i = 0; i < NodePairSet.size(); ++i) {
    TubeLayerPairProperty const &prop = NodePairSet[i];
    if(prop.isValid) {
      GridNode  &firstNode  = Ingrid[ prop.firstNodeIndex ];
      GridNode  &secondNode = Ingrid[ prop.secondNodeIndex ];
      // Create dummy virtual
      GridNode Dummy_coord;
      // Find intersection point.
      if( IntersectionPoint_NeigborList(hitMap, firstNode, secondNode, Dummy_coord) ) {
	// Modify node ID
	Dummy_coord.m_detID      = StartVirtualID + NumTubesAdded;
	Dummy_coord.m_Orig_detID = Dummy_coord.m_detID;
	/*
	  Dummy_coord.m_halfLength = 0;
	  Dummy_coord.m_type = GridNode::VIRTUAL_NODE;
	  // Add parents to the neigboring list
	  (Dummy_coord.m_neighbors).clear();
	  (Dummy_coord.m_neighbors).push_back(firstNode.m_detID);
	  (Dummy_coord.m_neighbors).push_back(secondNode.m_detID);
	*/
	/*
	 * We need to fix neighboring relations later using using
	 * fixneighbouring function. By adding here the virtual
	 * neighbours of the actual nodes are not added to the grid
	 * yet.
	 */
	// Add to output.
	VNodes.push_back(Dummy_coord);
	NumTubesAdded++;
      }// If intersect
    }// END if the pair is valid
  }// End for NodePairSet
  // Reset visiting variable for all nodes in the input graph
  for(size_t i = 0; i < Ingrid.size(); ++i) {
    Ingrid[i].m_visited = false;
  }
  // Print Info
  std::cout << "<INFO> Determined "<< NumTubesAdded
	    << " virtual tubes (neighborList)."
	    << " Num possible V_nodes = " << NodePairSet.size()
	    << '\n';
}
//____________________ END Compute_Add_VirtualNodes_Neigbor ______________
//____________________ Compute_Add_VirtualNodesNeigbor2 _________________
void Compute_Add_VirtualNodes_Neigbor2(CoordGrid &hitMap, std::vector < GridNode > &VNodes)
{
  std::cout << "<INFO> Computing grid virtual tubes (Neigbor based OLD).\n";
  int NumTubesAdded  = 0;
  size_t numPossibleVirtuals = 0;
  int StartVirtualID = START_VIRTUAL_ID;
  /////////////////
  /*
    std::vector< TubeLayerPairProperty > NodePairSet;
    FindNodeBetweenLayerNodePairs(hitMap, NodePairSet);
  */
  ////////////////
  // List of all nodes available in the image (detector map)
  std::vector< GridNode > &Ingrid = hitMap.m_grid;
  // Sort nodes by layer
  std::cout << "\t<-I-> Sorting nodes, increasing layer number.\n";
  std::sort(Ingrid.begin(), Ingrid.end(), LessThanLayer);

  // Loop through the tubes and compute virtual nodes.
  for(size_t i = 0; i < Ingrid.size(); ++i) {
    GridNode &current_Node = Ingrid[i];
    // List of neighbors
    std::vector<int> &neighList = current_Node.GetNeighbors();
    // Fetch neigbors
    for( size_t j = 0; j < neighList.size(); ++j) {
      int neighbourID = neighList[j];
      int neigh_Index = hitMap.Find(neighbourID);
      // Neighbor node
      GridNode &neigh_node = Ingrid[neigh_Index];
      /*
       * Check if this is correct. Maybe we can reduce the number of
       * virtuals or we are calculating too many virtuals.
       */
      if( ( (current_Node.m_Slope > neigh_node.m_Slope) ||// Not the same direction
	    (current_Node.m_Slope < neigh_node.m_Slope) ) &&
	  ( labs(current_Node.m_Layer - neigh_node.m_Layer) == 1) ) {// Not the same layer
	current_Node.m_visited = true;
	neigh_node.m_visited   = true;
	numPossibleVirtuals++;
	// Create dummy virtual
	GridNode Dummy_coord;
	// Find intersection point.
	if( IntersectionPoint_NeigborList(hitMap, current_Node, neigh_node, Dummy_coord) ) {
	  //////if( IntersectionPoint_SectorList(hitMap, current_Node, neigh_node, Dummy_coord) ) {
	  // Modify node ID
	  Dummy_coord.m_detID      = StartVirtualID + NumTubesAdded;
	  Dummy_coord.m_Orig_detID = Dummy_coord.m_detID;
	  /*
	   * We need to fix neighboring relations later using
	   * fixneighbouring function. By adding here the virtual
	   * neighbours of the actual nodes are not added to the grid
	   * yet.
	   */
	  // Add to output.
	  VNodes.push_back(Dummy_coord);
	  NumTubesAdded++;
	}// If intersect
      }// If diff slope and one layer diff
    }// Loop Neighbor list
  }// Grid loop
  // Reset visiting variable for all nodes in the input graph
  for(size_t i = 0; i < Ingrid.size(); ++i) {
    Ingrid[i].m_visited = false;
  }
  // Print Info
  std::cout << "<INFO> Determined "<< NumTubesAdded
	    << " virtual tubes (neighborList)."
	    << " Num possible V_nodes = " << numPossibleVirtuals
	    << '\n';
}
//____________________ END Compute_Add_VirtualNodes_Neigbor2 ______________
//__________________________ IntersectionPoint_neigborList__________________
/*
 * FIXME; This method computes wrong coordinates. Needs to be
 * recalculated and checked !!!!!?????????. We compute wrong linear
 * translation, by searching for the nearest tube in the other
 * plane. This way, the translation is not orthogonal to the plane and
 * it contains shifts. This shift will affect the values computed for
 * the translated point in the x,y plane and consequently it will
 * affect the values for z and determination of, if the lines
 * intersect at all. This issue needs to be fixed for better tracking
 * and also vor the computed z-coordinates.(VERY IMPORTANT FIXME
 * FIXME).
 */
bool IntersectionPoint_NeigborList(CoordGrid const &hitMap,
				   GridNode &tubeA, GridNode &tubeB,
				   GridNode &out)
{
  // List of all nodes available in the image (detector map)
   std::vector< GridNode > const &Ingrid = hitMap.m_grid;

  // Find the colsest point in tubeB's layer and sector to tube A
  size_t indexMinDist = 0;
  float minDist = std::numeric_limits<float>::max();
  float tmpDist = 0;

  // Neighbor list of the tubeB
  std::vector<int> const& NeighborsB = tubeB.GetNeighbors();
  // Find nearest tube to A ( Do we really need to do this?)
  for(size_t i = 0; i < NeighborsB.size(); ++i) {
    GridNode const &current_Node = Ingrid[ hitMap.Find(NeighborsB[i]) ];
    if( (current_Node.m_Sector == tubeB.m_Sector) &&// Same sector as B
	(current_Node.m_Layer  == tubeB.m_Layer )// Same layer as B
      ){
      tmpDist = (current_Node.m_x - tubeA.m_x) * (current_Node.m_x - tubeA.m_x) +
	        (current_Node.m_y - tubeA.m_y) * (current_Node.m_y - tubeA.m_y) +
	        (current_Node.m_z - tubeA.m_z) * (current_Node.m_z - tubeA.m_z);
      if(tmpDist < minDist) {
       	minDist = tmpDist;
       	indexMinDist = hitMap.Find(NeighborsB[i]);
      }
    }
  }
  // Found the nearest tube to A in tubeB.sector and layer
  float dx = Ingrid[indexMinDist].m_x - tubeA.m_x;
  float dy = Ingrid[indexMinDist].m_y - tubeA.m_y;

  // translate A to the tubeB.sector (xy-Plane)
  GridNode TransA(tubeA);// Dit kan beter maar voor nu .... FIXME Later
  TransA.m_x = tubeA.m_x + dx;
  TransA.m_y = tubeA.m_y + dy;

  // Compose the line equations for (X,Y,Z) for each center point and
  // direction vector of the tube.
  TVector3 directionA = TVector3(tubeA.m_WireDirection);// Direction vector of tube A
  TVector3 directionB = TVector3(tubeB.m_WireDirection);// Direction vector of tube B
  /*
   * x = x_a + r(x_v)
   * y = y_a + r(y_v)
   * z = z_a + r(z_v)
   * Use z_v, we know it is never zero.
   */
  float A = tubeB.m_y - TransA.m_y + TransA.m_z * (directionA.Y()/directionA.Z()) - 
            tubeB.m_z * (directionB.Y()/directionB.Z());
  float B = (directionA.Y()/directionA.Z()) - (directionB.Y()/directionB.Z());
  float Z_intersect = A/B;
  float r = (Z_intersect - TransA.m_z)/directionA.Z();
  float s = (Z_intersect - tubeB.m_z)/directionB.Z();
  float X_intersect = tubeB.m_x +  s * directionB.X();
  float Y_intersect = TransA.m_x + r * directionA.X();

  // Maybe better to use only the z-Difference???
  float DistA = sqrt( (TransA.m_x - X_intersect) * (TransA.m_x - X_intersect) +
		      (TransA.m_y - Y_intersect) * (TransA.m_y - Y_intersect) +
		      (TransA.m_z - Z_intersect) * (TransA.m_z - Z_intersect)
		      );
  float DistB = sqrt( (tubeB.m_x - X_intersect) * (tubeB.m_x - X_intersect) +
		      (tubeB.m_y - Y_intersect) * (tubeB.m_y - Y_intersect) +
		      (tubeB.m_z - Z_intersect) * (tubeB.m_z - Z_intersect)
		      );
  // !(>) => (<=, using <= is not safe for floats)
  if( !( DistA > tubeA.m_halfLength ) &&
      !( DistB > tubeB.m_halfLength ) 
      ) {
    // FIXME FIXME Not really optimal (correct??)
    TransA.m_x     = (tubeA.m_x + tubeB.m_x)/2.0;
    TransA.m_y     = (tubeA.m_y + tubeB.m_y)/2.0;
    TransA.m_z     = Z_intersect;
    TransA.m_z_Det = Z_intersect;
    // Set node as virtual and weight becomes 0 (no added value for
    // the length or area size).
    TransA.m_type        = GridNode::VIRTUAL_NODE;
    TransA.m_weight      = 0;
    TransA.m_SectorLimit = 0;
    TransA.m_LayerLimit  = false;
    TransA.m_halfLength  = 0;// A point. No half length
    // Add parents to the neigboring list
    (TransA.m_neighbors).clear();
    (TransA.m_neighbors).push_back(tubeA.m_detID);
    (TransA.m_neighbors).push_back(tubeB.m_detID);


    (tubeA.m_neighbors).erase(std::remove((tubeA.m_neighbors).begin(), (tubeA.m_neighbors).end(), tubeB.m_detID), (tubeA.m_neighbors).end());
    (tubeB.m_neighbors).erase(std::remove((tubeB.m_neighbors).begin(), (tubeB.m_neighbors).end(), tubeA.m_detID), (tubeB.m_neighbors).end());
    // Fill the output parameter.
    out = TransA;
    //// DEBUG DEBUG
    // std::cout << "TubeA: (" << tubeA.m_x << ", " << tubeA.m_y << ", " << tubeA.m_z << ")\n"
    // 	      << "TubeB: (" << tubeB.m_x << ", " << tubeB.m_y << ", " << tubeB.m_z << ")\n"
    // 	      << "Out: ("   << out.m_x << ", " << out.m_y << ", " << out.m_z << ")\n\n";
    ////////////
    return true;
  }// END IF
  return false;
}
//__________________________ END IntersectionPoint_neigborList _____________
//__________________________ IntersectionPoint_SectorList___________________
bool IntersectionPoint_SectorList(CoordGrid const &hitMap, GridNode &tubeA, GridNode &tubeB,
				  GridNode &out)
{
  // List of all nodes available in the image (detector map)
   std::vector< GridNode > const &Ingrid = hitMap.m_grid;

  // Find the colsest point in tubeB layers and sector to tube A
  size_t indexMinDist = 0;
  float minDist = std::numeric_limits<float>::max();
  float tmpDist = 0;
  for(size_t i = 0; i < Ingrid.size(); ++i) {
    GridNode const &current_Node = Ingrid[i];
    if( (current_Node.m_Sector == tubeB.m_Sector) &&// Same sector as B
	(current_Node.m_Layer  == tubeB.m_Layer )// Same layer as B
      ) {
      tmpDist = (current_Node.m_x - tubeA.m_x) * (current_Node.m_x - tubeA.m_x) +
	        (current_Node.m_y - tubeA.m_y) * (current_Node.m_y - tubeA.m_y) +
	        (current_Node.m_z - tubeA.m_z) * (current_Node.m_z - tubeA.m_z);
      if(tmpDist < minDist) {
       	minDist = tmpDist;
       	indexMinDist = i;
      }
    }
  }
  // Found the nearest tube to A in tubeB.sector and layer
  float dx = Ingrid[indexMinDist].m_x - tubeA.m_x;
  float dy = Ingrid[indexMinDist].m_y - tubeA.m_y;

  // translate A to the tubeB.sector (xy-Plane)
  GridNode TransA(tubeA);// Dit kan beter maar voor nu .... FIXME Later
  TransA.m_x = tubeA.m_x + dx;
  TransA.m_y = tubeA.m_y + dy;

  // Compose the line equations for (X,Y,Z) for each center point and
  // direction vector of the tube.
  TVector3 directionA = TVector3(tubeA.m_WireDirection);// Direction vector of tube A
  TVector3 directionB = TVector3(tubeB.m_WireDirection);// Direction vector of tube B
  /*
   * x = x_a + r(x_v)
   * y = y_a + r(y_v)
   * z = z_a + r(z_v)
   * Use z_v, we know it is never zero.
   */
  float A = tubeB.m_y - TransA.m_y + TransA.m_z * (directionA.Y()/directionA.Z()) - 
            tubeB.m_z * (directionB.Y()/directionB.Z());
  float B = (directionA.Y()/directionA.Z()) - (directionB.Y()/directionB.Z());
  float Z_intersect = A/B;
  float r = (Z_intersect - TransA.m_z)/directionA.Z();
  float s = (Z_intersect - tubeB.m_z)/directionB.Z();
  float X_intersect = tubeB.m_x +  s * directionB.X();
  float Y_intersect = TransA.m_x + r * directionA.X();

  float DistA = sqrt( (TransA.m_x - X_intersect) * (TransA.m_x - X_intersect) +
		      (TransA.m_y - Y_intersect) * (TransA.m_y - Y_intersect) +
		      (TransA.m_z - Z_intersect) * (TransA.m_z - Z_intersect)
		      );
  float DistB = sqrt( (tubeB.m_x - X_intersect) * (tubeB.m_x - X_intersect) +
		      (tubeB.m_y - Y_intersect) * (tubeB.m_y - Y_intersect) +
		      (tubeB.m_z - Z_intersect) * (tubeB.m_z - Z_intersect)
		      );
  
  if( !( DistA > TransA.m_halfLength ) && 
      !( DistB > tubeB.m_halfLength  ) 
      ) {
    // FIXME FIXME Not really optimal (correct??)
    TransA.m_x = (tubeA.m_x + tubeB.m_x)/2.0;
    TransA.m_y = (tubeA.m_y + tubeB.m_y)/2.0;
    TransA.m_z = Z_intersect;
    TransA.m_z_Det = Z_intersect;
    // Set node as virtual and weight becomes 0 (no added value for
    // the length or area size).
    TransA.m_type = GridNode::VIRTUAL_NODE;
    TransA.m_weight = 0;
    TransA.m_SectorLimit = 0;

    out = TransA;
    //// DEBUG DEBUG
    // std::cout << "TubeA: (" << tubeA.m_x << ", " << tubeA.m_y << ", " << tubeA.m_z << ")\n"
    // 	      << "TubeB: (" << tubeB.m_x << ", " << tubeB.m_y << ", " << tubeB.m_z << ")\n"
    // 	      << "Out: ("   << out.m_x << ", " << out.m_y << ", " << out.m_z << ")\n\n";
    ////////////
    return true;
  }// END IF
  return false;
}
//__________________________ END IntersectionPoint_SectorList _______________
//______________________________ GridToNtuple    ____________________________
TNtuple* GridToNtuple(std::vector < GridNode > const &VNodes, std::string const &name)
{
  TNtuple* out = new TNtuple(name.c_str(),"Grid To Ntuple","x:y:det_z:z");
  for(size_t i = 0; i < VNodes.size(); ++i) {
    GridNode const &tube = VNodes[i];
    out->Fill(tube.m_x, tube.m_y, tube.m_z_Det, tube.m_z);
  }
  return out;
}
//_________________________ END GridToNtuple ________________________________
//__________________________ fixNeighboring _________________________________
void fixNeighboring(CoordGrid &hitMap)
{
  std::cerr << "<WARNING> Correcting the missing neighbor relations.\n";
  size_t NumFixed = 0;
  // List of all nodes available in the image (detector map)
  //std::vector< GridNode > &Ingrid = hitMap.GetGrid();
  std::vector< GridNode > &Ingrid = hitMap.m_grid;

  // Tubes loop
  for(size_t i = 0; i < Ingrid.size(); ++i) {
    GridNode &first = Ingrid[i];
    // neighbor tubes loop
    for(size_t j = 0; j < Ingrid.size(); ++j) {
      GridNode &second = Ingrid[j];
      // A neig B and 
      if( ( first.IsNeighboring(second.m_detID) ) &&
	  ( !second.IsNeighboring(first.m_detID) )
	  ) {
	// Add first to the second
	second.m_neighbors.push_back(first.m_detID);
	NumFixed++;
      }
      // Otherwise
      if( ( second.IsNeighboring(first.m_detID) ) &&
          ( !first.IsNeighboring(second.m_detID) ) ){
	first.m_neighbors.push_back(second.m_detID);
	NumFixed++;
      }
    }// END For (j )
  }// END for(i)
  std::cerr << "\t  Total number of corrected neighboring relations = "
	    << NumFixed << '\n';
}
//__________________________ END fixNeighboring _________________________
//_________________________ isolateSectorAndLayerLimits _________________
void isolateSectorAndLayerLimits(CoordGrid const &hitMap, TNtuple &Sections, TNtuple &Layers)
{
  size_t countL, countS;
  countL = countS = 0;

  // Fetch the list of all detectors.
  //std::vector< GridNode > const &Ingrid = hitMap.GetGrid();
  std::vector< GridNode > const &Ingrid = hitMap.m_grid;

  for(size_t i = 0; i < Ingrid.size(); ++i) {
    GridNode const &tube = Ingrid[i];
    
    if( tube.m_LayerLimit &&
	(tube.m_type != GridNode::VIRTUAL_NODE)
      ) {
      Layers.Fill(tube.m_x, tube.m_y, tube.m_z_Det, tube.m_z);
      countL++;
    }
    // it returns -1 if the tube is at the front, 1 if is at the back
    // of the sector border and 0 elsewhere.
    //( (tube.m_SectorLimit == 1) || (tube.m_SectorLimit == -1) || (tube.m_SectorLimit == 0) ) &&
    if( (tube.m_SectorLimit != 0)
	&& (tube.m_type != GridNode::VIRTUAL_NODE)
      ) {
      Sections.Fill(tube.m_x, tube.m_y, tube.m_z_Det, tube.m_z);
      countS++;
    }
  }
  // Set draw attributes for layer limits
  Layers.SetMarkerStyle(8);
  Layers.SetMarkerSize(0.5);
  Layers.SetMarkerColor(kBlack);
  // Set draw attributes for Sector limits
  Sections.SetMarkerStyle(8);
  Sections.SetMarkerSize(0.5);
  Sections.SetMarkerColor(kBlue);

  // Report numbers.
  std::cout << "<DEBUG> CountL = " << countL
	    << " countS = " << countS
	    << '\n';
}
//_________________________ END isolateSectorAndLayerLimits ____________
//________________________ Compute_Virtual_InterSector_Nodes ___________
void Compute_Virtual_InterSector_Nodes(CoordGrid &hitMap, size_t const numSectors,
				       std::vector < GridNode > &VNodes)
{
  // Fetch the list of all STT tubes.
  std::vector< GridNode > &Ingrid = hitMap.m_grid;
 
  // Sector Counting starts from zero.
  // List of per sector tube indices in the map
  std::vector< std::vector<size_t> > SectorLeft;
  std::vector< std::vector<size_t> > SectorRight;

  for(size_t i = 0; i < numSectors; ++i) {
    SectorLeft.push_back(std::vector<size_t>());
    SectorRight.push_back(std::vector<size_t>());
  }
  // DEBUG
  TNtuple limL ("SelectedSectorLeft","SELECTED SKEWED BETWEEN SECTORS","x:y:z");
  TNtuple limR ("SelectedSectorRight","SELECTED SKEWED BETWEEN SECTORS","x:y:z");
  // DEBUG
  // Find sector boundaries per sector and per side
  for(size_t i = 0; i < Ingrid.size(); ++i) {
    GridNode const &tube = Ingrid[i];
    if( (tube.m_type != GridNode::VIRTUAL_NODE) ) {
      // Left or right boundary limit
      if(tube.m_SectorLimit == -1) {
	(SectorLeft[tube.m_Sector]).push_back(i);
	limL.Fill(tube.m_x, tube.m_y, tube.m_z);
      }
      // Right limit
      else if(tube.m_SectorLimit == 1) {
	(SectorRight[tube.m_Sector]).push_back(i);
	limR.Fill(tube.m_x, tube.m_y, tube.m_z);
      }
    }
  }
  //_______
  std::cout << "<DEBUG> Total number of sectors = " << SectorLeft.size()
	    << "\n\t<-I-> Number of left boundary limits = " << SectorLeft.size()
	    << "\n\t<-I-> Number of right boundary limits = " << SectorRight.size() 
	    << '\n';
  /*
   *   0 / \ 5
   *   1 | | 4
   *   2 \ / 3
   */
  std::vector < std::pair<size_t, size_t> > InterSectorPairs;
  float minDist = std::numeric_limits<float>::max();
  float currDist = 0;
  // Index of the tube with the shortest dstance
  size_t shortestIndex = std::numeric_limits<size_t>::max();
  size_t R_Index = 0;
  // Loop all the sectors.
  for(size_t s = 0; s < numSectors; ++s) {
    std::vector<size_t> const &currentLeftSector = SectorLeft[s];
    // The counter sector is SectorRight[s - 1]. The right limts of
    // the previous sector.
    R_Index = ( s == 0 ) ? (SectorRight.size() - 1) : (s - 1);
    std::vector<size_t> const &currentRightSector = SectorRight[R_Index];
    // For each tube find its nearest neighbor in the other sector.
    for(size_t l = 0; l < currentLeftSector.size(); ++l) {
      GridNode const &CurLftNode = Ingrid[currentLeftSector[l]];
      minDist = std::numeric_limits<float>::max();
      for(size_t r = 0; r < currentRightSector.size(); ++r) {
	GridNode const &CurRgtNode = Ingrid[currentRightSector[r]];
	// If in the same layer
	if( (CurLftNode.m_Layer == CurRgtNode.m_Layer) &&
	    (CurLftNode.m_Sector != CurRgtNode.m_Sector) &&
	    (CurLftNode.m_type == GridNode::STT_TYPE_SKEW && CurRgtNode.m_type == GridNode::STT_TYPE_SKEW) &&
	    (CurLftNode.m_detID != CurRgtNode.m_detID)
	  ){
	  currDist = (CurLftNode.m_x - CurRgtNode.m_x) * (CurLftNode.m_x - CurRgtNode.m_x) +
	             (CurLftNode.m_y - CurRgtNode.m_y) * (CurLftNode.m_y - CurRgtNode.m_y) +
	             (CurLftNode.m_z - CurRgtNode.m_z) * (CurLftNode.m_z - CurRgtNode.m_z);
	  // Update shortest distance
	  if(minDist > currDist) {
	    minDist = currDist;
	    shortestIndex = currentRightSector[r];
	  }
	  //limL.Fill(CurLftNode.m_x, CurLftNode.m_y, CurLftNode.m_z);
	  //limR.Fill(CurRgtNode.m_x, CurRgtNode.m_y, CurRgtNode.m_z);
	}// If same layer And not same sector
      }// Right limit Loop
      // Found the tube with the shortest distance.
      InterSectorPairs.push_back( std::make_pair(currentLeftSector[l], shortestIndex) );
    }// Left limit loop
  }// Sectors loop
  // Now we have all the tube pairs for which we can add virtual nodes
  // between the sectors.
  // DEBUG NTUPLE
  limL.SetMarkerStyle(8);
  limL.SetMarkerSize(0.7);
  limL.SetMarkerColor( kRed + 2);
  limL.Write();

  limR.SetMarkerStyle(8);
  limR.SetMarkerSize(0.7);
  limR.SetMarkerColor( kCyan + 2);
  limR.Write();
  // END NTUPLE DEBUG
  //________________________________ DEBUG _________________________________________
  std::cout << "<DEBUG> total number of pairs = " << InterSectorPairs.size() << '\n';
  //____________________________________ DEBUG _____________________________________

  for(size_t p = 0; p < InterSectorPairs.size(); ++p) {
    std::pair<size_t, size_t> const& indexPair = InterSectorPairs[p];
    GridNode out;
    GridNode &first_tube  = Ingrid[indexPair.first];
    GridNode &second_tube = Ingrid[indexPair.second];
    InterSectorPoints(first_tube, second_tube, out);

    out.m_detID = 50000 + p;
    out.m_type = GridNode::VIRTUAL_NODE;
    out.m_weight = 0;
    out.m_neighbors.clear();
    (out.m_neighbors).push_back(first_tube.m_detID);
    (out.m_neighbors).push_back(second_tube.m_detID);
    out.m_SectorLimit = 0;
    
    // (first_tube.m_neighbors).push_back(out.m_detID);
    // (second_tube.m_neighbors).push_back(out.m_detID);
    
    VNodes.push_back(out);
  }
}
//______________________ END Compute_Virtual_InterSector_Nodes ___________
//__________________________ InterSectorPoints ___________________________
/**
 */
void InterSectorPoints(GridNode const& tubeA, GridNode const& tubeB, GridNode &out)
{
  /*
   * The line equations are of the form:
   * x = x_a + r * x_v
   * y = y_a + r * y_v
   * z = z_a + r * z_v
   * Distance(D) is half length in each direcion.
   * r = D/sqrt(x_v^2 + y_v^2 + z_v^2)
   * For the other direction just reverse the direction of
   * the vector V.
   */
  // Find the end points for both input tubes.
  TVector3 dirA = tubeA.m_WireDirection;
  float R_A =  tubeA.m_halfLength / sqrt( (dirA[0]*dirA[0]) + (dirA[1]*dirA[1]) + (dirA[2]*dirA[2]) );
  TVector3 dirB = tubeB.m_WireDirection;
  float R_B =  tubeB.m_halfLength / sqrt( (dirB[0]*dirB[0]) + (dirB[1]*dirB[1]) + (dirB[2]*dirB[2]) );
  // Determine the end points
  point3D A_left;
  A_left.m_x = tubeA.m_x - R_A * dirA[0];
  A_left.m_y = tubeA.m_y - R_A * dirA[1];
  A_left.m_z = tubeA.m_z - R_A * dirA[2];
  point3D A_right;
  A_right.m_x = tubeA.m_x + R_A * dirA[0];
  A_right.m_y = tubeA.m_y + R_A * dirA[1];
  A_right.m_z = tubeA.m_z + R_A * dirA[2];
  point3D B_left;
  B_left.m_x = tubeB.m_x - R_B * dirB[0];
  B_left.m_y = tubeB.m_y - R_B * dirB[1];
  B_left.m_z = tubeB.m_z - R_B * dirB[2];
  point3D B_right;
  B_right.m_x = tubeB.m_x + R_B * dirB[0];
  B_right.m_y = tubeB.m_y + R_B * dirB[1];
  B_right.m_z = tubeB.m_z + R_B * dirB[2];
  // Compute shortest distance
  float Dist_LL = PointDistance(A_left,  B_left);
  float Dist_LR = PointDistance(A_left,  B_right);
  float Dist_RL = PointDistance(A_right, B_left);
  float Dist_RR = PointDistance(A_right, B_right);
  // Construct the new virtual tube.
  if( (Dist_LL < Dist_LR) && (Dist_LL < Dist_RL) && ( Dist_LL < Dist_RR)
    ) {
    out.m_x = (A_left.m_x + B_left.m_x) / 2.00;
    out.m_y = (A_left.m_y + B_left.m_y) / 2.00;
    out.m_z = (A_left.m_z + B_left.m_z) / 2.00;
  }
  else if( ( Dist_LR < Dist_LL ) && (Dist_LR < Dist_RL) && ( Dist_LR < Dist_RR)
	 ) {
    out.m_x = (A_left.m_x + B_right.m_x) / 2.00;
    out.m_y = (A_left.m_y + B_right.m_y) / 2.00;
    out.m_z = (A_left.m_z + B_right.m_z) / 2.00;
  }
  else if( (Dist_RL < Dist_LL) && (Dist_RL < Dist_LR) && (Dist_RL < Dist_RR)
	 ) {
    out.m_x = (A_right.m_x + B_left.m_x) / 2.00;
    out.m_y = (A_right.m_y + B_left.m_y) / 2.00;
    out.m_z = (A_right.m_z + B_left.m_z) / 2.00;
  }
  else if( (Dist_RR < Dist_LL) && (Dist_RR < Dist_LR) && (Dist_RR < Dist_RL)
	 ) {
    out.m_x = (A_right.m_x + B_right.m_x) / 2.00;
    out.m_y = (A_right.m_y + B_right.m_y) / 2.00;
    out.m_z = (A_right.m_z + B_right.m_z) / 2.00;
  }
}
//__________________________ END InterSectorPoints _______________________
//________________________ PointDistance _________________________________
float PointDistance(point3D const &a, point3D const &b)
{
  return sqrt( (a.m_x - b.m_x) * (a.m_x - b.m_x) +
	       (a.m_y - b.m_y) * (a.m_y - b.m_y) +
	       (a.m_z - b.m_z) * (a.m_z - b.m_z)
	      );
}
//________________________ END PointDistance _____________________________
//____________________________ CircleFit _________________________________
int CircleFit(CoordGrid const &hitMap, std::vector<int> const &MemberVector, CurvatureParameters &curvature)
{
  /*
   * Least squares fit of circle to set of points. This is not
   * efficient enough but for now its OK. We need to find an additive
   * algorithm.
   *
   Input:  (x_i,y_i), 1 <= i <= N, where N >= 3 and not all points
   are collinear
   Output:  circle center (a,b) and radius r
   
   Energy function to be minimized is
   
   E(a,b,r) = sum_{i=1}^N (L_i-r)^2
   
   where L_i = |(x_i-a,y_i-b)|, the length of the specified vector.
   Taking partial derivatives and setting equal to zero yield the
   three nonlinear equations
   
   E_r = 0:  r = Average(L_i)
   E_a = 0:  a = Average(x_i) + r * Average(dL_i/da)
   E_b = 0:  b = Average(y_i) + r * Average(dL_i/db)
   
   Replacing r in the last two equations yields
   
   a = Average(x_i) + Average(L_i) * Average(dL_i/da) = F(a,b)
   b = Average(y_i) + Average(L_i) * Average(dL_i/db) = G(a,b)
   
   which can possibly be solved by fixed point iteration as
   
   a_{n+1} = F(a_n,b_n),  b_{n+a} = G(a_n,b_n)
   
   with initial guess a_0 = Average(x_i) and b_0 = Average(y_i).
   Derivative calculations show that
   
   dL_i/da = (a-x_i)/L_i,  dL_i/db = (b-y_i)/L_i.
  */
  // std::cout << "<INFO> Determine curvature using circle least squares fit.\n";
  // Fetch the entire grid
  std::vector< GridNode > const &Ingrid = hitMap.m_grid;
   /* Local constants */
   //256;//1000, 10000;//512;//256;
  const size_t maxIterations = MAX_NUMBER_OF_FIT_ITERATION;
  //1e-04;//1e-06;
  const double tolerance = LOCAL_CIRCLE_FIT_TOLERANCE;
  double a, b, r;

  /* compute the average of the data points */
  double xAvr = 0.0;
  double yAvr = 0.0;
  for(size_t k = 0; k < MemberVector.size(); ++k) {
    int mem_Id = MemberVector[k];
    int mem_Idx = hitMap.Find(mem_Id);
    GridNode const &node = Ingrid[mem_Idx];
    xAvr += node.m_xDet;
    yAvr += node.m_yDet;
  }
  xAvr /= static_cast<double> (MemberVector.size());
  yAvr /= static_cast<double> (MemberVector.size());

  // Initial guess
  a = xAvr;
  b = yAvr;

  size_t stepCounter = 0;
  
  for (size_t j = 0; j < maxIterations; j++) {
    /* update the iterates */
    double a0 = a;
    double b0 = b;

    /* compute average L, dL/da, dL/db */
    double LAvr  = 0.0;
    double LaAvr = 0.0;
    double LbAvr = 0.0;

    for (size_t i = 0; i < MemberVector.size(); i++) {
      int mem_Id  = MemberVector[i];
      int mem_Idx = hitMap.Find(mem_Id);
      GridNode const &node = Ingrid[mem_Idx];
      double dx = node.m_xDet - a;
      double dy = node.m_yDet - b;
      double L = sqrt(dx * dx + dy * dy);
      if (fabs(L) > tolerance) {
        LAvr  += L;
        LaAvr -= dx / L;
        LbAvr -= dy / L;
      }
    }
    LAvr  /= static_cast<double>(MemberVector.size());
    LaAvr /= static_cast<double>(MemberVector.size());
    LbAvr /= static_cast<double>(MemberVector.size());

    a = xAvr + LAvr * LaAvr;
    b = yAvr + LAvr * LbAvr;
    r = LAvr;
    // Exit, if (a,b) is not changing a lot
    if( (fabs(a - a0) <= tolerance) && (fabs(b - b0) <= tolerance) ){
      break;
    }
    stepCounter = j;
  }

  // Now we need to determine the value for E(a,b,r)
  double Li   = 0;
  double Eabr = 0;
  for(size_t i = 0; i < MemberVector.size(); ++i){
    int node_Id = MemberVector[i];
    int node_Idx = hitMap.Find(node_Id);
    GridNode const &node = Ingrid[node_Idx];
    // Li = sqrt( ((node.m_xDet - a) * (node.m_xDet - a)) + ((node.m_yDet - b) * (node.m_yDet - b)) );
    Li = sqrt( ((node.m_x - a) * (node.m_x - a)) + ((node.m_y - b) * (node.m_y - b)) );
    Eabr += ( (Li - r) * (Li - r) );
  }
  // Prepare the output
  curvature.m_a = a;
  curvature.m_b = b;
  curvature.m_ra = r;
  curvature.m_r = (1.00/ r);
  curvature.m_E = Eabr;

  return ( (stepCounter < maxIterations) ? stepCounter : -1);
}
//________________________ END CircleFit _________________________________
//_____________________ CIrcleFit Corrected for MC Data ___________
int CircleFit( std::vector<point3D> const &points, CurvatureParameters &curvature) {
  /* Local constants */
  const size_t maxIterations = MAX_NUMBER_OF_FIT_ITERATION;
  //1e-04;//1e-06;
  const double tolerance = LOCAL_CIRCLE_FIT_TOLERANCE;
  double a, b, r;

  /* compute the average of the data points */
  double xAvr = 0.0;
  double yAvr = 0.0;
  for(size_t k = 0; k < points.size(); ++k) {
    point3D const &curPoint = points[k];
    xAvr += curPoint.m_x;
    yAvr += curPoint.m_y;
  }
  xAvr /= static_cast<double> (points.size());
  yAvr /= static_cast<double> (points.size());

  // Initial guess
  a = xAvr;
  b = yAvr;

  size_t stepCounter = 0;
  
  for (size_t j = 0; j < maxIterations; j++) {
    /* update the iterates */
    double a0 = a;
    double b0 = b;

    /* compute average L, dL/da, dL/db */
    double LAvr  = 0.0;
    double LaAvr = 0.0;
    double LbAvr = 0.0;

    for (size_t i = 0; i < points.size(); i++) {
      point3D const &pt  = points[i];
      double dx = pt.m_x - a;
      double dy = pt.m_y - b;
      double L = sqrt(dx * dx + dy * dy);
      if (fabs(L) > tolerance) {
        LAvr  += L;
        LaAvr -= dx / L;
        LbAvr -= dy / L;
      }
    }
    LAvr  /= static_cast<double>(points.size());
    LaAvr /= static_cast<double>(points.size());
    LbAvr /= static_cast<double>(points.size());

    a = xAvr + LAvr * LaAvr;
    b = yAvr + LAvr * LbAvr;
    r = LAvr;
    // Exit, if (a,b) is not changing a lot
    if( (fabs(a - a0) <= tolerance) && (fabs(b - b0) <= tolerance) ){
      break;
    }
    stepCounter = j;
  }

  // Now we need to determine the value for E(a,b,r)
  double Li   = 0;
  double Eabr = 0;
  for(size_t i = 0; i < points.size(); ++i){
    point3D const &po = points[i];
    // Li = sqrt( ((node.m_xDet - a) * (node.m_xDet - a)) + ((node.m_yDet - b) * (node.m_yDet - b)) );
    Li = sqrt( ((po.m_x - a) * (po.m_x - a)) + ((po.m_y - b) * (po.m_y - b)) );
    Eabr += ( (Li - r) * (Li - r) );
  }
  // Prepare the output
  curvature.m_a = a;
  curvature.m_b = b;
  curvature.m_ra = r;
  curvature.m_r = (1.00/ r);
  curvature.m_E = Eabr;

  return ( (stepCounter < maxIterations) ? stepCounter : -1);
}
//_____________________ END CircleFit Corrected for MC data _______

//________________ updateCandidateNodesFitVal ____________________________
void updateCandidateNodesFitVal (CoordGrid  &hitMap, PathCandidate const &Candidate)
{
  // Set of node id's of all member nodes.
  std::set<int> *SetOfMembers = Candidate.m_memberIdSet;
  double a = (Candidate.m_CurV_par).m_a;
  double b = (Candidate.m_CurV_par).m_b;
  double Li = 0.00;
  std::set<int>::iterator it;
  for( it = SetOfMembers->begin(); it != SetOfMembers->end(); ++it) {
    int node_Id = *it;
    GridNode &node = hitMap.GetNodeByID(node_Id);
    // Li = sqrt( ((node.m_xDet - a) * (node.m_xDet - a)) + ((node.m_yDet - b) * (node.m_yDet - b)) );
    Li = sqrt( ((node.m_x - a) * (node.m_x - a)) + ((node.m_y - b) * (node.m_y - b)) );
    node.m_fitValue = Li;
  }
}
//________________________________________________________________________
//_____________________________ updateNodeFitValForAllPaths ______________
void updateNodeFitValForAllPaths(CoordGrid &hitMap, std::vector<PathCandidate*> const &Candidates)
{
  for(size_t i = 0; i < Candidates.size(); ++i) {
    PathCandidate const *candidate = Candidates[i];
    updateCandidateNodesFitVal(hitMap,(*candidate));
  }
}
//________________________________________________________________________
//___________________________ ComputeCurvatureForListOfnodes _______________
void ComputeCurvatureForListOfnodes(CoordGrid const &hitMap, std::vector<int> const &nodes, CurvatureParameters &OutPars)
{
  std::vector< GridNode > const &Ingrid = hitMap.m_grid;
  if(nodes.size() < 3) {
    std::cout << "<INFO> List has " << nodes.size()
              << " nodes. too few members. Skip.\n";
    return;
  }
  CurvatureParameters CUVPars;
  int cf_out = CircleFit(hitMap, nodes, CUVPars);
  if( cf_out == -1) {
    std::cerr << "<ERROR> Finding the best circle fit failed.\n"
              << "\t object will not be updated.\n";
  }
  else{
    OutPars.m_a = CUVPars.m_a;
    OutPars.m_b = CUVPars.m_b;
    OutPars.m_ra = CUVPars.m_ra;
    OutPars.m_r = CUVPars.m_r;
    OutPars.m_E = CUVPars.m_E;
  }
}
//____________________________ computePathCurvature ____________________
void computePathCurvature(CoordGrid const &hitMap, PathCandidate &Candidate)
{
  std::cout << "<INFO> Determining the curvature.\n";
  // Fetch the entire grid
  std::vector< GridNode > const &Ingrid = hitMap.m_grid;
  // Set of node id's of all member nodes.
  std::set<int> *SetOfMembers = Candidate.m_memberIdSet;
  if( SetOfMembers->size() < 3 ) {
    std::cout << "<INFO> Candidate contains " << SetOfMembers->size()
              << " nodes. too few members. Skip.\n";
    // FIXME kan dit???
    (Candidate.m_CurV_par).m_a = 1000000;
    (Candidate.m_CurV_par).m_b = 1000000;
    (Candidate.m_CurV_par).m_ra = 1000000;
    (Candidate.m_CurV_par).m_r = 1000000;
    (Candidate.m_CurV_par).m_E = 1000000;
    Candidate.m_hasCurvature   = false;
    return;
  }
  // Not efficient but for now OK
  std::vector<int> MemberVector(SetOfMembers->begin(), SetOfMembers->end());
  //DEBUG INFO  
  std::cout << "<INFO> Candidate contains " << MemberVector.size()
            << " Handeling track candidate.\n";
  // DEBUG INFO
  CurvatureParameters curvaturePars;
  int cf_out = CircleFit(hitMap, MemberVector, curvaturePars);
  if( cf_out == -1 ) {
    std::cerr << "<ERROR> Finding the best circle fit failed.\n"
              << "\t object will not be updated.\n";
    Candidate.m_hasCurvature   = false;
  }
  else{
    (Candidate.m_CurV_par).m_a = curvaturePars.m_a;
    (Candidate.m_CurV_par).m_b = curvaturePars.m_b;
    (Candidate.m_CurV_par).m_ra = curvaturePars.m_ra;
    (Candidate.m_CurV_par).m_r = curvaturePars.m_r;
    (Candidate.m_CurV_par).m_E = curvaturePars.m_E;
    Candidate.m_hasCurvature   = true;
  }
}
//________________________ END computePathCurvature ____________________
//____________________________ computeCurvatureForAllPaths _______________
void computeCurvatureForAllPaths(CoordGrid const &hitMap, std::vector<PathCandidate*> &Candidates)
{
  for(size_t i = 0; i < Candidates.size(); ++i){
    PathCandidate *Cand = Candidates[i];
    if( Cand->m_isValid){
      // Update curvature data.
      computePathCurvature(hitMap, (*Cand));
    }
  }
}
//________________________ END computeCurvatureForAllPaths _______________
//____________________________ cleanupDuplicateSingletons ________________
void cleanupDuplicateSingletons(std::vector<PathCandidate*> &Candidates)
{
  std::set<int>::iterator findIt;
  for(size_t i = 0; i < Candidates.size(); ++i) {
    PathCandidate *singletone = Candidates[i];
    std::set<int> *memberList = singletone->m_memberIdSet;
    unsigned int singletonID  = singletone->m_id;
    if( (memberList->size() == 1) && (singletone->m_isValid) ){// It is a singletone
      // If it appears in any other list
      for(size_t j = 0; j < Candidates.size(); ++j) {
        PathCandidate *otherCand = Candidates[j];
        if( otherCand->m_id != singletonID ) {// Not the same tracklet
          std::set<int> *otherMembers = otherCand->m_memberIdSet;
          findIt = std::find(otherMembers->begin(), otherMembers->end(), singletonID);
          if( findIt != otherMembers->end()) {// found in the list
            singletone->m_isValid = false;
          }
        }// Not the same candidate
      }// END for J
    }// Sigleton and valid
  }// For I
}
//________________________ END cleanupDuplicateSingletons ________________
//____________________________ mergeSingletons ___________________________
void mergeSingletons( CoordGrid const &hitMap, std::set<int> &SingleOfNodesID,
                      std::vector<PathCandidate*> &Candidates)
{
  std::cout << "\t<-I-> Merge Singletons. Num input ids = " << SingleOfNodesID.size()
            << " Num input candidates = " << Candidates.size() << '\n';
  // Fetch the entire grid
  std::vector< GridNode > const &Ingrid = hitMap.m_grid;
  // Local variables
  std::set<int>::iterator FindIt;
  std::set<int>::iterator elementIt;
  std::vector<int> ListOfNotMerged;
  float first_dist, second_dist, tmpDist;
  first_dist = second_dist = tmpDist = 0;
  size_t bestMatchCandIndex = 0;
  float  shortestDist = std::numeric_limits<float>::max();
  int    nodeId_ShortestDist = -1;
  std::vector<int>::iterator neighbourIt;
  /* Nodes are sorted in increasing order, from inner most layer to
     outer most layer, spiral form. */
  std::vector<int> singletons(SingleOfNodesID.begin(), SingleOfNodesID.end());
  std::sort(singletons.begin(), singletons.end());
  //bool matched = false;
  // Loop singleton node ID list
  while (!singletons.empty()) {
    int Singlenode_ID = singletons.back();
    singletons.pop_back();
    GridNode const &S_node = hitMap.GetNodeByID(Singlenode_ID);
    std::vector<int> ListOfActiveNeighbours;
    // Get the list of first and second-order active neighbours.
    hitMap.ActiveNeighboursAndSecondOrderNeighbours(Singlenode_ID, ListOfActiveNeighbours);
    
    for(size_t i = 0; i < Candidates.size(); ++i) {
      PathCandidate *cand = Candidates[i];
      std::set<int> *CandMembers = cand->m_memberIdSet;
      FindIt = std::find(CandMembers->begin(), CandMembers->end(), Singlenode_ID);
      // If not added before, check if it fits to the current
      // candidate
      if( FindIt == CandMembers->end() ){
        elementIt = CandMembers->begin();
        int firstNodeID = *elementIt;
        GridNode const &First_CandNode = hitMap.GetNodeByID(firstNodeID);
        // first_dist = sqrt( (S_node.m_xDet - First_CandNode.m_xDet) * (S_node.m_xDet - First_CandNode.m_xDet) +
        //                    (S_node.m_yDet - First_CandNode.m_yDet) * (S_node.m_yDet - First_CandNode.m_yDet) );
        first_dist = sqrt( (S_node.m_x - First_CandNode.m_x) * (S_node.m_x - First_CandNode.m_x) +
                           (S_node.m_y - First_CandNode.m_y) * (S_node.m_y - First_CandNode.m_y) );
        elementIt = CandMembers->end();
        elementIt--;
        int lastNodeID = *elementIt;
        // To avoid virtuals.
        while(lastNodeID > 5000) {
          elementIt--;
          lastNodeID = *elementIt;
        }
        GridNode const &Second_CandNode = hitMap.GetNodeByID(lastNodeID);
        // second_dist = sqrt( (S_node.m_xDet - Second_CandNode.m_xDet) * (S_node.m_xDet - Second_CandNode.m_xDet) +
        //                     (S_node.m_yDet - Second_CandNode.m_yDet) * (S_node.m_yDet - Second_CandNode.m_yDet) );
        second_dist = sqrt( (S_node.m_x - Second_CandNode.m_x) * (S_node.m_x - Second_CandNode.m_x) +
                            (S_node.m_y - Second_CandNode.m_y) * (S_node.m_y - Second_CandNode.m_y) );

        tmpDist = ( (first_dist < second_dist) ? first_dist : second_dist);
        if(tmpDist < shortestDist){
          shortestDist = tmpDist;
          bestMatchCandIndex = i;
          nodeId_ShortestDist = ( (first_dist < second_dist) ? firstNodeID : lastNodeID);
        }
      }
    }// END candidate list loop
    // Candidates with the minimum Euclidean distance is found. Check
    // if we can merge.
    PathCandidate *candToMergeWith  = Candidates[bestMatchCandIndex];
    std::set<int> *BestMatchMembers = candToMergeWith->m_memberIdSet;
    // If the node with the shortest distance in the neighbour list.
    neighbourIt = std::find(ListOfActiveNeighbours.begin(), ListOfActiveNeighbours.end(), nodeId_ShortestDist);
    // If one of the neighbours or a second-order neighbour in the
    // list
    if( neighbourIt != ListOfActiveNeighbours.end() ) {
      BestMatchMembers->insert(Singlenode_ID);
    }
    else{// Could not merge
      /*FIXME: Dit klopt niet helemaal. Op deze manier wordt de node
       * aan een van de paden toegevoegd, zelfs als de afstand te
       * groot is en geen van zijn buren of over(over) buren in een
       * pad voorkomen.  Dit moet nog verbeterd worden.*/
      ListOfNotMerged.push_back(Singlenode_ID);

      //BestMatchMembers->insert(Singlenode_ID);
    }
    // Reset, go to next singleton
    bestMatchCandIndex = 0;
    shortestDist = std::numeric_limits<float>::max();
  }// END While not empty sigletons
  // Clean input list and replace with the not merged node ids
  SingleOfNodesID.clear();
  SingleOfNodesID.insert(ListOfNotMerged.begin(), ListOfNotMerged.end());
}
//________________________ END mergeSingletons ___________________________
//____________________________ mergeShortPaths ___________________________
void mergeShortPaths(CoordGrid const &hitMap, std::vector < PathCandidate* > &CandidesContainer, size_t minL)
{
  /* Local variables */
  float distHeadHead, distHeadTail, distTailHead, distTailTail;
  distHeadHead = distHeadTail = distTailHead = distTailTail = 0.0;
  // Find possible merge candidates for all short tracks.
  for(size_t i = 0; i < CandidesContainer.size(); ++i) {
    PathCandidate *Shortcandidate = CandidesContainer[i];
    if(Shortcandidate->m_isShort){
      // List of distances of the current tracklet to all other
      // tracklets.
      std::vector<candidateDistobject> MergePairs;
      GridNode const &F_head = hitMap.GetNodeByID(Shortcandidate->m_maxlayerNodeId);
      GridNode const &F_tail = hitMap.GetNodeByID(Shortcandidate->m_minlayerNodeId);
      for(size_t j = 0; j < CandidesContainer.size(); ++j) {
        if(j!= i){
          PathCandidate *mergeCand = CandidesContainer[j];
          GridNode const &S_head = hitMap.GetNodeByID(mergeCand->m_maxlayerNodeId);
          GridNode const &S_tail = hitMap.GetNodeByID(mergeCand->m_minlayerNodeId);
          // Compute 2D distances
          distHeadHead = F_head.GetDistance2D(S_head);
          distHeadTail = F_head.GetDistance2D(S_tail);
          distTailHead = F_tail.GetDistance2D(S_head);
          distTailTail = F_tail.GetDistance2D(S_tail);
          candidateDistobject candDistObject;
          float minDist = std::numeric_limits<float>::max();
          // Determine minimum distance
          if( distHeadHead < minDist) {
            minDist = distHeadHead;
            candDistObject.distance = minDist;
            candDistObject.FirstCandNodeId  = Shortcandidate->m_maxlayerNodeId;
            candDistObject.SecondCandNodeId = mergeCand->m_maxlayerNodeId;
            candDistObject.indexFirstCand  = i;
            candDistObject.indexSecondCand = j;
          }
          if( distHeadTail < minDist) {
            minDist = distHeadTail;
            candDistObject.distance = minDist;
            candDistObject.FirstCandNodeId  = Shortcandidate->m_maxlayerNodeId;
            candDistObject.SecondCandNodeId = mergeCand->m_minlayerNodeId;
            candDistObject.indexFirstCand  = i;
            candDistObject.indexSecondCand = j;
          }
          if( distTailHead < minDist) {
            minDist = distTailHead;
            candDistObject.distance = minDist;
            candDistObject.FirstCandNodeId  = Shortcandidate->m_minlayerNodeId;
            candDistObject.SecondCandNodeId = mergeCand->m_maxlayerNodeId;
            candDistObject.indexFirstCand  = i;
            candDistObject.indexSecondCand = j;
          }
          if( distTailTail < minDist){
            minDist = distTailTail;
            candDistObject.distance = minDist;
            candDistObject.FirstCandNodeId  = Shortcandidate->m_minlayerNodeId;
            candDistObject.SecondCandNodeId = mergeCand->m_minlayerNodeId;
            candDistObject.indexFirstCand  = i;
            candDistObject.indexSecondCand = j;
          }
          MergePairs.push_back(candDistObject);          
        }// End if( j != i)
      }// Other candidates loop (possible merge Candidates)
      // Sort distances
      std::sort(MergePairs.begin(), MergePairs.end());
      // Try to merge
      for(size_t k = 0; k < MergePairs.size(); ++k) {
        candidateDistobject const &dist_obj = MergePairs[k];
        // Dit moet nog beter FIXME FIXME
        if( dist_obj.distance < SHORT_CANDIDATES_MERGE_DISTANCE ) {
          PathCandidate *First  = CandidesContainer[dist_obj.indexFirstCand];
          std::set<int> *FirstMembersSet = First->m_memberIdSet;
          PathCandidate *Second = CandidesContainer[dist_obj.indexSecondCand];
          std::set<int> *SecondMembersSet = Second->m_memberIdSet;
          if( First->m_isValid && Second->m_isValid) {
            SecondMembersSet->insert(FirstMembersSet->begin(), FirstMembersSet->end());
            Second->m_isShort = (SecondMembersSet->size() < minL);
            FirstMembersSet->insert(SecondMembersSet->begin(), SecondMembersSet->end());
            First->m_isShort = ( FirstMembersSet->size() < minL);
          }
        }
      }
    }// If current is short
  }// First candidates loop
}
//____________________________ END mergeShortPaths _______________________
//___________________________ updatePathCandidateHeadAndTailNodes ________
void updatePathCandidateHeadAndTailNodes(PathCandidate &cand)
{
  std::set<int> *Members = cand.m_memberIdSet;
  std::set<int>::iterator it = Members->begin();
  if(Members->size() == 1){// Singleton
    cand.m_minlayerNodeId = *it;
    cand.m_maxlayerNodeId = *it;
  }
  else{
    cand.m_minlayerNodeId = *it;
    it = Members->end();
    --it;
    cand.m_maxlayerNodeId = *it;
    // Skip virtuals. Virtuals don't have well defined layer numbers.
    while( (cand.m_maxlayerNodeId >= START_VIRTUAL_ID) && (it != Members->begin()) ) {
      //--it;
      //cand.m_maxlayerNodeId = *it;
      cand.m_maxlayerNodeId = *(--it);
    }
  }
}
//_______________________ END updatePathCandidateHeadAndTailNodes ________
void updateAllCandidateHeadAndTailNodes(std::vector<PathCandidate*> &Candidates)
{
  for(size_t i = 0; i < Candidates.size(); ++i){
    PathCandidate *Cand = Candidates[i];
    if( Cand->m_isValid){
      // Update head and tail data.
      updatePathCandidateHeadAndTailNodes(*Cand);
    }
  }
}
//_______________________ Path validation and remove invalid paths _________
bool isInValidPath(PathCandidate const *cand)
{
  return (!cand->m_isValid);
}

void removeInvalidSubPaths(std::vector < PathCandidate* > &CandidesContainer)
{
  std::cout << "<INFO> Input has "<< CandidesContainer.size();
  CandidesContainer.erase( std::remove_if(CandidesContainer.begin(), CandidesContainer.end(), isInValidPath),
                           CandidesContainer.end());
  std::cout << " After removing, container has " << CandidesContainer.size()
            << " tracklets.\n";
}
//___________________ END Path validation and remove invalid paths _________
//_______________________ pathNumOfEndStops _______________
size_t pathNumOfEndStops(CoordGrid const &hitMap, std::set<int> const &nodeList)
{
  /*The list is sorted(ascending), virtuals have large ids; this way the appear
    at the end of the list.*/
  /* FIXME FIXME Dit klopt natuurlijk niet.*/
  size_t numEnds   = 10;
  return numEnds;
}
//___________________ END pathNumOfEndStops _______________
//_______________________ candidateNumOfEndStops __________
size_t candidateNumOfEndStops(CoordGrid const &hitMap, PathCandidate const &inputCand)
{
  std::set<int> const *nodeList = inputCand.m_memberIdSet;
  size_t numEnds = pathNumOfEndStops(hitMap, (*nodeList));
  return numEnds;
}
//___________________ END candidateNumOfEndStops __________
//___________________________ numEndStopsIfMerged __________________________
size_t numEndStopsIfMerged(CoordGrid const &hitMap, PathCandidate const &p1, PathCandidate const &p2)
{
  size_t numNstops = 3;
  std::set<int> const *firstMembers  = p1.m_memberIdSet;
  std::set<int> const *secondMembers = p2.m_memberIdSet;
  int min1_id = *std::min_element(firstMembers->begin(), firstMembers->end());
  int max1_id = *std::max_element(firstMembers->begin(), firstMembers->end());
  int min2_id = *std::min_element(secondMembers->begin(), secondMembers->end());
  int max2_id = *std::max_element(secondMembers->begin(), secondMembers->end());
  GridNode const &firstMin  = hitMap.GetNodeByID(p1.m_minlayerNodeId);
  GridNode const &firstMax  = hitMap.GetNodeByID(p1.m_maxlayerNodeId);
  GridNode const &secondMin = hitMap.GetNodeByID(p2.m_minlayerNodeId);
  GridNode const &secondMax = hitMap.GetNodeByID(p2.m_maxlayerNodeId);
  /* If one of the candidates is a singleton or a short one*/
  if( (firstMembers->size() < 3) || (secondMembers->size() < 3) ){
    std::cout << "Length less than 3. F is " << firstMembers->size()
              << " S = " << secondMembers->size() << '\n';
    numNstops = 2;
  }
  else {// Not a singletone (determine intersection)
    std::vector<int> *intersect = findIntersectionElements(p1, p2);
    if( intersect->size() == 0) {
      if( (firstMin.m_Layer >= secondMax.m_Layer) || (firstMax.m_Layer <= secondMin.m_Layer) ) {
        std::cout << "intersection 0. G 1 & 2.\n";
        numNstops = 2;
      }
      else if( (firstMin.m_Layer <= secondMin.m_Layer) && (firstMax.m_Layer >= secondMin.m_Layer) ){
        std::cout << "intersection 0. G 6.\n";
        numNstops = 3;
      }
      else if( (firstMin.m_Layer <= secondMax.m_Layer) && (firstMin.m_Layer > secondMin.m_Layer) ) {
        int maxlayerDiff = ( (firstMax.m_Layer > secondMax.m_Layer) ? (firstMax.m_Layer - secondMax.m_Layer) :
                             (secondMax.m_Layer - firstMax.m_Layer) );
        int minlayerDiff = ( (firstMin.m_Layer > secondMin.m_Layer) ? (firstMin.m_Layer - secondMin.m_Layer) :
                             (secondMin.m_Layer - firstMin.m_Layer) );
        if( (maxlayerDiff < 4) && (minlayerDiff < 4) ) {
        numNstops = 3;// FIXME is this correct?
        }
        else{
          numNstops = 2;// FIXME is this correct?
        }
        std::cout << "Haha raar 2\n";
      }
      else {
        numNstops = 3;
        std::cout << "\nfirstMin.m_Layer = " << firstMin.m_Layer
                  << " firstMax.m_Layer = " << firstMax.m_Layer
                  << " secondMin.m_Layer = " << secondMin.m_Layer
                  << " secondMax.m_Layer = " << secondMax.m_Layer
                  << '\n';
      }
    }// End no intersection
    else if(intersect->size() == 1) {// One of the three possiblilities
      if( (min1_id == max2_id) || (max1_id == min2_id) ){
        numNstops = 2;
        std::cout << "intersection min1 max2 equal.\n";
      }
      else if( (min1_id == min2_id) || (max1_id == max2_id) ){
        numNstops = 3;
        std::cout << "intersection min1 min2\n";
      }
      else{
        int maxNodeIdDiff = std::abs(max1_id - max2_id);
        int minNodeDiff   = std::abs(min1_id - min2_id);
        if( (maxNodeIdDiff > 2) && (minNodeDiff > 2) ){
          numNstops = 2;
        }
        else{
          std::cout << " Intersection element is " << intersect->at(0) << " niet min of max gelijk\n";
          numNstops = 3;
          std::cout << "intersection 1_3\n";
        }
      }
    }
    else {// intersect->size() > 1
      std::cout << " Size is " << intersect->size() << '\n';
      numNstops = 2;
    }
    // Clean memory
    delete intersect;
  }// End Not Singleton
  // numNstops = 2;
  // Return
  return numNstops;
}
//_______________________ END numEndStopsIfMerged __________________________
//_______________________________ selectcandidatesMatchNodes _____________
void selectcandidatesMatchNodes( CoordGrid const &hitMap, PathCandidate const &FirstCand,
                                 PathCandidate const &SecondCand, candidateDistobject &output)
{
  // Min and max layer nodes of input tracklets
  GridNode const &firstmaxLayer  = hitMap.GetNodeByID(FirstCand.m_maxlayerNodeId);
  GridNode const &firstminLayer  = hitMap.GetNodeByID(FirstCand.m_minlayerNodeId);
  GridNode const &secondmaxLayer = hitMap.GetNodeByID(SecondCand.m_maxlayerNodeId);
  GridNode const &secondminLayer = hitMap.GetNodeByID(SecondCand.m_minlayerNodeId);
  bool notSkewedLayer = (firstmaxLayer.m_type != GridNode::STT_TYPE_SKEW) &&
    (firstminLayer.m_type  != GridNode::STT_TYPE_SKEW) &&
    (secondmaxLayer.m_type != GridNode::STT_TYPE_SKEW) &&
    (secondminLayer.m_type != GridNode::STT_TYPE_SKEW);
  /* Layer differences between max-odes and min-nodes */
  size_t layerDiffMaxMax = ((firstmaxLayer.m_Layer > secondmaxLayer.m_Layer)? (firstmaxLayer.m_Layer - secondmaxLayer.m_Layer):
                            (secondmaxLayer.m_Layer - firstmaxLayer.m_Layer) );
  size_t layerDiffMinMin = ((firstminLayer.m_Layer > secondminLayer.m_Layer)? (firstminLayer.m_Layer - secondminLayer.m_Layer):
                            (secondminLayer.m_Layer - firstminLayer.m_Layer));
  size_t layerDiffMinMax = ((firstminLayer.m_Layer > secondmaxLayer.m_Layer)? (firstminLayer.m_Layer - secondmaxLayer.m_Layer):
                            (secondmaxLayer.m_Layer - firstminLayer.m_Layer));
  size_t layerDiffMaxMin = ((firstmaxLayer.m_Layer > secondminLayer.m_Layer)? (firstmaxLayer.m_Layer - secondminLayer.m_Layer):
                            (secondminLayer.m_Layer - firstmaxLayer.m_Layer));
  size_t layerDif = 0;
  // Select path ends  
  if( (firstminLayer.m_Layer >= secondmaxLayer.m_Layer) ) {
    output.FirstCandNodeId  = FirstCand.m_minlayerNodeId;
    output.SecondCandNodeId = SecondCand.m_maxlayerNodeId;
    layerDif = ((firstminLayer.m_Layer > secondmaxLayer.m_Layer)? (firstminLayer.m_Layer - secondmaxLayer.m_Layer):
                (secondmaxLayer.m_Layer - firstminLayer.m_Layer));
    if( layerDif > MAX_ALLOWED_LAYER_DIFF_HEAD_TAIL){
      output.isValid = false;
    }
  }
  else if ( (secondminLayer.m_Layer >= firstmaxLayer.m_Layer) ) {
    output.FirstCandNodeId  = FirstCand.m_maxlayerNodeId;
    output.SecondCandNodeId = SecondCand.m_minlayerNodeId;
    layerDif = ((secondminLayer.m_Layer > firstmaxLayer.m_Layer)? (secondminLayer.m_Layer - firstmaxLayer.m_Layer):
                (firstmaxLayer.m_Layer - secondminLayer.m_Layer));
    if( layerDif > MAX_ALLOWED_LAYER_DIFF_HEAD_TAIL){
      output.isValid = false;
    }
  }
  else if( (firstminLayer.m_Layer < secondmaxLayer.m_Layer) &&
           (firstminLayer.m_Layer > secondminLayer.m_Layer) ) {
    if( (firstminLayer.GetDistance2D(secondmaxLayer)) > (firstmaxLayer.GetDistance2D(secondminLayer)) ){
      output.FirstCandNodeId  = FirstCand.m_minlayerNodeId;
      output.SecondCandNodeId = SecondCand.m_maxlayerNodeId;
    }
    else {
      output.FirstCandNodeId  = FirstCand.m_maxlayerNodeId;
      output.SecondCandNodeId = SecondCand.m_minlayerNodeId;
    }
  }
  // FIXME twijfel geval
  else if( (firstminLayer.m_Layer <= secondminLayer.m_Layer) &&
           (firstmaxLayer.m_Layer <= secondmaxLayer.m_Layer) ) {
    output.FirstCandNodeId  = FirstCand.m_minlayerNodeId;
    output.SecondCandNodeId = SecondCand.m_maxlayerNodeId;
  }
  else{
    /* Mark pair as invalid.*/
    //output.isValid = false;
  }
  if( ( (layerDiffMaxMax == 0) || (layerDiffMinMin == 0) ||
        (layerDiffMinMax == 0) || (layerDiffMaxMin == 0) ) &&
      (notSkewedLayer)){
    //output.isValid = false;
  }
}
//___________________________ END selectcandidatesMatchNodes _____________
//_________________________ candidateShortestdistPairs _____________________
void candidateShortestdistPairs( CoordGrid const &hitMap, PathCandidate const& inputCand,
                                 std::vector < PathCandidate* > const &CandContainer,
                                 std::vector<candidateDistobject> &OutCandPairs,
                                 float maxDistTh, bool useMahalanobis)
{
  // Clear output container
  OutCandPairs.clear();
  // Local variables
  size_t indexOfInput = 0;
  float minDist = std::numeric_limits<float>::max();
  float distHeadHead, distHeadTail, distTailHead, distTailTail;
  distHeadHead = distHeadTail = distTailHead = distTailTail = 0.0;
  // Find the index of input candidate in the list of all candidates
  while(inputCand.m_id != (CandContainer[indexOfInput])->m_id) {
    indexOfInput++;
  }
  // Nodes with min and max layer-number of the input
  // candidate(tracklet)
  GridNode const &inputHead = hitMap.GetNodeByID(inputCand.m_maxlayerNodeId);
  GridNode const &inputTail = hitMap.GetNodeByID(inputCand.m_minlayerNodeId);
  // Process all candidates
  for(size_t i = 0; i < CandContainer.size(); ++i) {
    PathCandidate const *secondCand = CandContainer[i];
    // They are not the same candidate!!
    if ( inputCand.m_id != secondCand->m_id){
      if( (inputCand.m_isValid) && (secondCand->m_isValid) ){
        // Reset to find minimum
        minDist = std::numeric_limits<float>::max();
        GridNode const &secondHead = hitMap.GetNodeByID(secondCand->m_maxlayerNodeId);
        GridNode const &secondTail = hitMap.GetNodeByID(secondCand->m_minlayerNodeId);
        distHeadHead = inputHead.GetDistance2D(secondHead);
        distHeadTail = inputHead.GetDistance2D(secondTail);
        distTailHead = inputTail.GetDistance2D(secondHead);
        distTailTail = inputTail.GetDistance2D(secondTail);
        // Object to hold the current pair of tracklets
        candidateDistobject dstObj;
        dstObj.indexFirstCand    = indexOfInput;
        dstObj.FirstCandidateId  = inputCand.m_id;
        dstObj.indexSecondCand   = i;
        dstObj.SecondCandidateId = secondCand->m_id;
        // Find number of end-stops if merged.
        size_t endsCount = numEndStopsIfMerged(hitMap, inputCand, (*secondCand));
        dstObj.isValid = (endsCount == 2);// ? true : false;
        // Determine which ends to use.
        selectcandidatesMatchNodes(hitMap, inputCand, (*secondCand), dstObj);
        // If an acceptable pair has been found
        if(dstObj.isValid) {
          // Find minimum distance
          if(distHeadHead < minDist){
            minDist = distHeadHead;
            dstObj.distance = minDist;
            //dstObj.FirstCandNodeId  = inputCand.m_maxlayerNodeId;
            //dstObj.SecondCandNodeId = secondCand->m_maxlayerNodeId;
          }
          if( distHeadTail < minDist){
            minDist = distHeadTail;
            dstObj.distance = minDist;
            //dstObj.FirstCandNodeId  = inputCand.m_maxlayerNodeId;
            //dstObj.SecondCandNodeId = secondCand->m_minlayerNodeId;
          }
          if( distTailHead < minDist){
            minDist = distTailHead;
            dstObj.distance = minDist;
            //dstObj.FirstCandNodeId  = inputCand.m_minlayerNodeId;
            //dstObj.SecondCandNodeId = secondCand->m_maxlayerNodeId;
          }
          if( distTailTail < minDist){
            minDist = distTailTail;
            dstObj.distance = minDist;
            //dstObj.FirstCandNodeId  = inputCand.m_minlayerNodeId;
            //dstObj.SecondCandNodeId = secondCand->m_minlayerNodeId;
          }
          // Fetch head/tail nodes with smallest Euclidean distance
          GridNode const &FirstCandHead  = hitMap.GetNodeByID(dstObj.FirstCandNodeId);
          GridNode const &SecondCandHead = hitMap.GetNodeByID(dstObj.SecondCandNodeId);
          // Determine mahalanobis distances
          dstObj.mahalanobis_FirstNode_SecondCand = ComputeMahalanobisDistNodeToCandidate((*secondCand), FirstCandHead);
          dstObj.mahalanobis_SecondNode_FirstCand = ComputeMahalanobisDistNodeToCandidate(inputCand, SecondCandHead);
          dstObj.mahalanobisDist = std::max(dstObj.mahalanobis_FirstNode_SecondCand, dstObj.mahalanobis_SecondNode_FirstCand);
          // If in the acceptable distance range, add to the output list
          if( (useMahalanobis) && (dstObj.mahalanobisDist < maxDistTh) ) {
            // Add to output list
            OutCandPairs.push_back(dstObj);
          }
          else if( !(useMahalanobis) && (dstObj.distance < maxDistTh) ) {
            // Add to output list
            OutCandPairs.push_back(dstObj);
          }
          std::cout << "min dist = " << dstObj.distance << " for " << inputCand.m_id
                    << " ";
        }// If the selected pair is valid
      }// If both valid
    }// END if not the same candidate
  }// END all candidates loop
}
//_______________________ END candidateShortestdistPairs _____________________
//_________________________ candidateShortestdistPairsGraphPathBased __________
void candidateShortestdistPairsGraphPathBased( CoordGrid const &hitMap, PathCandidate const &inputCand,
                                               std::vector < PathCandidate* > const &CandContainer,
                                               std::vector<candidateDistobject> &OutCandPairs,
                                               float maxDistThreshold, bool useMahalanobis,
                                               unsigned int maxNeighbourDepth)
{
  // Clear output container
  OutCandPairs.clear();
  size_t indexOfInput = 0;
  size_t inputCandId  = inputCand.m_id;
  // Find the index of input candidate in the list
  while(inputCand.m_id != (CandContainer[indexOfInput])->m_id) {
    indexOfInput++;
  }
  // Nodes with min and max layer-number of the input
  // candidate(tracklet)
  GridNode const &inputHead = hitMap.GetNodeByID(inputCand.m_maxlayerNodeId);
  GridNode const &inputTail = hitMap.GetNodeByID(inputCand.m_minlayerNodeId);
  // Find the shortest distance candidates to both nodes, head and tail.
  std::vector<int> neighboursHeadNode;
  std::vector<int> neighboursTailNode;
  // Populate the list of neighbours in the range of allowed depth for
  // both head and tail.
  hitMap.ActiveNeighboursList(inputHead.m_detID, neighboursHeadNode, maxNeighbourDepth);
  hitMap.ActiveNeighboursList(inputTail.m_detID, neighboursTailNode, maxNeighbourDepth);
  // DEBUG info may be removed.
  std::cout << "<DEBUG> NeighboursHeadNode " << neighboursHeadNode.size()
            << " neighboursTailNode " << neighboursTailNode.size() << '\n';
  // END DEBUG
  std::vector<int> allNeighbourNodes(neighboursHeadNode);
  allNeighbourNodes.insert(allNeighbourNodes.end(), neighboursTailNode.begin(), neighboursTailNode.end());
  // Index of valid candidates containing the nodes.
  std::set<size_t> setOfValidCandIndices;
  // Loop through all candidates
  for(size_t i = 0; i < CandContainer.size(); ++i) {
    PathCandidate const *currentCand = CandContainer[i];
    // If valid and not the same candidate
    if( (currentCand->m_isValid) && (currentCand->m_id != inputCandId) ) {
      // Selected neighbour nodes loop
      for(size_t j = 0; j < allNeighbourNodes.size(); ++j) {
        int neighbourID = allNeighbourNodes[j];
        if( currentCand->isInCandidate(neighbourID) ) {
          setOfValidCandIndices.insert(i);
        }
      }// End neighbour node loop
    }// End if valid and not the same
  }// End CandidateContainer loop
  // All valid candidates containing at least one of the neighbour ids
  // are collected (process list).

  /// DEBUG MAG WEG
  std::cout << "<DEBUG> Found " << setOfValidCandIndices.size() << " possible merge candidates.\n";
  //// DEBUG MAG WEG

  float minDist = std::numeric_limits<float>::max();
  float distHeadHead, distHeadTail, distTailHead, distTailTail;
  distHeadHead = distHeadTail = distTailHead = distTailTail = 0.0;
  for(std::set<size_t>::iterator it = setOfValidCandIndices.begin(); it != setOfValidCandIndices.end(); ++it) {
    size_t candIndex = *it;
    PathCandidate const *secondCnd = CandContainer[candIndex];
    minDist = std::numeric_limits<float>::max();
    GridNode const &secondHead = hitMap.GetNodeByID(secondCnd->m_maxlayerNodeId);
    GridNode const &secondTail = hitMap.GetNodeByID(secondCnd->m_minlayerNodeId);
    distHeadHead = inputHead.GetDistance2D(secondHead);
    distHeadTail = inputHead.GetDistance2D(secondTail);
    distTailHead = inputTail.GetDistance2D(secondHead);
    distTailTail = inputTail.GetDistance2D(secondTail);
    // Object to hold the current pair of tracklets
    candidateDistobject dstObj;
    dstObj.indexFirstCand    = indexOfInput;
    dstObj.FirstCandidateId  = inputCand.m_id;
    dstObj.indexSecondCand   = candIndex;
    dstObj.SecondCandidateId = secondCnd->m_id;
    // Find number of end-stops if merged.
    size_t endsCount = numEndStopsIfMerged(hitMap, inputCand, (*secondCnd));
    dstObj.isValid = (endsCount == 2);// ? true : false;
    // Determine which ends to use.
    selectcandidatesMatchNodes(hitMap, inputCand, (*secondCnd), dstObj);
    // If an acceptable pair has been found
    if(dstObj.isValid) {
      // Find minimum distance
      if(distHeadHead < minDist){
        minDist = distHeadHead;
        dstObj.distance = minDist;
      }
      if( distHeadTail < minDist){
        minDist = distHeadTail;
        dstObj.distance = minDist;
      }
      if( distTailHead < minDist){
        minDist = distTailHead;
        dstObj.distance = minDist;
      }
      if( distTailTail < minDist){
        minDist = distTailTail;
        dstObj.distance = minDist;
      }
      // Fetch head/tail nodes with smallest Euclidean distance
      GridNode const &FirstCandHead  = hitMap.GetNodeByID(dstObj.FirstCandNodeId);
      GridNode const &SecondCandHead = hitMap.GetNodeByID(dstObj.SecondCandNodeId);
      // Determine mahalanobis distances
      dstObj.mahalanobis_FirstNode_SecondCand = ComputeMahalanobisDistNodeToCandidate((*secondCnd), FirstCandHead);
      dstObj.mahalanobis_SecondNode_FirstCand = ComputeMahalanobisDistNodeToCandidate(inputCand, SecondCandHead);
      dstObj.mahalanobisDist = std::max(dstObj.mahalanobis_FirstNode_SecondCand, dstObj.mahalanobis_SecondNode_FirstCand);
      // If in the acceptable distance range, add to the output list
      if( (useMahalanobis) && (dstObj.mahalanobisDist < maxDistThreshold) ) {
        // Add to output list
        OutCandPairs.push_back(dstObj);
      }
      else if( !(useMahalanobis) && (dstObj.distance < maxDistThreshold) ) {
        // Add to output list
        OutCandPairs.push_back(dstObj);
      }
      std::cout << "min dist = " << dstObj.distance << " for " << inputCand.m_id
                << " ";
    }// If valid object
  }//END for it(set iterator)
}
//_________________________ END candidateShortestdistPairsGraphPathBased ________
//_________________________________ mergeCandidatesMahalanobisDist _____________________
unsigned int mergeCandidatesMahalanobisDist( CoordGrid const &hitMap,
                                             std::vector < PathCandidate* > &AllCandides,
                                             float distThreshold, unsigned int firstOutputPathId,
                                             bool graphBasedSelect, bool useMahalanobisNodeDist,
                                             size_t maxNumLayerDepth)
{

  std::cout << "<INFO> Try to merge Mahalanobis distance based. input has "
            << AllCandides.size() << " members. \n";
  unsigned int lastInputId = firstOutputPathId;
  std::vector<candidateDistobject> mergePairs;
  size_t index = 0;
  size_t numInputCand = AllCandides.size();
  // To store new merged path candidates
  std::vector<PathCandidate*> mahalanobisMerged;
  // Process input list
  while( index < numInputCand ) {
    PathCandidate* currentCandidate = AllCandides[index];
    if(currentCandidate->m_isValid) {
      // Compile the list of possible merge pairs(Euclidean pre-selection)
      if(graphBasedSelect){//O(N) Graph path based
        candidateShortestdistPairsGraphPathBased(hitMap,(*currentCandidate), AllCandides,
                                                 mergePairs, distThreshold, useMahalanobisNodeDist,
                                                 maxNumLayerDepth);
      }
      else {// Euclidean dist between all tracklets(O(N^2))
        candidateShortestdistPairs(hitMap,(*currentCandidate), AllCandides,
                                   mergePairs, distThreshold, useMahalanobisNodeDist);
      }
      // If there are any merge pairs
      if(!mergePairs.empty()) {
        // Sort to find the minimum
        if( mergePairs.size() > 1) {
          std::sort(mergePairs.begin(), mergePairs.end(), lessthanMahalanobis);
        }
        // Get candidate with the smallest distance
        candidateDistobject &bestMatchPair = mergePairs[0];
        PathCandidate *LeftCandidate   = AllCandides[bestMatchPair.indexFirstCand];
        std::set<int> *LeftMembersSet  = LeftCandidate->m_memberIdSet;
        PathCandidate *RightCandidate  = AllCandides[bestMatchPair.indexSecondCand];
        std::set<int> *RightMembersSet = RightCandidate->m_memberIdSet;
        // Create merged tracklet
        PathCandidate *mergedCand = new PathCandidate();
        mergedCand->m_id = lastInputId;
        lastInputId++;
        mergedCand->m_level = std::max(LeftCandidate->m_level, RightCandidate->m_level) + 1;
        mergedCand->m_orientation = -1;
        (mergedCand->m_memberIdSet)->insert(LeftMembersSet->begin(), LeftMembersSet->end());
        (mergedCand->m_memberIdSet)->insert(RightMembersSet->begin(), RightMembersSet->end());
        mergedCand->m_isMerged = true;
        (mergedCand->m_parents).push_back(LeftCandidate->m_id);
        (mergedCand->m_parents).push_back(RightCandidate->m_id);
        (LeftCandidate->m_childeren).push_back(mergedCand->m_id);
        (RightCandidate->m_childeren).push_back(mergedCand->m_id);
        // Update Head and tail nodes
        updatePathCandidateHeadAndTailNodes((*mergedCand));
        // Update end nodes min and may layers
        updatePathCandidateHeadAndTailLayer(hitMap, (*mergedCand));
        // Add cov, mean, invCov
        computeCandCovariance(hitMap, (*mergedCand));
        // Add to output list
        mahalanobisMerged.push_back(mergedCand);
        // Marked: hier ben je bezig. Markeer beide ouders invalid
        // Zie of dit klopt.
        //if(mergePairs.size() == 1) {
        LeftCandidate->m_isValid  = false;
        RightCandidate->m_isValid = false;
        //}
      }
      else{
        std::cout << "\t<Warning> Mahalanobis. empty list for"
                  << " Cand ID = " << currentCandidate->m_id << " dist = " << distThreshold
                  << '\n'; 
      }
    }// If (currentCandidate is valid)
    index++;
  }// END of While(index < numInputCand)
  AllCandides.insert(AllCandides.end(), mahalanobisMerged.begin(), mahalanobisMerged.end());
  std::cout << "\t<-I-> Added " << mahalanobisMerged.size()
            << " input is exteded to " << AllCandides.size() << '\n';
  // Return output par(id number for the next merged candidate)
  return lastInputId;
}
//_______________________ END mergeCandidatesMahalanobisDist ____________________
//__________________________ candidateCurvCompatPairList ________________________
void candidateCurvCompatPairList(CoordGrid const &hitMap, PathCandidate const& inputCand,
                                 std::vector < PathCandidate* > const &CandContainer,
                                 std::vector<candidateDistobject> &OutCandPairs, float threshold)
{
  // Clean the output list, to avoid mistakes.
  OutCandPairs.clear();
  // If input hasn't a valid curvature, break.
  if(!inputCand.m_hasCurvature){
    std::cout << "<DEBUG> Not processing candidate " << inputCand.m_id
              << " no curvature.\n";
    return;
  }
  // Local variables
  CurvatureParameters const &firstCurv = inputCand.m_CurV_par;
  float curvatureDiff = 0.00;
  float centreDist    = 0.00;
  float radiusdiff    = 0.00;
  size_t indexOfInput = 0;
  // Find the index of input candidate
  while(inputCand.m_id != (CandContainer[indexOfInput])->m_id) {
    indexOfInput++;
  }
  // Visit all candidates to compile the list.
  for(size_t i = 0; i < CandContainer.size(); ++i) {
    PathCandidate const *secondCand = CandContainer[i];
    if(inputCand.m_id != secondCand->m_id){
      if( (inputCand.m_isValid) && (secondCand->m_isValid) ){
        if( secondCand->m_hasCurvature ){
          CurvatureParameters const &secondCurv = secondCand->m_CurV_par;
          radiusdiff    = fabs(firstCurv.m_ra - secondCurv.m_ra);
          curvatureDiff = fabs(firstCurv.m_r - secondCurv.m_r);
          centreDist = sqrt( ((firstCurv.m_a - secondCurv.m_a) * (firstCurv.m_a - secondCurv.m_a)) +
                             ((firstCurv.m_b - secondCurv.m_b) * (firstCurv.m_b - secondCurv.m_b)) );
          if( curvatureDiff < threshold){
            candidateDistobject mergePair;
            mergePair.indexFirstCand    = indexOfInput;
            mergePair.FirstCandidateId  = inputCand.m_id;
            mergePair.indexSecondCand   = i;
            mergePair.SecondCandidateId = secondCand->m_id;
            mergePair.curvDiff   = curvatureDiff;
            mergePair.radiusdiff = radiusdiff;
            mergePair.centreDist = centreDist;
            // Add to output list
            OutCandPairs.push_back(mergePair);
            std::cout << "\t<-D-> curvatureDiff = " << curvatureDiff
                      << " centreDist = " << centreDist
                      << " radiusdiff = " << radiusdiff << '\n';
          }// END if (curvatureDiff < ...)
        }// End second has a valid curvature
      }// If both valid
    }// If not the same candidate
  }//END list of candidates loop
}
//__________________________ END candidateCurvCompatPairList ___________________
//________________________ mergeCandidatesCurvatureCompat _______________________
unsigned int mergeCandidatesCurvatureCompat( CoordGrid const &hitMap,
                                             std::vector < PathCandidate* > &CandidatesContainer,
                                             float differenceThreshold, unsigned int firstOutputCandidateId)
{
  std::cout << "<INFO> Try to merge curvature compatible paths.\n\t Input has "
            << CandidatesContainer.size()
            << " members." << " differenceThreshold = " << differenceThreshold << '\n';
  unsigned int outputId = firstOutputCandidateId;
  std::vector<candidateDistobject> mergePairs;
  /* First we need to determine curvature for all available and valid candidates.*/
  // FIXME_LATER: Dit is niet efficient, lage urgentie want het klopt
  computeCurvatureForAllPaths(hitMap, CandidatesContainer);
  // HIERBEN JE BEZIG. FIXME FIXME
  /* Compile the list of possible merge candidates, possible iff curvature is valid.*/
  int teller = 0;
  for(size_t i = 0; i < CandidatesContainer.size(); ++i) {
    PathCandidate *ca = CandidatesContainer[i];
    if( ca->m_isValid ) {
      candidateCurvCompatPairList(hitMap, (*ca), CandidatesContainer, mergePairs, differenceThreshold);
      teller++;
      //+++++++ DEBUG DELETE ME
      std::cout << "<DEBUG> Found " << mergePairs.size()
                << " possible pairs to merge for threshold = " << differenceThreshold << '\n';
      //++++ DEBUG DELETE ME
    }
  }
  std::cout << "Aantal valide paden was " << teller << '\n';
  return outputId;
}
//________________________ END mergeCandidatesCurvatureCompat ___________________
//___________________________ mergeCandidatesCurvaturePriority __________________
unsigned int mergeCandidatesCurvaturePriority( CoordGrid const &hitMap,
                                               std::vector < PathCandidate* > &CandidesList,
                                               float distThreshold, unsigned int firstOutputCandidateId,
                                               bool useMahalanobisSelection)
{
  std::cout << "<INFO> Try to merge curvature based. input has " << CandidesList.size()
            << " members. \n";
  unsigned int lastInputId = firstOutputCandidateId;
  std::vector<unsigned int> hasCurvature;
  std::vector<unsigned int>::iterator findIt;
  std::vector<candidateDistobject> mergePairs;
  size_t index = 0;
  size_t numInputCand = CandidesList.size();
  // To store new merged path candidates
  std::vector<PathCandidate*> curvaturMerged;
  // Process input list
  while(index < numInputCand){
    PathCandidate *currentCandidate      = CandidesList[index];
    std::set<int> *currentCandMembersSet = currentCandidate->m_memberIdSet;
    // Compile list of possible merge pairs for current sub-path( Euclidean dist).
    if( !useMahalanobisSelection) {
      candidateShortestdistPairs(hitMap,(*currentCandidate), CandidesList, mergePairs, distThreshold, false);
    }
    else {
      // Compile list of possible merge pairs for current sub-path( Mahalanobis dist).
      candidateShortestdistPairs(hitMap,(*currentCandidate), CandidesList, mergePairs, distThreshold, true);
    }
    // If there are no pairs to merge
    if(mergePairs.empty()){
      std::cout << "\t<Warning> Curvatur. Could not merge index " << index
                << " Cand ID = " << currentCandidate->m_id
                << " empty list for dist = " << distThreshold << '\n';
    }
    else {
      // Determine curvature for all members of selected pairs
      // (current is always the first)
      findIt = std::find(hasCurvature.begin(), hasCurvature.end(), currentCandidate->m_id);
      if(findIt == hasCurvature.end()) {
        computePathCurvature(hitMap, (*currentCandidate));
        hasCurvature.push_back(currentCandidate->m_id);
        currentCandidate->m_hasCurvature = true;
      }
      for(size_t j = 0; j < mergePairs.size(); ++j) {
        candidateDistobject &pair = mergePairs[j];
        PathCandidate *thisCand = CandidesList[pair.indexSecondCand];
        std::set<int> *thisCandMembersSet = thisCand->m_memberIdSet;
        // If curvature was not calculated for this path
        findIt = std::find(hasCurvature.begin(), hasCurvature.end(), thisCand->m_id);
        if(findIt == hasCurvature.end()) {
          computePathCurvature(hitMap, (*thisCand));
          hasCurvature.push_back(thisCand->m_id);
          thisCand->m_hasCurvature = true;
        }
        // Compute curvature for combination
        std::vector<int> allMergedNodes(currentCandMembersSet->begin(), currentCandMembersSet->end());
        allMergedNodes.insert(allMergedNodes.end(), thisCandMembersSet->begin(), thisCandMembersSet->end());
        ComputeCurvatureForListOfnodes(hitMap, allMergedNodes, pair.m_curvaturePar);
      }// END for all mergelist elements
      // Try to merge
      double meanE = 0.00;
      for(size_t k = 0; k < mergePairs.size(); ++k) {
        candidateDistobject &pairTomerge = mergePairs[k];
        PathCandidate *othCand  = CandidesList[pairTomerge.indexSecondCand];
        std::set<int> *othCandMembersSet = othCand->m_memberIdSet;
        // Compute mean E(a,b,R) for both input candidates
        meanE = ( (currentCandidate->m_CurV_par).m_E + (othCand->m_CurV_par).m_E )/ 2.00;
        // If (E_Merged < E_mean), then merge
        if( (pairTomerge.m_curvaturePar).m_E <= meanE ) {
          PathCandidate *newCandidate = new PathCandidate();
          (newCandidate->m_memberIdSet)->insert(currentCandMembersSet->begin(), currentCandMembersSet->end());
          (newCandidate->m_memberIdSet)->insert(othCandMembersSet->begin(), othCandMembersSet->end());
          newCandidate->m_CurV_par = pairTomerge.m_curvaturePar;
          newCandidate->m_hasCurvature = true;
          newCandidate->m_id = lastInputId;
          lastInputId++;
          newCandidate->m_level = std::max(currentCandidate->m_level, othCand->m_level) + 1;
          newCandidate->m_orientation = -1;
          newCandidate->m_isMerged = true;
          (newCandidate->m_parents).push_back(currentCandidate->m_id);
          (newCandidate->m_parents).push_back(othCand->m_id);
          (currentCandidate->m_childeren).push_back(newCandidate->m_id);
          (othCand->m_childeren).push_back(newCandidate->m_id);
          // Update Head and tail nodes
          updatePathCandidateHeadAndTailNodes((*newCandidate));
          // Add cov, mean, invCov
          computeCandCovariance(hitMap, (*newCandidate));
          // Add to list
          curvaturMerged.push_back(newCandidate);
        }// If smaller than mean
        else {
          std::cerr << "<Warning> Could not merge " << currentCandidate->m_id
                    << " with E = " << (currentCandidate->m_CurV_par).m_E << " and "
                    << othCand->m_id << " with E = " << (othCand->m_CurV_par).m_E
                    << " meanE = " << meanE
                    << " merged E = " << (pairTomerge.m_curvaturePar).m_E
                    << '\n';
        }
      }
    }//END else
    index++;
  }// END while
  // Add to output list
  CandidesList.insert(CandidesList.end(), curvaturMerged.begin(), curvaturMerged.end());
  return lastInputId;
}
//___________________________ END mergeCandidatesCurvaturePriority ______________
//_______________________  computePathCovariance _____________________________
void computePathCovariance(CoordGrid const &hitMap, std::vector<int> const &nodeIds,
                           std::vector< float > &mean, std::vector< std::vector<float> > &covMatrix,
                           std::vector< std::vector<float> > &inveCovMatrix)
{
  float meanX, meanY, axRotFact;
  meanX = meanY = 0.00;
  // To avoid that the matrix becomes zero, we may rotate the axis
  axRotFact = 1.00/12.00;
  // Compute mean
  for(size_t i = 0; i < nodeIds.size(); ++i) {
    int nodeID = nodeIds[i];
    GridNode const &node = hitMap.GetNodeByID(nodeID);
    meanX += node.m_x;
    meanY += node.m_y;
  }
  // Normalise means
  meanX /= static_cast<float>(nodeIds.size());
  meanY /= static_cast<float>(nodeIds.size());
  // Fill output par
  mean[0] = meanX;
  mean[1] = meanY;
  // Compute var(X), var(Y) and coVar(X,Y)
  float varX, varY, covarXY;
  varX = varY = covarXY = 0.00;
  for(size_t i = 0; i < nodeIds.size(); ++i) {
    int nodeID = nodeIds[i];
    GridNode const &node = hitMap.GetNodeByID(nodeID);
    varX    += (((node.m_x - meanX) * (node.m_x - meanX)) + axRotFact);
    varY    += (((node.m_y - meanY) * (node.m_y - meanY)) + axRotFact);
    covarXY += ((node.m_x - meanX) * (node.m_y - meanY));
  }
  // Normalise Var and coVar
  varX    /= static_cast<float>(nodeIds.size());
  varY    /= static_cast<float>(nodeIds.size());
  covarXY /= static_cast<float>(nodeIds.size());
  // Fill output par( Covar)
  covMatrix[0][0] = varX;//a
  covMatrix[0][1] = covarXY;//b
  covMatrix[1][0] = covarXY;//c
  covMatrix[1][1] = varY;//d
  // Fill invers coVar matrix
  float determinant = (varX * varY) - ( covarXY * covarXY);
  if( (determinant > 0.00) || (determinant < 0.00) ){
    determinant = 1.00 / determinant;
    inveCovMatrix[0][0] = (varY * determinant);
    inveCovMatrix[0][1] = (-1.0 * covarXY * determinant);
    inveCovMatrix[1][0] = (-1.0 * covarXY * determinant);
    inveCovMatrix[1][1] = (varX * determinant);
  }
}
//_______________________  END computePathCovariance _________________________
//________________________ computeCandCovariance ___________________________
void computeCandCovariance(CoordGrid const &hitMap, PathCandidate &cand)
{
  std::set<int> const *Members = cand.m_memberIdSet;
  std::vector<int> memVector(Members->begin(), Members->end());
  computePathCovariance(hitMap, memVector, cand.m_meanVector, cand.m_covMatrix, cand.m_inveCovMatrix);
}
//________________________ END computeCandCovariance _______________________
//___________________________ computeCovMatForAllCandidates __________________
void computeCovMatForAllCandidates(CoordGrid const &hitMap, std::vector < PathCandidate* > &CandidesContainer)
{
  for(size_t i = 0; i < CandidesContainer.size(); ++i) {
    PathCandidate *cand = CandidesContainer[i];
    computeCandCovariance( hitMap, (*cand) );
  }
}
//___________________________ END computeCovMatForAllCandidates ______________
//___________________________ computeCandidateEigenValVectors _________________
void computeCandidateEigenValVectors(PathCandidate &candidate)
{
  //| a  b |
  //| c  d |
  std::vector< std::vector<float> > &invCovarMat = candidate.m_inveCovMatrix;
  float trace = (invCovarMat[0][0] + invCovarMat[1][1]);
  float det   = (invCovarMat[0][0] * invCovarMat[1][1]) - (invCovarMat[0][1] + invCovarMat[1][0]);
  float L1 = (trace/2) + sqrt( ((trace * trace)/4.0) - det );
  float L2 = (trace/2) - sqrt( ((trace * trace)/4.0) - det );
  candidate.m_eigevValue1 = L1;
  candidate.m_eigevValue2 = L2;
  // if b not zero
  if( (invCovarMat[0][1] < 0.0) || (invCovarMat[0][1] > 0.0)) {
    candidate.m_eigenVector1[0] = invCovarMat[0][1];//b
    candidate.m_eigenVector1[1] = L1 - invCovarMat[0][0];//L1 - a
    candidate.m_eigenVector2[0] = invCovarMat[0][1];//b
    candidate.m_eigenVector1[1] = L2 - invCovarMat[0][0];// L2 -a
  }
  // if c is not zero
  else if( (invCovarMat[1][0] < 0.0) || (invCovarMat[1][0] > 0.0) ) {
    candidate.m_eigenVector1[0] = L1 - invCovarMat[1][1];// L1 - d
    candidate.m_eigenVector1[1] = invCovarMat[1][0];//c
    candidate.m_eigenVector2[0] = L2 - invCovarMat[1][1];//L2 - d
    candidate.m_eigenVector1[1] = invCovarMat[1][0];//c;
  }
  //If both b and c are zero
  else if( (std::fabs(invCovarMat[0][1]) < LOCAL_ZERO_COMPARE_EPSILON) &&
           (std::fabs(invCovarMat[1][0]) < LOCAL_ZERO_COMPARE_EPSILON) ){
    candidate.m_eigenVector1[0] = 1.0;
    candidate.m_eigenVector1[1] = 0.0;
    candidate.m_eigenVector2[0] = 0.0;
    candidate.m_eigenVector1[1] = 1.0;
  }
}
//__________________________ END computeCandidateEigenValVectors _______________
//______________________________ computeEigenValVectorsAllCandidates ____________
void computeEigenValVectorsAllCandidates(std::vector < PathCandidate* > &allCandiddates)
{
  for(size_t i = 0; i < allCandiddates.size(); ++i) {
    PathCandidate *candidate = allCandiddates[i];
    computeCandidateEigenValVectors((*candidate));
  }
}
//__________________________ END computeEigenValVectorsAllCandidates ____________
//___________________________ ComputeMahalanobisDist _________________________
void ComputeMahalanobisDist(CoordGrid &hitMap, PathCandidate const &cand)
{
  std::vector< float > const &mean = cand.m_meanVector;
  std::vector< std::vector<float> > const &inCov = cand.m_inveCovMatrix;
  std::set<int> const *Members = cand.m_memberIdSet;
  std::set<int>::iterator it;
  std::vector<float> nodeDistToMean(2);
  std::vector<float> mahaTmp(2);
  float mahalanobisDist = 0.00;
  // Visit all members of this tracklet.
  for(it = Members->begin(); it != Members->end(); ++it) {
    int nodeId = *it;
    GridNode &node = hitMap.GetNodeByID(nodeId);
    // (x-m)
    nodeDistToMean[0] = node.m_x - mean[0];
    nodeDistToMean[1] = node.m_y - mean[1];
    // (X-m)^t * inCov * (x - m)
    mahaTmp[0] = (nodeDistToMean[0] * inCov[0][0]) + (nodeDistToMean[1] * inCov[1][0]);
    mahaTmp[1] = (nodeDistToMean[0] * inCov[0][1]) + (nodeDistToMean[1] * inCov[1][1]);
    mahalanobisDist = (mahaTmp[0] * nodeDistToMean[0]) + (mahaTmp[1] * nodeDistToMean[1]);
    node.m_mahalanobisDist = mahalanobisDist;
    //std::cout << " Maho dist = " << mahalanobisDist << std::endl;
  }
}
//___________________________ END ComputeMahalanobisDist _____________________
//_______________________ ComputeMahalanobisDistForAllCandidates ______________
void ComputeMahalanobisDistForAllCandidates( CoordGrid &hitMap,
                                             std::vector < PathCandidate* > const &CandidesContainer)
{
  for(size_t i = 0; i < CandidesContainer.size(); ++i){
    PathCandidate const *cand = CandidesContainer[i];
    ComputeMahalanobisDist(hitMap, (*cand));
  }
}
//_______________________ END ComputeMahalanobisDistForAllCandidates __________
//_______________________ ComputeMahalanobisDistNodeToCandidate _______________
float ComputeMahalanobisDistNodeToCandidate(PathCandidate const &cand, GridNode const &node)
{
  std::vector< float > const &mean = cand.m_meanVector;
  std::vector< std::vector<float> > const &inCov = cand.m_inveCovMatrix;
  std::vector<float> nodeDistToMean(2);
  std::vector<float> mahaTmp(2);
  float mahalanobisDist = 0.00;
  //(x-m)
  nodeDistToMean[0] = node.m_x - mean[0];
  nodeDistToMean[1] = node.m_y - mean[1];
  //(X-m)^t * ivcov * (x - m)
  mahaTmp[0] = (nodeDistToMean[0] * inCov[0][0]) + (nodeDistToMean[1] * inCov[1][0]);
  mahaTmp[1] = (nodeDistToMean[0] * inCov[0][1]) + (nodeDistToMean[1] * inCov[1][1]);
  mahalanobisDist = (mahaTmp[0] * nodeDistToMean[0]) + (mahaTmp[1] * nodeDistToMean[1]);
  return mahalanobisDist;
}
//_______________________ END ComputeMahalanobisDistNodeToCandidate ___________
//____________________________ markAllShortPaths _____________________
void markAllShortPaths(std::vector<PathCandidate*> &Candidates, size_t minL)
{
  std::cout << "<INFO> Mark All paths shorter than " << minL << '\n';
  size_t count = 0;
  for(size_t i = 0; i < Candidates.size(); ++i) {
    PathCandidate* cand = Candidates[i];
    std::set<int> *candMemSet = cand->m_memberIdSet;
    if( candMemSet->size() < minL) {
      cand->m_isShort = true;
      count++;
    }
  }
  std::cout << "\t<-I-> Marked " << count << " out of "
            << Candidates.size() << " as short.\n";
}
//________________________ END markAllShortPaths _____________________
//___________________________ markSingletonsAsInvalid ________________
void markAsInvalidLengthBased(std::vector < PathCandidate* > &CandidesContainer, size_t minLength)
{
  size_t length = 0;
  for(size_t i = 0; i < CandidesContainer.size(); ++i) {
    PathCandidate *cand = CandidesContainer[i];
    length = (cand->m_memberIdSet)->size();
    if(length < minLength) {
      cand->m_isValid = false;
    }
  }
}
//_______________________ END markSingletonsAsInvalid ________________
void markMergedAsInvalid(std::vector < PathCandidate* > &AllCandides, unsigned int minlevel)
{
  for(size_t i = 0; i < AllCandides.size(); ++i) {
    PathCandidate *cand = AllCandides[i];
    if( cand->m_level < minlevel) {
      cand->m_isValid = false;
    }
  }
}

void markInvalidLevelBased(std::vector < PathCandidate* > &AllCandides, unsigned int level)
{
  for(size_t i = 0; i < AllCandides.size(); ++i) {
    PathCandidate *cand = AllCandides[i];
    if( cand->m_level != level) {
      cand->m_isValid = false;
    }
  }
}

void resetValidy(std::vector < PathCandidate* > &AllCandides)
{
  for(size_t i = 0; i < AllCandides.size(); ++i) {
    PathCandidate *cand = AllCandides[i];
    cand->m_isValid = true;
  }
}

void MarkInvalidWithChilderen(std::vector < PathCandidate* > &AllCandides, unsigned int numChilderen)
{
  for(size_t i = 0; i < AllCandides.size(); ++i) {
    PathCandidate *cand = AllCandides[i];
    if( (cand->m_childeren).size() >= numChilderen ){
      cand->m_isValid = false;
    }
  }
}
//________________________ mergeMarkOverlappingAsInvalid _____________
/* Dit moet nog veranders worden. In plaats van toevoegen aan een van
   de sub-paths, er moet een nieuwe aangemaakt worden en toevoegen aan
   de lijst.*/
void mergeMarkOverlappingAsInvalid(std::vector < PathCandidate* > &AllCandides, int minOverlap)
{
  // Merge and mark overlapping candidates.
  size_t counter = 0;
  for(size_t i = 0; i < AllCandides.size(); ++i) {
    std::vector<candidateDistobject> possibleMergePairs;
    PathCandidate *FirstCand = AllCandides[i];
    for(size_t j = 0; j < AllCandides.size(); ++j) {
      PathCandidate *SecondCand = AllCandides[j];
      if( (FirstCand->m_id != SecondCand->m_id ) &&
          (FirstCand->m_isValid && SecondCand->m_isValid) ) {
        if(minOverlap > 0) {
          // Determine intesection
          std::vector<int> *intersect = findIntersectionElements((*FirstCand), (*SecondCand));
          int numIntersectionElement = static_cast<int>(intersect->size());
          if( numIntersectionElement >= minOverlap ){//
            candidateDistobject dst;
            dst.indexFirstCand    = i;
            dst.FirstCandidateId  = FirstCand->m_id;
            dst.indexSecondCand   = j;
            dst.SecondCandidateId = SecondCand->m_id;
            possibleMergePairs.push_back(dst);
          }
          delete intersect;
        }
        else{//minOverlap <= 0
          // If one tracklet containes the other.
          std::set<int> const *F_set = FirstCand->m_memberIdSet;
          std::set<int> const *S_set = SecondCand->m_memberIdSet;
          bool contained = (F_set->size() > S_set->size()) ? (std::includes(F_set->begin(), F_set->end(),
                                                                            S_set->begin(),S_set->end())) :
            (std::includes( S_set->begin(), S_set->end(), F_set->begin(), F_set->end()));
          if(contained) {
            candidateDistobject mpair;
            mpair.indexFirstCand    = i;
            mpair.FirstCandidateId  = FirstCand->m_id;
            mpair.indexSecondCand   = j;
            mpair.SecondCandidateId = SecondCand->m_id;
            possibleMergePairs.push_back(mpair);
          }
        }
      }// If not the same && both valid
    }// END for(j = 0..) Second candidates loop
    // This is not correct. FIXME
    if( possibleMergePairs.size() == 1 ) {
      candidateDistobject &candPairs = possibleMergePairs[0];
      PathCandidate *candOne = AllCandides[candPairs.indexFirstCand];
      std::set<int> *candOneMembers = candOne->m_memberIdSet;
      PathCandidate *candTwo = AllCandides[candPairs.indexSecondCand];
      std::set<int> *candTwoMembers = candTwo->m_memberIdSet;
      candTwoMembers->insert(candOneMembers->begin(), candOneMembers->end());
      candOne->m_isValid = false;
      counter++;
    }
  }// End for (i = 0 ...)
  std::cout << "<INFO> Number of merged = " << counter << '\n';
}
//___________________ END mergeMarkOverlappingAsInvalid ____________________
//_______________________ findIntersectionElements _________________________
std::vector<int>* findIntersectionElements(PathCandidate const &F_cand, PathCandidate const &S_cand)
{
  std::vector<int>::iterator findIt;
  std::set<int> const *FirstMemberSet  = F_cand.m_memberIdSet;
  std::set<int> const *SecondMemberSet = S_cand.m_memberIdSet;
  std::vector<int> *intersection = new std::vector<int>(FirstMemberSet->size() + SecondMemberSet->size());
  // Determine set intersection
  findIt = std::set_intersection( FirstMemberSet->begin(), FirstMemberSet->end(),
                                  SecondMemberSet->begin(), SecondMemberSet->end(),
                                  intersection->begin());
  // Resize to the actual intersection set
  intersection->resize( (findIt - intersection->begin()) );
  return intersection;
}
//___________________ END findIntersectionElements _________________________
//_______________________ updatePathCandidateHeadAndTailLayer ______________
void updatePathCandidateHeadAndTailLayer(CoordGrid const &hitMap, PathCandidate &cand)
{
  if( (cand.m_maxlayerNodeId < 0) || (cand.m_minlayerNodeId < 0) ){
    cand.updateHeadAndTailNodes();
  }
  GridNode const &maxNode  = hitMap.GetNodeByID(cand.m_maxlayerNodeId);
  cand.m_maxLayernodeLayer = maxNode.m_Layer;
  GridNode const &minNode  = hitMap.GetNodeByID(cand.m_minlayerNodeId);
  cand.m_minlayerNodeLayer = minNode.m_Layer;
}
//____________________ END updatePathCandidateHeadAndTailLayer _____________
//__________________________ updateAllCandidateHeadAndTailLayers ___________
void updateAllCandidateHeadAndTailLayers(CoordGrid const &hitMap, std::vector<PathCandidate*> &Candidates)
{
  for(size_t i = 0; i < Candidates.size(); ++i) {
    PathCandidate *cand = Candidates[i];
    updatePathCandidateHeadAndTailLayer(hitMap, (*cand));
  }
}
//______________________ END updateAllCandidateHeadAndTailLayers ___________