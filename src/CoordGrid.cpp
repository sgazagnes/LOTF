/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
// Standard C && C++ headers.
#include <iostream>
#include <algorithm>
#include <set>

// Local includes.
#include "auxfunctions.h"
#include "hitcoordinate.h"
#include "gridNode.h"
#include "CoordGrid.h"
#include "logc.h"

#define COORD_GRID_DEBUG 0
#if(COORD_GRID_DEBUG > 0)
#include <cassert>
#endif


//_____ DEBUG
#if (WRITE_GRID_TO_ASCII > 0)
#include <fstream>
#include "TNtuple.h"
#endif
//_____ DEBUG

struct ordering {
  bool operator ()(pair<int, double> const& a, 
		   pair<int, double> const& b) {
        return a.second < b.second;
    }
};

// Constructor
CoordGrid::CoordGrid()
  : m_grid(std::vector< GridNode >()),
    m_MVD_grid(std::vector< GridNode >()),
    m_STT_idx(std::set< int >())
{}

// Destructor
CoordGrid::~CoordGrid()
{
  m_grid.clear();
  m_MVD_grid.clear();
  m_STT_idx.clear();
}

// Init Grid data structure.
void CoordGrid::Initialize(std::vector < GridNode > const &detNodes)
{
  dbggrid("Init Grid grid for %lu detectors.",detNodes.size());

  // Clean structure:
   m_grid.clear();
   m_MVD_grid.clear();
   m_STT_idx.clear();
   // Fill the Detector coordinates into the grid 
   for(size_t i = 0; i < detNodes.size(); ++i) {
     // Copy current node
     if( (detNodes[i].m_type == GridNode::STT_TYPE_PARA) ||
	 (detNodes[i].m_type == GridNode::STT_TYPE_SKEW) ) {
       m_grid.push_back(detNodes[i]);
       // error("%d, %d", i, detNodes[i].m_detID);
     }
     else if((detNodes[i].m_type == GridNode::MVD_TYPE_PIXEL) ||
	     (detNodes[i].m_type == GridNode::MVD_TYPE_STRIP) ) {
       m_MVD_grid.push_back(detNodes[i]);
     }
   }
   firstVirtIdx = (int) m_grid.size()+1;
   
   dbggrid("Total number of inserts CoordGrid: %lu, first id: %d, last id: %d", m_grid.size(), m_grid[0].m_detID, m_grid[ (m_grid.size() - 1) ].m_detID);
}

void CoordGrid::CorrectLayerLimit()
{
  dbggrid("Correcting the layer Limit  variable");
  for(size_t i = 0; i < m_grid.size(); ++i) {
    GridNode &node = m_grid[i];
    if(node.m_LayerLimit) continue;
    bool found =0;
    for(size_t j =0; j < node.m_neighbors.size(); j++){
      GridNode &neigh = m_grid[Find(node.m_neighbors[j])];
      if(neigh.m_Layer > node.m_Layer)
	found = 1;
    }
    if(!found)
      node.m_LayerLimit = 1;
  }
  
}

	 
/***/
void CoordGrid::ExtendedGrid(std::vector < GridNode > const &nodes)
{
  for(size_t i = 0; i < nodes.size(); ++i) {
    m_grid.push_back(nodes[i]);
  }
}

void CoordGrid::AddNodeToGrid(GridNode const &node)
{
  m_grid.push_back(GridNode(node));
}


void  CoordGrid::FindNodeBetweenLayerNodePairs(std::vector< TubeLayerPairProperty > &NodePairSet)
{
  dbggrid("Finding tube pairs to compute virtual nodes.");

  // List of all nodes available in the image (detector map)
  std::vector< GridNode > const &Ingrid = m_grid;
  size_t layerDiff = 0;
  
  // Grid node loop to find the possible pairs.
  for(size_t i = 0; i < Ingrid.size(); ++i) {
    GridNode const &current_Node = Ingrid[i];
    // List of neighbors
    std::vector<int> const &neighList = current_Node.GetNeighbors();
    // Neighbors loop
    for( size_t j = 0; j < neighList.size(); ++j) {
      int neighbourID = neighList[j];// Neighbour ID
      int neigh_Index = Find(neighbourID);// Index in the map
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
  dbggrid("Found %lu pairs and %lu valid ones", NodePairSet.size(), cnt);
}
//____________________ END FindNodeBetweenLayerNodePairs _____________
bool compareByID(const GridNode &a, const GridNode &b)
{
    return a.m_detID < b.m_detID;
}

void CoordGrid::AddVirtualNodes(std::vector < GridNode > &VNodesLayer,  std::vector < GridNode > &VNodesSector)
{
  dbggrid("Computing grid virtual tubes (Neigbor based).");
  // Dummy local counters.
  int NumTubesAdded  = 0;
  // Start node ID for virtual nodes.
  int StartVirtualID = (int) firstVirtIdx;//START_VIRTUAL_ID;//6000;

  // List of all nodes available in the detector map.
  std::vector< GridNode > &Ingrid = m_grid;
  // Sort all nodes by their layer number
  dbggrid("Sorting nodes, increasing layer number");
  std::sort(Ingrid.begin(), Ingrid.end(), LessThanLayer);
  // Determine the pairs of nodes between the layers for which we want
  // to compute virtual nodes and discard duplicates.
  std::vector< TubeLayerPairProperty > NodePairSet;
  FindNodeBetweenLayerNodePairs(NodePairSet);

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
      //     if(firstNode.m_type == GridNode::STT_TYPE_SKEW && secondNode.m_type == GridNode::STT_TYPE_SKEW)
      //     //	IntersectionPoint(hitMap, firstNode, secondNode, Dummy_coord) ;
      //     else
      // if(firstNode.m_type == GridNode::STT_TYPE_SKEW && secondNode.m_type == GridNode::STT_TYPE_SKEW)
	IntersectionPointTubes(firstNode, secondNode, Dummy_coord);
//else
//	IntersectionPointSkeSke(firstNode, secondNode, Dummy_coord);

	// Modify node ID
	Dummy_coord.m_detID      = StartVirtualID + NumTubesAdded;
	Dummy_coord.m_Orig_detID = Dummy_coord.m_detID;
	std::pair<float, float> r_Theta;
	float theta_deg = Cartesian_To_Polar(Dummy_coord.m_xDet, Dummy_coord.m_yDet, r_Theta);
	Dummy_coord.m_r = r_Theta.first;
	Dummy_coord.m_thetaDeg = theta_deg;

	VNodesLayer.push_back(Dummy_coord);
	NumTubesAdded++;
      // If intersect
    }// END if the pair is valid
  }// End for NodePairSet
  // Reset visiting variable for all nodes in the input graph

  std::vector< std::vector<size_t> > SectorLeft;
  std::vector< std::vector<size_t> > SectorRight;

  for(size_t i = 0; i < 6; ++i) {
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
  dbggrid(" Total number of sectors = %d", SectorLeft.size());
  dbggrid(" Number of left boundary limits = %d" , SectorLeft.size());
  dbggrid(" Number of right boundary limits = %d", SectorRight.size());
  
  double currDist = 0;
  size_t R_Index = 0;
  // Loop all the sectors.
  for(size_t s = 0; s < 6; ++s) { // Numsectors = 6
    std::vector<size_t> const &currentLeftSector = SectorLeft[s];
    // The counter sector is SectorRight[s - 1]. The right limts of
    // the previous sector.
    R_Index = ( s == 0 ) ? (SectorRight.size() - 1) : (s - 1);
    std::vector<size_t> const &currentRightSector = SectorRight[R_Index];
    // For each tube find its nearest neighbor in the other sector.
    for(size_t l = 0; l < currentLeftSector.size(); ++l) {
      GridNode  &CurLftNode = Ingrid[currentLeftSector[l]];
      for(size_t r = 0; r < currentRightSector.size(); ++r) {
	GridNode  &CurRgtNode = Ingrid[currentRightSector[r]];
	// If in the same layer
	if( (labs(CurLftNode.m_Layer - CurRgtNode.m_Layer) <=1 ) &&
	    (CurLftNode.m_Sector != CurRgtNode.m_Sector) &&
	    (CurLftNode.m_type == GridNode::STT_TYPE_SKEW && CurRgtNode.m_type == GridNode::STT_TYPE_SKEW) &&
	    (CurLftNode.m_detID != CurRgtNode.m_detID)
	    ){
	  currDist = distanceBetweenTube(CurLftNode, CurRgtNode);
	  //sqrt((CurLftNode.m_x - CurRgtNode.m_x) * (CurLftNode.m_x - CurRgtNode.m_x) + (CurLftNode.m_y - CurRgtNode.m_y) * (CurLftNode.m_y - CurRgtNode.m_y));
	  
	  //if( CurLftNode.m_detID == 1318 || CurRgtNode.m_detID == 1318)
	  //    /info("%d, %d, %lf", firstNode.m_detID, secondNode.m_detID, distanceBetweenTube(firstNode, secondNode));
	  //info("Node %d and %d, dist %f", CurLftNode.m_detID, CurRgtNode.m_detID, currDist); 
	  if(currDist <3.5){
	    (CurLftNode.m_neighbors).push_back(CurRgtNode.m_detID);
	    /*	    GridNode Dummy_coord;
	    // Find intersection point.
	    IntersectionPointSkeSke(hitMap, CurLftNode, CurRgtNode, Dummy_coord);

	    Dummy_coord.m_detID      = StartVirtualID + NumTubesAdded;
	     Dummy_coord.m_Orig_detID = Dummy_coord.m_detID;
	      std::pair<float, float> r_Theta;
	      float theta_deg = Cartesian_To_Polar(Dummy_coord.m_xDet, Dummy_coord.m_yDet, r_Theta);
	      Dummy_coord.m_r = r_Theta.first;
	      Dummy_coord.m_thetaDeg = theta_deg;

	     VNodesSector.push_back(Dummy_coord);
	     NumTubesAdded++;*/
	  }
	}
      }
    }
  }
  std::sort(Ingrid.begin(), Ingrid.end(), compareByID);
  dbggrid("Determined %d virtual tubes (neighborList)", NumTubesAdded);
  // error("%lu", VNodesLayer[10].m_detID);
}

void CoordGrid::fixNeighboring()
{
  dbggrid("Correcting the missing neighbor relations.");
  size_t NumFixed = 0;
  // List of all nodes available in the image (detector map)
  //std::vector< GridNode > &Ingrid = hitMap.GetGrid();
  std::vector< GridNode > &Ingrid = m_grid;

  // Tubes loop
  for(size_t i = 0; i < Ingrid.size()-1; ++i) {

    GridNode &first = Ingrid[i];
    // neighbor tubes loop
    for(size_t j = i+1; j < Ingrid.size(); ++j) {
      GridNode &second = Ingrid[j];
      // A neig B and
      // if(first.m_detID == 1768 || second.m_detID == 1768)
      //	info("%d, %d", first.m_detID,second.m_detID);

	
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
  dbggrid("Total number of corrected neighboring relations = %lu", NumFixed);
}



void CoordGrid::isolateSectorAndLayerLimits(TNtuple &Sections, TNtuple &Layers)
{
  size_t countL, countS;
  countL = countS = 0;

  // Fetch the list of all detectors.
  //std::vector< GridNode > const &Ingrid = hitMap.GetGrid();
  std::vector< GridNode > const &Ingrid = m_grid;

  for(size_t i = 0; i < Ingrid.size(); ++i) {
    GridNode const &tube = Ingrid[i];

    if( tube.m_LayerLimit &&
	(tube.m_type != GridNode::VIRTUAL_NODE)
      ) {
      Layers.Fill(tube.m_x, tube.m_y, tube.m_zDet, tube.m_z);
      countL++;
    }
    // it returns -1 if the tube is at the front, 1 if is at the back
    // of the sector border and 0 elsewhere.
    //( (tube.m_SectorLimit == 1) || (tube.m_SectorLimit == -1) || (tube.m_SectorLimit == 0) ) &&
    if( (tube.m_SectorLimit != 0)
	&& (tube.m_type != GridNode::VIRTUAL_NODE)
      ) {
      Sections.Fill(tube.m_x, tube.m_y, tube.m_zDet, tube.m_z);
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
  dbggrid("Count layer limit = %lu, sector limit: %lu", countL, countS);
}


/*
 * Note: Current implementation is not optimal and might be very
 * slow. As it is just linear searching in a matrix.
 */
// FIXME (NOTE): We use Det_ID. This might not work for MVD and the
// other ones.
void CoordGrid::FillGrid(std::vector < HitCoordinate* > const& hitcoords)
{
  dbggrid("Fill grid with current tracks. Num input points %lu" , hitcoords.size());
            
  unsigned int  sttCnt, mvdCnt;
  sttCnt = 0;
  mvdCnt = 0;
  // FIXME for now we are adding MVD points dynamically to the
  // geometry. This needs to be fixed. WE NEED TO REMOVE THIS
  m_MVD_grid.clear();

  for(size_t i = 0; i < hitcoords.size(); ++i) {

    HitCoordinate const *hc = hitcoords[i];
    if( hc->type == HitCoordinate::STT_TYPE ) {
      int detID = hc->m_detID;
      //int idx = Find(detID);
      m_STT_idx.insert(detID);
      (m_grid[detID-1]).m_active = true;
      sttCnt++;
    }// IF STT type
    
    // For now it can only be MVD type
    else if( hc->type == HitCoordinate::MVD_TYPE ) {
      GridNode node;
      node.m_detID  = hc->m_detID;
      node.m_active = true;
      node.m_x      = hc->x;
      node.m_y      = hc->y;
      node.m_z      = hc->z;
      node.m_zDet   = hc->z;
      node.m_weight = 1;
      m_MVD_grid.push_back(node);
      mvdCnt++;
    }// if MVD type
  }// FOR hitcoordinates

  int vcnt = 0;
  /* for(size_t i = (size_t) firstVirtIdx; i < m_grid.size(); i++){
    GridNode &virt = m_grid[i];
    int id1 =  virt.m_neighbors[0];
    int id2 =  virt.m_neighbors[1];
   
    GridNode &virtN1 = m_grid[Find(id1)];
    GridNode &virtN2 = m_grid[Find(id2)];
    
    virt.m_active = (virtN1.m_active + virtN2.m_active)/2;
    vcnt += virt.m_active;
  }*/
  // All real tubes are activated.  Now we need to activate virtual
  // tubes as well. Note that every virtual tube has ONLY two
  // neighbors. These are the parent neigbours.
  dbggrid("Stored STT points = %u (indiv %u) and %u MVD points", sttCnt, m_STT_idx.size(), mvdCnt);
  dbggrid("Virt poiints %d", vcnt);

}

bool Is_STT_SplitSkewedNode( GridNode const &node)
{
  return node.IsSTTSplitSkewedNode();
}
bool IsVirtualSplitNode( GridNode const &node)
{
  return node.IsVirtualSplitNode();
}
//___________________________ ResetGrid __________________________
void CoordGrid::ResetGrid()
{
  for(size_t i = 0; i < m_grid.size(); ++i) {
    //if (i < 15) info("%d, %d", i, m_grid[i].m_detID);
    (m_grid[i]).m_active       = false;
    // The z values for virtual tubes are computed and stored in the grid,
    // so we do not need to reset them
    if( (m_grid[i]).m_type == GridNode::STT_TYPE_PARA ) {
      (m_grid[i]).m_xDet  = (m_grid[i]).m_x;
      (m_grid[i]).m_yDet  = (m_grid[i]).m_y;
      (m_grid[i]).m_zDet  = 0.00;
    }
    else if( (m_grid[i]).m_type == GridNode::STT_TYPE_SKEW ) {
      (m_grid[i]).m_xDet  = (m_grid[i]).m_x;
      (m_grid[i]).m_yDet  = (m_grid[i]).m_y;
      (m_grid[i]).m_zDet  = (m_grid[i]).m_z;
    }
    else if( (m_grid[i]).m_type == GridNode::VIRTUAL_NODE ) {
      (m_grid[i]).m_zDet  = (m_grid[i]).m_z;
    }
  }
  m_MVD_grid.clear();
  m_STT_idx.clear();
  // Set tollerance to zero
  
  /* The following code segment needs to be replaced by a more
     efficient implementation of the grid. This is very slow as the
     nodes are created and removed for each event.*/

  // Remove Skewed-Split nodes and resize container.
  std::vector< GridNode >::iterator it = std::remove_if(m_grid.begin(), m_grid.end(), Is_STT_SplitSkewedNode);
  m_grid.erase(it, m_grid.end());
  
  // Remove split virtuals and resize.
  it = std::remove_if(m_grid.begin(), m_grid.end(), IsVirtualSplitNode);
  m_grid.erase(it, m_grid.end());
}
//___________________________ END ResetGrid ______________________


int CoordGrid::Find(int id) const
{
  size_t index = 0;
  while( ( index < m_grid.size() ) &&
         ( m_grid[index]).m_detID != id ) {
    index++;
  }
  if( index >= m_grid.size() ) {
    return -1;
  }
  return index;
}

GridNode& CoordGrid::GetNodeByID(int DetId)
{
  int nodeIndex = CoordGrid::Find(DetId);
  return (m_grid[nodeIndex]);
}

GridNode const& CoordGrid::GetNodeByID(int DetId) const
{
  int nodeIndex = CoordGrid::Find(DetId);
  return (m_grid[nodeIndex]);
}

int CoordGrid::FindMVD(int mvdID) const
{
  size_t index = 0;
  while( ( index < m_MVD_grid.size() ) &&
         ( m_MVD_grid[index]).m_detID != mvdID ) {
    index++;
  }
  if( index >= m_MVD_grid.size() ) {
    return -1;
  }
  return index;
}

void CoordGrid::ActiveNeighboursAndSecondOrderNeighbours(int NodeID, std::vector<int> &ListOfNeighbours) const
{
  // Clear output parameter.
  ListOfNeighbours.clear();
  std::vector<int>::iterator findIt;
  // Fetch node and its neighbours.
  GridNode const &inputNode = this->GetNodeByID(NodeID);
  std::vector<int> const &neighbourList = inputNode.m_neighbors;
  // Collect all active direct and second order neighbours
  for(size_t i = 0; i < neighbourList.size(); ++i) {
    int n_id = neighbourList[i];
    GridNode const &direct_neighbour = this->GetNodeByID(n_id);
    std::vector<int> const &secondNeighbourList = direct_neighbour.m_neighbors;
    if(direct_neighbour.m_active){
      findIt = std::find(ListOfNeighbours.begin(), ListOfNeighbours.end(), n_id);
      // Was not added before
      if( findIt == ListOfNeighbours.end() ){
        ListOfNeighbours.push_back(n_id);
      }
    }
    for(size_t j = 0; j < secondNeighbourList.size(); ++j){
      int sec_id = secondNeighbourList[j];
      GridNode const &second_neighbour = this->GetNodeByID(sec_id);
      if( (sec_id != NodeID) && (second_neighbour.m_active) ){
        findIt = std::find(ListOfNeighbours.begin(), ListOfNeighbours.end(), sec_id);
        // Was not added before
        if( findIt == ListOfNeighbours.end() ){
          ListOfNeighbours.push_back(sec_id);
        }
      }
    }
  }
}
//______________________________ ActiveNeighboursList ___________________
void CoordGrid::ActiveNeighboursList(int NodeID, std::vector<int> &ListOfNeighbours, unsigned int depth) const
{
  // Temporary structure to store ids.
  std::set<int> activeNeighbourSet;
  size_t layerDiff = 0;
  // Fetch the node and its neigbours
  GridNode const &inputNode = this->GetNodeByID(NodeID);
  // List of direct neighbours
  std::vector<int> const &neighbourList = inputNode.m_neighbors;
  // Select active direct neighbours.
  for(size_t i = 0; i < neighbourList.size(); ++i) {
    int n_id = neighbourList[i];// node id
    GridNode const &direct_neighbour = this->GetNodeByID(n_id);//fetch node
    // Active?
    if(direct_neighbour.m_active) {
      activeNeighbourSet.insert(n_id);
    }
  }//END for(neigbourlist)
  // If go deeper!?
  unsigned int depthCounter = 0;
  if(depthCounter < depth) {
    std::vector<int> nodeIdCollection(activeNeighbourSet.begin(), activeNeighbourSet.end());
    // Store only active nodes
    std::vector<int> tmpNdIdSet(activeNeighbourSet.begin(), activeNeighbourSet.end());
    // Store all nodes on the search path
    std::vector<int> secondOrderIdList;
    // No active neighbours? use all neighbours
    if( nodeIdCollection.empty() ) {
      nodeIdCollection = inputNode.m_neighbors;
    }
    // Visit deeper layers and collect active nodes.
    while( depthCounter < depth) {
      depthCounter++;
      for(size_t j = 0; j < nodeIdCollection.size(); ++j) {
        int d_Id = nodeIdCollection[j];
        // Fetch node
        GridNode const &d_ND = this->GetNodeByID(d_Id);
        // Fetch neigbours of neighbour node, add the active ones
        std::vector<int> const &d_neighbours = d_ND.m_neighbors;
        for(size_t k = 0; k < d_neighbours.size(); ++k) {
          int nb_id = d_neighbours[k];
          // Fetch node
          GridNode const &a_node = this->GetNodeByID(nb_id);
          // Add all to the follower list
          secondOrderIdList.push_back(nb_id);
          // Active?add:nothing
          if( (a_node.m_active) ) {
            tmpNdIdSet.push_back(nb_id);
          }
        }//End neighbours loop(second order neighbours)
      }//END node collection loop
      // Found some active neighbours
      if(!tmpNdIdSet.empty()) {
        // insert in output
        activeNeighbourSet.insert(tmpNdIdSet.begin(), tmpNdIdSet.end());
        // Swap lists
        nodeIdCollection = tmpNdIdSet;
        tmpNdIdSet.clear();// Clear for next round
      }
      else{//empty active neighbours list(large gap)
        nodeIdCollection = secondOrderIdList;
      }
      // Clear the list(it might contain inactive nodes)
      secondOrderIdList.clear();
    }//END while()
  }// END if (depth)
  //________________________________________
  // Clear output list and insert collected members
  /* Not the most efficient way, but for now clear and correct.*/
  ListOfNeighbours = std::vector<int>(activeNeighbourSet.begin(), activeNeighbourSet.end());
  //DEBUG
  std::cout << "\t<DEBUG NUM NODES> Total number of depth: " << depth << " is " << ListOfNeighbours.size()
            << std::endl;
}
//__________________________ END ActiveNeighboursList ___________________


//_________________________ END isolateSectorAndLayerLimits __________


//_____________ Debug function to dump in a file
#if (WRITE_GRID_TO_ASCII > 0)
void CoordGrid::WriteGrid(std::string const& fileName)const
{
  size_t tubeCnt = 0;
  int actv = 0;

  std::string format = "########\n";
  format += "#\n# ID, sector, Layer, active, x, y, neighbour list (IDs of the neighbors)\n";
  format = format + "########\n";

  //// DEBUG DEBUG  Open ROOT File to write coordinates
  TNtuple coords ("Coords" , "Coordinates in x y plane", "x:y");
  //////// DEBUG DEBUG

  // Open File to write the coordinates ( text file)
  std::ofstream OutTxtFile;
  OutTxtFile.open (fileName.c_str());

  // Perform actual writing
  if (OutTxtFile.is_open()){
    OutTxtFile << format;
    // Write
    for(size_t j = 0; j < m_grid.size(); ++j) {
      // Current node
      GridNode const& curPoint = m_grid[j];
      // If active
      actv = curPoint.m_active ? 1 : 0;
      tubeCnt++;
      
      OutTxtFile << curPoint.m_detID  << ", "
		 << curPoint.m_Sector << ", "
		 << curPoint.m_Layer  << ", "
		 << actv << ", "
		 << curPoint.m_x  << ", "
		 << curPoint.m_y  << ", ";
      for(size_t n = 0; n < (curPoint.m_neighbors).size(); ++n) {
	OutTxtFile << (curPoint.m_neighbors)[n] << ' ';
      }
      OutTxtFile << '\n';
      coords.Fill(curPoint.m_x, curPoint.m_y);
    }// Grid loop
  }// IF(OutTxtFile....)
  OutTxtFile.close();

  ///////////
  coords.SaveAs("DUMP_GridCoordinaten.root");
  ///////////
  std::cout << "<DEBUG_INFO> The total number of written grid nodes = "
            << tubeCnt << '\n';
}
#endif
