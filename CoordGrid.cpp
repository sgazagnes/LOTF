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
#include "hitcoordinate.h"
#include "gridNode.h"
#include "CoordGrid.h"
#include "logc.h"

#define COORD_GRID_DEBUG 1
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
  : m_AttSpaceTolerance(0.0),
    m_grid(std::vector< GridNode >()),
    m_MVD_grid(std::vector< GridNode >())
{}

// Destructor
CoordGrid::~CoordGrid()
{
  m_grid.clear();
  m_MVD_grid.clear();
}

// Init Grid data structure.
void CoordGrid::Initialize(std::vector < GridNode > const &detNodes)
{
  dbggrid("Init Grid grid for %lu detectors.",detNodes.size());

  // Clean structure:
   m_grid.clear();
   m_MVD_grid.clear();

   // Fill the Detector coordinates into the grid 
   for(size_t i = 0; i < detNodes.size(); ++i) {
     // Copy current node
     if( (detNodes[i].m_type == GridNode::STT_TYPE_PARA) ||
	 (detNodes[i].m_type == GridNode::STT_TYPE_SKEW) ) {
       m_grid.push_back(detNodes[i]);

     }
     else if((detNodes[i].m_type == GridNode::MVD_TYPE_PIXEL) ||
	     (detNodes[i].m_type == GridNode::MVD_TYPE_STRIP) ) {
       m_MVD_grid.push_back(detNodes[i]);
     }
   }
   firstVirtIdx = m_grid.size();

   dbggrid("Total number of inserts CoordGrid: %lu, first id: %d, last id: %d", m_grid.size(), m_grid[0].m_detID, m_grid[ (m_grid.size() - 1) ].m_detID);
}

void CoordGrid::CorrectLayerLimit()
{
  dbggrid("Correcting the layer Limit  variable");
  int nLayerAdd = 0;
  bool found =0;
  for(size_t i = 0; i < m_grid.size(); ++i) {
    GridNode &node = m_grid[i];
    if(node.m_LayerLimit) continue;
    int numNeigh =  node.m_neighbors.size();
    for(size_t j =0; j < numNeigh; j++){
      GridNode &neigh = m_grid[Find(node.m_neighbors[j])];
      if(neigh.m_Layer > node.m_Layer)
	found = 1;
    }
    if(!found)
      node.m_LayerLimit = 1;
  
  }
}


	 
/*
 * Modify the node id numbering in order to include the virtual tubes
 * in the correct positions between the layes based on their
 * id's. Geometrically their coordinates are correct.
 *
 * "m_Orig_detID": Corresponds with the original ID's as provided by
 * the detector maps.
 *
 * "m_detID": Is the id after modifications.
 */
void CoordGrid::Rearange_nodeIds()
{
  std::cout << "<Warning> Rearange_nodeIds is called but not implemented yet.\n";
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

/*
 * Note: Current implementation is not optimal and might be very
 * slow. As it is just linear searching in a matrix.
 */
// FIXME (NOTE): We use Det_ID. This might not work for MVD and the
// other ones.
void CoordGrid::FillGrid(std::vector < HitCoordinate* > const& hitcoords)
{
  dbggrid("Fill grid with current tracks. Num input points %lu" , hitcoords.size());
            
  unsigned int  sttCnt, vCnt, mvdCnt;
  sttCnt = 0;
  mvdCnt = 0;
  // FIXME for now we are adding MVD points dynamically to the
  // geometry. This needs to be fixed. WE NEED TO REMOVE THIS
  m_MVD_grid.clear();

  for(size_t i = 0; i < hitcoords.size(); ++i) {
    HitCoordinate const *hc = hitcoords[i];
    if( hc->type == HitCoordinate::STT_TYPE ) {
      int detID = hc->m_detID;
      int idx = Find(detID);
      //m_grid[index]
	(m_grid[idx]).m_active = true;
      /* for(size_t index = 0; index < m_grid.size(); ++index) {
       	if( ( (m_grid[index]).m_Orig_detID == detID ) &&
            ( (m_grid[index]).m_type != GridNode::VIRTUAL_NODE) ) {
	  (m_grid[index]).m_active = true;
      	  sttCnt++;
       	}
	}*/
    }// IF STT type
    // For now it can only be MVD type
    else if( hc->type == HitCoordinate::MVD_TYPE ) {
      GridNode node;
      node.m_detID  = hc->m_detID;
      node.m_active = true;
      node.m_x      = hc->x;
      node.m_y      = hc->y;
      node.m_z      = hc->z;
      node.m_z_Det  = hc->z;
      node.m_weight = 1;
      m_MVD_grid.push_back(node);
      mvdCnt++;
    }// if MVD type
  }// FOR hitcoordinates

  vCnt = 0;
  // All real tubes are activated.  Now we need to activate virtual
  // tubes as well. Note that every virtual tube has ONLY two
  // neighbors. These are the parent neigbours.
   std::vector<int> tocheck;
  char *visited 	= (char *) calloc(firstVirtIdx, sizeof(char));
  /*
  for(size_t j = firstVirtIdx; j < m_grid.size(); j++) {
    GridNode &Virtual_tube = m_grid[j];
    if( Virtual_tube.m_type == GridNode::VIRTUAL_NODE) {
      std::vector<int> const &Neig_List = Virtual_tube.m_neighbors;
      //assert (Neig_List.size() == 2);
      // Note that each virtual tube can only have 2 neigbors.
      int nID = Neig_List[0];
      int nIndex = Find(nID);
      GridNode const &First  = m_grid[nIndex];
      nID = Neig_List[1];
      nIndex = Find(nID);
      GridNode const &Second = m_grid[nIndex];

      if( (First.m_active) && (Second.m_active) ) {
	if(!visited[First.m_detID]){
	  tocheck.push_back(First.m_detID);
	  visited[First.m_detID] = 1;
	}
	if(!visited[Second.m_detID]){
	  tocheck.push_back(Second.m_detID);
	  visited[Second.m_detID] = 1;
	}

	Virtual_tube.m_active = true;
	vCnt++;
      }
    }
  }


  //  printf("need to check %d \n", tocheck.size());

  for(size_t j = 0; j < tocheck.size(); j++) {
    int nID = tocheck[j];

    int nIndex = Find(nID);
    GridNode  &First  = m_grid[nIndex];
    if(First.m_neighbors.size() < 4)
      continue;

    char *neighVis 	= (char *) calloc(First.m_neighbors.size(), sizeof(char));

    for(int k = 0; k < First.m_neighbors.size(); k++){
      int neigh = First.m_neighbors[k];
      int neighIdx = Find(neigh);
      GridNode const &neighNode  = m_grid[neighIdx];
      
      if( neighNode.m_type != GridNode::VIRTUAL_NODE || neighVis[k] || !neighNode.m_active){
	neighVis[k] = 1;
	continue;
      }

      std::vector<int> inDir;
      std::vector<double> distDir;

      int cnt = 0;
      int xdir = First.m_x - neighNode.m_x > 0? 1: 0;
      int ydir = First.m_y - neighNode.m_y > 0? 1: 0;
      double dist = sqrt((First.m_x - neighNode.m_x)*(First.m_x - neighNode.m_x) +
			 (First.m_y - neighNode.m_y)*(First.m_y - neighNode.m_y));
      inDir.push_back(neigh);
      distDir.push_back(dist);
	
      for(int l = 0; l < First.m_neighbors.size(); l++){
	int neigh2 = First.m_neighbors[l];
	if(neigh == neigh2 || neighVis[l])
	  continue;
	
	int neighIdx2 = Find(neigh2);
	GridNode const &neighNode2  = m_grid[neighIdx2];
	if( neighNode2.m_type != GridNode::VIRTUAL_NODE || !neighNode2.m_active) {
	  neighVis[l] = 1;
	  continue;
	}
 
	int neighxdir = First.m_x - neighNode2.m_x > 0? 1: 0;
	int neighydir = First.m_y - neighNode2.m_y > 0? 1: 0;
	if(neighxdir != xdir || neighydir != ydir)
	  continue;
	
	double neighdist =  sqrt((First.m_x - neighNode2.m_x)*(First.m_x - neighNode2.m_x) +
				 (First.m_y - neighNode2.m_y)*(First.m_y - neighNode2.m_y));
	inDir.push_back(neigh2);
	distDir.push_back(neighdist);
	neighVis[l] = 1;	
      }
      if(inDir.size() > 2){
	dbggrid("ID %d has several virtual nodes in the same direction, we can clean ", nID);
	vector<pair<int, double>> order(inDir.size());
	for (int i=0; i<inDir.size(); ++i)
	  order[i] = make_pair(inDir[i], distDir[i]);
	sort(order.begin(), order.end(), ordering());
	
	for (int i=2; i<inDir.size(); ++i){
	  int id = order[i].first;
	  int idx = Find(id);
	  GridNode  &remove  = m_grid[idx];
	  // remove.m_active = false;
	  (First.m_neighbors).erase(std::remove((First.m_neighbors).begin(), (First.m_neighbors).end(),id), (First.m_neighbors).end());
	}
      }
      neighVis[k] = 1;
    }
  }

	//	int skeID = neighNode.m_neighbors[0] == nID? neighNode.m_neighbors[1]: neighNode.m_neighbors[0];
	//	int skeIdx = Find(skeID);
	//	GridNode const &skeNode  = m_grid[neighIdx];
	
	*/
  /*int cnt = 0;
    int xdir = First.m_x - Virtual_tube.m_x > 0? 1: 0;
    int ydir = First.m_y - Virtual_tube.m_y > 0? 1: 0;
    double dist = sqrt((First.m_x - Virtual_tube.m_x)*(First.m_x - Virtual_tube.m_x) +
    (First.m_y - Virtual_tube.m_y)*(First.m_y - Virtual_tube.m_y));

    for(int k = 0; k < First.m_neighbors.size(); k++){
    int neigh = First.m_neighbors[k];
    if(neigh == Virtual_tube.m_detID)
    continue;
    int neighIdx = Find(nID);
    GridNode const &neighNode  = m_grid[nIndex];
    if( neighNode.m_type == GridNode::VIRTUAL_NODE) {
    int neighxdir = First.m_x - neighNode.m_x > 0? 1: 0;
    int neighydir = First.m_y - neighNode.m_y > 0? 1: 0;
    if(neighxdir != xdir || neighydir != ydir)
    continue;
    double neighdist =  sqrt((First.m_x - neighNode.m_x)*(First.m_x - neighNode.m_x) +
    (First.m_y - neighNode.m_y)*(First.m_y - neighNode.m_y));
    if(neighdist < dist)
    cnt++;
    }
    }
    if(cnt > 1)
    continue;

    cnt = 0;

    for(int k = 0; k < Second.m_neighbors.size(); k++){
    int neigh = Second.m_neighbors[k];
    if(neigh == Virtual_tube.m_detID)
    continue;
    int neighIdx = Find(nID);
    GridNode const &neighNode  = m_grid[nIndex];
    if( neighNode.m_type == GridNode::VIRTUAL_NODE) {
    int neighxdir = Second.m_x - neighNode.m_x > 0? 1: 0;
    int neighydir = Second.m_y - neighNode.m_y > 0? 1: 0;
    if(neighxdir != xdir || neighydir != ydir)
    continue;
    double neighdist =  sqrt((Second.m_x - neighNode.m_x)*(First.m_x - neighNode.m_x) +
    (Second.m_y - neighNode.m_y)*(First.m_y - neighNode.m_y));
    if(neighdist < dist)
    cnt++;
    }
    }

    if(cnt > 1)
    continue;*/
  free(visited);
  dbggrid("Stored real STT points = %u,  %u virtual nodes, and %u MVD points", sttCnt,vCnt, mvdCnt);
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
    (m_grid[i]).m_active          = false;
    (m_grid[i]).m_length          = 0;
    (m_grid[i]).m_lengthFW        = 0;
    (m_grid[i]).m_lengthBW        = 0;
    (m_grid[i]).m_forwardVisited  = false;
    (m_grid[i]).m_backwardVisited = false;
    (m_grid[i]).m_area            = 0;
    (m_grid[i]).m_visited         = false;
    (m_grid[i]).m_maxPathVisited  = false;
    ((m_grid[i]).m_orientations).clear();
    (m_grid[i]).m_MVDAssigned  = false;
    (m_grid[i]).m_maxOrientVal = 1;
    (m_grid[i]).m_minOrientVal = 1;
    (m_grid[i]).m_orintVisited = false;
    (m_grid[i]).m_maxOrientIndex = 0;
    (m_grid[i]).m_minOrientIndex = 0;
    (m_grid[i]).m_zestiVisited = false;
    (m_grid[i]).m_times_visited = 0;
    // Fit value of the current node in a tracklet.
    (m_grid[i]).m_fitValue = 0.00;
    (m_grid[i]).m_mahalanobisDist = 0.00;
    // The z values for virtual tubes are computed and stored in the grid,
    // so we do not need to reset them
    if( (m_grid[i]).m_type == GridNode::STT_TYPE_PARA ) {
      (m_grid[i]).m_xDet  = (m_grid[i]).m_x;
      (m_grid[i]).m_yDet  = (m_grid[i]).m_y;
      (m_grid[i]).m_z_Det = 0.00;
    }
    else if( (m_grid[i]).m_type == GridNode::STT_TYPE_SKEW ) {
      (m_grid[i]).m_xDet  = (m_grid[i]).m_x;
      (m_grid[i]).m_yDet  = (m_grid[i]).m_y;
      (m_grid[i]).m_z_Det = (m_grid[i]).m_z;
    }
    else if( (m_grid[i]).m_type == GridNode::VIRTUAL_NODE ) {
      (m_grid[i]).m_z_Det = (m_grid[i]).m_z;
    }
  }
  m_MVD_grid.clear();
  // Set tollerance to zero
  m_AttSpaceTolerance = 0.0;

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

void CoordGrid::ResetGridOrientations()
{
  for(size_t i = 0; i < m_grid.size(); ++i) {
    (m_grid[i]).resetNodeOrientation();
  }
  m_AttSpaceTolerance = 0.0;
}

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
