/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#pragma once
#ifndef COORD_GRID_CLASS_H
#define COORD_GRID_CLASS_H

// C & C++ headers.
#include <vector>

// Local headers
struct HitCoordinate;
struct GridNode;

// Write the enire grid to a ascii file
#define WRITE_GRID_TO_ASCII 0
#define START_VIRTUAL_ID 6000

bool Is_STT_SplitSkewedNode( GridNode const &node);
bool IsVirtualSplitNode( GridNode const &node);

struct CoordGrid {
  
  CoordGrid();
  
  virtual ~CoordGrid();
  CoordGrid( CoordGrid const& ot);
  CoordGrid& operator=(CoordGrid const& ot);

// Get the total number of available nodes in grid.  
  inline size_t GetNumNodes() const;

  /*
   * Init Grid data structure.
   */
  void Initialize(std::vector < GridNode > const &detNodes);
  void CorrectLayerLimit();

  /* Add node(s) to grid */
  void ExtendedGrid(std::vector < GridNode > const &Nodes);
  void AddNodeToGrid(GridNode const &Node);

  void Rearange_nodeIds();
  void FillGrid(std::vector < HitCoordinate* > const& hitcoords);
  
  // Functions to rest grid node properties
  void ResetGrid();
  void ResetGridOrientations();
  
  // Find a node with a given node ID
  int Find(int DetId) const;
  GridNode& GetNodeByID(int DetId);
  GridNode const& GetNodeByID(int DetId) const;

  /**
   * Compiles a list of ids of active neighbours and neighbours of
   * neighbours.
   *@param NodeId: Id of the input node.
   *@param ListOfNeighbours: (output) List containing the active
   * direct neighbours and neighbours of the direct neigbours. Note:
   * This list is cleared at each call to this function.
   */
  void ActiveNeighboursAndSecondOrderNeighbours(int NodeID, std::vector<int> &ListOfNeighbours) const;
  void ActiveNeighboursList(int NodeID, std::vector<int> &ListOfNeighbours,
                            unsigned int depth = 0) const;

  // MVD subsystem of the grid (not fully included yet)
  int FindMVD(int ID) const;
   
  //___________ Debug function to dump in a file
#if (WRITE_GRID_TO_ASCII > 0)
  void WriteGrid(std::string const& fileName="GridDumpFile.txt") const;
  //___________ END DEBUG FUNCTIONS     
#endif
  // Tolerance for computing the orientation based Att-spcace
  float m_AttSpaceTolerance;
  size_t firstVirtIdx;
  // Members (currently only STT tubes.
  std::vector< GridNode > m_grid;
  /*
   * FIXME FIXME We need to include MVD into the static grid geometry.
   *
   * Store MVD points separately This is for now until we find out how
   * to store full MVD and handle detector indices.
   */
  std::vector< GridNode > m_MVD_grid;
  
  //_______________________ PROTECTED members ____________
  // protected:
  
  //______________________ PRIVATE members ________________
 private:
  bool operator==(CoordGrid const& ot) const;
  bool operator>(CoordGrid const& ot) const;
  bool operator<(CoordGrid const& ot) const;
};
// Inline functions implementations
inline size_t CoordGrid::GetNumNodes() const
{
  return (m_grid.size() + m_MVD_grid.size());
}
#endif// END interface
