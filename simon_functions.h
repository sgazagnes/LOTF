
#pragma once
#ifndef SIMON_FUNCTIONS_H
#define SIMON_FUNCTIONS_H



// Local headers

#include "CoordGrid.h"
#include "gridNode.h"
#include "pathCandidate.h"

#define MIN(a,b)  ((a<=b) ? (a) : (b))
#define MAX(a,b)  ((a>=b) ? (a) : (b))

#define pdd pair<double, double> ;


typedef enum {
  DOWN = 1,
  SAME = 2,
  UP = 4
} Direction;

typedef enum {
  ONGOING = 1,
  TOMERGE = 2,
  FINISHED = 3
} Status;

bool sortbysec(const pair<int,unsigned short> &a, 
	       const pair<int,unsigned short> &b);
//bool sortNeighbors(CoordGrid &gr, GridNode *currentNode, std::vector<int> &prev, std::vector<int> &same, std::vector<int> &next, std::vector<int> &virt, char *visited, int *dir);
bool sortNeighbors(CoordGrid &gr, GridNode *currentNode,  PathCandidate &cand, std::vector<int> &prev, std::vector<int> &same, std::vector<int> &next, std::vector<int> &virt, char *visited, int *dir);
int determineSkewed_XYPlane_new( CoordGrid &hitMap, GridNode const &VNode,
                            std::vector<int> &ListOfSkewedNodesIndex,
				 std::vector<int> &ListOfVirtualNodesIndex, char *visited);

void resetLists(char *visited, std::vector<int> &prev, std::vector<int> &same, std::vector<int> &next);

void removeIdFromNeigh(GridNode *neighNode, std::vector<int> *prevNodes, int curId);

//bool areAdjacent(CoordGrid &gr, std::vector<int> *v);
bool areAdjacent(CoordGrid &gr,  std::vector< GridNode > &Ingrid, std::vector<int> *v);

double returnAngle(double x1, double x2, double x3, double y1, double y2, double y3);

double *polyFit(std::vector<double>  x, std::vector<double>  y, int n);

int returnDirection(double prev, double cur);


double distanceBetweenTube(GridNode & tubeA, GridNode & tubeB);

bool PointsLineIntersectLive( GridNode &tubeC, float x1, float x2, float y1, float y2); //Output 
bool PointsLineIntersectFinal( GridNode &tubeC, float x1, float x2, float y1, float y2); //Output 


void Fix_InterSector_Nodes(CoordGrid &hitMap, size_t const numSectors);

double IntersectionPointSkePar(CoordGrid const &hitMap,
				   GridNode &tubeA, GridNode &tubeB,  GridNode &out);
double IntersectionPointSkeSke(CoordGrid const &hitMap,
				   GridNode &tubeA, GridNode &tubeB,  GridNode &out);
void IntersectionPointCoord(GridNode &tubeA, GridNode &tubeB);
void Add_VirtualNodes(CoordGrid &hitMap, std::vector < GridNode > &VNodesLayer,  std::vector < GridNode > &VNodesSector);
double IntersectionXY(double startX1, double endX1, double startY1, double endY1, double startX2, double endX2, double startY2, double endY2);
bool LineLineIntersect( GridNode &tubeA, GridNode &tubeB, GridNode &tubeC, float &ixOut, float &iyOut, float &izOut); //Output 
#endif
