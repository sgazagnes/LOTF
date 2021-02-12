
#pragma once
#ifndef SIMON_FUNCTIONS_H
#define SIMON_FUNCTIONS_H



// Local headers

#include "CoordGrid.h"
#include "gridNode.h"
#include "pathCandidate.h"


/*#include "CollectSttMvdPoints.h"
#include "SttMVDEventDataReader.h"
#include "pathopen.h"
#include "hitcoordinate.h"
#include "floodingFilter.h"
#include "utilfunctions.h"
#include "trackObject.h"
#include "logc.h"
#include "queue.h"
//#include "performFilter.h"
#include "path_queue.h"
*/

typedef enum {
  DOWN = 1,
  SAME = 2,
  UP = 4
} Direction;


bool sortNeighbors(CoordGrid &gr, GridNode *currentNode, std::vector<int> &prev, std::vector<int> &same, std::vector<int> &next, std::vector<int> &virt, char *visited, int *dir);

int determineSkewed_XYPlane_new( CoordGrid &hitMap, GridNode const &VNode,
                            std::vector<int> &ListOfSkewedNodesIndex,
				 std::vector<int> &ListOfVirtualNodesIndex, char *visited);

void resetLists(char *visited, std::vector<int> &prev, std::vector<int> &same, std::vector<int> &next);

void removeIdFromNeigh(GridNode *neighNode, std::vector<int> *prevNodes, int curId);

bool areAdjacent(CoordGrid &gr, std::vector<int> *v);

double returnAngle(double x1, double x2, double x3, double y1, double y2, double y3);

double *polyFit(std::vector<double>  x, std::vector<double>  y, int n);

int returnDirection(double prev, double cur);

int fitNextId(CoordGrid &gr, PathCandidate &cand, std::vector<int> &next, int k);

double distanceBetweenTube(GridNode & tubeA, GridNode & tubeB);



void Fix_InterSector_Nodes(CoordGrid &hitMap, size_t const numSectors);

double IntersectionPointSkePar(CoordGrid const &hitMap,
				   GridNode &tubeA, GridNode &tubeB,  GridNode &out);
double IntersectionPointSkeSke(CoordGrid const &hitMap,
				   GridNode &tubeA, GridNode &tubeB,  GridNode &out);
void IntersectionPointCoord(GridNode &tubeA, GridNode &tubeB);
void addTracklets (CoordGrid &gr, PathCandidate *newCand, PathCandidate &mergeCand,  int curdir, int mergedir);
#endif
void Add_VirtualNodes(CoordGrid &hitMap, std::vector < GridNode > &VNodesLayer,  std::vector < GridNode > &VNodesSector);
double IntersectionXY(double startX1, double endX1, double startY1, double endY1, double startX2, double endX2, double startY2, double endY2);
bool LineLineIntersect( GridNode &tubeA, GridNode &tubeB, GridNode &tubeC, float &ixOut, float &iyOut, float &izOut); //Output 
