
#pragma once
#ifndef SIMON_FUNCTIONS_H
#define SIMON_FUNCTIONS_H



// Local headers
#include "auxiliaryfunctions.h"
#include "CollectSttMvdPoints.h"
#include "SttMVDEventDataReader.h"
#include "pathopen.h"
#include "hitcoordinate.h"
#include "floodingFilter.h"
#include "utilfunctions.h"
#include "trackObject.h"
#include "logc.h"
#include "queue.h"
//#include "performFilter.h"
#include "pathCandidate.h"
#include "path_queue.h"


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


double *polyFit(std::vector<double>  x, std::vector<double>  y);

int returnDirection(double prev, double cur);

int fitNextId(CoordGrid &gr, std::vector<double> &x, std::vector<double> &y, std::vector<int> &next, int method);

double distanceBetweenTube(GridNode & tubeA, GridNode & tubeB);

void Fix_InterSector_Nodes(CoordGrid &hitMap, size_t const numSectors);


#endif
