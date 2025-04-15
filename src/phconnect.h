#pragma once
#ifndef PHCONNECT_H
#define PHCONNECT_H


// Local headers
#include "pathCandidate.h"
#include "CoordGrid.h"
#include "trackObject.h"

void addNodesToCand (CoordGrid &gr, std::vector< GridNode > &Ingrid,  PathCandidate &cand, char *visited, std::vector<int> &v);
void findEasyTracks (CoordGrid &gr, std::vector< GridNode > &Ingrid, std::vector < PathCandidate* > &tracklets,     std::vector<pair<int, unsigned short>> idToProcess, char *visited, int *candidateId);


#endif
