#pragma once
#ifndef PHCONNECT_H
#define PHCONNECT_H


// Local headers
#include "pathCandidate.h"
#include "CoordGrid.h"
#include "trackObject.h"


void findEasyTracks (CoordGrid &gr, std::vector< GridNode > &Ingrid, std::vector < PathCandidate* > &tracklets,     std::vector<pair<int, unsigned short>> idToProcess, char *visited, int *candidateId);


#endif
