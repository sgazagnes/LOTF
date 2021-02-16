#pragma once
#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H


// Local headers
#include "pathCandidate.h"
#include "CoordGrid.h"
#include "trackObject.h"


void findEasyTracks (CoordGrid &gr, std::vector< GridNode > &Ingrid, std::vector < PathCandidate* > &tracklets, std::vector<int> &activeId, char *visited, int *candidateId);
//void fitZCoordinates(CoordGrid &hitMap, PathCandidate *trk);
void TrackZ_CoordinatesDistNorm(CoordGrid &hitMap, std::vector<TrackObject*>* TrackList = 0);
void CompZCoordinates(CoordGrid &hitMap, PathCandidate *trk);


#endif
