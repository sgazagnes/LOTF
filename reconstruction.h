#pragma once
#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H


// Local headers
#include "pathCandidate.h"
#include "CoordGrid.h"
#include "trackObject.h"

void fittingPhase(CoordGrid &gr, std::vector< GridNode > &Ingrid, std::vector < PathCandidate* > &tracklets, std::vector<pair<int, unsigned short>> idToProcess, char *visited, int **sayYes);

//void fitZCoordinates(CoordGrid &hitMap, PathCandidate *trk);
void TrackZ_CoordinatesDistNorm(CoordGrid &hitMap, std::vector<TrackObject*>* TrackList = 0);
void CompZCoordinates(CoordGrid &hitMap, PathCandidate *trk);


#endif
