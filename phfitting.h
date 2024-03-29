#pragma once
#ifndef PHFITTING_H
#define PHFITTING_H

#include "pathCandidate.h"
#include "CoordGrid.h"


int fitNextId(CoordGrid &gr, std::vector< GridNode > &Ingrid, PathCandidate &cand, std::vector<int> &next, int k);
void fittingPhase(CoordGrid &gr, std::vector< GridNode > &Ingrid, std::vector < PathCandidate* > &tracklets, std::vector<pair<int, unsigned short>> idToProcess, char *visited, int **toMergeWith);

#endif
