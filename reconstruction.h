#pragma once
#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H


// Local headers
#include "simon_functions.h"


void findEasyTracks (CoordGrid &gr, std::vector < PathCandidate* > &tracklets, std::vector<int> &activeId, char *visited, int *candidateId);

#endif
