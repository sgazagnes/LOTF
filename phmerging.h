#pragma once
#ifndef PHMERGING_H
#define PHMERGING_H

#include "pathCandidate.h"
#include "CoordGrid.h"
#include "gridNode.h"

void mergeTracks (CoordGrid &gr, std::vector< GridNode > &Ingrid, std::vector < PathCandidate* > &tracklets,  int *candidateId );

#endif
