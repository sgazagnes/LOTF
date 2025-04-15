#pragma once
#ifndef PHZINTERP_H
#define PHZINTERP_H


// Local headers
#include "gridNode.h"
#include "pathCandidate.h"

void ZCoordinates(CoordGrid &gr, std::vector< GridNode > &Ingrid,std::vector < PathCandidate* > &tracklets);

#endif
