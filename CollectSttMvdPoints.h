/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#pragma once
#ifndef COLLECT_STT_MVD_POINTS_H
#define COLLECT_STT_MVD_POINTS_H

// C and C++ standard headers
#include <vector>

// Root headers
#include "TNtuple.h"
#include "TClonesArray.h"

// Local headers
#include "gridNode.h"
#include "CoordGrid.h"

// Local defines
#define HIT_EXCLUSION      -10000.0
#define NOT_AVAILABLE  -10000.0

void CollectSttDetecorCoords(TClonesArray const &TubeArray, std::vector < GridNode > &detNodes);

void CollectGridToTree( CoordGrid const &gr, TNtuple &out);

std::vector < std::vector<HitCoordinate*>* >*
CollectSttMvdPoints( std::vector < GridNode >& detNodes, TFile &OutFile,
		     int firstEvt = -1, int lastEvt = -1);

void WriteEventPlotsToFile(std::vector < std::vector<HitCoordinate*>* > const &evtData, TFile &OutFile);
#endif
