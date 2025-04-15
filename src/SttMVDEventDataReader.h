/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#pragma once
#ifndef Stt_MVD_EVENT_DATA_READER_H
#define Stt_MVD_EVENT_DATA_READER_H

// C and C++ standard headers
#include <string>
#include <vector>

// Local headers
class HitCoordinate;
#include "gridNode.h"

typedef enum GROUP_HITS {
  BY_EVENT     = 0,
  BY_TIMESTAMP = 1,
} GroupingOrder;

void SttMVDParametersTreeFileRead( std::string const &InFileName= "Infile.root",// Inputfile name
				   std::string const &EvenTree = "Treeinput", // Input event tree name
				   int firstEvt = -1, int lastEvt = -1,
				   GroupingOrder grouping = BY_EVENT,
				   std::vector < std::vector<HitCoordinate*>* > *OutputContainer = 0// Output parameter.
				   );

void GroupeByEvent( int firstEvt, int lastEvt, 
		    std::vector <HitCoordinate*> const &HitContainer,
		    std::vector < std::vector<HitCoordinate*>* > &Output);

void GroupByTimeStamp( int firstEvt, int lastEvt, float tZero, float offset,
		       std::vector <HitCoordinate*> const &HitContainer,
		       std::vector < std::vector<HitCoordinate*>* > &Output);

std::vector <GridNode>* ReadGridGeometry( std::string const &InFileName,// Inputfile name
					  std::string const &GeoTree // Input node geometry
					  );
#endif
