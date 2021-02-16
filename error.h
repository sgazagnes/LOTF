#pragma once
#ifndef ERROR_H
#define ERROR_H

#include <vector>

#include "CoordGrid.h"
#include "trackObject.h"

//_____ Structure to hold error evaluation parameters _____
typedef struct MCMatchingErrorStruct{
  // public:
  // Constructor
MCMatchingErrorStruct()
:   isNotmatched(0),
    NumberOfMCTracks(0),
    NumberOfTracklets(0),
    Error_underMerge(0.0),
    Error_overMerge(0.0),
    TotalError(0.0),
    Error_underMergeNorm(0.0),
    Error_overMergeNorm(0.0),
    TotalErrorNorm(0.0),
    BestMatchMCLength(0),
    CurrentTrackLength(0),
    MCMinCurrentLength(0),
    CurrentMinMCLength(0),
    MC_a(0.0),
    MC_b(0.0),
    MC_r(0.0),
    MC_E(0.0),
    tr_a(0.0),
    tr_b(0.0),
    tr_r(0.0),
    tr_E(0.0){};
  
  // Destructor
  virtual ~MCMatchingErrorStruct(){};

  // Local variables(Per image values)
  float  isNotmatched;//if it is a no matching error. Total error.
  size_t NumberOfMCTracks;// MC-Tracks
  size_t NumberOfTracklets;// Reconstructed
  int Complex;
  float  Error_underMerge;
  float  Error_overMerge;
  float  TotalError;
  // Pre-normalized values (total error does not need to be
  // normalized)
  float  Error_underMergeNorm;
  float  Error_overMergeNorm;
  float  TotalErrorNorm;
  float Jacardsingle;
  float Jacardaverage;
  // Variables for per track error evaluations.
  size_t BestMatchMCLength;// Length of best MC match
  size_t CurrentTrackLength;// Length of current tack
  size_t MCMinCurrentLength;//R_j - T_j
  size_t CurrentMinMCLength;//T_j - R_k
  /*Curvature parameters for members.*/
  // MC tracks
  double MC_a;
  double MC_b;
  double MC_r;
  double MC_E;
  // Reco tracklets
  double tr_a;
  double tr_b;
  double tr_r;
  double tr_E;
} MCMatchingError;
////////////////////////////////
/**
 *@param
 */

void complexSectors(CoordGrid &gr,  std::vector< int > &activeId, std::vector< int >* sectorTC);
void complexTracks(CoordGrid &gr,  std::vector< MCTrackObject* >  const *MCTracks, std::vector< int > *idComplex);

MCMatchingError* MatchMCTracksWithConnectedComponents(std::vector< MCTrackObject* >  const *MCTracks,
						      std::vector< std::set<int>* >  const *connectedComp);

MCMatchingError* MatchComplexMCTracks(CoordGrid &gr,
				      std::vector< MCTrackObject* >  const *MCTracks,
				      std::vector< std::set<int>* >  const *connectedComp,
				      std::vector<int >  &sectorList );

std::vector< MCMatchingError* >* MatchPerTrackWithMCTracks(CoordGrid const &hitMap,
                                                           std::vector < MCTrackObject* > const *MCTracks,
                                                           std::vector< std::set<int>* > const *connectedComp,
							   std::vector<int >  &idComplex);

#endif
