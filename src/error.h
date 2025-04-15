#pragma once
#ifndef ERROR_H
#define ERROR_H

#include <vector>

#include "CoordGrid.h"
#include "trackObject.h"
#include "pathCandidate.h"

// Structure holding per event errors
typedef struct EvtErrorStruct{
  // public:
  // Constructor
EvtErrorStruct()
  : nMCTracks(0),
    nRecoTracks(0),
    UnderMerge(0.0),
    OverMerge(0.0),
    TotalError(0.0),
    UnderMergeNorm(0.0),
    OverMergeNorm(0.0),
    TotalErrorNorm(0.0){};
  
  // Destructor
  virtual ~EvtErrorStruct(){};

  // Local variables(Per image values)
  size_t nMCTracks;// MC-Tracks
  size_t nRecoTracks;// Reconstructed
  float  UnderMerge;
  float  OverMerge;
  float  TotalError;
  // Pre-normalized values 
  float  UnderMergeNorm;
  float  OverMergeNorm;
  float  TotalErrorNorm;
  // double tr_E;
} EvtErrorStruct;
////////////////////////////////

//_____ Structure to hold error evaluation parameters _____
typedef struct TrackErrorStruct{
  // public:
  // Constructor
TrackErrorStruct()
  : isComplex(0),
    isNotmatched(0),
    MCTrackLength(0),
    RecoTrackLength(0),
    IntersectionLength(0),
    UnionLength(0),
    F1score(0),
    MC_px(0.0),
    MC_py(0.0),
    MC_pz(0.0),
    MC_a(0.0),
    MC_b(0.0),
    MC_r1(0.0), // Using MCSTTPoints
    MC_r2(0.0), // Using MCHitsCoords
    tr_a(0.0),
    tr_b(0.0),
    tr_rIsoRand(0.0),
    tr_rIsoRand16(0.0),
    tr_rIsoRand84(0.0),
    tr_rAnc(0.0),
    tr_rPts(0.0),
    tr_scattAngle(-1.),
    tr_rank(0),
    tr_isClone(0){};
  
  // Destructor
  virtual ~TrackErrorStruct(){};

  // Local variables
  int    isComplex;
  float  isNotmatched;
  float  F1score;
  size_t MatchIndex;
  size_t MCTrackLength;
  size_t RecoTrackLength;
  size_t IntersectionLength;
  size_t UnionLength;
  
  /*Curvature parameters for members.*/
  double MC_px;
  double MC_py;
  double MC_pz;
  double MC_a;
  double MC_b;
  double MC_r1;
  double MC_r2;
  
  // Reco tracklets
  double tr_a;
  double tr_b;
  double tr_rIsoRand;
  double tr_rIsoRand16;
  double tr_rIsoRand84;
  double tr_rAnc;
  double tr_rPts;
  double tr_scattAngle;
  int tr_rank;
  int tr_isClone;

  // Displacements of coordinates average
  double MeanDiffX;
  double MeanDiffY;
  double MeanDiffZ;

  // Displacements vectors
  std::vector<float> DiffX;
  std::vector<float> DiffY;
  std::vector<float> DiffZ;

} TrackErrorStruct;
////////////////////////////////


/**
 *@param
 */
bool sortbysec2(const pair<int,int> &a, 
              const pair<int,int> &b) 
{ 
    return (a.second > b.second); 
} 

//void complexSectors(CoordGrid &gr,  std::vector< int > &activeId, std::vector< int >* sectorTC);
void ComplexTracks(CoordGrid &gr,
		   std::vector< MCTrackObject* >  const *MCTracks,
		   std::vector< int > *ListIDComplex);

EvtErrorStruct* ComputeGlobalEvtErrors(CoordGrid &gr,
				       std::vector< MCTrackObject* >  const *MCTracks,
				       std::vector< std::set<int>* >  const *RecoTracks);

std::vector< TrackErrorStruct* >* ComputeErrorPerRecoTrack(CoordGrid const &hitMap,
                                                           std::vector < MCTrackObject* > const *MCTracks,
                                                           std::vector < PathCandidate* > const *RecoTracks,
							   std::vector<int >  &ListIDComplex,
							   std::vector<int >  &IDMatchesMCReco);

std::vector< TrackErrorStruct* >* PandaErrorMetric(CoordGrid const &hitMap,
						  std::vector < MCTrackObject* > const *MCTracks,
						  std::vector < PathCandidate* > const *RecoTracks);

std::vector< int > MatchBestRecoToMC( CoordGrid const &hitMap,
				      std::vector < MCTrackObject* > const *MCTracks,
				      std::vector < PathCandidate* > const *tracklets );

#endif
