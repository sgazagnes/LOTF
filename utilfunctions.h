/**************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#pragma once
#ifndef PATH_UTIL_FUNCTIONS_H
#define PATH_UTIL_FUNCTIONS_H

#include <utility>
#include <vector>
#include <string>
#include <algorithm>

// Local include and declarations.
#include "auxiliaryfunctions.h"
struct CoordGrid;
struct GridNode;
struct TNtuple;
struct PathCandidate;

// Local zero compare constants defines
#define LOCAL_ZERO_COMPARE_EPSILON  1.0e-8
#define LOCAL_PND_TRACKING_EPSILON  0.001
// When fitting a circle to data points in a path
#define LOCAL_CIRCLE_FIT_TOLERANCE 1e-04//1e-06
#define MAX_NUMBER_OF_FIT_ITERATION 5000//1000, 10000;//512;//256;

//------------------------------------------------
// Used to merge short candidates(max acceptable distance)
#define SHORT_CANDIDATES_MERGE_DISTANCE 3.0 //4.0

// Used for pre-selection of possible merge candidates.
#define CAND_MERGE_DISTANCE_SELECTION_TOLERANCE 2.0 //3.0, 4.0~6.0

// Use in "selectcandidatesMatchNodes" to allow max laye diff
#define MAX_ALLOWED_LAYER_DIFF_HEAD_TAIL 5

// Use to pre-select candidates for curvature based merging
#define CANDIDATE_DISTANCE_CURVATURE_MERGE  6.0 //6.0
#define CANDIDATE_MAHALANOBIS_DISTANCE_CURVATURE_MERGE 200//200

#define MIN(a,b)  ((a<=b) ? (a) : (b))
#define MAX(a,b)  ((a>=b) ? (a) : (b))

//________________ candidateDistobject __________________________________
struct candidateDistobject{
  /** Constructors */
candidateDistobject()
: distance(-1.0), FirstCandNodeId(-1), SecondCandNodeId(-1),
    indexFirstCand(0), indexSecondCand(0), FirstCandidateId(0), SecondCandidateId(0),
    isValid(true), m_curvaturePar(CurvatureParameters()),
    mahalanobis_FirstNode_SecondCand(0.0), mahalanobis_SecondNode_FirstCand(0.0),
    mahalanobisDist(0.0), curvDiff(1000.0), radiusdiff(1000.0), centreDist(1000.0)
  {};
  explicit candidateDistobject(float d, int FCandnId, int SCandnId, size_t idxFCand, size_t idxSCand)
    : distance(d), FirstCandNodeId(FCandnId),SecondCandNodeId(SCandnId),
    indexFirstCand(idxFCand),indexSecondCand(idxSCand),
    FirstCandidateId(0), SecondCandidateId(0), isValid(true),
    m_curvaturePar(CurvatureParameters()),
    mahalanobis_FirstNode_SecondCand(0.0), mahalanobis_SecondNode_FirstCand(0.0),
    mahalanobisDist(0.0), curvDiff(1000.0), radiusdiff(1000.0), centreDist(1000.0)
  {};
  /** Copy const.*/
candidateDistobject(candidateDistobject const &ot)
: distance(ot.distance),FirstCandNodeId(ot.FirstCandNodeId),
    SecondCandNodeId(ot.SecondCandNodeId), indexFirstCand(ot.indexFirstCand),
    indexSecondCand(ot.indexSecondCand),
    FirstCandidateId(ot.FirstCandidateId), SecondCandidateId(ot.SecondCandidateId),
    isValid(ot.isValid), m_curvaturePar(ot.m_curvaturePar),
    mahalanobis_FirstNode_SecondCand(ot.mahalanobis_FirstNode_SecondCand),
    mahalanobis_SecondNode_FirstCand(ot.mahalanobis_SecondNode_FirstCand),
    mahalanobisDist(ot.mahalanobisDist), curvDiff(ot.curvDiff),
    radiusdiff(ot.radiusdiff), centreDist(ot.centreDist)
    
  {};
  // Assignment
  candidateDistobject& operator=(candidateDistobject const &ot)
  {
    if(this != &ot) {
      this->distance = ot.distance;
      this->FirstCandNodeId = ot.FirstCandNodeId;
      this->SecondCandNodeId = ot.SecondCandNodeId;
      this->indexFirstCand = ot.indexFirstCand;
      this->indexSecondCand = ot.indexSecondCand;
      this->isValid = ot.isValid;
      this->m_curvaturePar = ot.m_curvaturePar;      
      this->FirstCandidateId = ot.FirstCandidateId;
      this->SecondCandidateId = ot.SecondCandidateId;
      this->mahalanobis_FirstNode_SecondCand = ot.mahalanobis_FirstNode_SecondCand;
      this->mahalanobis_SecondNode_FirstCand = ot.mahalanobis_SecondNode_FirstCand;
      this->mahalanobisDist = ot.mahalanobisDist;
      this->curvDiff = ot.curvDiff;
      this->radiusdiff = ot.radiusdiff;
      this->centreDist = ot.centreDist;
    }
    return (*this);
  };
  /** Destructor */
  virtual ~candidateDistobject(){};

  /** Compare members functions .(ordening) */
  bool operator>(candidateDistobject const &ot) const
  {
    return (this->distance > ot.distance);
  };
  bool operator<(candidateDistobject const &ot) const
  {
    return (this->distance < ot.distance);
  };
  
  inline friend bool lessThanCurve(candidateDistobject const &lc, candidateDistobject const &rc);
  inline friend bool greaterThanCurve(candidateDistobject const &lc, candidateDistobject const &rc);
  inline friend bool lessthanMahalanobis(candidateDistobject const &lc, candidateDistobject const &rc);
  inline friend bool greaterThanMahalanobis(candidateDistobject const &lc, candidateDistobject const &rc);
  inline friend bool lessThanFirstCandId(candidateDistobject const &lc, candidateDistobject const &rc);
  inline friend bool greaterThanFirstCandId(candidateDistobject const &lc, candidateDistobject const &rc);
  
  // Members (all public, but we don't care, :) )
  float  distance;
  int    FirstCandNodeId;
  int    SecondCandNodeId;
  size_t indexFirstCand;
  size_t indexSecondCand;
  unsigned int FirstCandidateId;
  unsigned int SecondCandidateId;
  bool   isValid;
  CurvatureParameters m_curvaturePar;
  float mahalanobis_FirstNode_SecondCand;
  float mahalanobis_SecondNode_FirstCand;
  float mahalanobisDist;
  float curvDiff;
  float radiusdiff;
  float centreDist;
};// END class definition
//____________ Function implementations ____
bool lessThanFirstCandId(candidateDistobject const &lc, candidateDistobject const &rc)
{
  return(lc.FirstCandidateId < rc.FirstCandidateId);
}
bool greaterThanFirstCandId(candidateDistobject const &lc, candidateDistobject const &rc)
{
    return(lc.FirstCandidateId > rc.FirstCandidateId);
}
bool lessthanMahalanobis(candidateDistobject const &lc, candidateDistobject const &rc)
{
  return (lc.mahalanobisDist < rc.mahalanobisDist);
}
bool greaterThanMahalanobis(candidateDistobject const &lc, candidateDistobject const &rc)
{
  return (lc.mahalanobisDist > rc.mahalanobisDist);
}
bool lessThanCurve(candidateDistobject const &lc, candidateDistobject const &rc)
{
  CurvatureParameters const &lpar = lc.m_curvaturePar;
  CurvatureParameters const &rpar = rc.m_curvaturePar;
  return (lpar.m_E < rpar.m_E);
}
bool greaterThanCurve(candidateDistobject const &lc, candidateDistobject const &rc)
{
  CurvatureParameters const &lpar = lc.m_curvaturePar;
  CurvatureParameters const &rpar = rc.m_curvaturePar;
  return (lpar.m_E > rpar.m_E);
}
//____________________ END candidateDistobject ___________________________
/**
 * Data structure to represent a point in 3D space
 */
struct point3D{
  // public:
point3D()
: m_x(0.000),
    m_y(0.000),
    m_z(0.000)
  {};
  
explicit point3D(float x, float y, float z)
  : m_x(x),
    m_y(y),
    m_z(z)
  {};
  
point3D(point3D const &ot)
  : m_x(ot.m_x),
    m_y(ot.m_y),
    m_z(ot.m_z)
  {};

  virtual ~point3D(){};

  point3D& operator=(point3D const &ot)
  {
    if(this != &ot) {
      this->m_x = ot.m_x;
      this->m_y = ot.m_y;
      this->m_z = ot.m_z;
    }
    return (*this);
  };

  // member parameters.
  float m_x;
  float m_y;
  float m_z;

private:
  bool operator==(point3D const &ot) const;
  bool operator>(point3D const &ot) const;
  bool operator<(point3D const &ot) const;
};
//---------------------------------------------------------
struct TubeLayerPairProperty{
  //public:
  // Constructor
TubeLayerPairProperty()
: firstNodeID(0),
    secondNodeID(0),
    firstNodeIndex(0),
    secondNodeIndex(0),
    isValid(true)
  {};
  // Copy constructor
  TubeLayerPairProperty(TubeLayerPairProperty const &ot)
  : firstNodeID(ot.firstNodeID),
    secondNodeID(ot.secondNodeID),
    firstNodeIndex(ot.firstNodeIndex),
    secondNodeIndex(ot.secondNodeIndex),
    isValid(ot.isValid)
  {};
  // Destructor
  virtual ~TubeLayerPairProperty(){};
  
  // Operators
  TubeLayerPairProperty& operator=(TubeLayerPairProperty const &ot)
  {
    // check for self-assignment
    if(&ot != this){
      this->firstNodeID     = ot.firstNodeID;
      this->secondNodeID    = ot.secondNodeID;
      this->firstNodeIndex  = ot.firstNodeIndex;
      this->secondNodeIndex = ot.secondNodeIndex;
      this->isValid         = ot.isValid;
    }
    return *this;
  };
  bool operator==(TubeLayerPairProperty const &ot) const
  {
    return  (( (firstNodeID == ot.firstNodeID)  && ( secondNodeID == ot.secondNodeID) ) ||
	     ( (firstNodeID == ot.secondNodeID) && ( secondNodeID == ot.firstNodeID) ));
  };

  friend bool equalTubePairProperty( TubeLayerPairProperty const &l,
				     TubeLayerPairProperty const &r)
  {
    return (( (l.firstNodeID == r.firstNodeID)  && ( l.secondNodeID == r.secondNodeID) ) ||
	    ( (l.firstNodeID == r.secondNodeID) && ( l.secondNodeID == r.firstNodeID) ));
  };

  // Local members
  int    firstNodeID;
  int    secondNodeID;
  size_t firstNodeIndex;
  size_t secondNodeIndex;
  bool   isValid;
  //protected:
  //private:
};
//_____________________________
/**
 */
float PointDistance(point3D const &a, point3D const &b);
// _____________________________________________________________
/**
 * Polarout <R, theta>: To compute the value, the function takes into
 * account the sign of both arguments in order to determine the
 * quadrant.
 */
double Cartesian_To_Polar(float const x, float const y,
                          std::pair<float,float>& polarOut,
                          bool useSign = true);

/**
 * Computes angle of a line between two given points in the space and
 * the x-axis. NOte that the output is in the range of [0, pi).
 *
 *@param x1,y1 define x and y coordinates of the first point.
 *@param x2,y2 define x and y coordinates of the second point.
 *
 *@param Rad_deg Holds the output in degrees and
 * radiants. Rad_deg.first = radiant, Rad_deg.second = degrees.
 */
void ComputeSlope(float const x1, float const y1, float const x2, float const y2,
                  std::pair<float,float>& Rad_deg, float epsilon = -100.00);
/**
 * Compute the virtual nodes between the layers with different slope.
 *
 *@param hitMap: Input grid
 *@param VNodes: Output, list of determined virtual nodes
 */
void Compute_Add_VirtualNodes_Neigbor(CoordGrid &hitMap, std::vector < GridNode > &VNodes);
// Old method
void Compute_Add_VirtualNodes_Neigbor2(CoordGrid &hitMap, std::vector < GridNode > &VNodes);

/**
 *@param hitMap, the STT tubes graph.
 *@param NodePairSet: The list of pairs of nodes for which a virtual
 * node can be defined. The duplicated paire (equivalt) are marked as
 * invalid.
 */
void FindNodeBetweenLayerNodePairs(CoordGrid const &hitMap,
				   std::vector< TubeLayerPairProperty > &NodePairSet);

/**
 * Determine the virtual intersection point between the two input
 * nodes (based on their neighbourhood relations).
 */
bool IntersectionPoint_NeigborList(CoordGrid const &hitMap, GridNode &a, GridNode &b,
				   GridNode &out);
/**
 * Determine the virtual intersection point between the two input
 * nodes (based on the sector where the are placed).
 */
bool IntersectionPoint_SectorList (CoordGrid const &hitMap, GridNode &a, GridNode &b,
				   GridNode &out);
/**
 * Compute virtual nodes between different sectors. This method is not
 * implemented correctly yet. Because it is not yet clear how to
 * compute the position of the nodes between two different
 * neighbouring sectors.
 *@param hitMap: Input grid
 *@param VNodes: List of virtual nodes between different sectors.
 */
void Compute_Virtual_InterSector_Nodes(CoordGrid &hitMap, size_t const numSectors,
				       std::vector < GridNode > &VNodes);

//void Fix_InterSector_Nodes(CoordGrid &hitMap, size_t const numSectors);
/**
 * Given two tubes first and second in two different neighbouring
 * sectors. Determinde their virtual intersection point.
*/
//double distanceBetweenTube(GridNode & tubeA, GridNode & tubeB) ;

void InterSectorPoints (GridNode const& first, GridNode const& second, GridNode &outPut);

/**
 * Stores a given grid of tubes into a ROOT TNtuple.
 *@param Nodes The vector containing input tupes.
 *@return TNtuple containing coordinates of tubes.
 */
TNtuple* GridToNtuple(std::vector < GridNode > const &Nodes, std::string const &name = "TestTuple");

//______________________________ Helper functions __________________________
/**
 * Temporary function to fix the neighboring issue in the current grid
 * implementation. If the implementation is fixed, we will not need
 * this one and can be remoeved.
 */
void fixNeighboring(CoordGrid &hitMap);

/**
 * Isolate nodes that are sector or layer boudaries.
 *@param hitMap Input graph.
 *@param sectorLimits List of section limits.
 *@param layerLimits  List of layer limits.
*/
void isolateSectorAndLayerLimits(CoordGrid const &hitMap, TNtuple &sectorLimits, TNtuple &layerLimits);

/**
 * Determine the curvature for all candidate paths.
 *@param hitMap Coordiantes graph(input)
 *@param Candidate  Candidates(input & output).
 */
void computePathCurvature(CoordGrid const &hitMap, PathCandidate &Candidate);
void computeCurvatureForAllPaths(CoordGrid const &hitMap, std::vector<PathCandidate*> &Candidates);
void ComputeCurvatureForListOfnodes( CoordGrid const &hitMap, std::vector<int> const &nodes,
                                     CurvatureParameters &OutPars);

void updateCandidateNodesFitVal (CoordGrid  &hitMap, PathCandidate const &Candidate);
void updateNodeFitValForAllPaths(CoordGrid &hitMap, std::vector<PathCandidate*> const &Candidates);

/** If a given path candidate is invalid. */
bool isInValidPath(PathCandidate const *candX);

void removeInvalidSubPaths(std::vector < PathCandidate* > &CandidesContainer);

void markAsInvalidLengthBased(std::vector < PathCandidate* > &CandidesContainer, size_t minLength = 0);

void mergeMarkOverlappingAsInvalid(std::vector < PathCandidate* > &AllCandides, int minOverlap = -1);

std::vector<int>* findIntersectionElements(PathCandidate const &cand1, PathCandidate const &cand2);

void markMergedAsInvalid(std::vector < PathCandidate* > &AllCandides, unsigned int minlevel = 1);

void markInvalidLevelBased(std::vector < PathCandidate* > &AllCandides, unsigned int level = 1);

void resetValidy(std::vector < PathCandidate* > &AllCandides);

void MarkInvalidWithChilderen(std::vector < PathCandidate* > &AllCandides, unsigned int numChilderen = 1);

/**
 * Fits a circle to the member points of the current candidate
 */
int CircleFit( CoordGrid const &hitMap, std::vector<int> const &MemberIdVector,
               CurvatureParameters &curvature);
/**
 * Given the coordinates of a set of points in 2D, fit a circle
 * through those points.
 */
int CircleFit( std::vector<point3D> const &points, CurvatureParameters &curvature);

/*
 * Path candidates merging functions and procedures.
 */
void updatePathCandidateHeadAndTailNodes(PathCandidate &cand);
void updateAllCandidateHeadAndTailNodes(std::vector<PathCandidate*> &Candidates);

void updatePathCandidateHeadAndTailLayer(CoordGrid const &hitMap, PathCandidate &cand);
void updateAllCandidateHeadAndTailLayers(CoordGrid const &hitMap, std::vector<PathCandidate*> &Candidates);


// Candidate CoVar and Mean
void computePathCovariance(CoordGrid const &hitMap, std::vector<int> const &nodeIds,
                           std::vector< float > &mean, std::vector< std::vector<float> > &covMatrix,
                           std::vector< std::vector<float> > &invCovMatrix);
void computeCandCovariance(CoordGrid const &hitMap, PathCandidate &cand);
void computeCovMatForAllCandidates(CoordGrid const &hitMap, std::vector < PathCandidate* > &CandidesContainer);
void computeCandidateEigenValVectors(PathCandidate &candidate);
void computeEigenValVectorsAllCandidates(std::vector < PathCandidate* > &allCandiddates);

/**
 * Find shortest distance between all candidate pairs (Euclidean distance based).
 */
void candidateShortestdistPairs( CoordGrid const &hitMap, PathCandidate const &cand,
                                 std::vector < PathCandidate* > const &CandContainer,
                                 std::vector<candidateDistobject> &OutPutCandidatePairs,
                                 float maxDistThreshold = 0.0, bool useMahalanobis = false);

/**
 * Compile a list of all candidates with a path_length <= to
 * maxDistThreshold from the input candidate. Possible nodes are
 * selected based on the available paths in the detector graph and are
 * selected based on either Euclidean or Mahalanobis
 * distance.
 *@param hitMap Input detector graph.
 *@param cand Input path candidate.
 *@param CandContainer The list of all available (merge)path-candidates.
 *@param OutPutCandidatePairs The list of pairs of candidates that
 * meet the selection criterion.
 *@param maxDistThreshold The selection distance threshold.
 *@param useMahalanobis  If use and evaluate Mahalanobis-distance to select.
 *@param maxNeighbourDepth Maximum allowed number of layers to skip if needed
 * to find the next active node in the graph
 */
void candidateShortestdistPairsGraphPathBased( CoordGrid const &hitMap, PathCandidate const &cand,
                                               std::vector < PathCandidate* > const &CandContainer,
                                               std::vector<candidateDistobject> &OutPutCandidatePairs,
                                               float maxDistThreshold = 0.0, bool useMahalanobis = false,
                                               unsigned int maxNeighbourDepth = 0);

void selectcandidatesMatchNodes( CoordGrid const &hitMap, PathCandidate const &FirstCand,
                                 PathCandidate const &SecondCand, candidateDistobject &output);

/**
 * Find the number of ends for the input candidate.
 */
size_t candidateNumOfEndStops(CoordGrid const &hitMap, PathCandidate const &inputCand);
size_t pathNumOfEndStops(CoordGrid const &hitMap, std::set<int> const &nodeList);
size_t numEndStopsIfMerged(CoordGrid const &hitMap, PathCandidate const &p1, PathCandidate const &p2);
/**
 * Mahalanobis distance based functions.
 */
// Determine Mahalanobis distance for all nodes in a given path
void ComputeMahalanobisDist(CoordGrid &hitMap, PathCandidate const &cand);
// Mahalanobis distance between all candidate pairs.
void ComputeMahalanobisDistForAllCandidates(CoordGrid &hitMap, std::vector < PathCandidate* > const &CandidesContainer);
// Mahalanobis distance between a candidate to a node (normally head
// or tail of another candidate)
float ComputeMahalanobisDistNodeToCandidate(PathCandidate const &cand, GridNode const &node);

// Merge sigletons
void mergeSingletons( CoordGrid const &hitMap, std::set<int> &SetOfNodesID, std::vector<PathCandidate*> &Candidates);
void cleanupDuplicateSingletons(std::vector<PathCandidate*> &Candidates);

// Merge short paths
void mergeShortPaths( CoordGrid const &hitMap, std::vector < PathCandidate* > &CandidesContainer, size_t minL);

// Mark all paths with (length < minL)
void markAllShortPaths(std::vector<PathCandidate*> &Candidates, size_t minL);

/**
 * Merge based on Mahalanobis distance
 */
unsigned int mergeCandidatesMahalanobisDist( CoordGrid const &hitMap, std::vector < PathCandidate* > &AllCandides,
                                             float distThreshold = 2.0, unsigned int firstOutputCandidateId = 0,
                                             bool graphBasedSelect = true, bool useMahalanobisNodeDist = false,
                                             size_t maxNumLayerDepth = 3);
/**
 * The curvature value has a higher priority in selecting pairs to
 * merge.
 */
unsigned int mergeCandidatesCurvaturePriority( CoordGrid const &hitMap, std::vector < PathCandidate* > &CandidesContainer,
                                               float distThreshold = 2.0, unsigned int firstOutputCandidateId = 0,
                                               bool useMahalanobisSelection = false);

/**
 * Merge candidates for which the curvature value-difference is below a given threshold.
 */
unsigned int mergeCandidatesCurvatureCompat( CoordGrid const &hitMap, std::vector < PathCandidate* > &CandidesContainer,
                                             float differenceThreshold = 1.0, unsigned int firstOutputCandidateId = 0);
/**
 * Compile list of all candidates with a compatible curvature ((1/R1 - 1/R2) < threshold).
 */
void candidateCurvCompatPairList(CoordGrid const &hitMap, PathCandidate const& inputCand,
                                 std::vector < PathCandidate* > const &CandContainer,
                                 std::vector<candidateDistobject> &OutCandPairs,
                                 float threshold = 1.0);
#endif// END of interface definition.
