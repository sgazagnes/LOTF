/**************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#pragma once
#ifndef PATH_OPEN_H
#define PATH_OPEN_H

#include <cstdlib>
#include <vector>
#include <utility>
#include <set>

struct CoordGrid;
struct GridNode;
struct TrackObject;

/**
 * Performs path opening on a give detector grid.
 *@param hitMap The detector (STT) map containing the detector hits.
 *@param length The minimum path length required. The shorter paths
 * are removed from the map.
*/
void pathopen( CoordGrid &hitMap, size_t length);

/**
 * Performs area opening on a give detector grid.
 *@param hitMap The detector (STT) map containing the detector hits.
 *@param area The minimum required area (count of pixels). The smaller
 * areas are removed from the map.
 *@return The output vector contains lists of participating detector
 * indices participating in different connected components with an
 * area >= to "area".
*/
std::vector< std::set<int>* >* areaOpen( CoordGrid &hitMap, size_t area);

/**
 *@param hitMap The map containing fired detectors.
 *@param node The tube(pixel) from which the path is grown.
 *@param orintation The angle of the path to the x-axis(rad).
 *@param tol Tollerance on angle deviation.
 *@param Leng_Radius
 *@param ListOfTubes
 */
void PathOpenTransform(CoordGrid &hitMap, GridNode &node,
		       float orintation, float  tol,
		       std::pair<size_t, float> &Leng_Radius,
		       std::vector<int>& ListOfTubes,
		       size_t gapSize = 0);

float DetermineAngleNthNeighDepth(CoordGrid &hitMap, GridNode &node,
				  size_t depthOflevel, std::vector<int> &NeighborsTubesList,
				  bool usePCA = false);

/**
 * Compute the atrribute space based on the number of slices. Assignes
 * one or more orientation value to each participating pixel in the
 * image (detector map).
 *@param hitMap The detector (STT) map containing the detector hits.
 *@param NumSlices Number of orientations in the space. The
 *@param lambda A tunable parameter. Orientation selectivity
 *@param tol Tolerance orientations are determined between [0, pi)
 * range.
 *@return number of orientations.
 */
// We need to allow gaps between the layers and neighbors in order to
// merge the segmented tracks.
size_t computeSoftOrienAttributes(CoordGrid &hitMap, float lambda = 1.1,
				  float tol = 0.001, size_t gapSize = 0,
				  std::vector<float> const *OrientationList = 0);

/** One layer at the time is followed.
 *@param tol: tol is tolerance used for finding (matching)
 * orientations in the orientations list of a given node.
 */
void compOrientAttLayerBased_Local(CoordGrid &hitMap, float lambda = 1.0,
				   float tol= 0.001, size_t gapSize = 0,
				   bool OuterToInner = true);
/**
 *@param tol: tol is tolerance used for finding (matching)
 * orientations in the orientations list of a given node.
 */
void compOrientAttLayerBased_Alt(CoordGrid &hitMap, float tol= 0.001, bool OuterToInner = true);

/**
 *  Process Skewed in XY-Plane. This procedure corrects the
 *  xy-coordinates of the skewed tubes based on the neighbouring
 *  virtual nodes in the graph. It is used in layer-based orientation
 *  space determination.
 *@param hitMap: The full detector graph.
 *@param VNode: The virtual node from which to start.
 *@param SkewedNodesIdx: Index list of skewed involved nodes(output).
 *@param VirtualNodesIdx: Index list of virtual involved nodes(output).
 *@param OuterToInner: Determines the direction of processing.
 *@return ID of the last visited virtual node.
 */
int determineSkewed_XYPlane ( CoordGrid &hitMap, GridNode const &VNode, std::vector<int> &SkewedNodesIdx,
                              std::vector<int> &VirtualNodesIdx, bool OuterToInner = true);

/**
 * Perform attribute space connected component analysis.
 *@param hitMap Detector map containing all available grid nodes
 * (detectors).
 *@param MinPathOrint Minimum size of the connected component (Threshold).
 *@return List of the sets of indices of all detectors participating
 * in each connected component.
 */
std::vector< std::set<int>* >* AttSpaceConnectedComp( CoordGrid &hitMap, size_t MinResponce = 2);

std::vector< std::set<int>* >* AttSpaceConCompLayerBasedHard( CoordGrid &hitMap, size_t MinResponce = 2);

std::vector< std::set<int>* >* ConCompLayerBasedLocalOrient( CoordGrid &hitMap, float tol = 0.001);

std::vector< std::set<int>* >* AttSpaceConCompLayerBasedSoft( CoordGrid &hitMap, size_t MinResponce = 2);
/***/
//std::vector< std::set<int>* >* AttSpaceConCompLayerBasedSoft( CoordGrid &hitMap);

/**
 * Union find based algorithm.
 */
std::vector< std::set<int>* >* AttSpaceConnectedComp_Union( CoordGrid &hitMap, size_t MinResponce = 2);

/**
 *
 */
void MergeConnectedComponents(CoordGrid const &hitMap, std::vector< std::set<int>* >* stt_component = 0);

/**
 *
 */
void FindCompatibleMVD(CoordGrid &hitMap, GridNode const &node, float angle,
		       std::vector<int> &MVD_indices);

/**
 */
void ComputeMVD_OrientationAttr(CoordGrid &hitMap);

/**
 *
 */
std::vector<TrackObject*>* MergeConnectedComponentsWithMVD(CoordGrid &hitMap,
							   std::vector< std::set<int>* >* STTcompList = 0);

/**
 * Determine the track z coordinates normalized by the distance
 * between tubes. In z direction the step size is constant as function
 * of the distance between the points (nodes). Again hier the virtual
 * tubes are not modified.
 */
void TrackZ_CoordinatesDistNorm(CoordGrid &hitMap, std::vector<TrackObject*>* TrackList = 0);

void SplitSharedNodes(CoordGrid &hitMap, size_t threshold,
		      std::vector < GridNode > &outNodes);
#endif// END of interface
