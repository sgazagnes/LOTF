/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#pragma once
#ifndef CONSTANTS_H
#define CONSTANTS_H
//+++++++++++++++++++++++++++++++++++++++++++++
//_________ Collect MVDSTT data
#define INCLUDE_STT_POINTS  1
#define INCLUDE_MVD_POINTS  0
#define DEBUG_PRINT     0
#define WRITE_TO_ASCII_FILE 1
#define HIT_EXCLUSION      -10000.0
#define NOT_AVAILABLE  -10000.0
//+++++++++++++++++++++++++++++++++++++++++++++
//___________________ Coordinates grid
#define HIT_EXCLUSION      -10000.0
#define NOT_AVAILABLE  -10000.0
//+++++++++++++++++++++++++++++++++++++++++++++
//______________ Pathopen
/// Local tollerance. This needs to be defined global???
#define SKEWED_VIRTUAL_ANGLE_TOLLERANCE 0.8 // ~60 degrees
// If we can use c++11 with Root, then use almost equal from utility
// function header
#define Local_ConComp_Epsilon 0.0001
/// DEBUG AND DEFINES
#define INCLUDE_LEFTLEFT_RIGHTRIGHT 0 // Hyper connectivity in orient space.
#define PATH_DEBUG_PRINT 1
#define ORIENT_DEBUG_INFO 1
#define MVD_MERGE_DEBUG_PRINT 1
#define Z_DETERMINATION_DEBUG 0
#define CONNECTED_COMPONENT_DEBUG 1
#define ADD_SHORTQUEUE_TO_TAIL 1
#define REORDER_AND_NORM_NODES_FOR_Z 1
#define LAYERBASED_DEBUG_PRINTS 1
//0: Outer->inner, 1: Inner->outer
#define LAYERBASED_OUTER_TO_INNER 0
//+++++++++++++++++++++++++++++++++++++++++++++
//_____________ PerformFilter___________
#define EXCLUDE_STTSKEWED_PLOT   0
#define EXCLUDE_VIRTUALS_IN_PLOT 0

#define WRITE_CONNECTED_COMPONENTS 1
#define EVALUATE_ERROR 1
#define PRINT_DEBUG_INFO 0
#define INCLUDE_MVD_INOUTPUT_TRACK 1
#define PRINT_TRACKS_DEBUG_INFO 1
#define PRINT_DEBUG_INFO_COMP_MATCH 0
//+++++++++++++++++++++++++++++++++++++++++++++
//_______________ Used in utilityfunctions
#define LOCAL_ZERO_COMPARE_EPSILON 0.00000001
#define LOCAL_PND_TRACKING_EPSILON 0.001
#define CURVATURE_TOLERANCE  0.2
#define CAND_MERGE_DISTANCE_SELECTION_TOLERANCE 15.0
#define CAND_MERGE_DISTANCE_TOLERANCE 6.0
#define SHORT_CANDIDATES_MERGE_DISTANCE 5.00
//+++++++++++++++++++++++++++++++++++++++++++++
//___________________ Defined in auxiliaryfunctions.h ____
//Start id for virtual nodes
#define START_VIRTUAL_ID 6000

#endif
