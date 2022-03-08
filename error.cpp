#include <iostream>
#include <set>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <vector>

#include "error.h"
#include "logc.h"
#include "gridNode.h"
#include "auxfunctions.h"


/* Find Complex tracks based on criterion defined in paper (very brutal approach but who cares?) */

void ComplexTracks(CoordGrid &gr,
		   std::vector< MCTrackObject* >  const *MCTracks,
		   std::vector< int > *ListIDComplex)
{
  std::vector< GridNode > &Ingrid = gr.m_grid;
  
  for(size_t j = 0; j < MCTracks->size(); ++j) {
    MCTrackObject const *MCtrack = MCTracks->at(j);           // Current MC track
    std::vector<int> MCSttComp(MCtrack->m_STT_Component);     // Stt component of the current MC track.

    for(size_t k = 0; k < MCSttComp.size(); k++){
      int CurID         = MCSttComp[k];
      GridNode &CurNode =  Ingrid[CurID-1];
      
      for (size_t i = j+1; i <  MCTracks->size(); i++){
	MCTrackObject const *MCtrack2 = MCTracks->at(i);     // Cross-checking all MC tracks
	std::vector<int> MCSttComp2(MCtrack2->m_STT_Component);
	
	for(size_t l = 0; l < MCSttComp2.size(); l++){
	  int CrossID = MCSttComp2[l];

	  // Check if node is shared between two tracks
	  if(CrossID == CurID){ 
	    
	    if(std::find(ListIDComplex->begin(), ListIDComplex->end(), MCtrack->m_trackID)==ListIDComplex->end()){
	      ListIDComplex->push_back(MCtrack->m_trackID);
	    }
	    
	    if(std::find(ListIDComplex->begin(), ListIDComplex->end(),MCtrack2->m_trackID)==ListIDComplex->end()){
	      ListIDComplex->push_back(MCtrack2->m_trackID);
	    }
	    
	    break;
	  }
	  
	  // Check if node is neighbor to other track
	  if(CurNode.IsNeighboring(CrossID)){ 

	    if(std::find(ListIDComplex->begin(), ListIDComplex->end(),MCtrack->m_trackID)==ListIDComplex->end()){
	      ListIDComplex->push_back(MCtrack->m_trackID);
	    }
	    
	    if(std::find(ListIDComplex->begin(), ListIDComplex->end(),MCtrack2->m_trackID)==ListIDComplex->end()){
	      ListIDComplex->push_back(MCtrack2->m_trackID);
	    }
	    
	    break;
	  }
	}
      }

      // Check if there is a gap such that node i is not neighbors of i-1
      if(k > 0 && !CurNode.IsNeighboring(MCSttComp[k-1])){ 

	if(std::find(ListIDComplex->begin(), ListIDComplex->end(),MCtrack->m_trackID)==ListIDComplex->end()){
	  ListIDComplex->push_back(MCtrack->m_trackID);
	}
	
	break;
      }
    }
  }
  
  return;
  
}


/* Function to compute global event errors */
EvtErrorStruct* ComputeGlobalEvtErrors(CoordGrid &gr,
				       std::vector< MCTrackObject* >  const *MCTracks,
				       std::vector< std::set<int>* >  const *RecoTracks)
{

  dbgtrkerror("Computing event errors based on Babai+16 paper");
    
  int nMC = 0;

  if( (MCTracks == 0) || (RecoTracks->size() == 0) ) {
    error("One of the input parameters for matching is empty.");
    return 0;
  }
  
  // Compute the total MC area (true area) for the current event.
  float TotalArea = 0;
  float Error_underMerge = 0;
  float Error_overMerge  = 0;
  float TotalError = 0;
  
  // Per error type normalized values
  float Error_underMergeNorm = 0;
  float Error_overMergeNorm  = 0;
  float TotalErrorNorm = 0;
  
  
  // Determine the total area of the all tracks in the current event.
  for(size_t i = 0; i < MCTracks->size(); ++i) {
    MCTrackObject const *MCtrack = MCTracks->at(i);
    std::set<int> MCSttComp((MCtrack->m_STT_Component).begin(), (MCtrack->m_STT_Component).end());
    
    if(MCSttComp.size()>5){ // Only considering MC tracks with >5 hits
      TotalArea += (MCtrack->m_STT_Component).size();
      nMC++;
    }
  }

  // Find for each reco'd track the best matched MC track
  std::vector<int>::iterator it;
  std::set<int>::iterator compIt;
  
  for(size_t i = 0; i < RecoTracks->size(); ++i) {
    std::set<int> const *Cur_comp = RecoTracks->at(i);
    std::set<int> Cur_Comp_list;
    
    for(compIt = Cur_comp->begin(); compIt != Cur_comp->end(); ++compIt) {
      int id = *compIt;
      
      // Only Real stt tubes
      if( (id < gr.firstVirtIdx ) && (id >= 0) ) 
	Cur_Comp_list.insert(id);
      
    }

    //   std::sort(Cur_Comp_list.begin(), Cur_Comp_list.end());

    if(Cur_Comp_list.size() < 5) // Discarding CC with less thhan 5 hits
      continue;

    float MatchValue = 0, UnionValue = 0, MCLength = 0;
    int   IDMatchedMC = -1;

    // Find which MC track has the largest overlap (R_k)
    for(size_t j = 0; j < MCTracks->size(); ++j) {
      MCTrackObject const *MCtrack = MCTracks->at(j);
      std::set<int> MCSttComp((MCtrack->m_STT_Component).begin(), (MCtrack->m_STT_Component).end());
      if(MCSttComp.size()>5){
	//	std::sort(MCSttComp.begin(), MCSttComp.end());
      
	// Intersection result output
	std::vector<int> IntersectionList( (Cur_Comp_list.size() + MCSttComp.size()), 0 );
	it = std::set_intersection( Cur_Comp_list.begin(), Cur_Comp_list.end(),
				    MCSttComp.begin(), MCSttComp.end(),
				    IntersectionList.begin());
	IntersectionList.resize(it - IntersectionList.begin());
      
	if(static_cast<float>(IntersectionList.size()) > MatchValue) {
	  MatchValue = IntersectionList.size();
	  MCLength   = static_cast<float>(MCSttComp.size());

	  std::vector<int> UnionList( (Cur_Comp_list.size() + MCSttComp.size()), 0 );
	  it=std::set_union ( Cur_Comp_list.begin(), Cur_Comp_list.end(),
			      MCSttComp.begin(), MCSttComp.end(),
			      UnionList.begin());
	  UnionList.resize(it - UnionList.begin());
	  UnionValue   = static_cast<float>(UnionList.size());
	  IDMatchedMC  =  j;
	}
      }
    }// MC tracks loop

    
    // We have found ( (T_j intersection R_k) is maximum )
    Error_underMerge     += ( ( ( MCLength - MatchValue ) * MatchValue )/ MCLength );
    Error_underMergeNorm += ( ( ( MCLength - MatchValue ) * MatchValue )/ MCLength ) / TotalArea;

    // (alpha_j = Cur_Comp_list.size())
    Error_overMerge     += (static_cast<float>(Cur_Comp_list.size()) - MatchValue);
    Error_overMergeNorm += (static_cast<float>(Cur_Comp_list.size()) - MatchValue) / TotalArea;
   
  }// Components loop
  

  TotalError = (Error_underMerge * Error_underMerge) +
               (Error_overMerge  * Error_overMerge);
  TotalError = sqrt(TotalError)/TotalArea;
  
  // Sub terms are already normalised
  TotalErrorNorm = sqrt( (Error_underMergeNorm * Error_underMergeNorm) +
                         (Error_overMergeNorm  * Error_overMergeNorm));
  
  dbgtrkerror("UnderMerge = %.3f, OverMerge = %.3f, TotalError = %.3f", Error_underMergeNorm, Error_overMergeNorm, TotalErrorNorm);

  // Prepare output
  EvtErrorStruct *outError = new EvtErrorStruct();
  outError->UnderMerge = Error_underMerge;
  outError->OverMerge  = Error_overMerge;
  outError->TotalError = TotalError;
  // Normalised version
  outError->UnderMergeNorm = Error_underMergeNorm;
  outError->OverMergeNorm  = Error_overMergeNorm;
  outError->TotalErrorNorm = TotalErrorNorm;
  outError->nMCTracks      = nMC;
  outError->nRecoTracks    = RecoTracks->size();

  
  dbgtrkerror("End of global error assessment per event (Babai+16) \n");

  return outError;
}
//__________________ END ComputeGlobalEvtErrors _______________



//___________________________ MatchPerTrackWithMCTracks ____________________
/* Function to evaluate error per track. MC is used as ground truth. */
std::vector< TrackErrorStruct* >* ComputeErrorPerRecoTrack(CoordGrid const &hitMap,
                                                           std::vector < MCTrackObject* > const *MCTracks,
                                                           std::vector < PathCandidate* > const *RecoTracks,
							   std::vector<int >  &ListIDComplex,
							   std::vector<int >  &IDMatchesMCReco)
{

  dbgtrkerror("Computing track errors using metrics defined in 2021 thesis");

  std::vector< GridNode > Ingrid = hitMap.m_grid;  

  if( (MCTracks == 0) || (RecoTracks->size() == 0) ) {
    error("One of the input parameters for matching is empty.");
    return 0;
  }
  
  // Create output parameter
  std::vector< TrackErrorStruct* >* outPutPar = 0;
  outPutPar = new std::vector< TrackErrorStruct* >();
  if( outPutPar == 0) {
    error("Could not allocate memory for output list in MatchPerTrackWithMCTracks.");
    delete outPutPar;
    return 0;
  }
 
  std::set<int>::iterator compIt;
  srand(time(NULL)); // Seed the time

  // MC tracks loop
  for(size_t i = 0; i < MCTracks->size(); ++i) {
    int TrkComplex = 0; // Track complex
    int IncRecZ    = 0; // Z info
    // Current MC track(R_k)
    MCTrackObject const *MCtrack = MCTracks->at(i);
    
    // Stt component of the current MC track.
    std::set<int>        MCSttComp((MCtrack->m_STT_Component).begin(), (MCtrack->m_STT_Component).end());
    std::vector<point3D> MC3DPt((MCtrack->m_pointSTTCoordList).begin(), (MCtrack->m_pointSTTCoordList).end());
    std::vector<point3D> MCmomentum((MCtrack->m_STT_Momentum).begin(), (MCtrack->m_STT_Momentum).end());
    std::vector<float>   isochrone((MCtrack->m_STT_Isochrone).begin(), (MCtrack->m_STT_Isochrone).end());
    std::vector<int>     MCSttCompVect(MCtrack->m_STT_Component);

    // We only care about tracks with > 5 hits
    if( (MCSttCompVect.size() > 5) ) { 
      //If track has been identified as complex
      if(std::find(ListIDComplex.begin(), ListIDComplex.end(), MCtrack->m_trackID) != ListIDComplex.end())
	TrkComplex = 1;

      std::vector<point3D> RecoAnchorsCoords, RecoHitsCoords, MCHitsCoords;

      for(compIt = MCSttComp.begin(); compIt != MCSttComp.end(); ++compIt) {
	int id = *compIt;
	GridNode &node = Ingrid[id-1];
	
	if(node.m_type == GridNode::STT_TYPE_PARA)
	  IncRecZ |= 1;
	else if(node.m_type == GridNode::STT_TYPE_SKEW)
	  IncRecZ |= 2;

	point3D pt(node.m_x,node.m_y,node.m_z);
       	MCHitsCoords.push_back(pt);	
      }

      int    MatchTrackIndex = -1;
      float  MatchLength     = 0;
      float  MCLength        = (float) MCSttComp.size();
      float  MatchValue      = std::numeric_limits<float>::min();
      float  UnionValue      = std::numeric_limits<float>::min();

      // Recover the id of the best reco'd track
      if(std::find(IDMatchesMCReco.begin(), IDMatchesMCReco.end(), MCtrack->m_trackID) != IDMatchesMCReco.end()){
	auto it = find(IDMatchesMCReco.begin(),IDMatchesMCReco.end(),  MCtrack->m_trackID);
	int index = it - IDMatchesMCReco.begin();
	MatchTrackIndex = index;
      }

      if(MatchTrackIndex >= 0) {	
        std::set<int> const *bestMatchComponent = (RecoTracks->at(MatchTrackIndex))->m_memberIdSet;
	std::vector<int> const *bestMemberList  = (RecoTracks->at(MatchTrackIndex))->m_memberList;

	std::vector<GridNode> *anchors = &(RecoTracks->at(MatchTrackIndex))->m_anchors;
	std::vector<double> const &x   = RecoTracks->at(MatchTrackIndex)->m_x;
	std::vector<double> const &y   = RecoTracks->at(MatchTrackIndex)->m_y;
	std::vector<double> const &z   = RecoTracks->at(MatchTrackIndex)->m_z;

        std::vector<int>   bestMatchlist;
        std::vector<float> ListIsochroneR;
	
        TrackErrorStruct *erroObject = new TrackErrorStruct();
	erroObject->isComplex = TrkComplex;

	
	// Storing  coordinates of all tubes + isochrone radius of parallel ones
	for(size_t j = 0; j < bestMemberList->size(); j++){
          int id = bestMemberList->at(j);
	  point3D pt(x[j],y[j],z[j]);
	  RecoHitsCoords.push_back(pt);
	  GridNode &node = Ingrid[id-1];
	  
          if( (id < hitMap.firstVirtIdx ) && (id >= 0) ) {
            bestMatchlist.push_back(id);
	    
	    //Finding isochrone radius
	    auto it = std::find(MCSttCompVect.begin(), MCSttCompVect.end(), node.m_detID);
	    
	    if (it != MCSttCompVect.end() && node.m_type != GridNode::STT_TYPE_SKEW) {
		int index = std::distance(MCSttCompVect.begin(), it);
		ListIsochroneR.push_back(isochrone[index]);
	    }  else
	      ListIsochroneR.push_back(0);
	  } else
	    ListIsochroneR.push_back(0);
	  
        }
	
        std::sort(bestMatchlist.begin(), bestMatchlist.end());
	
	// Storing coordinates of anchors
	for(size_t j = 0; j <  anchors->size(); j ++){
	  GridNode &anc = anchors->at(j);
	  //  float randomX = 0.;//((float) rand()) / (float) RAND_MAX/2.;
	  //  float randomY = 0.;//((float) rand()) / (float) RAND_MAX/2.;
	  float randomX = 0;//((float) rand()) / (float) RAND_MAX - 0.5;
	  float randomY = 0;//((float) rand()) / (float) RAND_MAX - 0.5;
	  point3D pt(anc.m_xDet+randomX, anc.m_yDet+randomY, anc.m_zDet);
	  RecoAnchorsCoords.push_back(pt);
	}

	// Determine the sequence intersection. Only real tubes.
	std::vector<int> IntersectionList( (bestMatchlist.size() + MCSttComp.size()), 0 );
	auto it = std::set_intersection(bestMatchlist.begin(), bestMatchlist.end(),
				   MCSttComp.begin(),     MCSttComp.end(),
				   IntersectionList.begin());
      	IntersectionList.resize(it - IntersectionList.begin());
	MatchValue = static_cast<float>(IntersectionList.size());
	MatchLength = static_cast<float>(bestMatchlist.size());

	// Determine the sequence union. Only real tubes.
	std::vector<int> UnionList( (bestMatchlist.size() + MCSttComp.size()), 0);
        it = std::set_union( bestMatchlist.begin(), bestMatchlist.end(),
			     MCSttComp.begin(),     MCSttComp.end(),
			     UnionList.begin());
	UnionList.resize(it - UnionList.begin());
	UnionValue = static_cast<float>(UnionList.size());


	dbgtrkerror("Reco'd length %.0f - MC length %.0f - Intersection length %.0f", MatchLength, MCLength, MatchValue);
	erroObject->MatchIndex         = static_cast<float>(MatchTrackIndex);
        erroObject->MCTrackLength      = MCLength ;
        erroObject->RecoTrackLength    = MatchLength;
        erroObject->IntersectionLength = MatchValue;
        erroObject->UnionLength        = UnionValue; 

	
        // Determine differences
        std::vector<int> MCdiffCurComp((bestMatchlist.size() + MCSttComp.size()), 0);
        it = std::set_difference( MCSttComp.begin(), MCSttComp.end(),
                                  bestMatchlist.begin(), bestMatchlist.end(),
                                  MCdiffCurComp.begin());
        MCdiffCurComp.resize(it - MCdiffCurComp.begin());
	
        // All x in Tracklet and not in MC
        std::vector<int> CurCompdiffMC((bestMatchlist.size() + MCSttComp.size()), 0);
        it = std::set_difference( bestMatchlist.begin(), bestMatchlist.end(),
                                  MCSttComp.begin(), MCSttComp.end(),
                                  CurCompdiffMC.begin());
        CurCompdiffMC.resize(it - CurCompdiffMC.begin());
	
	// Computing the F1 score
	float TP =  MatchValue;
	float FP =  static_cast<float>(CurCompdiffMC.size());
	float FN =  static_cast<float>(MCdiffCurComp.size());
	float F1 = 2*TP/(2*TP+FP+FN);
        erroObject->F1score  = F1;

	dbgtrkerror("F1 score %.3f", F1);

	
	// Compute coordinates errors
	std::vector<double> RecoMatchHitsX, RecoMatchHitsY, RecoMatchHitsZ;
	
	// Getting coordinates of matched tubes 
	for(size_t j = 0; j < (size_t) MatchValue; j++){
	  int id  = IntersectionList[j];
	  auto it = find((RecoTracks->at(MatchTrackIndex))->m_memberList->begin(),
			 (RecoTracks->at(MatchTrackIndex))->m_memberList->end(),
			 id);
	  int index = it - (RecoTracks->at(MatchTrackIndex))->m_memberList->begin();
	  RecoMatchHitsX.push_back(x[index]);
	  RecoMatchHitsY.push_back(y[index]);
	  RecoMatchHitsZ.push_back(z[index]);
	}

       	float DiffHitsCoordsX = 0, DiffHitsCoordsY = 0, DiffHitsCoordsZ = 0;
	float MeanDiffX = 0, MeanDiffY = 0, MeanDiffZ = 0;	  
	std::vector<float> ListDiffHitsCoordsX, ListDiffHitsCoordsY, ListDiffHitsCoordsZ;
	
	for(size_t k = 0; k < RecoMatchHitsX.size(); k++){
	  double Crit  = std::numeric_limits<double>::max();
	  double DiffX = 0, DiffY = 0, DiffZ = 0;

	  for(size_t l = 0; l < MC3DPt.size(); l++){
	    point3D MCHitsCoords = MC3DPt[l];
	    double  Dist = sqrt( pow(RecoMatchHitsX[k] - MCHitsCoords.m_x, 2) +
				 pow(RecoMatchHitsY[k] - MCHitsCoords.m_y, 2) );
	    
	    if(Dist < Crit) {
	      DiffX  = fabs(RecoMatchHitsX[k]-MCHitsCoords.m_x);
	      DiffY  = fabs(RecoMatchHitsY[k]-MCHitsCoords.m_y);
	      DiffZ  = fabs(RecoMatchHitsZ[k]-MCHitsCoords.m_z);
	      
	      DiffHitsCoordsX = (float) (RecoMatchHitsX[k] - MCHitsCoords.m_x);
	      DiffHitsCoordsY = (float) (RecoMatchHitsY[k] - MCHitsCoords.m_y);
	      DiffHitsCoordsZ = (float) (RecoMatchHitsZ[k] - MCHitsCoords.m_z);
	      
	      Crit = Dist;
	    }
	    
	  }
	  
	  MeanDiffX += DiffX;
	  MeanDiffY += DiffY;
	  MeanDiffZ += IncRecZ == 3 ? DiffZ : 200;

	  ListDiffHitsCoordsX.push_back(DiffHitsCoordsX);
	  ListDiffHitsCoordsY.push_back(DiffHitsCoordsY);
	  
	  if(IncRecZ == 3)
	    ListDiffHitsCoordsZ.push_back(DiffHitsCoordsZ);
	  else
	    ListDiffHitsCoordsZ.push_back(200);

	}
	
	MeanDiffX /= (double)RecoMatchHitsX.size();
	MeanDiffY /= (double)RecoMatchHitsY.size();
	MeanDiffZ /= (double)RecoMatchHitsZ.size();

	erroObject->MeanDiffX = MeanDiffX;
	erroObject->MeanDiffY = MeanDiffY;
	erroObject->MeanDiffZ = MeanDiffZ;
	
	erroObject->DiffX = ListDiffHitsCoordsX;
	erroObject->DiffY = ListDiffHitsCoordsY;
	erroObject->DiffZ = ListDiffHitsCoordsZ;

	

	/* Compute Curvature parameters */
	
	std::vector<point3D>  MCSTTpoints = MCtrack->m_pointSTTCoordList;
        CurvatureParameters   MCSTTCurv, MCHitsCurv, RecoAnchorsCurv, RecoHitsCurv;
	
	fit_circle(MCSTTpoints,  MCSTTCurv);
	fit_circle(MCHitsCoords, MCHitsCurv);
	fit_circle(RecoAnchorsCoords, RecoAnchorsCurv);
	fit_circle(RecoHitsCoords,    RecoHitsCurv);

	// Randomize coordinates using isochroone radius?
	CurvatureParameters RecoIsoCurv;
	std::vector<float> ListFitR;
	
	for(int j = 0; j < 100; j++){
	  RecoHitsCoords.clear();
	  
	  for(size_t k = 0; k < ListIsochroneR.size(); k++) {
	    float rx = 0, ry = 0;
	    
	    if(ListIsochroneR[k] >0){
	      float ri = ListIsochroneR[k];
	      rx = ((float) rand()) / (float) RAND_MAX * ri * 2 -ri;
	      ry =sqrt(pow(ri,2) - pow(rx,2));
	      float randNeg = ((float) rand()) / (float) RAND_MAX - 0.5;
	      // ry = ((float) rand()) / (float) RAND_MAX * ry * 2 -ry;
	      ry = randNeg < 0? -1*ry:ry;	
	    }
	    
	    point3D pt(x[k]+rx,y[k]+ry,z[k]);
	    RecoHitsCoords.push_back(pt);
	  }
	  fit_circle(RecoHitsCoords, RecoIsoCurv);
	  ListFitR.push_back(RecoIsoCurv.m_r);
	}
	
	sort(ListFitR.begin(), ListFitR.end() );
	
	if(F1 == 1 && MCLength >= 20 && 0.3*2.*RecoAnchorsCurv.m_r*0.01>2){
	  error("MCtrack %d, MC pt %f, real %f", i, 0.3*2.*MCSTTCurv.m_r*0.01,
		sqrt(pow(MCmomentum[0].m_x,2) +pow(MCmomentum[0].m_y,2))*2);
	  error("Anc rec pt %f",0.3*2.*1./RecoAnchorsCurv.m_r*0.01);
	  error("CC rec pt %f, %f, %f\n", 0.3*2.*0.01*ListFitR[ListFitR.size()/2],
		0.3*2.*0.01*ListFitR[16]-0.3*2.*0.01*ListFitR[ListFitR.size()/2],
		0.3*2.*0.01*ListFitR[ListFitR.size()/2]-0.3*2.*0.01*ListFitR[84]);
	}
	
    
        // Store curvature parameters. MC
	erroObject->MC_px = (double) MCmomentum[0].m_x;
	erroObject->MC_py = (double) MCmomentum[0].m_y;
	erroObject->MC_pz = (double) MCmomentum[0].m_z;
        erroObject->MC_a  = MCSTTCurv.m_a;
        erroObject->MC_b  = MCSTTCurv.m_b;
        erroObject->MC_r1 = MCSTTCurv.m_r;
        erroObject->MC_r2 = MCHitsCurv.m_r;
		
        // Matched tracklet
	erroObject->tr_a    = RecoAnchorsCurv.m_a;
	erroObject->tr_b    = RecoAnchorsCurv.m_b;
	erroObject->tr_rAnc = RecoAnchorsCurv.m_r;
	erroObject->tr_rPts = RecoHitsCurv.m_r;
	erroObject->tr_rIsoRand   = ListFitR[ListFitR.size()/2];
        erroObject->tr_rIsoRand16 = ListFitR[84]-ListFitR[ListFitR.size()/2];
        erroObject->tr_rIsoRand84 = ListFitR[ListFitR.size()/2]-ListFitR[16];
	erroObject->tr_scattAngle = (RecoTracks->at(MatchTrackIndex))->m_scattAngle;// is 1

        // Add to output list
        outPutPar->push_back(erroObject);
      }
      else{// not (MatchTrackIndex >= 0)
        /* There was not a matching tracklet. Totally mis matched.*/	
        TrackErrorStruct *noMatchError = new TrackErrorStruct();
        noMatchError->isNotmatched       = 1.0;
        noMatchError->MCTrackLength      = MCLength;
        noMatchError->RecoTrackLength    = 0;
        noMatchError->IntersectionLength = 0;
        noMatchError->UnionLength        = 0; 
        noMatchError->MatchIndex         = -1; 

	noMatchError->MeanDiffX = -1 ;
	noMatchError->MeanDiffY = -1 ;
	noMatchError->MeanDiffZ = -1 ;
	noMatchError->isComplex = TrkComplex;
	noMatchError->F1score  = 0;// Not Computed here
	
        outPutPar->push_back(noMatchError);
      }
    }//END if(MCSttComp.size() != 0) 
  }//MC tracks loop R_j
  
  dbgtrkerror("End of track error assessment based on 2021 PhD thesis \n");

  return outPutPar;
}
//_______________________ END ComputeErrorPerRecoTrack ___________________




//___________________________ PandaErrorMetric ____________________
std::vector< TrackErrorStruct* >* PandaErrorMetric(CoordGrid const &hitMap,
						  std::vector < MCTrackObject* > const *MCTracks,
						  std::vector < PathCandidate* > const *RecoTracks)
{
  dbgtrkerror("Computing track errors based on PANDA QA");

  if( (MCTracks == 0) || (RecoTracks->size() == 0) ) {
    error("PANDA QA: One of the input parameters for matching is empty.");
    return 0;
  }
  
  // Create output parameter
  std::vector< TrackErrorStruct* >* outPutPar = 0;
  outPutPar = new std::vector< TrackErrorStruct* >();
  if( outPutPar == 0) {
    error("Could not allocate memory for output list in MatchPerTrackWithMCTracks.");
    delete outPutPar;
    return 0;
  }
 

  std::set<int>::iterator compIt;

  int *MCfound      = (int*) calloc(MCTracks->size(), sizeof(int));
  int *MCBestID     = (int*) calloc(MCTracks->size(), sizeof(int));
  int *MCBestRank   = (int*) calloc(MCTracks->size(), sizeof(int));
  int *RecoRank     = (int*) calloc(RecoTracks->size(), sizeof(int));
  int *RecoIDMatch  = (int*) calloc(RecoTracks->size(), sizeof(int));

  int FullyPure = 0, FullyImpure = 0, PartiallyPure = 0, PartiallyImpure = 0, Ghosts = 0, Clones = 0;

  for(size_t i = 0; i < RecoTracks->size(); ++i) {
    int   TrkRank = 0, isTrackClone = 0;
    float MatchValue = 0,  UnionValue = 0, MCLength = 0;

    // Current connected component to analyse
    std::set<int> const *Cur_comp = (RecoTracks->at(i))->m_memberIdSet;
    std::vector<int> Cur_Comp_list;
    
    for(compIt = Cur_comp->begin(); compIt != Cur_comp->end(); ++compIt) {
      int id = *compIt;
      // Only Real stt tubes
      if( (id < hitMap.firstVirtIdx ) && (id >= 0) ) {
	Cur_Comp_list.push_back(id);
      }
    }
    
    if(Cur_Comp_list.size() < 5){
      RecoIDMatch[i] = -1;
      continue;
    }
    
    float RecoLength = Cur_Comp_list.size();
    std::sort(Cur_Comp_list.begin(), Cur_Comp_list.end());
    
    // Find which MC track has the largest overlap
    int IDMatchedMC = -1;
    for(size_t j = 0; j < MCTracks->size(); ++j) {
      // Current MC track
      MCTrackObject const *MCtrack = MCTracks->at(j);
      std::set<int> MCSttComp((MCtrack->m_STT_Component).begin(), (MCtrack->m_STT_Component).end());
      
      // Intersection result output
      std::vector<int> IntersectionList( (Cur_Comp_list.size() + MCSttComp.size()), 0 );
      auto it = std::set_intersection( Cur_Comp_list.begin(), Cur_Comp_list.end(),
				  MCSttComp.begin(), MCSttComp.end(),
				  IntersectionList.begin());
      IntersectionList.resize(it - IntersectionList.begin());
      
      if(static_cast<float>(IntersectionList.size()) > MatchValue) {
	MatchValue = IntersectionList.size();
	// A_k
	MCLength = static_cast<float>(MCSttComp.size());

	std::vector<int> UnionList( (Cur_Comp_list.size() + MCSttComp.size()), 0 );
	it=std::set_union ( Cur_Comp_list.begin(), Cur_Comp_list.end(),
				  MCSttComp.begin(), MCSttComp.end(),
			     UnionList.begin());
	UnionList.resize(it - UnionList.begin());
	UnionValue = static_cast<float>(UnionList.size());
	IDMatchedMC  =  j;
      }
    }// MC tracks loop

    dbgtrkerror("Matched track %d with MC %d", i, IDMatchedMC);
    dbgtrkerror("MatchValue = %.0f, MCLength = %.0f, RecoLength %.0f", MatchValue,MCLength, RecoLength);
    RecoIDMatch[i] = IDMatchedMC;
    
    if(MatchValue == MCLength  && MatchValue == RecoLength){
      //dbgtrkerror("Track %d is fully pure", i);
      FullyPure++;
      TrkRank = 1;
    } else if(MatchValue == MCLength && MatchValue/RecoLength >= 0.7){
      //dbgtrkerror("Track %d is fully impure", i);
      FullyImpure++;
      TrkRank = 2;
    } else if(MatchValue == RecoLength && MatchValue/MCLength > 0.7){
      // dbgtrkerror("Track %d is partially pure", i);
      PartiallyPure++;
      TrkRank = 3;
    } else if(MatchValue/RecoLength >= 0.7){
      // dbgtrkerror("Track %d is partially impure", i);
      PartiallyImpure++;
      TrkRank = 4;
    } else if(MatchValue/RecoLength <0.7){
      //dbgtrkerror("Track %d is a ghost", i);
      Ghosts++;
      TrkRank = 5;
    } else
      error("Not assigned");

    RecoRank[i]= TrkRank;
    
    /*  if(MCfound[IDMatchedMC]>=1){
      dbgtrkerror("MC track %d already matched, we might have a clone here!");
      isTrackClone = 1;
      Clones++;
      }*/
    
    MCfound[IDMatchedMC]++;
    if(TrkRank < MCBestRank[IDMatchedMC] || MCBestRank[IDMatchedMC] == 0){
      MCBestRank[IDMatchedMC] =  TrkRank;
      MCBestID[IDMatchedMC]   = i;
    }

    dbgtrkerror(" ");

    
  }


    
  for(size_t j = 0; j < MCTracks->size(); ++j) {
    //dbgtrkerror("Track %d has been found %d times", j, MCfound[j]);
    if(MCfound[j]>1){
      Clones += MCfound[j]-1;
      for(size_t i = 0; i < RecoTracks->size(); ++i) {
	if( RecoIDMatch[i] == j && i != MCBestID[j] && RecoRank[i] != 5){
	  // dbgtrkerror("Track %d updated to Clone", i);
	  RecoRank[i]= 6;
	}
      }
    }
  }

  for(size_t i = 0; i < RecoTracks->size(); ++i) {
    if( RecoIDMatch[i] != -1){
      TrackErrorStruct *erroObject = new TrackErrorStruct();
      erroObject->tr_rank = RecoRank[i];
      dbgtrkerror("Track %d rank is %d", i, RecoRank[i]);

      //erroObject->tr_isClone = isTrackClone;
      outPutPar->push_back(erroObject);
    }
  }
  
  dbgtrkerror("End of PANDA Error metric function \n");
  delete(MCfound);
  delete(MCBestRank);
  delete(MCBestID);
  delete(RecoRank);
  delete(RecoIDMatch);
    

  return outPutPar;
}
//__________________


std::vector< int > MatchBestRecoToMC( CoordGrid const &hitMap,
				      std::vector < MCTrackObject* > const *MCTracks,
				      std::vector < PathCandidate* > const *RecoTracks)
{

  std::vector<int> IDMatchesMCReco;
  for(size_t i = 0;i< RecoTracks->size(); i++)
    IDMatchesMCReco.push_back(-1);

  std::vector<int> matchedCCId;
 
  if( (MCTracks == 0) || (RecoTracks->size() == 0) ) {
    error("One of the input parameters for matching is empty.");
    return IDMatchesMCReco;
  }

  using trMatch = std::pair<int, int>;
  std::vector<std::vector<trMatch>> TrkMatches(MCTracks->size());
  std::set<int>::iterator compIt;

  // MC tracks loop
  for(size_t i = 0; i < MCTracks->size(); ++i) {
    MCTrackObject const *MCtrack = MCTracks->at(i);
    std::set<int> MCSttComp((MCtrack->m_STT_Component).begin(), (MCtrack->m_STT_Component).end());
    std::vector<int> MCSttCompVect(MCtrack->m_STT_Component);

    if( MCSttCompVect.size() > 5 ) {
 
      int    MatchTrackIndex = -1;
      float  MatchLength  = 0;
      float  MCLength     = (float) MCSttComp.size();
      float  MatchValue   = 0;

      for(size_t j = 0; j < RecoTracks->size(); ++j) {
          std::set<int> const *Cur_comp = (RecoTracks->at(j))->m_memberIdSet;
          std::set<int> Cur_Comp_list;

          for(compIt = Cur_comp->begin(); compIt != Cur_comp->end(); ++compIt) {
            int id = *compIt;
            if( (id < hitMap.firstVirtIdx ) && (id >= 0) ) {
              Cur_Comp_list.insert(id);
            }
          }//

          std::vector<int> IntersectionList( (Cur_Comp_list.size() + MCSttComp.size()), 0 );
	  auto it = std::set_intersection( Cur_Comp_list.begin(), Cur_Comp_list.end(),
                                      MCSttComp.begin(), MCSttComp.end(),
                                      IntersectionList.begin());
          IntersectionList.resize(it - IntersectionList.begin());

	  TrkMatches[i].push_back(make_pair(j, (int) IntersectionList.size()));
	
           if(static_cast<float>(IntersectionList.size()) > MatchValue) {
            MatchValue = static_cast<float>(IntersectionList.size());
            MatchLength = static_cast<float>(Cur_Comp_list.size());
            MatchTrackIndex = j;
          }
	  
      }// Components loop (T_j)

      int matchIdx = -1;
      sort(TrkMatches[i].begin(), TrkMatches[i].end(), sortbysec2);
      
      for(size_t p = 0; p < RecoTracks->size(); p++){
	if(std::find(matchedCCId.begin(), matchedCCId.end(),
		     TrkMatches[i][p].first) == matchedCCId.end() && TrkMatches[i][p].second != 0 ){
	  IDMatchesMCReco[TrkMatches[i][p].first] = MCtrack->m_trackID;
	  matchIdx = TrkMatches[i][p].first;
	  matchedCCId.push_back(TrkMatches[i][p].first);
	  break;
	}
      }
     
      dbgtrkerror("MC track %d has best match index %d",MCtrack->m_trackID, matchIdx);
      
    }//END if(MCSttComp.size() != 0) 
    //Go to next MC-Track
  }//MC tracks loop R_j
  // Debug info before return
  return IDMatchesMCReco;
}
//____________
